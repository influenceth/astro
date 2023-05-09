import * as math from 'mathjs';
import { E_to_nu, F_to_nu } from '../src/angles.js';
import { modulo, rotation_matrix } from '../src/utils.js';

/**
 * Converts from classical orbital elements to state vectors
 *
 * @param mu Standard gravitational parameter (km^3 / s^2)
 * @param p Semi latus rectum or parameter
 * @param ecc Eccentricity
 * @param inc Inclination (rad)
 * @param raan Longitude of ascending node (rad)
 * @param argp Argument of periapsis (rad)
 * @param nu True anomaly (rad)
 */
export const coe2rv = (mu, p, ecc, inc, raan, argp, nu) => {
  const pqw = math.dotMultiply(
    [[ Math.cos(nu), Math.sin(nu), 0 ], [ -Math.sin(nu), ecc + Math.cos(nu), 0 ]],
    [[0, 0, 0].fill(p / (1 + ecc * Math.cos(nu)), 0, 3), [0, 0, 0].fill(Math.sqrt(mu / p), 0, 3)]
  );

  let rm = rotation_matrix(raan, 2);
  rm = math.multiply(rm, rotation_matrix(inc, 0));
  rm = math.multiply(rm, rotation_matrix(argp, 2));

  return math.multiply(pqw, math.transpose(rm));
};

/**
 * Converts from state vectors to classical orbital elements
 *
 * @param mu Standard gravitational parameter (km^3 / s^2)
 * @param r Position vector (km)
 * @param v Velocity vector (km / s)
 * @param tol Tolerance for eccentricity and inclination checks
 * @returns {Object} Classical orbital elements { p, ecc, inc, raan, argp, nu }
 */
export const rv2coe = (mu, r, v, tol = 1e-8) => {
  let raan, argp, nu;

  const h = math.cross(r, v);
  const n = math.cross([0, 0, 1], h);
  const e = math.divide(
    math.subtract(
      math.multiply(math.dot(v, v) - mu / math.norm(r), r),
      math.multiply(math.dot(r, v), v)
    ),
    mu
  );

  const ecc = math.norm(e);
  const p = math.dot(h, h) / mu;
  const inc = Math.acos(h[2] / math.norm(h));

  const circular = ecc < tol;
  const equatorial = Math.abs(inc) < tol;

  if (equatorial && !circular) {
    // Equatorial elliptical orbit
    raan = 0;
    argp = Math.atan2(e[1], e[0]) % (2 * Math.PI);
    nu = Math.atan2(math.dot(h, math.cross(e, r)) / math.norm(h), math.dot(r, e));
  } else if (!equatorial && circular) {
    // Non-equatorial circular orbit
    raan = modulo(Math.atan2(n[1], n[0]), (2 * Math.PI));
    argp = 0;
    nu = Math.atan2(math.dot(r, math.cross(h, n)) / math.norm(h), math.dot(r, n));
  } else if (equatorial && circular) {
    // Equatorial circular orbit
    raan = 0;
    argp = 0;
    nu = modulo(Math.atan2(r[1], r[0]), (2 * Math.PI));
  } else {
    const a = p / (1 - (ecc ** 2));
    const mu_a = mu * a;

    if (a > 0) {
      // Elliptical orbit
      const e_se = math.dot(r, v) / Math.sqrt(mu_a);
      const e_ce = math.norm(r) * math.dot(v, v) / mu - 1;
      const E = Math.atan2(e_se, e_ce);
      nu = E_to_nu(E, ecc);
    } else {
      // Hyperbolic orbit
      const e_sh = math.dot(r, v) / Math.sqrt(-mu_a);
      const e_ch = math.norm(r) * (math.norm(v) ** 2) / mu - 1;
      const F = Math.log((e_ch + e_sh) / (e_ch - e_sh)) / 2;
      nu = F_to_nu(F, ecc);
    }

    raan = modulo(Math.atan2(n[1], n[0]), (2 * Math.PI));
    const px = math.dot(r, n);
    const py = math.dot(r, math.cross(h, n)) / math.norm(h);
    argp = modulo((Math.atan2(py, px) - nu), (2 * Math.PI));
  }

  // Shift true anomaly into range of -pi to pi
  nu = ((nu % (2 * Math.PI) + 3 * Math.PI) % (2 * Math.PI)) - Math.PI;

  return { p, ecc, inc, raan, argp, nu };
}

export default {
  coe2rv,
  rv2coe
};
