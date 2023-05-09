import { modulo } from './utils.js';
import { rv2coe, coe2rv } from './elements.js';
import {
  nu_to_E,
  E_to_M,
  M_to_E,
  E_to_nu,
  nu_to_D,
  D_to_M,
  D_to_M_near_parabolic,
  M_to_D,
  M_to_D_near_parabolic,
  D_to_nu,
  nu_to_F,
  F_to_M,
  M_to_F,
  F_to_nu
} from './angles.js';

// See: https://github.com/poliastro/poliastro/blob/main/src/poliastro/core/propagation/farnocchia.py

/**
 * Propagates an orbit using mean motion and returns new position and velocity vectors
 *
 * @param {number} mu Gravitational parameter (km^3 / s^2)
 * @param {number[]} r0 Initial position vector (km)
 * @param {number[]} v0 Initial velocity vector (km / s)
 * @param {number} tof Time of flight (sec)
 * @returns {number[][]} Final position and velocity vectors
 */
export const farnocchia_rv = (mu, r0, v0, tof) => {
  const { p, ecc, inc, raan, argp, nu: nu0 } = rv2coe(mu, r0, v0);
  const nu = farnocchia_coe(mu, p, ecc, inc, raan, argp, nu0, tof);
  return coe2rv(mu, p, ecc, inc, raan, argp, nu);
}

/**
 * Returns the new true anomaly for a given time of flight
 *
 * @param {number} mu Gravitation parameter (km^3 / s^2)
 * @param {number} p Semi-latus rectum (km)
 * @param {number} ecc Eccentricity (rad)
 * @param {number} inc Inclination (rad)
 * @param {number} raan Longitude of ascending node (rad)
 * @param {number} argp Argument of periapsis (rad)
 * @param {number} nu True anomaly (rad)
 * @param {number} tof Time of flight (sec)
 */
export const farnocchia_coe = (mu, p, ecc, inc, raan, argp, nu, tof) => {
  const q = p / (1 + ecc);
  const delta_t0 = delta_t_from_nu(nu, ecc, mu, q);
  const delta_t = delta_t0 + tof;
  return nu_from_delta_t(delta_t, ecc, mu, q);
};

/**
 * Time elapsed since periapsis for given true anomaly.
 *
 * @param {number} nu True anomaly (rad)
 * @param {number} ecc Eccentricity
 * @param {number} mu Gravity parameter
 * @param {number} q Periapsis distance
 * @param {number} delta Parameter that controls the size of the near parabolic region
 * @returns {number} Time elapsed since periapsis
 */
export const delta_t_from_nu = (nu, ecc, mu = 1, q = 1, delta = 1e-2) => {
  if (ecc < 0) throw new Error('ecc must be in [0, ∞)');
  if ( nu >= Math.PI || nu < -Math.PI ) throw new Error('nu must be in [-pi, pi)');
  let M, n;

  if (ecc < 1 - delta) {
    // Strong elliptic
    const E = nu_to_E(nu, ecc); // (-pi, pi]
    M = E_to_M(E, ecc); // (-pi, pi]
    n = Math.sqrt(mu * (1 - ecc) ** 3 / q ** 3);
  } else if (ecc < 1) {
    const E = nu_to_E(nu, ecc); // (-pi, pi]

    if (delta <= 1 - ecc * Math.cos(E)) {
      // Strong elliptic
      M = E_to_M(E, ecc); // (-pi, pi]
      n = Math.sqrt(mu * (1 - ecc) ** 3 / q ** 3)
    } else {
      // Near parabolic
      const D = nu_to_D(nu); // (-∞, ∞)
      // If |nu| is far from pi this result is bounded
      // because the near parabolic region shrinks in its vicinity,
      // otherwise the eccentricity is very close to 1
      // and we are really far away
      M = D_to_M_near_parabolic(D, ecc);
      n = Math.sqrt(mu / (2 * q ** 3));
    }
  } else if (ecc == 1) {
    // Parabolic
    const D = nu_to_D(nu); // (-∞, ∞)
    M = D_to_M(D); // (-∞, ∞)
    n = Math.sqrt(mu / (2 * q ** 3));
  } else if (1 + ecc * Math.cos(nu) < 0) {
    // Unfeasible region
    return NaN;
  } else if (ecc <= 1 + delta) {
    // NOTE: Do we need to wrap nu here?
    // For hyperbolic orbits, it should anyway be in
    // (-arccos(-1 / ecc), +arccos(-1 / ecc))
    const F = nu_to_F(nu, ecc); // (-∞, ∞)

    if (delta <= ecc * Math.cosh(F) - 1) {
      // Strong hyperbolic
      M = F_to_M(F, ecc); // (-∞, ∞)
      n = Math.sqrt(mu * (ecc - 1) ** 3 / q ** 3);
    } else {
      // Near parabolic
      const D = nu_to_D(nu); // (-∞, ∞)
      M = D_to_M_near_parabolic(D, ecc); // (-∞, ∞)
      n = Math.sqrt(mu / (2 * q ** 3));
    }
  } else {
    // Strong hyperbolic
    const F = nu_to_F(nu, ecc); // (-∞, ∞)
    M = F_to_M(F, ecc); // (-∞, ∞)
    n = Math.sqrt(mu * (ecc - 1) ** 3 / q ** 3);
  }

  return M / n;
};

/**
 * True anomaly for given elapsed time since periapsis.
 *
 * @param {number} delta_t TIme elapsed since periapsis
 * @param {number} ecc Eccentricity
 * @param {number} mu Gravity parameter
 * @param {number} q Periapsis distance
 * @param {number} delta Parameter that controls the size of the near parabolic region
 */
export const nu_from_delta_t = (delta_t, ecc, mu = 1, q = 1, delta = 1e-2) => {
  let nu;

  if (ecc < 1 - delta) {
    // Strong elliptic
    const n = Math.sqrt(mu * (1 - ecc) ** 3 / q ** 3);
    const M = n * delta_t;
    // This might represent several revolutions,
    // so we wrap the true anomaly
    const E = M_to_E(modulo((M + Math.PI), (2 * Math.PI)) - Math.PI, ecc);
    nu = E_to_nu(E, ecc);
  } else if (ecc < 1) {
    const E_delta = Math.acos((1 - delta) / ecc);
    // Compute M assuming we are in the strong elliptic case and verify later
    const n = Math.sqrt(mu * (1 - ecc) ** 3 / q ** 3);
    const M = n * delta_t;

    // Check against abs(M) because E_delta could also be negative
    if (E_to_M(E_delta, ecc) <= abs(M)) {
      // Strong elliptic. This might represent several revolutions, so wrap the true anomaly
      const E = M_to_E(modulo((M + Math.PI), (2 * Math.PI)) - Math.PI, ecc);
      nu = E_to_nu(E, ecc);
    } else {
      // Near parabolic, recompute M
      const n = Math.sqrt(mu / (2 * q ** 3));
      const M = n * delta_t;
      const D = M_to_D_near_parabolic(M, ecc);
      nu = D_to_nu(D);
    }
  } else if (ecc == 1) {
    // Parabolic
    const n = Math.sqrt(mu / (2 * q ** 3));
    const M = n * delta_t;
    const D = M_to_D(M);
    nu = D_to_nu(D);
  } else if (ecc <= 1 + delta) {
    const F_delta = Math.acosh((1 + delta) / ecc);
    // Compute M assuming we are in the strong hyperbolic case and verify later
    const n = Math.sqrt(mu * (ecc - 1) ** 3 / q ** 3);
    const M = n * delta_t;

    // Check against abs(M) because F_delta could also be negative
    if (F_to_M(F_delta, ecc) <= abs(M)) {
      // Strong hyperbolic, proceed
      const F = M_to_F(M, ecc);
      nu = F_to_nu(F, ecc);
    } else {
      // Near parabolic, recompute M
      const n = Math.sqrt(mu / (2 * q ** 3));
      const M = n * delta_t;
      const D = M_to_D_near_parabolic(M, ecc);
      nu = D_to_nu(D);
    }
  } else {
    // Strong hyperbolic
    const n = Math.sqrt(mu * (ecc - 1) ** 3 / q ** 3);
    const M = n * delta_t;
    const F = M_to_F(M, ecc);
    nu = F_to_nu(F, ecc);
  }

  return nu;
}

export default {
  delta_t_from_nu,
  nu_from_delta_t,
  farnocchia_coe,
  farnocchia_rv
};
