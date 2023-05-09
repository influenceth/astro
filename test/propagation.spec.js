import { expect } from 'chai';
import almostEqual from 'almost-equal';
import * as math from 'mathjs';
import elements from '../src/elements.js';
import propagation from '../src/propagation.js';

// Test fixtures [ k, p, ecc, inc, raan, argp, nu ]
const elliptical = [ 1.1368823e11 , 2.9328214e8, 0.325, 0.002443461, 3.410897, 5.2838098, 0.94806285 ];
const hyperbolic = [ 398600.44, 48848.56334147761, 1.7311, 0.122138, 1.00681, 3.10686, 0.12741601769795755 ];
const parabolic_low = [ 398600.44, 265608.0, 0.995, 0.122138, 1.00681, 3.10686, 0.5 ];
const parabolic_high = [ 398600.44, 265608.0, 1.005, 0.122138, 1.00681, 3.10686, 1.0 ];

describe('Propagation with Farnocchia', function () {
  it('should calculate time past periapsis for elliptical', function () {
    const nu = elliptical[6];
    const ecc = elliptical[2];
    const mu = elliptical[0];
    const q = elliptical[1] / (1 + ecc);
    const t0 = propagation.delta_t_from_nu(nu, ecc, mu, q);
    expect(almostEqual(t0, 8667649.6863514, 0, 1e-7)).to.be.true;
  });

  it('should calculate time past periapsis for hyperbolic', function () {
    const nu = hyperbolic[6];
    const ecc = hyperbolic[2];
    const mu = hyperbolic[0];
    const q = hyperbolic[1] / (1 + ecc);
    const t0 = propagation.delta_t_from_nu(nu, ecc, mu, q);
    expect(almostEqual(t0, 293.1233793, 0, 1e-7)).to.be.true;
  });

  it('should calculate time past periapsis for close to parabolic (low)', function () {
    const nu = parabolic_low[6];
    const ecc = parabolic_low[2];
    const mu = parabolic_low[0];
    const q = parabolic_low[1] / (1 + ecc);
    const t0 = propagation.delta_t_from_nu(nu, ecc, mu, q);
    expect(almostEqual(t0, 28421.60085025905, 0, 1e-7)).to.be.true;
  });

  it('should calculate time past periapsis for close to parabolic (high)', function () {
    const nu = parabolic_high[6];
    const ecc = parabolic_high[2];
    const mu = parabolic_high[0];
    const q = parabolic_high[1] / (1 + ecc);
    const t0 = propagation.delta_t_from_nu(nu, ecc, mu, q);
    expect(almostEqual(t0, 64825.61097546987, 0, 1e-7)).to.be.true;
  });

  it('should propogate one full period', function () {
    const a = elliptical[1] / (1 - elliptical[2] ** 2);
    const period = 2 * Math.PI * Math.sqrt(a ** 3 / elliptical[0]);
    const nu = propagation.farnocchia_coe(...elliptical, period);
    expect(almostEqual(nu, elliptical[6], 0, 1e-7)).to.be.true;
  });

  it('should propogate with rv vectors', function () {
    const [ r0, v0 ] = elements.coe2rv(...elliptical);
    const a = elliptical[1] / (1 - elliptical[2] ** 2);
    const period = 2 * Math.PI * Math.sqrt(a ** 3 / elliptical[0]);
    const [ r, v ] = propagation.farnocchia_rv(elliptical[0], r0, v0, period);
    expect(almostEqual(math.norm(r0), math.norm(r), 0, 1e-7)).to.be.true;
    expect(almostEqual(math.norm(v0), math.norm(v), 0, 1e-7)).to.be.true;
  });
});
