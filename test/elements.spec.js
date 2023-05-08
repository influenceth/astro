import { expect } from 'chai';
import elements from '../src/elements.js';

// Test fixtures
const elliptical = [ 398600.4418, 11067.790, 0.83285, 1.5336208137, 3.9774308324, 0.9316567547, 1.6115497648 ];
const circular = [ 3.9860047e14, 24464560.0, 0.0, 0.122138, 1.00681, 0.0, 0.048363 ];
const hyperbolic = [ 3.9860047e14, 4.884856334147761e7, 1.7311, 0.122138, 1.00681, 3.10686, 0.12741601769795755 ];
const equatorial = [ 3.9860047e14, 1.13880762905224e7, 0.7311, 0.0, 0.0, 3.10686, 0.44369564302687126 ];
const circular_equatorial = [ 3.9860047e14, 1.13880762905224e7, 0.0, 0.0, 0.0, 0.0, 0.44369564302687126 ];

describe('Elements conversions', function () {
  it('should generate classical elements from state vectors', function () {
    const r = [ 6.52536812e3,  6.86153183e3,  6.44911861e3];
    const v = [ 4.90227865e0,  5.53313957e0, -1.97571010e0];
    const { p, ecc, inc, raan, argp, nu } = elements.rv2coe(elliptical[0], r, v);

    expect(p.toFixed(3)).to.equal((elliptical[1]).toFixed(3));
    expect(ecc.toFixed(7)).to.equal((elliptical[2]).toFixed(7));
    expect(inc.toFixed(7)).to.equal((elliptical[3]).toFixed(7));
    expect(raan.toFixed(7)).to.equal((elliptical[4]).toFixed(7));
    expect(argp.toFixed(7)).to.equal((elliptical[5]).toFixed(7));
    expect(nu.toFixed(7)).to.equal((elliptical[6]).toFixed(7));
  });

  it('should generate classical elements from state vectors for hyperbolic orbits', function () {
    const r = [ -8173532.704990985, -16011738.121722244, -202722.41788328032 ];
    const v = [ 6580.954362290776, -4052.9254479420347, -948.6406516239671 ];
    const { p, ecc, inc, raan, argp, nu } = elements.rv2coe(hyperbolic[0], r, v);

    expect(p.toFixed(3)).to.equal((hyperbolic[1]).toFixed(3));
    expect(ecc.toFixed(7)).to.equal((hyperbolic[2]).toFixed(7));
    expect(inc.toFixed(7)).to.equal((hyperbolic[3]).toFixed(7));
    expect(raan.toFixed(7)).to.equal((hyperbolic[4]).toFixed(7));
    expect(argp.toFixed(7)).to.equal((hyperbolic[5]).toFixed(7));
    expect(nu.toFixed(7)).to.equal((hyperbolic[6]).toFixed(7));
  });

  it('should generate state vectors from classical elements', function () {
    const [ r, v ] = elements.coe2rv(...elliptical);

    expect(r[0].toFixed(5)).to.equal((6.52536812e3).toFixed(5));
    expect(r[1].toFixed(5)).to.equal((6.86153183e3).toFixed(5));
    expect(r[2].toFixed(5)).to.equal((6.44911861e3).toFixed(5));
    expect(v[0].toFixed(5)).to.equal((4.90227865e0).toFixed(5));
    expect(v[1].toFixed(5)).to.equal((5.53313957e0).toFixed(5));
    expect(v[2].toFixed(5)).to.equal((-1.97571010e0).toFixed(5));
  });

  it('should generate state vectors from classical elements for hyperbolic orbits', function () {
    const [ r, v ] = elements.coe2rv(...hyperbolic);

    expect(r[0].toFixed(5)).to.equal((-8173532.70499).toFixed(5));
    expect(r[1].toFixed(5)).to.equal((-16011738.12172).toFixed(5));
    expect(r[2].toFixed(5)).to.equal((-202722.41788).toFixed(5));
    expect(v[0].toFixed(5)).to.equal((6580.95436).toFixed(5));
    expect(v[1].toFixed(5)).to.equal((-4052.92545).toFixed(5));
    expect(v[2].toFixed(5)).to.equal((-948.64065).toFixed(5));
  });
});
