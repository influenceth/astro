import { expect } from 'chai';
import almostEqual from 'almost-equal';
import Orbit from '../src/orbit.js';

const elliptical = [ 1.1368823e11 , 2.9328214e8, 0.325, 0.002443461, 3.410897, 5.2838098, 0.94806285 ];

describe('Orbit model', function () {
  it('should be able to create an orbit from classical elements', function () {
    const orbit = Orbit.fromClassicElements(...elliptical);
    expect(orbit).to.be.an.instanceof(Orbit);
  });

  it('should be able to create an orbit from state vectors', function () {
    let rv = Orbit.fromClassicElements(...elliptical).rv;
    const orbit = Orbit.fromStateVectors(elliptical[0], rv.r, rv.v);
    expect(orbit).to.be.an.instanceof(Orbit);
  });

  it('should return ephem samples', function () {
    const orbit = Orbit.fromClassicElements(...elliptical);
    const samples = orbit.ephem(10);
    expect(samples).to.be.an('array');
    expect(samples.length).to.equal(10);
  });

  it('should return the semi major axis', function () {
    const orbit = Orbit.fromClassicElements(...elliptical);
    expect(almostEqual(orbit.a, 3.2791853e8, 0, 1e-7)).to.be.true;
  });

  it('should return the period', function () {
    const orbit = Orbit.fromClassicElements(...elliptical);
    expect(almostEqual(orbit.period, 1.10655e8, 0, 1e-7)).to.be.true;
  });

  it('should propagate by tof', function () {
    const orbit = Orbit.fromClassicElements(...elliptical);
    orbit.propagateFor(orbit.period / 2);
    expect(almostEqual(orbit.epoch, orbit.period / 2, 0, 1e-7)).to.be.true;
  });
});
