import { expect } from 'chai';
import almostEqual from 'almost-equal';
import angles from '../src/angles.js';

describe('Angle conversions', function () {
  it('should calculate eccentric anomaly from mean anomaly', function () {
    expect(angles.M_to_E(1, 0.5).toFixed(7)).to.equal('1.4987011');
    expect(angles.M_to_E(2, 0.25).toFixed(7)).to.equal('2.2018514');
    expect(angles.M_to_E(3, 0.125).toFixed(7)).to.equal('3.0156956');
    expect(angles.M_to_E(4, 0.9).toFixed(7)).to.equal('3.6009583');
    expect(angles.M_to_E(5, 0.99).toFixed(7)).to.equal('4.1581841');
    expect(angles.M_to_E(-4, 0.5).toFixed(7)).to.equal('-3.7246928');
  });

  it('should calculate true eccentric anomaly from mean anomaly', function () {
    expect(angles.M_to_F(1, 1.1).toFixed(7)).to.equal('1.5928117');
    expect(angles.M_to_F(2, 1.5).toFixed(7)).to.equal('1.6126858');
    expect(angles.M_to_F(3, 2.0).toFixed(7)).to.equal('1.5628462');
    expect(angles.M_to_F(4, 5.0).toFixed(7)).to.equal('0.8616757');
    expect(angles.M_to_F(-4, 5.0).toFixed(7)).to.equal('-0.8616757');
  });

  it('should calculate true anomaly from eccentric anomaly', function () {
    expect(angles.E_to_nu(-4, 0.125).toFixed(7)).to.equal('2.3743532');
    expect(angles.E_to_nu(-2, 0.5).toFixed(7)).to.equal('-2.4315800');
    expect(angles.E_to_nu(0, 0.75).toFixed(7)).to.equal('0.0000000');
    expect(angles.E_to_nu(2, 0.9).toFixed(7)).to.equal('2.8490840');
    expect(angles.E_to_nu(4, 0.99).toFixed(7)).to.equal('-3.0767304');
  });

  it('should calculate true anomaly from hyperbolic eccentric anomaly', function () {
    expect(angles.F_to_nu(-5, 1.125).toFixed(7)).to.equal('-2.6594996');
    expect(angles.F_to_nu(-3, 1.5).toFixed(7)).to.equal('-2.2237955');
    expect(angles.F_to_nu(-1, 2.75).toFixed(7)).to.equal('-1.1895181');
    expect(angles.F_to_nu(1, 2.9).toFixed(7)).to.equal('1.1696349');
    expect(angles.F_to_nu(3, 3.99).toFixed(7)).to.equal('1.7265870');
  });

  it('should calculate true anomaly from mean anomaly', function () {
    // Elliptical
    expect(angles.M_to_nu(0.9480628496833199, 0.325).toFixed(7)).to.equal('1.5891713');
    expect(angles.M_to_nu(2.042558823608964, 0.05).toFixed(7)).to.equal('2.1290652');
    expect(angles.M_to_nu(27, 0.5).toFixed(7)).to.equal('2.6063312');

    // Hyperbolic
    expect(angles.M_to_nu(27, 2.5).toFixed(7)).to.equal('1.9053011');
    expect(angles.M_to_nu(-14, 5.5).toFixed(7)).to.equal('-1.4133834');
  });

  it('should calculate M from D when near parabolic', function () {
    expect(almostEqual(angles.D_to_M_near_parabolic(0.5, 0.995), 0.542104676492085, 0, 1e-7)).to.be.true;
    expect(almostEqual(angles.D_to_M_near_parabolic(0.5, 1.005), 0.5412296752578174, 0, 1e-7)).to.be.true;
  });
});
