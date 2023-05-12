import angles from './angles.js';
import { coe2rv, rv2coe } from './elements.js';
import { farnocchia_coe } from './propagation.js';

class Orbit {
  /**
   * Creates a new Orbit object
   *
   * @param {number} mu Gravitational parameter (km^3 / s^2)
   * @param {number} p Semi-latus rectum (km)
   * @param {number} ecc Eccentricity (rad)
   * @param {number} inc Inclination (rad)
   * @param {number} raan Longitude of ascending node (rad)
   * @param {number} argp Argument of periapsis (rad)
   * @param {number} nu True anomaly (rad)
   * @param {number} epoch Time of true anomaly (sec)
   */
  constructor(mu, p, ecc, inc, raan, argp, nu, epoch = 0) {
    this.mu = mu;
    this.p = p;
    this.ecc = ecc;
    this.inc = inc;
    this.raan = raan;
    this.argp = argp;
    this.nu = nu;
    this.epoch = epoch;
  }

  /**
   * Creates a new Orbit object from a set of state vectors
   *
   * @param {number} mu Gravitational parameter (km^3 / s^2)
   * @param {number[]} r Position vector (km)
   * @param {number[]} v Velocity vector (km / s)
   * @param {number} epoch Time of true anomaly (sec)
   * @returns
   */
  static fromStateVectors(mu, r, v, epoch = 0) {
    const { p, ecc, inc, raan, argp, nu } = rv2coe(mu, r, v);
    return new Orbit(mu, p, ecc, inc, raan, argp, nu, epoch);
  }

  // See constructor for signature
  static fromClassicElements(mu, p, ecc, inc, raan, argp, nu, epoch = 0) {
    return new Orbit(mu, p, ecc, inc, raan, argp, nu, epoch);
  }

  // Returns the semi-major axis in km
  get a() {
    return this.p / (1 - this.ecc ** 2);
  }

  // Returns the classical orbital elements
  get coe() {
    return { p: this.p, ecc: this.ecc, inc: this.inc, raan: this.raan, argp: this.argp, nu: this.nu };
  }

  // Returns the period in seconds
  get period() {
    return 2 * Math.PI * Math.sqrt(this.a ** 3 / this.mu);
  }

  // Returns the orbital radius in km
  get radius() {
    return this.p / (1 + this.ecc * Math.cos(this.nu));
  }

  // Returns the state vectors in km / sec and km
  get rv() {
    const [ r, v ] = coe2rv(this.mu, this.p, this.ecc, this.inc, this.raan, this.argp, this.nu);
    return { r, v };
  }

  /**
   * Propagates the orbit by a given time of flight. Updates epoch and true anomaly
   *
   * @param {number} tof Time of flight (sec)
   */
  propagateFor(tof) {
    this.nu = farnocchia_coe(this.mu, this.p, this.ecc, this.inc, this.raan, this.argp, this.nu, tof);
    this.epoch += tof;
  }

  /**
   * Propagates the orbit to a specific time. Updates epoch and true anomaly
   *
   * @param {number} epoch Epoch to propagate to (sec)
   */
  propagateTo(epoch) {
    const tof = epoch - this.epoch;
    this.propagateFor(tof);
  }

  /**
   * Samples the orbit at a given time and returns state vectors.
   *
   * @param {number} epoch Time of true anomaly (sec)
   * @returns
   */
  sampleAtEpoch(epoch) {
    const tof = epoch - this.epoch;
    const nu = farnocchia_coe(this.mu, this.p, this.ecc, this.inc, this.raan, this.argp, this.nu, tof);
    return this.sampleAtAngle(nu);
  }

  /**
   * Samples the orbit at a given true anomaly and returns state vectors.
   *
   * @param {number} nu True anomaly (rad)
   * @returns
   */
  sampleAtAngle(nu) {
    const [ r, v ] = coe2rv(this.mu, this.p, this.ecc, this.inc, this.raan, this.argp, nu);
    return { r, v };
  }

  /**
   * Samples the orbit at set of times and returns position vectors.
   *
   * @param {number} samples Number of samples to take
   * @param {number} tof Time of flight between samples (sec)
   * @param {number} start Epoch of first sample (sec
   * @returns {number[]{}} Array of position vectors (km / s & km)
   */
  ephem(samples = 100, tof = null, start = null) {
    if (!tof && this.ecc >= 1) throw new Error('tof must be specified for non elliptical orbits');
    if (!tof) tof = this.period;
    if (!start) start = this.epoch;

    const dt = tof / samples;
    const times = Array.from({ length: samples }, (_, i) => start + i * dt);
    return times.map(t => this.sampleAtEpoch(t));
  }
}

export default Orbit;
