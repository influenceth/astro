import * as math from 'mathjs';

/**
 * Implements modulo (as distinct from the remainder operator)
 *
 * @param {number} x
 * @param {number} y
 * @returns {number}
 */
export const modulo = (x, y) => ((x % y) + y) % y;

/**
 * Generates a rotation matrix
 *
 * @param {number} angle
 * @param {number} axis
 */
export const rotation_matrix = (angle, axis) => {
  const c = Math.cos(angle);
  const s = Math.sin(angle);
  const a1 = modulo((axis + 1), 3);
  const a2 = modulo((axis + 2), 3);

  const R = math.zeros([3, 3]);
  R[axis][axis] = 1;
  R[a1][a1] = c;
  R[a1][a2] = -s;
  R[a2][a1] = s;
  R[a2][a2] = c;

  return R;
};

export default {
  modulo,
  rotation_matrix
};
