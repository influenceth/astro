import { modulo } from './utils.js';

// See: https://github.com/poliastro/poliastro/blob/main/src/poliastro/core/angles.py

// Converts eccentric anomaly to mean anomaly
export const E_to_M = (E, ecc) => {
  return E - ecc * Math.sin(E);
};

// Converts mean anomaly to eccentric anomaly
export const M_to_E = (M, ecc) => {
  let E1, fVal, fDer, step;
  let E = M < 0 ? M - ecc : M + ecc; // Initial guess for elliptical eccentric anomaly

  for (let i = 0; i < 50; i++) {
    fVal = E_to_M(E, ecc) - M;
    fDer = 1 - ecc * Math.cos(E);
    step = fVal / fDer;
    E1 = E - step;

    if (Math.abs(E1 - E) < 1e-7) {
      break;
    } else {
      E = E1;
    }
  }

  return E1;
}

// Converts eccentric anomaly to true anomaly
export const E_to_nu = (E, ecc) => {
  return 2 * Math.atan(Math.sqrt((1 + ecc) / (1 - ecc)) * Math.tan(E / 2));
}

// Converts true anomaly to eccentric anomaly
export const nu_to_E = (nu, ecc) => {
  return 2 * Math.atan(Math.sqrt((1 - ecc) / (1 + ecc)) * Math.tan(nu / 2));
}

// Converts true anomaly to hyperbolic eccentric anomaly
export const F_to_M = (F, ecc) => {
  return ecc * Math.sinh(F) - F;
};

// Converts mean anomaly to hyperbolic eccentric anomaly
export const M_to_F = (M, ecc) => {
  let F1, fVal, fDer, step;
  let F = Math.asinh(M / ecc); // Initial guess for hyperbolic eccentric anomaly

  for (let i = 0; i < 50; i++) {
    fVal = F_to_M(F, ecc) - M;
    fDer = ecc * Math.cosh(F) - 1;
    step = fVal / fDer;
    F1 = F - step;

    if (Math.abs(F1 - F) < 1e-7) {
      break;
    } else {
      F = F1;
    }
  }

  return F;
}

// Converts hyperbolic eccentric anomaly to true anomaly
export const F_to_nu = (F, ecc) => {
  return 2 * Math.atan(Math.sqrt((ecc + 1) / (ecc - 1)) * Math.tanh(F / 2));
}

// Converts true anomaly to hyperbolic eccentric anomaly
export const nu_to_F = (nu, ecc) => {
  return Math.acosh((ecc + Math.cos(nu)) / (1 + ecc * Math.cos(nu)));
}

// Converts mean anomaly to true anomaly
export const M_to_nu = (M, ecc, delta = 1e-2) => {
  if (ecc < 1 - delta) {
    // Eccentric
    M = modulo((M + Math.PI), (2 * Math.PI)) - Math.PI;
    return E_to_nu(M_to_E(M, ecc), ecc);
  } else if (ecc < 1) {
    // Near parabolic low
    return D_to_nu(M_to_D_near_parabolic(M, ecc));
  } else if (ecc == 1) {
    // Parabolic
    return D_to_nu(M_to_D(M));
  } else if (ecc < 1 + delta) {
    // Near parabolic high
    return D_to_nu(M_to_D_near_parabolic(M, ecc));
  } else {
    // Hyperbolic
    return F_to_nu(M_to_F(M, ecc), ecc);
  }
}

// Converts true anomaly to mean anomaly
export const nu_to_M = (nu, ecc, delta = 1e-2) => {
  if (ecc < 1 - delta) {
    // Elliptic
    return E_to_M(nu_to_E(nu, ecc), ecc);
  } else if (ecc < 1) {
    // Near parabolic low
    return D_to_M_near_parabolic(nu_to_D(nu), ecc);
  } else if (ecc == 1) {
    // Parabolic
    return D_to_M(nu_to_D(nu));
  } else if (ecc < 1 + delta) {
    // Near parabolic high
    return D_to_M_near_parabolic(nu_to_D(nu), ecc);
  } else {
    // Hyperbolic
    return F_to_M(nu_to_F(nu, ecc), ecc);
  }
}

// Converts parabolic anomaly to mean anomaly
export const D_to_M = (D) => {
  return D + D ** 3 / 3;
}

export const D_to_M_near_parabolic = (D, ecc) => {
  const x = (ecc - 1) / (ecc + 1) * (D ** 2);
  if (Math.abs(x) >= 1) throw new Error('abs(x) must be less than 1');
  const S = _S_x(ecc, x);
  return Math.sqrt(2 / (1 + ecc)) * D + Math.sqrt(2 / (1 + ecc) ** 3) * (D**3) * S;
}

// Parabolic eccentric anomaly from mean anomaly, near parabolic case.
export const M_to_D_near_parabolic = (M, ecc) => {
  let D1, fVal, fDer, step;
  let D = M_to_D(M); // initial guess for parabolic eccentric anomaly

  for (let i = 0; i < 50; i++) {
    fVal = D_to_M_near_parabolic(D, ecc) - M;
    fDer = _kepler_equation_prime_near_parabolic(D, ecc);
    step = fVal / fDer;
    D1 = D - step;

    if (Math.abs(D1 - D) < 1.48e-8) {
      break;
    } else {
      D = D1;
    }
  }

  return D;
}

// Converts mean anomaly to parabolic anomaly
export const M_to_D = (M) => {
  const B = 3 * M / 2;
  const A = (B + (1 + B ** 2) ** 0.5) ** (2 / 3);
  return 2 * A * B / (1 + A + A ** 2);
}

// Converts parabolic anomaly to true anomaly
export const D_to_nu = (D) => {
  return 2 * Math.atan(D)
}

// Converts true anomaly to parabolic anomaly
export const nu_to_D = (nu) => {
  return Math.tan(nu / 2);
}

const _kepler_equation_prime_near_parabolic = (D, ecc) => {
  const x = (ecc - 1) / (ecc + 1) * (D ** 2);
  if (Math.abs(x) >= 1) throw new Error('abs(x) must be less than 1');
  S = _dS_x_alt(ecc, x);
  return Math.sqrt(2 / (1 + ecc)) + Math.sqrt(2 / (1 + ecc) ** 3) * (D ** 2) * S;
};

const _S_x = (ecc, x, atol = 1e-12) => {
  if (Math.abs(x) >= 1) throw new Error('abs(x) must be less than 1');
  let S = 0;
  let k = 0;

  while (true) {
    const S_old = S;
    S += (ecc - 1 / (2 * k + 3)) * x ** k;
    k += 1;

    if (Math.abs(S - S_old) < atol) return S;
  }
};

const _dS_x_alt = (ecc, x, atol = 1e-12) => {
  if (Math.abs(x) >= 1) throw new Error('abs(x) must be less than 1');
  let S = 0;
  let k = 0;

  while (true) {
    S_old = S
    S += (ecc - 1 / (2 * k + 3)) * (2 * k + 3) * x**k
    k += 1

    if (Math.abs(S - S_old) < atol) return S;
  }
};

export default {
  E_to_M,
  M_to_E,
  E_to_nu,
  nu_to_E,
  F_to_M,
  M_to_F,
  F_to_nu,
  nu_to_F,
  M_to_nu,
  nu_to_M,
  D_to_M,
  M_to_D,
  D_to_nu,
  nu_to_D,
  M_to_D_near_parabolic,
  D_to_M_near_parabolic
};
