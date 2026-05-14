// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "CemMod.h"
#include "mat_fun.h"
#include "utils.h"

#include <cmath>

CemMod::CemMod()
{
}

// -----------------------------------------------------------------------
// Private helpers
// -----------------------------------------------------------------------

/// Gompertz-type transition rate: eps(xi) = eps_0 + (eps_i-eps_0)*exp(-exp(-xi_T*(xi-xi_crit)))
double CemMod::eps_func(const double xi) const
{
  return eps_0 + (eps_i - eps_0) * exp(-exp(-xi_T * (xi - xi_crit)));
}

/// Sarcomere force-length relationship F_SL(I4f) as a Fourier series.
/// Returns 0 outside [SLmin, SLmax].
double CemMod::force_length(const double I4f) const
{
  const double SL = I4f * SL0;
  if (SL < SLmin || SL > SLmax) {
    return 0.0;
  }
  return 0.5*f0
       + fc1*cos(SL)       + fs1*sin(SL)
       + fc2*cos(2.0*SL)   + fs2*sin(2.0*SL)
       + fc3*cos(3.0*SL)   + fs3*sin(3.0*SL);
}

/// Right-hand side of the fiber stretch ODE:
///   d(gf)/dt = [ F_active(c, I4f) + 2*I4f*(1/(1+gf)^3 - 1) ] / (mu_C * c^2)
double CemMod::actv_strn_rhs(const double c_Ca, const double I4f, const double gf) const
{
  const double FL   = force_length(I4f);
  const double Fact = alFa * (c_Ca - c0) * (c_Ca - c0) * FL;
  const double rtmp = 2.0 * I4f * (1.0 / pow(1.0 + gf, 3.0) - 1.0);
  return (Fact + rtmp) / (mu_C * c_Ca * c_Ca);
}

// -----------------------------------------------------------------------
// Active stress — time integration
// -----------------------------------------------------------------------

/// Semi-implicit (Rush-Larsen) step for dT/dt = eps*(eta_T*(xi-xi_rest) - T).
/// This is equivalent to exact integration for the linear ODE at fixed xi
/// and matches the original actv_strs implementation in CepModAp/Bo/Ttp.
void CemMod::actv_strs_si(const double xi, const double dt, double& Tact) const
{
  const double e = eps_func(xi);
  Tact = (Tact + e * dt * eta_T * (xi - xi_rest)) / (1.0 + e * dt);
}

/// Forward-Euler step for dT/dt = eps*(eta_T*(xi-xi_rest) - T).
void CemMod::actv_strs_fe(const double xi, const double dt, double& Tact) const
{
  const double e = eps_func(xi);
  Tact = Tact + dt * e * (eta_T * (xi - xi_rest) - Tact);
}

/// Fourth-order Runge-Kutta step for dT/dt = eps*(eta_T*(xi-xi_rest) - T).
/// eps(xi) is treated as constant over the substep (xi is the EP driver value).
void CemMod::actv_strs_rk(const double xi, const double dt, double& Tact) const
{
  const double e = eps_func(xi);
  const auto f = [&](const double T) { return e * (eta_T * (xi - xi_rest) - T); };

  const double k1 = f(Tact);
  const double k2 = f(Tact + 0.5*dt*k1);
  const double k3 = f(Tact + 0.5*dt*k2);
  const double k4 = f(Tact +     dt*k3);

  Tact = Tact + (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
}

// -----------------------------------------------------------------------
// Active strain — time integration
// -----------------------------------------------------------------------

/// Forward-Euler step for the fiber stretch ODE — matches the original
/// actv_strn implementation in CepModBo/Ttp.
void CemMod::actv_strn_fe(const double c_Ca, const double I4f,
                            const double dt, double& gf) const
{
  gf = gf + dt * actv_strn_rhs(c_Ca, I4f, gf);
}

/// Fourth-order Runge-Kutta step for the fiber stretch ODE.
/// c_Ca and I4f are treated as constant over the substep.
void CemMod::actv_strn_rk(const double c_Ca, const double I4f,
                            const double dt, double& gf) const
{
  const double k1 = actv_strn_rhs(c_Ca, I4f, gf);
  const double k2 = actv_strn_rhs(c_Ca, I4f, gf + 0.5*dt*k1);
  const double k3 = actv_strn_rhs(c_Ca, I4f, gf + 0.5*dt*k2);
  const double k4 = actv_strn_rhs(c_Ca, I4f, gf +     dt*k3);

  gf = gf + (dt / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
}

// -----------------------------------------------------------------------
// Active deformation gradient
// -----------------------------------------------------------------------

/// Compute the active deformation gradient F_a given the fiber stretch
/// increment gamma_f and the fiber/sheet direction array fl.
///
/// F_a = I + gamma_f*(f x f) + gamma_s*(s x s) + gamma_n*(n x n)
///
/// Transverse constraints (Nobile, Quarteroni & Ruiz-Baier 2012):
///   gamma_n = 4 * gamma_f
///   gamma_s = 1 / ((1 + gamma_f)*(1 + gamma_n)) - 1   (det F_a = 1)
///
/// Replaces the former free function actv_strain() in mat_models.
void CemMod::compute_Fa(const int nsd, const double gf,
                         const Array<double>& fl, Array<double>& Fa) const
{
  using namespace mat_fun;

  const auto af = fl.col(0);
  const auto as = fl.col(1);
  const auto an = utils::cross(fl);

  const double gn = 4.0 * gf;
  const double gs = 1.0 / ((1.0 + gf) * (1.0 + gn)) - 1.0;

  const auto IDm = mat_id(nsd);
  const auto Hf  = mat_dyad_prod(af, af, nsd);
  const auto Hs  = mat_dyad_prod(as, as, nsd);
  const auto Hn  = mat_dyad_prod(an, an, nsd);

  Fa = IDm + gf*Hf + gs*Hs + gn*Hn;
}
