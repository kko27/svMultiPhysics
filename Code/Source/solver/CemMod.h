// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// Cardiac electromechanics (CEM) coupling model.
//
// Consolidates the former cemModelType (from CepMod.h) and all electromechanics
// parameters and time-integration routines that were previously scattered across
// CepModAp, CepModBo, and CepModTtp.

#ifndef CEM_MOD_H
#define CEM_MOD_H

#include "Array.h"
#include "Vector.h"

/// @brief Cardiac electromechanics coupling model.
///
/// Holds the coupling flags, the per-node activation variable Ya, all
/// electromechanics parameters, and explicit time-stepping routines for the
/// active stress and active strain formulations.
///
/// Active stress ODE (Nash & Panfilov 2004; Keldermann et al. 2010):
///   dT_act/dt = eps(xi) * [ eta_T*(xi - xi_rest) - T_act ]
/// where xi is the driving variable (voltage for AP/BO; Ca_i for TTP) and
///   eps(xi) = eps_0 + (eps_i - eps_0)*exp(-exp(-xi_T*(xi - xi_crit)))
///
/// Active strain ODE (Nobile, Quarteroni & Ruiz-Baier 2012):
///   mu_C * c^2 * d(gamma_f)/dt = F_active(c, I4f) + 2*I4f*(1/(1+gamma_f)^3 - 1)
class CemMod
{
  public:
    CemMod();

    // -----------------------------------------------------------------------
    // Coupling flags (formerly cemModelType)
    // -----------------------------------------------------------------------

    /// Whether electrophysiology and mechanics are coupled
    bool cpld = false;

    /// Whether active stress formulation is employed
    bool aStress = false;

    /// Whether active strain formulation is employed
    bool aStrain = false;

    /// Per-node activation variable (size = tnNo).
    ///   Active stress:  T_act [MPa] — fiber activation force
    ///   Active strain:  gamma_f    — fiber stretch increment
    Vector<double> Ya;

    // -----------------------------------------------------------------------
    // Active stress parameters
    // -----------------------------------------------------------------------

    /// Resting value of driving variable (V_rest [mV] for AP/BO; Ca_rest [mM] for TTP)
    double xi_rest = 0.0;

    /// Critical value of driving variable for the Gompertz transition
    double xi_crit = 0.0;

    /// Activation-stress scale factor [MPa/mM or ms^{-1}mV^{-1}]
    double eta_T = 0.0;

    /// Minimum transition rate [ms^{-1}]
    double eps_0 = 0.1;

    /// Maximum transition rate [ms^{-1}]
    double eps_i = 1.0;

    /// Gompertz sharpness [mM^{-1} or mV^{-1}]
    double xi_T = 0.0;

    // -----------------------------------------------------------------------
    // Active strain parameters
    // -----------------------------------------------------------------------

    /// Active sarcomere force constant [mM^{-2}]
    double alFa = -4.0e6;

    /// Resting Ca^{2+} concentration [mM]
    double c0 = 2.155e-4;

    /// Viscous-type constant [ms mM^{-2}]
    double mu_C = 5.0e6;

    /// Initial sarcomere length [um]
    double SL0 = 1.95;

    /// Minimum sarcomere length [um]
    double SLmin = 1.70;

    /// Maximum sarcomere length [um]
    double SLmax = 2.60;

    // Fourier coefficients for the sarcomere force-length relationship
    double f0  = -4333.618335582119;
    double fc1 =  2570.395355352195;
    double fs1 = -2051.827278991976;
    double fc2 =  1329.53611689133;
    double fs2 =   302.216784558222;
    double fc3 =   104.943770305116;
    double fs3 =   218.375174229422;

    // -----------------------------------------------------------------------
    // Active stress time integration
    // -----------------------------------------------------------------------

    /// Semi-implicit (Rush-Larsen) step — exact for the linear activation ODE
    /// at fixed xi; this matches the original actv_strs behaviour in all ionic models.
    void actv_strs_si(double xi, double dt, double& Tact) const;

    /// Forward-Euler step for the activation ODE.
    void actv_strs_fe(double xi, double dt, double& Tact) const;

    /// Fourth-order Runge-Kutta step for the activation ODE.
    void actv_strs_rk(double xi, double dt, double& Tact) const;

    // -----------------------------------------------------------------------
    // Active strain time integration
    // -----------------------------------------------------------------------

    /// Forward-Euler step for the fiber stretch ODE — matches original behaviour.
    void actv_strn_fe(double c_Ca, double I4f, double dt, double& gf) const;

    /// Fourth-order Runge-Kutta step for the fiber stretch ODE.
    void actv_strn_rk(double c_Ca, double I4f, double dt, double& gf) const;

    // -----------------------------------------------------------------------
    // Active deformation gradient
    // -----------------------------------------------------------------------

    /// Compute the active deformation gradient tensor F_a from gamma_f.
    /// Replaces the former free function actv_strain() in mat_models.
    ///
    /// F_a = I + gamma_f*(f x f) + gamma_s*(s x s) + gamma_n*(n x n)
    /// with  gamma_n = 4*gamma_f  and  det(F_a) = 1 => gamma_s from constraint.
    void compute_Fa(int nsd, double gf, const Array<double>& fl, Array<double>& Fa) const;

  private:

    /// Gompertz-type transition rate epsilon(xi).
    double eps_func(double xi) const;

    /// Sarcomere force-length function F_SL(I4f).
    double force_length(double I4f) const;

    /// Right-hand side of the fiber stretch ODE df(gf)/dt.
    double actv_strn_rhs(double c_Ca, double I4f, double gf) const;
};

#endif
