// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef IONIC_TTP_SLOW_NA_H
#define IONIC_TTP_SLOW_NA_H

#include "IonicModelTTP.h"

/**
 * @brief Ten Tusscher-Panfilov ionic model with an ultra-slow I_Na
 * inactivation gate.
 *
 * This model is identical to @ref TTP except for one extra gate, @f$k@f$, on
 * the fast sodium current:
 * @f[
 *   I_\text{Na} = G_\text{Na} m^3 h j k^3 (V - E_\text{Na})\;,
 * @f]
 * where @f$k@f$ follows
 * @f$\frac{\text{d}k}{\text{d}t} = \frac{k_\infty(V) - k}{\tau_k(V)}@f$, with
 * @f[ \begin{aligned}
 *   k_\infty(V) &= \frac{1}{1 + \exp\left((V - V_{1/2,k}) / s_k\right)}\;, \\
 *   \tau_k(V) &= 5 + \frac{\tau_\text{usi}}{1 + \exp\left((V + 15) / 4\right)}
 *     \;\text{[ms]}.
 * \end{aligned} @f]
 * The defaults are @f$V_{1/2,k} = 10@f$ mV, @f$s_k = 5@f$ mV and
 * @f$\tau_\text{usi} = 1050@f$ ms. With these values @f$\tau_k@f$ is about
 * 5 ms at plateau potentials and about 1055 ms at rest. The gate therefore
 * inactivates quickly during the action potential and recovers very slowly
 * during diastole, giving a rate-dependent reduction of @f$I_\text{Na}@f$.
 * Like the other TTP gates, @f$k@f$ is advanced with Rush-Larsen.
 *
 * The XML parameters are those of @ref TTP plus the optional parameters
 * @c usi_param (@f$\tau_\text{usi}@f$), @c V_half_k (@f$V_{1/2,k}@f$) and
 * @c slope_k (@f$s_k@f$). The gate's initial value is the gating variable
 * @c k_ultraslow_Na (default 1).
 *
 * Source: the openCARP model @c tenTusscherPanfilov_slow_Na (the variables
 * @c k, @c k_inf, @c tau_k and @c usi_param in its @c .model file).
 */
class TTPSlowNa : public TTP {
public:
  /// Model label.
  static inline const std::string label = "TTP_slowNa";

  /// Gating variables: the TTP gates followed by the ultra-slow I_Na gate.
  static inline const InitialStates initial_Xg = [] {
    InitialStates states = TTP::initial_Xg;
    states.emplace_back("k_ultraslow_Na", 1.0);
    return states;
  }();

  /// Index of the ultra-slow I_Na gate in the gating-variable vector.
  static constexpr unsigned int k_index = 12;

  /// Model parameters class.
  class Parameters : public TTP::Parameters {
  public:
    Parameters()
        : TTP::Parameters(TTPSlowNa::label, TTP::initial_X,
                          TTPSlowNa::initial_Xg) {
      constexpr bool required = true;

      add_parameter("usi_param", 1050.0, !required);
      add_parameter("V_half_k", 10.0, !required);
      add_parameter("slope_k", 5.0, !required);
    }
  };

  /// Constructor.
  TTPSlowNa() : TTP(TTP::initial_X, TTPSlowNa::initial_Xg) {}

  /// Construct an instance of model parameters.
  virtual std::unique_ptr<IonicModelParameters>
  get_parameters() const override {
    return std::make_unique<Parameters>();
  }

  /// Read model parameters from a parameter object.
  virtual void read_parameters(const IonicModelParameters &params) override;

  /// Distribute model parameters to all parallel processes.
  virtual void distribute_parameters(const CmMod &cm_mod,
                                     const cmType &cm) override;

protected:
  /// @name Ultra-slow I_Na gate parameters
  /// @{

  /// Amplitude of the voltage-dependent part of tau_k [ms]
  double usi_param = 1050.0;

  /// Half-inactivation voltage of k_inf [mV]
  double V_half_k = 10.0;

  /// Slope factor of k_inf [mV]
  double slope_k = 5.0;

  /// @}

  /// Update the TTP gates, then the ultra-slow I_Na gate (Rush-Larsen).
  virtual void update_g(const unsigned int zone_id, const double dt,
                        const Vector<double> &X,
                        Vector<double> &Xg) const override;

  /// Returns @f$k^3@f$.
  virtual double
  sodium_current_gate_factor(const Vector<double> &Xg) const override;
};

#endif
