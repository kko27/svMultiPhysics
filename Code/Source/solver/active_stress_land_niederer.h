// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_LAND_NIEDERER_H
#define ACTIVE_STRESS_LAND_NIEDERER_H

#include "ActiveStressODE.h"

class LandNiederer : public ActiveStressODE {
public:
  /// Model label, used for factory registration and XML selection. 
  static inline const std::string label = "LandNiederer";

  /*
   * @brief Model parameters class.
   * Declares the parameters required by the model. All parameters are 
   * marked as required, and omitting a parameter will cause a parse error.
  */ 
  class Parameters : public ActiveStressODE::Parameters {
  public:
    Parameters() : ActiveStressODE::Parameters(label) {
      constexpr bool required = true;

      // Reference values: Land et al. 2017 
      add_parameter("CaRef", 0.805, required);
      add_parameter("eta_Tm", 5.0, required);
      add_parameter("k_uw", 0.182, required);
      add_parameter("k_ws", 0.012, required);
      add_parameter("Tref", 120.0, required);
      add_parameter("k_TRPN", 0.1, required);
      add_parameter("eta_TRPN", 2.0, required);
      add_parameter("k_u", 1.0, required); 
      add_parameter("TRPN50", 0.35, required);
      add_parameter("rw", 0.5, required);
      add_parameter("rs", 0.25, required);
      add_parameter("gamma_s", 0.0085, required);
      add_parameter("gamma_w", 0.615, required);
      add_parameter("phi", 2.23, required);
      add_parameter("Aeff", 25.0, required);
      add_parameter("beta0", 2.3, required);
      add_parameter("beta1", -2.4, required);      
    }
  };

  /**
   * @brief Constructor.
   */
  LandNiederer() : ActiveStressODE(7) {}

  /**
   * @brief Construct an instance of model parameters.
   */
  virtual std::unique_ptr<ActiveStressModelParameters>
  get_parameters() const override {
    return std::make_unique<Parameters>();
  }

protected:
  /**
   * @brief Read model parameters from a parameter object.
   */
  virtual void read_model_specific_parameters(
      const ActiveStressModelParameters &params) override;

  /**
   * @brief Distribute model parameters to all parallel processes.
   */
  virtual void distribute_model_specific_parameters(const CmMod &cm_mod,
                                                    const cmType &cm) override;

  /**
   * @brief Initialize the state vector for a single node.
   *
   * @param[out] state State vector for a single node, to be initialized by
   *   this function.
   */
  virtual void init_local(Vector<double> &state) const override;

  virtual void advance_time_step_local(const double t, const double dt,
                                       const double calcium,
                                       const double fiber_stretch,
                                       const double fiber_stretch_rate,
                                       Vector<double> &state) const override;

  /**
   * @brief Compute the rate of change in the state variables.
   */
  virtual Vector<double> getf(const double t, const Vector<double> &state,
                              const double calcium, const double fiber_stretch,
                              const double fiber_stretch_rate) const override;

  /**
   * @brief Compute the active tension for a single node.
   */
  virtual double
  compute_active_tension_local(const Vector<double> &state) const override;

  /**
   * @brief Compute the active stiffness for a single node.
   */
  virtual double 
  compute_active_stiffness_local(const Vector<double> &state, 
                                 const double fiber_stretch, 
                                 const double fiber_stretch_rate) const override;

  /// @name Model parameters.
  /// @{
  /// Reference intracellular calcium concentration giving half-maximal
  /// troponin C saturation, @f$[Ca^{2+}]_{T50,ref}@f$, used (with length
  /// dependence via beta1) in the CaTRPN binding ODE [uM]
  double CaRef;

  /// Cooperativity (Hill) exponent @f$n_{Tm}@f$ for the tropomyosin
  /// blocked/unblocked (B/U) transition; sets the steepness of thin-filament
  /// activation by CaTRPN [-]
  double eta_Tm;

  /// Rate constant for transition from unbound (U) to pre-power stroke (W) state [1/ms]
  double k_uw;

  /// Rate constant for transition from pre-power stroke (W) to post-power stroke (S) state [1/ms]
  double k_ws;

  /// Maximum observable tension @f$\Tref@f$ at resting length [kPa]
  double Tref;

  /// Rate constant @f$k_{TRPN}@f$ governing calcium binding/unbinding
  /// kinetics of troponin C (the "pace" of CaTRPN relaxation to its
  /// steady-state value) [1/ms]
  double k_TRPN;

  /// Cooperativity (Hill) exponent @f$n_{TRPN}@f$ of the calcium-troponin C
  /// binding rate [-]
  double eta_TRPN;

  /// Tropomyosin unblocking rate constant @f$k_u@f$ (rate at which blocked
  /// binding sites become unblocked); together with rw, rs, and TRPN50 it
  /// fixes the blocking rate kb [1/ms]
  double k_u;

  /// The value of CaTRPN where bounded state B = 0.5 in steady-state [-]
  double TRPN50;

  /// Steady-state ratio @f$r_w@f$ between the pre-powerstroke (W) population
  /// and the non-strongly-bound (U+W) population at equilibrium; sets the
  /// reverse rate constant kwu [-]
  double rw;

  /// Steady-state duty ratio @f$r_s@f$: the fraction of cross-bridges in the
  /// strongly-bound (post-powerstroke, S) state at equilibrium; sets the
  /// reverse rate constant ksu [-]
  double rs;

  /// Strain-dependent detachment-rate coefficient @f$\gamma_s@f$ for the
  /// strongly-bound (post-powerstroke) S state; scales how fast cross-bridges
  /// are pulled off by distortion in that state [-]
  double gamma_s;

  /// Strain-dependent detachment-rate coefficient @f$\gamma_w@f$ for the
  /// weakly-bound (pre-powerstroke) W state; scales how fast cross-bridges
  /// are pulled off by distortion in that state [-]
  double gamma_w;

  /// Proportionality factor @f$\phi@f$ relating the cross-bridge distortion
  /// decay rates (cw, cs) to the forward cycling rates k_uw, k_ws [-]
  double phi;

  /// Rescaled magnitude @f$A_{eff}@f$ of the immediate cross-bridge
  /// distortion response to fiber stretch rate; sets the distortion
  /// sensitivities As, Aw [-]
  double Aeff;

  /// Coefficient @f$\beta_0@f$ controlling the length-dependence of maximal
  /// tension through changes in filament overlap (Frank-Starling effect) [-]
  double beta0;

  /// Coefficient @f$\beta_1@f$ controlling the length-dependence of calcium
  /// sensitivity (shifts CaRef/[Ca2+]T50 with sarcomere stretch) [-]
  double beta1;
};



#endif