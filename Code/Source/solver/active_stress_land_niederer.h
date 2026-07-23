// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_LAND_NIEDERER_H
#define ACTIVE_STRESS_LAND_NIEDERER_H

#include "active_stress_ode.h"

class LandNiederer : public ActiveStressODE {
public:
  /// Model label.
  static inline const std::string label = "LandNiederer";

  /// Model parameters class.
  class Parameters : public ActiveStressODE::Parameters {
  public:
    Parameters() : ActiveStressODE::Parameters(label) {
      constexpr bool required = true;

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

  double CaRef;

  double eta_Tm;

  double k_uw;

  double k_ws;

  double Tref;

  double k_TRPN;

  double eta_TRPN; 

  double k_u; 

  double TRPN50; 

  double rw; 

  double rs; 

  double gamma_s; 

  double gamma_w; 

  double phi; 

  double Aeff; 

  double beta0; 

  double beta1; 

  /// @}
};



#endif