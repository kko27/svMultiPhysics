// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_TANH_FIBER_STRETCH_H
#define ACTIVE_STRESS_TANH_FIBER_STRETCH_H

#include "ActiveStress.h"

#include "FourierInterpolation.h"

/**
 * @brief Active stress model driven by a time-dependent activation signal
 * and a tanh^2 fiber-stretch switch.
 *
 * Defines an active tension
 * @f[
 *   \Tact(t, \calcium, \fiberstretch, \fiberstretchrate, \astressstate) =
 *     \eta(t) \tanh^2\left(a_1 (\fiberstretch - a_2)\right)\;,
 * @f]
 * where @f$\eta(t)@f$ is a user-defined function of time (constant ramp or an
 * arbitrary time series, as in @ref ActiveStressUniformUnsteady), and
 * @f$a_1@f$, @f$a_2@f$ are scalar model parameters. @f$\fiberstretch@f$ is
 * the fiber stretch, i.e. @f$\sqrt{I_{4f}}@f$.
 */
class ActiveStressTanhFiberStretch : public ActiveStress {
public:
  /// Model label.
  static inline const std::string label = "TanhFiberStretch";

  /// Model parameters class.
  class Parameters : public ActiveStressModelParameters {
  public:
    Parameters() : ActiveStressModelParameters(label) {
      constexpr bool required = true;

      add_parameter("a1", 0.0, required);
      add_parameter("a2", 0.0, required);
      add_parameter("Ramp", false, required);
      add_parameter("Temporal_values_file_path", std::string(""), required);
    }
  };

  /**
   * @brief Constructor.
   */
  ActiveStressTanhFiberStretch()
      : ActiveStress(/* n_states = */ 0,
                     /* needs_fiber_stretch = */ true,
                     /* needs_fiber_stretch_rate = */ false) {}

  /**
   * @brief Construct an instance of model parameters.
   */
  virtual std::unique_ptr<ActiveStressModelParameters>
  get_parameters() const override {
    return std::make_unique<Parameters>();
  }

  /**
   * @brief Initialization.
   *
   * Calls the parent class initialization method, and reads eta(t) from file.
   */
  virtual void init(const unsigned int tnNo) override;

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
   * This model has no states, so this function does nothing.
   */
  virtual void init_local(Vector<double> &state) const override {}

  /**
   * @brief Advance in time for a single node.
   *
   * This model has no states, so this function does nothing.
   */
  virtual void advance_time_step_local(const double t, const double dt,
                                       const double calcium,
                                       const double fiber_stretch,
                                       const double fiber_stretch_rate,
                                       Vector<double> &state) const override {}

  /**
   * @brief Compute the active tension for a single node.
   */
  virtual double
  compute_active_tension_local(const Vector<double> &state,
                               const double fiber_stretch) const override;

  /// Coefficient controlling the steepness of the tanh^2 switch.
  double a1;

  /// Fiber stretch threshold around which the tanh^2 switch is centered.
  double a2;

  /// Toggle between ramp or Fourier transform for eta(t).
  bool ramp;

  /// Name of the file containing the temporal values for eta(t).
  std::string temporal_values_file_path;

  /// Fourier interpolation of the time dependent activation signal eta(t).
  FourierInterpolation fourier_interpolation;
};

#endif
