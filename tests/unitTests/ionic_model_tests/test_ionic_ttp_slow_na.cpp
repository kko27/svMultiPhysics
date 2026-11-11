// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

/// @file

#include "ionic_model_test_helpers.h"
#include "IonicModelTTP.h"
#include "IonicModelTTPSlowNa.h"
#include "gtest/gtest.h"

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace {

/// Stimulus of the canonical TTP reference protocol.
double canonical_stimulus(double time)
{
  return 10.0 <= time && time < 11.0 ? -52.0 : 0.0;
}

void run_ttp_slow_na_trajectory(
    const TTPSlowNa::Parameters &parameters,
    int zone_id,
    const std::string &reference_csv_filename,
    const std::vector<double> &initial_X_override = {},
    const std::vector<double> &initial_Xg_override = {})
{
  IonicModelTrajectoryConfiguration configuration;
  configuration.integration_type = TimeIntegrationType::FE;
  configuration.zone_id = zone_id;
  configuration.time_step = 0.005;
  configuration.update_count = 120000;
  if (!initial_X_override.empty())
    configuration.initial_X_override = initial_X_override;
  if (!initial_Xg_override.empty())
    configuration.initial_Xg_override = initial_Xg_override;
  configuration.reference_csv_filename = reference_csv_filename;
  configuration.state_reference_columns = {
      "V_mV", "Ki",  "Nai", "Cai", "Ca_ss", "Ca_SR", "R_prime",
      "Xr1",  "Xr2", "Xs",  "m",   "h",     "j",     "d",
      "f",    "f2",  "fCass", "s",  "r",    "k"};
  configuration.tolerance = 1.0e-10;
  configuration.stimulus_at_time = canonical_stimulus;
  configuration.sac_coefficient = 0.0;

  IonicModelTrajectoryTest<TTPSlowNa> trajectory(parameters, configuration);
  trajectory.run();
}

/// Construct a model and read its parameters.
template <class Model>
std::unique_ptr<IonicModel> make_model(const typename Model::Parameters &params)
{
  std::unique_ptr<IonicModel> model = std::make_unique<Model>();
  model->read_parameters(params);
  return model;
}

/// Summary of a paced simulation, measured on the last beat.
struct PacingSummary {
  /// Ultra-slow gate value just before the last stimulus (NaN without gate).
  double k_before_last_stimulus = std::numeric_limits<double>::quiet_NaN();

  /// Maximum upstroke velocity dV/dt of the last beat [mV/ms].
  double max_dvdt_last_beat = 0.0;
};

/// Pace a model with 1 ms, -52 pA/pF stimuli every @p cycle_length ms.
PacingSummary pace(const IonicModel &model, double cycle_length, int beats)
{
  constexpr double dt = 0.01;
  constexpr double stimulus_duration = 1.0;
  const auto steps_per_beat =
      static_cast<std::size_t>(std::llround(cycle_length / dt));

  odeType solver;
  solver.tIntType = TimeIntegrationType::FE;

  Vector<double> X(model.nX());
  Vector<double> Xg(model.nG());
  model.init(X, Xg);

  PacingSummary summary;
  for (int beat = 0; beat < beats; ++beat) {
    const bool last_beat = beat == beats - 1;
    if (last_beat && model.nG() > TTPSlowNa::k_index)
      summary.k_before_last_stimulus = Xg(TTPSlowNa::k_index);

    for (std::size_t step = 0; step < steps_per_beat; ++step) {
      const double time_in_beat = static_cast<double>(step) * dt;
      const double time = beat * cycle_length + time_in_beat;
      const double stimulus =
          time_in_beat < stimulus_duration - 1.0e-12 ? -52.0 : 0.0;
      const double V_old = X(0);
      model.integ(solver, 1, time, dt, stimulus, 0.0, X, Xg);
      if (last_beat)
        summary.max_dvdt_last_beat =
            std::max(summary.max_dvdt_last_beat, (X(0) - V_old) / dt);
    }
  }
  return summary;
}

} // namespace

/**
 * @test Advance the EPI TTP_slowNa model (zone 1) from the CellML EPI initial
 * state, with @f$k=1@f$, for 600 ms using @f$\Delta t=0.005\,\mathrm{ms}@f$.
 *
 * The protocol is the same as @c IonicModelTrajectory.TTPEpi. The reference
 * driver adds the ultra-slow gate to the CellML TP06 core; see
 * `generate_ttp_slow_na.py` in the [generator README](../reference_generators/ionic_model/ten_tusscher_panfilov/README.md).
 * See @ref TTPSlowNa for the model formulation.
 */
TEST(IonicModelTrajectory, TTPSlowNaEpi)
{
  TTPSlowNa::Parameters parameters;
  run_ttp_slow_na_trajectory(parameters, 1,
                             "ionic_ttp_slow_na_epi_trajectory.csv");
}

/**
 * @test Advance the ENDO TTP_slowNa model (zone 2) from the CellML ENDO
 * initial state, with @f$G_\mathrm{to}=0.073@f$ and @f$k=1@f$, using the
 * canonical TTP protocol.
 */
TEST(IonicModelTrajectory, TTPSlowNaEndo)
{
  TTPSlowNa::Parameters parameters;
  parameters.set_scalar("G_to", 0.073);
  run_ttp_slow_na_trajectory(
      parameters, 2, "ionic_ttp_slow_na_endo_trajectory.csv",
      {-86.709, 138.4, 10.355, 1.3e-4, 3.6e-4, 3.715, 0.9068},
      {0.00448, 0.476, 0.0087, 0.00155, 0.7573, 0.7225,
       3.164e-5, 0.8009, 0.9778, 0.9953, 0.3212, 2.235e-8, 1.0});
}

/**
 * @test Advance the M-cell TTP_slowNa model (zone 3) from the CellML M-cell
 * initial state, with @f$G_\mathrm{Ks}=0.098@f$ and @f$k=1@f$, using the
 * canonical TTP protocol.
 */
TEST(IonicModelTrajectory, TTPSlowNaM)
{
  TTPSlowNa::Parameters parameters;
  parameters.set_scalar("G_Ks", 0.098);
  run_ttp_slow_na_trajectory(
      parameters, 3, "ionic_ttp_slow_na_m_trajectory.csv",
      {-85.423, 138.52, 10.132, 1.53e-4, 4.2e-4, 4.272, 0.8978},
      {0.0165, 0.473, 0.0174, 0.00165, 0.749, 0.6788,
       3.288e-5, 0.7026, 0.9526, 0.9942, 0.999998, 2.347e-8, 1.0});
}

/**
 * @test With @f$V_{1/2,k}@f$ so large that @f$k_\infty \equiv 1@f$ in
 * floating point, the k gate stays at exactly 1. TTP_slowNa must then
 * reproduce TTP bit for bit, for both FE and RK4.
 */
TEST(IonicModelTTPSlowNa, ReducesToTTPWhenGateIsPinned)
{
  TTP::Parameters ttp_parameters;
  TTPSlowNa::Parameters slow_parameters;
  slow_parameters.set_scalar("V_half_k", 1.0e3);

  const auto ttp = make_model<TTP>(ttp_parameters);
  const auto slow = make_model<TTPSlowNa>(slow_parameters);
  ASSERT_EQ(slow->nX(), ttp->nX());
  ASSERT_EQ(slow->nG(), ttp->nG() + 1);

  for (const auto method : {TimeIntegrationType::FE, TimeIntegrationType::RK4}) {
    odeType solver;
    solver.tIntType = method;

    Vector<double> X_ttp(ttp->nX()), Xg_ttp(ttp->nG());
    Vector<double> X_slow(slow->nX()), Xg_slow(slow->nG());
    ttp->init(X_ttp, Xg_ttp);
    slow->init(X_slow, Xg_slow);

    constexpr double dt = 0.01;
    for (int step = 0; step < 60000; ++step) {
      const double time = step * dt;
      const double stimulus = canonical_stimulus(time);
      ttp->integ(solver, 1, time, dt, stimulus, 0.0, X_ttp, Xg_ttp);
      slow->integ(solver, 1, time, dt, stimulus, 0.0, X_slow, Xg_slow);
    }

    for (int i = 0; i < X_ttp.size(); ++i)
      EXPECT_EQ(X_slow(i), X_ttp(i)) << "X[" << i << "], method " << method;
    for (int i = 0; i < Xg_ttp.size(); ++i)
      EXPECT_EQ(Xg_slow(i), Xg_ttp(i)) << "Xg[" << i << "], method " << method;
    EXPECT_EQ(Xg_slow(TTPSlowNa::k_index), 1.0) << "method " << method;
  }
}

/**
 * @test The ultra-slow gate makes I_Na availability rate dependent. With
 * faster pacing, k has less time to recover before the next beat, and the
 * upstroke velocity is lower than in TTP under the same pacing.
 */
TEST(IonicModelTTPSlowNa, RateDependentSodiumAvailability)
{
  TTP::Parameters ttp_parameters;
  TTPSlowNa::Parameters slow_parameters;
  const auto ttp = make_model<TTP>(ttp_parameters);
  const auto slow = make_model<TTPSlowNa>(slow_parameters);

  const PacingSummary fast = pace(*slow, 400.0, 2);
  const PacingSummary slow_rate = pace(*slow, 1000.0, 2);
  const PacingSummary ttp_fast = pace(*ttp, 400.0, 2);

  EXPECT_LT(fast.k_before_last_stimulus, 1.0);
  EXPECT_LT(fast.k_before_last_stimulus, slow_rate.k_before_last_stimulus);
  EXPECT_LT(fast.max_dvdt_last_beat, ttp_fast.max_dvdt_last_beat);
  EXPECT_LT(fast.max_dvdt_last_beat, slow_rate.max_dvdt_last_beat);
}
