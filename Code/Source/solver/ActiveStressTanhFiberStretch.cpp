// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "ActiveStressTanhFiberStretch.h"

#include <cmath>

void ActiveStressTanhFiberStretch::init(const unsigned int tnNo) {
  ActiveStress::init(tnNo);

  fourier_interpolation = FourierInterpolation::from_time_series_file(
      temporal_values_file_path, /* n_components = */ 1, ramp);
}

void ActiveStressTanhFiberStretch::read_model_specific_parameters(
    const ActiveStressModelParameters &params) {
  a1 = params.get_scalar("a1");
  a2 = params.get_scalar("a2");
  ramp = params.get_bool("Ramp");
  temporal_values_file_path = params.get_string("Temporal_values_file_path");
}

void ActiveStressTanhFiberStretch::distribute_model_specific_parameters(
    const CmMod &cm_mod, const cmType &cm) {
  cm.bcast(cm_mod, &a1);
  cm.bcast(cm_mod, &a2);
  cm.bcast(cm_mod, &ramp);
  cm.bcast(cm_mod, temporal_values_file_path);
}

double ActiveStressTanhFiberStretch::compute_active_tension_local(
    const Vector<double> &state, const double fiber_stretch) const {
  const double eta = fourier_interpolation.value(time)[0];
  const double phi = std::tanh(a1 * (fiber_stretch - a2));
  return eta * phi * phi;
}

REGISTER_ACTIVE_STRESS_MODEL("TanhFiberStretch", ActiveStressTanhFiberStretch);
