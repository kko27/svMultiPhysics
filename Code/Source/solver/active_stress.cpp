// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "active_stress.h"

bool supports_active_stress(const consts::EquationType eq_type) {
  return eq_type == consts::EquationType::phys_struct ||
         eq_type == consts::EquationType::phys_ustruct ||
         eq_type == consts::EquationType::phys_FSI;
}

void ActiveStress::read_parameters(const ActiveStressParameters &params) {
  eta_f = params.get_eta_f();
  eta_s = params.get_eta_s();
  eta_n = params.get_eta_n();

  use_stabilization = params.get_use_stabilization(); 

  read_model_specific_parameters(
      params.get_parameters(params.get_model_name()));
}

void ActiveStress::distribute_parameters(const CmMod &cm_mod,
                                         const cmType &cm) {
  cm.bcast(cm_mod, &eta_f);
  cm.bcast(cm_mod, &eta_s);
  cm.bcast(cm_mod, &eta_n);
  cm.bcast(cm_mod, &use_stabilization);

  distribute_model_specific_parameters(cm_mod, cm);
}

void ActiveStress::init(const unsigned int tnNo) {
  states.resize(n_states, tnNo);

  if (n_states > 0) {
    Vector<double> state_loc(n_states);
    init_local(state_loc);

    for (unsigned int i = 0; i < tnNo; ++i)
      for (unsigned int j = 0; j < n_states; ++j)
        states(j, i) = state_loc(j);
  }

  active_tension.resize(tnNo);
  raw_active_tension.resize(tnNo); 
  previous_fiber_stretch.resize(tnNo); 
  has_previous_fiber_stretch.resize(tnNo); 

  for (unsigned int i = 0; i < tnNo; ++i) {
    active_tension[i] = 0.0;  
    raw_active_tension[i] = 0.0;
    previous_fiber_stretch[i] = 1.0;
    has_previous_fiber_stretch[i] = 0;
  }
}

void ActiveStress::advance_time_step(const double t, const double dt,
                                     const Vector<double> &calcium,
                                     const Vector<double> &fiber_stretch,
                                     const Vector<double> &fiber_stretch_rate) {
  time = t;

  for (unsigned int i = 0; i < states.ncols(); ++i) {
    Vector<double> state_loc = states.col(i);
    advance_time_step_local(t, dt, calcium[i], fiber_stretch[i],
                            fiber_stretch_rate[i], state_loc);
    states.set_col(i, state_loc);

    const double Ta = compute_active_tension_local(state_loc); 
    raw_active_tension[i] = Ta; 

    if (!has_previous_fiber_stretch[i]) {
      previous_fiber_stretch[i] = fiber_stretch[i];
      has_previous_fiber_stretch[i] = 1;
    }

    if (use_stabilization) {
      const double Ka = 
        compute_active_stiffness_local(state_loc, fiber_stretch[i], 
                                       fiber_stretch_rate[i]);

      active_tension[i] = apply_stabilization_local(Ta, Ka, fiber_stretch[i], 
                                                previous_fiber_stretch[i]);
    } else { 
      active_tension[i] = Ta; 
    }

    previous_fiber_stretch[i] = fiber_stretch[i];
  }
}