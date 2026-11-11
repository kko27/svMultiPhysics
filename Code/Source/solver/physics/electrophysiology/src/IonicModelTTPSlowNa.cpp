// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "IonicModelTTPSlowNa.h"

void TTPSlowNa::read_parameters(const IonicModelParameters &params) {
  TTP::read_parameters(params);

  usi_param = params.get_scalar("usi_param");
  V_half_k = params.get_scalar("V_half_k");
  slope_k = params.get_scalar("slope_k");
}

void TTPSlowNa::distribute_parameters(const CmMod &cm_mod, const cmType &cm) {
  TTP::distribute_parameters(cm_mod, cm);

  cm.bcast(cm_mod, &usi_param);
  cm.bcast(cm_mod, &V_half_k);
  cm.bcast(cm_mod, &slope_k);
}

void TTPSlowNa::update_g(const unsigned int zone_id, const double dt,
                         const Vector<double> &X, Vector<double> &Xg) const {
  // TTP::update_g only reads and writes Xg(0)..Xg(11), so Xg(k_index) still
  // holds the old value here.
  TTP::update_g(zone_id, dt, X, Xg);

  // k: ultra-slow inactivation gate for I_Na
  const double V = X(0);
  const double k = Xg(k_index);
  const double ki = 1.0 / (1.0 + exp((V - V_half_k) / slope_k));
  const double tau = 5.0 + usi_param / (1.0 + exp((V + 15.0) / 4.0));

  Xg(k_index) = ki - (ki - k) * exp(-dt / tau);
}

double TTPSlowNa::sodium_current_gate_factor(const Vector<double> &Xg) const {
  return pow(Xg(k_index), 3.0);
}

REGISTER_IONIC_MODEL("TTP_slowNa", TTPSlowNa);
