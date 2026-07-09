// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "RigidPlaneBoundaryCondition.h"

#include <cmath>
#include <stdexcept>

RigidPlaneBoundaryCondition::RigidPlaneBoundaryCondition(double penalty_factor,
    const std::vector<double>& plane_point, const std::vector<double>& plane_normal,
    SimulationLogger& logger)
    : penalty_factor_(penalty_factor)
    , plane_point_(plane_point.size())
    , plane_normal_(plane_normal.size())
    , logger_(&logger)
{
    if (penalty_factor_ <= 0.0) {
        throw std::runtime_error("[RigidPlaneBoundaryCondition] Penalty_factor must be positive.");
    }

    if (plane_point.size() == 0 || plane_point.size() != plane_normal.size()) {
        throw std::runtime_error("[RigidPlaneBoundaryCondition] Plane_point and Plane_normal must have the same non-zero size.");
    }

    double normal_norm_sq = 0.0;
    for (int i = 0; i < plane_point_.size(); i++) {
        plane_point_(i) = plane_point[i];
        plane_normal_(i) = plane_normal[i];
        normal_norm_sq += plane_normal_(i) * plane_normal_(i);
    }

    double normal_norm = std::sqrt(normal_norm_sq);
    if (normal_norm <= 0.0) {
        throw std::runtime_error("[RigidPlaneBoundaryCondition] Plane_normal must have non-zero length.");
    }

    for (int i = 0; i < plane_normal_.size(); i++) {
        plane_normal_(i) /= normal_norm;
    }
    initialized_ = true;

    logger_->log_message("[RigidPlaneBoundaryCondition] Initialized rigid-plane penalty BC");
    logger_->log_message("\t Penalty factor:", penalty_factor_);
}

void RigidPlaneBoundaryCondition::distribute(const CmMod& cm_mod, const cmType& cm, int nsd)
{
    cm.bcast(cm_mod, &initialized_);
    if (!initialized_) {
        return;
    }

    cm.bcast(cm_mod, &penalty_factor_);

    if (plane_point_.size() == 0) {
        plane_point_.resize(nsd);
    }
    if (plane_normal_.size() == 0) {
        plane_normal_.resize(nsd);
    }

    cm.bcast(cm_mod, plane_point_);
    cm.bcast(cm_mod, plane_normal_);
}
