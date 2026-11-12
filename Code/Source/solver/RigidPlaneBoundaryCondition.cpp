// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "RigidPlaneBoundaryCondition.h"

#include "Core/Exception.h"

#include <cmath>

namespace {
using svmp::CoreException;
constexpr auto kInvalidArg = svmp::StatusCode::InvalidArgument;
}

RigidPlaneBoundaryCondition::RigidPlaneBoundaryCondition(double penalty_factor,
    const std::vector<double>& plane_point, const std::vector<double>& plane_normal,
    int nsd, SimulationLogger& logger)
    : penalty_factor_(penalty_factor)
    , plane_point_(plane_point.size())
    , plane_normal_(plane_normal.size())
    , logger_(&logger)
{
    svmp::check<CoreException>(penalty_factor_ > 0.0,
        "Penalty_factor must be positive.", kInvalidArg);

    svmp::check<CoreException>(plane_point.size() == static_cast<size_t>(nsd) &&
        plane_normal.size() == static_cast<size_t>(nsd),
        "Plane_point and Plane_normal must match the number of spatial dimensions (" +
        std::to_string(nsd) + ").", kInvalidArg);

    double normal_norm_sq = 0.0;
    for (int i = 0; i < plane_point_.size(); i++) {
        plane_point_(i) = plane_point[i];
        plane_normal_(i) = plane_normal[i];
        normal_norm_sq += plane_normal_(i) * plane_normal_(i);
    }

    double normal_norm = std::sqrt(normal_norm_sq);
    svmp::throw_if<CoreException>(normal_norm <= 0.0,
        "Plane_normal must have non-zero length.", kInvalidArg);

    for (int i = 0; i < plane_normal_.size(); i++) {
        plane_normal_(i) /= normal_norm;
    }
    initialized_ = true;

    logger_->log_message("[RigidPlaneBoundaryCondition] Initialized rigid-plane penalty BC");
    logger_->log_message("\t Penalty factor:", penalty_factor_);
}

double RigidPlaneBoundaryCondition::gap(const Vector<double>& y, double plane_displacement) const
{
    double g = -plane_displacement;
    for (int i = 0; i < plane_normal_.size(); i++) {
        g += plane_normal_(i) * (y(i) - plane_point_(i));
    }
    return g;
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
