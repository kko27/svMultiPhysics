// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef RIGID_PLANE_BOUNDARY_CONDITION_H
#define RIGID_PLANE_BOUNDARY_CONDITION_H

#include "CmMod.h"
#include "SimulationLogger.h"
#include "Vector.h"

#include <vector>

/// @brief Boundary-condition data for unilateral rigid-plane penalty contact.
class RigidPlaneBoundaryCondition
{
public:
    RigidPlaneBoundaryCondition() = default;

    RigidPlaneBoundaryCondition(double penalty_factor, const std::vector<double>& plane_point,
        const std::vector<double>& plane_normal, int nsd, SimulationLogger& logger);

    bool is_initialized() const noexcept { return initialized_; }

    double get_penalty_factor() const noexcept { return penalty_factor_; }

    const Vector<double>& get_plane_point() const noexcept { return plane_point_; }

    const Vector<double>& get_plane_normal() const noexcept { return plane_normal_; }

    /// @brief Signed gap g = n·(y − p) − d of a point y against the plane
    /// translated by d along its unit normal n. g > 0 means penetration.
    double gap(const Vector<double>& y, double plane_displacement) const;

    /// @brief Distribute model parameters to all parallel processes.
    void distribute(const CmMod& cm_mod, const cmType& cm, int nsd);

private:
    bool initialized_ = false;
    double penalty_factor_ = 0.0;
    Vector<double> plane_point_;
    Vector<double> plane_normal_;
    const SimulationLogger* logger_ = nullptr;
};

#endif
