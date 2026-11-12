// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the 
// University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef RIGID_PLANE_BOUNDARY_CONDITION_H
#define RIGID_PLANE_BOUNDARY_CONDITION_H

#include "CmMod.h"
#include "SimulationLogger.h"
#include "Vector.h"

#include <vector>

/// @brief Boundary-condition data for unilateral rigid-plane penalty contact.
///
/// A face is kept on one side of a rigid plane that passes through the
/// point p with unit normal n. The plane can move along n by a prescribed scalar
/// displacement d(t), read from the BC's temporal values file. The rigid obstacle
/// occupies the half-space on the +n side of the plane. Points on that side have
/// a positive gap g (see gap()) and are penetrating.
///
/// At each face quadrature point where g > 0, the face gets a restoring traction
/// of magnitude penalty_factor * g, directed along -n (pushed back out of the
/// obstacle). Points with g <= 0 feel no force, so the contact is unilateral:
/// the plane can push the face but never pull it. Assembly is done in
/// set_bc_rigid_plane_l() (set_bc.cpp).
///
/// This class stores only the plane geometry and the penalty parameter. The
/// time history d(t) is kept in the owning bcType (bcType::gt).
///
/// XML parameters: <Penalty_factor>, <Plane_point>, <Plane_normal>.
class RigidPlaneBoundaryCondition
{
public:
    /// @brief Construct an inactive BC; is_initialized() returns false.
    RigidPlaneBoundaryCondition() = default;

    /// @brief Construct and validate a rigid-plane BC.
    ///
    /// @param penalty_factor Penalty stiffness, as traction per unit gap
    ///        (force / length^3). Must be > 0.
    /// @param plane_point    A point p on the plane at d = 0. Must have nsd components.
    /// @param plane_normal   Plane normal pointing into the rigid obstacle. Need
    ///                       not be unit length; it is normalized here. Must have
    ///                       nsd components and non-zero length.
    /// @param nsd            Number of spatial dimensions.
    /// @param logger         Logger for the initialization summary.
    ///
    /// @throws svmp::CoreException (StatusCode::InvalidArgument) if any argument
    ///         is invalid.
    RigidPlaneBoundaryCondition(double penalty_factor, const std::vector<double>& plane_point,
        const std::vector<double>& plane_normal, int nsd, SimulationLogger& logger);

    /// @brief True if this BC was constructed from parameters, or received them
    /// through distribute(). A default-constructed BC returns false.
    bool is_initialized() const noexcept { return initialized_; }

    double get_penalty_factor() const noexcept { return penalty_factor_; }

    /// @brief Point p on the plane at zero plane displacement.
    const Vector<double>& get_plane_point() const noexcept { return plane_point_; }

    /// @brief Unit normal n, pointing into the rigid obstacle.
    const Vector<double>& get_plane_normal() const noexcept { return plane_normal_; }

    /// @brief Signed gap g = n·(y − p) − d of a point y against the plane
    /// translated by d along its unit normal n. g > 0 means penetration.
    ///
    /// @param y                  Current (deformed) position of the point.
    /// @param plane_displacement Plane displacement d along n at the current time.
    double gap(const Vector<double>& y, double plane_displacement) const;

    /// @brief Distribute model parameters to all parallel processes.
    void distribute(const CmMod& cm_mod, const cmType& cm, int nsd);

private:
    bool initialized_ = false;
    double penalty_factor_ = 0.0;
    Vector<double> plane_point_;   ///< Point p on the plane at d = 0.
    Vector<double> plane_normal_;  ///< Unit normal n, pointing into the obstacle.
    const SimulationLogger* logger_ = nullptr;  ///< Not distributed; null on non-master ranks.
};

#endif
