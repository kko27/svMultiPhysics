Example structural case for the `RigidPlanePenalty` boundary condition.

It reuses the `block_compression` mesh and applies unilateral penalty contact on face `Z1`
against a rigid plane defined by `Plane_point` and `Plane_normal`, with prescribed plane
translation from `plate_motion.dat`.
