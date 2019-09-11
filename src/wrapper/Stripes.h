#pragma once

#include <Eigen/Dense>

namespace DDG {

///
/// Compute a stripe pattern aligned to an input direction field.
///
/// @param[in]  V                       #V x 3 input mesh vertices.
/// @param[in]  F                       #F x 3 input mesh faces.
/// @param[in]  directionField          #V x 3 input direction field. If empty, compute a curvature
///                                     aligned direction field.
/// @param[in]  frequency               Frequency of output stripes.
/// @param[in]  numCoordinateFunctions  Number of parametrization to compute (D=1 or 2).
/// @param[out] branchIndex             #F x 1 branch index of each parametrization. A zero
///                                     indicates a regular face.
/// @param[out] parameterization        3*#F x D parametrization at each face corner.
/// @param[out] zeroIndex               #F x D zero index of each parametrization.
/// @param[out] isBorder                #F x 1 array indicating if a face is on the border of the
///                                     mesh.
///
void computeStripePatterns(const Eigen::MatrixXd &V,
                           const Eigen::MatrixXi &F,
                           const Eigen::MatrixXd &directionField,
                           double frequency,
                           int numCoordinateFunctions,
                           Eigen::VectorXi &branchIndex,
                           Eigen::MatrixXd &parameterization,
                           Eigen::MatrixXi &zeroIndex,
                           std::vector<bool> &isBorder);

}  // namespace DDG
