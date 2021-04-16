/**
 * @file Fourier.h
 * @brief Fourier decompositions
 * @author Frank Dellaert
 * @date November 23, 2014
 */

#pragma once

#include <gtsam/base/Matrix.h>

/*
 *  Concept of a Basis:
 *    - CalculateWeights
 *
 *  Concept of SythesizingBasis refines Basis
 *    - Parameters
 *    - EvaluationFunctor(double x)(Parameters c) = weights(x) * c
 *      where c are coefficients for each of the basis function
 *
 *  Concept of InterpolatingBasis refines Basis
 *    - Points, Values
 *    - point(size_t j)
 *    - EvaluationFunctor(double x)(Values f) = weights(x) * f
 *      where f are function values at the point x
 */

namespace gtsam {

/// CRTP Base class for function bases
template <typename Derived>
struct Basis {

  /// Call weights for all x in vector
  static Matrix WeightMatrix(const Vector& X) {
    const size_t m = X.size();
    Matrix row0 = Derived::CalculateWeights(X(0));
    Matrix W(m,row0.cols());
    W.row(0) = row0;
    for (size_t i = 1;i<m;i++)
      W.row(i) = Derived::CalculateWeights(X(i));
    return W;
  }

};

}
