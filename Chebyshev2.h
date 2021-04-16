/**
 * @file Chebyshev2.h
 * @brief Chebyshev parameterizations on Chebyshev points of second kind
 * @author Frank Dellaert
 * @date November 23, 2014
 */

#pragma once

#include "Basis.h"
#include <gtsam/base/Manifold.h>
#include <gtsam/base/OptionalJacobian.h>
#include "Eigen/KroneckerProduct"
// #include <gtsam/3rdparty/Eigen/unsupported/Eigen/KroneckerProduct>
#include <boost/function.hpp>

namespace gtsam
{

/**
 * Chebyshev Interpolation on Chebyshev points of the second kind
 * Note that N here, the #points, is one less than N from Trefethen00book (pg.42)
 */
template <size_t N>
struct Chebyshev2 : Basis<Chebyshev2<N>>
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    typedef Eigen::Matrix<double, 1, N> Weights;
    typedef Eigen::Matrix<double, N, 1> Parameters;

    /**
   * A matrix of M*N values at the Chebyshev points, where M is the dimension of T
   * template argument T: the type you would like to EvaluationFunctor using polynomial interpolation
   * template argument N: the number of Chebyshev points of the second kind
   */
    template <typename T>
    struct ParameterMatrix
    {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        typedef Eigen::Matrix<double, traits<T>::dimension, N> type;
    };

    /// Specific Chebyshev point
    static double Point(int j)
    {
        assert(j >= 0 && size_t(j) < N);
        static double const dtheta = M_PI / (double)(N - 1);
        return cos(-M_PI + dtheta * j);
        // sin(- M_PI_2 + dtheta*j); also works
    }

    /// Specific Chebyshev point within [a,b] interval
    static double Point(int j, double a, double b)
    {
        assert(j >= 0 && size_t(j) < N);
        static double const dtheta = M_PI / (double)(N - 1);
        return a + (b - a) * (1. + cos(-M_PI + dtheta * j)) / 2;
    }

    /// All Chebyshev points
    static Eigen::Matrix<double, N, 1> Points()
    {
        Eigen::Matrix<double, N, 1> points;
        for (size_t j = 0; j < N; j++)
            points(j) = Point(j);
        return points;
    }

    /// All Chebyshev points within [a,b] interval
    static Eigen::Matrix<double, N, 1> Points(double a, double b)
    {
        Eigen::Matrix<double, N, 1> points = Points();
        const double T1 = (a + b) / 2, T2 = (b - a) / 2;
        for (size_t j = 0; j < N; j++)
            points(j) = T1 + T2 * points(j);
        return points;
    }

    /**
   * Evaluate Chebyshev Weights on [-1,1] at any x up to order N-1 (N values)
   * These weights implement barycentric inerpolation at a specific x.
   * More precisely, f(x) ~ [w0;...;wN] * [f0;...;fN], where the fj are the values
   * of the function f at the Chebyshev points. As such, for a given x
   * we obtain a linear map from parameter vectors f to interpolated values f(x).
   * Optional [a,b] interval can be specified as well.
   */
    static Weights CalculateWeights(double x, double a = -1, double b = 1)
    {
        // Trefethen13chapters p 34, formula 5.11

        // Allocate space for weights
        Weights weights;

        // We start by getting distances from x to all Chebyshev points
        // as well as getting smallest distance
        Weights distances;
        for (size_t j = 0; j < N; j++)
        {
            const double dj = x - Point(j, a, b); // only thing that depends on [a,b]
            if (std::abs(dj) < 1e-10)
            {
                // exceptional case: x coincides with a Chebyshev point
                weights.setZero();
                weights(j) = 1;
                return weights;
            }
            distances(j) = dj;
        }

        // Beginning of interval, j = 0, x0 = a
        weights(0) = 0.5 / distances(0);

        // All intermediate points j=1:N-2
        double d = weights(0), s = -1; // changes sign s at every iteration
        for (size_t j = 1; j < N - 1; j++, s = -s)
        {
            weights(j) = s / distances(j);
            d += weights(j);
        }

        // End of interval, j = N-1, x0 = 1.0
        weights(N - 1) = 0.5 * s / distances(N - 1);
        d += weights(N - 1);

        // normalize
        return weights / d;
    }

    /**
   *  Evaluate derivative weights.
   *  NOTE(duy): This would be nice to have but doesn't work yet!
   */
    static Weights DerivativeWeights(double x, double a = -1, double b = 1)
    {
        // Trefethen13chapters p 34, formula 5.11

        // We start by getting distances from x to all Chebyshev points
        // as well as getting smallest distance
        Weights distances;
        for (size_t j = 0; j < N; j++)
            distances(j) = x - Point(j, a, b);

        // Allocate space for weights
        Weights weights;

        // Beginning of interval, j = 0, x0 = -1.0
        weights(0) = -0.5 / pow(distances(0), 2);

        // All intermediate points j=1:N-2
        double d = weights(0), s = -1; // changes sign s at every iteration
        for (size_t j = 1; j < N - 1; j++, s = -s)
        {
            weights(j) = -s / pow(distances(j), 2);
            d += weights(j);
        }

        // End of interval, j = N-1, x0 = 1.0
        weights(N - 1) = -0.5 * s / pow(distances(N - 1), 2);
        d += weights(N - 1);

        // normalize
        return weights / d;
    }

    /**
   *  Evaluate Clenshaw-Curtis integration weights.
   *  Trefethen00book, pg 128, clencurt.m
   *  Note that N in clencurt.m is 1 less than our N
   *  K = N-1;
      theta = pi*(0:K)'/K;
      w = zeros(1,N); ii = 2:K; v = ones(K-1, 1);
      if mod(K,2) == 0
          w(1) = 1/(K^2-1); w(N) = w(1);
          for k=1:K/2-1, v = v-2*cos(2*k*theta(ii))/(4*k^2-1); end
          v = v - cos(K*theta(ii))/(K^2-1);
      else
          w(1) = 1/K^2; w(N) = w(1);
          for k=1:K/2, v = v-2*cos(2*k*theta(ii))/(4*k^2-1); end
      end
      w(ii) = 2*v/K;

   */
    static Weights IntegrationWeights(double a = -1, double b = 1)
    {
        // Allocate space for weights
        Weights weights;
        size_t K = N - 1, // number of intervals between N points
            K2 = K * K;
        weights(0) = 0.5 * (b - a) / (K2 + K % 2 - 1);
        weights(N - 1) = weights(0);
        size_t last_k = K / 2 + K % 2 - 1;
        for (size_t ii = 1; ii <= N - 2; ++ii)
        {
            double theta = ii * M_PI / K;
            weights(ii) = (K % 2 == 0) ? 1 - cos(K * theta) / (K2 - 1) : 1;
            for (size_t k = 1; k <= last_k; ++k)
                weights(ii) -= 2 * cos(2 * k * theta) / (4 * k * k - 1);
            weights(ii) *= (b - a) / K;
        }
        return weights;
    }

    /// Call weights for all x in vector, with interval
    /// This is mainly for help with plotting trajectories in MATLAB
    static Matrix WeightMatrix(const Vector &X, double a = -1, double b = 1)
    {
        Matrix W(X.size(), N);
        for (int i = 0; i < X.size(); i++)
            W.row(i) = CalculateWeights(X(i), a, b);
        return W;
    }

    /// Given values f at the Chebyshev points, predict value f(x) at x
    /// Optional [a,b] interval
    class EvaluationFunctor
    {
    protected:
        Weights weights_;

    public:
        /// Constructor with standard interval [-1,1]
        EvaluationFunctor(double x) : weights_(CalculateWeights(x))
        {
        }
        /// Constructor with non-standard interval [a,b]
        EvaluationFunctor(double x, double a, double b) : weights_(CalculateWeights(x, a, b))
        {
        }
        /// Regular 1D interpolation
        double apply(const Parameters &f, //
                     OptionalJacobian<1, N> H = boost::none) const
        {
            if (H)
                *H = weights_;
            return weights_ * f;
        }
        /// c++ sugar
        double operator()(const Parameters &f, //
                          OptionalJacobian<1, N> H = boost::none) const
        {
            return apply(f, H); // might call apply in derived
        }
    };

    /// Vector interpolation, e.g. given 3*N matrix yields 3-vector
    template <int M>
    class VectorEvaluationFunctor
    {
    protected:
        typedef Eigen::Matrix<double, M, 1> VectorM;
        typedef Eigen::Matrix<double, M, N> MatrixMN;
        typedef Eigen::Matrix<double, M, M * N> Jacobian;
        Jacobian H_;
        Weights weights_;
        void calculateJacobian()
        {
            typedef Eigen::Matrix<double, M, M> MatrixM;
            H_ = Eigen::kroneckerProduct(weights_, MatrixM::Identity());
        }

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        typedef MatrixMN Parameters;

        /// Constructor, standard interval
        VectorEvaluationFunctor(double x) : weights_(CalculateWeights(x))
        {
            calculateJacobian();
        }
        /// Constructor, with non-standard interval [a,b]
        VectorEvaluationFunctor(double x, double a, double b) :

                                                                weights_(CalculateWeights(x, a, b))
        {
            calculateJacobian();
        }
        VectorM apply(const MatrixMN &F, //
                      OptionalJacobian<M, M *N> H = boost::none) const
        {
            if (H)
                *H = H_;
            return F * weights_.transpose();
        }
        /// c++ sugar
        VectorM operator()(const MatrixMN &F, //
                           OptionalJacobian<M, M *N> H = boost::none) const
        {
            return apply(F, H);
        }
    };

    /// Manifold interpolation
    template <class T>
    class type : public VectorEvaluationFunctor<traits<T>::dimension>
    {
        enum
        {
            M = traits<T>::dimension
        };
        typedef VectorEvaluationFunctor<M> Base;

    public:
        /// Constructor, possibly with standard interval
        type(double x) : Base(x)
        {
        }
        /// Constructor, with non-standard interval [a,b]
        type(double x, double a, double b) : Base(x, a, b)
        {
        }
        /// Manifold interpolation
        T apply(const typename Base::MatrixMN &F, //
                OptionalJacobian<M, M *N> H = boost::none) const
        {
            // Interpolate the M-dimensional vector to yield a vector in tangent space
            Eigen::Matrix<double, M, 1> xi = Base::operator()(F, H);
            // Now call retract with this M-vector, possibly with derivatives
            Eigen::Matrix<double, M, M> D_result_xi;
            T result = T::ChartAtOrigin::Retract(xi, H ? &D_result_xi : 0);
            // Finally, if derivatives are asked, apply chain rule
            // where H is Mx(M*N) derivative of interpolation
            // and D_result_xi is MxM derivative of retract
            if (H)
                *H = D_result_xi * (*H);
            // and return a T
            return result;
        }
        /// c++ sugar
        T operator()(const typename Base::MatrixMN &F, //
                     OptionalJacobian<M, M *N> H = boost::none) const
        {
            return apply(F, H); // might call apply in derived
        }
    };

    /// Given M*N Matrix at Chebyshev points, predict component for given row
    template <int M>
    struct one_component : public EvaluationFunctor
    {
    protected:
        typedef Eigen::Matrix<double, M, N> MatrixMN;
        typedef Eigen::Matrix<double, 1, M * N> Jacobian;
        size_t row_;
        Jacobian H_;
        void calculateJacobian()
        {
            H_.setZero();
            for (int j = 0; j < EvaluationFunctor::weights_.size(); j++)
                H_(0, row_ + j * M) = EvaluationFunctor::weights_(j);
        }

    public:
        /// Construct with row index
        one_component(size_t i, double x, double a = -1, double b = 1) : EvaluationFunctor(x, a, b), row_(i)
        {
            calculateJacobian();
        }
        /// Calculate component of component row_ of F
        double apply(const MatrixMN &F, //
                     OptionalJacobian<1, M *N> H = boost::none) const
        {
            if (H)
                *H = H_;
            return F.row(row_) * EvaluationFunctor::weights_.transpose();
        }
        /// c++ sugar
        double operator()(const MatrixMN &F, //
                          OptionalJacobian<1, M *N> H = boost::none) const
        {
            return apply(F, H);
        }
    };

    /// compute D = differentiation matrix, Trefethen00book p.53
    /// when given a parameter vector f of function values at the Chebyshev points,
    /// D*f are the values of f'.
    typedef Eigen::Matrix<double, N, N> DiffMatrix;
    static DiffMatrix DifferentiationMatrix(double a = -1, double b = 1)
    {
        DiffMatrix D;
        if (N == 1)
        {
            D(0, 0) = 1;
            return D;
        }
        for (size_t i = 0; i < N; i++)
        {
            double xi = Point(i, a, b);
            double ci = (i == 0 || i == N - 1) ? 2. : 1.;
            int s = pow(-1, i);
            for (size_t j = 0; j < N; j++, s = -s)
            {
                if (i != j)
                {
                    double xj = Point(j, a, b);
                    double cj = (j == 0 || j == N - 1) ? 2. : 1.;
                    D(i, j) = (ci / cj) * s / (xi - xj);
                }
            }
            D(i, i) = 0.0;
            for (size_t j = 0; j < N; j++)
                if (i != j)
                    D(i, i) -= D(i, j);
        }
        return D;
    }

    /// Base class for functors below that calculates weights
    class DerivativeFunctorBase
    {
    protected:
        Weights weights_;

    public:
        DerivativeFunctorBase(double x) :
#ifdef CHEBYSHEV_WHAT_WE_SHOULD_BE_USING
                                          weights_(DerivativeWeights(x))
#else
                                          weights_(CalculateWeights(x) * DifferentiationMatrix())
#endif
        {
        }

        DerivativeFunctorBase(double x, double a, double b) :
#ifdef CHEBYSHEV_WHAT_WE_SHOULD_BE_USING
                                                              weights_(DerivativeWeights(x, a, b))
#else
                                                              weights_(CalculateWeights(x, a, b) * DifferentiationMatrix(a, b))
#endif
        {
        }
    };

    /// Given values f at the Chebyshev points, predict derivative at x
    struct DerivativeFunctor : protected DerivativeFunctorBase
    {
        DerivativeFunctor(double x) : DerivativeFunctorBase(x)
        {
        }
        DerivativeFunctor(double x, double a, double b) : DerivativeFunctorBase(x, a, b)
        {
        }
        double apply(const Parameters &f, //
                     OptionalJacobian<1, N> H = boost::none) const
        {
            if (H)
                *H = this->weights_;
            return this->weights_ * f;
        }
        /// c++ sugar
        double operator()(const Parameters &f, //
                          OptionalJacobian<1, N> H = boost::none) const
        {
            return apply(f, H); // might call apply in derived
        }
    };

    /// Vector interpolation, e.g. given 3*N matrix yields 3-vector
    template <int M>
    class vec_derivative : protected DerivativeFunctorBase
    {
    protected:
        typedef Eigen::Matrix<double, 1, M> VectorM;
        typedef Eigen::Matrix<double, M, N> MatrixMN;
        typedef Eigen::Matrix<double, M, M * N> Jacobian;
        Jacobian H_;
        void calculateJacobian()
        {
            typedef Eigen::Matrix<double, M, M> MatrixM;
            H_ = Eigen::kroneckerProduct(this->weights_, MatrixM::Identity());
        }

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        /// Constructor, standard interval
        vec_derivative(double x) : DerivativeFunctorBase(x)
        {
            calculateJacobian();
        }
        /// Constructor, with non-standard interval [a,b]
        vec_derivative(double x, double a, double b) : DerivativeFunctorBase(x, a, b)
        {
            calculateJacobian();
        }
        VectorM apply(const MatrixMN &F, //
                      OptionalJacobian<M, M *N> H = boost::none) const
        {
            if (H)
                *H = H_;
            return F * this->weights_.transpose();
        }
        /// c++ sugar
        VectorM operator()(const MatrixMN &F, //
                           OptionalJacobian<M, M *N> H = boost::none) const
        {
            return apply(F, H);
        }
    };

    /// Given M*N Matrix at Chebyshev points, predict derivative for given row
    template <int M>
    struct one_derivative : protected DerivativeFunctorBase
    {
    protected:
        typedef Eigen::Matrix<double, M, N> MatrixMN;
        typedef Eigen::Matrix<double, 1, M * N> Jacobian;
        size_t row_;
        Jacobian H_;
        void calculateJacobian()
        {
            H_.setZero();
            for (int j = 0; j < this->weights_.size(); j++)
                H_(0, row_ + j * M) = this->weights_(j);
        }

    public:
        /// Construct with row index
        one_derivative(size_t i, double x, double a = -1, double b = 1) : DerivativeFunctorBase(x, a, b), row_(i)
        {
            calculateJacobian();
        }
        /// Calculate derivative of component row_ of F
        double apply(const MatrixMN &F, //
                     OptionalJacobian<1, M *N> H = boost::none) const
        {
            if (H)
                *H = H_;
            return F.row(row_) * this->weights_.transpose();
        }
        /// c++ sugar
        double operator()(const MatrixMN &F, //
                          OptionalJacobian<1, M *N> H = boost::none) const
        {
            return apply(F, H);
        }
    };

    // Vector version for MATLAB :-(
    static double Derivative(double x, const Vector &f, //
                             OptionalJacobian<1, N> H = boost::none)
    {
        return DerivativeFunctor(x)(f.transpose(), H);
    }

    /**
   * Create matrix of values at Chebyshev points given vector-valued function.
   */
    template <size_t M>
    static Eigen::Matrix<double, M, N> matrix(
        boost::function<Eigen::Matrix<double, M, 1>(double)> f, //
        double a = -1, double b = 1)
    {
        Eigen::Matrix<double, 12, N> Xmat;
        for (size_t j = 0; j < N; j++)
            Xmat.col(j) = f(Point(j, a, b));
        return Xmat;
    }
};
// \ Chebyshev2

typedef Chebyshev2<7> Chebyshev27;    // TODO delete
typedef Chebyshev2<3> Chebyshev23;    // TODO delete
typedef Chebyshev2<16> Chebyshev2_16; // TODO delete
typedef Chebyshev2<32> Chebyshev2_32; // TODO delete
typedef Chebyshev2<64> Chebyshev2_64; // TODO delete

} // namespace gtsam
