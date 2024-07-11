#ifndef FIVO_TIME_H
#define FIVO_TIME_H

#include <functional>
#include <cmath>

namespace fivo {

namespace time {

struct RK {
  auto stages() const { return c.size(); }
  auto butcher() const { return std::make_tuple(std::cref(a), std::cref(b), std::cref(c)); }

  template<typename State, typename RHS>
  void operator()(State& X, double const& t, double const& dt, RHS const& f) const {
    std::vector<State> k(stages());
    for (std::size_t i = 0; i < stages(); ++i) {
      State tmp(X.size());
      for (std::size_t j = 0; j < i; ++j) { tmp += a[i * stages() + j] * k[j]; }
      k[i] = f(t + c[i] * dt, X + dt * tmp);
    }
    for (std::size_t i = 0; i < stages(); ++i) { X += dt * b[i] * k[i]; }
  }

protected:
  explicit RK(std::size_t stages) : a(stages * stages), b(stages), c(stages) {}
  std::vector<double> a, b, c;
};

/* FIRST ORDER METHODS */
struct RK1 : RK {
  RK1() : RK(1) {a = {0}, b = {1}, c = {0}; }
};

/* SECOND ORDER METHODS */
struct RK2Heun : RK {
  RK2Heun() : RK(2) {
    a = {0 , 0,
         1 , 0};
    b = {0.5, 0.5};
    c = {0, 1};
  }
};

struct RK2Midpoint : RK {
  RK2Midpoint() : RK(2) {
    a = {0   , 0,
         0.5 , 0};
    b = {0, 1};
    c = {0, 0.5};
  }
};

struct RK2Ralston : RK {
  RK2Ralston() : RK(2) {
    a = {0     , 0,
         2./3. , 0};
    b = {0.25, 0.75};
    c = {0, 2./3.};
  }
};

struct RK2Generic : RK {
  explicit RK2Generic(double alpha) : RK(2) {
    a = {0     , 0,
         alpha , 0};
    b = {1 - 1 / (2 * alpha), 1 / (2 * alpha)};
    c = {0, alpha};
  }
};

/* THIRD ORDER METHODS */
struct RK3Kutta : RK {
  RK3Kutta() : RK(3) {
    a = {0   , 0 , 0,
         0.5 , 0 , 0,
         -1  , 2 , 0};
    b = {1./6., 2./3., 1./6.};
    c = {0, 0.5, 1};
  }
};

struct RK3Heun : RK {
  RK3Heun() : RK(3) {
    a = {0     , 0     , 0,
         1./3. , 0     , 0,
         0     , 2./3. , 0};
    b = {1./4., 0, 3./4.};
    c = {0, 1./3., 2./3.};
  }
};

struct RK3VHW : RK {
  RK3VHW() : RK(3) {
    a = {0      , 0      , 0,
         8./15. , 0      , 0,
         1./4.  , 5./12. , 0};
    b = {1./4., 0, 3./4.};
    c = {0, 8./15., 2./3.};
  }
};

struct RK3Ralston : RK {
  RK3Ralston() : RK(3) {
    a = {0     , 0     , 0,
         1./2. , 0     , 0,
         0     , 3./4. , 0};
    b = {2./9., 1/3., 4./9.};
    c = {0, 1./2., 3./4.};
  }
};

struct RK3SSP : RK {
  RK3SSP() : RK(3) {
    a = {0     , 0     , 0,
         1.    , 0     , 0,
         1./4. , 1./4. , 0};
    b = {1./6., 1./6., 2./3.};
    c = {0, 1., 1./2.};
  }
};

struct RK3Generic : RK {
  explicit RK3Generic(double alpha) : RK(3) {
    auto const f1 = (1-alpha)/(alpha*(3*alpha-2));
    auto const f2 = 1/(6*(1-alpha));
    a = {0      , 0   , 0,
         alpha  , 0   , 0,
         1 + f1 , -f1 , 0};
    b = {1./2. - 1/(6*alpha), f2/alpha, (2 - 3 * alpha) * f2};
    c = {0, alpha, 1};
  }
};

/* FOURTH ORDER METHODS */
struct RK4 : RK {
  RK4() : RK(4) {
    a = {0   , 0   , 0 , 0,
         0.5 , 0   , 0 , 0,
         0   , 0.5 , 0 , 0,
         0   , 0   , 1 , 0};
    b = {1./6., 1./3., 1./3., 1./6.};
    c = {0, 1./2., 1./2., 1};
  }
};

struct RK438 : RK {
  RK438() : RK(4) {
    a = {0      , 0  , 0 , 0,
         1./3.  , 0  , 0 , 0,
         -1./3. , 1  , 0 , 0,
         1      , -1 , 1 , 0};
    b = {1./8., 3./8., 3./8., 1./8.};
    c = {0, 1./3., 2./3., 1};
  }
};

// struct EmbeddedRK {
//   template<typename State>
//   auto error(State const& X, State const& Xh, State const& X0) const {
//     return 0;
//   }

//   void adapt(double err, int p, int ph, double& dt) const {
//     auto const q = std::min(p, ph);
//     dt = dt * std::min(facmax, std::max(facmin, fac * std::pow(1 / err, 1. / (q + 1))));
//   }

//   template<typename State, typename RHS>
//   void operator()(State& X, double const& t, double& dt, RHS const& f) const;

// protected:
//   EmbeddedRK(double abs, double rel, double fac, double facmax, double facmin) :
//     abs(abs), rel(rel), fac(fac), facmax(facmax), facmin(facmin)
//   {}
//   double abs, rel, fac, facmax, facmin;
// };

// struct RK2Heun_1 : EmbeddedRK {
//   template<typename State, typename RHS>
//   void operator()(State& X, double const& t, double& dt, RHS const& f) const {
//     State X0 = X;
//     double err, maxsc;
//     do {
//       // High order
//       auto const k1 = f(t, X);
//       auto const k2 = f(t + dt, X + dt * k1);
//       X += dt * 0.5 * (k1 + k2);
//       // Low order
//       auto Xh = X + dt * k1;
//       // Estimate the error and compute new step size
//       err = error(X, Xh, X0);
//       adapt(err, 2, 1, dt);
//     } while (err > maxsc);
//   }
// };

} // namespace time

} // namespace fivo

#endif // FIVO_TIME_H
