---
title: "Description of the fivo finite volume solver"
author: "Gabriel Suau"
date: "05/07/2024"
---

# fivo
Simple FInite VOlumes solver for 1D hyperbolic systems of conservation laws, i.e. systems that can be written in the form

$$ \partial_t U + \nabla\cdot F(U) = S(U) $$

where $U$ is the vector of conservative variables, $F$ is the flux function and $S$ is a source term.

## Highlights
* Header-only library
* Requires C++17
* Extensible

## Systems

Several systems are already implemented in fivo, they are available in the `fivo::system` namespace :
- `LinearAdvection`: advection with constant velocity $a$
  ```math
  U = c \qquad F(U) = ac \qquad S(U) = 0 $$
  ```

- `LWRTrafficFlow`: LWR model for traffic flow
  ```math
  U = \rho \qquad F(U) = \rho u(\rho) \qquad S(U) = 0
  ```
  with EOS $u(\rho) = u_m(1 - \frac{\rho}{\rho_m}$ where $\rho_m$ and $u_m$ are the maximum density and velocity.

- `Burgers`: inviscid Burgers equation
  ```math
  U = u \qquad F(U) = \frac{u^2}{2} \qquad S(U) = 0
  ```

- `LinearAcousticsPressure`: pressure-velocity formulation of the linear acoustics equations
  ```math
  U = \begin{pmatrix} p \\ u \end{pmatrix} \qquad F(U) = \begin{pmatrix} u_0 & \rho_0 c^2 \\ \frac{1}{\rho_0} & u_0 \end{pmatrix} U \qquad S(U) = \begin{pmatrix} 0 \\ 0 \end{pmatrix}
  ```
  with EOS $p(\rho) = \rho c^2$ where $c$ is the speed of sound.

- `LinearAcousticsDensity`: density-velocity formulation of the linear acoustics equations
  ```math
  U = \begin{pmatrix} p \\ u \end{pmatrix} \qquad F(U) = \begin{pmatrix} u_0 & \rho_0 \\ c^2 & u_0 \end{pmatrix} U \qquad S(U) = \begin{pmatrix} 0 \\ 0 \end{pmatrix}
  ```
  with EOS $p(\rho) = \rho c^2$ where $c$ is the speed of sound.

- `SWE`: shallow-water equations
  ```math
  U = \begin{pmatrix} h \\ hu \end{pmatrix} \qquad F(U) = \begin{pmatrix} hu \\ hu^2 + \frac{gh^2}{2} \end{pmatrix} \qquad S(U) = \begin{pmatrix} 0 \\ -gh\partial_x z \end{pmatrix}
  ```
  where $g$ is the gravity acceleration and $z$ is the bed elevation.

- `IsentropicEuler`: simplification of Euler's equations for gas dynamics for isentropic systems
  ```math
  U = \begin{pmatrix} \rho \\ \rho u \end{pmatrix} \qquad F(U) = \begin{pmatrix} \rho u \\ \rho u^2 + p \end{pmatrix} \qquad S(U) = \begin{pmatrix} 0 \\ 0 \end{pmatrix}
  ```
  with EOS $p(\rho) = \kappa \rho^\gamma$, where $\gamma = \frac{c_p}{c_v}$ is the adiabatic index and $\kappa$ a constant.

- `IsothermalEuler`: simplification of Euler's equations for gas dynamics for isothermal systems
  ```math
  U = \begin{pmatrix} \rho \\ \rho u \end{pmatrix} \qquad F(U) = \begin{pmatrix} \rho u \\ \rho u^2 + p \end{pmatrix} \qquad S(U) = \begin{pmatrix} 0 \\ 0 \end{pmatrix}
  ```
  with EOS $p(\rho) = \kappa \rho^\gamma$, where $\gamma = \frac{c_p}{c_v}$ is the adiabatic index and $\kappa$ a constant.

- `Euler`: Euler's equations for gas dynamics
  ```math
  U = \begin{pmatrix} \rho \\ \rho u \\ \epsilon \end{pmatrix} \qquad F(U) = \begin{pmatrix} \rho u \\ \rho u^2 + p \\ u(\epsilon + p) \end{pmatrix} \qquad S(U) = \begin{pmatrix} 0 \\ 0 \end{pmatrix}
  ```
  This is an abstract class that has a pure virtual method to compute the pressure (EOS). At the moment, two derived classes are implemented :
  - `IdealGasEuler`: $p(\rho) = e(\gamma - 1)$, where $\gamma = \frac{c_p}{c_v}$ is the adiabatic index and $e = \epsilon - \frac{\rho u^2}{2}$ is the internal energy per unit volume ($\epsilon$ is the total energy per unit volume)
  - `StiffenedGasEuler`: $p(\rho) = e(\gamma - 1) - \gamma p_0$, where $p_0$ is a pressure correction term that accounts for the additional internal pressure due to molecular interactions.

## Adding a new system

There are two main syntaxes to create a new system : deriving it from one of the above or creating it from scratch.

### Deriving a new system from an existing one

If you only want to change small things to a system implemented in fivo, you can just derive your new system struct from it and only change your application specific implementation. For instance, creating a Euler system with a custom EOS for the pressure is easily done with
```c++
struct EulerCustomEOS : fivo::system::Euler {
  EulerCustomEOS(Mesh const& mesh,
                 std::shared_ptr<Euler::BC> const& left_bc,
                 std::shared_ptr<Euler::BC> const& right_bc,
                 typename fivo::system::Euler::value_type const& gamma,
                 /* optional arguments (needed by your EOS for instance) */)
    : fivo::System::Euler(mesh, left_bc, right_bc, gamma), ... 
  {}
  
  // Define your EOS
  value_type pressure(state_type const& s) const override {
    /*...*/;
  }
};
```

### Add a completely new system
To add a completely new system that can't be derived from on of the existing one, the syntax is a bit more verbose, and is easier to illustrate with an example.

Let's say you'd like to code a new system that has 2 conservative variables, and a source term. The minimal syntax is the following
```c++
using state_type = fivo::state<double, 2>
struct NewSystem : fivo::system::System<state_type> {
  using state_type = fivo::state<double, 2>;
  using global_state_type = fivo::global_state<state_type>;
  using value_type = typename state_type::value_type;

  struct MyBCType : fivo::system::System<state_type>::BC {
    static inline auto make() { return std::make_shared<MyBCType>(); }
    state_type compute(fivo::Mesh const& mesh, value_type const& t, 
                       state_type const& in, state_type const& in_opbound) const final {
      // impl
    }
  };

  // Other BCs...

  /* Define a constructor */
  NewSystem(fivo::Mesh const& mesh, 
            std::shared_ptr<BC> const& left_bc,
            std::shared_ptr<BC> const& right_bc,
            /* optional arguments (needed internally for instance) */)
    : fivo::system::System(mesh, left_bc, right_bc) /*, ... */
  {}

  /* Define the physical flux, wave speeds and the admissible state space */
  state_type flux(state_type const& s) const override { /*...*/; }
  state_type wave_speeds(state_type const& s) const override { /*...*/; }
  bool admissible(state_type const& s) const override { /*...*/; }

  /* Define the global source function */
  global_state_type source(Mesh const& mesh, 
                           value_type const& t, 
                           global_state_type const& X) const override {
    /*...*/;
  }
};
```

Lets break it down. First, your new struct must always derive from the `fivo::system::System` template struct with the appropriate state type. Also, the three typedefs in the struct are needed for fivo solver to work correctly. 

Next, you have to define at least one boundary condition type for your system. This is done by declaring structs with a `compute` method (with the exact same signature as above). You can define them inside or outside your system struct. A minimal constructor shoud have three argument (a mesh and two boundary conditions), and should call the `fivo::system::System` constructor with them.

A new system should at least implement the `flux`, `wave_speeds`, `admissible`, `prim_to_cons` and `cons_to_prim` functions that respectively compute the physical flux, the wave speeds (i.e. the eigenvalues of the flux jacobian), the set of admissible states, and the conversion between primitive and conservative variables. Optionnaly, to add a source, it should also override the `source` function (base implementation sets a zero source everywhere).

## Numerical fluxes

Several numerical fluxes are already implemented in fivo, they are available in the `fivo::flux` namespace :

| System name               | Upwind | Godunov | Rusanov | HLL | HLLC |
|:--------------------------|:------:|---------|---------|-----|------|
| `LinearAdvection`         | Yes    | Yes     | Yes     | Yes | No   |
| `LWRTrafficflow`          | N/A    | Yes     | Yes     | Yes | No   |
| `Burgers`                 | N/A    | Yes     | Yes     | Yes | No   |
| `LinearAcousticsDensity`  | N/A    | Yes     | Yes     | Yes | No   |
| `LinearAcousticsPressure` | N/A    | Yes     | Yes     | Yes | No   |
| `SWE`                     | N/A    | No      | Yes     | Yes | No   |
| `IsothermalEuler`         | N/A    | No      | Yes     | Yes | No   |
| `IsentropicEuler`         | N/A    | No      | Yes     | Yes | No   |
| `Euler` and derived       | N/A    | No      | Yes     | Yes | Yes  |

## Adding a new numerical flux

Numerical fluxes follow the CRTP pattern for static inheritance. To create a new numerical flux, the sytax is the following
```c++
struct NewFlux : fivo::flux::NumericalFlux<NewFlux> {
  // Generic implementation
  template<typename System>
  auto compute(System const& sys,
               typename System::state_type const& left,
               typename System::state_type const& right) const {
    return typename System::state_type{/*...*/};
  }
  
  // Specific implementation for SWE
  auto compute(fivo::system::SWE const& sys,
               typename fivo::system::SWE::state_type const& left,
               typename fivo::system::SWE::state_type const& right) const {
    return typename fivo::system::SWE::state_type{/*...*/};
  }

  // Add as many specialisations as you need...
};
```
