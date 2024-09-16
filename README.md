# Ball-Majumdar Singular Potential 

## Basic Info

One method of studying [non-equilibrium](https://en.wikipedia.org/wiki/Non-equilibrium_thermodynamics) [nematic liquid crystal](https://en.wikipedia.org/wiki/Liquid_crystal) systems is by using a [tensorial order parameter](https://arxiv.org/abs/1409.3542) field to model diffusive time evolution.
The standard way of writing down a time evolution equation for such a tensor field, called the [Landau-de Gennes theory](https://en.wikipedia.org/wiki/Landau-de_Gennes_theory), runs into problems when one wishes to distinguish between different [modes of distortion](https://en.wikipedia.org/wiki/Distortion_free_energy_density) in the nematic system.
A different model, proposed by Ball and Majumdar [here](https://doi.org/10.1080/15421401003795555) rectifies some of these problems at the cost of having a time evolution equation which cannot be directly calculated from the tensorial order parameter.
However, it may be numerically calculated and was first computationally implemented [here](https://doi.org/10.1103/physreve.101.032702) by Schimming and Viñals. 

This repository is meant to encapsulate the necessary numerical details in an easily-usable pure C++ implementation of this singular potential theory.
Additionally, it extends the calculation to the purely 2-dimensional case, which is calculable without numerical integration, and in the future will implement a more efficient process which reduces the number of operations by 100 times. 
The hope is that this will increase the adoption of this theory to computational scientists working on mesoscale nematic systems, because C++ code is easily portable to whichever system or language one prefers.
Indeed, a python port is in the works using the wonderful [pybind11](https://pybind11.readthedocs.io/en/stable/index.html).

## Basic Theory

The idea is to write down a free energy in the usual way, use a mean-field approximation to define an energy density, and then locally maximize the entropy subject to a constraint on the probability distribution function associated to molecular direction.
The result is a free energy given by:
```math
F_\text{bulk}[\mathbf Q]
=
E[\mathbf Q]
- T \Delta \mathcal S[\mathbf Q]
```
with bulk energy $E$ given by:
```math
E[\mathbf Q]
=
\int_{\Omega} -\kappa \, \text{tr}\left[\mathbf Q^2\right]
```
where $\kappa$ is a positive constant specifying the alignment strength of the nematic and $\Omega$ is the domain.
The entropy is given by:
```math
\Delta \mathcal S[\mathbf Q]
=
-n k_B \int_{\Omega}
\left[
    \ln 4 \pi
    - \ln Z[\mathbf Q]
    + \boldsymbol \Lambda[\mathbf Q] : \left( \mathbf Q + \frac13 \mathbf I \right)
\right]
```
with $n$ the number density of molecules and $k_B$ Boltzmann's constant.
The Lagrange multiplier $\boldsymbol \Lambda[\mathbf Q]$ is given by:
```math
\mathbf Q
=
\frac{\partial \ln Z}{\partial \boldsymbol \Lambda}
- \frac13 \mathbf I
```
with the partition function $Z[\boldsymbol \Lambda]$ given by:
```math
Z[\boldsymbol \Lambda]
=
\int_{S^2}
\exp\left(\mathbf p^T \boldsymbol \Lambda \mathbf p\right) d \sigma
```
.
$S^2$ is the 2-sphere and $d\sigma$ is the area measure, with $\mathbf p$ the integration variable.
We remark that in two dimensions, these become:
```math
\mathbf Q
=
\frac{\partial \ln Z}{\partial \boldsymbol \Lambda}
- \frac12 \mathbf I
```
and
```math
Z[\boldsymbol \Lambda]
=
\int_{S^1}
\exp\left(\mathbf p^T \boldsymbol \Lambda \mathbf p\right) d \sigma
```
respectively, with the tensors and vectors changing dimension appropriately.
With these relations, one may numerically solve for $\boldsymbol \Lambda$ in terms of $\mathbf Q$.
See the [theory section](theory/README.md) for more details.

## Installation and usage

To use, one constructs a `SingularPotential` object which is templated by a `NematicDimension` enum.
`NematicDimension` can take three values, `full_2D`, `quasi_2D`, and `full_3D`.
`full_3D` corresponds to the well-known $3\times 3$ $Q$-tensor with 5 degrees of freedom.
`full_2D` corresponds to the $2\times 2$ $Q$-tensor with 2 degrees of freedom.
`quasi_2D` is like `full_3D` except the $x$-$z$ and $y$-$z$ components are taken to be zero, giving 3 degrees of freedom.
Physically this is taken to model a thin-film nematic system in which it is energetically unfavorable for the nematic to point out of the plane, but for which in-plane biaxiality is still allowed.

A brief set-up then looks like:
``` cpp
#include "ball_majumdar_singular_potential.hpp"

#include <Eigen/Dense>
#include <iostream>

namespace bmsp = ball_majumdar_singular_potential;

constexpr auto dim = bmsp::NematicDimension::full_3D;

using vec = Eigen::Vector<double, bmsp::QTensorShape<dim>::n_degrees_of_freedom>;

int main()
{
    const unsigned int lebedev_order = 590;
    const double damping_parameter = 1.0;
    const double tolerance = 1e-9;
    const unsigned int maximum_iterations = 100;

    bmsp::SingularPotential<dim> singular_potential(lebedev_order,
                                                    damping_parameter,
                                                    tolerance,
                                                    maximum_iterations);

    vec Q = {1.0 / 3.0, -1.0 / 6.0, 0.0, 0.0, 0.0};

    singular_potential.invert_Q(Q);
    auto Lambda = singular_potential.return_Lambda();

    std::cout << Lambda << std::endl;

    return 0;
}
```

## Installation

### Header only

1. Clone the repository from GitHub
```
git clone --recurse-submodules https://github.com/lucasmyers97/ball-majumdar-singular-potential.git
```
2. Add `ball-majumdar-singular-potential/src` to your include paths
3. `#include <ball_majumdar_singular_potential.hpp>`

### Compiled library

If you include the headers from this library in a large number of compilation units it may slow down compilation time. To bypass this, you must compile the library to a shared or static library within your project.

1. Write a global header in your project which includes `ball_majumdar_singular_potential`, and set the `SINGULAR_POTENTIAL_HEADER_ONLY` macro to `0`.
```cpp
// singular_potential_implementation.hpp
// this is a global header in your project which you will include in many 
// compilation units

#ifndef SINGULAR_POTENTIAL_IMPLEMENTATION_HPP
#define SINGULAR_POTENTIAL_IMPLEMENTATION_HPP

#define SINGULAR_POTENTIAL_HEADER_ONLY 0
#include "ball_majumdar_singular_potential.hpp"

#endif
```

2. Define `#SINGULAR_POTENTIAL_IMPLEMENTATION` before including  `singular_potential_implementation.hpp`
```cpp
// singular_potential_implementation.cpp

#define SINGULAR_POTENTIAL_IMPLEMENTATION 1
#include "singular_potential_implementation.hpp"
```

If you are using CMake for your project, you can create an internal library that you can link to:
```cmake
add_library(singular_potential_implementation
    singular_potential_implementation.cpp)
target_include_directories(singular_potential_implementation
    PRIVATE
    path/to/ball_majumdar_singular_potential/src)
target_include_directories(singular_potential_implementation
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR})
```
