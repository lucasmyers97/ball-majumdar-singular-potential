#ifndef SINGULAR_POTENTIAL_QUADRATURE_HPP
#define SINGULAR_POTENTIAL_QUADRATURE_HPP

#ifndef SINGULAR_POTENTIAL_HEADER_ONLY
#define SINGULAR_POTENTIAL_HEADER_ONLY 1
#endif

#if !SINGULAR_POTENTIAL_HEADER_ONLY
    #define LEBEDEV_HEADER_ONLY 0
#endif
#if SINGULAR_POTENTIAL_IMPLEMENTATION
    #define LEBEDEV_IMPLEMENTATION 1
#endif 

#include "preprocessor.hpp"
#include "Q_tensor_shape.hpp"
#include "singular_potential.hpp"

#if SINGULAR_POTENTIAL_HEADER_ONLY || SINGULAR_POTENTIAL_IMPLEMENTATION

#include "singular_potential.inl"

#endif

#endif
