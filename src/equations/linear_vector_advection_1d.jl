# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin


@doc raw"""
    LinearVectorAdvectionEquation1D

The linear Vector advection equation
```math
\partial_t \vec{u} + a \partial_1 \vec{u}  = 0
```
in one space dimension with constant velocity `a`.
"""
struct LinearVectorAdvectionEquation1D{NVARS, RealT<:Real} <: AbstractLinearVectorAdvectionEquation{1, NVARS}
  advection_velocity::SVector{1, RealT}
  function LinearVectorAdvectionEquation1D{NVARS, RealT}(a) where {NVARS, RealT <: Real}
    return new(a)
  end
end

function LinearVectorAdvectionEquation1D(a::SVector{1, RealT}, nvars::Integer) where {RealT <: Real}
  LinearVectorAdvectionEquation1D{nvars,RealT}(a)
end

function LinearVectorAdvectionEquation1D(a::RealT, nvars::Integer) where {RealT <: Real}
  LinearVectorAdvectionEquation1D{nvars,RealT}(SVector(a))
end


function varnames(::typeof(cons2cons), equations::LinearVectorAdvectionEquation1D)
  vars = ntuple(n -> "u_" * string(n), Val(ncomponents(equations)))
  return vars
end

function varnames(::typeof(cons2prim), equations::LinearVectorAdvectionEquation1D)
  vars = ntuple(n -> "u_" * string(n), Val(ncomponents(equations)))
  return vars
end

# Set initial conditions at physical location `x` for time `t`
"""
    initial_condition_constant(x, t, equations::LinearVectorAdvectionEquation1D)

A constant initial condition to test free-stream preservation.
"""
function initial_condition_constant(x, t, equations::LinearVectorAdvectionEquation1D)
  # Store translated coordinate for easy use of exact solution
  x_trans = x - equations.advection_velocity * t

  # return SVector(2.0)
  return SVector{ncomponents(equations)}(fill(2., ncomponents(equations)))
end


"""
    initial_condition_convergence_test(x, t, equations::LinearVectorAdvectionEquation1D)

A smooth initial condition used for convergence tests
(in combination with [`BoundaryConditionDirichlet(initial_condition_convergence_test)`](@ref)
in non-periodic domains).
"""
function initial_condition_convergence_test(x, t, equations::LinearVectorAdvectionEquation1D)
  # Store translated coordinate for easy use of exact solution
  x_trans = x - equations.advection_velocity * t

  c = 1.0
  A = 0.5
  L = 2
  f = 1/L
  omega = 2 * pi * f
  Vector = c + A * sin(omega * sum(x_trans))
  # return SVector(Vector)
  return SVector{ncomponents(equations)}(fill(Vector, ncomponents(equations)))
end


"""
    initial_condition_gauss(x, t, equations::LinearVectorAdvectionEquation1D)

A Gaussian pulse used together with
[`BoundaryConditionDirichlet(initial_condition_gauss)`](@ref).
"""
function initial_condition_gauss(x, t, equations::LinearVectorAdvectionEquation1D)
  # Store translated coordinate for easy use of exact solution
  x_trans = x - equations.advection_velocity * t

  Vector = exp(-(x_trans[1]^2))
  # return SVector(Vector)
  return SVector{ncomponents(equations)}(fill(Vector, ncomponents(equations)))
end


"""
    initial_condition_sin(x, t, equations::LinearVectorAdvectionEquation1D)

A sine wave in the conserved variable.
"""
function initial_condition_sin(x, t, equations::LinearVectorAdvectionEquation1D)
  # Store translated coordinate for easy use of exact solution
  x_trans = x - equations.advection_velocity * t

  Vector = sinpi(2 * x_trans[1])
  # return SVector(Vector)
  return SVector{ncomponents(equations)}(fill(Vector, ncomponents(equations)))
end


"""
    initial_condition_linear_x(x, t, equations::LinearVectorAdvectionEquation1D)

A linear function of `x[1]` used together with
[`boundary_condition_linear_x`](@ref).
"""
function initial_condition_linear_x(x, t, equations::LinearVectorAdvectionEquation1D)
  # Store translated coordinate for easy use of exact solution
  x_trans = x - equations.advection_velocity * t

  # return SVector(x_trans[1])
  return SVector{ncomponents(equations)}(fill(xtrans[1], ncomponents(equations)))
end

"""
    boundary_condition_linear_x(u_inner, orientation, direction, x, t,
                                surface_flux_function,
                                equation::LinearVectorAdvectionEquation1D)

Boundary conditions for
[`initial_condition_linear_x`](@ref).
"""
function boundary_condition_linear_x(u_inner, orientation, direction, x, t,
                                     surface_flux_function,
                                     equations::LinearVectorAdvectionEquation1D)
  u_boundary = initial_condition_linear_x(x, t, equations)

  # Calculate boundary flux
  if direction == 2  # u_inner is "left" of boundary, u_boundary is "right" of boundary
    flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
  else # u_boundary is "left" of boundary, u_inner is "right" of boundary
    flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
  end

  return flux
end


# Pre-defined source terms should be implemented as
# function source_terms_WHATEVER(u, x, t, equations::LinearVectorAdvectionEquation1D)


# Calculate 1D flux in for a single point
@inline function flux(u, orientation::Integer, equations::LinearVectorAdvectionEquation1D)
  a = equations.advection_velocity[orientation]
  return a * u
end


# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
@inline function max_abs_speed_naive(u_ll, u_rr, orientation, equations::LinearVectorAdvectionEquation1D)
  λ_max = abs(equations.advection_velocity[orientation])
end


@inline have_constant_speed(::LinearVectorAdvectionEquation1D) = Val(true)

@inline function max_abs_speeds(equations::LinearVectorAdvectionEquation1D)
  return abs.(equations.advection_velocity)
end


# Convert conservative variables to primitive
@inline cons2prim(u, equations::LinearVectorAdvectionEquation1D) = u

# Convert conservative variables to entropy variables
@inline cons2entropy(u, equations::LinearVectorAdvectionEquation1D) = u


# Calculate entropy for a conservative state `cons`
@inline entropy(u::Real, ::LinearVectorAdvectionEquation1D) = 0.5 * sum(u.^2)
@inline entropy(u, equations::LinearVectorAdvectionEquation1D) = entropy(sum(u), equations)


# Calculate total energy for a conservative state `cons`
@inline energy_total(u::Real, ::LinearVectorAdvectionEquation1D) = 0.5 * sum(u.^2)
@inline energy_total(u, equations::LinearVectorAdvectionEquation1D) = energy_total(sum(u), equations)


end # @muladd
