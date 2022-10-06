# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin


@doc raw"""
    LinearVectorAdvectionEquation2D

The linear Vector advection equations
```math
\partial_t \vec{u} + a_1 \partial_1 \vec{u} + a_2 \partial_2 \vec{u} = 0
```
in two space dimensions with constant velocity `a`.
"""
struct LinearVectorAdvectionEquation2D{NVARS, RealT<:Real} <: AbstractLinearVectorAdvectionEquation{2, NVARS}
  advection_velocity::SVector{2, RealT}
  function LinearVectorAdvectionEquation2D{NVARS, RealT}(a) where {NVARS, RealT <: Real}
    return new(a)
  end
end

function LinearVectorAdvectionEquation2D(a::NTuple{2,RealT}, nvars::Integer) where {RealT <: Real}
  LinearVectorAdvectionEquation2D{nvars,RealT}(SVector(a[1], a[2]))
end

function LinearVectorAdvectionEquation2D(a1::RealT, a2::RealT, nvars::Integer) where {RealT <: Real}
  LinearVectorAdvectionEquation2D{nvars,RealT}(SVector(a1, a2))
end

function varnames(::typeof(cons2cons), equations::LinearVectorAdvectionEquation2D)
  vars = ntuple(n -> "u_" * string(n), Val(ncomponents(equations)))
  return vars
end

function varnames(::typeof(cons2prim), equations::LinearVectorAdvectionEquation2D)
  vars = ntuple(n -> "u_" * string(n), Val(ncomponents(equations)))
  return vars
end

# Calculates translated coordinates `x` for a periodic domain
# function x_trans_periodic_2d(x, domain_length = SVector(10, 10), center = SVector(0, 0))
#   x_normalized = x .- center
#   x_shifted = x_normalized .% domain_length
#   x_offset = ((x_shifted .< -0.5*domain_length) - (x_shifted .> 0.5*domain_length)) .* domain_length
#   return center + x_shifted + x_offset
# end

# Set initial conditions at physical location `x` for time `t`
# """
#     initial_condition_constant(x, t, equations::LinearVectorAdvectionEquation2D)

# A constant initial condition to test free-stream preservation.
# """
# function initial_condition_constant(x, t, equations::LinearVectorAdvectionEquation2D)
#   # Store translated coordinate for easy use of exact solution
#   x_trans = x_trans_periodic_2d(x - equations.advection_velocity * t)

#   return SVector{ncomponents(equations)}(fill(2.0, ncomponents(equations)))
# end


# """
#     initial_condition_convergence_test(x, t, equations::LinearVectorAdvectionEquation2D)

# A smooth initial condition used for convergence tests.
# """
# function initial_condition_convergence_test(x, t, equations::LinearVectorAdvectionEquation2D)
#   # Store translated coordinate for easy use of exact solution
#   x_trans = x - equations.advection_velocity * t

#   c = 1.0
#   A = 0.5
#   L = 2
#   f = 1/L
#   omega = 2 * pi * f
#   Vector = c + A * sin(omega * sum(x_trans))
#   return SVector{ncomponents(equations)}(fill(Vector, ncomponents(equations)))
# end


# """
#     initial_condition_gauss(x, t, equations::LinearVectorAdvectionEquation2D)

# A Gaussian pulse used together with
# [`BoundaryConditionDirichlet(initial_condition_gauss)`](@ref).
# """
# function initial_condition_gauss(x, t, equations::LinearVectorAdvectionEquation2D)
#   # Store translated coordinate for easy use of exact solution
#   x_trans = x_trans_periodic_2d(x - equations.advection_velocity * t)

#   Vector = exp(-(x_trans[1]^2 + x_trans[2]^2))
#   return SVector{ncomponents(equations)}(fill(Vector, ncomponents(equations)))
# end


# """
#     initial_condition_sin_sin(x, t, equations::LinearVectorAdvectionEquation2D)

# A sine wave in the conserved variable.
# """
# function initial_condition_sin_sin(x, t, equations::LinearVectorAdvectionEquation2D)
#   # Store translated coordinate for easy use of exact solution
#   x_trans = x - equations.advection_velocity * t

#   Vector = sinpi(2 * x_trans[1]) * sinpi(2 * x_trans[2])
#   return SVector{ncomponents(equations)}(fill(Vector, ncomponents(equations)))
# end


# """
#     initial_condition_linear_x_y(x, t, equations::LinearVectorAdvectionEquation2D)

# A linear function of `x[1] + x[2]` used together with
# [`boundary_condition_linear_x_y`](@ref).
# """
# function initial_condition_linear_x_y(x, t, equations::LinearVectorAdvectionEquation2D)
#   # Store translated coordinate for easy use of exact solution
#   x_trans = x - equations.advection_velocity * t

#   return SVector{ncomponents(equations)}(fill(sum(x_trans), ncomponents(equations)))
# end

# """
#     boundary_condition_linear_x_y(u_inner, orientation, direction, x, t,
#                                   surface_flux_function,
#                                   equations::LinearVectorAdvectionEquation2D)

# Boundary conditions for
# [`initial_condition_linear_x_y`](@ref).
# """
# function boundary_condition_linear_x_y(u_inner, orientation, direction, x, t,
#                                        surface_flux_function,
#                                        equations::LinearVectorAdvectionEquation2D)
#   u_boundary = initial_condition_linear_x_y(x, t, equations)

#   # Calculate boundary flux
#   if direction in (2, 4) # u_inner is "left" of boundary, u_boundary is "right" of boundary
#     flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
#   else # u_boundary is "left" of boundary, u_inner is "right" of boundary
#     flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
#   end

#   return flux
# end


# """
#     initial_condition_linear_x(x, t, equations::LinearVectorAdvectionEquation2D)

# A linear function of `x[1]` used together with
# [`boundary_condition_linear_x`](@ref).
# """
# function initial_condition_linear_x(x, t, equations::LinearVectorAdvectionEquation2D)
#   # Store translated coordinate for easy use of exact solution
#   x_trans = x - equations.advection_velocity * t

#   return SVector{ncomponents(equations)}(fill(x_trans[1], ncomponents(equations)))
# end

# """
#     boundary_condition_linear_x(u_inner, orientation, direction, x, t,
#                                 surface_flux_function,
#                                 equations::LinearVectorAdvectionEquation2D)

# Boundary conditions for
# [`initial_condition_linear_x`](@ref).
# """
# function boundary_condition_linear_x(u_inner, orientation, direction, x, t,
#                                      surface_flux_function,
#                                      equations::LinearVectorAdvectionEquation2D)
#   u_boundary = initial_condition_linear_x(x, t, equations)

#   # Calculate boundary flux
#   if direction in (2, 4) # u_inner is "left" of boundary, u_boundary is "right" of boundary
#     flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
#   else # u_boundary is "left" of boundary, u_inner is "right" of boundary
#     flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
#   end

#   return flux
# end


# """
#     initial_condition_linear_y(x, t, equations::LinearVectorAdvectionEquation2D)

# A linear function of `x[1]` used together with
# [`boundary_condition_linear_y`](@ref).
# """
# function initial_condition_linear_y(x, t, equations::LinearVectorAdvectionEquation2D)
#   # Store translated coordinate for easy use of exact solution
#   x_trans = x - equations.advection_velocity * t

#   return SVector{ncomponents(equations)}(fill(x_trans[2], ncomponents(equations)))
# end

# """
#     boundary_condition_linear_y(u_inner, orientation, direction, x, t,
#                                 surface_flux_function,
#                                 equations::LinearVectorAdvectionEquation2D)

# Boundary conditions for
# [`initial_condition_linear_y`](@ref).
# """
# function boundary_condition_linear_y(u_inner, orientation, direction, x, t,
#                                      surface_flux_function,
#                                      equations::LinearVectorAdvectionEquation2D)
#   u_boundary = initial_condition_linear_y(x, t, equations)

#   # Calculate boundary flux
#   if direction in (2, 4) # u_inner is "left" of boundary, u_boundary is "right" of boundary
#     flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
#   else # u_boundary is "left" of boundary, u_inner is "right" of boundary
#     flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
#   end

#   return flux
# end


# Pre-defined source terms should be implemented as
# function source_terms_WHATEVER(u, x, t, equations::LinearVectorAdvectionEquation2D)


# Calculate 1D flux for a single point
@inline function flux(u, orientation::Integer, equations::LinearVectorAdvectionEquation2D)
  a = equations.advection_velocity[orientation]
  return a * u
end


# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation
@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer, equations::LinearVectorAdvectionEquation2D)
  λ_max = abs(equations.advection_velocity[orientation])
end


# Calculate 1D flux for a single point in the normal direction
# Note, this directional vector is not normalized
@inline function flux(u, normal_direction::AbstractVector, equations::LinearVectorAdvectionEquation2D)
  a = dot(equations.advection_velocity, normal_direction) # velocity in normal direction
  return a * u
end


# Calculate maximum wave speed in the normal direction for local Lax-Friedrichs-type dissipation
@inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector, equations::LinearVectorAdvectionEquation2D)
  a = dot(equations.advection_velocity, normal_direction) # velocity in normal direction
  return abs(a)
end


@inline have_constant_speed(::LinearVectorAdvectionEquation2D) = Val(true)

@inline function max_abs_speeds(equations::LinearVectorAdvectionEquation2D)
  return abs.(equations.advection_velocity)
end


# Convert conservative variables to primitive
@inline cons2prim(u, equations::LinearVectorAdvectionEquation2D) = u

# Convert conservative variables to entropy variables
@inline cons2entropy(u, equations::LinearVectorAdvectionEquation2D) = u


# Calculate entropy for a conservative state `cons`
@inline entropy(u::Real, ::LinearVectorAdvectionEquation2D) = 0.5 * sum(u.^2)
@inline entropy(u, equations::LinearVectorAdvectionEquation2D) = entropy(sum(u), equations)


# Calculate total energy for a conservative state `cons`
@inline energy_total(u::Real, ::LinearVectorAdvectionEquation2D) = 0.5 * sum(u.^2)
@inline energy_total(u, equations::LinearVectorAdvectionEquation2D) = energy_total(sum(u), equations)


end # @muladd
