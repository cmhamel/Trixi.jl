
using OrdinaryDiffEq
using Plots
using Trixi

###############################################################################
# semidiscretization of the linear advection equation

advection_velocity = (0.2, -0.7)
# advection_velocity = SVector{2, Float64}([0.2, -0.7])

equations = LinearVectorAdvectionEquation2D(advection_velocity, 2)
# equations = LinearVectorAdvectionEquation2D(0.2, -0.7, 2)

# Create DG solver with polynomial degree = 3 and (local) Lax-Friedrichs/Rusanov flux as surface flux
solver = DGSEM(polydeg=3, surface_flux=flux_lax_friedrichs)

coordinates_min = (-1.0, -1.0) # minimum coordinates (min(x), min(y))
coordinates_max = ( 1.0,  1.0) # maximum coordinates (max(x), max(y))

# Create a uniformly refined mesh with periodic boundaries
mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level=4,
                n_cells_max=30_000) # set maximum capacity of tree data structure

function Trixi.initial_condition_convergence_test(x, t, equations::LinearVectorAdvectionEquation2D)
    # Store translated coordinate for easy use of exact solution
    x_trans = x - equations.advection_velocity * t
    
    c = 1.0
    A = 0.5
    L = 2
    f = 1/L
    omega = 2 * pi * f
    Vector = c + A * sin(omega * sum(x_trans))
    return SVector{ncomponents(equations)}(fill(Vector, ncomponents(equations)))
end

# A semidiscretization collects data structures and functions for the spatial discretization
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_convergence_test, solver)


###############################################################################
# ODE solvers, callbacks etc.

# Create ODE problem with time span from 0.0 to 1.0
ode = semidiscretize(semi, (0.0, 10.0));

# At the beginning of the main loop, the SummaryCallback prints a summary of the simulation setup
# and resets the timers
summary_callback = SummaryCallback()

# The AnalysisCallback allows to analyse the solution in regular intervals and prints the results
analysis_callback = AnalysisCallback(semi, interval=100)

# The SaveSolutionCallback allows to save the solution to a file in regular intervals
save_solution = SaveSolutionCallback(interval=100,
                                     solution_variables=cons2prim)

# The StepsizeCallback handles the re-calculcation of the maximum Δt after each time step
stepsize_callback = StepsizeCallback(cfl=1.6)

# Create a CallbackSet to collect all callbacks such that they can be passed to the ODE solver
callbacks = CallbackSet(summary_callback, analysis_callback, save_solution, stepsize_callback)


###############################################################################
# run the simulation

# OrdinaryDiffEq's `solve` method evolves the solution in time and executes the passed callbacks
sol = solve(ode, CarpenterKennedy2N54(williamson_condition=false),
            dt=1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
            save_everystep=false, callback=callbacks);

# Print the timer summary
summary_callback()
