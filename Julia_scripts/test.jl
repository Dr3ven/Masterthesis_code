# #using CairoMakie
# #=
# function av_x(B)
#     A = 0.5 .* (B[2:end] .+ B[1:end-1])
# end

# # Define parameters
# L = 1.0         # Length of the domain
# Nx = 100        # Number of grid points
# T = 0.5         # Total simulation time
# dt = 0.001      # Time step size
# dx = L / Nx     # Grid spacing
# γ = 1.4         # Ratio of specific heats

# # Initialize arrays
# x_v = range(0, stop=L, length=Nx+1)   # Grid points
# x = av_x(x_v)   # Grid points
# Q = zeros(3, Nx)                    # Conserved variables: [ρ, ρ*u, E]

# # Initial condition: square wave
# ρ_l = 1.0       # Density on the left
# ρ_r = 0.1       # Density on the right
# u_l = 0.0       # Velocity on the left
# u_r = 0.0       # Velocity on the right
# p_l = ρ_l       # Pressure on the left (assuming isentropic relation)
# p_r = ρ_r       # Pressure on the right (assuming isentropic relation)
# E_l = p_l / (γ - 1) + 0.5 * ρ_l * u_l^2  # Total energy on the left
# E_r = p_r / (γ - 1) + 0.5 * ρ_r * u_r^2  # Total energy on the right

# # Set initial conditions inside the square wave
# Q[1, 1:49] .= ρ_l
# Q[2, 1:49] .= ρ_l * u_l
# Q[3, 1:49] .= E_l
# Q[1, 50:end] .= ρ_r
# Q[2, 50:end] .= ρ_r * u_r
# Q[3, 50:end] .= E_r

# # Function to compute flux
# function compute_flux(Q)
#     ρ = Q[1, :]
#     u = Q[2, :] ./ ρ
#     p = (γ - 1) .* (Q[3, :] .- 0.5 .* ρ .* u.^2)
#     return hcat(Q[2, :], Q[2, :] .* u .+ p, (Q[3, :] .+ p) .* u)
# end

# # Godunov method with forward Euler time-stepping
# function godunov_forward_euler(Q, dt, dx, γ, num_steps)
#     for n in 1:num_steps
#         # Compute fluxes at cell interfaces
#         F = compute_flux(Q)
#         @info F

#         # Update solution using forward Euler
#         Q[:, 2:end-1] .-= (dt / dx) * (F[:, 2:end-1] - F[:, 1:end-2])

#         # Apply periodic boundary conditions
#         Q[:, 1] .= Q[:, end-1]
#         Q[:, end] .= Q[:, 2]
#     end
#     return Q
# end

# # Perform simulation
# num_steps = Int(T / dt)
# Q = godunov_forward_euler(Q, dt, dx, γ, num_steps)

# # Plot results
# ρ = Q[1, :]
# u = Q[2, :] ./ ρ
# p = (γ - 1) .* (Q[3, :] .- 0.5 .* ρ .* u.^2)
# fig = Figure(size=(1000,800))
# ax1 = Axis(fig[1,1], title="Density", xlabel="x", ylabel="ρ")
# scatter!(x, ρ, label="Density", xlabel="x", ylabel="ρ", legend=:topright)
# lines!(x, ρ, label="Density", xlabel="x", ylabel="ρ", legend=:topright)
# display(fig)
# =#
# using CairoMakie
# using OrdinaryDiffEq
# using Trixi

# ###############################################################################
# # semidiscretization of the compressible Euler equations
# gamma = 1.4
# equations = CompressibleEulerEquations2D(gamma)

# initial_condition = initial_condition_density_wave

# solver = DGSEM(polydeg = 5, surface_flux = flux_central)

# coordinates_min = (-1.0, -1.0)
# coordinates_max = (1.0, 1.0)
# mesh = TreeMesh(coordinates_min, coordinates_max,
#                 initial_refinement_level = 2,
#                 n_cells_max = 30_000)

# semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

# ###############################################################################
# # ODE solvers, callbacks etc.

# tspan = (0.0, 2.0)
# ode = semidiscretize(semi, tspan)

# summary_callback = SummaryCallback()

# analysis_interval = 100
# analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

# alive_callback = AliveCallback(analysis_interval = analysis_interval)

# save_solution = SaveSolutionCallback(interval = 100,
#                                      save_initial_solution = true,
#                                      save_final_solution = true,
#                                      solution_variables = cons2prim)

# stepsize_callback = StepsizeCallback(cfl = 1.6)

# callbacks = CallbackSet(summary_callback,
#                         analysis_callback, alive_callback,
#                         save_solution,
#                         stepsize_callback)

# ###############################################################################
# # run the simulation

# sol = solve(ode, CarpenterKennedy2N54(williamson_condition = false),
#             dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
#             save_everystep = false, callback = callbacks);
# summary_callback() # print the timer summary

# fig = Figure(size=(1000,800); fontsize=26)
# Label(fig[1,1:2], "Time = 0.14", tellwidth=false, font=:bold ,fontsize=30)
# ax1 = Axis(fig[2,1], title="Density", xticklabelsvisible=false, xticksvisible=false)
# ax2 = Axis(fig[2,2], title="Velocity", xticklabelsvisible=false, xticksvisible=false)
# ax3 = Axis(fig[3,1], title="Pressure")#, limits=(nothing, nothing, P0, P_max))
# ax4 = Axis(fig[3,2], title="Energy")
# lines!(ax1, x_c, sol.u[1]; opts...)
# lines!(ax2, x_c, sol.u[2]; opts...)
# #lines!(ax3, x_c, sol.; opts...)
# #lines!(ax4, x_c, sol.; opts...)

A = 1.0
xc = range(start=10.0e6, stop=10.0e8, length=6)
σ = 100.0
test = A .* 1.0 ./ exp.(.-1.0 .* (xc ./ σ).^2.0)

fig = Figure(size=(1000,800); fontsize=26)
ax1 = Axis(fig[1,1], title="test")
lines!(ax1, xc, test)
display(fig)
