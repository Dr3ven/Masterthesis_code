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

using CairoMakie

function av_x(B)
    A = 0.5 .* (B[2:end,:] .+ B[1:end-1,:])
end

function ac_wave_1D_noadvect()
    # Physics
    Lx = 1.0                           # domain
    K = 1.0e10                             # shear modulus
    ρ0 = 3000.0                          # initial density at all points
    P0 = 1.0e6                          # initial pressure at all points
    Vx0 = 0.0                          # initial velocity in x-direction for all points
    # Gaussian parameters
    A = 1.0e7                          # gaussian maximum amplitude
    σ = Lx * 0.04                            # standard deviation of the initial pressure distribution
    
    # Plotting parameters
    divisor = 3000 

    # Numerics
    nx = 200                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 10000                             # number of time steps

    # Grid definition
    xc = -(Lx - dx) / 2:dx:(Lx - dx) / 2        # grid nodes in x-direction
    xv = - Lx       / 2:dx: Lx       / 2        # grid vertices in x-direction

    # Allocations
    ρ = ones(nx) .* ρ0
    ρ_old = zeros(nx)
    P = ones(nx) #.* P0
    Vx = zeros(nx + 1)
    Vxdρdx = zeros(nx - 2)
    ρdVxdt = zeros(nx - 2)
    Vxdρdt = zeros(nx - 1)
    ρdPdx = zeros(nx - 1)
    ρVxdVxdx = zeros(nx - 1)
    VxdρVxdx = zeros(nx - 1)
    e = zeros(nx)

    # Initial conditions
    dP = A .* exp.(.- 1.0 .* (xc ./ σ).^2.0)       # initial pressure distribution
    P .= P0 .+ dP
    c = sqrt(K / ρ0)                # speed of sound
    dρ = dP ./ c^2.0                                 # initial density distribution
    ρ .+= dρ

    dt = 1.0e-8 #dx / (c * 4.0)                      # time step size
    t = 0.0                                         # initial time

    xc_vec = Array(xc)
    xv_vec = Array(xv)

    linplots = []

    # Initial plotting
    fig = Figure(size=(1000, 800), fontsize=20)
    #ax1 = Axis(fig[1,1], title="Pressure", ylabel="Pressure", xlabel="Domain",)# limits=(nothing, nothing, P0-(A*3/5), P0+A))
    #ax2 = Axis(fig[3,1], title="Velocity", ylabel="Velocity", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
    ax1 = Axis(fig[1,1], title="Density" , ylabel="Density", xticklabelsvisible=false, xticksvisible=false)
    ax2 = Axis(fig[1,2], title="Velocity", ylabel="Velocity", xticklabelsvisible=false, xticksvisible=false)
    ax3 = Axis(fig[2,1], title="Pressure", xlabel="Domain", ylabel="Pressure")
    ax4 = Axis(fig[2,2], title="Energy"  , xlabel="Domain", ylabel="Energy")
    l0 = lines!(ax1, xc_vec, ρ)
    push!(linplots, l0)
    lines!(ax2, xv_vec, Vx)
    lines!(ax3, xc_vec, P)
    lines!(ax4, xc_vec, e)
    #save("../Plots/Navier-Stokes_acoustic_wave/with_advection_more_realistic_params/$(0).png", fig)
    #display(fig)

    for i = 1:nt
        ρ_old .= ρ
        #c = sqrt(K / maximum(ρ))                # speed of sound
        #@show c

        #Vxdρdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(ρ), dims=1) ./ dx
        ρdVxdt .= .-ρ[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        dρ = ρdVxdt .* dt
        ρ[2:end-1] .= ρ[2:end-1] .+ dρ

        dP = c.^2.0 .* dρ
        P[2:end-1] .= P[2:end-1] .+ dP

        Vxdρdt .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* (av_x(ρ) .- av_x(ρ_old))
        ρdPdx .= .-(1.0 ./ av_x(ρ)) .* diff(P, dims=1) ./ dx
        ρVxdVxdx .= .-(1.0 ./ av_x(ρ)) .* (av_x(ρ) .* Vx[2:end-1]) .* diff(av_x(Vx), dims=1) ./ dx
        #VxdρVxdx .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* diff(ρ .* av_x(Vx), dims=1) ./ dx
        Vx[2:end-1] .= Vx[2:end-1] .+ ρdPdx .* dt .+ ρVxdVxdx .* dt .+ Vxdρdt .* dt

        t += dt
        if i % divisor == 0
            #fig2 = Figure(size=(1000, 800))
            #ax1 = Axis(fig2[1,1], title="Pressure, time = $t")#, limits=(nothing, nothing, -0.25, 0.25))
            #ax2 = Axis(fig2[2,1], title="Velocity")#, limits=(nothing, nothing, -0.25, 0.25))
            #li = lines!(ax1, xc_vec, P, label="time = $t")
            #push!(linplots, li)
            # ax1 = Axis(fig2[1,1], title="Density, time = $t")
            # ax2 = Axis(fig2[1,2], title="Velocity")
            # ax3 = Axis(fig2[2,1], title="Pressure")#, limits=(nothing, nothing, P0, P_max))
            # ax4 = Axis(fig2[2,2], title="Energy")
            l0 = lines!(ax1, xc_vec, ρ)
            push!(linplots, l0)
            lines!(ax2, xv_vec, Vx)
            lines!(ax3, xc_vec, P)
            lines!(ax4, xc_vec, e)
            #lines!(ax2, xv_vec[2:end-1], Vx[2:end-1])
            #save("../Plots/Navier-Stokes_acoustic_wave/with_advection_more_realistic_params/$(i).png", fig2)
            #display(fig2)
        end
    end
    Legend(fig[3,:], linplots, string.(0:dt*divisor:dt*nt), "Total time", nbanks=Int(floor((nt/divisor)+1)), tellhight = false, tellwidth = false)
    rowsize!(fig.layout, 3, 40)
    #save("/home/nils/Masterthesis_code/Plots/Nils_noncons_Euler-equations_acoustic_wave/no_advection/2/4_in_one_acoustic.png", fig)
    #display(fig)
    return ρ, Vx, P, e
end


function ac_wave1D()
    # Physics
    Lx = 1.0                           # domain
    K = 1.0e10                             # shear modulus
    ρ0 = 3000.0                          # initial density at all points
    P0 = 1.0e6                          # initial pressure at all points
    Vx0 = 0.0                          # initial velocity in x-direction for all points
    # Gaussian parameters
    A = 1.0e7                          # gaussian maximum amplitude
    σ = Lx * 0.04                            # standard deviation of the initial pressure distribution
    
    # Plotting parameters
    divisor = 3000 

    # Numerics
    nx = 200                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 10000                             # number of time steps

    # Grid definition
    xc = -(Lx - dx) / 2:dx:(Lx - dx) / 2        # grid nodes in x-direction
    xv = - Lx       / 2:dx: Lx       / 2        # grid vertices in x-direction

    # Allocations
    ρ = ones(nx) .* ρ0
    ρ_old = zeros(nx)
    P = ones(nx) #.* P0
    Vx = zeros(nx + 1)
    Vxdρdx = zeros(nx - 2)
    ρdVxdt = zeros(nx - 2)
    Vxdρdt = zeros(nx - 1)
    ρdPdx = zeros(nx - 1)
    ρVxdVxdx = zeros(nx - 1)
    VxdρVxdx = zeros(nx - 1)
    e = zeros(nx)

    # Initial conditions
    dP = A .* exp.(.- 1.0 .* (xc ./ σ).^2.0)       # initial pressure distribution
    P .= P0 .+ dP
    c = sqrt(K / ρ0)                # speed of sound
    dρ = dP ./ c^2.0                                 # initial density distribution
    ρ .+= dρ

    dt = 1.0e-8 #dx / (c * 4.0)                      # time step size
    t = 0.0                                         # initial time

    xc_vec = Array(xc)
    xv_vec = Array(xv)

    linplots = []

    # Initial plotting
    fig = Figure(size=(1000, 800), fontsize=20)
    #ax1 = Axis(fig[1,1], title="Pressure", ylabel="Pressure", xlabel="Domain",)# limits=(nothing, nothing, P0-(A*3/5), P0+A))
    #ax2 = Axis(fig[3,1], title="Velocity", ylabel="Velocity", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
    ax1 = Axis(fig[1,1], title="Density" , ylabel="Density", xticklabelsvisible=false, xticksvisible=false)
    ax2 = Axis(fig[1,2], title="Velocity", ylabel="Velocity", xticklabelsvisible=false, xticksvisible=false)
    ax3 = Axis(fig[2,1], title="Pressure", xlabel="Domain", ylabel="Pressure")
    ax4 = Axis(fig[2,2], title="Energy"  , xlabel="Domain", ylabel="Energy")
    l0 = lines!(ax1, xc_vec, ρ)
    push!(linplots, l0)
    lines!(ax2, xv_vec, Vx)
    lines!(ax3, xc_vec, P)
    lines!(ax4, xc_vec, e)
    #save("../Plots/Navier-Stokes_acoustic_wave/with_advection_more_realistic_params/$(0).png", fig)
    #display(fig)

    for i = 1:nt
        ρ_old .= ρ
        #c = sqrt(K / maximum(ρ))                # speed of sound
        #@show c

        Vxdρdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(ρ), dims=1) ./ dx
        ρdVxdt .= .-ρ[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        dρ = Vxdρdx .* dt .+ ρdVxdt .* dt
        ρ[2:end-1] .= ρ[2:end-1] .+ dρ

        dP = c.^2.0 .* dρ
        P[2:end-1] .= P[2:end-1] .+ dP

        Vxdρdt .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* (av_x(ρ) .- av_x(ρ_old))
        ρdPdx .= .-(1.0 ./ av_x(ρ)) .* diff(P, dims=1) ./ dx
        ρVxdVxdx .= .-(1.0 ./ av_x(ρ)) .* (av_x(ρ) .* Vx[2:end-1]) .* diff(av_x(Vx), dims=1) ./ dx
        VxdρVxdx .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* diff(ρ .* av_x(Vx), dims=1) ./ dx
        Vx[2:end-1] .= Vx[2:end-1] .+ ρdPdx .* dt .+ ρVxdVxdx .* dt .+ VxdρVxdx .* dt .+ Vxdρdt .* dt

        t += dt
        if i % divisor == 0
            #fig2 = Figure(size=(1000, 800))
            #ax1 = Axis(fig2[1,1], title="Pressure, time = $t")#, limits=(nothing, nothing, -0.25, 0.25))
            #ax2 = Axis(fig2[2,1], title="Velocity")#, limits=(nothing, nothing, -0.25, 0.25))
            #li = lines!(ax1, xc_vec, P, label="time = $t")
            #push!(linplots, li)
            # ax1 = Axis(fig2[1,1], title="Density, time = $t")
            # ax2 = Axis(fig2[1,2], title="Velocity")
            # ax3 = Axis(fig2[2,1], title="Pressure")#, limits=(nothing, nothing, P0, P_max))
            # ax4 = Axis(fig2[2,2], title="Energy")
            l0 = lines!(ax1, xc_vec, ρ)
            push!(linplots, l0)
            lines!(ax2, xv_vec, Vx)
            lines!(ax3, xc_vec, P)
            lines!(ax4, xc_vec, e)
            #lines!(ax2, xv_vec[2:end-1], Vx[2:end-1])
            #save("../Plots/Navier-Stokes_acoustic_wave/with_advection_more_realistic_params/$(i).png", fig2)
            #display(fig2)
        end
    end
    Legend(fig[3,:], linplots, string.(0:dt*divisor:dt*nt), "Total time", nbanks=Int(floor((nt/divisor)+1)), tellhight = false, tellwidth = false)
    rowsize!(fig.layout, 3, 40)
    #save("/home/nils/Masterthesis_code/Plots/Nils_noncons_Euler-equations_acoustic_wave/with_realistic_parameters/4_in_one_acoustic.png", fig)
    #display(fig)
    return ρ, Vx, P, e, xc, xv
end

ρ, Vx, P, e, xc, xv = ac_wave1D()
ρ_noadvect, Vx_noadvect, P_noadvect, e_noadvect = ac_wave_1D_noadvect()

error_ρ = abs.(ρ .- ρ_noadvect)
error_Vx = abs.(Vx .- Vx_noadvect)
error_P = abs.(P .- P_noadvect)
error_e = abs.(e .- e_noadvect)

errorsum_ρ = sum(error_ρ)
errorsum_Vx= sum(error_Vx)
errorsum_P = sum(error_P)
errorsum_e = sum(error_e)

# fig = Figure(size=(1000, 800), fontsize=20)
# ax1 = Axis(fig[1,1], title="Density" , ylabel="Error for density", xticklabelsvisible=false, xticksvisible=false)
# ax2 = Axis(fig[1,2], title="Velocity", ylabel="Error for velocity", xticklabelsvisible=false, xticksvisible=false)
# ax3 = Axis(fig[2,1], title="Pressure", xlabel="Grid points", ylabel="Error for pressure")
# ax4 = Axis(fig[2,2], title="Energy"  , xlabel="Grid points", ylabel="Error for energy")
# lines!(ax1, xc, error_ρ)
# scatter!(ax1, xc, error_ρ)
# lines!(ax2, xv, error_Vx)
# scatter!(ax2, xv, error_Vx)
# lines!(ax3, xc, error_P)
# scatter!(ax3, xc, error_P)
# lines!(ax4, xc, error_e)
# scatter!(ax4, xc, error_e)
# display(fig)

###--------------------------- Different style to plot the data ---------------------------###

function replace_inf(x)
    for i in eachindex(x)
        if x[i] == Inf || x[i] == -Inf
            x[i] = 0.0
        end
    end
    return x
end


logerror_ρ = log10.(error_ρ)
logerror_Vx = log10.(error_Vx[1:end-1])
logerror_P = log10.(error_P)
logerror_e = log10.(error_e)

logerror_ρ = replace_inf(logerror_ρ)
logerror_Vx = replace_inf(logerror_Vx)
logerror_P = replace_inf(logerror_P)
logerror_e = replace_inf(logerror_e)

# f = Figure(size=(1000, 800), fontsize=20)
# ax = Axis(f[1,1], title="Logarithmic absolute difference between advection and no advection", xlabel="Domain", ylabel="log10(error)", yticks=(-24:2:6, string.(-24:2:6)))
# scatterlines!(ax, xc, logerror_ρ, label="Density")
# scatterlines!(ax, xc, logerror_Vx, label="Velocity")
# scatterlines!(ax, xc, logerror_P, label="Pressure")
# lines!(ax, xc, logerror_e, label="Energy", linestyle=:dash, color=:purple)
# axislegend(;position=:cb)
# #save("/home/nils/Masterthesis_code/Plots/Nils_noncons_Euler-equations_acoustic_wave/difference_advection_vs_noadvection.png", f)
# display(f)

###--------------------------- Normalized difference to both models ---------------------------###

norm_error_ρ = error_ρ ./ ρ
norm_error_Vx = error_Vx ./ Vx
norm_error_P = error_P ./ P
norm_error_e = error_e ./ e

norm_error_ρ_noadv = error_ρ ./ ρ_noadvect
norm_error_Vx_noadv = error_Vx ./ Vx_noadvect
norm_error_P_noadv = error_P ./ P_noadvect
norm_error_e_noadv = error_e ./ e_noadvect

f = Figure(size=(1000, 800), fontsize=20)
ax = Axis(f[1,1], title="Relative difference between advection and no advection compared with advection", xlabel="Domain", ylabel="Relative error")
scatterlines!(ax, xc, norm_error_ρ, label="Density")
scatterlines!(ax, xc, norm_error_Vx[1:end-1], label="Velocity")
scatterlines!(ax, xc, norm_error_P, label="Pressure")
lines!(ax, xc, norm_error_e, label="Energy", linestyle=:dash, color=:purple)
axislegend(;position=:rt)
#save("/home/nils/Masterthesis_code/Plots/Nils_noncons_Euler-equations_acoustic_wave/difference_advection_vs_noadvection.png", f)
#display(f)

f = Figure(size=(1000, 800), fontsize=20)
ax = Axis(f[1,1], title="Relative difference between advection and no advection compared with no advection", xlabel="Domain", ylabel="Relative error")
scatterlines!(ax, xc, norm_error_ρ_noadv, label="Density")
scatterlines!(ax, xc, norm_error_Vx_noadv[1:end-1], label="Velocity")
scatterlines!(ax, xc, norm_error_P_noadv, label="Pressure")
lines!(ax, xc, norm_error_e_noadv, label="Energy", linestyle=:dash, color=:purple)
axislegend(;position=:rt)
#save("/home/nils/Masterthesis_code/Plots/Nils_noncons_Euler-equations_acoustic_wave/difference_advection_vs_noadvection.png", f)
#display(f)

###--------------------------- Difference of both relative differences ---------------------------###

diff_norms_ρ = abs.(norm_error_ρ .- norm_error_ρ_noadv)
diff_norms_Vx = abs.(norm_error_Vx .- norm_error_Vx_noadv)
diff_norms_P = abs.(norm_error_P .- norm_error_P_noadv)
diff_norms_e = abs.(norm_error_e .- norm_error_e_noadv)

f = Figure(size=(1000, 800), fontsize=20)
ax = Axis(f[1,1], title="(((|adv - no adv|) / adv) - ((|adv - no adv|) / no adv))", xlabel="Domain", ylabel="Difference")
scatterlines!(ax, xc, diff_norms_ρ, label="Density")
scatterlines!(ax, xc, diff_norms_Vx[2:end], label="Velocity")
scatterlines!(ax, xc, diff_norms_P, label="Pressure")
lines!(ax, xc, diff_norms_e, label="Energy", linestyle=:dash, color=:purple)
axislegend(;position=:rt)
save("/home/nils/Masterthesis_code/Plots/Nils_noncons_Euler-equations_acoustic_wave/difference_of_relative_differences_advection_vs_noadvection.png", f)
display(f)