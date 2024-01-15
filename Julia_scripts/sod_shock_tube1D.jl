using CairoMakie
using Infiltrator
using SodShockTube

function sod_shock_tube(ρ, u, p, γ)
    # Compute derived quantities
    e = p / ((γ - 1)* ρ) + 0.5 * u^2  # Specific internal energy
    #c = sqrt(γ * p / ρ)  # Speed of sound
    
    # Conservative form of the equations
    ρ_flux = ρ * u
    momentum_flux = ρ * u^2 + p
    energy_flux = u * (ρ * e + p)

    return [ρ_flux, momentum_flux, energy_flux]
end

function sod_shock_tube_upwind!(ρ, u, p, Δt, Δx, γ)
    N = length(ρ)

    # Create arrays to store fluxes at cell interfaces
    fluxes = zeros(N-1, 3)

    # Compute fluxes using upwind scheme
    for i in 1:N-1
        if u[i] >= 0
            fluxes[i, :] = sod_shock_tube(ρ[i], u[i], p[i], γ)
        else
            fluxes[i, :] = sod_shock_tube(ρ[i+1], u[i+1], p[i+1], γ)
        end
    end

    # Update solution
    for i in 2:N-1
        ρ[i-1] -= Δt / Δx * (fluxes[i-1, 1] - fluxes[i, 1])
        u[i-1] -= Δt / Δx * (fluxes[i-1, 2] - fluxes[i, 2])
        p[i-1] -= Δt / Δx * (fluxes[i-1, 3] - fluxes[i, 3])
    end

    # Reflective boundary condition at the right end
    ρ[end] = ρ[end-1]
    u[end] = -u[end-1]
    p[end] = p[end-1]
end

function sod_shock_tube_maccormack!(ρ, u, p, Δt, Δx, γ)
    N = length(ρ)

    # Predictor step
    ρ_star = copy(ρ)
    u_star = copy(u)
    p_star = copy(p)

    for i in 2:N-1
        ρ_star[i] = ρ[i] - Δt/Δx * (ρ[i]*u[i+1] - ρ[i-1]*u[i])
        u_star[i] = u[i] - Δt/Δx * (ρ[i]*u[i+1]*u[i+1] + p[i+1] - ρ[i-1]*u[i]*u[i] - p[i]) / (ρ[i+1] - ρ[i-1])
        p_star[i] = p[i] - Δt/Δx * (u[i+1]*(ρ[i+1]*u[i+1] + p[i+1]) - u[i]*(ρ[i-1]*u[i] + p[i-1])) / (ρ[i+1] - ρ[i-1])
    end

    # Corrector step
    for i in 2:N-1
        ρ[i] = 0.5 * (ρ[i] + ρ_star[i] - Δt/Δx * (ρ_star[i]*u_star[i] - ρ_star[i-1]*u_star[i-1]))
        @infiltrate
        u[i] = 0.5 * (u[i] + u_star[i] - Δt/Δx * (ρ_star[i]*u_star[i]*u_star[i] + p_star[i] - ρ_star[i-1]*u_star[i-1]*u_star[i-1] - p_star[i-1]) / (ρ_star[i] - ρ_star[i-1]))
        p[i] = 0.5 * (p[i] + p_star[i] - Δt/Δx * (u_star[i]*(ρ_star[i]*u_star[i] + p_star[i]) - u_star[i-1]*(ρ_star[i-1]*u_star[i-1] + p_star[i-1])) / (ρ_star[i] - ρ_star[i-1]))
        if any(isnan, ρ)
            println("NaN found in ρ at iteration $i")
            @infiltrate
        end
    end

    # Reflective boundary condition at the right end
    ρ[end] = ρ[end-1]
    u[end] = -u[end-1]
    p[end] = p[end-1]
end

function sod_shock_tube_laxwendroff!(ρ, u, p, Δt, Δx, γ)
    N = length(ρ)

    # Predictor step
    ρ_star = copy(ρ)
    u_star = copy(u)
    p_star = copy(p)

    for i in 2:N-1
        ρ_star[i] = ρ[i] - 0.5 * Δt/Δx * (ρ[i+1]*u[i+1] - ρ[i-1]*u[i-1])
        u_star[i] = u[i] - 0.5 * Δt/Δx * ((ρ[i+1]*u[i+1]^2 + p[i+1]) - (ρ[i-1]*u[i-1]^2 + p[i-1]))
        p_star[i] = p[i] - 0.5 * Δt/Δx * (u[i+1]*(ρ[i+1]*u[i+1] + p[i+1]) - u[i-1]*(ρ[i-1]*u[i-1] + p[i-1]))
    end

    # Corrector step
    for i in 2:N-1
        ρ[i] = ρ[i] - Δt/Δx * (ρ_star[i]*u_star[i] - ρ_star[i-1]*u_star[i-1])
        u[i] = u[i] - Δt/Δx * ((ρ_star[i]*u_star[i]^2 + p_star[i]) - (ρ_star[i-1]*u_star[i-1]^2 + p_star[i-1]))
        p[i] = p[i] - Δt/Δx * (u_star[i]*(ρ_star[i]*u_star[i] + p_star[i]) - u_star[i-1]*(ρ_star[i-1]*u_star[i-1] + p_star[i-1]))
    end

    # Reflective boundary condition at the right end
    ρ[end] = ρ[end-1]
    u[end] = -u[end-1]
    p[end] = p[end-1]
end

function run()
    # Simulation parameters
    γ = 1.4  # Ratio of specific heats for air
    N = 501   # Number of nodes
    Δx = 0.01
    Δt = 0.0001
    t_final = 0.2

    x = Array(range(0.0, stop=1.0, length=N))
    # Initial conditions for the Sod Shock Tube
    ρ = ones(N)
    u = zeros(N)
    p = ones(N)
    ρ[div(N, 2):end] .= 0.125
    p[div(N, 2):end] .= 0.1

    # Main time-stepping loop
    t = 0.0
    while t < t_final
        sod_shock_tube_upwind!(ρ, u, p, Δt, Δx, γ)
        #sod_shock_tube_maccormack!(ρ, u, p, Δt, Δx, γ)
        #sod_shock_tube_laxwendroff!(ρ, u, p, Δt, Δx, γ)
        t += Δt
    end

    # Display the final results
    println("Density: ", ρ)
    println("Velocity: ", u)
    println("Pressure: ", p)

    fig = Figure(size = (800, 800))
    ax1 = Axis(fig[1,1], title="Density", xlabel="x", ylabel="ρ")
    ax2 = Axis(fig[1,2], title="Velocity", xlabel="x", ylabel="u")
    ax3 = Axis(fig[2,1], title="Pressure", xlabel="x", ylabel="p")

    lines!(ax1, x, ρ)
    lines!(ax2, x, u)
    lines!(ax3, x, p)
    #save("./sod_shock_tube.png", fig)
    display(fig)
end

###-------------------------------------------------------------------------###

function shock_tube_benchmark()
    problem = ShockTubeProblem(
            geometry = (0.0, 1.0, 0.5),
            left_state = (ρ = 1.0, u = 0.0, p = 1.0),
            right_state = (ρ = 0.125, u = 0.0, p = 0.1),
            t = 0.2,
            γ = 1.4
        );

    xs = LinRange(0.0, 1.0, 500);
    positions, regions, values = solve(problem, xs)

    f = Figure(size = (800, 800))
    ax_ρ = Axis(f[1,1], xlabel = "x", ylabel = "ρ", title = "Density")
    ax_u = Axis(f[2,1], xlabel = "x", ylabel = "u", title = "Velocity")
    ax_p = Axis(f[1,2], xlabel = "x", ylabel = "p", title = "Pressure")
    ax_E = Axis(f[2,2], xlabel = "x", ylabel = "E", title = "Stagnation Energy")

    opts = (;linewidth = 4)

    lines!(ax_ρ, values.x, values.ρ; opts...)
    lines!(ax_u, values.x, values.u; opts...)
    lines!(ax_p, values.x, values.p; opts...)
    lines!(ax_E, values.x, values.e; opts...)

    display(f)
    save("./sod_shock_tube_package.png", f)
end