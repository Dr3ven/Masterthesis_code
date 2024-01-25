using CairoMakie
using Infiltrator
using SodShockTube

function av_x(B)
    A = 0.5 .* (B[2:end] .+ B[1:end-1])
end

function check_upwind(u)
    if u >= 0.0
        return -1
    elseif u < 0.0 
        return 1
    end
end

function check_sign(u)
    if u >= 0.0
        return 1
    elseif u < 0.0 
        return -1
    end
end

function sod_shock_tube(ρ, u, p, γ)
    # Compute derived quantities
    m = ρ * u   
    e = p / ((γ - 1.0) * ρ)
    E = ρ * e + 0.5 * m * u

    # Conservative form of the equations
    ρ_flux = m
    momentum_flux = m^2 / ρ
    energy_flux = (m / ρ) * (E + p)

    return [ρ_flux, momentum_flux, energy_flux], m, e, E ,p
end

function sod_staggered()
    # Simulation parameters
    γ = 1.4                             # Ratio of specific heats for air
    N = 1000                             # Number of nodes
    L = 1.0                             # Length of the domain              
    x = LinRange(0.0, L, N)             # Spatial discretization
    Δx = x[2] - x[1]                    # Spatial step  
    t = 0.0                             # Initial time    
    Δt = 0.0001                         # Time step      
    t_final = 0.14                      # Final time

    # Initial conditions for the Sod Shock Tube
    fluxes = zeros(eltype(γ), N + 1, 3)
    c = zeros(eltype(γ), N)
    ρ = ones(eltype(γ), N)
    u = zeros(eltype(ρ), N + 1)
    p = ones(eltype(ρ), N)
    m = zeros(eltype(ρ), N + 1)
    e = ones(eltype(ρ), N)
    e_ana = zeros(eltype(ρ), N)
    E = zeros(eltype(ρ), N)

    ρ[div(N, 2):end] .= 0.125
    p[div(N, 2):end] .= 0.1
    #e .= p ./ ((γ .- 1.0) .* ρ)

    # Main time-stepping loop
    t = 0.0

    fig = Figure(size = (800, 600))
    ax1 = Axis(fig[1,1], title="Density", xlabel="x", ylabel="ρ")
    ax2 = Axis(fig[1,2], title="Velocity", xlabel="x", ylabel="u")
    ax3 = Axis(fig[2,1], title="Pressure", xlabel="x", ylabel="p")
    ax4 = Axis(fig[2,2], title="Energy", xlabel="x", ylabel="e")

    scatter!(ax1, x, ρ)
    scatter!(ax2, x, av_x(m) ./ ρ)
    scatter!(ax3, x, p)
    scatter!(ax4, x, E)
    display(fig)
    counter = 0
    problem = ShockTubeProblem(
                geometry = (0.0, 1.0, 0.5),
                left_state = (ρ = 1.0, u = 0.0, p = 1.0),
                right_state = (ρ = 0.125, u = 0.0, p = 0.1),
                t = t_final,
                γ = 1.4
           );
    positions, regions, values = solve(problem, x)
    ρ_ana = values.ρ
    u_ana = values.u
    p_ana = values.p
    e_ana .= p_ana ./ ((γ .- 1.0) .* ρ_ana)
    
    while t < t_final
        counter += 1
        #c .= sqrt.((1.0 ./ p) .* ρ)
        # σ = maximum(abs.(u) + c) .* (Δt ./ Δx) # Stability criterion σ ≤ 1
        for i in 2:N-1
                # Check sign of velocity and set parameters accordingly
                s_u = check_upwind(u[i])
                sgn = check_sign(u[i])
                
                # Compute fluxes using upwind scheme
                if u[i] >= 0.0
                    fluxes[i, :], m[i] , e[i], E[i], p[i] = sod_shock_tube(ρ[i], u[i], p[i], γ)
                elseif u[i] < 0.0
                    fluxes[i, :], m[i] , e[i], E[i], p[i] = sod_shock_tube(ρ[i+1], u[i+1], p[i+1], γ)
                end

                # Update solution
                ρ[i] = ρ[i] - sgn * (Δt / Δx) * (fluxes[i,1] - fluxes[i+s_u,1])
                m[i] = m[i] - sgn * (Δt / Δx) * (fluxes[i,2] - fluxes[i+s_u,2]) - (Δt / (2 * Δx)) * (p[i+1] - p[i-1])
                E[i] = E[i] - sgn * (Δt / Δx) * (fluxes[i,3] - fluxes[i+s_u,3])
                
                # Calculate pressure and velocity
                p[i] = (γ - 1.0) * (E[i] - 0.5 * m[i] * u[i])
                u[i] = m[i] / ρ[i]
        end

        t += Δt

        if mod(counter, t_final / Δt) == 0.0
            # Display the final results
            fig = Figure(size = (800, 600))
            ax1 = Axis(fig[1,1], title="Density, t=$t", xlabel="x", ylabel="ρ")
            ax2 = Axis(fig[1,2], title="Velocity", xlabel="x", ylabel="u")
            ax3 = Axis(fig[2,1], title="Pressure", xlabel="x", ylabel="p")
            ax4 = Axis(fig[2,2], title="Energy", xlabel="x", ylabel="e")

            opts = (;linewidth = 2, color = :red)
            scatter!(ax1, x, ρ)
            scatter!(ax2, x, av_x(m) ./ ρ)
            scatter!(ax3, x, p)
            scatter!(ax4, x, e)
            lines!(ax1, x, ρ_ana; opts...)
            lines!(ax2, x, u_ana; opts...)
            lines!(ax3, x, p_ana; opts...)
            lines!(ax4, x, e_ana; opts...)
            #save("./sod_shock_tube.png", fig)
            display(fig)
        end
    end
end