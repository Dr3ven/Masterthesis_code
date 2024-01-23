using CairoMakie
using Infiltrator
using SodShockTube

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

function sod_collocated()
    # Simulation parameters
    γ = 1.4  # Ratio of specific heats for air
    N = 101   # Number of nodes
    Δx = 0.01
    Δt = 0.0001
    t_final = 0.14

    x = Array(range(0.0, stop=1.0, length=N))
    # Initial conditions for the Sod Shock Tube
    fluxes = zeros(eltype(γ), N + 1, 3)
    ρ = ones(eltype(γ), N)
    u = zeros(eltype(ρ), N)
    p = ones(eltype(ρ), N)
    m = zeros(eltype(ρ), N)
    e = ones(eltype(ρ), N) .* 2.5
    E = zeros(eltype(ρ), N)

    ρ[div(N, 2):end] .= 0.125
    p[div(N, 2):end] .= 0.1
    e[div(N, 2):end] .= 2.0

    # Main time-stepping loop
    t = 0.0

    fig = Figure(size = (800, 600))
    ax1 = Axis(fig[1,1], title="Density", xlabel="x", ylabel="ρ")
    ax2 = Axis(fig[1,2], title="Velocity", xlabel="x", ylabel="u")
    ax3 = Axis(fig[2,1], title="Pressure", xlabel="x", ylabel="p")
    ax4 = Axis(fig[2,2], title="Energy", xlabel="x", ylabel="E")

    scatter!(ax1, x, ρ)
    scatter!(ax2, x, m ./ ρ)
    scatter!(ax3, x, p)
    scatter!(ax4, x, E)
    display(fig)
    counter = 0
    while t < t_final
        counter += 1
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

        if mod(counter, 20) == 0.0
            # Display the final results
            fig = Figure(size = (800, 600))
            ax1 = Axis(fig[1,1], title="Density, t=$t", xlabel="x", ylabel="ρ")
            ax2 = Axis(fig[1,2], title="Velocity", xlabel="x", ylabel="u")
            ax3 = Axis(fig[2,1], title="Pressure", xlabel="x", ylabel="p")
            ax4 = Axis(fig[2,2], title="Energy", xlabel="x", ylabel="E")

            scatter!(ax1, x, ρ)
            scatter!(ax2, x, m ./ ρ)
            scatter!(ax3, x, p)
            scatter!(ax4, x, e)
            #save("./sod_shock_tube.png", fig)
            display(fig)
        end
    end
end