using CairoMakie, SodShockTube

average(x) = (x[2:end] + x[1:end-1])/2
extend_vertices(x) = [x[1]; x; x[end]];

function flux_upwind(u, G, dx)
    G_p = (G[2:end-1] - G[1:end-2])/dx
    G_m = (G[2:end-1] + G[3:end])/dx
    Gflux = sign.(u[2:end-1]).*(G_p.*(u[2:end-1] .> 0) + G_m.*(u[2:end-1] .< 0))

    return Gflux
end

function flux_upwind_center(u, G, dx)
    G_p = (G[2:end-1] - G[1:end-2])/dx
    G_m = (G[2:end-1] + G[3:end])/dx
    u_c = average(u)
    Gflux = sign.(u_c[2:end-1]).*(G_p.*(u_c[2:end-1] .> 0) + G_m.*(u_c[2:end-1] .< 0))

    return Gflux
end

function sod()
    # Parameters
    Lx = 1.0                # Length of the domain
    Nx = 200                # Number of spatial grid points
    Nt = 1400               # Number of time steps
    dx = Lx / (Nx - 1)      # Spatial grid spacing
    dt = 1e-4               # Temporal grid spacing
    divisor = 350           # Plotting divisor
    γ  = 1.4                # Adiabatic index
    
    # Allocations
    x   = 0:dx:Lx
    x_c = average(x)
    u   = zeros(Nx)
    ρ   = ones(Nx-1)
    ρ_v = ones(Nx)
    P   = ones(Nx-1)
    e   = zeros(Nx-1)
    ρ[x_c .> Lx / 2.0] .= 0.125
    P[x_c .> Lx / 2.0] .= 0.1
    E = P./((γ - 1.0)) + 0.5.* ρ .* average(u).^2
     
    # Setup conservative variables
    m = average(ρ).*u[2:end-1]
    m = [m[1]; m; m[end]]
    t = 0

    linplots = []

    # Analytical solution calculated by the SodShockTube.jl package
    problem = ShockTubeProblem(
                    geometry = (0.0, Lx, Lx / 2.0), # left edge, right edge, initial shock location
                    left_state = (ρ = 1.0, u = 0.0, p = 1.0),#1.0e6),
                    right_state = (ρ = 0.125, u = 0.0, p = 0.1),# 1.0e5),
                    t = 0.14, γ = γ)
    positions, regions, vals = solve(problem, x_c);    
    e .= P./(γ-1)./ρ
    e_anal = vals.p./((γ-1).*vals.ρ)

    # Initial plotting
    fig = Figure(size=(1000,800); fontsize=20)
    ax1 = Axis(fig[1,1], title="Density, time = 0.0", ylabel="Density", xticklabelsvisible=false, xticksvisible=false)
    ax2 = Axis(fig[1,2], title="Velocity", ylabel="Velocity", xticklabelsvisible=false, xticksvisible=false)
    ax3 = Axis(fig[2,1], title="Pressure", ylabel="Pressure", xlabel="Domain")
    ax4 = Axis(fig[2,2], title="Energy", ylabel="Energy", xlabel="Domain")
    #opts = (;linewidth = 2, color = :red, )  
    # lines!(ax1, x_c, vals.ρ; opts...)
    # lines!(ax2, x_c, vals.u; opts...)
    # lines!(ax3, x_c, vals.p; opts...)
    # lines!(ax4, x_c, e_anal; opts...)
    l0 = lines!(ax1, x_c, ρ)
    lines!(ax2, x_c, average(u))
    lines!(ax3, x_c, P)
    lines!(ax4, x_c, e)
    push!(linplots, l0)
    display(fig)

    # Solver
    for n in 1:Nt

        # Conservation of mass
        ρ_v = extend_vertices(average(ρ))
        ρ[2:end-1]      -=   dt * flux_upwind_center(u, ρ.*average(u), dx)
 
        ρ[end] = ρ[end-1]

        # Conservation of momentum
        m[2:end-1]      -=   dt * flux_upwind(u, (m.^2)./ρ_v, dx) + dt/dx.*(P[2:end] - P[1:end-1])

        # Conservation of energy
        E[2:end-1]      -=  dt* flux_upwind_center(u, average(u).*(E + P), dx)

        # Update velocity
        ρ_v = extend_vertices(average(ρ))
        u   = m./ρ_v
        u_c = average(u)
        
        # Equation of state for pressure of an isentropic gas
        P = (γ - 1.0) .* (E - 0.5*ρ.* u_c.^2)
        
        # Internal energy calculation
        e .= P./(γ-1)./ρ

        t += dt

        # Plotting
        if n % divisor == 0
            # fig = Figure(size=(1000,800); fontsize=20)
            # ax1 = Axis(fig[1,1], title="Density, time = 0.0", ylabel="Density", xticklabelsvisible=false, xticksvisible=false)
            # ax2 = Axis(fig[1,2], title="Velocity", ylabel="Velocity", xticklabelsvisible=false, xticksvisible=false)
            # ax3 = Axis(fig[2,1], title="Pressure", ylabel="Pressure", xlabel="Domain")
            # ax4 = Axis(fig[2,2], title="Energy", ylabel="Energy", xlabel="Domain")
            # opts = (;linewidth = 2, color = :red)
            li = lines!(ax1, x_c, ρ)
            lines!(ax2, x_c, average(u))
            lines!(ax3, x_c, P)
            lines!(ax4, x_c, e)
            push!(linplots, li)
            # lines!(ax1, x_c, vals.ρ; opts...)
            # lines!(ax2, x_c, vals.u; opts...)
            # lines!(ax3, x_c, vals.p; opts...)
            # lines!(ax4, x_c, e_anal; opts...)

            display(fig)
        end
    end
    #Legend(fig[3,:], linplots, string.(round.(0:dt*divisor:dt*Nt, digits=5)), "Total time", tellwidth = false, nbanks=Int(floor((Nt/divisor)+1)))#, tellhight = false, tellwidth = false)#, orientation=:horizontal, tellhight = false, tellwidth = false)
    #rowsize!(fig.layout, 3, 60)
    #save("/home/nils/Masterthesis_code/Plots/Boris_sod_shock_code/Time_evolution_sod_shock_staggered.png", fig)
    #display(fig)
end
