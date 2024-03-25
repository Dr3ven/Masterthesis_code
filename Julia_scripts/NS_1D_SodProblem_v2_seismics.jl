# One-dimensional Navier-Stokes solver using finite differences
# This code creates the setup for the Sod Tube model problem, which is a typical test case for compressible flow solvers
#
# luckily someone implemented the analytical solution for this problem in julia:
# https://github.com/archermarx/SodShockTube.jl
#
#
# Note: here we use a staggered grid with P & ρ at the center
using CairoMakie, SodShockTube

average(x) = (x[2:end] + x[1:end-1])/2
extend_vertices(x) = [x[1]; x; x[end]];

"""
Compute the flux using a 1D upwind scheme
"""
function flux_upwind(u, G, dx)
    G_p = (G[2:end-1] .- G[1:end-2])./dx
    G_m = (G[3:end] .- G[2:end-1])./dx
    Gflux = #=sign.(u[2:end-1]).*=#(G_p.*(u[2:end-1] .> 0) .+ G_m.*(u[2:end-1] .< 0))

    return Gflux
end

# attempt to do upwind in staggered formulation
function flux_upwind_center(u, G, dx)
    G_p = (G[2:end-1] .- G[1:end-2])./dx
    G_m = (G[3:end] .- G[2:end-1])./dx
    u_c = average(u)
    Gflux = #=sign.(u_c[2:end-1]).*=#(G_p.*(u_c[2:end-1] .> 0) .+ G_m.*(u_c[2:end-1] .< 0))

    return Gflux
end

function wave_b(; wave=true)
    # Parameters
    L = 1.0             # Length of the domain
    Nx = 100           # Number of spatial grid points
    Nt = 100000; #1400           # Number of time steps
    dx = L / (Nx - 1)   # Spatial grid spacing
    dt = 1e-8           # Temporal grid spacing
    σ = L * 0.04
    A = 10.0

    # Initial conditions
    ρ0 = 3000.0     # Initial density
    P0 = 1.0e6      # Initial pressure
    K  = 1.0e10     # Bulk modulus
    μ  = 1e-3
    γ  = 1.4
    c = sqrt(K / ρ0)                # speed of sound
    #c = 5000.0 # velocity of sound
    
    # Arrays to store results
    x   = 0:dx:L
    x_c = average(x)
    u   = zeros(Nx)           # initial u at rest
    ρ   = ones(Nx-1) .* ρ0
    ρ_v = ones(Nx) .* ρ0             # @ vertices
    P   = ones(Nx-1) 

    if wave == false
        ρ[x_c .> 0.5] .= 0.125
        P[x_c .> 0.5] .= 0.1
    else 
        P .= P0 .+ A .*exp.(.-1.0 .* ((x_c .- 0.5*L) ./ σ).^2.0)        # initial pressure distribution
        ρ .= P./c.^2.0  
    end

     
    # compute conservative variables
    m = average(ρ).*u[2:end-1]
    m = [m[1]; m; m[end]]
    E = P ./ (γ .- 1.0) .+ 0.5.*average(u).^2
    t = 0

    # Time-stepping loop
    for n in 1:Nt
        t += dt

        # Update density using the continuity equation
        # Note that this specified at the center of a control volume
        # ∂ρ/∂t + ∂(ρu)/∂x = 0
        ρ_v = extend_vertices(average(ρ))
        ρ[2:end-1]      .-=   dt .* flux_upwind_center(u, ρ.*average(u), dx)

        # update ρu (momentum) using the momentum equation:
        # This is formulated around the vertices of the control volume  
        # ∂(m)/∂t + ∂(u*m + P)/∂x = 0
        #P_v = extend_vertices(average(P))
        m[2:end-1]      .-=   dt .* flux_upwind(u, (m.^2)./ρ_v, dx) .+ (dt./dx) .* diff(P, dims=1)

        # update energy 
        # ∂e/∂t + ∂/∂x( (ρu/ρ)(E + P) ) = 0
        E[2:end-1]      .-=  dt .* flux_upwind_center(u, average(u).*(E .+ P), dx)
        
        ρ_v = extend_vertices(average(ρ))
        u   = m./ρ_v
        u_c = average(u)

        if wave == false
            P .= (γ .- 1.0) .* (E .- 0.5.*ρ.* u_c.^2)
        else
            P .= ρ.*c.^2.0;
        end

        if mod(n,1) == 0
            c_anal = sqrt.(abs.(γ*P./ρ))     # shouldnt go negative but sometimes does
            #dt_courant = dx/(maximum(abs.(u_c) + c_anal))
            #@show dt_courant, t

            problem = ShockTubeProblem(
                            geometry = (0.0, 1.0, 0.5), # left edge, right edge, initial shock location
                            left_state = (ρ = 1.0, u = 0.0, p = 1.0),
                            right_state = (ρ = 0.125, u = 0.0, p = 0.1),
                            t = t, γ = γ)
            positions, regions, values = solve(problem, x_c);
            #e = (E .- 0.5*ρ.*u.^2)
            
            
            e = P./(γ-1)./ρ
            e_anal = values.p./((γ-1).*values.ρ)      # it seems the e calculation is a bit different in the analytical code
        end

        if n % 1000 == 0
            fig2 = Figure()
            ax1 = Axis(fig2[1,1], title="Density, time = $t")
            ax2 = Axis(fig2[1,2], title="Velocity")
            ax3 = Axis(fig2[2,1], title="Pressure")#, limits=(nothing, nothing, P0, P_max))
            ax4 = Axis(fig2[2,2], title="Energy")

            #ylims!(ax, -1.2, 1.2)
            opts = (;linewidth = 2, color = :red)
            lines!(ax1, x_c, ρ)
            lines!(ax2, x_c, average(u))
            lines!(ax3, x_c, P)
            lines!(ax4, x_c, e)
            #lines!(ax1, x_c, values.ρ; opts...)
            #lines!(ax2, x_c, values.u; opts...)
            #lines!(ax3, x_c, values.p; opts...)
            #lines!(ax4, x_c, e_anal; opts...)
            display(fig2)
        end
    end

   return ρ, u, P, E
end

#ρ_b, u_b, P_b, E_b = run()