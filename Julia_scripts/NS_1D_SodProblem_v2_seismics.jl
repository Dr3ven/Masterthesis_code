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
    Gflux = sign.(u[2:end-1]).*(G_p.*(u[2:end-1] .> 0) .+ G_m.*(u[2:end-1] .< 0))

    return Gflux
end

# attempt to do upwind in staggered formulation
function flux_upwind_center(u, G, dx)
    G_p = (G[2:end-1] .- G[1:end-2])./dx
    G_m = (G[3:end] .- G[2:end-1])./dx
    u_c = average(u)
    Gflux = sign.(u_c[2:end-1]).*(G_p.*(u_c[2:end-1] .> 0) .+ G_m.*(u_c[2:end-1] .< 0))

    return Gflux
end

function wave_b(; wave=true)
    # Parameters
    L = 1.0             # Length of the domain
    Nx = 200           # Number of spatial grid points
    Nt = 4#14000; #1400           # Number of time steps
    dx = L / (Nx - 1)   # Spatial grid spacing
    dt = 5e-6           # Temporal grid spacing
    σ = L * 0.04
    A = 1.0e7
    divisor = 1

    if wave == false
        # Initial conditions
        ρ0 = 1.0     # Initial density
        γ  = 1.4
    else
        # Initial conditions
        ρ0 = 3000.0     # Initial density
        P0 = 1.0e6      # Initial pressure
        K  = 1.0e10     # Bulk modulus
        μ  = 1e-3
        γ  = 1.4
        c = sqrt(K / ρ0)                # speed of sound
    end
    # Arrays to store results
    x   = 0:dx:L
    x_c = average(x)
    u   = zeros(Nx)           # initial u at rest
    ρ   = ones(Nx-1) .* ρ0
    ρ_v = ones(Nx) .* ρ0             # @ vertices
    P   = ones(Nx-1) 
    E   = zeros(Nx-1)
    e   = zeros(Nx-1)

    if wave == false
        ρ[x_c .> 0.5] .= 0.125
        P[x_c .> 0.5] .= 0.1
        E .= P ./ (γ .- 1.0) .+ 0.5.*average(u).^2
        e .= P./(γ-1)./ρ
    else 
        dP = A .* exp.(.-1.0 .* ((x_c .- 0.5*L) ./ σ).^2.0)        # initial pressure distribution
        P .= P0 .+ dP
        dρ = dP ./ c.^2.0 
        ρ .= ρ0 .+ dρ 
    end
     
    # compute conservative variables
    m = average(ρ).*u[2:end-1]
    m = [m[1]; m; m[end]]
    t = 0

    linplots = []

    # Initial plotting
    fig = Figure(size=(1000,800), fontsize=20)
    ax1 = Axis(fig[1,1], title="Density", ylabel="Density", xticklabelsvisible=false, xticksvisible=false)
    ax2 = Axis(fig[1,2], title="Velocity", ylabel="Velocity", xticklabelsvisible=false, xticksvisible=false)
    ax3 = Axis(fig[2,1], title="Pressure", xlabel="Domain", ylabel="Pressure")#, limits=(nothing, nothing, P0, P_max))
    ax4 = Axis(fig[2,2], title="Energy", xlabel="Domain", ylabel="Energy")
    l0 = lines!(ax1, x_c, ρ)
    push!(linplots, l0)
    lines!(ax2, x_c, average(u))
    lines!(ax3, x_c, P)
    lines!(ax4, x_c, e)
    display(fig)

    # Time-stepping loop
    for n in 1:Nt
        t += dt

        # Update density using the continuity equation
        # Note that this specified at the center of a control volume
        # ∂ρ/∂t + ∂(ρu)/∂x = 0
        ρ_v = extend_vertices(average(ρ))
        dρ = dt .* flux_upwind_center(u, ρ.*average(u), dx)
        ρ[2:end-1]      .-=  dρ

        if wave 
            dP = dρ .* c.^2.0;
            P[2:end-1] .+= dP
        end

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
            e .= P./(γ-1)./ρ
            P .= (γ .- 1.0) .* (E .- 0.5.*ρ.* u_c.^2)
        end

        if mod(n,Nt) == 0
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
            
            
            e_anal = values.p./((γ-1).*values.ρ)      # it seems the e calculation is a bit different in the analytical code
        end

        if n % divisor == 0
            # fig2 = Figure(size=(1000,800), fontsize=20)
            # ax1 = Axis(fig2[1,1], title="Density", ylabel="Density", xticklabelsvisible=false, xticksvisible=false)
            # ax2 = Axis(fig2[1,2], title="Velocity", ylabel="Velocity", xticklabelsvisible=false, xticksvisible=false)
            # ax3 = Axis(fig2[2,1], title="Pressure", xlabel="Domain", ylabel="Pressure")#, limits=(nothing, nothing, P0, P_max))
            # ax4 = Axis(fig2[2,2], title="Energy", xlabel="Domain", ylabel="Energy")
            opts = (;linewidth = 2, color = :red)
            li = lines!(ax1, x_c, ρ)
            push!(linplots, li)
            lines!(ax2, x_c, u_c)
            lines!(ax3, x_c, P)
            lines!(ax4, x_c, e)
            if wave == false && n % Nt == 0
                lines!(ax1, x_c, values.ρ; opts...)
                lines!(ax2, x_c, values.u; opts...)
                lines!(ax3, x_c, values.p; opts...)
                lines!(ax4, x_c, e_anal; opts...)
            end
            display(fig)
        end
    end
    Legend(fig[3,:], linplots, string.(0:dt*divisor:dt*Nt), "Total time", nbanks=Int(floor((Nt/divisor)+1)), tellhight = false, tellwidth = false)
    rowsize!(fig.layout, 3, 50)
    #save("C:\\Users\\Nils\\Desktop\\Masterarbeit_plots\\Boris_cons_Euler_acoustic_wave_upwind\\4_in_one_cons_acousticwave_upwind.png", fig)
    display(fig)
end
#ρ_b, u_b, P_b, E_b = run()