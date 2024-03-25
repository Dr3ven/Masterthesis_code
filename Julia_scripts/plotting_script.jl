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
    G_p = (G[2:end-1] - G[1:end-2])/dx
    G_m = (G[2:end-1] + G[3:end])/dx
    Gflux = sign.(u[2:end-1]).*(G_p.*(u[2:end-1] .> 0) + G_m.*(u[2:end-1] .< 0))

    return Gflux
end

# attempt to do upwind in staggered formulation
function flux_upwind_center(u, G, dx)
    
    G_p = (G[2:end-1] - G[1:end-2])/dx
    G_m = (G[2:end-1] + G[3:end])/dx
    u_c = average(u)
    Gflux = sign.(u_c[2:end-1]).*(G_p.*(u_c[2:end-1] .> 0) + G_m.*(u_c[2:end-1] .< 0))

    return Gflux
end

function sod()
    # Parameters
    Lx = 1.0 #1000.0             # Length of the domain
    Nx = 100           # Number of spatial grid points
    Nt = 1400; #1400           # Number of time steps
    dx = Lx / (Nx - 1)   # Spatial grid spacing
    dt = 1e-4           # Temporal grid spacing
    divisor = 280
   
    # Initial conditions
    ρ0 = 1.0    # Initial density
    μ  = 1e-3
    γ  = 1.4
    c = 1.0 # velocity of sound
    
    # Arrays to store results
    x   = 0:dx:Lx
    x_c = average(x)
    u   = zeros(Nx)           # initial u at rest
    ρ   = ones(Nx-1)
    ρ_v = ones(Nx)            # @ vertices
    P   = ones(Nx-1) #.* 1.0e6
    e   = zeros(Nx-1)
    ρ[x_c .> Lx / 2.0] .= 0.125
    P[x_c .> Lx / 2.0] .= 0.1 #1.0e5
     
    # compute conservative variables
    m = average(ρ).*u[2:end-1]
    m = [m[1]; m; m[end]]
    E = P./((γ - 1.0)) + 0.5.*average(u).^2
    t = 0

    linplots = []

    #c = sqrt.(abs.(γ*P./ρ))     # shouldnt go negative but sometimes does
    #dt_courant = dx/(maximum(abs.(u_c) + c))
    #@show dt_courant, t

    problem = ShockTubeProblem(
                    geometry = (0.0, Lx, Lx / 2.0), # left edge, right edge, initial shock location
                    left_state = (ρ = 1.0, u = 0.0, p = 1.0),#1.0e6),
                    right_state = (ρ = 0.125, u = 0.0, p = 0.1),# 1.0e5),
                    t = 0.14, γ = γ)
    positions, regions, vals = solve(problem, x_c);    
    e .= P./(γ-1)./ρ
    e_anal = vals.p./((γ-1).*vals.ρ)

    opts = (;linewidth = 2, color = :red, )  
    fig = Figure(size=(1000,800); fontsize=26)
    Label(fig[1,1:2], "Time = 0.14", tellwidth=false, font=:bold ,fontsize=30)
    ax1 = Axis(fig[2,1], title="Density", ylabel="Density", xticklabelsvisible=false, xticksvisible=false)
    ax2 = Axis(fig[2,2], title="Velocity", ylabel="Velocity", xticklabelsvisible=false, xticksvisible=false)
    ax3 = Axis(fig[3,1], title="Pressure", xlabel="Domain", ylabel="Pressure")#, limits=(nothing, nothing, P0, P_max))
    ax4 = Axis(fig[3,2], title="Energy", xlabel="Domain", ylabel="Energy")
    lines!(ax1, x_c, vals.ρ; opts...)
    lines!(ax2, x_c, vals.u; opts...)
    lines!(ax3, x_c, vals.p; opts...)
    lines!(ax4, x_c, e_anal; opts...)
    # l0 = lines!(ax1, x_c, ρ)
    # lines!(ax2, x_c, average(u))
    # lines!(ax3, x_c, P)
    # lines!(ax4, x_c, e)
    # push!(linplots, l0)
    #vs1a = vspan!(ax1, 0.33, 0.49, color=(:lightskyblue1, 0.3))
    #vs1b = vspan!(ax1, [0.33, 0.61, 0.7325], [0.49, 0.64, 0.7625], color=[(:lightskyblue1, 0.3), (:orange, 0.3), (:purple, 0.3)])
    vs1 = vspan!(ax1, [0.33, 0.61, 0.7325], [0.49, 0.64, 0.7625], color=[(:lightskyblue1, 0.3), (:orange, 0.3), (:purple, 0.3)])
    vs2 = vspan!(ax2, [0.33, 0.61, 0.7325], [0.49, 0.64, 0.7625], color=[(:lightskyblue1, 0.3), (:orange, 0.3), (:purple, 0.3)])
    vs3 = vspan!(ax3, [0.33, 0.61, 0.7325], [0.49, 0.64, 0.7625], color=[(:lightskyblue1, 0.3), (:orange, 0.3), (:purple, 0.3)])
    vs4 = vspan!(ax4, [0.33, 0.61, 0.7325], [0.49, 0.64, 0.7625], color=[(:lightskyblue1, 0.3), (:orange, 0.3), (:purple, 0.3)])
    
    #myblue = colorant"rgba(124,240,255,0.3)"
    #myorange = colorant"rgba(191,0,255,0.3)"
    #mypurple = colorant"rgba(224,166,108,0.3)"
    #box1 = PolyElement(color=RGBAf(224,166,108,0.3))
    #box2 = PolyElement(color=:orange, opacity=0.3)
    #box3 = PolyElement(color=:purple, opacity=0.3)
    #coloredboxes = [box1, box2, box3]
    #Legend(fig[3,:], coloredboxes, ["Rarefaction", "Contact", "Shock"], "Wave types", tellwidth = false, nbanks=3)#, tellhight = false, tellwidth = false)#, orientation=:horizontal, tellhight = false, tellwidth = false)
    #rowsize!(fig.layout, 3, 60)
    
    #save("/home/nils/Masterthesis_code/Plots/Boris_sod_shock_code/Analytical_solution_colored.png", fig)
    display(fig)

    # Time-stepping loop
    #=for n in 1:Nt
        t += dt

        # Update density using the continuity equation
        # Note that this specified at the center of a control volume
        # ∂ρ/∂t + ∂(ρu)/∂x = 0
        ρ_v = extend_vertices(average(ρ))
        ρ[2:end-1]      -=   dt * flux_upwind_center(u, ρ.*average(u), dx)

        # update ρu (momentum) using the momentum equation:
        # This is formulated around the vertices of the control volume  

        # ∂(m)/∂t + ∂(u*m + P)/∂x = 0
        m[2:end-1]      -=   dt * flux_upwind(u, (m.^2)./ρ_v, dx) + dt/dx.*(P[2:end] - P[1:end-1])

        # update energy 
        # ∂e/∂t + ∂/∂x( (ρu/ρ)(E + P) ) = 0
        E[2:end-1]      -=  dt* flux_upwind_center(u, average(u).*(E + P), dx)

        # 
        ρ_v = extend_vertices(average(ρ))
        u   = m./ρ_v
        u_c = average(u)
        
        P = (γ - 1.0) .* (E - 0.5*ρ.* u_c.^2)
        e .= P./(γ-1)./ρ

        #P = ρ.*c^2.0;
        #P = ρ

        #=if n % 280 == 0
            c = sqrt.(abs.(γ*P./ρ))     # shouldnt go negative but sometimes does
            dt_courant = dx/(maximum(abs.(u_c) + c))
            #@show dt_courant, t

            problem = ShockTubeProblem(
                            geometry = (0.0, Lx, Lx / 2.0), # left edge, right edge, initial shock location
                            left_state = (ρ = 1.0, u = 0.0, p = 1.0),#1.0e6),
                            right_state = (ρ = 0.125, u = 0.0, p = 0.1),# 1.0e5),
                            t = t, γ = γ)
            positions, regions, vals = solve(problem, x_c);
            #e = (E .- 0.5*ρ.*u.^2)
            
            
            e_anal = vals.p./((γ-1).*vals.ρ)      # it seems the e calculation is a bit different in the analytical code
        end=#

        if n % divisor == 0
            #fig = Figure(size=(1000,800); fontsize=20)
            #ax1 = Axis(fig[1,1], title="Density, time = $t")
            #ax2 = Axis(fig[1,2], title="Velocity")
            #ax3 = Axis(fig[2,1], title="Pressure")#, limits=(nothing, nothing, P0, P_max))
            #ax4 = Axis(fig[2,2], title="Energy")
            opts = (;linewidth = 2, color = :red)
            li = lines!(ax1, x_c, ρ)
            lines!(ax2, x_c, average(u))
            lines!(ax3, x_c, P)
            lines!(ax4, x_c, e)
            push!(linplots, li)
            #lines!(ax1, x_c, vals.ρ; opts...)
            #lines!(ax2, x_c, vals.u; opts...)
            #lines!(ax3, x_c, vals.p; opts...)
            #lines!(ax4, x_c, e_anal; opts...)
            display(fig)
        end
    end
    Legend(fig[3,:], linplots, string.(round.(0:dt*divisor:dt*Nt, digits=5)), "Total time", tellwidth = false, nbanks=Int(floor((Nt/divisor)+1)))#, tellhight = false, tellwidth = false)#, orientation=:horizontal, tellhight = false, tellwidth = false)
    rowsize!(fig.layout, 3, 60)
    #save("/home/nils/Masterthesis_code/Plots/Boris_sod_shock_code/Time_evolution_sod_shock_staggered.png", fig)
    display(fig)
    #return ρ, u, P, E=#
end

#ρ_b, u_b, P_b, E_b = sod()