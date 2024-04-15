using CairoMakie
using SodShockTube

average(x) = (x[2:end] + x[1:end-1])/2
extend_vertices(x) = [x[1]; x; x[end]];

function flux_upwind(u, G, dx)
    G_p = (G[2:end-1] .- G[1:end-2])./dx
    G_m = (G[3:end] .- G[2:end-1])./dx
    Gflux = sign.(u[2:end-1]).*(G_p.*(u[2:end-1] .> 0) .+ G_m.*(u[2:end-1] .< 0))

    return Gflux
end

function flux_upwind_center(u, G, dx)
    G_p = (G[2:end-1] .- G[1:end-2])./dx
    G_m = (G[3:end] .- G[2:end-1])./dx
    u_c = average(u)
    Gflux = sign.(u_c[2:end-1]).*(G_p.*(u_c[2:end-1] .> 0) .+ G_m.*(u_c[2:end-1] .< 0))

    return Gflux
end

function ac_wave_cons(; wave=true)
    # Parameters
    L = 1.0                         # Length of the domain
    Nx = 1000                       # Number of spatial grid points
    Nt = 500                        # Number of time steps
    dx = L / (Nx - 1)               # Spatial grid spacing
    dt = 1e-8                       # Temporal grid spacing
    σ = L * 0.04                    # With of the Gaussian
    A = 1.0e7                       # Maximum Gaussian amplitude
    divisor = 125                   # Plotting divisor
    ρ0 = 3000.0                     # Initial density
    P0 = 1.0e6                      # Initial pressure
    K  = 1.0e10                     # Bulk modulus
    γ  = 1.4                        # Adiabatic index
    c = sqrt(K / ρ0)                # speed of sound

    # Allocations
    x   = 0:dx:L
    x_c = average(x)
    u   = zeros(Nx)
    ρ   = ones(Nx-1) .* ρ0
    ρ_v = ones(Nx) .* ρ0
    P   = ones(Nx-1) 
    E   = zeros(Nx-1)
    e   = zeros(Nx-1)

    # Initial conditions
    dP = A .* exp.(.-1.0 .* ((x_c .- 0.5*L) ./ σ).^2.0) 
    P .= P0 .+ dP
    dρ = dP ./ c.^2.0 
    ρ .= ρ0 .+ dρ 
     
    # Setup conservative variables
    m = average(ρ).*u[2:end-1]
    m = [m[1]; m; m[end]]
    t = 0

    linplots = []

    # Initial plotting
    fig = Figure(size=(1000,800), fontsize=20)
    ax1 = Axis(fig[1,1], title="Density", ylabel="Density", xticklabelsvisible=false, xticksvisible=false)
    ax2 = Axis(fig[1,2], title="Velocity", ylabel="Velocity", xticklabelsvisible=false, xticksvisible=false)
    ax3 = Axis(fig[2,1], title="Pressure", xlabel="Domain", ylabel="Pressure")
    ax4 = Axis(fig[2,2], title="Energy", xlabel="Domain", ylabel="Energy")
    l0 = lines!(ax1, x_c, ρ)
    push!(linplots, l0)
    lines!(ax2, x_c, average(u))
    lines!(ax3, x_c, P)
    lines!(ax4, x_c, e)
    display(fig)

    # Solver
    for n in 1:Nt

        # Conservation of mass
        ρ_v = extend_vertices(average(ρ))
        dρ = dt .* flux_upwind_center(u, ρ.*average(u), dx)
        ρ[2:end-1]      .-=  dρ

        # Conservation of momentum
        m[2:end-1]      .-=   dt .* flux_upwind(u, (m.^2)./ρ_v, dx) .+ (dt./dx) .* diff(P, dims=1)

        # Conservation of energy 
        E[2:end-1]      .-=  dt .* flux_upwind_center(u, average(u).*(E .+ P), dx)
        
        # Calculate velocity
        ρ_v = extend_vertices(average(ρ))
        u   = m./ρ_v
        u_c = average(u)

        # Equation of state for pressure difference and the pressure update
        dP = dρ .* c.^2.0;
        P[2:end-1] .= P[2:end-1] .+ dP

        # Internal energy calculation
        e .= P./(γ-1)./ρ

        t += dt

        if n % divisor == 0
            # fig2 = Figure(size=(1000,800), fontsize=20)
            # ax1 = Axis(fig2[1,1], title="Density", ylabel="Density", xticklabelsvisible=false, xticksvisible=false)
            # ax2 = Axis(fig2[1,2], title="Velocity", ylabel="Velocity", xticklabelsvisible=false, xticksvisible=false)
            # ax3 = Axis(fig2[2,1], title="Pressure", xlabel="Domain", ylabel="Pressure")
            # ax4 = Axis(fig2[2,2], title="Energy", xlabel="Domain", ylabel="Energy")
            # opts = (;linewidth = 2, color = :red)
            li = lines!(ax1, x_c, ρ)
            push!(linplots, li)
            lines!(ax2, x_c, u_c)
            lines!(ax3, x_c, P)
            lines!(ax4, x_c, e)
            display(fig)
        end
    end
    Legend(fig[3,:], linplots, string.(0:dt*divisor:dt*Nt), "Total time", nbanks=Int(floor((Nt/divisor)+1)), tellhight = false, tellwidth = false)
    rowsize!(fig.layout, 3, 50)
    #save("Pictures/4_in_one_cons_acousticwave_upwind_withenergy.png", fig)
    display(fig)
end
