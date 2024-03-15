using CairoMakie
using Infiltrator


function av_x(B)
    A = 0.5 .* (B[2:end] .+ B[1:end-1])
end

function flux_upwind(u, G, dx)
    G_p = (G[2:end-1] .- G[1:end-2])./dx
    G_m = (G[2:end-1] .+ G[3:end])./dx
    Gflux = sign.(u[2:end-1]).*(G_p.*(u[2:end-1] .> 0) .+ G_m.*(u[2:end-1] .< 0))

    return Gflux
end

function flux_upwind_center(u, G, dx)
    
    G_p = (G[2:end-1] - G[1:end-2])/dx
    G_m = (G[2:end-1] + G[3:end])/dx
    u_c = av_x(u)
    Gflux = sign.(u_c[2:end-1]).*(G_p.*(u_c[2:end-1] .> 0) + G_m.*(u_c[2:end-1] .< 0))

    return Gflux
end

extend_vertices(x) = [x[1]; x; x[end]];

function wave_1D_b(;wave=false)
    # Physics
    Lx = 1.0                           # domain
    P0 = 0.0                          # initial pressure at all points
    γ  = 1.4                            # adiabatic index
    σ  = 0.1                            # gaussian amplitude
    A  = 10.0                            # gaussian standard deviation  
    divisor = 100

    # Numerics
    nx = 100                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 14000                             # number of time steps

    # Grid definition
    xv = 0:dx:1        # grid vertices in x-direction
    xc = av_x(xv)        # grid nodes in x-direction

    # Allocations
    e = zeros(nx)
    E = zeros(nx)
    ρ = ones(nx)
    ρ_v = ones(nx + 1)
    Mx = zeros(nx + 1)
    P = ones(nx) 
    Vx = zeros(nx + 1)
    dEdt = zeros(nx)
    dρdt = zeros(nx)
    dMxdt = zeros(nx - 1)

    # Initial conditions
    if wave == false
        ρ[div(nx, 2):end] .= 0.125
        P[div(nx, 2):end] .= 0.1
        E .= P ./ (γ .- 1.0) + 0.5 .* av_x(Vx).^2
    else    
        P .= P0 .+ A .* exp.(.-1.0 .* (xc./ σ).^2.0)        # initial pressure distribution
    end

    P_max = maximum(P)
    t = 0.0                              # initial time
    dt = 1.0e-5

    xc_vec = Vector(xc)
    xv_vec = Vector(xv)

    # Initial plotting
    fig = Figure()
    ax = Axis(fig[1,1], title="Pressure at t = $t")
    lines!(ax, xc_vec, P)
    display(fig)

    for i = 1:nt
        ρ_v = extend_vertices(av_x(ρ))

        dρdt[2:end-1] .= flux_upwind_center(Vx, ρ .* av_x(Vx), dx)#-(ρ[2:end-1] .* diff(Vx[2:end-1], dims=1)) ./ dx
        ρ[2:end-1] .= ρ[2:end-1] .- dρdt[2:end-1] .* dt
        
        dMxdt .= flux_upwind(Vx, (Mx.^2.0) ./ ρ_v, dx) .+ (P[2:end] - P[1:end-1]) ./ dx #-dVxdx[2:end-1] .- diff(P[2:end-1], dims=1) ./ dx
        Mx[2:end-1] .= Mx[2:end-1] .- dMxdt .* dt

        dEdt[2:end-1] .= flux_upwind_center(Vx, av_x(Vx) .* (E .+ P), dx) #-dEudx[2:end-1] .- dPudx[2:end-1]
        E[2:end-1] .= E[2:end-1] .- dEdt[2:end-1] .* dt
        
        Vx[2:end-1] .= Mx[2:end-1] ./ av_x(ρ)
        if wave == false
            P .= (γ .- 1.0) .* (E .- 0.5 .* ρ .* av_x(Vx).^2.0)
        else
            P .= c.^2.0 .* ρ
        end
        e = P ./ (γ .- 1) ./ ρ

        t += dt
        if i % divisor == 0
            fig2 = Figure()
            ax1 = Axis(fig2[1,1], title="Density, time = $t")#, limits=(nothing, nothing, P0, P_max))
            ax2 = Axis(fig2[1,2], title="Velocity")
            ax3 = Axis(fig2[2,1], title="Pressure")
            ax4 = Axis(fig2[2,2], title="Energy")
            #ylims!(ax, -1.2, 1.2)
            scatter!(ax1, xc_vec[2:end-1], ρ[2:end-1])
            lines!(ax1, xc_vec[2:end-1], ρ[2:end-1])
            scatter!(ax2, xv_vec[2:end-1], Vx[2:end-1])
            lines!(ax2, xv_vec[2:end-1], Vx[2:end-1])
            scatter!(ax3, xc_vec[2:end-1], P[2:end-1])
            lines!(ax3, xc_vec[2:end-1], P[2:end-1])
            scatter!(ax4, xc_vec[2:end-1], e[2:end-1])
            lines!(ax4, xc_vec[2:end-1], e[2:end-1])
            display(fig2)
        end
    end
end