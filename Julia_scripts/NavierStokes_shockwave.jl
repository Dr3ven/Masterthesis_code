using CairoMakie
using Infiltrator


function av_x(B)
    A = 0.5 .* (B[2:end] .+ B[1:end-1])
end

function wave_1D()
    # Physics
    Lx = 10.0                           # domain
    P0 = 0.0                          # initial pressure at all points
    γ  = 1.4                            # adiabatic index

    # Numerics
    nx = 100                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 100000                             # number of time steps

    # Grid definition
    xc = -(Lx - dx) / 2:dx:(Lx - dx) / 2        # grid nodes in x-direction
    xv = - Lx       / 2:dx: Lx       / 2        # grid vertices in x-direction

    # Allocations
    E = zeros(nx)
    c = ones(nx)
    ρ = ones(nx)
    Mx = zeros(nx + 1)
    P = ones(nx) 
    Vx = zeros(nx + 1)
    dEudx = zeros(nx)
    dPudx = zeros(nx)
    dEdt = zeros(nx)
    dVxdx = zeros(nx - 1)
    dρdt = zeros(nx)
    dMxdt = zeros(nx - 1)

    # Initial conditions
    P .= P0 .+ exp.(.- 1.0 .* xc.^2.0)        # initial pressure distribution
    E .= P ./ (γ .- 1.0) + (1.0 ./ 2.0) .* av_x(Vx).^2

    P_max = maximum(P)
    t = 0.0                              # initial time
    dt = 0.0001

    xc_vec = Vector(xc)
    xv_vec = Vector(xv)

    # Initial plotting
    fig = Figure()
    ax = Axis(fig[1,1], title="t = $t")
    lines!(ax, xc_vec, P)
    display(fig)

    for i = 1:nt
        #c .= sqrt.(P ./ ρ)                            # speed of sound

        dρdt[2:end-1] .= -(ρ[2:end-1] .* diff(Vx[2:end-1], dims=1)) ./ dx
        ρ[2:end-1] .= ρ[2:end-1] .+ dρdt[2:end-1] .* dt
        
        dVxdx[2:end-1] .= diff(ρ[2:end-1] .* (av_x(Vx[2:end-1]) .> 0.0) .* Vx[2:end-2].^2.0, dims=1) ./ dx .+ diff(ρ[2:end-1] .* (av_x(Vx[2:end-1]) .< 0.0) .* Vx[3:end-1].^2.0, dims=1) ./ dx
        dMxdt[2:end-1] .= -dVxdx[2:end-1] .- diff(P[2:end-1], dims=1) ./ dx
        Mx[2:end-1] .= Mx[2:end-1] .+ dMxdt .* dt

        dEudx[2:end-1] .= diff((Vx[2:end-1] .> 0.0) .* Vx[2:end-1] .* av_x(E), dims=1) ./ dx .+ diff((Vx[2:end-1] .< 0.0) .* Vx[2:end-1] .* av_x(E), dims=1) ./ dx
        dPudx[2:end-1] .= diff((Vx[2:end-1] .> 0.0) .* Vx[2:end-1] .* av_x(P), dims=1) ./ dx .+ diff((Vx[2:end-1] .< 0.0) .* Vx[2:end-1] .* av_x(P), dims=1) ./ dx
        dEdt[2:end-1] .= -dEudx[2:end-1] .- dPudx[2:end-1]
        E[2:end-1] .= E[2:end-1] .+ dEdt[2:end-1] .* dt
        
        Vx[2:end-1] .= Mx[2:end-1] ./ av_x(ρ)
        P .= (av_x(Vx) .> 0.0) .* (γ .- 1.0) .* (E .- 0.5 .* Mx[1:end-1] .* Vx[1:end-1]) .+ (av_x(Vx) .< 0.0) .* (γ .- 1.0) .* (E .- 0.5 .* Mx[2:end] .* Vx[2:end]) #c.^2.0 .* ρ

        #Vx[1] = Vx[2]
        #Vx[end] = Vx[end-1]

        t += dt
        if i % 100 == 0
            fig2 = Figure()
            ax1 = Axis(fig2[1,1], title="Pressure, time = $t", limits=(nothing, nothing, P0, P_max))
            ax2 = Axis(fig2[1,2], title="Velocity")
            ax3 = Axis(fig2[2,1], title="Energy")
            ax4 = Axis(fig2[2,2], title="Density")
            #ylims!(ax, -1.2, 1.2)
            lines!(ax1, xc_vec, P)
            lines!(ax2, xv_vec, Vx)
            lines!(ax3, xc_vec, E)
            lines!(ax4, xc_vec, ρ)
            display(fig2)
        end
    end
end