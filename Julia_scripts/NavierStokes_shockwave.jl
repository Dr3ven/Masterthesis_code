using CairoMakie
using Infiltrator


function av_x(B)
    A = 0.5 .* (B[2:end] .+ B[1:end-1])
end

function wave_1D_v1(;wave=false)
    # Physics
    Lx = 1.0                           # domain
    P0 = 0.0                          # initial pressure at all points
    γ  = 1.4                            # adiabatic index
    σ  = 1.0                            # gaussian amplitude
    A  = 10.0                            # gaussian standard deviations
    c  = 1.0

    divisor = 100

    # Numerics
    nx = 200                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 14000                             # number of time steps

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

    ρ_old = zeros(nx)
    Mx_old = zeros(nx + 1)
    E_old = zeros(nx)

    # Initial conditions
    if wave
        P .= P0 .+ A .* exp.(.-(1.0 ./ σ) .* xc.^2.0)        # initial pressure distribution
        ρ .= P./ c.^2.0
    else
        ρ[div(nx, 2):end] .= 0.125
        P[div(nx, 2):end] .= 0.1
        E .= P ./ (γ .- 1.0) + (1.0 ./ 2.0) .* av_x(Vx).^2
    end

    P_max = maximum(P)
    t = 0.0                              # initial time
    dt = 1.0e-5

    xc_vec = Vector(xc)
    xv_vec = Vector(xv)

    # Initial plotting
    fig = Figure(size=(1000,800))
    ax = Axis(fig[1,1], title="t = $t")
    lines!(ax, xc_vec[2:end-1], P[2:end-1])
    display(fig)

    for i = 1:nt
        ρ_old .= ρ
        Mx_old .= Mx
        E_old .= E
        #c .= sqrt.(P ./ ρ)                            # speed of sound

        upwind_ρ = ((ρ[2:end-2] .* Vx[3:end-2] .* (Vx[3:end-2] .> 0.0)).+ (ρ[3:end-1] .* Vx[3:end-2] .* (Vx[3:end-2] .< 0.0))) ./ dx
        dρdt[3:end-2] .= .-(diff(upwind_ρ, dims=1) ./ dx)
        ρ[2:end-1] .= ρ_old[2:end-1] .+ dρdt[2:end-1] .* dt
        
        # Mx_new = Mx_old + (-∂_x(ρVx^2) - ∂_x(P)) * ∂t
        upwind_mx = ((Mx[2:end-2] .* av_x(Vx[2:end-1]) .* (av_x(Vx[2:end-1]) .> 0.0)) .+ (Mx[3:end-1] .* av_x(Vx[2:end-1]) .* (av_x(Vx[2:end-1]) .< 0.0))) ./ dx
        dVxdx[2:end-1] .= .-(diff(upwind_mx .+ P[2:end-1], dims=1) ./ dx)
        Mx[2:end-1] .= Mx_old[2:end-1] .+ dMxdt .* dt
        Vx[2:end-1] .= Mx[2:end-1] ./ av_x(ρ)

        #dEudx[2:end-1] .= diff((Vx[2:end-1] .> 0.0) .* Vx[2:end-1] .* E[1:end-1] .+ (Vx[2:end-1] .< 0.0) .* Vx[2:end-1] .* E[2:end], dims=1) ./ dx
        #dPudx[2:end-1] .= diff((Vx[2:end-1] .> 0.0) .* Vx[2:end-1] .* P[1:end-1] .+ (Vx[2:end-1] .< 0.0) .* Vx[2:end-1] .* P[2:end], dims=1) ./ dx
        upwind_E = (Vx[3:end-2] .> 0.0) .* Vx[3:end-2] .* (P[2:end-2] + E[2:end-2]) .+ (Vx[3:end-2] .< 0.0) .* Vx[3:end-2] .* (P[3:end-1] + E[3:end-1])
        #dEdt[2:end-1] .= .-dEudx[2:end-1] .- dPudx[2:end-1]
        dEdt[3:end-2] .= .-(diff(upwind_E, dims=1) ./ dx)
        E[2:end-1] .= E_old[2:end-1] .+ dEdt[2:end-1] .* dt
        
        if wave 
            P .= c.^2.0 .* ρ
        else
            P .= (γ .- 1.0) .* (E .- 0.5 .* ρ .* av_x(Vx).^2.0)
        end
        e = P ./ (γ .- 1) ./ ρ

        t += dt
        if i % divisor == 0
            fig2 = Figure(size=(1000,800))
            ax1 = Axis(fig2[2,1], title="Pressure, time = $t")#, limits=(nothing, nothing, 0.0, 2.0))
            ax2 = Axis(fig2[1,2], title="Velocity")
            ax3 = Axis(fig2[2,2], title="Energy")
            ax4 = Axis(fig2[1,1], title="Density")
            #ylims!(ax, -1.2, 1.2)
            lines!(ax1, xc_vec[2:end-1], P[2:end-1])
            scatter!(ax1, xc_vec[2:end-1], P[2:end-1])
            lines!(ax2, xv_vec[2:end-1], Vx[2:end-1])
            scatter!(ax2, xv_vec[2:end-1], Vx[2:end-1])
            lines!(ax3, xc_vec[2:end-1], e[2:end-1])
            scatter!(ax3, xc_vec[2:end-1], e[2:end-1])
            lines!(ax4, xc_vec[2:end-1], ρ[2:end-1])
            scatter!(ax4, xc_vec[2:end-1], ρ[2:end-1])
            display(fig2)
        end
    end
end