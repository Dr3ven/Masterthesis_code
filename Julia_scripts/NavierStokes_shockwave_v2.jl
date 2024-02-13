using CairoMakie
using Infiltrator


function av_x(B)
    A = 0.5 .* (B[2:end] .+ B[1:end-1])
end

extend_vertices(x) = [x[1]; x; x[end]];

function wave_1D()
    # Physics
    Lx = 1.0                           # domain
    P0 = 0.0                          # initial pressure at all points
    γ  = 1.4                            # adiabatic index
    σ  = 1.0                            # gaussian amplitude

    # Numerics
    nx = 100                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 1400                             # number of time steps

    # Grid definition
    xv = 0:dx:1        # grid vertices in x-direction
    xc = av_x(xv)        # grid nodes in x-direction


    # Allocations
    e = zeros(nx)
    E = zeros(nx)
    c = ones(nx)
    ρ = ones(nx)
    ρ_v = ones(nx)
    Mx = zeros(nx + 1)
    P = ones(nx) 
    Vx = zeros(nx + 1)
    #dEudx = zeros(nx)
    #dPudx = zeros(nx)
    dEdt = zeros(nx)
    dVxdx = zeros(nx + 1)
    dρdt = zeros(nx)
    dMxdt = zeros(nx + 1)

    # Initial conditions
    ρ[div(nx, 2):end] .= 0.125
    P[div(nx, 2):end] .= 0.1

    # P .= P0 .+ exp.(.-(1.0 ./ σ) .* xc.^2.0)        # initial pressure distribution
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
        
        ρ_v = extend_vertices(av_x(ρ))
        Vx_c = av_x(Vx[2:end-1])
        dρdt[2:end-1] .= sign.(Vx_c) .* ((((ρ[2:end-1] .* Vx_c) .- (ρ[1:end-2] .* av_x(Vx[1:end-2]))) ./ dx) .* (Vx_c .> 0.0)) .+ (((ρ[2:end-1] .* Vx_c) .+ (ρ[3:end] .* av_x(Vx[3:end]))./ dx) .* (Vx_c .< 0.0)) 
        ρ[2:end-1] .= ρ[2:end-1] .- dρdt[2:end-1] .* dt
        
        dVxdx[2:end-1] .= sign.(Vx[2:end-1]) .* (((((ρ_v[2:end-1] .* Vx[2:end-1].^2.0) .- (ρ_v[1:end-2] .* Vx[1:end-2].^2.0)) ./ dx) .* (Vx[2:end-1] .> 0.0)) .+ ((((ρ_v[2:end-1] .* Vx[2:end-1].^2.0) .+ (ρ_v[3:end] .* Vx[3:end].^2.0)) ./ dx) .* (Vx[2:end-1] .< 0.0)))
        dMxdt[2:end-1] .= .-dVxdx[2:end-1] .- (P[2:end] .- P[1:end-1]) ./ dx
        Mx[2:end-1] .= Mx[2:end-1] .+ dMxdt[2:end-1] .* dt

        dEdt[2:end-1] .= .-sign.(Vx_c) .* (((((Vx_c .* (E[2:end-1] .+ P[2:end-1])) .- (av_x(Vx[1:end-2]) .* (E[1:end-2] .+ P[1:end-2])))./ dx) .* (Vx_c .> 0.0)) .+ ((((Vx_c .* (E[2:end-1] .+ P[2:end-1])) .+ (av_x(Vx[3:end]) .* (E[3:end] .+ P[3:end])))./ dx) .* (Vx_c .< 0.0)))
        E[2:end-1] .= E[2:end-1] .+ dEdt[2:end-1] .* dt
        
        Vx[2:end-1] .= Mx[2:end-1] ./ av_x(ρ)
        P .= (γ .- 1.0) .* (E .- 0.5 .* ρ .* av_x(Vx))  # c.^2.0 .* ρ 

        e .= P ./ (γ .- 1) ./ ρ

        t += dt
        if i % 100 == 0
            fig2 = Figure()
            ax1 = Axis(fig2[1,1], title="Density, time = $t")#, limits=(nothing, nothing, 0.0, 2.0))
            ax2 = Axis(fig2[1,2], title="Velocity")
            ax3 = Axis(fig2[2,1], title="Pressure")
            ax4 = Axis(fig2[2,2], title="Energy")
            #ylims!(ax, -1.2, 1.2)
            scatter!(ax1, xc_vec, ρ)
            scatter!(ax2, xv_vec, Vx)
            scatter!(ax3, xc_vec, P)
            scatter!(ax4, xc_vec, e)
            display(fig2)
        end
    end
end