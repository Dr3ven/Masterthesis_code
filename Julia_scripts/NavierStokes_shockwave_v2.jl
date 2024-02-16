using CairoMakie
using Infiltrator

function av_x(B)
    A = 0.5 .* (B[2:end] .+ B[1:end-1])
end

extend_vertices(x) = [x[1]; x; x[end]];

function wave_1D(; wave=true)
    # Physics
    γ  = 1.4                            # adiabatic index
    if wave == true
        Lx = 1.0
        P0 = 1.0                          # initial pressure at all points
        σ  = Lx / 10                            # gaussian amplitude width
        A  = 1.0                          # gaussian maximum amplitude   
        c  = 1.0
        β  = 1.0
        ρ0 = 1.0
    else
        Lx = 1.0
    end

    # Numerics
    nx = 100                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 14000                             # number of time steps

    # Grid definition
    xv = 0:dx:Lx        # grid vertices in x-direction
    xc = av_x(xv)        # grid nodes in x-direction

    # Allocations
    if wave == true
        ρ = ones(nx) .* ρ0
        ρ_v = ones(nx) .* ρ0
        P = ones(nx)
    else
        ρ = ones(nx)
        ρ_v = ones(nx)
        P = ones(nx) .* P0
    end
    e = zeros(nx)
    E = zeros(nx)
    #c = ones(nx)
    Mx = zeros(nx + 1)
    Vx = zeros(nx + 1)
    #dEudx = zeros(nx)
    #dPudx = zeros(nx)
    dEdt = zeros(nx)
    dVxdx = zeros(nx + 1)
    dρdt = zeros(nx)
    dMxdt = zeros(nx + 1)

    # Initial conditions
    if wave == false
        ρ[div(nx, 2):end] .= 0.125
        P[div(nx, 2):end] .= 0.1
        E .= P ./ (γ .- 1.0) + (1.0 ./ 2.0) .* av_x(Vx).^2
    else
        P .= P0 .+ A .* exp.(.-1.0 .* ((xc .- (Lx .* 0.5)) ./ σ).^2.0)        # initial pressure distribution
        #P[Int(div(nx, 2.7)):Int(div(nx, 1.6))] .= 2.0 .* P0 
        #P[div(nx, 2):end] .= 0.1
        ρ .= P ./ c.^2.0
    end
    
    P_max = maximum(P)
    t = 0.0                              # initial time
    dt = 0.0001

    xc_vec = Vector(xc)
    xv_vec = Vector(xv)

    # Initial plotting
    fig = Figure()
    ax = Axis(fig[1,1], title="Pressure at t = $t")
    scatter!(ax, xc_vec, P)
    display(fig)

    for i = 1:nt
        #c = maximum(sqrt.(P ./ ρ))                            # speed of sound
        β = 1.0 / maximum(P)

        ρ_v = extend_vertices(av_x(ρ))
        Vx_c = av_x(Vx[2:end-1])
        dρdt[2:end-1] .= .-sign.(Vx_c) .* (((((ρ[2:end-1] .* Vx_c) .- (ρ[1:end-2] .* av_x(Vx[1:end-2]))) ./ dx) .* (Vx_c .> 0.0)) .+ (((ρ[2:end-1] .* Vx_c) .+ (ρ[3:end] .* av_x(Vx[3:end]))./ dx) .* (Vx_c .< 0.0))) 
        ρ[2:end-1] .= ρ[2:end-1] .+ dρdt[2:end-1] .* dt

        #ρ[1] = ρ[2]
        #ρ[end] = ρ[end-1]

        dVxdx[2:end-1] .= .-sign.(Vx[2:end-1]) .* (((((ρ_v[2:end-1] .* Vx[2:end-1].^2.0) .- (ρ_v[1:end-2] .* Vx[1:end-2].^2.0)) ./ dx) .* (Vx[2:end-1] .> 0.0)) .+ ((((ρ_v[2:end-1] .* Vx[2:end-1].^2.0) .+ (ρ_v[3:end] .* Vx[3:end].^2.0)) ./ dx) .* (Vx[2:end-1] .< 0.0)))
        dMxdt[2:end-1] .= dVxdx[2:end-1] .- (P[2:end] .- P[1:end-1]) ./ dx
        Mx[2:end-1] .= Mx[2:end-1] .+ dMxdt[2:end-1] .* dt

        #Mx[1] = Mx[2]
        #Mx[end] = Mx[end-1]

        dEdt[2:end-1] .= .-sign.(Vx_c) .* (((((Vx_c .* (E[2:end-1] .+ P[2:end-1])) .- (av_x(Vx[1:end-2]) .* (E[1:end-2] .+ P[1:end-2])))./ dx) .* (Vx_c .> 0.0)) .+ ((((Vx_c .* (E[2:end-1] .+ P[2:end-1])) .+ (av_x(Vx[3:end]) .* (E[3:end] .+ P[3:end])))./ dx) .* (Vx_c .< 0.0)))
        E[2:end-1] .= E[2:end-1] .+ dEdt[2:end-1] .* dt
        
        Vx[2:end-1] .= Mx[2:end-1] ./ av_x(ρ)

        #Vx[1] = Vx[2]
        #Vx[end] = Vx[end-1]

        if wave == false
            P .= (γ .- 1.0) .* (E .- 0.5 .* ρ .* av_x(Vx).^2.0)  
        else 
            P .= P .- (1.0 ./ β) .* diff(Vx, dims=1) .* (dt ./ dx) #c.^2.0 .* ρ # P0 .+ (1.0 ./ β) .* log.(ρ ./ ρ0)
        end

        e .= P ./ (γ .- 1) ./ ρ
        #@infiltrate

        t += dt
        if i % 10 == 0
            val, ind = findmax(dVxdx)

            fig2 = Figure()
            ax1 = Axis(fig2[1,1], title="Density, time = $t")#, limits=(nothing, nothing, -0.1, 2.5))
            ax2 = Axis(fig2[1,2], title="Velocity")
            ax3 = Axis(fig2[2,1], title="Pressure")#, limits=(nothing, nothing, -0.1, 2.5))
            ax4 = Axis(fig2[2,2], title="Energy")
            #ylims!(ax, -1.2, 1.2)
            scatter!(ax1, xc_vec, ρ)
            scatter!(ax1, xc_vec[ind], ρ[ind], color=:red)
            scatter!(ax2, xv_vec, Vx)
            scatter!(ax3, xc_vec, P)
            scatter!(ax4, xc_vec, e)
            save("./Pictures/sod_shock_tube_$i.png", fig2)
            display(fig2)
            @infiltrate
        end
    end
end