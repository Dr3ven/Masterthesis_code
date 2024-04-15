using CairoMakie


function av_x(B)
    A = 0.5 .* (B[2:end] .+ B[1:end-1])
end

function upwind(u, G, dx)
    G_p = (G[2:end-1] .- G[1:end-2])./dx
    G_m = (G[3:end] .- G[2:end-1])./dx
    Gflux = (G_p.*(u[2:end-1] .> 0) .+ G_m.*(u[2:end-1] .< 0))

    return Gflux
end

function upwind_center(u, G, dx)
    G_p = (G[2:end-1] .- G[1:end-2])./dx
    G_m = (G[3:end] .- G[2:end-1])./dx
    u_c = av_x(u)
    Gflux = (G_p.*(u_c[2:end-1] .> 0) .+ G_m.*(u_c[2:end-1] .< 0))

    return Gflux
end

extend_vertices(x) = [x[1]; x; x[end]];

function coupled_wave1D_up()
    # Physics
    Lx = 1000.0                                 # domain
    γ = 1.4                                     # adiabatic index/ratio of specific heats
    K = 1.0e10                                  # shear modulus
    ρ0 = 2800.0                                 # initial density at all points
    P0_a = 1.0e5                                # initial pressure for the air part
    P0_s = 1.0e6                                # initial pressure for the solid part
    P0_c = P0_s                                 # initial pressure for the chamber part
    c = sqrt(K / ρ0)                            # speed of sound
    # Gaussian parameters
    A = 1.0e5                                   # gaussian maximum amplitude
    σ = Lx * 0.04                               # standard deviation of the initial pressure distribution

    # Numerics
    nx = 1000                                   # number of nodes in x
    dx = Lx / nx                                # step size in x
    dt = 1.0e-8 #dx / (c * 4.0)                 # time step size
    nt = 700000000                              # number of time steps

    # Plotting parameters
    divisor = 1000000

    # Grid definition
    xc = 0:dx:Lx-dx                             # grid nodes in x-direction
    xv = 0:dx:Lx                                # grid vertices in x-direction

    # Allocations
    ρ = ones(nx)
    ρ_old = zeros(nx)
    ρ_t_av = zeros(nx)
    P = ones(nx)
    P_old = ones(nx)
    Vx = zeros(nx + 1)
    Mx = zeros(nx + 1)
    Mx_old = zeros(nx + 1)
    E = zeros(nx)
    e = zeros(nx)
    Vxdρdx = zeros(nx - 2)
    ρdVxdx = zeros(nx - 2)
    Vxdρdt = zeros(nx - 1)
    ρdPdx = zeros(nx - 1)
    ρVxdVxdx = zeros(nx - 1)
    VxdρVxdx = zeros(nx - 1)

    EdρVxdx = zeros(nx - 2)
    VxdPdx = zeros(nx - 2)
    PdVxdx = zeros(nx - 2)
    ρVxdEdx = zeros(nx - 2)

    # Indices for borders of the three phases
    idx_solid = Int(floor((33/100)*nx))
    idx_solidVx = Int(floor((33/100)*(nx+1)))
    idx_chamber = Int(idx_solid+1)
    idx_chamberVx = Int(idx_solidVx+1)
    idx_chamber2 = Int(floor((65/100)*(nx)))
    idx_chamber2Vx = Int(floor((65/100)*(nx+1)))
    idx_air = Int(idx_chamber2+1)
    idx_airVx = Int(idx_chamber2Vx+1)

    # Mask arrays
    masksolid = ones(nx)
    masksolidVx = ones(nx+1)
    masksolid[idx_chamber:end] .= 0.0
    masksolidVx[idx_chamber:end] .= 0.0
    maskair = masksolid .- 1.0
    maskair = abs.(maskair)
    maskairVx = masksolidVx .- 1.0
    maskairVx = abs.(maskairVx)

    # IC for solid
    ρ[1:idx_solid] .= ρ0                                                    # initial density of the solid 
    P[1:idx_solid] .= P0_s                                                  # initial pressure of the solid
    
    # IC for chamber
    dPc = A .* exp.(.- 1.0 .* ((xc .- 0.5*Lx)./ σ).^2.0)                    # initial Gaussian pressure difference distribution
    P[idx_chamber:idx_chamber2] .= P0_c .+ dPc[idx_chamber:idx_chamber2]    # initial pressure of the chamber
    dρc = dPc[idx_chamber:idx_chamber2] ./ c.^2.0                           # initial Gaussian density difference
    ρ[idx_chamber:idx_chamber2] .= ρ0 .+ dρc                                # initial density of the chamber
    
    # IC for air
    ρ[idx_air:end] .= ρ0                                                    # initial density of the air
    P[idx_air:end] .= P0_a                                                  # initial pressure of the air
    E .= P ./ ((γ - 1.0)) + 0.5 .* ρ .* av_x(Vx).^2                         # initial total energy of the whole domain
    e .= P ./ (γ - 1.0) ./ ρ                                                # initial internal energy of the whole domain

    # Total time
    t = 0.0                                                                 # initial time

    # Domain arrays for plotting
    xc_vec = Array(xc)
    xv_vec = Array(xv)

    # Initial plotting
    fig = Figure(size=(1000, 900), fontsize=20)
    Label(fig[1,1:2], "Time = 0.0", tellwidth=false, font=:bold ,fontsize=26)
    ax1 = Axis(fig[3,1], title="Pressure", ylabel="Pressure", xlabel="Domain",)
    ax2 = Axis(fig[2,2], title="Velocity", ylabel="Velocity", xticklabelsvisible=false, xticksvisible=false)
    ax3 = Axis(fig[3,2], title="Energy", ylabel="Energy", xlabel="Domain")
    ax4 = Axis(fig[2,1], title="Density", ylabel="Density", xticklabelsvisible=false, xticksvisible=false)
    lines!(ax1, xc_vec, P, label="time = 0")
    lines!(ax2, xv_vec, Vx)
    lines!(ax3, xc_vec, e)
    lines!(ax4, xc_vec, ρ)

    # Borders between the three parts in all intitial plots
    vlines!(ax1, Lx*(33/100), color=:orange)
    vlines!(ax2, Lx*(33.5/100), color=:orange)
    vlines!(ax3, Lx*(33/100), color=:orange)
    vlines!(ax4, Lx*(33/100), color=:orange)
    vlines!(ax1, Lx*(65/100), color=:cyan)
    vlines!(ax2, Lx*(65.5/100), color=:cyan)
    vlines!(ax3, Lx*(65/100), color=:cyan)
    vlines!(ax4, Lx*(65/100), color=:cyan)

    save("iteration_0_lines.png", fig)
    # display(fig)

    for i = 1:nt
        #β .= 1.0 ./ P
        ρ_old .= ρ
        Mx_old .= Mx
        P_old .= P

        Vxdρdx .= .-av_x(Vx[2:end-1]) .* upwind_center(Vx, ρ, dx)
        ρdVxdx .= .-ρ[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        dρ = Vxdρdx .* dt .+ ρdVxdx .* dt
        ρ[2:end-1] .= ρ_old[2:end-1] .+ dρ
        ρ_t_av .= (ρ .+ ρ_old) .* 0.5
        
        # Velocity formulation
        dρ_t = av_x(ρ .- ρ_old)
        ρ_v = extend_vertices(av_x(ρ))
        Vxdρdt .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* (dρ_t ./ dt)
        ρdPdx .= .-(1.0 ./ av_x(ρ)) .* diff(P, dims=1) ./ dx
        ρVxdVxdx .= .-(1.0 ./ av_x(ρ)) .* (av_x(ρ) .* Vx[2:end-1]) .* diff(av_x(Vx), dims=1) ./ dx
        VxdρVxdx .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* upwind(Vx, ρ_v .* Vx, dx)
        Vx[2:end-1] .= Vx[2:end-1] .+ ρdPdx .* dt .+ ρVxdVxdx .* dt .+ VxdρVxdx .* dt .+ Vxdρdt .* dt

        # Velocity formulation
        EdρVxdx.= .-E[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        VxdPdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(P), dims=1) ./ dx
        PdVxdx .= .-P[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        ρVxdEdx .= .-av_x(Vx[2:end-1]) .* upwind_center(Vx, E, dx)
        E[2:end-1] .= E[2:end-1] .+ EdρVxdx .* dt .+ VxdPdx .* dt .+ PdVxdx .* dt .+ ρVxdEdx .* dt1

        # Internal energy calculations
        e = P ./ (γ - 1.0) ./ ρ

        # Equation of state for the pressure difference and the updated pressure
        dP = dρ .* c.^2.0#
        P[2:end-1] .+= dP       

        t = round(t + dt, digits=Int(abs(log10(dt))))

        # Plotting of the progression in time
        if i % divisor == 0 || i <= 10
            println("#------- Iteration $i-------#")

            # Plots domain vs. density, velocity, pressure and interal energy
            fig2 = Figure(size=(1000, 900), fontsize=20)
            Label(fig2[1,1:2], "Time = $(round(t, digits=8))", tellwidth=false, font=:bold ,fontsize=26)
            ax1 = Axis(fig2[3,1], title="Pressure", ylabel="Pressure",xlabel="Domain")
            ax2 = Axis(fig2[2,2], title="Velocity", ylabel="Velocity", xticklabelsvisible=false, xticksvisible=false)
            ax3 = Axis(fig2[3,2], title="Energy", ylabel="Energy", xlabel="Domain")
            ax4 = Axis(fig2[2,1], title="Density", ylabel="Density", xticklabelsvisible=false, xticksvisible=false)
            lines!(ax1, xc_vec, P, label="time = $t")
            lines!(ax2, xv_vec, Vx)
            lines!(ax3, xc_vec, e)
            lines!(ax4, xc_vec, ρ)

            # Borders of between the three phases in all plots
            vlines!(ax1, Lx*(33/100), color=:orange)
            vlines!(ax2, Lx*(33.5/100), color=:orange)
            vlines!(ax3, Lx*(33/100), color=:orange)
            vlines!(ax4, Lx*(33/100), color=:orange)
            vlines!(ax1, Lx*(65/100), color=:cyan)
            vlines!(ax2, Lx*(65.5/100), color=:cyan)
            vlines!(ax3, Lx*(65/100), color=:cyan)
            vlines!(ax4, Lx*(65/100), color=:cyan)
            
            save("iteration_$(i)_lines.png", fig2)
            # display(fig2)
            end
        end
    end
end