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

function ac_wave1D_up()
    # Physics
    Lx = 1.0                                # domain
    γ = 1.4                                 # ratio of specific heats
    K = 1.0e10                              # bulk modulus
    ρ0 = 3.0e3                              # initial density at all points
    P0 = 1.0e6                              # initial pressure at all points
    # Gaussian parameters
    A = 1.0e7                               # maximum gaussian amplitude
    σ = Lx * 0.04                           # standard deviation of the initial pressure distribution
    
    # Plotting parameters
    divisor = 3000                          # after how much iterations plots shall be shown

    # Numerics
    nx = 200                                # number of nodes in x
    dx = Lx / nx                            # step size in x
    nt = 9000                               # number of time steps
    dt = 1.0e-8                             # time step size

    # Grid definition
    xv = 0.0:dx:Lx                          # grid nodes in x-direction
    xc = av_x(xv)                           # grid vertices in x-direction

    # Allocations
    ρ               = ones(nx) .* ρ0
    ρ_old           = zeros(nx)
    P               = ones(nx) 
    P_old           = ones(nx) 
    Vx              = zeros(nx + 1)
    Vx_old          = zeros(nx + 1)

    Vxdρdx          = zeros(nx - 2)
    ρdVxdt          = zeros(nx - 2)

    Vxdρdt          = zeros(nx - 1)
    ρdPdx           = zeros(nx - 1)
    ρVxdVxdx        = zeros(nx - 1)
    VxdρVxdx        = zeros(nx - 1)

    EdρVxdx         = zeros(nx - 2)
    VxdPdx          = zeros(nx - 2)
    PdVxdx          = zeros(nx - 2)
    ρVxdEdx         = zeros(nx - 2)

    E               = zeros(nx)
    e               = zeros(nx)

    idx_solid = Int(floor((33/100)*nx))
    idx_solidVx = Int(floor((33/100)*(nx+1)))
    idx_chamber = Int(idx_solid+1)
    idx_chamberVx = Int(idx_solid+1)
    idx_chamber2 = Int(floor((65/100)*(nx+1)))
    idx_chamber2Vx = Int(floor((65/100)*(nx+1)))
    idx_air = Int(idx_chamber2+1)
    idx_airVx = Int(idx_chamber2+1)

    # Initial conditions
    c = sqrt(K / ρ0)                                        # speed of sound
    dP = A .* exp.(.- 1.0 .* ((xc .- 0.5.*Lx) ./ σ).^2.0)   # initial pressure difference
    P .= P0 .+ dP                                           # initial total pressure 
    dρ = dP ./ c^2.0                                        # initial density difference
    ρ .+= dρ                                                # initial total density
    e = P ./ (γ - 1.0) ./ ρ                                 # initial internal energy

    t = 0.0                                                 # initial time

    # Domain for plotting
    xc_vec = Array(xc)
    xv_vec = Array(xv)

    # Vector for Legend entries
    linplots = []

    # Initial plotting
    fig = Figure(size=(1000, 800), fontsize=20)
    # Label(fig[1,1:2], "Time = 0.0", tellwidth=false, font=:bold ,fontsize=30)
    #ax1 = Axis(fig[1,1], title="Pressure", ylabel="Pressure", xlabel="Domain")
    #ax2 = Axis(fig[3,1], title="Velocity", ylabel="Velocity", xlabel="Domain")
    #l0 = lines!(ax1, xc_vec, P, label="time = 0")
    #push!(linplots, l0)
    #lines!(ax2, xv_vec, Vx)
    ax1 = Axis(fig[1,1], title="Density", ylabel="Density", xticklabelsvisible=false, xticksvisible=false)
    ax2 = Axis(fig[1,2], title="Velocity",ylabel="Velocity", xticklabelsvisible=false, xticksvisible=false)
    ax3 = Axis(fig[2,1], title="Pressure", xlabel="Domain", ylabel="Pressure")
    ax4 = Axis(fig[2,2], title="Energy", xlabel="Domain", ylabel="Energy")
    l0 = lines!(ax1, xc_vec, ρ)
    push!(linplots, l0)
    lines!(ax2, xv_vec, Vx)
    lines!(ax3, xc_vec, P)
    lines!(ax4, xc_vec, e)
    #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(0).png", fig)
    display(fig)

    for i = 1:nt
        # Memonry variables
        ρ_old .= ρ
        Vx_old .= Vx
        P_old .= P

        # Mass conservation
        Vxdρdx .= .-av_x(Vx[2:end-1]) .* upwind_center(Vx, ρ, dx)
        ρdVxdt .= .-ρ[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        dρ = Vxdρdx .* dt .+ ρdVxdt .* dt
        ρ[2:end-1] .= ρ_old[2:end-1] .+ dρ

        # Momentum conservation
        ρ_v = extend_vertices(av_x(ρ))
        Vxdρdt .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* ((av_x(ρ) .- av_x(ρ_old)) ./ dt)
        ρdPdx .= .-(1.0 ./ av_x(ρ)) .* diff(P, dims=1) ./ dx
        ρVxdVxdx .= .-(1.0 ./ av_x(ρ)) .* (av_x(ρ) .* Vx[2:end-1]) .* diff(av_x(Vx), dims=1) ./ dx
        VxdρVxdx .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* upwind(Vx, ρ_v .* Vx, dx)
        Vx[2:end-1] .= Vx[2:end-1] .+ ρdPdx .* dt .+ ρVxdVxdx .* dt .+ VxdρVxdx .* dt .+ Vxdρdt .* dt

        # Energy conservation
        EdρVxdx.= .-E[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        VxdPdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(P), dims=1) ./ dx
        PdVxdx .= .-P[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        ρVxdEdx .= .-av_x(Vx[2:end-1]) .* upwind_center(Vx, E, dx)
        E[2:end-1] .= E[2:end-1] .+ EdρVxdx .* dt .+ VxdPdx .* dt .+ PdVxdx .* dt .+ ρVxdEdx .* dt

        # Equation of state for pressure difference and total pressure update
        dP = c.^2.0 .* dρ #.- c.^2.0 .* ρ0 .+ P0
        P[2:end-1] .+= dP
        
        # Internal energy caluclation
        e = P ./ (γ - 1.0) ./ ρ

        t += dt
        # Plotting for time progression
        if i % divisor == 0
            # fig2 = Figure(size=(1000, 800), fontsize=20)
            # Label(fig2[1,1:2], "Time = $(round(t, digits=8))", tellwidth=false, font=:bold ,fontsize=30)
            #ax1 = Axis(fig2[1,1], title="Pressure, time = $t")
            #ax2 = Axis(fig2[2,1], title="Velocity")
            #li = lines!(ax1, xc_vec, P, label="time = $t")
            #push!(linplots, li)
            #lines!(ax2, xv_vec[2:end-1], Vx[2:end-1])
            # ax1 = Axis(fig2[1,1], title="Density, time = ")
            # ax2 = Axis(fig2[1,2], title="Velocity")
            # ax3 = Axis(fig2[2,1], title="Pressure")
            # ax4 = Axis(fig2[2,2], title="Energy")
            l0 = lines!(ax1, xc_vec, ρ)
            push!(linplots, l0)
            lines!(ax2, xv_vec, Vx)
            lines!(ax3, xc_vec, P)
            lines!(ax4, xc_vec, e)
            #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(i).png", fig2)
            display(fig)
        end
    end
    # Plotting legend after all iterations
    Legend(fig[3,:], linplots, string.(round.(0:dt*divisor:dt*nt, digits=8)), "Total time", nbanks=Int(floor((nt/divisor)+1)), tellhight = false, tellwidth = false)
    rowsize!(fig.layout, 3, 60)
    #save("Pictures/4_in_one_noncons_acoustic_upwind_withenergy.png", fig)
    display(fig)
end