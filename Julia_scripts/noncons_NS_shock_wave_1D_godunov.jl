using CairoMakie
using SodShockTube

function av_x(B)
    A = 0.5 .* (B[2:end,:] .+ B[1:end-1,:])
end

function upwind(u, G, dx)
    G_p = (G[2:end-1] .- G[1:end-2])./dx
    G_m = (G[2:end-1] .+ G[3:end])./dx
    #@infiltrate
    Gflux = (G_p.*(u[2:end-1] .> 0) .+ G_m.*(u[2:end-1] .< 0))

    return Gflux
end

function upwind_center(u, G, dx)
    G_p = (G[2:end-1] .- G[1:end-2])./dx
    G_m = (G[2:end-1] .+ G[3:end])./dx
    u_c = av_x(u)
    Gflux = (G_p.*(u_c[2:end-1] .> 0) .+ G_m.*(u_c[2:end-1] .< 0))

    return Gflux
end

function godunov_flux_bar(u, F, dt, dx)
    u_n = 0.5 .* (u[2:end] .+ u[1:end-1]) .- (dt ./ dx) .* (F[2:end] .- F[1:end-1])
end

extend_vertices(x) = [x[1]; x; x[end]];

# Sod shock paper 
function shock_wave1D_god_v1()
    # Physics
    Lx = 1.0                           # domain
    γ = 1.4                                # adiabatic index/ratio of specific heats
    K = 1.0e10                             # shear modulus
    ρ0 = 1.0                          # initial density at all points
    P0 = 1.0#e5                          # initial pressure at all points
    Vx0 = 0.0                          # initial velocity in x-direction for all points
    β = 1.0 / K
    # Gaussian parameters
    A = 10.0                          # gaussian maximum amplitude
    σ = Lx * 0.04                            # standard deviation of the initial pressure distribution
    
    # Plotting parameters
    divisor = 10

    # Numerics
    nx = 100                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 1400                             # number of time steps

    # Grid definition
    xc = -(Lx - dx) / 2:dx:(Lx - dx) / 2        # grid nodes in x-direction
    xv = - Lx       / 2:dx: Lx       / 2        # grid vertices in x-direction

    # Allocations
    ρ = ones(nx + 2) .* ρ0
    ρ_old = zeros(nx + 2)
    P = ones(nx + 2) .* P0
    Vx = zeros(nx + 3)
    Vx_old = zeros(nx + 3)
    Mx = zeros(nx + 3)
    Mx_old = zeros(nx + 3)
    E = zeros(nx + 2)
    Vxdρdx = zeros(nx+1)
    ρdVxdt = zeros(nx+1)
    Vxdρdt = zeros(nx+2)
    ρdPdx  = zeros(nx+2)
    ρVxdVxdx = zeros(nx+2)
    VxdρVxdx = zeros(nx+2)
    VxdPdx = zeros(nx + 1)
    PdVxdx = zeros(nx + 1)
    VxdEdx = zeros(nx + 1)
    EdVxdx = zeros(nx + 1)

    fρ  = zeros(nx+1)
    fMx = zeros(nx+2)
    fE  = zeros(nx+1)

    # Initial conditions
    #P .= P0 .+ A .* exp.(.- 1.0 .* (xc ./ σ).^2.0)       # initial pressure distribution
    P[Int((50/100)*nx):end] .= 0.1
    ρ[Int((50/100)*nx):end] .= 0.125
    c = sqrt(K / ρ0)                # speed of sound
    E .= P./((γ - 1.0)) + 0.5 .* av_x(Vx).^2
    e = P ./ (γ - 1.0) ./ ρ
    Mx[2:end-1] .= av_x(ρ) .* Vx[2:end-1]

    dt = 1.0e-4#8 #dx / (c * 4.0)                      # time step size
    t = 0.0                                         # initial time

    xc_vec = Array(xc)
    xv_vec = Array(xv)

    # Analytical solution 
    problem = ShockTubeProblem(
                geometry = (-(Lx - dx) / 2, (Lx - dx) / 2, 0.0), # left edge, right edge, initial shock location
                left_state = (ρ = 1.0, u = 0.0, p = 1.0),
                right_state = (ρ = 0.125, u = 0.0, p = 0.1),
                t = 0.14, γ = γ)
    positions, regions, values = solve(problem, xc);
    e_anal = values.p./((γ-1).*values.ρ)      # it seems the e calculation is a bit different in the analytical code

    linplots = []

    # Initial plotting
    fig = Figure(size=(1000, 800))
    ax1 = Axis(fig[2,1], title="Pressure", ylabel="Pressure", xlabel="Domain",)# limits=(nothing, nothing, P0-(A*3/5), P0+A))
    ax2 = Axis(fig[1,2], title="Velocity", ylabel="Velocity", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
    ax3 = Axis(fig[2,2], title="Energy", ylabel="Energy", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
    ax4 = Axis(fig[1,1], title="Density", ylabel="Density", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
    l0 = lines!(ax1, xc_vec, P[2:end-1], label="time = 0")
    push!(linplots, l0)
    lines!(ax2, xv_vec, Vx[2:end-1])
    lines!(ax3, xc_vec, e[2:end-1])
    lines!(ax4, xc_vec, ρ[2:end-1])
    #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(0).png", fig)
    #display(fig)

    for i = 1:nt
        ρ_old .= ρ
        Mx_old .= Mx
        Vx_old .= Vx

        # Fluxes F
        Vxdρdx .= .-Vx[2:end-1] .* diff(ρ, dims=1) ./ dx
        ρdVxdt .= .-av_x(ρ) .* diff(av_x(Vx), dims=1) ./ dx
        # Difference of fluxes F
        Fρ_p = (ρdVxdt[3:end]) .- (ρdVxdt[2:end-1])
        Fρ_m = (ρdVxdt[2:end-1]) .- (ρdVxdt[1:end-2])
        # Fluxes F_bar
        Fρ_bar_p = godunov_flux_bar(av_x(ρ[2:end-1]), Fρ_p, dt, dx)
        Fρ_bar_m = godunov_flux_bar(av_x(ρ[2:end-1]), Fρ_m, dt, dx)
        # Update ρ
        ρ[3:end-2] .= ρ[3:end-2] .- (dt ./ dx) .* (Fρ_bar_p .- Fρ_bar_m)

        P .= (γ .- 1.0) .* (E .- 0.5 .* ρ.* av_x(Vx).^2)

        ## Velocity formulation
        # # Prerequisites
        # dρ = (ρ .- ρ_old)
        # ρ_v = extend_vertices(av_x(ρ))
        # P_v = extend_vertices(av_x(P))
        # # Fluxes F
        # Vxdρdt .= .-(1.0 ./ ρ) .* av_x(Vx) .* (dρ ./ dt)
        # ρdPdx .= .-(1.0 ./ ρ) .* diff(P_v, dims=1) ./ dx
        # ρVxdVxdx .= .-(1.0 ./ ρ) .* (ρ .* av_x(Vx)) .* diff(Vx, dims=1) ./ dx
        # VxdρVxdx .= .-(1.0 ./ ρ) .* av_x(Vx) .* diff(ρ_v .* Vx, dims=1) ./ dx
        # # Difference of fluxes F
        # FVx_p = (Vxdρdt[3:end] .+ ρdPdx[3:end] .+ ρVxdVxdx[3:end] .+ VxdρVxdx[3:end]) .- (Vxdρdt[2:end-1] .+ ρdPdx[2:end-1] .+ ρVxdVxdx[2:end-1] .+ VxdρVxdx[2:end-1])
        # FVx_m = (Vxdρdt[2:end-1] .+ ρdPdx[2:end-1] .+ ρVxdVxdx[2:end-1] .+ VxdρVxdx[2:end-1]) .- (Vxdρdt[1:end-2] .+ ρdPdx[1:end-2] .+ ρVxdVxdx[1:end-2] .+ VxdρVxdx[1:end-2])
        # # Fluxes F_bar
        # FVx_bar_p = godunov_flux_bar(av_x(Vx[2:end-1]), FVx_p, dt, dx)
        # FVx_bar_m = godunov_flux_bar(av_x(Vx[2:end-1]), FVx_m, dt, dx)
        # # Update Vx
        # Vx[3:end-2] .= Vx[3:end-2] .- (dt ./ dx) .* (FVx_bar_p .- FVx_bar_m)

        ## Momentum formulation
        # Prerequisites
        ρ_v = extend_vertices(av_x(ρ))
        P_v = extend_vertices(av_x(P))
        # Fluxes F
        ρdPdx .= .-diff(P_v, dims=1) ./ dx
        ρVxdVxdx .= .-av_x(Mx) .* diff(Vx, dims=1) ./ dx
        VxdρVxdx .= .-av_x(Vx) .* diff(Mx, dims=1) ./ dx
        # Difference of fluxes F
        FMx_p = (ρdPdx[3:end] .+ ρVxdVxdx[3:end] .+ VxdρVxdx[3:end]) .- (ρdPdx[2:end-1] .+ ρVxdVxdx[2:end-1] .+ VxdρVxdx[2:end-1])
        FMx_m = (ρdPdx[2:end-1] .+ ρVxdVxdx[2:end-1] .+ VxdρVxdx[2:end-1]) .- (ρdPdx[1:end-2] .+ ρVxdVxdx[1:end-2] .+ VxdρVxdx[1:end-2])
        # Fluxes F_bar
        FMx_bar_p = godunov_flux_bar(av_x(Mx[2:end-1]), FMx_p, dt, dx)
        FMx_bar_m = godunov_flux_bar(av_x(Mx[2:end-1]), FMx_m, dt, dx)
        # Update Mx and Vx
        Mx[3:end-2] .= Mx[3:end-2] .- (dt ./ dx) .* (FMx_bar_p .- FMx_bar_m)
        Vx[2:end-1] .= Mx[2:end-1] ./ av_x(ρ)
        
        # Fluxes F
        VxdPdx .= .-Vx[2:end-1] .* diff(P, dims=1) ./ dx
        PdVxdx .= .-av_x(P) .* diff(av_x(Vx), dims=1) ./ dx
        VxdEdx .= .-Vx[2:end-1] .* diff(E, dims=1) ./ dx
        EdVxdx .= .-av_x(E) .* diff(av_x(Vx), dims=1) ./ dx
        # Difference of fluxes F
        FE_p = (VxdPdx[3:end] .+ PdVxdx[3:end] .+ VxdEdx[3:end] .+ EdVxdx[3:end]) .- (VxdPdx[2:end-1] .+ PdVxdx[2:end-1] .+ VxdEdx[2:end-1] .+ EdVxdx[2:end-1])
        FE_m = (VxdPdx[2:end-1] .+ PdVxdx[2:end-1] .+ VxdEdx[2:end-1] .+ EdVxdx[2:end-1]) .- (VxdPdx[1:end-2] .+ PdVxdx[1:end-2] .+ VxdEdx[1:end-2] .+ EdVxdx[1:end-2])
        # Fluxes F_bar
        FE_bar_p = godunov_flux_bar(av_x(E[2:end-1]), FE_p, dt, dx)
        FE_bar_m = godunov_flux_bar(av_x(E[2:end-1]), FE_m, dt, dx)
        # Update E
        E[3:end-2] .= E[3:end-2] .- (dt ./ dx) .* (FE_bar_p .- FE_bar_m)

        # Update e (internal energy)
        e .= P ./ (γ - 1.0) ./ ρ

        t += dt
        if i % divisor == 0
            fig2 = Figure(size=(1000, 800))
            ax1 = Axis(fig2[2,1], title="Pressure, time = $(round(t, digits=4))", ylabel="Pressure", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            ax2 = Axis(fig2[1,2], title="Velocity", ylabel="Velocity", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            ax3 = Axis(fig2[2,2], title="Energy", ylabel="Energy", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            ax4 = Axis(fig2[1,1], title="Density", ylabel="Density", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            opts = (;linewidth = 2, color = :red)
            lines!(ax4, xc, values.ρ; opts...)
            lines!(ax2, xc, values.u; opts...)
            lines!(ax1, xc, values.p; opts...)
            lines!(ax3, xc, e_anal; opts...)
            li = lines!(ax1, xc_vec, P[2:end-1], label="time = $t")
            push!(linplots, li)
            lines!(ax2, xv_vec, Vx[2:end-1])
            lines!(ax3, xc_vec, e[2:end-1])
            lines!(ax4, xc_vec, ρ[2:end-1])
            
            #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(i).png", fig2)
            display(fig2)
            if i % nt == 0
                Legend(fig2[3,:], linplots, string.(round.(0:dt*divisor:dt*nt, digits=8)), "Total time", tellwidth = false, nbanks=Int((nt/divisor)+1))#, tellhight = false, tellwidth = false)#, orientation=:horizontal, tellhight = false, tellwidth = false)
                rowsize!(fig2.layout, 3, 40)
                #save("C:\\Users\\Nils\\Desktop\\Masterthesis_code\\Plots\\Navier-Stokes_shock_wave\\nonconservative\\Shock_upwind_vs_analytical.png", fig)
                display(fig2)
            end
        end
    end
end

# normal fluxes taken
function shock_wave1D_god_v2()
    # Physics
    Lx = 1.0                           # domain
    γ = 1.4                                # adiabatic index/ratio of specific heats
    K = 1.0e10                             # shear modulus
    ρ0 = 1.0                          # initial density at all points
    P0 = 1.0#e5                          # initial pressure at all points
    Vx0 = 0.0                          # initial velocity in x-direction for all points
    β = 1.0 / K
    # Gaussian parameters
    A = 10.0                          # gaussian maximum amplitude
    σ = Lx * 0.04                            # standard deviation of the initial pressure distribution
    
    # Plotting parameters
    divisor = 10

    # Numerics
    nx = 100                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 1400                             # number of time steps

    # Grid definition
    xc = -(Lx - dx) / 2:dx:(Lx - dx) / 2        # grid nodes in x-direction
    xv = - Lx       / 2:dx: Lx       / 2        # grid vertices in x-direction

    # Allocations
    ρ = ones(nx + 2) .* ρ0
    ρ_old = zeros(nx + 2)
    P = ones(nx + 2) .* P0
    Vx = zeros(nx + 3)
    Vx_old = zeros(nx + 3)
    Mx = zeros(nx + 3)
    Mx_old = zeros(nx + 3)
    E = zeros(nx + 2)
    Vxdρdx = zeros(nx+1)
    ρdVxdt = zeros(nx+1)
    Vxdρdt = zeros(nx+2)
    ρdPdx  = zeros(nx+2)
    ρVxdVxdx = zeros(nx+2)
    VxdρVxdx = zeros(nx+2)
    VxdPdx = zeros(nx + 1)
    PdVxdx = zeros(nx + 1)
    VxdEdx = zeros(nx + 1)
    EdVxdx = zeros(nx + 1)

    fρ  = zeros(nx+1)
    fMx = zeros(nx+2)
    fE  = zeros(nx+1)

    # Initial conditions
    #P .= P0 .+ A .* exp.(.- 1.0 .* (xc ./ σ).^2.0)       # initial pressure distribution
    P[Int((50/100)*nx):end] .= 0.1
    ρ[Int((50/100)*nx):end] .= 0.125
    c = sqrt(K / ρ0)                # speed of sound
    E .= P./((γ - 1.0)) + 0.5 .* av_x(Vx).^2
    e = P ./ (γ - 1.0) ./ ρ
    Mx[2:end-1] .= av_x(ρ) .* Vx[2:end-1]

    dt = 1.0e-4#8 #dx / (c * 4.0)                      # time step size
    t = 0.0                                         # initial time

    xc_vec = Array(xc)
    xv_vec = Array(xv)

    # Analytical solution 
    problem = ShockTubeProblem(
                geometry = (-(Lx - dx) / 2, (Lx - dx) / 2, 0.0), # left edge, right edge, initial shock location
                left_state = (ρ = 1.0, u = 0.0, p = 1.0),
                right_state = (ρ = 0.125, u = 0.0, p = 0.1),
                t = 0.14, γ = γ)
    positions, regions, values = solve(problem, xc);
    e_anal = values.p./((γ-1).*values.ρ)      # it seems the e calculation is a bit different in the analytical code

    linplots = []

    # Initial plotting
    fig = Figure(size=(1000, 800))
    ax1 = Axis(fig[2,1], title="Pressure", ylabel="Pressure", xlabel="Domain",)# limits=(nothing, nothing, P0-(A*3/5), P0+A))
    ax2 = Axis(fig[1,2], title="Velocity", ylabel="Velocity", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
    ax3 = Axis(fig[2,2], title="Energy", ylabel="Energy", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
    ax4 = Axis(fig[1,1], title="Density", ylabel="Density", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
    l0 = lines!(ax1, xc_vec, P[2:end-1], label="time = 0")
    push!(linplots, l0)
    lines!(ax2, xv_vec, Vx[2:end-1])
    lines!(ax3, xc_vec, e[2:end-1])
    lines!(ax4, xc_vec, ρ[2:end-1])
    #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(0).png", fig)
    #display(fig)

    for i = 1:nt
        ρ_old .= ρ
        Mx_old .= Mx
        Vx_old .= Vx

        # Fluxes F
        Vxdρdx .= .-Vx[2:end-1] .* diff(ρ, dims=1) ./ dx
        ρdVxdt .= .-av_x(ρ) .* diff(av_x(Vx), dims=1) ./ dx
        # Difference of fluxes F
        Fρ_p = (Vxdρdx[3:end] .+ ρdVxdt[3:end]) .- (Vxdρdx[2:end-1] .+ ρdVxdt[2:end-1])
        Fρ_m = (Vxdρdx[2:end-1] .+ ρdVxdt[2:end-1]) .- (Vxdρdx[1:end-2] .+ ρdVxdt[1:end-2])
        # Fluxes F_bar
        Fρ_bar_p = godunov_flux_bar(av_x(ρ[2:end-1]), Fρ_p, dt, dx)
        Fρ_bar_m = godunov_flux_bar(av_x(ρ[2:end-1]), Fρ_m, dt, dx)
        # Update ρ
        ρ[3:end-2] .= ρ[3:end-2] .- (dt ./ dx) .* (Fρ_bar_p .- Fρ_bar_m)

        P .= (γ .- 1.0) .* (E .- 0.5 .* ρ.* av_x(Vx).^2)

        ## Velocity formulation
        # # Prerequisites
        # dρ = (ρ .- ρ_old)
        # ρ_v = extend_vertices(av_x(ρ))
        # P_v = extend_vertices(av_x(P))
        # # Fluxes F
        # Vxdρdt .= .-(1.0 ./ ρ) .* av_x(Vx) .* (dρ ./ dt)
        # ρdPdx .= .-(1.0 ./ ρ) .* diff(P_v, dims=1) ./ dx
        # ρVxdVxdx .= .-(1.0 ./ ρ) .* (ρ .* av_x(Vx)) .* diff(Vx, dims=1) ./ dx
        # VxdρVxdx .= .-(1.0 ./ ρ) .* av_x(Vx) .* diff(ρ_v .* Vx, dims=1) ./ dx
        # # Difference of fluxes F
        # FVx_p = (Vxdρdt[3:end] .+ ρdPdx[3:end] .+ ρVxdVxdx[3:end] .+ VxdρVxdx[3:end]) .- (Vxdρdt[2:end-1] .+ ρdPdx[2:end-1] .+ ρVxdVxdx[2:end-1] .+ VxdρVxdx[2:end-1])
        # FVx_m = (Vxdρdt[2:end-1] .+ ρdPdx[2:end-1] .+ ρVxdVxdx[2:end-1] .+ VxdρVxdx[2:end-1]) .- (Vxdρdt[1:end-2] .+ ρdPdx[1:end-2] .+ ρVxdVxdx[1:end-2] .+ VxdρVxdx[1:end-2])
        # # Fluxes F_bar
        # FVx_bar_p = godunov_flux_bar(av_x(Vx[2:end-1]), FVx_p, dt, dx)
        # FVx_bar_m = godunov_flux_bar(av_x(Vx[2:end-1]), FVx_m, dt, dx)
        # # Update Vx
        # Vx[3:end-2] .= Vx[3:end-2] .- (dt ./ dx) .* (FVx_bar_p .- FVx_bar_m)

        ## Momentum formulation
        # Prerequisites
        ρ_v = extend_vertices(av_x(ρ))
        P_v = extend_vertices(av_x(P))
        # Fluxes F
        ρdPdx .= .-diff(P_v, dims=1) ./ dx
        ρVxdVxdx .= .-av_x(Mx) .* diff(Vx, dims=1) ./ dx
        VxdρVxdx .= .-av_x(Vx) .* diff(Mx, dims=1) ./ dx
        # Difference of fluxes F
        FMx_p = (ρdPdx[3:end] .+ ρVxdVxdx[3:end] .+ VxdρVxdx[3:end]) .- (ρdPdx[2:end-1] .+ ρVxdVxdx[2:end-1] .+ VxdρVxdx[2:end-1])
        FMx_m = (ρdPdx[2:end-1] .+ ρVxdVxdx[2:end-1] .+ VxdρVxdx[2:end-1]) .- (ρdPdx[1:end-2] .+ ρVxdVxdx[1:end-2] .+ VxdρVxdx[1:end-2])
        # Fluxes F_bar
        FMx_bar_p = godunov_flux_bar(av_x(Mx[2:end-1]), FMx_p, dt, dx)
        FMx_bar_m = godunov_flux_bar(av_x(Mx[2:end-1]), FMx_m, dt, dx)
        # Update Mx and Vx
        Mx[3:end-2] .= Mx[3:end-2] .- (dt ./ dx) .* (FMx_bar_p .- FMx_bar_m)
        Vx[2:end-1] .= Mx[2:end-1] ./ av_x(ρ)
        
        # Fluxes F
        VxdPdx .= .-Vx[2:end-1] .* diff(P, dims=1) ./ dx
        PdVxdx .= .-av_x(P) .* diff(av_x(Vx), dims=1) ./ dx
        VxdEdx .= .-Vx[2:end-1] .* diff(E, dims=1) ./ dx
        EdVxdx .= .-av_x(E) .* diff(av_x(Vx), dims=1) ./ dx
        # Difference of fluxes F
        FE_p = (VxdPdx[3:end] .+ PdVxdx[3:end] .+ VxdEdx[3:end] .+ EdVxdx[3:end]) .- (VxdPdx[2:end-1] .+ PdVxdx[2:end-1] .+ VxdEdx[2:end-1] .+ EdVxdx[2:end-1])
        FE_m = (VxdPdx[2:end-1] .+ PdVxdx[2:end-1] .+ VxdEdx[2:end-1] .+ EdVxdx[2:end-1]) .- (VxdPdx[1:end-2] .+ PdVxdx[1:end-2] .+ VxdEdx[1:end-2] .+ EdVxdx[1:end-2])
        # Fluxes F_bar
        FE_bar_p = godunov_flux_bar(av_x(E[2:end-1]), FE_p, dt, dx)
        FE_bar_m = godunov_flux_bar(av_x(E[2:end-1]), FE_m, dt, dx)
        # Update E
        E[3:end-2] .= E[3:end-2] .- (dt ./ dx) .* (FE_bar_p .- FE_bar_m)

        # Update e (internal energy)
        e .= P ./ (γ - 1.0) ./ ρ

        t += dt
        if i % divisor == 0
            fig2 = Figure(size=(1000, 800))
            ax1 = Axis(fig2[2,1], title="Pressure, time = $(round(t, digits=4))", ylabel="Pressure", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            ax2 = Axis(fig2[1,2], title="Velocity", ylabel="Velocity", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            ax3 = Axis(fig2[2,2], title="Energy", ylabel="Energy", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            ax4 = Axis(fig2[1,1], title="Density", ylabel="Density", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            opts = (;linewidth = 2, color = :red, label="Analytical")
            lines!(ax4, xc, values.ρ; opts...)
            lines!(ax2, xc, values.u; opts...)
            lines!(ax1, xc, values.p; opts...)
            lines!(ax3, xc, e_anal; opts...)
            li = lines!(ax1, xc_vec, P[2:end-1], label="Godunov")#"time = $t")
            push!(linplots, li)
            lines!(ax2, xv_vec, Vx[2:end-1])
            lines!(ax3, xc_vec, e[2:end-1])
            lines!(ax4, xc_vec, ρ[2:end-1])
            
            save("/home/nils/Masterthesis_code/Plots/Nils_Euler-equations_shock_wave/godunov_try/$(i).png", fig2)
            display(fig2)
            if i % nt == 0
                Legend(fig2[3,:], linplots, string.(round.(0:dt*divisor:dt*nt, digits=8)), "Total time", tellwidth = false, nbanks=Int((nt/divisor)+1))#, tellhight = false, tellwidth = false)#, orientation=:horizontal, tellhight = false, tellwidth = false)
                rowsize!(fig2.layout, 3, 40)
                #save("C:\\Users\\Nils\\Desktop\\Masterthesis_code\\Plots\\Navier-Stokes_shock_wave\\nonconservative\\Shock_upwind_vs_analytical.png", fig)
                display(fig2)
            end
        end
    end
end