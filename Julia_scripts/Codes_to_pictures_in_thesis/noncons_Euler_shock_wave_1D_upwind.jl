using CairoMakie
using SodShockTube

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

function shock_wave1D_up()
    # Physics
    Lx = 1.0                                    # domain
    γ = 1.4                                     # adiabatic index/ratio of specific heats
    ρ0 = 1.0                                    # initial density
    P0 = 1.0                                    # initial pressure
    
    # Plotting parameters
    divisor = 35000 

    # Numerics
    nx = 200                                    # number of nodes in x
    dx = Lx / nx                                # step size in x
    nt = 140000                                 # number of time steps

    # Grid definition
    xv_vec = 0:dx:Lx                            # grid face nodes in x-direction
    xc_vec = av_x(xv_vec)                       # grid center nodes in x-direction

    # Allocations
    ρ = ones(nx) .* ρ0
    ρ_old = zeros(nx)
    P = ones(nx) .* P0
    Vx = zeros(nx + 1)
    Vx_old = zeros(nx + 1)
    Mx = zeros(nx + 1)
    E = zeros(nx)
    E_old = zeros(nx)
    Vxdρdx = zeros(nx - 2)
    ρdVxdt = zeros(nx - 2)
    Vxdρdt = zeros(nx - 1)
    ρdPdx = zeros(nx - 1)
    ρVxdVxdx = zeros(nx - 1)
    VxdρVxdx = zeros(nx - 1)
    VxdPdx = zeros(nx - 2)
    PdVxdx = zeros(nx - 2)
    VxdEdx = zeros(nx - 2)
    EdVxdx = zeros(nx - 2)

    # Initial conditions
    P[Int((50/100)*nx):end] .= 0.1
    ρ[Int((50/100)*nx):end] .= 0.125
    E .= P./((γ - 1.0)) + 0.5 .* ρ .* av_x(Vx).^2
    e = P ./ (γ - 1.0) ./ ρ

    dt = 1.0e-6#8 #dx / (c * 4.0)                                   # time step size
    t = 0.0                                                         # initial time

    # Analytical solution 
    problem = ShockTubeProblem(
        geometry = (0, Lx, 0.5), # left edge, right edge, initial shock location
        left_state = (ρ = 1.0, u = 0.0, p = 1.0),
        right_state = (ρ = 0.125, u = 0.0, p = 0.1),
        t = 0.14, γ = γ)
    positions, regions, values = solve(problem, xc_vec);
    e_anal = values.p./((γ-1).*values.ρ)      
    linplots = []

    # Initial plotting
    fig = Figure(size=(1000, 800), fontsize=20)
    ax1 = Axis(fig[2,1], title="Pressure", ylabel="Pressure", xlabel="Domain")
    ax2 = Axis(fig[1,2], title="Velocity", ylabel="Velocity", xticklabelsvisible=false, xticksvisible=false)
    ax3 = Axis(fig[2,2], title="Energy", ylabel="Energy", xlabel="Domain")
    ax4 = Axis(fig[1,1], title="Density", ylabel="Density", xticklabelsvisible=false, xticksvisible=false)
    display(fig)
    li2 = nothing
    for i = 1:nt
        ρ_old .= ρ
        Vx_old .= Vx
        E_old .= E
        
        # Conservation of mass
        Vxdρdx .= .-av_x(Vx[2:end-1]) .* upwind_center(Vx, ρ, dx)
        ρdVxdt .= .-ρ[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        ρ[2:end-1] .= ρ_old[2:end-1] .+ Vxdρdx .* dt .+ ρdVxdt .* dt

        # Conservation of momentum
        dρ = av_x(ρ .- ρ_old)
        ρ_v = extend_vertices(av_x(ρ))
        Vxdρdt .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* (dρ ./ dt)
        ρdPdx .= .-(1.0 ./ av_x(ρ)) .* diff(P, dims=1) ./ dx
        ρVxdVxdx .= .-(1.0 ./ av_x(ρ)) .* (av_x(ρ) .* Vx[2:end-1]) .* diff(av_x(Vx), dims=1) ./ dx
        VxdρVxdx .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* upwind(Vx, ρ_v .* Vx, dx)
        Vx[2:end-1] .= Vx_old[2:end-1] .+ ρdPdx .* dt .+ ρVxdVxdx .* dt .+ VxdρVxdx .* dt .+ Vxdρdt .* dt
        
        # Conservation of energy
        VxdPdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(P), dims=1) ./ dx
        PdVxdx .= .-P[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        VxdEdx .= .-av_x(Vx[2:end-1]) .* upwind_center(Vx, E, dx)
        EdVxdx .= .-E[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        E[2:end-1] .= E_old[2:end-1] .+ VxdPdx .* dt .+ PdVxdx .* dt .+ VxdEdx .* dt .+ EdVxdx .* dt
        
        # Equation of state for pressure of an isentropic gas
        P .= (γ .- 1.0) .* (E .- 0.5 .* ρ.* av_x(Vx).^2)

        # Internal energy calculation
        e = P ./ (γ - 1.0) ./ ρ

        t += dt

        # Plotting
        if i % divisor == 0
            # fig2 = Figure(size=(1000, 800), fontsize=20)
            # Label(fig2[1,1:2], "Time =  $(round(t, digits=4))
            # ax1 = Axis(fig2[3,1], title="Pressure", ylabel="Pressure", xlabel="Domain")
            # ax2 = Axis(fig2[2,2], title="Velocity", ylabel="Velocity", xlabel="Domain")
            # ax3 = Axis(fig2[3,2], title="Energy", ylabel="Energy", xlabel="Domain")
            # ax4 = Axis(fig2[2,1], title="Density", ylabel="Density", xlabel="Domain")
            opts = (;linewidth = 2, color = :red, label="analyt. @ 0.14")
            li2 = lines!(ax4, xc_vec, values.ρ; opts...)
            lines!(ax2, xc_vec, values.u; opts...)
            lines!(ax1, xc_vec, values.p; opts...)
            lines!(ax3, xc_vec, e_anal; opts...)
            li = lines!(ax1, xc_vec, P, label="time = $t")
            push!(linplots, li)
            lines!(ax2, xv_vec, Vx)
            lines!(ax3, xc_vec, e)
            lines!(ax4, xc_vec, ρ)
            
            #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(i).png", fig2)
            #display(fig2)
        end
    end
    # Legend plotting
    linplots = append!(linplots, [li2])
    text = string.(round.(0:dt*divisor:dt*nt, digits=8))
    text = append!(text, ["0.14 (analytical)"])
    Legend(fig[3,:], linplots, text, "Total time", tellwidth = false, nbanks=Int(floor((nt/divisor)+2)))
    rowsize!(fig.layout, 3, 50)
    #save("/local/home/nimeding/Ma2/Julia_scripts/Shock_upwind_time_evolution_corrected_n200_new.png", fig)
    display(fig)
end