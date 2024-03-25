using CairoMakie
using SodShockTube

function av_x(B)
    A = 0.5 .* (B[2:end,:] .+ B[1:end-1,:])
end

function upwind(u, G, dx)
    G_p = (G[2:end-1] .- G[1:end-2])./dx
    G_m = (G[2:end-1] .+ G[3:end])./dx
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

extend_vertices(x) = [x[1]; x; x[end]];

function shock_wave1D_up()
    # Physics
    Lx = 1.0                           # domain
    γ = 1.4                                # adiabatic index/ratio of specific heats
    K = 1.0e10                             # shear modulus
    ρ0 = 1.0                          # initial density at all points
    P0 = 1.0#e5                          # initial pressure at all points
    Vx0 = 0.0                          # initial velocity in x-direction for all points
    # Gaussian parameters
    A = 10.0                          # gaussian maximum amplitude
    σ = Lx * 0.04                            # standard deviation of the initial pressure distribution
    
    # Plotting parameters
    divisor = 280 

    # Numerics
    nx = 1000                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 1400                             # number of time steps

    # Grid definition
    xc = -(Lx - dx) / 2:dx:(Lx - dx) / 2        # grid nodes in x-direction
    xv = - Lx       / 2:dx: Lx       / 2        # grid vertices in x-direction

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
    #P .= P0 .+ A .* exp.(.- 1.0 .* (xc ./ σ).^2.0)       # initial pressure distribution
    P[Int((50/100)*nx):end] .= 0.1
    ρ[Int((50/100)*nx):end] .= 0.125
    c = sqrt(K / ρ0)                # speed of sound
    E .= P./((γ - 1.0)) + 0.5 .* av_x(Vx).^2
    e = P ./ (γ - 1.0) ./ ρ

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
    l0 = lines!(ax1, xc_vec, P, label="time = 0")
    push!(linplots, l0)
    lines!(ax2, xv_vec, Vx)
    lines!(ax3, xc_vec, e)
    lines!(ax4, xc_vec, ρ)
    #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(0).png", fig)
    #display(fig)

    for i = 1:nt
        t += dt

        ρ_old .= ρ
        Vx_old .= Vx
        E_old .= E
        
        Vxdρdx .= .-av_x(Vx[2:end-1]) .* upwind_center(Vx, ρ, dx)
        ρdVxdt .= .-ρ[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx # hier
        ρ[2:end-1] .= ρ_old[2:end-1] .+ Vxdρdx .* dt .+ ρdVxdt .* dt

        dρ = av_x(ρ .- ρ_old)
        ρ_v = extend_vertices(av_x(ρ))
        Vxdρdt .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* (dρ ./ dt)
        ρdPdx .= .-(1.0 ./ av_x(ρ)) .* diff(P, dims=1) ./ dx
        ρVxdVxdx .= .-(1.0 ./ av_x(ρ)) .* (av_x(ρ) .* Vx[2:end-1]) .* diff(av_x(Vx), dims=1) ./ dx # hier
        VxdρVxdx .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* upwind(Vx, ρ_v .* Vx, dx)
        Vx[2:end-1] .= Vx_old[2:end-1] .+ ρdPdx .* dt .+ ρVxdVxdx .* dt .+ VxdρVxdx .* dt .+ Vxdρdt .* dt
        
        VxdPdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(P), dims=1) ./ dx # hier
        PdVxdx .= .-P[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        VxdEdx .= .-av_x(Vx[2:end-1]) .* upwind_center(Vx, E, dx)
        EdVxdx .= .-E[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        E[2:end-1] .= E_old[2:end-1] .+ VxdPdx .* dt .+ PdVxdx .* dt .+ VxdEdx .* dt .+ EdVxdx .* dt
        
        P .= (γ .- 1.0) .* (E .- 0.5 .* ρ.* av_x(Vx).^2)

        e = P ./ (γ - 1.0) ./ ρ

        # Analytical solution 
        # problem = ShockTubeProblem(
        #     geometry = (-(Lx - dx) / 2, (Lx - dx) / 2, 0.0), # left edge, right edge, initial shock location
        #     left_state = (ρ = 1.0, u = 0.0, p = 1.0),
        #     right_state = (ρ = 0.125, u = 0.0, p = 0.1),
        #     t = t, γ = γ)
        # positions, regions, values = solve(problem, xc);
        # e_anal = values.p./((γ-1).*values.ρ)      # it seems the e calculation is a bit different in the analytical code


        if i % divisor == 0
            fig2 = Figure(size=(1000, 800), fontsize=20)
            Label(fig2[1,1:2], "Time =  $(round(t, digits=4))", tellwidth=false, font=:bold ,fontsize=26)
            ax1 = Axis(fig2[3,1], title="Pressure", ylabel="Pressure", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            ax2 = Axis(fig2[2,2], title="Velocity", ylabel="Velocity", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            ax3 = Axis(fig2[3,2], title="Energy", ylabel="Energy", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            ax4 = Axis(fig2[2,1], title="Density", ylabel="Density", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            opts = (;linewidth = 2, color = :red)
            li2 = lines!(ax4, xc, values.ρ; opts...)
            lines!(ax2, xc, values.u; opts...)
            lines!(ax1, xc, values.p; opts...)
            lines!(ax3, xc, e_anal; opts...)
            li = lines!(ax1, xc_vec, P, label="time = $t")
            #scatter!(ax1, xc_vec, P, label="time = $t")
            push!(linplots, li)
            lines!(ax2, xv_vec, Vx)
            #scatter!(ax2, xv_vec[2:end-1], Vx[2:end-1])
            lines!(ax3, xc_vec, e)
            #scatter!(ax3, xc_vec, e)
            lines!(ax4, xc_vec, ρ)
            #scatter!(ax4, xc_vec, ρ)
            
            #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(i).png", fig2)
            #display(fig2)
            if i % nt == 0
                #Legend(fig[3,:], linplots, string.(round.(0:dt*divisor:dt*nt, digits=8)), "Total time", tellwidth = false, nbanks=Int((nt/divisor)+1))#, tellhight = false, tellwidth = false)#, orientation=:horizontal, tellhight = false, tellwidth = false)
                Legend(fig2[4,:], [li2, linplots[end]], ["analytical", "numerical"], "Total time", tellwidth = false, nbanks=Int((nt/divisor)+1))#, tellhight = false, tellwidth = false)#, orientation=:horizontal, tellhight = false, tellwidth = false)
                rowsize!(fig2.layout, 4, 50)
                #save("/home/nils/Masterthesis_code/Plots/Nils_Euler-equations_shock_wave/all_terms_upwind_sod_shock_setup/Shock_upwind_vs_analytical_corrected.png", fig2) #time_evolution
                display(fig2)
            end
        end
    end
    #Legend(fig[3,:], linplots, string.(round.(0:dt*divisor:dt*nt, digits=8)), "Total time", tellwidth = false, nbanks=Int((nt/divisor)+1))#, tellhight = false, tellwidth = false)#, orientation=:horizontal, tellhight = false, tellwidth = false)
    #rowsize!(fig.layout, 3, 40)
    #save("/home/nils/Masterthesis_code/Plots/Navier-Stokes_shock_wave/all_terms_upwind_sod_shock_setup/Shock_upwind_time_evolution.png", fig)
    #display(fig)
end