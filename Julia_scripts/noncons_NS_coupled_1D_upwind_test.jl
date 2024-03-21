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

extend_vertices(x) = [x[1]; x; x[end]];

function coupled_wave1D_up()
    # Physics
    Lx = 1000.0                           # domain
    γ = 1.4                                # adiabatic index/ratio of specific heats
    K = 1.0e10                             # shear modulus
    ρ0 = 3000.0                          # initial density at all points
    P0 = 1.0e5                          # initial pressure at all points
    Vx0 = 0.0                          # initial velocity in x-direction for all points
    # Gaussian parameters
    A = 10.0                          # gaussian maximum amplitude
    σ = Lx * 0.04                            # standard deviation of the initial pressure distribution
    
    # Plotting parameters
    divisor = 10 
    plotlegend = false

    # Numerics
    nx = 100                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 500000                             # number of time steps

    # Grid definition
    xc = 0:dx:Lx-dx        # grid nodes in x-direction
    xv = 0:dx:Lx        # grid vertices in x-direction

    # Allocations
    β = zeros(nx)
    ρ = ones(nx)
    ρ_old = zeros(nx)
    ρ_t_av = zeros(nx)
    P = ones(nx) .* P0
    P_old = ones(nx)
    Vx = zeros(nx + 1)
    #Vx_old = zeros(nx + 1)
    #Vx_t_av = zeros(nx + 1)
    Mx = zeros(nx + 1)
    Mx_old = zeros(nx + 1)
    E = zeros(nx)
    e = zeros(nx)
    Vxdρdx = zeros(nx - 2)
    ρdVxdt = zeros(nx - 2)
    Vxdρdt = zeros(nx - 1)
    ρdPdx = zeros(nx - 1)
    ρVxdVxdx = zeros(nx - 1)
    VxdρVxdx = zeros(nx - 1)

    Edρdt  = zeros(nx - 2)
    EdρVxdx = zeros(nx - 2)
    VxdPdx = zeros(nx - 2)
    PdVxdx = zeros(nx - 2)
    ρVxdEdx = zeros(nx - 2)

    masksolid = ones(nx)
    maskair = masksolid
    maskair[1:Int((33/100)*nx)] .= 0.0
    masksolid[Int((34/100)*nx):end] .= 0.0
    
    # Initial conditions
    #P .= P0 .+ A .* exp.(.- 1.0 .* (xc ./ σ).^2.0)       # initial pressure distribution
    depth_solid = 33:-1:1
    depth_air = 1:1:67
    # IC for solid
    P[1:Int((33/100)*nx)] .= P0 .+ (ρ[1:Int((33/100)*nx)] .* 9.81 .* depth_solid)
    ρ[1:Int((33/100)*nx)] .= 3000.0
    # IC for chamber
    P[Int((34/100)*nx):end] .= P0 .* exp.(-9.81 .* depth_air .* 0.028 ./ 288.0 ./ 8.314)
    ρ[Int((34/100)*nx):65] .= 2800.0 ./ P0 .* P[Int((34/100)*nx):65]
    # IC for air
    ρ[Int((66/100)*nx):end] .= 1.225 ./ P0 .* P[Int((66/100)*nx):end]
    P[Int((66/100)*nx):end] .= (P[Int((66/100)*nx):end] .* ρ[Int((66/100)*nx):end]) ./ 1.225
    c = sqrt(K / ρ0)                # speed of sound
    E .= P ./ ((γ - 1.0)) + 0.5 .* av_x(Vx).^2
    e = P ./ (γ - 1.0) ./ ρ

    dt = 1.0e-6#8 #dx / (c * 4.0)                      # time step size
    t = 0.0                                         # initial time

    xc_vec = Array(xc)
    xv_vec = Array(xv)

    # Analytical solution 
    #problem = ShockTubeProblem(
    #            geometry = (-(Lx - dx) / 2, (Lx - dx) / 2, 0.0), # left edge, right edge, initial shock location
    #            left_state = (ρ = 1.0, u = 0.0, p = 1.0),
    #            right_state = (ρ = 0.125, u = 0.0, p = 0.1),
    #            t = 0.14, γ = γ)
    #positions, regions, values = solve(problem, xc);
    #e_anal = values.p./((γ-1).*values.ρ)      # it seems the e calculation is a bit different in the analytical code

    linplots = []

    # Initial plotting
    fig = Figure(size=(1000, 800))
    ax1 = Axis(fig[2,1], title="Pressure", ylabel="Pressure", xlabel="Domain",)# limits=(nothing, nothing, P0-(A*3/5), P0+A))
    ax2 = Axis(fig[1,2], title="Velocity", ylabel="Velocity", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
    ax3 = Axis(fig[2,2], title="Energy", ylabel="Energy", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
    ax4 = Axis(fig[1,1], title="Density", ylabel="Density", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
    l0 = scatter!(ax1, xc_vec, P, label="time = 0")
    push!(linplots, l0)
    scatter!(ax2, xv_vec, Vx)
    scatter!(ax3, xc_vec, e)
    scatter!(ax4, xc_vec, ρ)
    #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(0).png", fig)
    display(fig)

    for i = 1:nt
        β .= 1.0 ./ P
        ρ_old .= ρ
        Mx_old .= Mx
        P_old .= P
        c_var_s = sqrt(K ./ maximum(ρ[1:33]))
        c_var_a = sqrt(K ./ maximum(ρ[34:end]))

        Vxdρdx .= .-av_x(Vx[2:end-1]) .* upwind_center(Vx, ρ, dx) .* maskair[2:end-1]
        ρdVxdt .= .-ρ[2:end-1] .* (diff(Vx[2:end-1], dims=1) ./ dx) .* maskair[2:end-1]
        ρ[2:end-1] .= (ρ[2:end-1] .+ Vxdρdx .* dt .+ ρdVxdt .* dt ) 
        ρ_t_av .= (ρ .+ ρ_old) .* 0.5

        ρ[1] = ρ[2]
        ρ[end] = ρ[end-1]
        
        # Velocity formulation
        dρ = av_x(ρ .- ρ_old)
        ρ_v = extend_vertices(av_x(ρ))
        Vxdρdt .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* (dρ ./ dt)
        ρdPdx .= .-(1.0 ./ av_x(ρ)) .* diff(P, dims=1) ./ dx
        ρVxdVxdx .= .-(1.0 ./ av_x(ρ)) .* (av_x(ρ) .* Vx[2:end-1]) .* diff(av_x(Vx), dims=1) ./ dx
        VxdρVxdx .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* upwind(Vx, ρ_v .* Vx, dx)
        Vx[2:end-1] .= Vx[2:end-1] .+ ρdPdx .* dt .+ ρVxdVxdx .* dt .+ VxdρVxdx .* dt .+ Vxdρdt .* dt
        
        Vx[1] = Vx[2]
        Vx[end] = Vx[end-1]

        # Momentum formulation
        # ρdPdx .= .-diff(P, dims=1) ./ dx
        # ρVxdVxdx .= .-Mx[2:end-1] .* diff(av_x(Vx), dims=1) ./ dx
        # VxdρVxdx .= .-Vx[2:end-1] .* upwind(Vx, Mx, dx)
        # Mx[2:end-1] .= Mx[2:end-1] .+ ρdPdx .* dt .+ ρVxdVxdx .* dt .+ VxdρVxdx .* dt
        # Vx[2:end-1] .= Mx[2:end-1] ./ av_x(ρ)
       
        # Velocity formulation
        EdρVxdx.= .-E[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        VxdPdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(P), dims=1) ./ dx
        PdVxdx .= .-P[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        ρVxdEdx .= .-av_x(Vx[2:end-1]) .* upwind_center(Vx, E, dx)
        E[2:end-1] .= E[2:end-1] .+ EdρVxdx .* dt .+ VxdPdx .* dt .+ PdVxdx .* dt .+ ρVxdEdx .* dt

        # E[1] = E[2]    # lassen die simulation explodieren wenn Rand erreicht wird
        # E[end] = E[end-1]

        # Momentum formulation
        # EdρVxdx.= .-E[2:end-1] .* upwind_center(Vx, av_x(Mx), dx)
        # VxdPdx .= .-av_x(Vx[2:end-1]) .* upwind_center(Vx, P, dx)
        # PdVxdx .= .-P[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        # ρVxdEdx .= .-av_x(Vx[2:end-1]) .* upwind_center(Vx, E, dx)
        # E[2:end-1] .= E[2:end-1] .+ EdρVxdx .* dt .+ VxdPdx .* dt .+ PdVxdx .* dt .+ ρVxdEdx .* dt

        P[Int((34/100)*nx):end] .= (γ .- 1.0) .* (E[Int((34/100)*nx):end] .- 0.5 .* ρ[Int((34/100)*nx):end] .* (av_x(Mx)[Int((34/100)*nx):end] ./ ρ[Int((34/100)*nx):end]).^2)
        #P[1:Int((33/100)*nx)] .= ρ[1:Int((33/100)*nx)] .* c.^2.0
        P[1:Int((33/100)*nx)] .= ((1.0 ./ β[1:Int((33/100)*nx)]) .* log.(ρ[1:Int((33/100)*nx)] ./ 2800.0) .+ P[1:Int((33/100)*nx)])

        e = P ./ (γ - 1.0) ./ ρ

        t += dt
        #if t > 0.0015
        #    divisor = 1
        #end
        if i % divisor == 0
            @show c_var_s
            @show c_var_a
            mach = maximum(abs.(Vx[34:end])) / 340.0
            @show mach
            fig2 = Figure(size=(1000, 800))
            ax1 = Axis(fig2[2,1], title="Pressure, time = $(round(t, digits=8))")#, limits=(nothing, nothing, -0.25, 0.25))
            ax2 = Axis(fig2[1,2], title="Velocity")#, limits=(nothing, nothing, -0.25, 0.25))
            ax3 = Axis(fig2[2,2], title="Energy", ylabel="Energy", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            ax4 = Axis(fig2[1,1], title="Density", ylabel="Density", xlabel="Domain")#, limits=(nothing, nothing,-50.0, 3050.0))
            #opts = (;linewidth = 2, color = :red)
            #lines!(ax4, xc, values.ρ; opts...)
            #lines!(ax2, xc, values.u; opts...)
            #lines!(ax1, xc, values.p; opts...)
            #lines!(ax3, xc, e_anal; opts...)
            li = scatter!(ax1, xc_vec, P, label="time = $t")
            li = lines!(ax1, xc_vec, P, label="time = $t")
            push!(linplots, li)
            scatter!(ax2, xv_vec, Vx)
            lines!(ax2, xv_vec, Vx)
            scatter!(ax3, xc_vec, e)
            lines!(ax3, xc_vec, e)
            scatter!(ax4, xc_vec, ρ)
            lines!(ax4, xc_vec, ρ)
            
            #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(i).png", fig2)
            display(fig2)
            if i % nt == 0 && plotlegend == true
                Legend(fig2[3,:], linplots, string.(round.(0:dt*divisor:dt*nt, digits=8)), "Total time", tellwidth = false, nbanks=Int((nt/divisor)+1))#, tellhight = false, tellwidth = false)#, orientation=:horizontal, tellhight = false, tellwidth = false)
                rowsize!(fig2.layout, 3, 40)
                #save("C:\\Users\\Nils\\Desktop\\Masterthesis_code\\Plots\\Navier-Stokes_shock_wave\\nonconservative\\Shock_upwind_vs_analytical.png", fig)
                display(fig2)
            end
        end
        
    end
    
end