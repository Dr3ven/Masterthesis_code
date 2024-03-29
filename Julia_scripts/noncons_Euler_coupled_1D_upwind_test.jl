using CairoMakie
using Infiltrator
using SodShockTube

function av_x(B)
    A = 0.5 .* (B[2:end] .+ B[1:end-1])
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

function coupled_wave1D_up()
    # Physics
    Lx = 1000.0                           # domain
    γ = 1.4                                # adiabatic index/ratio of specific heats
    K = 1.0e10                             # shear modulus
    β = 1.0 / K
    ρ0 = 2800.0                          # initial density at all points
    P0_a = 1.0e5                          # initial pressure at all points
    P0_s = 1.0e6
    P0_c = P0_s + 1.0e7
    Vx0 = 0.0                          # initial velocity in x-direction for all points
    # Gaussian parameters
    A = 10.0                          # gaussian maximum amplitude
    σ = Lx * 0.04                            # standard deviation of the initial pressure distribution
    
    # Plotting parameters
    divisor = 10
    plotlegend = false

    # Numerics
    nx = 200                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 700000                             # number of time steps

    # Grid definition
    xc = 0:dx:Lx-dx        # grid nodes in x-direction
    xv = 0:dx:Lx        # grid vertices in x-direction

    # Allocations
    #β = zeros(nx)
    ρ = ones(nx)
    ρ_old = zeros(nx)
    ρ_t_av = zeros(nx)
    P = ones(nx)
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

    idx_solid = Int(floor((33/100)*nx))
    idx_chamber = idx_solid+1
    idx_chamber2 = Int(floor((65/100)*nx))
    idx_air = idx_chamber2+1

    masksolid = ones(nx)
    masksolidVx = ones(nx+1)
    maskair = masksolid
    maskairVx = masksolidVx
    maskair[1:idx_solid] .= 0.0
    maskairVx[1:idx_solid] .= 0.0
    masksolid[idx_chamber:end] .= 0.0
    masksolidVx[idx_chamber:end] .= 0.0
    
    # Initial conditions
    #P .= P0 .+ A .* exp.(.- 1.0 .* (xc ./ σ).^2.0)       # initial pressure distribution
    #depth_solid = 67:-1:1
    depth_air = 1:1:length(idx_air:nx)

    # IC for solid
    # Constant solid density ρ = 3000.0
    # Lithosstatic pressure P = P0 + ρgh
    ρ[1:idx_solid] .= ρ0
    P[1:idx_solid] .= P0_s #.+ (ρ[1:idx_solid] .* 9.81 .* reverse(Vector(34:1:66)))
    # IC for chamber
    # Constant chamber ideal gas density ρ = 2800.0 
    # Ideal gas pressure P = P0 *  ρ / ρ0
    ρ[idx_chamber:idx_chamber2] .= ρ0 .+ 3.0
    P[idx_chamber:idx_chamber2] .= P0_c #.* ρ[idx_chamber:idx_chamber2] ./ 2800 
    E[idx_chamber:idx_chamber2] .= P[idx_chamber:idx_chamber2] ./ ((γ - 1.0)) .+ 0.5 .* av_x(Vx)[idx_chamber:idx_chamber2].^2
    e[idx_chamber:idx_chamber2] .= P[idx_chamber:idx_chamber2] ./ (γ - 1.0) ./ ρ[idx_chamber:idx_chamber2]
    #P[idx_chamber:idx_chamber2] .= (γ .- 1.0) .* (E[idx_chamber:idx_chamber2] .- 0.5 .* ρ[idx_chamber:idx_chamber2] .* av_x(Vx)[idx_chamber:idx_chamber2].^2)#    P0_c #.* ρ[idx_chamber:idx_chamber2] ./ 2800 
    # IC for air
    # Constant ideal gas density ρ = 1.225
    # Ideal gas pressure P = P0 * exp(-ghM/RT)
    ρ[idx_air:end] .= 1.225
    P[idx_air:end] .= P0_a .* exp.(-9.81 .* depth_air .* 0.028 ./ 288.0 ./ 8.314)
    c = sqrt(K / ρ0)                # speed of sound
    E .= P ./ ((γ - 1.0)) + 0.5 .* av_x(Vx).^2
    e .= P ./ (γ - 1.0) ./ ρ

    dt = 1.0e-4#8 #dx / (c * 4.0)                      # time step size
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
        #β .= 1.0 ./ P
        ρ_old .= ρ
        Mx_old .= Mx
        P_old .= P
        c_var_s = sqrt(K ./ maximum(ρ[1:idx_solid]))
        c_var_c = sqrt(K ./ maximum(ρ[idx_chamber:end]))

        Vxdρdx .= .-av_x(Vx[2:end-1]) .* upwind_center(Vx, ρ, dx) #.* maskair[2:end-1]
        ρdVxdt .= .-ρ[2:end-1] .* (diff(Vx[2:end-1], dims=1) ./ dx) #.* maskair[2:end-1]
        dρ = Vxdρdx .* dt .+ ρdVxdt .* dt
        ρ[2:end-1] .= ρ[2:end-1] .+ dρ
        ρ_t_av .= (ρ .+ ρ_old) .* 0.5

        ρ[1] = ρ[2]
        ρ[end] = ρ[end-1]
        
        # Velocity formulation
        dρ_t = av_x(ρ .- ρ_old)
        ρ_v = extend_vertices(av_x(ρ))
        Vxdρdt .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* (dρ_t ./ dt)
        ρdPdx .= .-(1.0 ./ av_x(ρ)) .* diff(P, dims=1) ./ dx
        ρVxdVxdx .= .-(1.0 ./ av_x(ρ)) .* (av_x(ρ) .* Vx[2:end-1]) .* diff(av_x(Vx), dims=1) ./ dx
        VxdρVxdx .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* upwind(Vx, ρ_v .* Vx, dx) .* maskairVx[2:end-1]
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
        E[end] = 0#E[end-1]

        # Momentum formulation
        # EdρVxdx.= .-E[2:end-1] .* upwind_center(Vx, av_x(Mx), dx)
        # VxdPdx .= .-av_x(Vx[2:end-1]) .* upwind_center(Vx, P, dx)
        # PdVxdx .= .-P[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        # ρVxdEdx .= .-av_x(Vx[2:end-1]) .* upwind_center(Vx, E, dx)
        # E[2:end-1] .= E[2:end-1] .+ EdρVxdx .* dt .+ VxdPdx .* dt .+ PdVxdx .* dt .+ ρVxdEdx .* dt

        if any(ρ .< 0.0)
            println("Negative density at iteration $i, time $t")
            findρ = findall(x-> x .< 0.0, ρ)
            #println("Index $findρ\nDensity: $(ρ[findρ])")
            break
        elseif any(P .< -1.0e6)
            findP = findall(x-> x .< -1.0e6, P)
            #println("Negative pressure at iteration $i, time $t")
            #println("Index $findP\nPressure: $(P[findP])")
            #break
        elseif any(abs.(Vx) .> 10000.0)
            println("Unrealistic velocity at iteration $i, time $t")
            findVx = findall(x-> x .> 10000.0, abs.(Vx))
            #println("Index $findVx\nVelocity: $(Vx[findVx])")
            #break
        #elseif any(e .< 0.0)
            #println("Negative energy at time $t")
            #E_ind = findall(x-> x .< 0.0, E)
            #@show E[E_ind]
            #return E
        end
        depth_cham_air = 1:1:length(idx_chamber:nx)
        P[idx_chamber+1:end] .= (γ .- 1.0) .* (E[idx_chamber+1:end] .- 0.5 .* ρ[idx_chamber+1:end] .* av_x(Vx)[idx_chamber+1:end].^2)
        P[idx_chamber] = P[idx_chamber+1]
        #L_t = E[idx_chamber-1:idx_chamber+1]
        #R_t = 0.5 .* ρ[idx_chamber-1:idx_chamber+1] .* av_x(Vx)[idx_chamber-1:idx_chamber+1].^2

        dP = dρ[1:idx_solid-1] .* c.^2.0
        P[2:idx_solid] .= P[2:idx_solid] .+ dP 
        #P[1:idx_solid] .= .-K .* log.(ρ[1:idx_solid] ./ ρ0) .+ P0_s

        #P[1] = 0#P[2]
        #P[end] = 0#P[end-1]

        e = P ./ (γ - 1.0) ./ ρ

        t += dt
        #if t > 0.0015
        #    divisor = 1
        #end
        if i % divisor == 0
            println("#------- Iteration $i-------#")
            @show c_var_s
            @show c_var_c
            mach = maximum(abs.(Vx[idx_chamber-1:end])) / 340.0
            @show mach
            @show av_x(Vx)[idx_chamber-1:idx_chamber+1].^2
            @show ρ[idx_chamber-1:idx_chamber+1]
            fig2 = Figure(size=(1000, 800))
            ax1 = Axis(fig2[2,1], title="Pressure, time = $(round(t, digits=8))")#, limits=(nothing, nothing, -0.25, 0.25))
            ax2 = Axis(fig2[1,2], title="Velocity")#, limits=(nothing, nothing, -0.25, 0.25))
            ax3 = Axis(fig2[2,2], title="Energy", ylabel="Energy", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            ax4 = Axis(fig2[1,1], title="Density", ylabel="Density", xlabel="Domain")#, limits=(nothing, nothing,-50.0, 3050.0))
            #opts = (;linewidth = 2, color = :red)
            # lines!(ax4, xc, values.ρ; opts...)
            # lines!(ax2, xc, values.u; opts...)
            # lines!(ax1, xc, values.p; opts...)
            # lines!(ax3, xc, e_anal; opts...)
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