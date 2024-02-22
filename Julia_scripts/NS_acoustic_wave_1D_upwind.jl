using CairoMakie
using Infiltrator

function av_x(B)
    A = 0.5 .* (B[2:end,:] .+ B[1:end-1,:])
end

function upwind(u, G, dx)
    G_p = (G[2:end-1] .- G[1:end-2])./dx
    G_m = (G[3:end] .- G[2:end-1])./dx
    #@infiltrate
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
    Lx = 1.0                           # domain
    K = 1.0e10                             # shear modulus
    ρ0 = 3.0e3                          # initial density at all points
    P0 = 1.0e6                          # initial pressure at all points
    Vx0 = 0.0                          # initial velocity in x-direction for all points
    # Gaussian parameters
    A = 10.0                          # gaussian maximum amplitude
    σ = Lx * 0.04                            # standard deviation of the initial pressure distribution
    
    # Plotting parameters
    divisor = 2000 

    # Numerics
    nx = 100                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 10000                             # number of time steps

    # Grid definition
    xc = -(Lx - dx) / 2:dx:(Lx - dx) / 2        # grid nodes in x-direction
    xv = - Lx       / 2:dx: Lx       / 2        # grid vertices in x-direction

    # Allocations
    ρ = ones(nx) .* ρ0
    P = ones(nx) #.* P0
    Vx = zeros(nx + 1)
    Vxdρdx = zeros(nx - 2)
    ρdVxdt = zeros(nx - 2)
    ρdPdx = zeros(nx - 1)
    ρVxdVxdx = zeros(nx - 1)
    VxdρVxdx = zeros(nx - 1)

    # Initial conditions
    P .= P0 .+ A .* exp.(.- 1.0 .* (xc ./ σ).^2.0)       # initial pressure distribution
    #P[Int((48/100)*nx):Int((48/100)*nx)] .= P0 .+ A
    c = sqrt(K / ρ0)                # speed of sound
    ρ .= P ./ c^2.0                                 # initial density distribution
    #P .= c.^2.0 .* ρ                                 # initial pressure distribution

    dt = 5.0e-9 #dx / (c * 4.0)                      # time step size
    t = 0.0                                         # initial time

    xc_vec = Array(xc)
    xv_vec = Array(xv)

    linplots = []

    # Initial plotting
    fig = Figure(size=(600, 800))
    ax1 = Axis(fig[1,1], title="Pressure", ylabel="Pressure", xlabel="Domain",)# limits=(nothing, nothing, P0-(A*3/5), P0+A))
    ax2 = Axis(fig[3,1], title="Velocity", ylabel="Velocity", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
    l0 = lines!(ax1, xc_vec, P, label="time = 0")
    push!(linplots, l0)
    lines!(ax2, xv_vec, Vx)
    #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(0).png", fig)
    display(fig)

    for i = 1:nt
        #c = sqrt(K / maximum(ρ))                # speed of sound
        #@show c

        Vxdρdx .= .-av_x(Vx[2:end-1]) .* upwind_center(Vx, ρ, dx)
        ρdVxdt .= .-ρ[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        ρ[2:end-1] .= ρ[2:end-1] .+ Vxdρdx .* dt .+ ρdVxdt .* dt

        P .= c.^2.0 .* ρ

        ρ_v = extend_vertices(av_x(ρ))
        ρdPdx .= .-(1.0 ./ av_x(ρ)) .* diff(P, dims=1) ./ dx
        ρVxdVxdx .= .-(1.0 ./ av_x(ρ)) .* (av_x(ρ) .* Vx[2:end-1]) .* diff(av_x(Vx), dims=1) ./ dx
        VxdρVxdx .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* upwind(Vx, ρ_v .* Vx, dx)
        Vx[2:end-1] .= Vx[2:end-1] .+ ρdPdx .* dt .+ ρVxdVxdx .* dt .+ VxdρVxdx .* dt

        t += dt
        if i % divisor == 0
            #fig2 = Figure(size=(600, 800))
            #ax1 = Axis(fig2[1,1], title="Pressure, time = $t")#, limits=(nothing, nothing, -0.25, 0.25))
            #ax2 = Axis(fig2[2,1], title="Velocity")#, limits=(nothing, nothing, -0.25, 0.25))
            li = lines!(ax1, xc_vec, P, label="time = $t")
            push!(linplots, li)
            lines!(ax2, xv_vec[2:end-1], Vx[2:end-1])
            #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(i).png", fig2)
            display(fig)
        end
    end
    #@infiltrate
    # +dt*divisor
    Legend(fig[2,1], linplots, string.(round.(0:dt*divisor:dt*nt, digits=8)), "Total time", nbanks=2, orientation=:horizontal, tellhight = false, tellwidth = false)
    save("../Plots/Navier-Stokes_acoustic_wave/with_realistic_parameters_upwind_v2/All_in_one.png", fig)
    display(fig)
end