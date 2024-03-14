using CairoMakie
using Infiltrator

function av_x(B)
    A = 0.5 .* (B[2:end,:] .+ B[1:end-1,:])
end

function shock_wave1D()
    # Physics
    Lx = 1.0                           # domain
    γ = 1.4
    K = 1.42e5                             # shear modulus
    ρ0 = 1.0                          # initial density at all points
    P0 = 1.0#e5                          # initial pressure at all points
    Vx0 = 0.0                          # initial velocity in x-direction for all points
    # Gaussian parameters
    A = 10.0                          # gaussian maximum amplitude
    σ = Lx * 0.04                            # standard deviation of the initial pressure distribution
    
    # Plotting parameters
    divisor = 280#400000 

    # Numerics
    nx = 100                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 14000#2000000                             # number of time steps

    # Grid definition
    xc = -(Lx - dx) / 2:dx:(Lx - dx) / 2        # grid nodes in x-direction
    xv = - Lx       / 2:dx: Lx       / 2        # grid vertices in x-direction

    # Allocations
    ρ = ones(nx) #.* ρ0
    ρ_old = zeros(nx)
    P = ones(nx) #.* P0
    Vx = zeros(nx + 1)
    E = zeros(nx)
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
    #P[1:Int((48/100)*nx)] .= P0 .+ A
    #ρ[1:Int((48/100)*nx)] .= ρ0 .+ A
    P[Int((50/100)*nx):end] .= 0.1
    ρ[Int((50/100)*nx):end] .= 0.125
    c = sqrt(K / ρ0)                # speed of sound
    E .= P./((γ - 1.0)) + 0.5.*av_x(Vx).^2
    e = P ./ (γ - 1.0) ./ ρ
    #ρ .= P ./ c^2.0                                 # initial density distribution

    dt = 1.0e-4#9 #dx / (c * 4.0)                      # time step size
    t = 0.0                                         # initial time

    xc_vec = Array(xc)
    xv_vec = Array(xv)

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
        #c = sqrt(K / maximum(ρ))                # speed of sound
        #@show c
        ρ_old .= ρ

        Vxdρdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(ρ), dims=1) ./ dx
        ρdVxdt .= .-ρ[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        ρ[2:end-1] .= ρ[2:end-1] .+ Vxdρdx .* dt .+ ρdVxdt .* dt

        P .= (γ .- 1.0) .* (E .- 0.5 .* ρ.* av_x(Vx).^2) #c.^2.0 .* ρ

        Vxdρdt .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* (av_x(ρ) .- av_x(ρ_old))
        ρdPdx .= .-(1.0 ./ av_x(ρ)) .* diff(P, dims=1) ./ dx
        ρVxdVxdx .= .-(1.0 ./ av_x(ρ)) .* (av_x(ρ) .* Vx[2:end-1]) .* diff(av_x(Vx), dims=1) ./ dx
        VxdρVxdx .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* diff(ρ .* av_x(Vx), dims=1) ./ dx
        Vx[2:end-1] .= Vx[2:end-1] .+ ρdPdx .* dt .+ ρVxdVxdx .* dt .+ VxdρVxdx .* dt .+ Vxdρdt 

        VxdPdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(P), dims=1) ./ dx
        PdVxdx .= .-P[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        VxdEdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(E), dims=1) ./ dx
        EdVxdx .= .-E[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        E[2:end-1] .= E[2:end-1] .+ VxdPdx .* dt .+ PdVxdx .* dt .+ VxdEdx .* dt .+ EdVxdx .* dt

        e .= P ./ (γ - 1.0) ./ ρ

        t += dt
        if i % divisor == 0
            fig2 = Figure(size=(1000, 800))
            ax1 = Axis(fig2[2,1], title="Pressure, time = $t")#, limits=(nothing, nothing, -0.25, 0.25))
            ax2 = Axis(fig2[1,2], title="Velocity")#, limits=(nothing, nothing, -0.25, 0.25))
            ax3 = Axis(fig2[2,2], title="Energy", ylabel="Energy", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            ax4 = Axis(fig2[1,1], title="Density", ylabel="Density", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            li = lines!(ax1, xc_vec, P, label="time = $t")
            push!(linplots, li)
            lines!(ax2, xv_vec[2:end-1], Vx[2:end-1])
            lines!(ax3, xc_vec, e)
            lines!(ax4, xc_vec, ρ)
            #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(i).png", fig2)
            display(fig2)
        end
    end
    #@infiltrate
    #Legend(fig[3,:], linplots, string.(round.(0:dt*divisor:dt*nt, digits=8)), "Total time", tellwidth = false, nbanks=Int(floor((nt/divisor)+1)))#, tellhight = false, tellwidth = false)#, orientation=:horizontal, tellhight = false, tellwidth = false)
    #rowsize!(fig.layout, 3, 40)
    #save("../Plots/Navier-Stokes_shock_wave/all_terms_sod_shock_setup/All_in_one.png", fig)
    #display(fig)
end