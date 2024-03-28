using CairoMakie

function av_x(B)
    A = 0.5 .* (B[2:end] .+ B[1:end-1])
end

function ac_wave1D()
    # Physics
    Lx = 1.0                           # domain
    K = 1.0e10                             # shear modulus
    ρ0 = 3000.0                          # initial density at all points
    P0 = 1.0e6                          # initial pressure at all points
    Vx0 = 0.0                          # initial velocity in x-direction for all points
    # Gaussian parameters
    A = 1.0e7                          # gaussian maximum amplitude
    σ = Lx * 0.04                            # standard deviation of the initial pressure distribution
    
    # Plotting parameters
    divisor = 3000 

    # Numerics
    nx = 200                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 10000                             # number of time steps

    # Grid definition
    xv = 0.0:dx:Lx                      # grid nodes in x-direction
    xc = av_x(xv)        # grid vertices in x-direction

    # Allocations
    ρ = ones(nx) .* ρ0
    ρ_old = zeros(nx)
    P = ones(nx) #.* P0
    Vx = zeros(nx + 1)
    Vxdρdx = zeros(nx - 2)
    ρdVxdt = zeros(nx - 2)
    Vxdρdt = zeros(nx - 1)
    ρdPdx = zeros(nx - 1)
    ρVxdVxdx = zeros(nx - 1)
    VxdρVxdx = zeros(nx - 1)
    e = zeros(nx)

    # Initial conditions
    dP = A .* exp.(.- 1.0 .* ((xc .- 0.5 * Lx)./ σ).^2.0)       # initial pressure distribution
    P .= P0 .+ dP
    c = sqrt(K / ρ0)                # speed of sound
    dρ = dP ./ c^2.0                                 # initial density distribution
    ρ .+= dρ

    dt = 1.0e-8 #dx / (c * 4.0)                      # time step size
    t = 0.0                                         # initial time

    xc_vec = Vector(xc)
    xv_vec = Vector(xv)

    linplots = []

    # Initial plotting
    fig = Figure(size=(1000, 800), fontsize=20)
    #ax1 = Axis(fig[1,1], title="Pressure", ylabel="Pressure", xlabel="Domain",)# limits=(nothing, nothing, P0-(A*3/5), P0+A))
    #ax2 = Axis(fig[3,1], title="Velocity", ylabel="Velocity", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
    ax1 = Axis(fig[1,1], title="Density" , ylabel="Density", xticklabelsvisible=false, xticksvisible=false)
    ax2 = Axis(fig[1,2], title="Velocity", ylabel="Velocity", xticklabelsvisible=false, xticksvisible=false)
    ax3 = Axis(fig[2,1], title="Pressure", xlabel="Domain", ylabel="Pressure")
    ax4 = Axis(fig[2,2], title="Energy"  , xlabel="Domain", ylabel="Energy")
    l0 = lines!(ax1, xc_vec, ρ)
    push!(linplots, l0)
    lines!(ax2, xv_vec, Vx)
    lines!(ax3, xc_vec, P)
    lines!(ax4, xc_vec, e)
    #save("../Plots/Navier-Stokes_acoustic_wave/with_advection_more_realistic_params/$(0).png", fig)
    display(fig)

    for i = 1:nt
        ρ_old .= ρ
        #c = sqrt(K / maximum(ρ))                # speed of sound
        #@show c

        Vxdρdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(ρ), dims=1) ./ dx
        ρdVxdt .= .-ρ[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        dρ = Vxdρdx .* dt .+ ρdVxdt .* dt
        ρ[2:end-1] .= ρ[2:end-1] .+ dρ

        dP = c.^2.0 .* dρ
        P[2:end-1] .= P[2:end-1] .+ dP

        Vxdρdt .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* (av_x(ρ) .- av_x(ρ_old))
        ρdPdx .= .-(1.0 ./ av_x(ρ)) .* diff(P, dims=1) ./ dx
        ρVxdVxdx .= .-(1.0 ./ av_x(ρ)) .* (av_x(ρ) .* Vx[2:end-1]) .* diff(av_x(Vx), dims=1) ./ dx
        VxdρVxdx .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* diff(ρ .* av_x(Vx), dims=1) ./ dx
        Vx[2:end-1] .= Vx[2:end-1] .+ ρdPdx .* dt .+ ρVxdVxdx .* dt .+ VxdρVxdx .* dt .+ Vxdρdt .* dt

        t += dt
        if i % divisor == 0
            #fig2 = Figure(size=(1000, 800))
            #ax1 = Axis(fig2[1,1], title="Pressure, time = $t")#, limits=(nothing, nothing, -0.25, 0.25))
            #ax2 = Axis(fig2[2,1], title="Velocity")#, limits=(nothing, nothing, -0.25, 0.25))
            #li = lines!(ax1, xc_vec, P, label="time = $t")
            #push!(linplots, li)
            # ax1 = Axis(fig2[1,1], title="Density, time = $t")
            # ax2 = Axis(fig2[1,2], title="Velocity")
            # ax3 = Axis(fig2[2,1], title="Pressure")#, limits=(nothing, nothing, P0, P_max))
            # ax4 = Axis(fig2[2,2], title="Energy")
            l0 = lines!(ax1, xc_vec, ρ)
            push!(linplots, l0)
            lines!(ax2, xv_vec, Vx)
            lines!(ax3, xc_vec, P)
            lines!(ax4, xc_vec, e)
            #lines!(ax2, xv_vec[2:end-1], Vx[2:end-1])
            #save("../Plots/Navier-Stokes_acoustic_wave/with_advection_more_realistic_params/$(i).png", fig2)
            #display(fig2)
        end
    end
    Legend(fig[3,:], linplots, string.(0:dt*divisor:dt*nt), "Total time", nbanks=Int(floor((nt/divisor)+1)), tellhight = false, tellwidth = false)
    rowsize!(fig.layout, 3, 40)
    #save("C:\\Users\\Nils\\Desktop\\Masterarbeit_plots\\Nils_noncons_Euler-equations_acoustic_wave\\cent_diff\\4_in_one_acoustic.png", fig)
    display(fig)
end