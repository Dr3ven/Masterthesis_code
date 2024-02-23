using CairoMakie

function av_x(B)
    A = 0.5 .* (B[2:end,:] .+ B[1:end-1,:])
end

function ac_wave1D()
    # Physics
    Lx = 1.0                           # domain
    K = 1.0                             # shear modulus
    ρ0 = 1.0                          # initial density at all points
    P0 = 1.0                          # initial pressure at all points
    Vx0 = 0.0                          # initial velocity in x-direction for all points
    # Gaussian parameters
    A = 10.0                          # gaussian maximum amplitude
    σ = Lx * 0.04                            # standard deviation of the initial pressure distribution

    # Numerics
    nx = 200                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 1000                             # number of time steps

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
    c = 1.800#sqrt(K / maximum(ρ))                # speed of sound
    P .= P0 .+ A .* exp.(.- 1.0 .* (xc ./ σ).^2.0)       # initial pressure distribution
    ρ .= P ./ c^2.0                                 # initial density distribution
    #P .= c.^2.0 .* ρ                                 # initial pressure distribution
    #ρ .= P ./ c^2.0                                 # initial density distribution

    dt = 0.0001 #dx / (c * 4.0)                      # time step size
    t = 0.0                                         # initial time

    xc_vec = Array(xc)
    xv_vec = Array(xv)

    # Initial plotting
    fig = Figure(size=(600, 800))
    ax1 = Axis(fig[1,1], title="Pressure, t = $t")
    ax2 = Axis(fig[2,1], title="Velocity")
    scatter!(ax1, xc_vec, P)
    scatter!(ax2, xv_vec, Vx)
    save("../Plots/Navier-Stokes_acoustic_wave/no_advection/$(0).png", fig)
    display(fig)

    for i = 1:nt
        #c = sqrt(K / maximum(ρ))                # speed of sound

        Vxdρdx .= 0.0#.-av_x(Vx[2:end-1]) .* diff(av_x(ρ), dims=1) ./ dx
        ρdVxdt .= .-ρ[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        ρ[2:end-1] .= ρ[2:end-1] .- (Vxdρdx .+ ρdVxdt) .* dt

        P .= c.^2.0 .* ρ

        ρdPdx .= .-(1.0 ./ av_x(ρ)) .* diff(P, dims=1) ./ dx
        ρVxdVxdx .= 0#.-(1.0 ./ av_x(ρ)) .* (av_x(ρ) .* Vx[2:end-1]) .* diff(av_x(Vx), dims=1) ./ dx
        VxdρVxdx .= 0#.-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* diff(ρ .* av_x(Vx), dims=1) ./ dx
        Vx[2:end-1] .= Vx[2:end-1] .- (ρdPdx .+ ρVxdVxdx .+ VxdρVxdx) .* dt

        t += dt
        if i % 10 == 0
            fig2 = Figure(size=(600, 800))
            ax1 = Axis(fig2[1,1], title="Pressure, time = $t")#, limits=(nothing, nothing, -0.25, 0.25))
            ax2 = Axis(fig2[2,1], title="Velocity")#, limits=(nothing, nothing, -0.25, 0.25))
            lines!(ax1, xc_vec, P)
            lines!(ax2, xv_vec, Vx)
            save("../Plots/Navier-Stokes_acoustic_wave/no_advection/$(i).png", fig2)
            display(fig2)
        end
    end
end