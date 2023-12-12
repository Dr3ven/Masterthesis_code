using CairoMakie
using TimerOutputs

@views function av_x(B)
    B = 0.5 .* (B[2:end,:] .+ B[1:end-1,:])
end

@views function av_y(B)
    B = 0.5 .* (B[:,2:end] .+ B[:,1:end-1])
end

@views function av_xy(B)
    B = 0.25 .* (B[2:end,2:end] .+ B[1:end-1,2:end] .+ B[2:end,1:end-1] .+ B[1:end-1,1:end-1])
end

@views function elastic_wave_2D()
    # Physics
    Lx = 100                           # domain in x
    Ly = Lx                             # domain in y
    ρ = 2800.0                             # density
    β = 3.0e-6                             # compressibility
    μ = 1.0                             # shear modulus
    λ = 1.0                             # Lamé parameter
    c = sqrt(1.0 / (β * ρ))             # speed of sound
    P0 = 0.0                          # initial pressure at all points

    # Numerics
    nx = 255                             # number of nodes in x
    ny = 255                             # number of nodes in y
    dx = Lx / nx                        # step size in x
    dy = Ly / ny                        # step size in y
    nt = 5000                             # number of time steps
    dt = min(dx, dy) / (c * sqrt(8.1))# (c * 6.1)                # time step size                       

    # Grid definition
    xc = -(Lx - dx) / 2.0:dx:(Lx - dx) / 2.0        # grid nodes in x-direction
    yc = -(Ly - dy) / 2.0:dy:(Ly - dy) / 2.0        # grid nodes in x-direction
    xv =   Lx       / 2.0:dx: Lx       / 2.0        # grid vertices in x-direction 
    yv =   Ly       / 2.0:dy: Ly       / 2.0        # grid vertices in x-direction

    # Allocations
    P = zeros(Float64, nx, ny)
    divV = zeros(Float64, nx, ny)
    Vx = zeros(Float64, nx + 1, ny)
    Vy = zeros(Float64, nx, ny + 1)
    εxx = zeros(Float64, nx, ny)
    εyy = zeros(Float64, nx, ny)
    εxy = zeros(Float64, nx + 1, ny + 1)
    τxx = zeros(Float64, nx, ny)
    τyy = zeros(Float64, nx, ny)
    τxy = zeros(Float64, nx + 1, ny + 1)
    dPdt = zeros(Float64, nx, ny)
    dVxdt = zeros(Float64, nx + 1, ny)
    dVydt = zeros(Float64, nx, ny + 1)
    dτxxdt = zeros(Float64, nx, ny)
    dτyydt = zeros(Float64, nx, ny)
    dτxydt = zeros(Float64, nx + 1, ny + 1)

    # Initial conditions
    P .= P0 .+ exp.((.-0.005 .* xc.^2.0) .+ (.-0.005 .* yc'.^2.0))          # initial pressure distribution
    global t = 0.0                                                      # initial time

    xc_vec = Vector(xc)
    yc_vec = Vector(yc)

    # Initial plotting
    fig = Figure()
    #ax = Axis(fig[1,1], title="t = $t")#, limits=(nothing, nothing, nothing, 1.1))
    ax3 = Axis3(fig[1,1], title="time = $t")
    limits!(ax3, -Lx / 2.0, Lx / 2.0, -Ly / 2.0, Ly / 2.0, -0.1, 1.0)
    sur2 = surface!(ax3, xc_vec, yc_vec, P)
    #lines!(ax, xc_vec, P[:, Int(nx/2)])
    #hm = heatmap!(ax, xc_vec, yc_vec, P)
    #Colorbar(fig[1,2], hm, label="Pressure")
    fig

    reset_timer!()
    for i = 1:nt
        @timeit "divV" divV .= diff(Vx, dims=1) ./ dx .+ diff(Vy, dims=2) ./ dy
        @timeit "dPdt" dPdt .= .-(1.0 ./ β) .* divV
        @timeit "P" P .= P .+ dPdt .* dt
        @timeit "εxx" εxx .= diff(Vx, dims=1) ./ dx .- (1.0 ./ 3.0) .* divV
        @timeit "εyy" εyy .= diff(Vy, dims=2) ./ dy .- (1.0 ./ 3.0) .* divV
        @timeit "εxy" εxy[2:end-1, 2:end-1] .= 0.5 .* (diff(Vx[2:end-1, :], dims=2) ./ dy .+ diff(Vy[:, 2:end-1], dims=1) ./ dx)
        @timeit "dτxxdt" dτxxdt .= (2.0 .* μ .+ λ) .* diff(Vx, dims=1) ./ dx  .+ λ .* diff(Vy, dims=2) ./ dy
        @timeit "dτyydt" dτyydt .= (2.0 .* μ .+ λ) .* diff(Vy, dims=2) ./ dy  .+ λ .* diff(Vx, dims=1) ./ dx
        @timeit "dτxydt" dτxydt[2:end-1, 2:end-1] .= μ .* (diff(Vx[2:end-1, :], dims=2) ./ dy .+ diff(Vy[:, 2:end-1], dims=1) ./ dx)
        @timeit "τxx" τxx .= τxx .+ dτxxdt .* dt
        @timeit "τyy" τyy .= τyy .+ dτyydt .* dt
        @timeit "τxy" τxy .= τxy .+ dτxydt .* dt
        @timeit "dVxdt" dVxdt[2:end-1, :] .= .-(1.0 ./ ρ) .* (diff(P, dims=1) ./ dx .+ diff(τxx, dims=1) ./ dx .+ diff(τxy[2:end-1, :], dims=2) ./ dy)
        @timeit "dVydt" dVydt[:, 2:end-1] .= .-(1.0 ./ ρ) .* (diff(P, dims=2) ./ dy .+ diff(τyy, dims=2) ./ dy .+ diff(τxy[:, 2:end-1], dims=1) ./ dx)
        @timeit "Vx" Vx[2:end-1, :] .= Vx[2:end-1, :] .+ dVxdt[2:end-1, :] .* dt 
        @timeit "Vy" Vy[:, 2:end-1] .= Vy[:, 2:end-1] .+ dVydt[: ,2:end-1] .* dt
        
        # Boundary conditions
        Vx[:, 1] .= 0.0 #.-Vx[:, 1]
        Vx[:, end] .= 0.0 #.-Vx[:, end]

        Vy[1, :] .= 0.0 #.-Vy[1, :]
        Vy[end, :] .= 0.0 #.-Vy[end, :]
        
        global t += dt

        if i % 5 == 0
            fig2 = Figure()
            ax2 = Axis(fig2[1,1], title="time = $t")#, limits=(nothing, nothing, nothing, 1.1))
            #ax3 = Axis3(fig2[1,1], title="time = $t")
            #limits!(ax3, -Lx / 2.0, Lx / 2.0, -Ly / 2.0, Ly / 2.0, -0.1, 1.0)
            #sur2 = surface!(ax3, xc_vec, yc_vec, P)
            #ylims!(ax, -1.2, 1.2)
            #lines!(ax2, xc_vec, P[:, Int(nx/2)])
            hm2 = heatmap!(ax2, xc_vec, yc_vec, P)#, colorrange=(0.0, 1.0))
            Colorbar(fig2[1,2], hm2, label="Pressure")
            display(fig2)
        end
    end
    print_timer()
end