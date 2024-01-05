using CairoMakie
using Infiltrator
#using GeoParams

# Helper functions
function av_x(B)
    A = 0.5 .* (B[2:end] .+ B[1:end-1])                 # try overriding the values of B, dont put them into a new array
end

function explicit_solver(dPdt, P, dτxxdt, τxx, dVxdt, Vx, ρ, β, μ, dx, dt)
    dPdt .= .-(1.0 ./ β) .* diff(Vx, dims=1) ./ dx
    P .= P .+ dPdt .* dt
    dτxxdt .= 2.0 .* μ .* diff(Vx, dims=1) ./ dx
    τxx .= τxx .+ dτxxdt .* dt
    dVxdt .= .-(1.0 ./ ρ) .* diff(P, dims=1) ./ dx .+ diff(τxx, dims=1) ./ dx
    Vx[2:end-1] .= Vx[2:end-1] .+ dVxdt .* dt
end

function PT_conservative_solver()
    
end

function conservative_seismic_wave_1D()
    # Physics
    L = 1000.0                              # length of the domain
    P0 = 1.0e5                              # initial pressure
    νs = 0.25                                # Poisson's ratio
    ρs = 2700.0                              # density
    μs = 1.0e11                              # shear modulus
    Ks = (2.0 *  μs * (1.0 + νs)) / (3.0 * (1.0 - 2.0 * νs))                # model poisson ratio of the solid - bulk modulus
    βs = 1.0 / Ks                             # compressibility
    cs = sqrt(1.0 / (ρs * βs))                         # speed of sound

    # Numerics
    nx = 300                                # number of grid points
    nt = 1000                               # number of time steps
    dt = 0.01                               # time step
    dx = L / nx                             # grid spacing
    CFL = 1.0 / 16.0                      # CFL number for velocity  
    ξ  = 0.95                               # relexation factor
    x = Array(range(-L / 2.0, stop=L / 2.0, length=nx))         # grid
    t = 0.0
    err_threshold = 1.0e-3

    # Allocations
    P = zeros(nx + 2)                           # pressure
    Vx = zeros(nx + 1)                          # velocity
    Fρx = zeros(nx + 1)                         # mass flux
    ρ = ones(nx + 2) .* ρs                      # solid density
    ρ_old = zeros(nx + 2)                       # solid density at previous time step
    dρdt = zeros(nx + 2)                        # rate of change of density
    ρRes = zeros(nx)                            # residual of density
    β = ones(nx + 2) .* βs                      # solid compressibility
    μ = ones(nx + 2) .* μs                      # solid shear modulus
    c = ones(nx) .* cs                          # solid speed of sound
    εxx = zeros(nx)                             # strain
    τxx = zeros(nx)                             # stress
    τxx_old = zeros(nx)                         # stress at previous time step
    divV = zeros(nx)                            # divergence of velocity
    Mx = zeros(nx + 1)                          # momentum
    Mx_old = zeros(nx + 1)                      # momentum at previous time step
    FMxx = zeros(nx)                            # flux of momentum
    dMxdt = zeros(nx + 1)                       # rate of change of momentum
    MxRes = zeros(nx - 1)                       # residual of momentum
    dMxdξ = zeros(nx + 1)                       # change of momentum with respect to ξ
    dtρ = zeros(nx)                             # time step for density 
    dtPT = zeros(nx)                            # pseudo transient time step
    
    # Initial conditions
    @. ρ[2:end-1] = ρs - 100.0 * exp(-0.0005 * x^2.0)
    @. P[2:end-1] = P0 + 100000.0 * exp(-0.0005 * x^2.0)

    # Initial plotting
    fig = Figure()
    ax = Axis(fig[1,1], title="t = $t")
    lines!(ax, x, P[2:end-1])
    display(fig)

    for i in 1:100
        err = 1.0
        pt_counter = 0
        it_counter = 0
        ρ_old .= ρ
        τxx_old .= τxx
        Mx_old .= Mx
        dMxdξ .= 0.0
        while err > err_threshold
            pt_counter += 1
            @. β = 1.0 / Ks
            dt = minimum([dx ./ maximum(abs.(Vx)), dx .* sqrt(maximum(ρ .* β))]).*4.0 # time step size
            @. c = sqrt(β[2:end-1] / ρ[2:end-1])                                                                                               # speed of sound
            dtPT .= min(dx ./ abs.(av_x(Vx)), dx ./ c) # time step size for pressure and temperature
            @. dtρ = 1.0 / (1.0 / dt + 1.0 / (dx / c / 4.1))
            
            @. Fρx = (Vx > 0.0) * Vx * ρ[1:end-1] + (Vx < 0.0) * Vx * ρ[2:end] # mass flux in x-direction
            @. dρdt = (ρ - ρ_old) / dt                                                                           # time derivative of density
            ρRes = .-dρdt[2:end-1] .- diff(Fρx, dims=1) ./ dx       # updating residual of density
            @. ρ[2:end-1] = ρ[2:end-1] + ρRes * dtρ * CFL

            @. P = P0 + (1.0 / β) * log(abs(ρ) / ρs)
            divV .= diff(Vx, dims=1) ./ dx
            @. εxx = (2.0 / 3.0) * divV
            @. τxx = (-P[2:end-1] + τxx_old) + 2.0 * μ[2:end-1] * εxx * dt

            Mx .= av_x(ρ) .* Vx
            FMxx .= (av_x(Vx) .> 0.0) .* av_x(Vx) .* Mx[1:end-1] .+ (av_x(Vx) .< 0.0) .* av_x(Vx) .* Mx[2:end]  # upwind advective momentum flux
            @. dMxdt = (Mx - Mx_old) / dt
            MxRes .= .-dMxdt[2:end-1] .- diff(FMxx .- τxx, dims=1) ./ dx
            @. dMxdξ[2:end-1] = MxRes + dMxdξ[2:end-1] * ξ
            Mx[2:end-1] .= Mx[2:end-1] .+ dMxdξ[2:end-1] .* av_x(dtPT) .* CFL

            # Boundary conditions
            Mx[1]      = 0.0 # Mx[end]
            Mx[end]    = 0.0 # Mx[1]

            Vx .= Mx ./ av_x(ρ)

            @infiltrate
            if mod(pt_counter, 1) == 0
                err = maximum(abs.([ρRes; MxRes]))                                                                      # error for time integration concerning density and momentum
                print("PT_iter = $pt_counter, err = $err, ρRes = $(maximum(abs.(ρRes[:]))), MxRes = $(maximum(abs.(MxRes[:])))\n")
            end

            if err <= err_threshold
                print("Iterationcount = $pt_counter")
            end
        end
        t = t + dt

        if mod(i-1, 1) == 0
            fig2 = Figure()
            ax2 = Axis(fig2[1,1], title="Pressure at time = $t")
            lines!(ax2, x, P[2:end-1])
            display(fig2)
        end
    end
end