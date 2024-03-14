using CairoMakie

function acoustic_wave()
    # Physics
    Lx = 1.0                           # domain
    ρ = 3000.0                             # density
    β = 1.0e-10                             # compressibility
    μ = 1.0                             # shear modulus
    c = sqrt(1.0 / (β * ρ))             # speed of sound
    P0 = 1.0e6                          # initial pressure at all points
    Vx0 = 0.0                          # initial velocity in x-direction for all points
    σ = Lx * 0.04                            # standard deviation of the initial pressure distribution
    A = 10.0                          # gaussian maximum amplitude
    divisor = 200

    # Numerics
    nx = 100                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 1000                             # number of time steps
    dt = 1.0e-7 #dx / (c * 4.0)                 # time step size

    # Grid definition
    xc = -(Lx - dx) / 2:dx:(Lx - dx) / 2        # grid nodes in x-direction
    xv = - Lx       / 2:dx: Lx       / 2        # grid vertices in x-direction

    # Allocations
    P = ones(nx)
    Vx = zeros(nx + 1)
    τxx = zeros(nx)
    dPdt = zeros(nx)
    dVxdt = zeros(nx - 1)

    # Initial conditions
    P .= P0 .+ A .* exp.(.- 1.0 .* (xc ./ σ).^2.0)        # initial pressure distribution
    t = 0.0                              # initial time

    xc_vec = Vector(xc)
    xv_vec = Vector(xv)
    linplots = []

    # Initial plotting
    linplots = []

    # Initial plotting
    fig = Figure(size=(600, 800))
    ax1 = Axis(fig[1,1], title="Pressure", ylabel="Pressure", xlabel="Domain",)# limits=(nothing, nothing, P0-(A*3/5), P0+A))
    ax2 = Axis(fig[3,1], title="Velocity", ylabel="Velocity", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
    l0 = lines!(ax1, xc_vec, P, label="time = 0")
    push!(linplots, l0)
    lines!(ax2, xv_vec, Vx)
    #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(0).png", fig)
    #display(fig)

    for i = 1:nt
        t += dt
        dPdt .= .-(1.0 ./ β) .* diff(Vx, dims=1) ./ dx
        P .= P .+ dPdt .* dt
        dVxdt .= .-(1.0 ./ ρ) .* diff(P, dims=1) ./ dx
        Vx[2:end-1] .= Vx[2:end-1] .+ dVxdt .* dt
        
        if i % divisor == 0
            #fig2 = Figure()
            #ax2 = Axis(fig2[1,1], title="time = $t", limits=(nothing, nothing, -0.01, 1.1))
            li = lines!(ax1, xc_vec, P, label="time = $t")
            push!(linplots, li)
            lines!(ax2, xv_vec, Vx)
            display(fig)
        end
    end
    Legend(fig[2,1], linplots, string.(round.(0:dt*divisor:dt*nt, digits=8)), "Total time", tellwidth = false, nbanks=2, tellhight = false, orientation=:horizontal)
    rowsize!(fig.layout, 2, 60)
    #save("/home/nils/Masterthesis_code/Plots/Navier-Stokes_acoustic_wave/normal_acoustic_wave/normal_acoustic_wave.png", fig)
    display(fig)
end