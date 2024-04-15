using CairoMakie
using SodShockTube

function av_x(B)
    A = 0.5 .* (B[2:end] .+ B[1:end-1])
end

function shock_wave1D()
    # Physics
    Lx = 1.0                                # domain
    γ = 1.4                                 # adiabatic index
    ρ0 = 1.0                                # initial density at all points
    
    # Plotting parameters
    divisor = 70000                         # plotting divisor

    # Numerics
    nx = 200                                # number of nodes in x
    dx = Lx / nx                            # step size in x
    nt = 140000                             # total number of time steps

    # Grid definition
    xv = 0:dx:Lx                            # grid vertices in x-direction
    xc = av_x(xv)                           # grid nodes in x-direction

    # Allocations
    ρ = ones(nx)
    ρ_old = zeros(nx)
    P = ones(nx)
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
    P[Int((50/100)*nx):end] .= 0.1
    ρ[Int((50/100)*nx):end] .= 0.125
    E .= P./((γ - 1.0)) + 0.5.*av_x(Vx).^2
    e = P ./ (γ - 1.0) ./ ρ

    dt = 1.0e-6                                         # time step size
    t = 0.0                                             # initial time

    xc_vec = Array(xc)
    xv_vec = Array(xv)

    # Analytical solution 
    problem = ShockTubeProblem(
        geometry = (0, Lx, 0.5), # left edge, right edge, initial shock location
        left_state = (ρ = 1.0, u = 0.0, p = 1.0),
        right_state = (ρ = 0.125, u = 0.0, p = 0.1),
        t = 0.14, γ = γ)
    positions, regions, values = solve(problem, xc);
    e_anal = values.p./((γ-1).*values.ρ)


    linplots = []

    # Initial plotting
    fig = Figure(size=(1000, 800), fontsize=20)
    ax1 = Axis(fig[2,1], title="Pressure", ylabel="Pressure", xlabel="Domain")
    ax2 = Axis(fig[1,2], title="Velocity", ylabel="Velocity", xticklabelsvisible=false, xticksvisible=false)
    ax3 = Axis(fig[2,2], title="Energy", ylabel="Energy", xlabel="Domain")
    ax4 = Axis(fig[1,1], title="Density", ylabel="Density", xticklabelsvisible=false, xticksvisible=false)
    l0 = lines!(ax1, xc_vec, P, label="time = 0")
    push!(linplots, l0)
    lines!(ax2, xv_vec, Vx)
    lines!(ax3, xc_vec, e)
    lines!(ax4, xc_vec, ρ)
    #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(0).png", fig)
    #display(fig)
    li2 = nothing
    for i = 1:nt
        ρ_old .= ρ

        # Conservation of mass
        Vxdρdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(ρ), dims=1) ./ dx
        ρdVxdt .= .-ρ[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        ρ[2:end-1] .= ρ[2:end-1] .+ Vxdρdx .* dt .+ ρdVxdt .* dt

        # Conservation of momentum
        Vxdρdt .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* (av_x(ρ) .- av_x(ρ_old))
        ρdPdx .= .-(1.0 ./ av_x(ρ)) .* diff(P, dims=1) ./ dx
        ρVxdVxdx .= .-(1.0 ./ av_x(ρ)) .* (av_x(ρ) .* Vx[2:end-1]) .* diff(av_x(Vx), dims=1) ./ dx
        VxdρVxdx .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* diff(ρ .* av_x(Vx), dims=1) ./ dx
        Vx[2:end-1] .= Vx[2:end-1] .+ ρdPdx .* dt .+ ρVxdVxdx .* dt .+ VxdρVxdx .* dt .+ Vxdρdt 

        # Conservation of energy
        VxdPdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(P), dims=1) ./ dx
        PdVxdx .= .-P[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        VxdEdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(E), dims=1) ./ dx
        EdVxdx .= .-E[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        E[2:end-1] .= E[2:end-1] .+ VxdPdx .* dt .+ PdVxdx .* dt .+ VxdEdx .* dt .+ EdVxdx .* dt

        # Equation of state for pressure of an isentropic gas
        P .= (γ .- 1.0) .* (E .- 0.5 .* ρ.* av_x(Vx).^2) 

        # Internal energy calculation
        e .= P ./ (γ - 1.0) ./ ρ

        t += dt

        # Plotting
        if i % divisor == 0
            # fig2 = Figure(size=(1000, 800), fontsize=20)
            # ax1 = Axis(fig2[2,1], title="Pressure", ylabel="Pressure", xlabel="Domain")
            # ax2 = Axis(fig2[1,2], title="Velocity", ylabel="Velocity", xlabel="Domain")
            # ax3 = Axis(fig2[2,2], title="Energy", ylabel="Energy", xlabel="Domain")
            # ax4 = Axis(fig2[1,1], title="Density", ylabel="Density", xlabel="Domain")
            opts = (;linewidth = 2, color = :red, label="analyt. @ 0.14")
            li2 = lines!(ax4, xc, values.ρ; opts...)
            lines!(ax2, xc, values.u; opts...)
            lines!(ax1, xc, values.p; opts...)
            lines!(ax3, xc, e_anal; opts...)
            li = lines!(ax1, xc_vec, P, label="time = $t")
            push!(linplots, li)
            lines!(ax2, xv_vec[2:end-1], Vx[2:end-1])
            lines!(ax3, xc_vec, e)
            lines!(ax4, xc_vec, ρ)
            #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(i).png", fig2)
            #display(fig2)
        end
    end
    # Legend plotting
    linplots[end] = li2
    text = string.(round.(0:dt*divisor:dt*nt, digits=8))
    text[end] = "0.14 (analytical)"
    Legend(fig[3,:], linplots, text, "Total time", tellwidth = false, nbanks=Int(floor((nt/divisor)+2)))
    rowsize!(fig.layout, 3, 50)
    #save("/local/home/nimeding/Ma2/Julia_scripts/3_in_one_shock_wave_vs_ana_n200_new.png", fig)
    display(fig)
end