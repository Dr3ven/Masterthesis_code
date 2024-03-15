# using CairoMakie

# function check_upwind(u)
#     s_u = zeros(length(u))
#     for i in 2:length(u)-1
#         if u[i] > 0.0
#             s_u[i] = -1
#         elseif u[i] < 0.0
#             s_u[i] = 1
#         else
#             s_u[i] = 0
#         end    end
#     return s_u
# end

# function sod_v1()
#     L = 1.0
#     γ = 1.4
#     N = 100
#     dx = L / N
#     dt = 0.0001
#     t = 0.0
#     t_final = 0.14
#     x = Array(range(0.0+dx, stop=1.0, length=N))


#     # Initial conditions for the Sod Shock Tube
#     ρ = ones(N)
#     u = zeros(N)
#     p = ones(N)
#     ρ_flux = zeros(N - 1)
#     m = zeros(N) 
#     m_flux = zeros(N - 1)
#     ε = zeros(N)
#     E_total = zeros(N)
#     E_flux = zeros(N - 1)

#     ρ[div(N, 3):end] .= 0.125
#     p[div(N, 3):end] .= 0.1

#     i = 0
#     while t_final > t
#         i += 1 

#         s_u = check_upwind(u)

#         m .= ρ .* u
#         ε .= p ./ (γ .- 1.0) .* ρ
#         E_total .= ρ .* ε .+ 0.5 .* m .* u  # Total energy
#         ρ_flux[2:end-1] .= m[2:end-2] .* (u[2:end-2] .> 0.0) .+ m[3:end-1] .* (u[2:end-2] .<= 0.0)
#         m_flux[2:end-1] .= (m[2:end-2].^2 ./ ρ[2:end-2]) .* (u[2:end-2] .> 0.0) .+ (m[3:end-1].^2 ./ ρ[3:end-1]) .* (u[2:end-2] .<= 0.0)
#         E_flux[2:end-1] .= ((m[2:end-2] ./ ρ[2:end-2]) .* (E_total[2:end-2] .+ p[2:end-2])) .* (u[2:end-2] .> 0.0) .+ ((m[3:end-1] ./ ρ[3:end-1]) .* (E_total[3:end-1] .+ p[3:end-1])) .* (u[2:end-2] .<= 0.0)

#         ρ[2:end-1] .-= s_u .* dt/dx .* (ρ_flux[2:end-1] .- ρ_flux[2+:end-1])
#         m[2:end-1] .-= s_u .* dt/dx .* (m_flux[2:end-1] .- m_flux[1:end-1]) .- dt/(2.0 .* dx) .* (p[3:end] .- p[1:end-2])
#         E_total[2:end-1] .-= s_u .* dt/dx .* (E_flux[2:end] .- E_flux[1:end-1])

#         p .= (γ .- 1.0) .* (ρ .* E_total .- 0.5 .* m .* u)
#         u .= m ./ ρ

#         t += dt
#         if mod(i, 20) == 0.0
#             fig = Figure(size=(800,800))
#             ax1 = Axis(fig[1,1], title="Density, t=$t")
#             ax2 = Axis(fig[1,2], title="Velocity, t=$t")
#             ax3 = Axis(fig[2,1], title="Pressure, t=$t")
#             ax4 = Axis(fig[2,2], title="Energy, t=$t")
#             scatter!(ax1, x, ρ)
#             scatter!(ax2, x, m ./ ρ)
#             scatter!(ax3, x, p)
#             scatter!(ax4, x, E_total)
#             display(fig)
#         end
#     end
# end

#=using CairoMakie

function run()
    # Simulation parameters
    γ = 1.4  # Ratio of specific heats for air
    g = 9.81 # Gravitational acceleration
    N = 501   # Number of nodes
    Δx = 0.01
    Δt = 0.0001
    t_final = 0.2

    x = Array(range(0.0, stop=1.0, length=N))
    # Initial conditions for the Sod Shock Tube
    ρ = ones(N + 2)
    Fρ = zeros(N)
    u = zeros(N + 1)
    p = ones(N + 2)
    Fp1 = zeros(N)
    Fp2 = zeros(N)
    m = zeros(N + 1)
    Fm1 = zeros(N)
    Fm2 = zeros(N)
    E = ones(N + 2) .* 2.5

    ρ[div(N, 2):end] .= 0.125
    p[div(N, 2):end] .= 0.1

    # Main time-stepping loop
    t = 0.0

    fig = Figure(size = (800, 600))
    ax1 = Axis(fig[1,1], title="Density", xlabel="x", ylabel="ρ")
    ax2 = Axis(fig[1,2], title="Velocity", xlabel="x", ylabel="u")
    ax3 = Axis(fig[2,1], title="Pressure", xlabel="x", ylabel="p")
    ax4 = Axis(fig[2,2], title="Energy", xlabel="x", ylabel="E")

    lines!(ax1, x, ρ)
    lines!(ax2, x, m ./ ρ)
    lines!(ax3, x, p)
    lines!(ax4, x, E)
    display(fig)

    while t < t_final
        Fρ .= (u .> 0.0) .* diff(ρ[1:end-1] .* u) ./ Δx .+ (u .< 0.0) .* diff(ρ[2:end] .* u) ./ Δx # density flux (upwind)
        @. ρ[2:end-1] = ρ[2:end-1] + Fρ * Δt                                            # density update

        #@. p = p[1] * (γ - 1) * ρ * e                           # pressure
        Fp1 .= (u .> 0.0) .* diff(p[1:end-1] .* u) ./ Δx .+ (u .< 0.0) .* diff(p[2:end] .* u) ./ Δx # pressure flux (upwind)
        Fp2 .= Fp1 .+ γ .* (diff(ρ[1:end-1] .* u) ./ Δx) .* (u > 0.0) .+ γ .* (diff(ρ[2:end] .* u) ./ Δx) .* (u .< 0.0) # pressure flux + source terms
        @. p = p + Fp2 * Δt                                           # pressure update

        @. m = (u > 0.0) * ρ[1:end-1] * u + (u < 0.0) * ρ[2:end] * u # momentum (upwind
        @. Fm1 = (u > 0.0) * m[1:end-1] * u[1:end-1] + (u < 0.0) * m[2:end] * u[2:end] # momentum flux (upwind)
        @. Fm2 = Fm1 + p - ρ * g                                     # momentum flux + source terms
        @. m = m + Fm2 * Δt                                           # momentum update

        @. e = p / ((γ - 1) * ρ)                          # specific internal energy
        @. E = ρ * e + 0.5 * m * u                                   # total energy
        #sod_shock_tube_upwind!(ρ, u, p, m, E, Δt, Δx, γ)
        t += Δt
    end

    # Display the final results
    println("Density: ", ρ)
    println("Velocity: ", u)
    println("Pressure: ", p)

    fig = Figure(size = (800, 600))
    ax1 = Axis(fig[1,1], title="Density", xlabel="x", ylabel="ρ")
    ax2 = Axis(fig[1,2], title="Velocity", xlabel="x", ylabel="u")
    ax3 = Axis(fig[2,1], title="Pressure", xlabel="x", ylabel="p")
    ax4 = Axis(fig[2,2], title="Energy", xlabel="x", ylabel="E")

    lines!(ax1, x, ρ)
    lines!(ax2, x, m ./ ρ)
    lines!(ax3, x, p)
    lines!(ax4, x, E)
    #save("./sod_shock_tube.png", fig)
    display(fig)
end



    # Time stepping loop
    while t_final > t
        # Compute fluxes
        m = ρ .* u
        E_total = p/(γ-1) + 0.5 * m .* u
        ρ_flux .= (u .< 0.0) m[1:end-1] .+ (u .> 0.0) .* m[2:end]
        m_flux = m .* u + p
        E_flux = u .* (E_total + p)

        # Update variables
        ρ[2:end-1] -= dt./dx .* (ρ_flux[2:end] .- ρ_flux[1:end-1])
        m[2:end-1] -= dt./dx .* (m_flux[2:end] .- m_flux[1:end-1])
        E_total[2:end-1] -= dt./dx .* (E_flux[2:end] .- E_flux[1:end-1])

        # Compute derived variables
        u .= m ./ ρ
        p .= (γ.-1) .* (E_total .- 0.5 .* m .* u)

        t += dt
    end

    # Plot density
    fig = Figure()
    ax = Axis(fig[1,1], title="Density")
    scatter!(ax, x, ρ, label="Density", xlabel="x", ylabel="Density")
end=#


#=using CairoMakie
using Infiltrator


function av_x(B)
    A = 0.5 .* (B[2:end] .+ B[1:end-1])
end

function wave_1D()
    # Physics
    Lx = 10.0                           # domain
    P0 = 0.0                          # initial pressure at all points
    γ  = 1.4                            # adiabatic index
    c = 3000.0

    # Numerics
    nx = 100                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 100000                             # number of time steps

    # Grid definition
    xc = -(Lx - dx) / 2:dx:(Lx - dx) / 2        # grid nodes in x-direction
    xv = - Lx       / 2:dx: Lx       / 2        # grid vertices in x-direction

    # Allocations
    E = zeros(nx)
    #c = ones(nx)
    ρ = ones(nx)
    Mx = zeros(nx + 1)
    P = ones(nx) 
    Vx = zeros(nx + 1)
    dEudx = zeros(nx)
    dPudx = zeros(nx)
    dEdt = zeros(nx)
    dVxdx = zeros(nx - 1)
    dρdt = zeros(nx)
    dMxdt = zeros(nx - 1)

    # Initial conditions
    P .= P0 .+ exp.(.- 1.0 .* xc.^2.0)        # initial pressure distribution
    #E .= P ./ (γ .- 1.0) + (1.0 ./ 2.0) .* av_x(Vx).^2

    P_max = maximum(P)
    t = 0.0                              # initial time
    dt = 0.0001

    xc_vec = Vector(xc)
    xv_vec = Vector(xv)

    # Initial plotting
    fig = Figure()
    ax = Axis(fig[1,1], title="t = $t")
    lines!(ax, xc_vec, P)
    display(fig)

    for i = 1:nt
        #c .= sqrt.(P ./ ρ)                            # speed of sound

        dρdt[2:end-1] .= -(ρ[2:end-1] .* diff(Vx[2:end-1], dims=1)) ./ dx
        ρ[2:end-1] .= ρ[2:end-1] .+ dρdt[2:end-1] .* dt
        
        dVxdx[2:end-1] .= diff(ρ[2:end-1] .* (av_x(Vx[2:end-1]) .> 0.0) .* Vx[2:end-2].^2.0, dims=1) ./ dx .+ diff(ρ[2:end-1] .* (av_x(Vx[2:end-1]) .< 0.0) .* Vx[3:end-1].^2.0, dims=1) ./ dx
        dMxdt[2:end-1] .= -dVxdx[2:end-1] .- diff(P[2:end-1], dims=1) ./ dx
        Mx[2:end-1] .= Mx[2:end-1] .+ dMxdt .* dt

        #dEudx[2:end-1] .= diff((Vx[2:end-1] .> 0.0) .* Vx[2:end-1] .* av_x(E), dims=1) ./ dx .+ diff((Vx[2:end-1] .< 0.0) .* Vx[2:end-1] .* av_x(E), dims=1) ./ dx
        #dPudx[2:end-1] .= diff((Vx[2:end-1] .> 0.0) .* Vx[2:end-1] .* av_x(P), dims=1) ./ dx .+ diff((Vx[2:end-1] .< 0.0) .* Vx[2:end-1] .* av_x(P), dims=1) ./ dx
        #dEdt[2:end-1] .= -dEudx[2:end-1] .- dPudx[2:end-1]
        #E[2:end-1] .= E[2:end-1] .+ dEdt[2:end-1] .* dt
        
        Vx[2:end-1] .= Mx[2:end-1] ./ av_x(ρ)
        P .= c.^2.0 .* ρ

        #Vx[1] = Vx[2]
        #Vx[end] = Vx[end-1]

        t += dt
        if i % 100 == 0
            fig2 = Figure()
            ax1 = Axis(fig2[1,1], title="Pressure, time = $t")#, limits=(nothing, nothing, P0, P_max))
            ax2 = Axis(fig2[1,2], title="Velocity")
            ax3 = Axis(fig2[2,1], title="Energy")
            ax4 = Axis(fig2[2,2], title="Density")
            #ylims!(ax, -1.2, 1.2)
            lines!(ax1, xc_vec, P)
            lines!(ax2, xv_vec, Vx)
            lines!(ax3, xc_vec, E)
            lines!(ax4, xc_vec, ρ)
            display(fig2)
        end
    end
end=#

using CairoMakie
using SodShockTube

# Average function
function av_x(B)
    A = 0.5 .* (B[2:end,:] .+ B[1:end-1,:])
end

# Riemann solver
function riemann_solver(ρL, uL, pL, EL, ρR, uR, pR, ER, γ)
    # Compute sound speeds
    cL = sqrt.(γ .* pL ./ ρL)
    cR = sqrt.(γ .* pR ./ ρR)

    # Compute wave speeds
    sL = min(minimum(av_x(uL) .- cL), minimum(av_x(uR) .- cR))
    sR = max(maximum(av_x(uL) .+ cL), maximum(av_x(uR) .+ cR))

    # Compute fluxes
    FL = ρL .* av_x(uL), ρL .* av_x(uL).^2 .+ pL, (ρL .* EL .+ pL) .* av_x(uL)
    FR = ρR .* av_x(uR), ρR .* av_x(uR).^2 .+ pR, (ρR .* ER .+ pR) .* av_x(uR)

    # Compute HLL flux
    if sL >= 0
        return FL
    elseif sR <= 0
        return FR
    else
        return ((sR .* FL .- sL .* FR) .+ sL .* sR .* (ρR .- ρL, av_x(uR .- uL), ER .- EL)) ./ (sR .- sL)
    end
end

# Godunov method
function godunov(Vx, ρ, p, E, γ, dt, dx)
    # Solve Riemann problem at boundaries
    ρL, ρR = ρ[1:end-1], ρ[2:end]
    VL, VR = Vx[1:end-1], Vx[2:end]
    pL, pR = p[1:end-1], p[2:end]
    EL, ER = E[1:end-1], E[2:end]
    ρ_star, V_star, p_star = riemann_solver(ρL, VL, pL, EL, ρR, VR, pR, ER, γ)

    # Compute fluxes
    F_ρ = ρ_star .* V_star
    F_Vx = ρ_star .* V_star.^2 + p_star
    F_E = (ρ_star .* EL + p_star) .* V_star

    # Update variables
    @show size(F_ρ)
    @show size(F_Vx)
    @show size(F_E)
    ρ_new = ρ[2:end-1] .- dt/dx .* (F_ρ[2:end] - F_ρ[1:end-1])
    Vx_new = av_x(Vx[2:end-1]).- dt/dx .* (F_Vx[2:end] - F_Vx[1:end-1])
    E_new = E[2:end-1] .- dt/dx .* (F_E[2:end] - F_E[1:end-1])

    return ρ_new, Vx_new, E_new
end

function godunov_test()
    # Constants
    γ = 1.4
    dx = 1.0 / 100
    dt = 0.0001
    t_end = 0.1
    x_c = 0:dx:1  # center points
    x_f = 0:dx:1+dx  # face points

    # Initial conditions (Sod's shock tube)
    ρ = [x < 0.5 ? 1.0 : 0.125 for x in x_c]
    Vx = zeros(length(x_f))
    E = [x < 0.5 ? 1.0 : 0.1 for x in x_c]
    P = (γ - 1.0) .* (E .- 0.5 .* ρ .* av_x(Vx).^2)

    # Time stepping
    t = 0.0
    while t < t_end
        ρ, Vx, E = godunov(Vx, ρ, P, E, γ, dt, dx)
        P .= (γ .- 1.0) .* (E .- 0.5 .* ρ .* av_x(Vx).^2)
        t += dt
    end

    # Plot results
    fig = Figure()
    ax1 = Axis(fig[1, 1], title="Density")
    lines!(ax1, x_c, ρ)
    ax2 = Axis(fig[2, 1], title="Velocity")
    lines!(ax2, x_f, Vx)
    ax3 = Axis(fig[3, 1], title="Energy")
    lines!(ax3, x_c, E)
    display(fig)
end