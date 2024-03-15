using CairoMakie
using SodShockTube
using Infiltrator

function av_x(B)
    A = 0.5 .* (B[2:end,:] .+ B[1:end-1,:])
end

function pvrs(W_L, W_R)
    # Heat Constant Ratio
    γ = 1.4

    # Obtaining parameters from primitive variable vectors 
    rho_L, u_L, p_L = W_L[1], W_L[2], W_L[3]
    rho_R, u_R, p_R = W_R[1], W_R[2], W_R[3]
    a_L, a_R = sqrt(γ * p_L / rho_L), sqrt(γ * p_R / rho_R)

    rho_avg = 0.5 * (rho_L + rho_R)
    a_avg = 0.5 * (a_L + a_R)

    # Calculating Star Region Properties
    p_star = 0.5 * (p_L + p_R) + 0.5 * (u_L - u_R) * (rho_avg * a_avg)
    u_star = 0.5 * (u_L + u_R) + 0.5 * (p_L - p_R) * (rho_avg * a_avg)
    rho_L_star = rho_L + ( u_L - u_star ) * (rho_avg/a_avg) 
    rho_R_star = rho_R + ( u_star - u_R ) * (rho_avg/a_avg)

    return p_star, u_star, rho_L_star, rho_R_star 
end

function q_k(p, p_star)
    if p_star <= p
        return 1
    else
        return ((1 + (γ+1)/(2*γ) * (p_star/p - 1))^0.5)
    end
end 

# HLLC Flux
function U_k_star(U_k, W_k, S_k, S_star)
    wave_speed_coeff = (S_k .- W_k) ./ (S_k .- S_star)
    star_region_comp_vec_energy_term = U_k[3]/W_k[1] + (S_star - W_k[2]) * 
                                        ( S_star + W_k[3] / (W_k[1]*(S_k - W_k[2])) )
    star_region_comp_vec = Array([1, S_star, star_region_comp_vec_energy_term])

    star_region_var_vec = W_k[1] .* wave_speed_coeff .* star_region_comp_vec 
    return star_region_var_vec 
end

function F_k_star(U_k, W_k, S_k, F_k)
    D_star = Array([0, 1, S_star])
    star_region_flux = ( S_star * (S_k * U_k - F_k) + 
                        S_k * (W_k[3] + rho_L*(S_k - W_k[2])*(S_star - W_k[2])) * D_star ) /
                        (S_k - S_star)
    return star_region_flux
end

function riemann_solver_hllc(uL, uR)
    γ = 1.4
    gm = γ - 1.0

    # Obtaining parameters from primitive variable vectors 
    rho_L, u_L, p_L = uL[1], uL[2] / uL[1], ( (γ-1)*uL[3] - 0.5*uL[2]^2 / uL[1] )
    rho_R, u_R, p_R = uR[1], uR[2] / uR[1], ( (γ-1)*uR[3] - 0.5*uR[2]^2 / uR[1] )
    a_L, a_R = sqrt(γ * p_L / rho_L), sqrt(γ * p_R / rho_R)
   
    # Initialising primitive variable vectors
    W_L = Array([rho_L, u_L, p_L]) 
    W_R = Array([rho_R, u_R, p_R]) 

    # Pressure Estimate 
    p_star, u_star, rho_L_star, rho_R_star = pvrs(W_L, W_R)
    p_star = max(0, p_star)

    # Wave Speed Estimates
    S_L = u_L - a_L * q_k(p_L, p_star)
    S_R = u_R + a_R * q_k(p_R, p_star)
    S_star = ( p_R - p_L + rho_L*u_L*(S_L - u_L) - rho_R*u_R*(S_R - u_R) ) / ( rho_L*(S_L - u_L) - rho_R*(S_R - u_R) )

    # Location of the required intercell flux 
    if 0 <= S_L 
        # Flux is calculated from the LEFT region
        f_1 = rho_L * u_L
        f_2 = rho_L * u_L^2 + p_L
        f_3 = u_L * (uL[3] + p_L)

        F_L = Array([f_1, f_2, f_3])
        F = F_L
    elseif ( (S_L <= 0) & (S_star >= 0) )
        # Flux is calculated from the LEFT STAR region
        f_1 = rho_L * u_L
        f_2 = rho_L * u_L^2 + p_L
        f_3 = u_L * (uL[3] + p_L)

        F_L = Array([f_1, f_2, f_3])
        F = F_L + S_L * (U_k_star(uL, W_L, S_L, S_star) - uL)
    elseif ( (S_star <= 0) & (S_R >= 0) )
        # Flux is calculated from the RIGHT STAR region
        f_1 = rho_R * u_R
        f_2 = rho_R * u_R^2 + p_R
        f_3 = u_R * (uR[3] + p_R)

        F_R = Array([f_1, f_2, f_3])
        F = F_R + S_R * (U_k_star(uR, W_R, S_R, S_star) - uR)
    else
        # Flux is calculated from the RIGHT region
        f_1 = rho_R * u_R
        f_2 = rho_R * u_R^2 + p_R
        f_3 = u_R * (uR[3] + p_R)

        F_R = Array([f_1, f_2, f_3])
        F = F_R
    end
    return F
end

function shock_riemann_godunov()
    # Simulation parameters
    Lx = 1.0                     # Length of the domain 
    nx = 100                    # Number of nodes
    x = LinRange(0.0, Lx, nx)    # Spatial discretization
    dx = x[2] - x[1]            # Spatial step
    γ = 1.4                     # Ratio of specific heats for air
    t = 0.0                     # Initial time
    t_final = 0.14               # Final time
    dt = 0.0001                 # Time step
    nt = t_final / dt            # Number of time steps
    divisor = 280               # Plotting divisor

    # Allocations
    U = zeros(eltype(γ), 3, nx + 2)
    A = zeros(eltype(γ), 3, 3)
    F = zeros(eltype(γ), 3, nx + 1)
    F_0 = zeros(3)


    # Initial conditions
    U[1, 1:div(nx, 2)] .= 1.0           # Density left side
    U[3, 1:div(nx, 2)] .= 1.0           # Pressure left side
    U[1, div(nx, 2):end] .= 0.125     # Density right side
    U[3, div(nx, 2):end] .= 0.1       # Pressure right side

    # Analytical solution 
    problem = ShockTubeProblem(
        geometry = (0.0, 1.0, 0.5), # left edge, right edge, initial shock location
        left_state = (ρ = 1.0, u = 0.0, p = 1.0),
        right_state = (ρ = 0.125, u = 0.0, p = 0.1),
        t = 0.14, γ = γ)
    positions, regions, vals = solve(problem, x);
    e_anal = vals.p./((γ-1).*vals.ρ)      # it seems the e calculation is a bit different in the analytical code
    e = U[3, :] ./ (γ - 1.0) ./ U[1, :]
    
    linplots = []

    # Initial plotting
    fig = Figure(size=(1000, 800))
    ax1 = Axis(fig[2,1], title="Pressure", ylabel="Pressure", xlabel="Domain",)# limits=(nothing, nothing, P0-(A*3/5), P0+A))
    ax2 = Axis(fig[1,2], title="Velocity", ylabel="Velocity", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
    ax3 = Axis(fig[2,2], title="Energy", ylabel="Energy", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
    ax4 = Axis(fig[1,1], title="Density", ylabel="Density", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
    l0 = lines!(ax1, x, U[3, 2:end-1], label="time = 0")
    push!(linplots, l0)
    lines!(ax2, x, U[2, 2:end-1])
    lines!(ax3, x, e[2:end-1])
    lines!(ax4, x, U[1, 2:end-1])
    #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(0).png", fig)
    display(fig)

    counter = 0
    while t <= t_final
        for i in eachindex(F)
            i = 1
            U_0 = riemann_solver_hllc(U[:, i], U[:, i+1])

            # F(U) = A(U) * U
            A[1, :] .= [U_0[2], U_0[1], 0.0]
            A[2, :] .= [0.0, U_0[2], 1.0 ./ U_0[1]]
            A[3, :] .= [0.0, γ .* U_0[3], U_0[2]]
            F_0[:, i] .= A * U_0

            F[:, i] .= F_0
        end

        # Godunov method
        U[:, 2:end-1] .+= (dt / dx) .* (F[:, 2:end] .- F[:, 1:end-1])
        U[:, 1]   .= U[:, 2]
        U[:, end] .= U[:, end-1]

        e = U[3, :] ./ (γ - 1.0) ./ U[1, :]

        counter += 1
        t += dt
        if counter % divisor == 0
            fig2 = Figure(size=(1000, 800))
            ax1 = Axis(fig2[2,1], title="Pressure, time = $(round(t, digits=4))", ylabel="Pressure", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            ax2 = Axis(fig2[1,2], title="Velocity", ylabel="Velocity", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            ax3 = Axis(fig2[2,2], title="Energy", ylabel="Energy", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            ax4 = Axis(fig2[1,1], title="Density", ylabel="Density", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            opts = (;linewidth = 2, color = :red)
            lines!(ax4, x, vals.ρ; opts...)
            lines!(ax2, x, vals.u; opts...)
            lines!(ax1, x, vals.p; opts...)
            lines!(ax3, x, e_anal; opts...)
            li = lines!(ax1, x, U[3, 2:end-1], label="time = $t")
            push!(linplots, li)
            lines!(ax2, x, U[2, 2:end-1])
            lines!(ax3, x, e[2:end-1])
            lines!(ax4, x, U[1, 2:end-1])
            display(fig2)
            if i % nt == 0
                Legend(fig2[3,:], linplots, string.(round.(0:dt*divisor:dt*nt, digits=8)), "Total time", tellwidth = false, nbanks=Int((nt/divisor)+1))#, tellhight = false, tellwidth = false)#, orientation=:horizontal, tellhight = false, tellwidth = false)
                rowsize!(fig2.layout, 3, 40)
                #save("C:\\Users\\Nils\\Desktop\\Masterthesis_code\\Plots\\Navier-Stokes_shock_wave\\nonconservative\\Shock_upwind_vs_analytical.png", fig)
                display(fig2)
            end
        end
    end
end



