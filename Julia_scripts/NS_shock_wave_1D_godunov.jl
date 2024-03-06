using CairoMakie
using SodShockTube

function av_x(B)
    A = 0.5 .* (B[2:end,:] .+ B[1:end-1,:])
end

function upwind(u, G, dx)
    G_p = (G[2:end-1] .- G[1:end-2])./dx
    G_m = (G[2:end-1] .+ G[3:end])./dx
    #@infiltrate
    Gflux = (G_p.*(u[2:end-1] .> 0) .+ G_m.*(u[2:end-1] .< 0))

    return Gflux
end

function upwind_center(u, G, dx)
    G_p = (G[2:end-1] .- G[1:end-2])./dx
    G_m = (G[2:end-1] .+ G[3:end])./dx
    u_c = av_x(u)
    Gflux = (G_p.*(u_c[2:end-1] .> 0) .+ G_m.*(u_c[2:end-1] .< 0))

    return Gflux
end

function riemann_solver_hllc(nx, γ, ρL, uL, EL, ρR, uR, ER, fρ, fρL, fρR, fu, fuL, fuR, fE, fEL, fER)
    gm = γ - 1.0
    Ds = Array{Float64}(undef, 3)
    Ds[1], Ds[2] = 0.0, 1.0

    for i in 1:nx+1
        #left scatter
        rhLL = ρL[i]
        uuLL = uL[i]
        eeLL = EL[i]
        pL = gm * (EL - 0.5 * ρL * uL^2)
        aL = sqrt(abs(γ*pL/ρL))

        # right scatter
        rhRR = ρR[i]
        uuRR = uR[i]
        eeRR = ER[i]
        pR = gm * (ER - 0.5 * ρR * uR^2)
        aR = sqrt(abs(γ*pR/ρR))
        
        # compute SL and sR
        SL = min(uL, uR) - max(aL, aR)
        SR = max(uL, uR) + max(aL, aR)

        # compute compound speed
        SP = (pR - pL + ρL*uL*(SL - uL) - ρR*uR*(SR - uR)) / (ρL*(SL - uL) - ρR*(SR - uR))

        # compute compound pressure
        PLR = 0.5 * (pL + pR + ρL * (SL - uL) * (SP - uL) + ρR * (SR - uR) * (SP - uR))

        Ds[3] = SP

        if SL >= 0.0
            fρ[i] = fρL[i]
            fu[i] = fuL[i]
            fE[i] = fEL[i]
        elseif SR <= 0.0
            fρ[i] = fρR[i]
            fu[i] = fuR[i]
            fE[i] = fER[i]
        elseif (SP >= 0.0) & (SL <= 0.0)
            fρ[i] = (SP * (SL*ρL[i] - fρL[i]) + SL * PLR * Ds[1]) / (SL - SP)
            fu[i] = (SP * (SL*uL[i] - fuL[i]) + SL * PLR * Ds[2]) / (SL - SP)
            fE[i] = (SP * (SL*EL[i] - fEL[i]) + SL * PLR * Ds[3]) / (SL - SP)
        elseif (SP <= 0.0) & (SR >= 0.0)
            fρ[i] = (SP * (SR*ρR[i] - fρR[i]) + SR * PLR * Ds[1]) / (SR - SP)
            fu[i] = (SP * (SR*uR[i] - fuR[i]) + SR * PLR * Ds[2]) / (SR - SP)
            fE[i] = (SP * (SR*ER[i] - fER[i]) + SR * PLR * Ds[3]) / (SR - SP)
        end
    end
    return fρ, fu, fE
end

function godunov_method(nx, γ, u, f)
    # Initialize left and right states and fluxes
    uL = zeros(nx+1, 3)
    uR = zeros(nx+1, 3)
    fL = zeros(nx+1, 3)
    fR = zeros(nx+1, 3)

    # Main loop
    for t in 0:dt:T
        # Compute left and right states and fluxes
        for i in 1:nx+1
            uL[i, :] = u[i, :]
            uR[i, :] = u[i+1, :]
            fL[i, :] = f(uL[i, :])
            fR[i, :] = f(uR[i, :])
        end

        # Call HLLC Riemann solver
        riemann_solver_hllc(nx, γ, uL, uR, f, fL, fR)

        # Update solution
        for i in 2:nx
            u[i, :] = u[i, :] - dt/dx * (f[i+1, :] - f[i, :])
        end
    end

    return u
end
# taken quantities from center nodes to compute flux for face nodes
function godunov_flux_faces(u, dx)
    f = (u[1:end-1] .- u[2:end]) ./ dx
    return f
end

# taken quantities from face nodes to compute flux for center nodes
function godunov_flux_center(u, dx)
    f = (u[2:end-2] .- u[3:end-1]) ./ dx
    return f
end

extend_vertices(x) = [x[1]; x; x[end]];

function shock_wave1D_god()
    # Physics
    Lx = 1.0                           # domain
    γ = 1.4                                # adiabatic index/ratio of specific heats
    K = 1.0e10                             # shear modulus
    ρ0 = 1.0                          # initial density at all points
    P0 = 1.0#e5                          # initial pressure at all points
    Vx0 = 0.0                          # initial velocity in x-direction for all points
    β = 1.0 / K
    # Gaussian parameters
    A = 10.0                          # gaussian maximum amplitude
    σ = Lx * 0.04                            # standard deviation of the initial pressure distribution
    
    # Plotting parameters
    divisor = 10

    # Numerics
    nx = 100                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 1400                             # number of time steps

    # Grid definition
    xc = -(Lx - dx) / 2:dx:(Lx - dx) / 2        # grid nodes in x-direction
    xv = - Lx       / 2:dx: Lx       / 2        # grid vertices in x-direction

    # Allocations
    ρ = ones(nx) .* ρ0
    ρ_old = zeros(nx)
    P = ones(nx) .* P0
    Vx = zeros(nx + 1)
    Vx_old = zeros(nx + 1)
    Mx = zeros(nx + 1)
    E = zeros(nx)
    Vxdρdx = zeros(nx - 1)
    ρdVxdt = zeros(nx - 1)
    Vxdρdt = zeros(nx - 1)
    ρdPdx = zeros(nx - 1)
    ρVxdVxdx = zeros(nx)
    VxdρVxdx = zeros(nx)
    VxdPdx = zeros(nx - 1)
    PdVxdx = zeros(nx - 1)
    VxdEdx = zeros(nx - 1)
    EdVxdx = zeros(nx - 1)

    # Initial conditions
    #P .= P0 .+ A .* exp.(.- 1.0 .* (xc ./ σ).^2.0)       # initial pressure distribution
    P[Int((50/100)*nx):end] .= 0.1
    ρ[Int((50/100)*nx):end] .= 0.125
    c = sqrt(K / ρ0)                # speed of sound
    E .= P./((γ - 1.0)) + 0.5 .* av_x(Vx).^2
    e = P ./ (γ - 1.0) ./ ρ

    dt = 1.0e-4#8 #dx / (c * 4.0)                      # time step size
    t = 0.0                                         # initial time

    xc_vec = Array(xc)
    xv_vec = Array(xv)

    # Analytical solution 
    problem = ShockTubeProblem(
                geometry = (-(Lx - dx) / 2, (Lx - dx) / 2, 0.0), # left edge, right edge, initial shock location
                left_state = (ρ = 1.0, u = 0.0, p = 1.0),
                right_state = (ρ = 0.125, u = 0.0, p = 0.1),
                t = 0.14, γ = γ)
    positions, regions, values = solve(problem, xc);
    e_anal = values.p./((γ-1).*values.ρ)      # it seems the e calculation is a bit different in the analytical code

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
        ρ_old .= ρ
        Vx_old .= Vx
        
        Vxdρdx .= .-Vx[2:end-1] .* godunov_flux_faces(ρ, dx)
        ρdVxdt .= .-av_x(ρ) .* godunov_flux_faces(av_x(Vx), dx)
        ρ[2:end-1] .= ρ[2:end-1] .+ ((Vxdρdx[1:end-1] .+ ρdVxdt[1:end-1]) .- (Vxdρdx[2:end] .+ ρdVxdt[2:end])) .* dt

        P .= (γ .- 1.0) .* (E .- 0.5 .* ρ.* av_x(Vx).^2)

        dρ = av_x(ρ .- ρ_old)
        ρ_v = extend_vertices(av_x(ρ))
        Vxdρdt .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* (dρ ./ dt)
        ρdPdx .= .-(1.0 ./ av_x(ρ)) .* diff(P, dims=1) ./ dx
        ρVxdVxdx .= .-(1.0 ./ ρ) .* (ρ .* av_x(Vx)) .* godunov_flux_faces(Vx, dx)
        VxdρVxdx .= .-(1.0 ./ ρ) .* av_x(Vx) .* godunov_flux_faces(ρ_v .* Vx, dx)
        Vx[2:end-1] .= Vx[2:end-1] .+ ((ρVxdVxdx[1:end-1] .+ VxdρVxdx[1:end-1]) .- (ρVxdVxdx[2:end] .+ VxdρVxdx[2:end])) .* dt .+ ρdPdx .* dt .+ Vxdρdt .* dt
        
        VxdPdx .= .-Vx[2:end-1] .* godunov_flux_faces(P, dx) # hier
        PdVxdx .= .-av_x(P) .* godunov_flux_faces(av_x(Vx), dx)
        VxdEdx .= .-Vx[2:end-1] .* godunov_flux_faces(E, dx)
        EdVxdx .= .-av_x(E) .* godunov_flux_faces(av_x(Vx), dx)
        E[2:end-1] .= E[2:end-1] .+ ((VxdPdx[1:end-1] .+ PdVxdx[1:end-1] .+ VxdEdx[1:end-1] .+ EdVxdx[1:end-1]) .- (VxdPdx[2:end] .+ PdVxdx[2:end] .+ VxdEdx[2:end] .+ EdVxdx[2:end])) .* dt

        #ρ, Vx, E = godunov(Vx, ρ, P, E, γ, β, dt, dx)
        e .= P ./ (γ - 1.0) ./ ρ

        t += dt
        if i % divisor == 0
            fig2 = Figure(size=(1000, 800))
            ax1 = Axis(fig2[2,1], title="Pressure, time = $(round(t, digits=4))", ylabel="Pressure", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            ax2 = Axis(fig2[1,2], title="Velocity", ylabel="Velocity", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            ax3 = Axis(fig2[2,2], title="Energy", ylabel="Energy", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            ax4 = Axis(fig2[1,1], title="Density", ylabel="Density", xlabel="Domain")#, limits=(nothing, nothing, -0.25, 0.25))
            opts = (;linewidth = 2, color = :red)
            lines!(ax4, xc, values.ρ; opts...)
            lines!(ax2, xc, values.u; opts...)
            lines!(ax1, xc, values.p; opts...)
            lines!(ax3, xc, e_anal; opts...)
            li = lines!(ax1, xc_vec, P, label="time = $t")
            push!(linplots, li)
            lines!(ax2, xv_vec[2:end-1], Vx[2:end-1])
            lines!(ax3, xc_vec, e)
            lines!(ax4, xc_vec, ρ)
            
            #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(i).png", fig2)
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