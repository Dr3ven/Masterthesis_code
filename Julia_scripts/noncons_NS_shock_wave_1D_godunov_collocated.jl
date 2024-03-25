using CairoMakie
using SodShockTube
using Infiltrator

function av_x(B)
    A = 0.5 .* (B[2:end,:] .+ B[1:end-1,:])
end

function riemann_solver_hllc(γ, ρLi, uLi, ELi, ρRi, uRi, ERi, fρL, fρR, fuL, fuR, fEL, fER)
    gm = γ - 1.0

    #left scatter
    ρL = ρLi
    uL = uLi
    EL = ELi
    pL = gm * (EL - 0.5 * ρL * uL^2)
    aL = sqrt(abs(γ*pL/ρL))

    # right scatter
    ρR = ρRi
    uR = uRi
    ER = ERi
    pR = gm * (ER - 0.5 * ρR * uR^2)
    aR = sqrt(abs(γ*pR/ρR))
    
    # middle regions
    aS = (1.0 / 2.0) * (aL + aR) + (1.0 / 4.0) * gm * (uR - uL)
    pS = (1.0 / 2.0) * (pL + pR) - (1.0 / 2.0) * (uR - uL) * ((1.0 / 2.0) * (ρL + ρR) * (aL + aR))
    uS = (1.0 / 2.0) * (uL + uR) - (1.0 / 2.0) * (pL - pR) / ((1.0 / 2.0) * (ρL + ρR) * (aL + aR)) # eq. 12 HLLC paper
    # uS = (1.0 / 2.0) * (uL + uR) + ((aL - aR) / gm) # eq. 10 HLLC paper
    ρSL = ρL + (uL - uS) * ((1.0 / 2.0) * (ρL + ρR)) / ((1.0 / 2.0) * (aL + aR))
    ρSR = ρR + (uS - uR) * ((1.0 / 2.0) * (ρL + ρR)) / ((1.0 / 2.0) * (aL + aR))

    # compute pressure ratios
    HSL = pS / pL
    HSR = pS / pR

    # compute qL and qR
    qL = HSL <= 1 ? 1 : sqrt(1 + (γ + 1) / (2 * γ) * (HSL - 1))
    qR = HSR <= 1 ? 1 : sqrt(1 + (γ + 1) / (2 * γ) * (HSR - 1))

    # compute SL and sR
    #SL = min(uL, uR) - max(aL, aR)
    SL = uL - aL * qL 
    #SR = max(uL, uR) + max(aL, aR)
    SR = uR - aR * qR

    # compute compound speed
    SM = (pR - pL + ρL*uL*(SL - uL) - ρR*uR*(SR - uR)) / (ρL*(SL - uL) - ρR*(SR - uR))

    # compute q1, q2, q3, q4, q5, q6
    q1 = SL * ρSL - ρSL * uS 
    q2 = SL * ρL * uS - ρL * uS^2 + pL
    q3 = SL * EL - uS * (EL + pL)

    q4 = SR * ρSR - ρSR * uS
    q5 = SR * ρSR * uS - ρSR * uS^2 + pR
    q6 = SR * ER - uS * (ER + pR)

    # use the q1-6 and SL, SM, SR to compute Uk
    ρSL = q1 / (SL - SM)
    pSL = SM * q1 - q2
    ESL = q3 + SM * pSL / (SL - SM)
    
    ρSR = q4 / (SR - SM)
    pSR = SM * q4 - q5
    ESR = q6 + SM * pSR / (SR - SM)

    # Fk is given as input parameters of this function (see fρL, fρR, fuL, fuR, fEL, fER)
    # use FSL = FL + SL * (USL - UL) and FSR = FR + SR * (USR - UR) to compute fluxes -> those are Fi+1/2 then!
    FρSL = fρL + SL * (ρL - ρSL)
    FuSL = fuL + SR * (ρR - ρSR)
    FESL = fEL + SL * (EL - ESL)

    FρSR = fρR + SR * (ρR - ρSR)
    FuSR = fuR + SR * (ρR - ρSR)
    FESR = fER + SR * (ER - ESR)

    # right supersonic flow
    if SL >= 0.0
        fρ = fρL
        fu = fuL
        fE = fEL
    # left supersonic flow
    elseif SR <= 0.0
        fρ = fρR
        fu = fuR
        fE = fER
    # subsonic flow
    elseif (SM >= 0.0) & (SL <= 0.0)
        fρ = FρSL
        fu = FuSL
        fE = FESL
    # what flow is this?
    elseif (SM <= 0.0) & (SR >= 0.0)
        fρ = FρSR
        fu = FuSR
        fE = FESR
    end
    return fρ, fu, fE
end

extend_vertices(x) = [x[1]; x; x[end]];

function shock_god_col()
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
    divisor = 280 

    # Numerics
    nx = 100                             # number of nodes in x
    dx = Lx / nx                        # step size in x
    nt = 14000                             # number of time steps

    # Grid definition
    xc = -(Lx - dx) / 2:dx:(Lx - dx) / 2        # grid nodes in x-direction
    xv = - Lx       / 2:dx: Lx       / 2        # grid vertices in x-direction

    # Allocations
    u = zeros(nx, 3)
    f = zeros(nx, 3)
    ρ = ones(nx) .* ρ0
    P = ones(nx) .* P0
    Mx = zeros(nx)
    Vx = zeros(nx)
    E = zeros(nx)
    fρ = zeros(nx-1)
    fMx= zeros(nx-1)
    fE = zeros(nx-1)

    # Initial conditions
    #P .= P0 .+ A .* exp.(.- 1.0 .* (xc ./ σ).^2.0)       # initial pressure distribution
    P[Int((50/100)*nx):end] .= 0.1
    ρ[Int((50/100)*nx):end] .= 0.125
    c = sqrt(K / ρ0)                # speed of sound
    E .= P./((γ - 1.0)) + 0.5 .* ρ .* Vx.^2
    e = P ./ (γ - 1.0) ./ ρ
    

    u[:, 1]   .= ρ
    u[:, 2]   .= Vx
    u[:, 3]   .= E

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
    lines!(ax2, xc_vec, Vx)
    lines!(ax3, xc_vec, e)
    lines!(ax4, xc_vec, ρ)
    #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(0).png", fig)
    #display(fig)

    for i = 1:nt
        # ρ_old .= ρ
        # Vx_old .= Vx
        
        # fρ1 = .-av_x(Vx) .* diff(ρ, dims=1) ./ dx
        # fρ2 = .-av_x(ρ) .* diff(Vx, dims=1) ./ dx 
        # fρ3 = fρ1 .+ fρ2
        fρ4 = ρ .* Vx

        fρR .= 0.5 * (ρ[2:end] .+ ρ[1:end-1]) .- (dt ./ dx) .* (fρ4[2:end] .- fρ4[1:end-1])
        fρL .= 0.5 * (ρ[2:end] .+ ρ[1:end-1]) .- (dt ./ dx) .* (fρ4[2:end] .- fρ4[1:end-1])
        ρ[2:end-1] .= ρ[2:end-1] .- (dt ./ dx) .* (fρ[2:end] .- fρ[1:end-1])

        # fMx1 = diff(P, dims=1) ./ dx
        # fMx2 = .-av_x(Mx) .* diff(Vx, dims=1) ./ dx
        # fMx3 = .-av_x(Vx) .* diff(Mx, dims=1) ./ dx
        # fMx4 = fMx1 .+ fMx2 .+ fMx3
        fMx5 = (Mx.^2.0 ./ ρ) + P
        fMxR .= 0.5 * (Mx[2:end] .+ Mx[1:end-1]) .- (dt ./ dx) .* (fMx5[2:end] .- fMx5[1:end-1])
        Mx[2:end-1] .= Mx[2:end-1] .- (dt ./ dx) .* (fMx[2:end] .- fMx[1:end-1])

        # fE1 = .-av_x(Vx) .* diff(P, dims=1) ./ dx
        # fE2 = .-av_x(P) .* diff(Vx, dims=1) ./ dx
        # fE3 = .-av_x(Vx) .* diff(E, dims=1) ./ dx
        # fE4 = .-av_x(E) .* diff(Vx, dims=1) ./ dx
        # fE5 = fE1 .+ fE2 .+ fE3 .+ fE4
        fE6 = (Mx ./ ρ) .* (E .+ P)
        fER .= 0.5 * (E[2:end] .+ E[1:end-1]) .- (dt ./ dx) .* (fE6[2:end] .- fE6[1:end-1])
        E[2:end-1] .= E[2:end-1] .- (dt ./ dx) .* (fE[2:end] .- fE[1:end-1])

        P .= (γ .- 1.0) .* (E .- 0.5 .* ρ .* Vx.^2.0)
        e = P ./ (γ - 1.0) ./ ρ


        #=for i in 1:nx-2
            fρi, fMxi, fEi = riemann_solver_hllc(γ, ρ[i], Vx[i], E[i], ρ[i+1], Vx[i+1], E[i+1], fρ3[i], fρ3[i+1], fMx4[i], fMx4[i+1], fE5[i], fE5[i+1])
            fρ[i] = fρi
            fMx[i] = fMxi
            fE[i]  = fEi
        end=#
        


        # ρ[2:end-1] .= 0.5 .* (ρ[2:end] .+ ρ[1:end-1]) .+ Vxdρdx .* dt .+ ρdVxdt .* dt

        # P .= (γ .- 1.0) .* (E .- 0.5 .* ρ.* av_x(Vx).^2)

        # dρ = av_x(ρ .- ρ_old)
        # ρ_v = extend_vertices(av_x(ρ))
        # Vxdρdt .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* (dρ ./ dt)
        # ρdPdx .= .-(1.0 ./ av_x(ρ)) .* diff(P, dims=1) ./ dx
        # ρVxdVxdx .= .-(1.0 ./ av_x(ρ)) .* (av_x(ρ) .* Vx[2:end-1]) .* diff(av_x(Vx), dims=1) ./ dx # hier
        # VxdρVxdx .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* upwind(Vx, ρ_v .* Vx, dx)
        # Vx[2:end-1] .= Vx[2:end-1] .+ ρdPdx .* dt .+ ρVxdVxdx .* dt .+ VxdρVxdx .* dt .+ Vxdρdt .* dt
        
        # VxdPdx .= .-av_x(Vx[2:end-1]) .* upwind_center(Vx, P, dx) # hier
        # PdVxdx .= .-P[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        # VxdEdx .= .-av_x(Vx[2:end-1]) .* upwind_center(Vx, E, dx)
        # EdVxdx .= .-E[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        # E[2:end-1] .= E[2:end-1] .+ VxdPdx .* dt .+ PdVxdx .* dt .+ VxdEdx .* dt .+ EdVxdx .* dt

        # e .= P ./ (γ - 1.0) ./ ρ

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
            lines!(ax2, xc_vec, Vx)
            lines!(ax3, xc_vec, e)
            lines!(ax4, xc_vec, ρ)
            
            #save("../Plots/Navier-Stokes_acoustic_wave/discontinous_initial_condition/$(i).png", fig2)
            display(fig2)
            if i % nt == 0
                #Legend(fig2[3,:], linplots, string.(round.(0:dt*divisor:dt*nt, digits=8)), "Total time", tellwidth = false, nbanks=Int((nt/divisor)+1))#, tellhight = false, tellwidth = false)#, orientation=:horizontal, tellhight = false, tellwidth = false)
                #rowsize!(fig2.layout, 3, 40)
                #save("C:\\Users\\Nils\\Desktop\\Masterthesis_code\\Plots\\Navier-Stokes_shock_wave\\nonconservative\\Shock_upwind_vs_analytical.png", fig)
                #display(fig2)
            end
        end
    end
end