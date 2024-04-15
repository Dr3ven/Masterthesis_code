using CairoMakie

function av_x(B)
    A = 0.5 .* (B[2:end] .+ B[1:end-1])
end

function ac_wave_1D_noadvect(nt)
    # Physics
    Lx = 1.0                                    # domain
    γ = 1.4
    K = 1.0e10                                  # shear modulus
    ρ0 = 3000.0                                 # initial density at all points
    P0 = 1.0e6                                  # initial pressure at all points
    # Gaussian parameters
    A = 1.0e7                                   # gaussian maximum amplitude
    σ = Lx * 0.04                               # standard deviation of the initial pressure distribution
    
    # Plotting parameters
    divisor = 3000 

    # Numerics
    nx = 200                                    # number of nodes in x
    dx = Lx / nx                                # step size in x
    dt = 1.0e-8                                 # time step size

    # Grid definition
    xv = 0:dx:Lx                                # grid vertices in x-direction
    xc = av_x(xv)                               # grid nodes in x-direction

    # Allocations
    ρ = ones(nx) .* ρ0
    ρ_old = zeros(nx)
    P = ones(nx)
    Vx = zeros(nx + 1)
    Vxdρdx = zeros(nx - 2)
    ρdVxdt = zeros(nx - 2)
    Vxdρdt = zeros(nx - 1)
    ρdPdx = zeros(nx - 1)
    ρVxdVxdx = zeros(nx - 1)
    VxdρVxdx = zeros(nx - 1)
    EdρVxdx = zeros(nx - 2)
    VxdPdx = zeros(nx - 2)
    PdVxdx = zeros(nx - 2)
    ρVxdEdx = zeros(nx - 2)
    E = zeros(nx)
    e = zeros(nx)

    # Initial conditions
    dP = A .* exp.(.- 1.0 .* ((xc .- 0.5*Lx) ./ σ).^2.0)        # initial pressure distribution
    P .= P0 .+ dP
    c = sqrt(K / ρ0)                                            # speed of sound
    dρ = dP ./ c^2.0                                            # initial density distribution
    ρ .+= dρ
    e = P ./ (γ - 1.0) ./ ρ


    t = 0.0                                                     # initial time

    xc_vec = Array(xc)
    xv_vec = Array(xv)

    linplots = []

    # Initial plotting
    fig = Figure(size=(1000, 800), fontsize=20)
    ax1 = Axis(fig[1,1], title="Density no advection" , ylabel="Density", xticklabelsvisible=false, xticksvisible=false)
    ax2 = Axis(fig[1,2], title="Velocity", ylabel="Velocity", xticklabelsvisible=false, xticksvisible=false)
    ax3 = Axis(fig[2,1], title="Pressure", xlabel="Domain", ylabel="Pressure")
    ax4 = Axis(fig[2,2], title="Energy"  , xlabel="Domain", ylabel="Energy")
    l0 = lines!(ax1, xc_vec, ρ)
    push!(linplots, l0)
    lines!(ax2, xv_vec, Vx)
    lines!(ax3, xc_vec, P)
    lines!(ax4, xc_vec, e)
    #save("../Plots/Navier-Stokes_acoustic_wave/with_advection_more_realistic_params/$(0).png", fig)
    #display(fig)

    for i = 1:nt
        ρ_old .= ρ

        #Vxdρdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(ρ), dims=1) ./ dx
        ρdVxdt .= .-ρ[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        dρ = ρdVxdt .* dt
        ρ[2:end-1] .= ρ[2:end-1] .+ dρ
        
        Vxdρdt .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* (av_x(ρ) .- av_x(ρ_old))
        ρdPdx .= .-(1.0 ./ av_x(ρ)) .* diff(P, dims=1) ./ dx
        ρVxdVxdx .= .-(1.0 ./ av_x(ρ)) .* (av_x(ρ) .* Vx[2:end-1]) .* diff(av_x(Vx), dims=1) ./ dx
        #VxdρVxdx .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* diff(ρ .* av_x(Vx), dims=1) ./ dx
        Vx[2:end-1] .= Vx[2:end-1] .+ ρdPdx .* dt .+ ρVxdVxdx .* dt .+ Vxdρdt .* dt

        EdρVxdx.= .-E[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        VxdPdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(P), dims=1) ./ dx
        PdVxdx .= .-P[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        #ρVxdEdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(E), dims=1) ./ dx
        E[2:end-1] .= E[2:end-1] .+ EdρVxdx .* dt .+ VxdPdx .* dt .+ PdVxdx .* dt 

        dP = c.^2.0 .* dρ
        P[2:end-1] .= P[2:end-1] .+ dP

        e = P ./ (γ - 1.0) ./ ρ

        t += dt
        if i % divisor == 0
            #fig2 = Figure(size=(1000, 800))
            #ax1 = Axis(fig2[1,1], title="Pressure, time = $t")
            #ax2 = Axis(fig2[2,1], title="Velocity")
            #li = lines!(ax1, xc_vec, P, label="time = $t")
            #push!(linplots, li)
            # ax1 = Axis(fig2[1,1], title="Density, time = $t")
            # ax2 = Axis(fig2[1,2], title="Velocity")
            # ax3 = Axis(fig2[2,1], title="Pressure")
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
    #Legend(fig[3,:], linplots, string.(0:dt*divisor:dt*nt), "Total time", nbanks=Int(floor((nt/divisor)+1)), tellhight = false, tellwidth = false)
    #rowsize!(fig.layout, 3, 40)
    #save("/home/nils/Masterthesis_code/Plots/Nils_noncons_Euler-equations_acoustic_wave/no_advection/2/4_in_one_acoustic.png", fig)
    #display(fig)
    return ρ, Vx, P, e
end


function ac_wave1D(nt)
    # Physics
    Lx = 1.0                                    # domain
    γ = 1.4
    K = 1.0e10                                  # shear modulus
    ρ0 = 3000.0                                 # initial density at all points
    P0 = 1.0e6                                  # initial pressure at all points
    Vx0 = 0.0                                   # initial velocity in x-direction for all points
    # Gaussian parameters
    A = 1.0e7                                   # gaussian maximum amplitude
    σ = Lx * 0.04                               # standard deviation of the initial pressure distribution
    
    # Plotting parameters
    divisor = 20000 

    # Numerics
    nx = 200                                    # number of nodes in x
    dx = Lx / nx                                # step size in x

    # Grid definition
    xv = 0:dx:Lx                                # grid vertices in x-direction
    xc = av_x(xv)                               # grid nodes in x-direction

    # Allocations
    ρ = ones(nx) .* ρ0
    ρ_old = zeros(nx)
    P = ones(nx)
    Vx = zeros(nx + 1)
    Vxdρdx = zeros(nx - 2)
    ρdVxdt = zeros(nx - 2)
    Vxdρdt = zeros(nx - 1)
    ρdPdx = zeros(nx - 1)
    ρVxdVxdx = zeros(nx - 1)
    VxdρVxdx = zeros(nx - 1)
    EdρVxdx = zeros(nx - 2)
    VxdPdx = zeros(nx - 2)
    PdVxdx = zeros(nx - 2)
    ρVxdEdx = zeros(nx - 2)
    E = zeros(nx)
    e = zeros(nx)

    # Initial conditions
    dP = A .* exp.(.- 1.0 .* ((xc .- 0.5*Lx)./ σ).^2.0)         # initial pressure distribution
    P .= P0 .+ dP
    c = sqrt(K / ρ0)                                            # speed of sound
    dρ = dP ./ c^2.0                                            # initial density distribution
    ρ .+= dρ
    e = P ./ (γ - 1.0) ./ ρ


    dt = 1.0e-8 #dx / (c * 4.0)                                 # time step size
    t = 0.0                                                     # initial time

    xc_vec = Array(xc)
    xv_vec = Array(xv)

    linplots = []

    # Initial plotting
    fig = Figure(size=(1000, 800), fontsize=20)
    #ax1 = Axis(fig[1,1], title="Pressure", ylabel="Pressure", xlabel="Domain",)
    #ax2 = Axis(fig[3,1], title="Velocity", ylabel="Velocity", xlabel="Domain")
    ax1 = Axis(fig[1,1], title="Density advection" , ylabel="Density", xticklabelsvisible=false, xticksvisible=false)
    ax2 = Axis(fig[1,2], title="Velocity", ylabel="Velocity", xticklabelsvisible=false, xticksvisible=false)
    ax3 = Axis(fig[2,1], title="Pressure", xlabel="Domain", ylabel="Pressure")
    ax4 = Axis(fig[2,2], title="Energy"  , xlabel="Domain", ylabel="Energy")
    # l0 = lines!(ax1, xc_vec, ρ)
    # push!(linplots, l0)
    # lines!(ax2, xv_vec, Vx)
    # lines!(ax3, xc_vec, P)
    # lines!(ax4, xc_vec, e)
    #save("../Plots/Navier-Stokes_acoustic_wave/with_advection_more_realistic_params/$(0).png", fig)
    #display(fig)

    for i = 1:nt
        ρ_old .= ρ

        # Conservation of mass
        Vxdρdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(ρ), dims=1) ./ dx
        ρdVxdt .= .-ρ[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        dρ = Vxdρdx .* dt .+ ρdVxdt .* dt
        ρ[2:end-1] .= ρ[2:end-1] .+ dρ

        # Conservation of momentum
        Vxdρdt .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* (av_x(ρ) .- av_x(ρ_old))
        ρdPdx .= .-(1.0 ./ av_x(ρ)) .* diff(P, dims=1) ./ dx
        ρVxdVxdx .= .-(1.0 ./ av_x(ρ)) .* (av_x(ρ) .* Vx[2:end-1]) .* diff(av_x(Vx), dims=1) ./ dx
        VxdρVxdx .= .-(1.0 ./ av_x(ρ)) .* Vx[2:end-1] .* diff(ρ .* av_x(Vx), dims=1) ./ dx
        Vx[2:end-1] .= Vx[2:end-1] .+ ρdPdx .* dt .+ ρVxdVxdx .* dt .+ VxdρVxdx .* dt .+ Vxdρdt .* dt

        # Conservation of energy
        EdρVxdx.= .-E[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        VxdPdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(P), dims=1) ./ dx
        PdVxdx .= .-P[2:end-1] .* diff(Vx[2:end-1], dims=1) ./ dx
        ρVxdEdx .= .-av_x(Vx[2:end-1]) .* diff(av_x(E), dims=1) ./ dx
        E[2:end-1] .= E[2:end-1] .+ EdρVxdx .* dt .+ VxdPdx .* dt .+ PdVxdx .* dt .+ ρVxdEdx .* dt

        # Equation of state for pressure difference
        dP = c.^2.0 .* dρ
        P[2:end-1] .= P[2:end-1] .+ dP

        # Internal energy calculation
        e = P ./ (γ - 1.0) ./ ρ
        
        t += dt
        if i % divisor == 0
            fig2 = Figure(size=(1000, 800))
            #ax1 = Axis(fig2[1,1], title="Pressure, time = $t")#, limits=(nothing, nothing, -0.25, 0.25))
            #ax2 = Axis(fig2[2,1], title="Velocity")#, limits=(nothing, nothing, -0.25, 0.25))
            #li = lines!(ax1, xc_vec, P, label="time = $t")
            #push!(linplots, li)
            ax1 = Axis(fig2[1,1], title="Density, time = $t")
            ax2 = Axis(fig2[1,2], title="Velocity")
            ax3 = Axis(fig2[2,1], title="Pressure")#, limits=(nothing, nothing, P0, P_max))
            ax4 = Axis(fig2[2,2], title="Energy")
            l0 = lines!(ax1, xc_vec, ρ)
            push!(linplots, l0)
            lines!(ax2, xv_vec, Vx)
            lines!(ax3, xc_vec, P)
            lines!(ax4, xc_vec, e)
            #lines!(ax2, xv_vec[2:end-1], Vx[2:end-1])
            #save("../Plots/Navier-Stokes_acoustic_wave/with_advection_more_realistic_params/$(i).png", fig2)
            # display(fig2)
        end
    end
    # Legend(fig[3,:], linplots, string.(0:dt*divisor:dt*nt), "Total time", nbanks=Int(floor((nt/divisor)+1)), tellhight = false, tellwidth = false)
    # rowsize!(fig.layout, 3, 40)
    # #save("/home/nils/Masterthesis_code/Plots/Nils_noncons_Euler-equations_acoustic_wave/with_realistic_parameters/4_in_one_acoustic.png", fig)
    # display(fig)
    return ρ, Vx, P, e, xc, xv
end

nt = [10000, 1000000]
ρ1, Vx1, P1, e1, xc1, xv1 = ac_wave1D(nt[1])
ρ_noadvect1, Vx_noadvect1, P_noadvect1, e_noadvect1 = ac_wave_1D_noadvect(nt[1])

error_ρ1  = abs.(ρ1 .- ρ_noadvect1)
error_Vx1 = abs.(Vx1 .- Vx_noadvect1)
error_P1  = abs.(P1 .- P_noadvect1)
error_e1  = abs.(e1 .- e_noadvect1)

ρ2, Vx2, P2, e2, xc2, xv2 = ac_wave1D(nt[2])
ρ_noadvect2, Vx_noadvect2, P_noadvect2, e_noadvect2 = ac_wave_1D_noadvect(nt[2])

error_ρ2  = abs.(ρ2 .- ρ_noadvect2)
error_Vx2 = abs.(Vx2 .- Vx_noadvect2)
error_P2  = abs.(P2 .- P_noadvect2)
error_e2  = abs.(e2 .- e_noadvect2)

###--------------------------- Normalized difference to both models / Absolute percentage error (APE) ---------------------------###

norm_error_ρ1  = error_ρ1  ./ abs.(ρ1) 
norm_error_Vx1 = error_Vx1 ./ abs.(Vx1)
norm_error_P1  = error_P1  ./ abs.(P1)
norm_error_e1  = error_e1  ./ abs.(e1)

norm_error_ρ2  = error_ρ2  ./ abs.(ρ2) 
norm_error_Vx2 = error_Vx2 ./ abs.(Vx2)
norm_error_P2  = error_P2  ./ abs.(P2)
norm_error_e2  = error_e2  ./ abs.(e2)

f = Figure(size=(1600, 800), fontsize=20)
Label(f[1,1:2], "APE of advection vs. no advection", tellwidth=false, font=:bold ,fontsize=30)
ax1 = Axis(f[2,1], title="Time = 0.0001", xlabel="Domain", ylabel="Absolute percentage error")
ax2 = Axis(f[2,2], title="Time = 0.01", xlabel="Domain", ylabel="Absolute percentage error")
scatterlines!(ax1, xc1, norm_error_ρ1 .* 100, label="Density")
scatterlines!(ax1, xc1, norm_error_Vx1[1:end-1] .* 100, label="Velocity")
scatterlines!(ax1, xc1, norm_error_P1 .* 100, label="Pressure")
scatterlines!(ax1, xc1, norm_error_e1 .* 100, label="Energy", linestyle=:dash, color=:purple)
scatterlines!(ax2, xc1, norm_error_ρ2 .* 100, label="Density")
scatterlines!(ax2, xc1, norm_error_Vx2[1:end-1] .* 100, label="Velocity")
scatterlines!(ax2, xc1, norm_error_P2 .* 100, label="Pressure")
scatterlines!(ax2, xc1, norm_error_e2 .* 100, label="Energy", linestyle=:dash, color=:purple)
axislegend(ax1 ;position=:rt)
axislegend(ax2 ;position=:rt)

save("APE_advection_vs_noadvection_2timesteps.png", f)
display(f)


