using CairoMakie
using GeoParams

@views function av_x(B)
    B = 0.5 .* (B[2:end,:] .+ B[1:end-1,:])
end

@views function av_y(B)
    B = 0.5 .* (B[:,2:end] .+ B[:,1:end-1])
end

@views function av_xy(B)
    B = 0.25 .* (B[2:end,2:end] .+ B[1:end-1,2:end] .+ B[2:end,1:end-1] .+ B[1:end-1,1:end-1])
end

@views function visco_elastic_wave_2D()
    # Physics
    Lx = 1000.0                          # domain in x
    Ly = Lx                             # domain in y
    ρ = 2800.0                             # density
    #β = 1.0e-10                             # compressibility
    η = 1.0e21                            # viscosity
    μ = 1.0e10                             # shear modulus
    g_y  = 9.81                            # gravity
    P0 = 1.0e6                            # initial pressure at all points
    ν  = 0.4                            # Poisson's ratio
    λ  = (2.0 * μ * ν) / (1.0 - 2.0 * ν)     # Lamé parameter 2
    K = (2.0 * μ * (1.0 + ν)) / (3.0 * (1.0 - 2.0 * ν))                # model poisson ratio of the solid
    β  = 1.0 / K                         # compressibility

    # Numerics
    nx = 255                            # number of nodes in x
    ny = 255                            # number of nodes in y
    dx = Lx / nx                        # step size in x
    dy = Ly / ny                        # step size in y
    nt = 5000                           # number of time steps
    t0 = 0.01                               # time of excitation      
    t = 0.0                             # initial time
    a = 40.0                              # duration of excitation

    # Grid definition
    xc = -(Lx - dx) / 2.0:dx:(Lx - dx) / 2.0        # grid nodes in x-direction
    yc = 0.0+dy:dy:Ly#-(Ly - dy) / 2.0:dy:(Ly - dy) / 2.0        # grid nodes in x-direction
    xv =   Lx       / 2.0:dx: Lx       / 2.0        # grid vertices in x-direction 
    yv =   Ly       / 2.0:dy: Ly       / 2.0        # grid vertices in x-direction

    # Nondimensionalization
    CharDim = GEO_units(length=1km, viscosity=1.0e19Pas, stress=100MPa)

    L       = 500.0                             # length of half the domain

    xc_geo = xc*m
    yc_geo = yc*m
    a_geo = a*m
    t_geo = t*s
    t0_geo = t0*s
    dx_geo = dx*m
    dy_geo = dy*m
    L_geo = L*m
    ρs_geo = ρ*kg/m^3
    βs_geo = β*1/Pa
    ηs_geo = η*Pa*s
    μs_geo = μ*Pa
    λs_geo = λ*Pa
    g_geo  = g_y*m/s^2
    P_geo  = P0*Pa

    xc_non = nondimensionalize(xc_geo, CharDim)
    yc_non = nondimensionalize(yc_geo, CharDim)
    a_non = nondimensionalize(a_geo, CharDim)
    t_non = nondimensionalize(t_geo, CharDim)
    t0_non = nondimensionalize(t0_geo, CharDim)
    dx_non = nondimensionalize(dx_geo, CharDim)
    dy_non = nondimensionalize(dy_geo, CharDim)
    L_non = nondimensionalize(L_geo, CharDim)
    ρs_non = nondimensionalize(ρs_geo, CharDim)
    βs_non = nondimensionalize(βs_geo, CharDim)
    ηs_non = nondimensionalize(ηs_geo, CharDim)
    μs_non = nondimensionalize(μs_geo, CharDim)
    λs_non = nondimensionalize(λs_geo, CharDim)
    g_non  = nondimensionalize(g_geo, CharDim)
    P_non  = nondimensionalize(P_geo, CharDim)

    # Allocations
    P = zeros(Float64, nx, ny)
    ρ_vec    = ones(Float64, nx, ny)         .* ρs_non
    ρ_vx_vec = ones(Float64, nx + 1, ny)     .* ρs_non
    ρ_vy_vec = ones(Float64, nx, ny + 1)     .* ρs_non
    β_vec    = ones(Float64, nx, ny)         .* βs_non
    η_vec    = ones(Float64, nx, ny)         .* ηs_non
    η_c_vec  = ones(Float64, nx + 1, ny + 1) .* ηs_non
    μ_vec    = ones(Float64, nx, ny)         .* μs_non
    μ_c_vec  = ones(Float64, nx + 1, ny + 1) .* μs_non
    divV = ones(Float64, nx, ny)
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
    τII_vec = []
    εII_vec = []

    # Initial conditions
    c  = sqrt(1.0/βs_non/ρs_non)                  # speed of sound / p-wave velocity
    dt = 0.1e-15 #min(min(dx_non, dy_non) / c / 4.5, min(dx_non^2.0, dy_non^2.0) / ((4.0 / 3.0) * ηs_non / ρs_non) / 4.5)        # time step size                       

    P .= P_non .+ exp.((.-5.0 .* (xc_non/ 0.01).^2.0) .+ (.-5.0 .* ((yc_non .- 1.0) / 0.01)'.^2.0))          # initial pressure distribution
    # initial equilibirum 
    for i in 1:ny
        τyy[:, (ny+1)-i] .= ρs_non .* g_non .* yc_non[i]
    end

    xc_vec = Vector(xc)
    yc_vec = Vector(yc)

    # Initial plotting
    P_dim  = dimensionalize(P, Pa, CharDim)
    τyy_dim = dimensionalize(τyy, Pa, CharDim)
    xc_dim = dimensionalize(xc_non, m, CharDim)
    yc_dim = dimensionalize(yc_non, m, CharDim)

    fig = Figure()
    ax = Axis(fig[1,1], title="t = $t")#, limits=(nothing, nothing, nothing, 1.1))
    #ax3 = Axis3(fig[1,1], title="time = $t")
    #limits!(ax3, -Lx / 2.0, Lx / 2.0, -Ly / 2.0, Ly / 2.0, -0.1, 1.0)
    #lines!(ax, xc_vec, P[:, Int(nx/2)])
    hm = heatmap!(ax, ustrip(xc_dim), ustrip(yc_dim), ustrip(τyy_dim))
    #sur2 = surface!(ax3, ustrip(xc_dim), ustrip(yc_dim), ustrip(P_dim))
    Colorbar(fig[1,2], hm, label="Pressure")
    fig

    for i = 1:nt
        divV .= diff(Vx, dims=1) ./ dx_non .+ diff(Vy, dims=2) ./ dy_non
        dPdt .= .-(1.0 ./ β_vec) .* divV
        P .= P .+ dPdt .* dt
        εxx .= diff(Vx, dims=1) ./ dx_non .- (1.0 ./ 3.0) .* divV
        εyy .= diff(Vy, dims=2) ./ dy_non .- (1.0 ./ 3.0) .* divV
        εxy[2:end-1, 2:end-1] .= 0.5 .* (diff(Vx[2:end-1, :], dims=2) ./ dy_non .+ diff(Vy[:, 2:end-1], dims=1) ./ dx_non)
        dτxxdt .= 2.0 .* μ_vec .* εxx .- (μ_vec ./ η_vec) .* τxx
        dτyydt .= 2.0 .* μ_vec .* εyy .- (μ_vec ./ η_vec) .* τyy         
        dτxydt[2:end-1, 2:end-1] .= 2.0 .* μ_c_vec[2:end-1, 2:end-1] .* εxy[2:end-1, 2:end-1] .- (μ_c_vec[2:end-1, 2:end-1] ./ η_c_vec[2:end-1, 2:end-1]) .* τxy[2:end-1, 2:end-1]             
        τxx .= τxx .+ dτxxdt .* dt
        τyy .= τyy .+ dτyydt .* dt
        τxy .= τxy .+ dτxydt .* dt
        #τxx[127, 255] = (t_non - t0_non >= 0.0) * (-2.0 * a_non * (t_non - t0_non) * exp(-a_non * (t_non - t0_non)^2.0))
        #τyy[127, 255] = (t_non - t0_non >= 0.0) * (-2.0 * a_non * (t_non - t0_non) * exp(-a_non * (t_non - t0_non)^2.0))
        #τxy[127, 255] = (t_non - t0_non >= 0.0) * (-2.0 * a_non * (t_non - t0_non) * exp(-a_non * (t_non - t0_non)^2.0))
        dVxdt[2:end-1, :] .= (1.0 ./ ρ_vx_vec[2:end-1, :]) .* (.-diff(P, dims=1) ./ dx_non .+ diff(τxx, dims=1) ./ dx_non .+ diff(τxy[2:end-1, :], dims=2) ./ dy_non)
        dVydt[:, 2:end-1] .= (1.0 ./ ρ_vy_vec[:, 2:end-1]) .* (.-diff(P, dims=2) ./ dy_non .+ diff(τyy, dims=2) ./ dy_non .+ diff(τxy[:, 2:end-1], dims=1) ./ dx_non .+ g_non .* ρ_vy_vec[:, 2:end-1])
        Vx[2:end-1, :] .= Vx[2:end-1, :] .+ dVxdt[2:end-1, :] .* dt 
        Vy[:, 2:end-1] .= Vy[:, 2:end-1] .+ dVydt[: ,2:end-1] .* dt
        
        t += dt

        if i % 50 ==0
            Vy_dim = dimensionalize(Vy, m/s, CharDim)
            Vx_dim = dimensionalize(Vx, m/s, CharDim)
            Vx_av = av_x(ustrip(Vx_dim))
            Vy_av = av_y(ustrip(Vy_dim))
            data_plt = sqrt.(Vx_av.^2.0 .+ Vy_av.^2.0)
            P_plt = dimensionalize(P, Pa, CharDim)
            #τxx_plt = ustrip(dimensionalize(τxx[127, 230], Pa, CharDim))
            #τyy_plt = ustrip(dimensionalize(τyy[127, 230], Pa, CharDim))
            #τxy_plt = ustrip(dimensionalize(τxy[127, 230], Pa, CharDim))
            #εxx_plt = ustrip(dimensionalize(εxx[127, 230], s^-1.0, CharDim))
            #εyy_plt = ustrip(dimensionalize(εyy[127, 230], s^-1.0, CharDim))
            #εxy_plt = ustrip(dimensionalize(εxy[127, 230], s^-1.0, CharDim))
            t_dim = dimensionalize(t, s, CharDim)
            #τII = 0.5 .* τxx_plt.*τyy_plt .- τxy_plt.^2.0
            #εII = 0.5 .* εxx_plt.*εyy_plt .- εxy_plt.^2.0
            #push!(τII_vec, τII)
            #push!(εII_vec, εII)

            fig2 = Figure()
            ax2 = Axis(fig2[1,1], title="time = $t_dim")#, limits=(nothing, nothing, nothing, 1.1))
            #ax3 = Axis3(fig2[1,1], title="time = $t")
            #limits!(ax3, -Lx / 2.0, Lx / 2.0, -Ly / 2.0, Ly / 2.0, -0.1, 1.0)
            #ylims!(ax, -1.2, 1.2)
            #lines!(ax2, εII, τII)
            #sur2 = surface!(ax3, xc_vec, yc_vec, P)
            hm2 = heatmap!(ax2, ustrip(xc_dim), ustrip(yc_dim), data_plt, colormap=Reverse(:roma))#, colorrange=(0.0, 1.0e-14))
            #hm2 = heatmap!(ax2,  ustrip(xc_dim), ustrip(yc_dim), ustrip(P_plt), colormap=Reverse(:roma))#, colorrange=(0.0, 1.0))
            #Colorbar(fig2[1,2], hm2, label="Pressure")
            Colorbar(fig2[1,2], hm2, label="Velocity [m/s]")
            display(fig2)
        end
    end
end