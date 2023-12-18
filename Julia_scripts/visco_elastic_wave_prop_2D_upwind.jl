using CairoMakie
using Infiltrator
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
    ρ_s = 2800.0                             # density
    #β = 1.0e-10                             # compressibility
    η = 1.0e21                            # viscosity
    μ = 1.0e10                             # shear modulus
    g_y  = 9.81                            # gravity
    P0 = 1.0e6                            # initial pressure at all points
    ν  = 0.4                            # Poisson's ratio
    λ  = (2.0 * μ * ν) / (1.0 - 2.0 * ν)     # Lamé parameter 2
    K_s = (2.0 * μ * (1.0 + ν)) / (3.0 * (1.0 - 2.0 * ν))                # model poisson ratio of the solid
    β  = 1.0 / K_s                         # compressibility
    E_s = 2.0 * μ * (1.0 + ν)                # Young's modulus
    V_p = sqrt(E_s / (ρ_s * ((1.0 + η) / E_s))) # P-wave velocity in Voigt model

    # Numerics
    nx = 101                            # number of nodes in x
    ny = 101                            # number of nodes in y
    dx = Lx / nx                        # step size in x
    dy = Ly / ny                        # step size in y
    nt = 5000                           # number of time steps
    t0 = 0.01                               # time of excitation      
    t = 0.0                             # initial time
    a = 40.0                              # duration of excitation

    # Grid definition
    xc = -(Lx - dx) / 2.0:dx:(Lx - dx) / 2.0        # grid nodes in x-direction
    yc = 0.0:dy:Ly-dy#-(Ly - dy) / 2.0:dy:(Ly - dy) / 2.0        # grid nodes in x-direction
    xv =   Lx       / 2.0:dx: Lx       / 2.0        # grid vertices in x-direction 
    yv =   Ly       / 2.0:dy: Ly       / 2.0        # grid vertices in x-direction

    # Nondimensionalization
    CharDim = GEO_units(length=1m, viscosity=1.0e19Pas, stress=100MPa)

    L       = 500.0                             # length of half the domain

    xc_geo = xc*m
    yc_geo = yc*m
    a_geo = a*m
    t_geo = t*s
    t0_geo = t0*s
    dx_geo = dx*m
    dy_geo = dy*m
    L_geo = L*m
    ρs_geo = ρ_s*kg/m^3
    βs_geo = β*1/Pa
    ηs_geo = η*Pa*s
    μs_geo = μ*Pa
    λs_geo = λ*Pa
    g_geo  = g_y*m/s^2
    P_geo  = P0*Pa
    E_geo  = E_s*Pa

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
    E_non  = nondimensionalize(E_geo, CharDim)

    Vp_non = sqrt(E_non / (ρs_non * ((1.0 + ηs_non) / E_non)))
    Vs_non = sqrt(μs_non / (ρs_non * ((1.0 + ηs_non) / μs_non)))


    # Allocations
    P = zeros(Float64, nx, ny)
    Pτxx = zeros(Float64, nx - 1, ny)
    Pτyy = zeros(Float64, nx, ny - 1)
    Frhox    = zeros(Float64, nx + 1, ny)
    Frhoy    = zeros(Float64, nx, ny + 1)
    drhodt   = zeros(Float64, nx, ny)
    rhoRes   = zeros(Float64, nx, ny)
    ρ_old    = ones(Float64, nx, ny) 
    ρ        = ones(Float64, nx, ny)         .* ρs_non
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
    c  = sqrt((1.0/βs_non)/ρs_non)                  # speed of sound / p-wave velocity
    dt = min(min(dx_non, dy_non) / c / 4.5, min(dx_non^2.0, dy_non^2.0) / ((4.0 / 3.0) * ηs_non / ρs_non) / 4.5)        # time step size                       

    P .= P_non .+ exp.((.-0.0000005 .* (xc_non/ 0.01).^2.0) .+ (.-0.0000005 .* ((yc_non .- 1000.0) / 0.01)'.^2.0))          # initial pressure distribution
    # initial equilibirum 
    for i in 1:ny
        τyy[:, (ny+1)-i] .= ρs_non .* g_non .* yc_non[i]
    end

    xc_vec = Vector(xc)
    yc_vec = Vector(yc)
    xticks = ["-500", "0", "500"]
    yticks = ["0", "500", "1000"]
    xtick_positions = range(start=-500, stop=500, step=500)
    ytick_positions = range(start=1000, stop=0, step=-500)

    # Initial plotting
    P_dim  = ustrip(dimensionalize(P, Pa, CharDim))
    τyy_dim = ustrip(dimensionalize(τyy, Pa, CharDim))
    xc_dim = ustrip(dimensionalize(xc_non, m, CharDim))
    yc_dim = ustrip(dimensionalize(yc_non, m, CharDim))

    fig = Figure(resolution=(2000,2000))
    ax = Axis(fig[1,1], title="t = $t")#, limits=(nothing, nothing, nothing, 1.1))
    #ax3 = Axis3(fig[1,1], title="time = $t")
    #limits!(ax3, -Lx / 2.0, Lx / 2.0, -Ly / 2.0, Ly / 2.0, -0.1, 1.0)
    #lines!(ax, xc_vec, P[:, Int(nx/2)])
    hm = heatmap!(ax, xc_dim, yc_dim, τyy_dim)
    #sur2 = surface!(ax3, ustrip(xc_dim), ustrip(yc_dim), ustrip(P_dim))
    Colorbar(fig[1,2], hm, label="Pressure")
    display(fig)

    for i = 1:nt
        β_vec .= 1.0 ./ P

        ρ_old .= ρ
        Frhox[2:end-1, :] .= (Vx[2:end-1, :] .>= 0.0).*Vx[2:end-1, :].*ρ[1:end-1, :] .+ (Vx[2:end-1, :] .< 0.0).*Vx[2:end-1, :].*ρ[2:end, :] # mass flux in x-direction (upwind scheme)
        Frhoy[:, 2:end-1] .= (Vy[:, 2:end-1] .>= 0.0).*Vy[:, 2:end-1].*ρ[:, 1:end-1] .+ (Vy[:, 2:end-1] .< 0.0).*Vy[:, 2:end-1].*ρ[:, 2:end] # mass flux in y-direction (upwind scheme)
        drhodt .= (ρ .- ρ_old)./dt                                                                           # time derivative of density
        rhoRes .= .-drhodt .- diff(Frhox, dims=1)./dx_non .- diff(Frhoy, dims=2)./dy_non        # updating residual of density
        ρ[2:end-1, 2:end-1] .= ρ[2:end-1, 2:end-1] .+ rhoRes[2:end-1, 2:end-1].* dt                          # updating density
        
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
        Pτxx .=.-diff(P, dims=1) ./ dx_non .+ diff(τxx, dims=1) ./ dx_non .+ diff(τxy[2:end-1, :], dims=2) ./ dy_non
        Pτyy .=.-diff(P, dims=2) ./ dy_non .+ diff(τyy, dims=2) ./ dy_non .+ diff(τxy[:, 2:end-1], dims=1) ./ dx_non .+ g_non .* (ρ[:, 1:end-1] .* (Vy[:, 2:end-1] .>= 0.0) .+ ρ[:, 2:end] .* (Vy[:, 2:end-1] .< 0.0))
        dVxdt[2:end-1, :] .= (1.0 ./ (ρ[1:end-1, :] .* (Vx[2:end-1, :] .>= 0.0) .+ ρ[2:end, :] .* (Vx[2:end-1, :] .< 0.0))) .* (Pτxx) .- Vx[2:end-1, :] .* diff(av_x(Vx), dims=1) ./ dx_non .+ av_x(av_y(Vy)) .* diff(av_x(Vx), dims=1) ./ dy_non
        dVydt[:, 2:end-1] .= (1.0 ./ (ρ[:, 1:end-1] .* (Vy[:, 2:end-1] .>= 0.0) .+ ρ[:, 2:end] .* (Vy[:, 2:end-1] .< 0.0))) .* (Pτyy) .- Vy[:, 2:end-1] .* diff(av_y(Vy), dims=2) ./ dy_non .+ av_y(av_x(Vx)) .* diff(av_y(Vy), dims=2) ./ dx_non
        Vx[2:end-1, :] .= Vx[2:end-1, :] .+ dVxdt[2:end-1, :] .* dt 
        Vy[:, 2:end-1] .= Vy[:, 2:end-1] .+ dVydt[: ,2:end-1] .* dt
        
        t += dt


        if i % 500 ==0
            fig = Figure(resolution=(2000,2000))
            Vy_dim = dimensionalize(Vy, m/s, CharDim)
            Vx_dim = dimensionalize(Vx, m/s, CharDim)
            Vx_av = av_x(ustrip(Vx_dim))
            Vy_av = av_y(ustrip(Vy_dim))
            data_plt = sqrt.(Vx_av.^2.0 .+ Vy_av.^2.0)
            P_plt = ustrip(dimensionalize(P, Pa, CharDim))
            t_dim = dimensionalize(t, s, CharDim)

            ax1 = Axis(fig[1,1], xticks=(xtick_positions, xticks), yticks=(ytick_positions, yticks),
                    yticklabelsize=25, xticklabelsize=25, xlabelsize=25, ylabelsize=25, title="time = $t_dim")
            ax2 = Axis(fig[2,1], xticks=(xtick_positions, xticks), yticks=(ytick_positions, yticks),
                    yticklabelsize=25, xticklabelsize=25, xlabelsize=25, ylabelsize=25, title="time = $t_dim")

            #lines!(ax2, xc_vec, P[:, Int(nx/2)])
            hm1 = heatmap!(ax1, xc_dim, yc_dim, P_plt, colormap=Reverse(:roma))#, colorrange=(P0, maximum(P)))
            hm2 = heatmap!(ax2, xc_dim, yc_dim, data_plt, colormap=Reverse(:roma))#, colorrange=(0, 8.0e-7))#, colorrange=(0.0, 1.0))
            #scatter!(ax2, x_circ, y_circ, color=:white, markersize=4.0)
            Colorbar(fig[1,2], hm1, label="Pressure [Pa]", labelsize=25, ticklabelsize=25)#, vertical=false)
            Colorbar(fig[2,2], hm2, label="Velocity [m/s]", labelsize=25, ticklabelsize=25)#, vertical=false)
            #ax3 = Axis3(fig2[1,2][1,1], title="time = $t")
            #limits!(ax3, -Lx / 2.0, Lx / 2.0, -Ly / 2.0, Ly / 2.0, -0.1, 1.0)
            #sur2 = surface!(ax3, xc_vec, yc_vec, P)
            display(fig)
            #save(".\\Plots\\Earthquake_2\\$(i).png", fig2)
            @infiltrate
        end
    end
end