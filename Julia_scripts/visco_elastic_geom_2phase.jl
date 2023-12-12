using CairoMakie
using GeoParams
using Infiltrator

function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

function inpolygon(X::Array{T,2}, Y::Array{T,2}, XV::Array{T,1}, YV::Array{T,1}) where T <: Real
    n = length(XV)
    inpts = zeros(Bool, size(X))
    onpts = zeros(Bool, size(X))
    
    for i in 1:size(X, 1), j in 1:size(X, 2)
        xi = X[i,j]
        yi = Y[i,j]
        wn = 0 
        
        for k in 1:n
            if k < n
                x1 = XV[k]
                y1 = YV[k]
                x2 = XV[k+1]
                y2 = YV[k+1]
            else
                x1 = XV[n]
                y1 = YV[n]
                x2 = XV[1]
                y2 = YV[1]
            end
            
            if y1 <= yi
                if y2 > yi && cross2d(x2-x1, y2-y1, xi-x1, yi-y1) >= 0
                    wn += 1
                end
            else
                if y2 <= yi && cross2d(x2-x1, y2-y1, xi-x1, yi-y1) <= 0
                    wn -= 1
                end
            end
            
            if y1 == yi && x1 == xi
                onpts[i,j] = true
            elseif (y2 == yi && x2 == xi) || (y1 < yi && y2 > yi) || (y2 < yi && y1 > yi)
                onpts[i,j] = true
            end
        end
        
        if wn != 0
            inpts[i,j] = true
        end
    end
    
    return inpts
end

@views function cross2d(a, b, c, d)
    return a*d - b*c
end

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
    Lx = 1.0e3                          # domain in x
    Ly = Lx                             # domain in y
    ρ_a = 1.225                           # density
    ρ_s = 2.8e3                           # density
    η_a = 1.0e-5                            # viscosity
    η_s = 1.0e21                             # viscosity
    μ_a = 1.0e-2                             # shear modulus
    μ_s = 1.0e10                             # shear modulus
    P0 = 1.0e5                            # initial pressure at all points
    g_y = 9.81                             # gravity
    ν_s  = 0.4                            # Poisson's ratio
    λ_s  = (2.0 * μ_s * ν_s) / (1.0 - 2.0 * ν_s)     # Lamé parameter 2
    K_s = (2.0 * μ_s * (1.0 + ν_s)) / (3.0 * (1.0 - 2.0 * ν_s))                # model poisson ratio of the solid
    β_a  = 1.0 / 141.0e3                         # compressibility
    β_s  = 1.0 / K_s                        # compressibility

    α = sqrt(μ_s / ρ_s)                  # shear wave velocity
    ξ = sqrt((λ_s + 2.0 * μ_s) / ρ_s)      # p-wave velocity

    # Numerics
    nx = 301                            # number of nodes in x
    ny = 301                            # number of nodes in y
    dx = Lx / nx                        # step size in x
    dy = Ly / ny                        # step size in y
    nt = 10000                           # number of time steps
    t = 0.0                                                      # initial time

    # Grid definition
    xc = -(Lx - dx) / 2.0:dx:(Lx - dx) / 2.0        # grid nodes in x-direction
    yc = -(Ly - dy) / 2.0:dy:(Ly - dy) / 2.0        # grid nodes in x-direction
    xv =  -Lx       / 2.0:dx: Lx       / 2.0        # grid vertices in x-direction 
    yv =  -Ly       / 2.0:dy: Ly       / 2.0        # grid vertices in x-direction

    # Nondimensionalization
    CharDim = GEO_units(length=1m, viscosity=1.0e19Pas, stress=100MPa)

    L       = 500.0                             # length of half the domain
    R       = 8.314                             # gas constant
    T       = 288.0                             # temperature
    M       = 0.028                             # molar mass of air

    xc_geo = xc*m
    yc_geo = yc*m
    t_geo = t*s
    dx_geo = dx*m
    dy_geo = dy*m
    L_geo = L*m
    R_geo = R*J/mol/K
    T_geo = T*K
    M_geo = M*kg
    ρa_geo = ρ_a*kg/m^3
    ρs_geo = ρ_s*kg/m^3
    βa_geo = β_a*1/Pa
    βs_geo = β_s*1/Pa
    ηa_geo = η_a*Pa*s
    ηs_geo = η_s*Pa*s
    μa_geo = μ_a*Pa
    μs_geo = μ_s*Pa
    λs_geo = λ_s*Pa
    g_geo  = g_y*m/s^2
    P_geo  = P0*Pa

    xc_non = nondimensionalize(xc_geo, CharDim)
    yc_non = nondimensionalize(yc_geo, CharDim)
    t_non  = nondimensionalize(t_geo, CharDim)
    dx_non = nondimensionalize(dx_geo, CharDim)
    dy_non = nondimensionalize(dy_geo, CharDim)
    L_non  = nondimensionalize(L_geo, CharDim)
    R_non  = nondimensionalize(R_geo, CharDim)
    T_non  = nondimensionalize(T_geo, CharDim)
    M_non  = nondimensionalize(M_geo, CharDim)
    ρa_non = nondimensionalize(ρa_geo, CharDim)
    ρs_non = nondimensionalize(ρs_geo, CharDim)
    βa_non = nondimensionalize(βa_geo, CharDim)
    βs_non = nondimensionalize(βs_geo, CharDim)
    ηa_non = nondimensionalize(ηa_geo, CharDim)
    ηs_non = nondimensionalize(ηs_geo, CharDim)
    μa_non = nondimensionalize(μa_geo, CharDim)
    μs_non = nondimensionalize(μs_geo, CharDim)
    λs_non = nondimensionalize(λs_geo, CharDim)
    g_non  = nondimensionalize(g_geo, CharDim)
    P_non  = nondimensionalize(P_geo, CharDim)

    # Allocations
    P = zeros(Float64, nx, ny)
    ρ    = ones(Float64, nx, ny)        .* ρa_non
    ρ_vx = ones(Float64, nx + 1, ny)    .* ρa_non
    ρ_vy = ones(Float64, nx, ny + 1)    .* ρa_non
    β    = ones(Float64, nx, ny)        .* βa_non
    η    = ones(Float64, nx, ny)        .* ηa_non
    η_c  = ones(Float64, nx + 1, ny + 1).* ηa_non
    μ    = ones(Float64, nx, ny)        .* μa_non
    μ_c  = ones(Float64, nx + 1, ny + 1).* μa_non
    divV = zeros(Float64, nx, ny)
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
    dτxxdt2 = zeros(Float64, nx, ny)
    dτyydt = zeros(Float64, nx, ny)
    dτyydt2 = zeros(Float64, nx, ny)
    dτxydt = zeros(Float64, nx + 1, ny + 1)
    dτxydt2 = zeros(Float64, nx + 1, ny + 1)
    τII_vec = []
    εII_vec = []

    maskrho_air = zeros(Float64, nx, ny)
    maskrho_solid = zeros(Float64, nx, ny)
    maskVx_air = zeros(Float64, nx + 1, ny)
    maskVx_solid = zeros(Float64, nx + 1, ny)
    maskVy_air = zeros(Float64, nx, ny + 1)
    maskVy_solid = zeros(Float64, nx, ny + 1)
    maskμc_air   = zeros(Float64, nx + 1, ny + 1)
    maskμc_solid = zeros(Float64, nx + 1, ny + 1)
    radrho = zeros(Float64, nx, ny)
    radVx = zeros(Float64, nx + 1, ny)
    radVy = zeros(Float64, nx, ny + 1)
    radμc = zeros(Float64, nx + 1, ny + 1)

    ### Geometrie part ---------------------------
    x2dc, y2dc = meshgrid(xc,yc)                                        # 2d mesh of x- and y-coordinates of density nodes
    x2dVx, y2dVx = meshgrid(xv,yc)                                      # 2d mesh of x- and y-coordinates of velocity nodes in x-direction
    x2dVy, y2dVy = meshgrid(xc,yv)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction
    x2dμc, y2dμc = meshgrid(xv,yv)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction

    x2dc_non, y2dc_non   = meshgrid(xc_non,yc_non)

    locX      = 0.0                                                     # x-coordinate of circle
    locY      = -0.35*Ly                                                # y-coordinate of circle
    diam      = 0.2*Lx                                                  # diameter of circle
    radrho    .= sqrt.((x2dc .-locX).^2.0 .+ (y2dc .-locY).^2.0)        # radius of circle
    radVx     .= sqrt.((x2dVx.-locX).^2.0 .+ (y2dVx.-locY).^2.0)        # radius of circle
    radVy     .= sqrt.((x2dVy.-locX).^2.0 .+ (y2dVy.-locY).^2.0)        # radius of circle
    radμc     .= sqrt.((x2dμc.-locX).^2.0 .+ (y2dμc.-locY).^2.0)        # radius of circle
    Xp        = [-1.0/2.0, -1.0/2.0, -0.3, -1.0/8.0, -0.01, -0.01, 0.01, 0.01, 1.0/8.0, 0.3, 1.0/2.0, 1.0/2.0].*Lx          # x-coordinates of polygon
    Yp        = [-1.0/2.0, -1.0/5.0, -1.0/5.0, -0.1, -0.15,  -0.28, -0.28, -0.15,  -0.1, -1.0/5.0, -1.0/5.0, -1.0/2.0].*Ly  # y-coordinates of polygon
    Xp2       = [-0.02, -0.02,  0.02,  0.02].*Lx                                                                            # x-coordinates of polygon
    Yp2       = [-0.15, -0.35, -0.35, -0.15].*Ly                                                                             # y-coordinates of polygon

    inpolyrho   = inpolygon(x2dc,y2dc,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon
    inpolyVx    = inpolygon(x2dVx,y2dVx,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon
    inpolyVy    = inpolygon(x2dVy,y2dVy,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon
    inpolyμc    = inpolygon(x2dμc,y2dμc,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon

    maskrho_air = 1.0 .- inpolyrho                                             # mask for density
    maskrho_air[findall(<(diam ./ 2.0), radrho)] .= 1.0                               # mask for density in circle
    maskVx_air  = 1.0 .- inpolyVx                               # mask for velocity in x-direction
    maskVx_air[findall(<(diam ./ 2.0), radVx)] .= 1.0                                 # mask for velocity in x-direction in circle
    maskVy_air  = 1.0 .- inpolyVy                               # mask for velocity in y-direction
    maskVy_air[findall(<(diam ./ 2.0), radVy)] .= 1.0                                 # mask for velocity in y-direction in circle
    maskμc_air  = 1.0 .- inpolyμc                               # mask for velocity in y-direction
    maskμc_air[findall(<(diam ./ 2.0), radμc)] .= 1.0                                 # mask for velocity in y-direction in circle

    maskrho_solid = 1.0 .- maskrho_air
    maskVx_solid  = 1.0 .- maskVx_air
    maskVy_solid  = 1.0 .- maskVy_air
    maskμc_solid  = 1.0 .- maskμc_air

    x_circ_ind = []
    y_circ_ind = []

    for y= 2:size(maskrho_air)[1] - 1
        for x = 2:size(maskrho_air)[2] - 1
            if (maskrho_air[x, y] == 0.0 && abs(maskrho_air[x, y] - maskrho_air[x - 1, y]) == 1.0)
                push!(x_circ_ind, x)
                push!(y_circ_ind, y)
            elseif (maskrho_air[x, y] == 0.0 && abs(maskrho_air[x, y] - maskrho_air[x + 1, y]) == 1.0 )
                push!(x_circ_ind, x)
                push!(y_circ_ind, y)
            elseif (maskrho_air[x, y] == 0.0 && abs(maskrho_air[x, y] - maskrho_air[x, y - 1]) == 1.0)
                push!(x_circ_ind, x)
                push!(y_circ_ind, y)
            elseif (maskrho_air[x, y] == 0.0 && abs(maskrho_air[x, y] - maskrho_air[x, y + 1]) == 1.0)
                push!(x_circ_ind, x)
                push!(y_circ_ind, y)
            end
        end
    end

    # Initial conditions
    c  = sqrt(1.0/βs_non/ρs_non)                  # speed of sound / p-wave velocity
    dt = 5.0e-16 #min(min(dx_non, dy_non) / c / 4.5, min(dx_non^2.0, dy_non^2.0) / ((4.0 / 3.0) * ηs_non / ρs_non) / 4.5)        # time step size                       
    # 1.0e-15 for seimic wave and 1.0e-23 for atmospheric waves


    #P .= P_non .* exp.(-g_non .* (y2dc_non .+ 1.0 ./ 5.0 .* L_non) .* M_non ./ T_non ./ R_non)
    P .+= P_non .+ exp.((.-0.5 .* ((xc_non)/ 10.0).^2.0) .+ (.-0.5 .* ((yc_non .+ 500.0) / 10.0)'.^2.0))          # initial pressure distribution 
    
    τyy .= reverse(cumsum(ρ .* g_non #=.* reverse(maskrho_solid, dims=2)=# .* dy_non, dims=2), dims=2)


    #@. ρ    += ρa_non * (maskrho_air == 1.0) + (maskrho_solid  == 1.0) * ρs_non                 # density
    #@. ρ_vx += ρa_non * (maskVx_air  == 1.0) + (maskVx_solid   == 1.0) * ρs_non                 # density
    #@. ρ_vy += ρa_non * (maskVy_air  == 1.0) + (maskVy_solid   == 1.0) * ρs_non                 # density
    #@. β    += βa_non * (maskrho_air == 1.0) + (maskrho_solid  == 1.0) * βs_non               # compressibility
    #@. η    += ηa_non * (maskrho_air == 1.0) + (maskrho_solid  == 1.0) * ηs_non               # viscosity
    #@. η_c  += ηa_non * (maskμc_air  == 1.0) + (maskμc_solid   == 1.0) * ηs_non               # viscosity
    #@. μ    += μa_non * (maskrho_air == 1.0) + (maskrho_solid  == 1.0) * μs_non               # shear modulus
    #@. μ_c  += μa_non * (maskμc_air  == 1.0) + (maskμc_solid   == 1.0) * μs_non               # shear modulus

    xc_vec = Vector(xc)
    yc_vec = Vector(yc)

    # Initial plotting
    x_circ = zeros(length(x_circ_ind))
    y_circ = zeros(length(y_circ_ind))
    x_circ = xc[x_circ_ind]
    y_circ = yc[y_circ_ind]

    P_dim  = ustrip(dimensionalize(P, Pa, CharDim))
    τyy_dim= ustrip(dimensionalize(τyy, Pa, CharDim))
    xc_dim = ustrip(dimensionalize(xc_non, m, CharDim))
    yc_dim = ustrip(dimensionalize(yc_non, m, CharDim))

    fig = Figure()
    ax = Axis(fig[1,1], title="t = $t")#, limits=(nothing, nothing, nothing, 1.1))
    #ax3 = Axis3(fig[1,1], title="time = $t")
    #limits!(ax3, -Lx / 2.0, Lx / 2.0, -Ly / 2.0, Ly / 2.0, -0.1, 1.0)
    #lines!(ax, xc_vec, P[:, Int(nx/2)])
    scatter!(ax, x_circ, y_circ, color=:white, markersize=4.0)
    hm = heatmap!(ax, xc_dim, yc_dim, P_dim, colormap=Reverse(:roma))
    #sur2 = surface!(ax3, ustrip(xc_dim), ustrip(yc_dim), ustrip(P_dim))
    Colorbar(fig[1,2], hm, label="Pressure")
    display(fig)

    for i = 1:nt
        divV .= diff(Vx, dims=1) ./ dx_non .+ diff(Vy, dims=2) ./ dy_non
        dPdt .= .-(1.0 ./ β) .* divV
        P .= P .+ dPdt .* dt
        εxx .= diff(Vx, dims=1) ./ dx_non .- (1.0 ./ 3.0) .* divV
        εyy .= diff(Vy, dims=2) ./ dy_non .- (1.0 ./ 3.0) .* divV
        εxy[2:end-1, 2:end-1] .= 0.5 .* (diff(Vx[2:end-1, :], dims=2) ./ dy_non .+ diff(Vy[:, 2:end-1], dims=1) ./ dx_non)
        #dτxxdt .=  2.0 .* (1.0 ./ ((1.0 ./ η) .+ (1.0 ./ (μ .* dt)))) .* (εxx .+ (τxx ./ (2.0 .* μ .* dt)))
        dτxxdt .= 2.0 .* μ .* εxx .- (μ ./ η) .* τxx
        #dτyydt .=  2.0 .* (1.0 ./ ((1.0 ./ η) .+ (1.0 ./ (μ .* dt)))) .* (εyy .+ (τyy ./ (2.0 .* μ .* dt)))         
        dτyydt .= 2.0 .* μ .* εyy .- (μ ./ η) .* τyy         
        #dτxydt[2:end-1, 2:end-1] .= 2.0 .* (1.0 ./ ((1.0 ./ η_c[2:end-1, 2:end-1]) .+ (1.0 ./ (μ_c[2:end-1, 2:end-1] .* dt)))) .* (εxy[2:end-1, 2:end-1] .+ (τxy[2:end-1, 2:end-1] ./ (2.0 .* μ_c[2:end-1, 2:end-1] .* dt)))             
        dτxydt[2:end-1, 2:end-1]  .= 2.0 .* μ_c[2:end-1, 2:end-1] .* εxy[2:end-1, 2:end-1] .- (μ_c[2:end-1, 2:end-1] ./ η_c[2:end-1, 2:end-1]) .* τxy[2:end-1, 2:end-1]             
        #@infiltrate
        τxx .= τxx .+ dτxxdt .* dt
        τyy .= τyy .+ dτyydt .* dt
        τxy .= τxy .+ dτxydt .* dt
        #@infiltrate
        #τxx[127, 255] = (t_non - t0_non >= 0.0) * (-2.0 * a_non * (t_non - t0_non) * exp(-a_non * (t_non - t0_non)^2.0))
        #τyy[127, 255] = (t_non - t0_non >= 0.0) * (-2.0 * a_non * (t_non - t0_non) * exp(-a_non * (t_non - t0_non)^2.0))
        #τxy[127, 255] = (t_non - t0_non >= 0.0) * (-2.0 * a_non * (t_non - t0_non) * exp(-a_non * (t_non - t0_non)^2.0))
        dVxdt[2:end-1, :] .= (1.0 ./ ρ_vx[2:end-1, :]) .* (.-diff(P, dims=1) ./ dx_non .+ diff(τxx, dims=1) ./ dx_non .+ diff(τxy[2:end-1, :], dims=2) ./ dy_non)                              .- Vx[2:end-1, :] .* diff(av_x(Vx), dims=1) ./ dx_non .+ av_x(av_y(Vy)) .* diff(av_x(Vx), dims=1) ./ dy_non
        dVydt[:, 2:end-1] .= (1.0 ./ ρ_vy[:, 2:end-1]) .* (.-diff(P, dims=2) ./ dy_non .+ diff(τyy, dims=2) ./ dy_non .+ diff(τxy[:, 2:end-1], dims=1) ./ dx_non .+ g_non .* ρ_vy[:, 2:end-1]) .- Vy[:, 2:end-1] .* diff(av_y(Vy), dims=2) ./ dy_non .+ av_y(av_x(Vx)) .* diff(av_y(Vy), dims=2) ./ dx_non
        Vx[2:end-1, :] .= Vx[2:end-1, :] .+ dVxdt[2:end-1, :] .* dt
        Vy[:, 2:end-1] .= Vy[:, 2:end-1] .+ dVydt[: ,2:end-1] .* dt
        
        t += dt

        if i % 10 ==0
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
            hm2 = heatmap!(ax2, xc_dim, yc_dim, data_plt, colormap=Reverse(:roma))#, colorrange=(0.0, 0.4))
            #hm2 = heatmap!(ax2,  ustrip(xc_dim), ustrip(yc_dim), ustrip(P_plt), colormap=Reverse(:roma))#, colorrange=(0.0, 1.0))
            #Colorbar(fig2[1,2], hm2, label="Pressure")
            Colorbar(fig2[1,2], hm2, label="Velocity [m/s]")
            display(fig2)
        end
    end
end