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

@views function visco_elastic_atmowave_2D_test()
    # Physics
    Lx = 1.0e3                          # domain in x
    Ly = Lx                             # domain in y
    ρ_a = 1.225                           # density
    ρ_c = 2.7e3                           # density chamber
    ρ_s = 2.8e3                           # density
    η_a = 1.0e-5                            # viscosity
    η_s = 1.0e19                             # viscosity
    μ_a = 1.0e-5                            # shear modulus
    μ_s = 1.0e11                             # shear modulus
    P0 = 1.0e5                            # initial pressure at all points
    g_y = 9.81                             # gravity
    ν  = 0.4                            # Poisson's ratio
    λ  = (2.0 * μ_s * ν) / (1.0 - 2.0 * ν)     # Lamé parameter 2
    K_s = (2.0 * μ_s * (1.0 + ν)) / (3.0 * (1.0 - 2.0 * ν))                # model poisson ratio of the solid
    β_a  = 1.0/141.0e3                         # compressibility
    β_s  = 1.0 / K_s                        # compressibility

    α = sqrt(μ_s / ρ_s)                  # shear wave velocity
    ξ = sqrt((λ + 2.0 * μ_s) / ρ_s)      # p-wave velocity

    # Numerics
    nx = 301                            # number of nodes in x
    ny = 301                            # number of nodes in y
    dx = Lx / nx                        # step size in x
    dy = Ly / ny                        # step size in y
    nt = 10000                          # number of time steps
    t = 0.0                             # initial time

    # Grid definition
    xc = -(Lx - dx) / 2.0:dx:(Lx - dx) / 2.0        # grid nodes in x-direction
    yc = -(Ly - dy) / 2.0:dy:(Ly - dy) / 2.0        # grid nodes in x-direction
    xv =  -Lx       / 2.0:dx: Lx       / 2.0        # grid vertices in x-direction 
    yv =  -Ly       / 2.0:dy: Ly       / 2.0        # grid vertices in x-direction

    # Nondimensionalization
    CharDim = GEO_units(length=1km, viscosity=1.0e19Pas, stress=100MPa)
    L       = 500.0                             # length of half the domain
    R       = 8.314                             # gas constant
    T       = 288.0                             # temperature
    M       = 0.028                             # molar mass of air

    # Making values GeoUnits for nondimensionalization
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
    λs_geo = λ*Pa
    g_geo  = g_y*m/s^2
    P_geo  = P0*Pa

    # Nondimensionalization
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
    P    = zeros(Float64, nx, ny)
    ρ    = zeros(Float64, nx, ny)         #.* ρa_non
    ρ_vx = zeros(Float64, nx + 1, ny)     #.* ρa_non
    ρ_vy = zeros(Float64, nx, ny + 1)     #.* ρa_non
    β    = zeros(Float64, nx, ny)         #.* βa_non
    η    = zeros(Float64, nx, ny)         #.* ηa_non
    η_c  = zeros(Float64, nx + 1, ny + 1) #.* ηa_non
    μ    = zeros(Float64, nx, ny)         #.* μa_non
    μ_c  = zeros(Float64, nx + 1, ny + 1) #.* μa_non
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
    dτyydt = zeros(Float64, nx, ny)
    dτxydt = zeros(Float64, nx + 1, ny + 1)
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
    Yp2       = [-0.15, -0.35, -0.35, -0.15].*Ly 

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
    
    ### End of geometrie part --------------------

    # Initial conditions
    cs_non  = sqrt(1.0/βs_non/ρs_non)                                                              # speed of sound / p-wave velocity in solid
    ca_non  = sqrt(1.0/βa_non/ρa_non)                                                              # speed of sound / p-wave velocity in air
    dt = 1.0e-15 #min(min(dx, dy) / c / 4.5, min(dx^2.0, dy^2.0) / ((4.0 / 3.0) * η / ρ) / 4.5)        # time step size                       

    # Equilibrium pressure distribution in the air
    #ρ[radrho .< diam ./ 2.0] .= ρa_non .+ ρc_non                                                          # initial density in the circle
    P .= P_non .* exp.(-g_non .* (y2dc_non .+ 1.0 ./ 5.0 .* L_non) .* M_non ./ T_non ./ R_non)
    # Pressure pertubation for the shockwave
    P .+= P_non .+ exp.((.-0.005 .* ((xc_non)/ 0.01).^2.0) .+ (.-0.005 .* ((yc_non .+ 0.35) / 0.01)'.^2.0))          # initial pressure distribution

    # Setting up matrices for density, compressibility, viscosity and shear modulus
    @. ρ    += ρs_non #* (maskrho_air == 1.0) #+ (maskrho_solid  == 1.0) * ρs_non                 # density
    @. ρ_vx += ρs_non #* (maskVx_air  == 1.0) #+ (maskVx_solid   == 1.0) * ρs_non                 # density
    @. ρ_vy += ρs_non #* (maskVy_air  == 1.0) #+ (maskVy_solid   == 1.0) * ρs_non                 # density
    @. β    += βs_non #* (maskrho_air == 1.0) #+ (maskrho_solid  == 1.0) * βs_non               # compressibility
    @. η    += ηs_non #* (maskrho_air == 1.0) #+ (maskrho_solid  == 1.0) * ηs_non               # viscosity
    @. η_c  += ηs_non #* (maskμc_air  == 1.0) #+ (maskμc_solid   == 1.0) * ηs_non               # viscosity
    @. μ    += μs_non #* (maskrho_air == 1.0) #+ (maskrho_solid  == 1.0) * μs_non               # shear modulus
    @. μ_c  += μs_non #* (maskμc_air  == 1.0) #+ (maskμc_solid   == 1.0) * μs_non               # shear modulus


    #xc_vec = Vector(xc)
    #yc_vec = Vector(yc)

    # Initial plotting
    # Points for the scatter plot showing the boundary between air and solid
    x_circ = zeros(length(x_circ_ind))
    y_circ = zeros(length(y_circ_ind))
    x_circ = xc[x_circ_ind]
    y_circ = yc[y_circ_ind]

    # Dimensionalization for plotting
    P_dim  = ustrip(dimensionalize(P, Pa, CharDim))
    τyy_dim = ustrip(dimensionalize(τyy, Pa, CharDim))
    xc_dim = ustrip(dimensionalize(xc_non, m, CharDim))
    yc_dim = ustrip(dimensionalize(yc_non, m, CharDim))

    # Ticks for the Axes
    xyticks = ["-400", "-200", "0", "200", "400"]
    xytick_positions = range(start=-400, stop=400, step=200)
    
    # Initial plotting
    fig = Figure()
    ax = Axis(fig[1,1][1,1], xticks=(xytick_positions, xyticks), yticks=(xytick_positions, xyticks),
                yticklabelsize=25, xticklabelsize=25, xlabelsize=25, ylabelsize=25)#, title="t = $t")#, aspect = DataAspect())#, limits=(nothing, nothing, nothing, 1.1))
    hm = heatmap!(ax, xc_dim, yc_dim, P_dim, colormap=Reverse(:roma))
    scatter!(ax, x_circ, y_circ, color=:white, markersize=4.0)
    Colorbar(fig[1,2][1,1], hm, label="Pressure", labelsize=25, ticklabelsize=25)#, vertical=false)
    display(fig)
    #save("../Plots/visco_elastic_shockwave_test/0.png", fig)

    # Time loop
    for i = 1:nt
        #β = 1.0 ./ P
        divV .= diff(Vx, dims=1) ./ dx_non .+ diff(Vy, dims=2) ./ dy_non
        dPdt .= .-(1.0 ./ β) .* divV
        P .= P .+ dPdt .* dt
        εxx .= diff(Vx, dims=1) ./ dx_non .- (1.0 ./ 3.0) .* divV
        εyy .= diff(Vy, dims=2) ./ dy_non .- (1.0 ./ 3.0) .* divV
        εxy[2:end-1, 2:end-1] .= 0.5 .* (diff(Vx[2:end-1, :], dims=2) ./ dy_non .+ diff(Vy[:, 2:end-1], dims=1) ./ dx_non)
        dτxxdt .= 2.0 .* μ .* εxx .- (μ ./ η) .* τxx
        dτyydt .= 2.0 .* μ .* εyy .- (μ ./ η) .* τyy
        dτxydt[2:end-1, 2:end-1] .= 2.0 .* μ_c[2:end-1, 2:end-1] .* εxy[2:end-1, 2:end-1] .- (μ_c[2:end-1, 2:end-1] ./ η_c[2:end-1, 2:end-1]) .* τxy[2:end-1, 2:end-1]
        τxx .= τxx .+ dτxxdt .* dt
        τyy .= τyy .+ dτyydt .* dt
        τxy .= τxy .+ dτxydt .* dt 
        dVxdt[2:end-1, :] .= (1.0 ./ av_x(ρ)) .* (.-diff(P, dims=1) ./ dx_non .+ diff(τxx, dims=1) ./ dx_non .+ diff(τxy[2:end-1, :], dims=2) ./ dy_non)                              .- Vx[2:end-1, :] .* diff(av_x(Vx), dims=1) ./ dx_non .+ av_x(av_y(Vy)) .* diff(av_x(Vx), dims=1) ./ dy_non
        dVydt[:, 2:end-1] .= (1.0 ./ av_y(ρ)) .* (.-diff(P, dims=2) ./ dy_non .+ diff(τyy, dims=2) ./ dy_non .+ diff(τxy[:, 2:end-1], dims=1) ./ dx_non .+ g_non .* av_y(ρ)) .- Vy[:, 2:end-1] .* diff(av_y(Vy), dims=2) ./ dy_non .+ av_y(av_x(Vx)) .* diff(av_y(Vy), dims=2) ./ dx_non
        #@infiltrate
        Vx[2:end-1, :] .= Vx[2:end-1, :] .+ dVxdt[2:end-1, :] .* dt .* maskVx_solid[2:end-1, :]
        Vy[:, 2:end-1] .= Vy[:, 2:end-1] .+ dVydt[: ,2:end-1] .* dt .* maskVy_solid[:, 2:end-1]
         
        Vx[1, :]   .= 0.0 #Vx[:, 2]
        Vx[end, :] .= 0.0 #Vx[:, end-1]
        Vy[1, :]   .= 0.0
        Vy[end, :] .= 0.0

        t += dt
        
        # Plotting for the time steps
        if i % 1 == 0
            # Dimensionalization for plotting
            Vy_dim = dimensionalize(Vy, m/s, CharDim)
            Vx_dim = dimensionalize(Vx, m/s, CharDim)
            Vx_av = av_x(ustrip(Vx_dim))
            Vy_av = av_y(ustrip(Vy_dim))
            data_plt = sqrt.(Vx_av.^2.0 .+ Vy_av.^2.0)
            P_plt = ustrip(dimensionalize(P, Pa, CharDim))
            t_dim = dimensionalize(t, s, CharDim)

            # Plotting
            fig1 = Figure()#size=(1000,1000))
            #ax1 = Axis(fig1[1,1], xticks=(xytick_positions, xyticks), yticks=(xytick_positions, xyticks),
            #        yticklabelsize=25, xticklabelsize=25, xlabelsize=25, ylabelsize=25, title="time = $t_dim")
            ax2 = Axis(fig1[1,1], xticks=(xytick_positions, xyticks), yticks=(xytick_positions, xyticks),
                    yticklabelsize=25, xticklabelsize=25, xlabelsize=25, ylabelsize=25, title="time = $t_dim")

            #lines!(ax2, xc_vec, P[:, Int(nx/2)])
            #hm1 = heatmap!(ax1, xc_dim, yc_dim, P_plt, colormap=Reverse(:roma))#, colorrange=(P0, maximum(P)))
            hm2 = heatmap!(ax2, xc_dim, yc_dim, data_plt, colormap=Reverse(:roma))#, colorrange=(0, 1500))#, colorrange=(0.0, 1.0))
            #scatter!(ax1, x_circ, y_circ, color=:white, markersize=4.0)
            scatter!(ax2, x_circ, y_circ, color=:white, markersize=4.0)
            #Colorbar(fig1[1,2], hm1, label="Pressure [Pa]", labelsize=25, ticklabelsize=25)#, vertical=false)
            Colorbar(fig1[1,2], hm2, label="Velocity [m/s]", labelsize=25, ticklabelsize=25)#, vertical=false)
            display(fig1)
            #save("../Plots/visco_elastic_shockwave_test/$i.png", fig1)

        end
    end
end