using CairoMakie
using Infiltrator
using GeoParams

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
    
    return inpts, onpts
end

function cross2d(a, b, c, d)
    return a*d - b*c
end

function av_x(B)
    B = 0.5 .* (B[2:end,:] .+ B[1:end-1,:])                 # try overriding the values of B, dont put them into a new array
end

function av_y(B)
    B = 0.5 .* (B[:,2:end] .+ B[:,1:end-1])
end

function av_xy(B)
    B = 0.25 .* (B[2:end,2:end] .+ B[1:end-1,2:end] .+ B[2:end,1:end-1] .+ B[1:end-1,1:end-1])
end

function advection2d(A, Vx, Vy, dx_non, dy_non, dt)
    dtadv = dt
    A[1:end-1,:] = A[1:end-1,:] .- (Vx[2:end-1,:] .< 0.0) .* Vx[2:end-1,:] .* diff(A, dims=1) ./ dx_non .* dtadv
    A[2:end,:] = A[2:end,:] .- (Vx[2:end-1,:] .> 0.0) .* Vx[2:end-1,:] .* diff(A, dims=1) ./ dx_non .* dtadv
    A[:,1:end-1] = A[:,1:end-1] .- (Vy[:,2:end-1] .< 0.0) .* Vy[:,2:end-1] .* diff(A, dims=2) ./ dy_non .* d
    tadv
    A[:,2:end] = A[:,2:end] .- (Vy[:,2:end-1] .> 0.0) .* Vy[:,2:end-1] .* diff(A, dims=2) ./ dy_non .* dtadv
    return A
end

function update_totalstresses_v(σxx, σyy, σxy, P, η, η_c, εxx, εyy, εxy)
    σxx .= .-P[2:end-1,2:end-1] .+ 2.0 .* η .* εxx                                                                          # total stress (dani class 5 equation)
    σyy .= .-P[2:end-1,2:end-1] .+ 2.0 .* η .* εyy                                                                          # total stress
    σxy .=                         2.0 .* η_c .* εxy                                                                            # total stress
    return σxx, σyy, σxy
end

function update_totalstresses_ve(σxx, σyy, σxy, σxx_old, σyy_old, σxy_old, P, P_old, η, η_c, μ, μ_c, εxx, εyy, εxy, dt)
    σxx .= (.-P[2:end-1,2:end-1] .+ 2.0 .* (1.0 ./ ((1.0 ./ η) .+ (1.0 ./ (μ .* dt)))) .* (εxx .+ ((σxx_old .+ P_old[2:end-1, 2:end-1]) ./ (2.0 .* μ .* dt)))) #.* maskrho_solid[2:end-1, 2:end-1] #.*μ.*  Exx .- (μ   ./ η) .*   (Sxx_old .+ P[2:end-1,2:end-1])                                                                          # stress in x-direction
    σyy .= (.-P[2:end-1,2:end-1] .+ 2.0 .* (1.0 ./ ((1.0 ./ η) .+ (1.0 ./ (μ .* dt)))) .* (εyy .+ ((σyy_old .+ P_old[2:end-1, 2:end-1]) ./ (2.0 .* μ .* dt)))) #.* maskrho_solid[2:end-1, 2:end-1] #.*μ.*  Eyy .- (μ   ./ η) .*   (Syy_old .+ P[2:end-1,2:end-1])                                             # stress in y-direction
    σxy .= (                        2.0 .* (1.0 ./ ((1.0 ./ η_c) .+ (1.0 ./ (μ_c .* dt)))) .* (εxy .+ (σxy_old                          ./ (2.0 .* μ_c .* dt)))) #.* maskμc_solid #.*μ_c.*Exy .- (μ_c ./ η_c) .* (Sxy_old .+ av_xy(P))                                                  # stress in xy-direction                                                                        # total stress
    return σxx, σyy, σxy
end

function conservative2D_ve()
    # Physics
    Lx      = 1.0e3                             # length of domain in x-direction
    Ly      = Lx                                # length of domain in y-direction
    rho0    = 1.225                               # density of air
    drho    = 2.7e3                               # density magma chamber
    rho_solid    = 2.8e3                        # density of solid
    Vx0     = 0.0                               # starting velocity in x-direction
    P0      = 1.0e5                               # pressure at rest
    beta    = 1.0/141.0e3                       # compressibility
    beta_air    = 1.0e-6                         # compressibility
    eta         = 1.81e-5                          # dynamic viscosity
    eta_air     = 1.0e-5                          # dynamic viscosity for air
    eta_solid   = 1.0e19                         # dynamic viscosity for solid
    #mu         = 1.81e-5                          # shear modulus
    mu_air      = 1.0e13                           # shear modulus for air
    mu_solid    = 1.0e11                           # shear modulus for solid
    g_y         = 9.81                             # gravitational acceleration
    a           = 40.0                             # wavelength content of excitation
    ν           = 0.4                            # Poisson's ratio
    K_s         = (2.0 * mu_solid * (1.0 + ν)) / (3.0 * (1.0 - 2.0 * ν))                # model poisson ratio of the solid
    beta_solid  = 1.0 / K_s                        # compressibility


    ν_solid = 0.4
    K_solid = (2.0 * mu_solid*(1.0+ν_solid)) / (3.0*(1.0 - 2.0 * ν_solid))                # model poisson ratio of the solid
    β_solid = 1.0 / K_solid                     # bulk modulus of solid
    Erde_ν  = 1.0e21 / (2.0 * (1.0e21 + 1.0e10)) # realistic poisson ratio of the earth


    # Numerics
    nx      = 101                               # number of nodes in x-direction
    ny      = 101                               # number of nodes in y-direction
    dx      = Lx/(nx)                           # grid spacing in x-direction
    dy      = Ly/(ny)                           # grid spacing in y-direction
    dt      = 1.0                               # time step size
    t    = 0.0                               # initial time
    t0      = 0.0                               # time of excitation      
    t_d     = 0.05                              # duration of excitation
    CFL_P   = 1.0/16.0                           # Courant number for pressure
    CFL_V   = 1.0/16.0                          # Courant number for velocity
    #psi     = 15.05 #5.25                               # dampening factor for pseudo transient iteration
    ksi     = 0.0 * 0.95 # 1.0 - (psi / nx)        # relaxation factor for stress when nx=ny
    save_plt    = false                             # save plots false or true
    display_err = false                         # display error false or true
    path    = "./Plots/test"                    # path to save plots

    # Nondimensionalization
    CharDim = GEO_units(length=1km, viscosity=1.0e19Pas, stress=100MPa)

    L       = 500.0                             # length of half the domain
    R       = 8.314                             # gas constant
    T       = 288.0                             # temperature
    M       = 0.028                             # molar mass of air

    a_geo = a*m
    t_geo = t*s
    t0_geo = t0*s
    dx_geo = dx*m
    dy_geo = dy*m
    L_geo = L*m
    R_geo = R*J/mol/K
    T_geo = T*K
    M_geo = M*kg
    ρa_geo = rho0*kg/m^3
    ρc_geo = drho*kg/m^3
    ρs_geo = rho_solid*kg/m^3
    βa_geo = beta_air*1/Pa
    βs_geo = beta_solid*1/Pa
    ηa_geo = eta_air*Pa*s
    ηs_geo = eta_solid*Pa*s
    μa_geo = mu_air*Pa
    μs_geo = mu_solid*Pa
    g_geo  = g_y*m/s^2
    P_geo  = P0*Pa

    a_non = nondimensionalize(a_geo, CharDim)
    t_non = nondimensionalize(t_geo, CharDim)
    t0_non = nondimensionalize(t0_geo, CharDim)
    dx_non = nondimensionalize(dx_geo, CharDim)
    dy_non = nondimensionalize(dy_geo, CharDim)
    L_non = nondimensionalize(L_geo, CharDim)
    R_non = nondimensionalize(R_geo, CharDim)
    T_non = nondimensionalize(T_geo, CharDim)
    M_non = nondimensionalize(M_geo, CharDim)
    ρa_non = nondimensionalize(ρa_geo, CharDim)
    ρc_non = nondimensionalize(ρc_geo, CharDim)
    ρs_non = nondimensionalize(ρs_geo, CharDim)
    βa_non = nondimensionalize(βa_geo, CharDim)
    βs_non = nondimensionalize(βs_geo, CharDim)
    ηa_non = nondimensionalize(ηa_geo, CharDim)
    ηs_non = nondimensionalize(ηs_geo, CharDim)
    μa_non = nondimensionalize(μa_geo, CharDim)
    μs_non = nondimensionalize(μs_geo, CharDim)
    g_non  = nondimensionalize(g_geo, CharDim)
    P_non  = nondimensionalize(P_geo, CharDim)

    # Reduce allocations
    β       = ones(Float64, nx + 2, ny + 2) #.* βs_non
    μ       = ones(Float64, nx, ny) #.* μs_non
    μ_c     = ones(Float64, nx + 1, ny + 1) #.* μs_non
    η       = ones(Float64, nx, ny) #.* ηs_non
    η_c     = ones(Float64, nx + 1, ny + 1) #.* ηs_non
    rho_old = zeros(Float64, nx, ny)
    Mx_old  = zeros(Float64, nx + 1, ny + 2)
    My_old  = zeros(Float64, nx + 2, ny + 1)
    c_loc   = zeros(Float64, nx, ny)
    dtPT    = zeros(Float64, nx, ny)
    dtrho   = zeros(Float64, nx, ny)
    Frhox   = zeros(Float64, nx + 1, ny + 2)
    Frhoy   = zeros(Float64, nx + 2, ny + 1)
    drhodt  = zeros(Float64, nx + 2, ny + 2)
    rhoRes  = zeros(Float64, nx, ny)
    rho     = ones(Float64, nx + 2, ny + 2) #.* ρs_non
    P       = zeros(Float64, nx + 2, ny + 2)
    P_dim   = zeros(Float64, nx + 2, ny + 2)
    P_old   = zeros(Float64, nx + 2, ny + 2)
    divV    = zeros(Float64, nx, ny)
    Exx     = zeros(Float64, nx, ny)
    Eyy     = zeros(Float64, nx, ny)
    Ezz     = zeros(Float64, nx, ny)
    Exy     = zeros(Float64, nx + 1, ny + 1)
    Sxx     = zeros(Float64, nx, ny)
    Sxx_old = zeros(Float64, nx, ny)
    Syy     = zeros(Float64, nx, ny)
    Syy_old = zeros(Float64, nx, ny)
    Szz     = zeros(Float64, nx, ny)
    Szz_old = zeros(Float64, nx, ny)
    Sxy     = zeros(Float64, nx + 1, ny + 1)
    Sxy_old = zeros(Float64, nx + 1, ny + 1)
    Mx      = zeros(Float64, nx + 1, ny + 2)
    FMxx    = zeros(Float64, nx, ny)
    FMxy    = zeros(Float64, nx - 1, ny + 1)
    dMxdt   = zeros(Float64, nx + 1, ny + 2)
    MxRes   = zeros(Float64, nx - 1, ny)
    Vx      = zeros(Float64, nx + 1, ny + 2)
    Vx_dim  = zeros(Float64, nx + 1, ny + 2)
    My      = zeros(Float64, nx + 2, ny + 1)
    FMyy    = zeros(Float64, nx, ny)
    FMyx    = zeros(Float64, nx + 1, ny - 1)
    dMydt   = zeros(Float64, nx + 2, ny + 1)
    MyRes   = zeros(Float64, nx, ny - 1)
    Vy      = zeros(Float64, nx + 2, ny + 1)
    Vy_dim  = zeros(Float64, nx + 2, ny + 1)
    rho_old = zeros(Float64, nx + 2, ny + 2)
    Mx_old  = zeros(Float64, nx + 1, ny + 2)
    My_old  = zeros(Float64, nx + 2, ny + 1)
    dMxdtau = zeros(Float64, nx - 1, ny)
    dMydtau = zeros(Float64, nx, ny - 1)
    radrho  = zeros(Float64, nx + 2, ny + 2)
    radVx   = zeros(Float64, nx + 1, ny + 2)
    radVy   = zeros(Float64, nx + 2, ny + 1)
    radμc   = zeros(Float64, nx + 1, ny + 1)
    beta_vec= zeros(Float64, nx + 2, ny + 2)
    maskrho_air   = zeros(Float64, nx, ny)
    maskVx_air    = zeros(Float64, nx, ny)
    maskVy_air    = zeros(Float64, nx, ny)
    maskμc_air    = zeros(Float64, nx + 1, ny + 1)

    maskrho_solid     = zeros(Float64, nx, ny)
    maskVx_solid      = zeros(Float64, nx, ny)
    maskVy_solid      = zeros(Float64, nx, ny)
    maskμc_solid      = zeros(Float64, nx + 1, ny + 1)

    maskrho_closedchamber_air = zeros(Float64, nx + 2, ny + 2) 
    maskrho_closedchamber_solid = zeros(Float64, nx + 2, ny + 2)
    maskμc_closedchamber_air = zeros(Float64, nx + 1, ny + 1) 
    maskμc_closedchamber_solid = zeros(Float64, nx + 1, ny + 1)

    # Initialization
    Xv      = range(-Lx/2.0, stop=Lx/2.0, step=dx)                      # x-coordinates of velocity nodes
    Xc      = range(-(Lx+dx)/2.0, stop=(Lx+dx)/2.0, step=dx)            # x-coordinates of density nodes
    Yv      = range(-Ly/2.0, stop=Ly/2.0, step=dy)                      # y-coordinates of velocity nodes
    Yc      = range(-(Ly+dy)/2.0, stop=(Ly+dy)/2.0, step=dy)            # y-coordinates of density nodes
    
    # Nondimensional coordinates
    Xv_non  = range(-L_non, stop=L_non, length=nx+1)                    # x-coordinates of velocity nodes
    Xc_non  = range(-(L_non+dx_non), stop=(L_non+dx_non), length=nx+2)  # x-coordinates of density nodes
    Yv_non  = range(-L_non, stop=L_non, length=ny+1)                    # y-coordinates of velocity nodes
    Yc_non  = range(-(L_non+dy_non), stop=(L_non+dy_non), length=ny+2)  # y-coordinates of density nodes
    
    x2dc, y2dc   = meshgrid(Xc,Yc)                                      # 2d mesh of x- and y-coordinates of density nodes
    x2dVx, y2dVx = meshgrid(Xv,Yc)                                      # 2d mesh of x- and y-coordinates of velocity nodes in x-direction
    x2dVy, y2dVy = meshgrid(Xc,Yv)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction
    x2dμc, y2dμc = meshgrid(Xv,Yv)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction

    x2dc_non, y2dc_non   = meshgrid(Xc_non,Yc_non)                                      # 2d mesh of x- and y-coordinates of density nodes
    x2dVx_non, y2dVx_non = meshgrid(Xv_non,Yc_non)                                      # 2d mesh of x- and y-coordinates of velocity nodes in x-direction
    x2dVy_non, y2dVy_non = meshgrid(Xc_non,Yv_non)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction
    x2dμc_non, y2dμc_non = meshgrid(Xv_non,Yv_non)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction

    Vx        .= Vx0.*ones(nx + 1,ny + 2)                               # initial velocity in x-direction
    Vx[1,:]   .= Vx0 
    Vx[end,:] .= Vx0
    locX      = 0.0                                                     # x-coordinate of circle
    locY      = -0.35*Ly                                                # y-coordinate of circle
    diam      = 0.2*Lx                                                  # diameter of circle
    radrho    .= sqrt.((x2dc .-locX).^2.0 .+ (y2dc .-locY).^2.0)        # radius of circle
    radVx     .= sqrt.((x2dVx.-locX).^2.0 .+ (y2dVx.-locY).^2.0)        # radius of circle
    radVy     .= sqrt.((x2dVy.-locX).^2.0 .+ (y2dVy.-locY).^2.0)        # radius of circle
    radμc     .= sqrt.((x2dμc.-locX).^2.0 .+ (y2dμc.-locY).^2.0)        # radius of circle
    Xp        = [-1.0/2.0, -1.0/2.0, -0.3, -1.0/8.0, -0.01, -0.01, 0.01, 0.01, 1.0/8.0, 0.3, 1.0/2.0, 1.0/2.0].*Lx          # x-coordinates of polygon
    Yp        = [-1.0/2.0, -1.0/5.0, -1.0/5.0, -0.1, -0.15,  -0.28, -0.28, -0.15,  -0.1, -1.0/5.0, -1.0/5.0, -1.0/2.0].*Ly  # y-coordinates of polygon
    Xp_flat   = [-1.0/1.9, -1.0/2.0,  1.0/2.0, 1.0/1.9].*Lx          # x-coordinates of polygon, 1.0/1.9 damit Ränder mit zum Polygon gehören und nicht Luft sind
    Yp_flat   = [-1.0/1.9, -1.0/5.0, -1.0/5.0, -1.0/1.9].*Ly  # y-coordinates of polygon, 1.0/1.9 damit Ränder mit zum Polygon gehören und nicht Luft sind
    Xp2       = [-0.02, -0.02,  0.02,  0.02].*Lx                                                                            # x-coordinates of polygon
    Yp2       = [-0.14,  -0.25, -0.25, -0.14].*Ly                                                                             # y-coordinates of polygon

    inpolyrho, onpolyrho   = inpolygon(x2dc,y2dc,Xp_flat,Yp_flat)                                 # 1: inside polygon, 0: outside polygon
    inpolyVx, onpolyVx    = inpolygon(x2dVx,y2dVx,Xp_flat,Yp_flat)                                 # 1: inside polygon, 0: outside polygon
    inpolyVy, onpolyVy    = inpolygon(x2dVy,y2dVy,Xp_flat,Yp_flat)                                 # 1: inside polygon, 0: outside polygon
    inpolyμc, onpolyμc    = inpolygon(x2dμc,y2dμc,Xp_flat,Yp_flat)                                 # 1: inside polygon, 0: outside polygon

    inpolyrho_cond, onpolyrho_cond = inpolygon(x2dc,y2dc,Xp2,Yp2)
    inpolyVx_cond, onpolyVx_cond  = inpolygon(x2dVx,y2dVx,Xp2,Yp2)
    inpolyVy_cond, onpolyVy_cond  = inpolygon(x2dVy,y2dVy,Xp2,Yp2)
    inpolyμc_cond, onpolyμc_cond  = inpolygon(x2dμc,y2dμc,Xp2,Yp2)

    maskrho_air = 1.0 .- inpolyrho                                             # mask for density
    maskμc_air = 1.0 .- inpolyμc                                             # mask for density
    maskVx_air  = 1.0 .- inpolyVx                               # mask for velocity in x-direction
    maskVy_air  = 1.0 .- inpolyVy                               # mask for velocity in y-direction

    maskrho_solid = 1.0 .- maskrho_air
    maskμc_solid  = 1.0 .- maskμc_air
    maskVx_solid  = 1.0 .- maskVx_air
    maskVy_solid  = 1.0 .- maskVy_air

    # mask for closed conduit and chamber completely surrounded by solid
    #maskrho_closedchamber_solid .= maskrho_solid
    #maskrho_closedchamber_solid[inpolyrho_cond .== 1.0] .= 1.0
    #maskrho_closedchamber_air = 1.0 .- maskrho_closedchamber_solid
    #maskμc_closedchamber_solid .= maskμc_solid
    #maskμc_closedchamber_solid[inpolyμc_cond .== 1.0] .= 1.0
    #maskμc_closedchamber_air = 1.0 .- maskμc_closedchamber_solid

    # mask for flat topography 
    maskrho_closedchamber_solid .= maskrho_solid
    maskrho_closedchamber_air = 1.0 .- maskrho_closedchamber_solid
    maskμc_closedchamber_solid .= maskμc_solid
    maskμc_closedchamber_air = 1.0 .- maskμc_closedchamber_solid

    x_circ_ind = []
    y_circ_ind = []

    for y in 2:size(maskrho_closedchamber_air)[1] - 1, x in 2:size(maskrho_closedchamber_air)[2] - 1
        if (maskrho_closedchamber_air[x, y] == 0.0 && abs(maskrho_closedchamber_air[x, y] - maskrho_closedchamber_air[x - 1, y]) == 1.0)
            push!(x_circ_ind, x)
            push!(y_circ_ind, y)
        elseif (maskrho_closedchamber_air[x, y] == 0.0 && abs(maskrho_closedchamber_air[x, y] - maskrho_closedchamber_air[x + 1, y]) == 1.0 )
            push!(x_circ_ind, x)
            push!(y_circ_ind, y)
        elseif (maskrho_closedchamber_air[x, y] == 0.0 && abs(maskrho_closedchamber_air[x, y] - maskrho_closedchamber_air[x, y - 1]) == 1.0)
            push!(x_circ_ind, x)
            push!(y_circ_ind, y)
        elseif (maskrho_closedchamber_air[x, y] == 0.0 && abs(maskrho_closedchamber_air[x, y] - maskrho_closedchamber_air[x, y + 1]) == 1.0)
            push!(x_circ_ind, x)
            push!(y_circ_ind, y)
        end
    end


    # Initial conditions
    @. β      = βa_non * (maskrho_closedchamber_air == 1.0) + βs_non * (maskrho_closedchamber_solid == 1.0)                  # initial viscosity distribution  

    P         .= (P_non.*exp.(-g_non.*(y2dc_non.+1.0./5.0.*L_non).*M_non./T_non./R_non)) .* maskrho_closedchamber_air .+ reverse(cumsum(ρs_non.*g_non.*reverse(maskrho_solid, dims=2).*dy_non,dims=2), dims=2)                                # barometric setting atmosphere: P = P0*exp(-(g*(h-h0)*M)/(T*R)) M: Mass density of air, T: Temperature, R: Gas constant, h: height, P0: Pressure at sea level
    rho       .= ρa_non .* maskrho_closedchamber_air .+ ρs_non .* maskrho_closedchamber_solid                         # equation of state for density depending on pressure
    rho[radrho .< diam ./ 2.0] .= ρa_non.+ρc_non                                                          # initial density in the circle
    #rho[inpolygon(x2dc,y2dc,Xp_flat,Yp_flat) .== 1.0 .&& y2dc.< locY.+diam./2.0] .= ρa_non.+ ρc_non          # initial density of the conduit overlapping with the circular chamber (so in the chamber)
    #rho[inpolygon(x2dc,y2dc,Xp_flat,Yp_flat) .== 1.0 .&& y2dc.>=locY.+diam./2.0] .= .-0.0 .* (y2dc[inpolygon(x2dc,y2dc,Xp2,Yp2).==1.0 .&& y2dc.>=locY.+diam./2.0].+0.15.*L_non).*ρc_non./(-0.15.-(locY.+diam./2.0)).+ρa_non    # initial density in the conduit

    #P .+= reverse(cumsum(ρs_non*g_non.*reverse(maskrho_solid, dims=2)*dx_non,dims=2), dims=2)
    #P .= P_non .+ exp.((.-50.0 .* Xc_non.^2.0) .+ (.-50.0 .* Yc_non'.^2.0))          # initial pressure distribution

    P[radrho .< diam ./ 2.0]         .= (P_non .* rho[radrho .< diam ./ 2.0]) ./ ρa_non #.+ (rho[radrho .< diam ./ 2.0] .* g .* depth[radrho .< diam ./ 2.0])#P0 ./ rho0 .* rho
    #P         .+=(P0 .+  (1.0 ./ beta) .* log.(rho ./ rho_solid)) .* maskrho_solid


    # Initial parameter matrices for both phases
    @. η += ηa_non * (maskrho_closedchamber_air[2:end-1, 2:end-1] == 1.0) + ηs_non * (maskrho_closedchamber_solid[2:end-1, 2:end-1] == 1.0)     # initial viscosity distribution
    @. μ += μa_non * (maskrho_closedchamber_air[2:end-1, 2:end-1] == 1.0) + μs_non * (maskrho_closedchamber_solid[2:end-1, 2:end-1] == 1.0)       # initial viscosity distribution

    @. η_c += ηa_non * (maskμc_closedchamber_air == 1.0) + ηs_non * (maskμc_closedchamber_solid == 1.0)   # initial viscosity distribution for corner nodes
    @. μ_c += μa_non * (maskμc_closedchamber_air == 1.0) + μs_non * (maskμc_closedchamber_solid == 1.0)     # initial viscosity distribution for corner nodes

    # Inital plot 
    P_dim = ustrip(dimensionalize(P, Pa, CharDim))
    
    fig = Figure()
    ax = Axis(fig[1, 1], xticks=([-500, 0, 500], ["-500", "0", "500"]), yticks=([-500, 0, 500], ["-500", "0", "500"]),
                yticklabelsize=25, xticklabelsize=25, xlabelsize=25, ylabelsize=25)#, title="0")

    x_circ = zeros(length(x_circ_ind))
    y_circ = zeros(length(y_circ_ind))
    x_circ = Xc[x_circ_ind]
    y_circ = Yc[y_circ_ind]

    U = av_y(Vx)
    V = av_x(Vy)
    data_plt = sqrt(U.^2.0 .+ V.^2.0)
    
    p1 = heatmap!(ax, Xc, Yc, P_dim, colormap=Reverse(:roma))
    #p1 = heatmap!(ax, x2dc, y2dc , data_plt, shading=false, colorrange=(0.0, 350.0), colormap=Reverse(:roma))
    #p1 = heatmap!(ax, x2dc, y2dc , rho, shading=false, colormap=Reverse(:roma))#, colorrange=(P0, P0*2))
    #Colorbar(fig[1, 2], p1, label="Velocity", labelsize=25, ticklabelsize=25)
    Colorbar(fig[1, 2], p1, label="Pressure ", labelsize=25, ticklabelsize=25)
    #points_XpYp = Point2f[]
    #=for i in eachindex(Xp)
        x = Xp[i]
        y = Yp[i]
        point = (x,y)
        push!(points_XpYp, point)
    end=#
    #poly!(ax,points_XpYp)
    #lines!(ax, Xp, Yp, color = :white) # only conduit
    #scatter!(ax, x_circ, y_circ, color = :white, markersize=4.0) # conduit and chamber
    display(fig)
    if save_plt
        mkdir(path)
        save(path * "/0.png", fig)
    end

    # Solver
    for it = 1:600                                                   # 60 iterations max for all equilibrated
        @show it
        P_old .= P
        rho_old .= rho                                              # save old density values for time integration
        Mx_old .= Mx                                                # save old momentum values in x-direction for time integration
        My_old .= My                                                # save old momentum values in y-direction for time integration
        Sxx_old .= Sxx
        Syy_old .= Syy
        Sxy_old .= Sxy
        Szz_old .= Szz
        err = 1.0                                                   # error for time integration
        iter = 0                                                    
        it_counter = 0
        dMxdtau .= 0.0                                              # stress derivative of momentum in x-direction
        dMydtau .= 0.0                                              # stress derivative of momentum in y-direction
        err_vec = zeros(Float64, 1000, 5) .* NaN
        while err > 1.0e-1 #&& iter < 25*360
            iter += 1
            #beta_vec .= 1.0 ./ P
            @. c_loc = 1.0/sqrt(rho[2:end-1,2:end-1] * β[2:end-1,2:end-1])                                                    # local speed of sound
            dt = minimum([dx_non./maximum(abs.(Vx)), dy_non./maximum(abs.(Vy)), min(dx_non, dy_non) .* sqrt(maximum(rho .* β))]).*4.1 # time step size  , , min(dx_non, dy_non).^2.0./maximum(η)
            dtPT .= min.(min.(min.(dx_non ./ abs.(av_x(Vx[:, 2:end-1])), dx_non ./ abs.(av_y(Vy[2:end-1, :]))), min(dx_non, dy_non).^2.0 ./ η), dx_non ./ c_loc) # time step size for pressure and temperature
            dtrho .= 1.0 ./ (1.0 ./ dt .+ 1.0 ./ (min(dx_non, dy_non) ./ c_loc ./ 4.1))                                                         # time step size for density

            # Conservation of mass
            Frhox .= (Vx .> 0.0).*Vx.*rho[1:end-1, :] .+ (Vx .< 0.0).*Vx.*rho[2:end, :] # mass flux in x-direction (upwind scheme)
            Frhoy .= (Vy .> 0.0).*Vy.*rho[:, 1:end-1] .+ (Vy .< 0.0).*Vy.*rho[:, 2:end] # mass flux in y-direction (upwind scheme)
            drhodt .= (rho .- rho_old)./dt                                                                           # time derivative of density
            rhoRes .= .-drhodt[2:end-1, 2:end-1] .- diff(Frhox[:, 2:end-1], dims=1)./dx_non .- diff(Frhoy[2:end-1, :], dims=2)./dy_non        # updating residual of density
            rho[2:end-1, 2:end-1] .= rho[2:end-1, 2:end-1] .+ rhoRes.*dtrho.*CFL_P                   # updating density

            #if iter >= 1
            #    @infiltrate
            #end

            # Boundary conditions Inflow and outflow densities
            #rho[1, :] .= rho[2, :]
            #rho[end, :] .= rho[end-1, :]
            # Boundary conditions impermeable walls
            #rho[:, 1] .= rho[:, 2]
            #rho[:, end] .= rho[:, end-1]

            # Strain-rates and stresses
            P = zeros(Float64, nx + 2, ny + 2)
            #@show findmin(rho)
            P              .+= ((1.0 ./ β) .* log.(rho ./ ρs_non) .+ P_non) #.* maskrho_solid             # equation of state for pressure depending on density
            #P              .+= ((P_non .* rho) ./ ρa_non ).* maskrho_air                                              # equation of state for pressure depending on density            
            divV         .= diff(Vx[:,2:end-1],dims=1)./dx_non .+ diff(Vy[2:end-1,:],dims=2)./dy_non                                             # divergence of velocity
            Exx          .= diff(Vx[:,2:end-1],dims=1)./dx_non .- 1.0./3.0.*divV                                                               # strain-rate in x-direction
            Eyy          .= diff(Vy[2:end-1,:],dims=2)./dy_non .- 1.0./3.0.*divV                                                               # strain-rate in y-direction
            Exy          .= 0.5.*(diff(Vy,dims=1)./dx_non .+ diff(Vx,dims=2)./dy_non) 
            Sxx, Syy, Sxy = update_totalstresses_v(Sxx, Syy, Sxy, P, η, η_c, Exx, Eyy, Exy)            
            #Sxx[10, 10]                        = (t_non - t0_non >= 0.0) * (-2.0 * a * (t_non - t0_non) * exp(-a * (t_non - t0_non)^2.0))
            #Syy[10, 10]                        = (t_non - t0_non >= 0.0) * (-2.0 * a * (t_non - t0_non) * exp(-a * (t_non - t0_non)^2.0))
            #Sxy[10, 10]                        = (t_non - t0_non >= 0.0) * (-2.0 * a * (t_non - t0_non) * exp(-a * (t_non - t0_non)^2.0))
            
            dtV           = 1.0 / (1.0 / dt + 1.0 / (min(dx_non,dy_non)^2.0 / maximum(η) / 4.0)) * CFL_V                                           # time step size for velocity
            
            # Conservation of the x-component of momentum
            Mx             .= av_x(rho).*Vx                                                             # momentum in x-direction
            FMxx         .= (av_x(Vx[ :     ,2:end-1]).> 0.0).*av_x(Vx[ :     ,2:end-1]).*Mx[1:end-1,2:end-1] .+ (av_x(Vx[ :     ,2:end-1]).< 0.0).*av_x(Vx[ :     ,2:end-1]).*Mx[2:end  ,2:end-1]  # mass flux in x-direction (upwind scheme)
            FMxy         .= (av_x(Vy[2:end-1, :     ]).> 0.0).*av_x(Vy[2:end-1, :     ]).*Mx[2:end-1,1:end-1] .+ (av_x(Vy[2:end-1, :     ]).< 0.0).*av_x(Vy[2:end-1, :     ]).*Mx[2:end-1,2:end  ]  # mass flux in y-direction (upwind scheme)
            dMxdt       .= (Mx.-Mx_old)./dt                                                                                             # time derivative of momentum in x-direction
            MxRes      .= .-dMxdt[2:end-1,2:end-1] .- diff((FMxx .- Sxx),dims=1)./dx_non .- diff(FMxy .- Sxy[2:end-1,:],dims=2)./dy_non       # updating residual of momentum in x-direction
            dMxdtau   .= MxRes .+ dMxdtau .* ksi                                                                                    # stress derivative of momentum in x-direction
            Mx[2:end-1,2:end-1]  .= Mx[2:end-1,2:end-1] .+ dMxdtau.*av_x(dtPT).*CFL_V                                                     # updating momentum in x-direction
            Vx           .= Mx./av_x(rho)                                                              # velocity in x-direction

            # BC fixed walls (normal velocity = 0)
            Mx[:,1]      .= Mx[:,2]
            Mx[:,end]    .= Mx[:,end-1]
            # BC no slip on vertical walls
            #Mx[1,:]      .= Mx[2,:]
            #Mx[end,:]    .= Mx[end-1,:]
            
            Vx[1  ,:] .= Vx0
            Vx[end,:] .= Vx0
            
            #if iter >= 1
            #    @infiltrate
            #end

            # Conservation of the y component of momentum
            My             .= av_y(rho).*Vy                                                              # momentum in y-direction
            FMyy         .= (av_y(Vy[2:end-1, :     ]).> 0.0).*av_y(Vy[2:end-1, :     ]).*My[2:end-1,1:end-1] .+ (av_y(Vy[2:end-1, :     ]).< 0.0).*av_y(Vy[2:end-1, :     ]).*My[2:end-1,2:end  ]  # mass flux in y-direction
            FMyx         .= (av_y(Vx[ :     ,2:end-1]).> 0.0).*av_y(Vx[ :     ,2:end-1]).*My[1:end-1,2:end-1] .+ (av_y(Vx[ :     ,2:end-1]).< 0.0).*av_y(Vx[ :     ,2:end-1]).*My[2:end  ,2:end-1]  # mass flux in x-direction
            dMydt       .= (My.-My_old)./dt                                                                                              # time derivative of momentum in y-direction
            MyRes      .= .-dMydt[2:end-1,2:end-1] .- diff(FMyy .- Syy,dims=2)./dy_non .- diff(FMyx .- Sxy[:,2:end-1],dims=1)./dx_non - g_non .* av_y(rho[2:end-1,2:end-1]) # drunken sailor chapter taras book: dt .* (av_xy(Vx[:, 2:end-1]) .* av_y(diff(rho[2:end,2:end-1], dims=1)) .+ Vy[2:end-1, 2:end-1] .* diff(rho[2:end-1,2:end-1], dims=2))
            dMydtau   .= MyRes .+ dMydtau.*ksi                                                                                      # stress derivative of momentum in y-direction
            My[2:end-1,2:end-1]  .= My[2:end-1,2:end-1] .+ dMydtau.*av_y(dtPT).*CFL_V                                                     # updating momentum in y-direction
            Vy             .= My./av_y(rho)                   # updating velocity in y-direction

            # BC fixed walls (normal velocity = 0)
            My[1,:]      .= -My[2,:]
            My[end,:]    .= -My[end-1,:]
            # BC no slip on horizontal walls
            My[:,1]      .= My[:,2]
            My[:,end]    .= My[:,end-1]

            if iter >= 17000
                @infiltrate
            end

            if mod(iter, 500) == 0
                it_counter += 1
                err = maximum(abs.([rhoRes[:]; MxRes[:]; MyRes[:]]))
                print("PT_iter = $iter, err = $err, rhoRes = $(maximum(abs.(rhoRes[:]))), MxRes = $(maximum(abs.(MxRes[:]))), MyRes = $(maximum(abs.(MyRes[:])))\n")
                if isnan(err) == 1 
                    break
                end
                err_vec[it_counter, 1] = it_counter
                err_vec[it_counter, 2] = err
                err_vec[it_counter, 3] = maximum(abs.(rhoRes[:]))
                err_vec[it_counter, 4] = maximum(abs.(MxRes[:]))
                err_vec[it_counter, 5] = maximum(abs.(MyRes[:]))
                figerr = Figure()
                axerr  = Axis(figerr[1,1], limits=(1, nothing, -6.1, 8.0), yticks=(Array(range(start=-6.0, stop=8.0, step=1.0)), ["-6.0","-5.0","-4.0","-3.0", "-2.0", "-1.0", "0.0", "1.0", "2.0", "3.0", "4.0", "5.0", "6.0", "7.0", "8.0"]),
                            yticklabelsize=25, xticklabelsize=25, xlabelsize=25, ylabelsize=25, title="PT #$iter error")
                sc2 = scatterlines!(axerr, err_vec[:, 1], log10.(err_vec[:, 3]), markersize=10.0, label="rhoRes")
                sc3 = scatterlines!(axerr, err_vec[:, 1], log10.(err_vec[:, 4]), markersize=10.0, label="MxRes")
                sc4 = scatterlines!(axerr, err_vec[:, 1], log10.(err_vec[:, 5]), markersize=10.0, label="MyRes")
                sc1 = scatterlines!(axerr, err_vec[:, 1], log10.(err_vec[:, 2]), markersize=10.0, label="err")
                l   = hlines!(axerr, -3.0, label="convergence", linestyle=:dash, color=:black)
                axislegend(axerr; position=:lb)
                if display_err
                    display(figerr)
                end

                if save_plt && err <= 1.0e-3
                    save(path * "/PT_error_$(it).png", figerr)
                elseif err <= 1.0e-3
                    #display(figerr)
                end
            end

            if err <= 1.0e-3
                @show it_counter
                print("------------------------------------------------------------------------------------\n")
            end
        end
        t = t + dt

        # Updating plot
        if mod(it-1, 1) == 0
            fig1 = Figure(resolution=(2000,2000))
            ax1 = Axis(fig1[1,1], xticks=([-500, 0, 500], ["-500", "0", "500"]), yticks=([-500, 0, 500], ["-500", "0", "500"]),
                    yticklabelsize=25, xticklabelsize=25, xlabelsize=25, ylabelsize=25, title="time = $t")
            ax2 = Axis(fig1[2,1], xticks=([-500, 0, 500], ["-500", "0", "500"]), yticks=([-500, 0, 500], ["-500", "0", "500"]),
                    yticklabelsize=25, xticklabelsize=25, xlabelsize=25, ylabelsize=25, title="time = $t")

            Vx_dim .= ustrip(dimensionalize(Vx, m/s, CharDim))
            Vy_dim .= ustrip(dimensionalize(Vy, m/s, CharDim))
            P_dim  .= ustrip(dimensionalize(P, Pa, CharDim))
            X = av_xy(x2dc)
            Y = av_xy(y2dc)
            U = av_y(Vx_dim)
            V = av_x(Vy_dim)

            data_plt = sqrt.(U.^2.0 .+ V.^2.0)

            hm = heatmap!(ax1, Xc, Yc, P_dim)# colorrange=(P0, P0*2))
            #hm = heatmap!(ax1, X[2:end-1, 2:end-1], Y[2:end-1, 2:end-1], data_plt[2:end-1, 2:end-1], shading=false, colormap=Reverse(:roma))#, colorrange=(0.0, 0.1),)
            hm2 = heatmap!(ax2, Xc[2:end-1], Yc[2:end-1], data_plt[2:end-1, 2:end-1], colormap=Reverse(:roma))#, colorrange=(0.0, 0.01),)
            Colorbar(fig1[1,2],  hm, label="Pressure", labelsize=25, ticklabelsize=25)
            #Colorbar(fig1[1,2],  hm, label="Velocity", labelsize=25, ticklabelsize=25)
            Colorbar(fig1[2,2],  hm2, label="Velocity", labelsize=25, ticklabelsize=25)
            #lines!(ax1, Xp, Yp, color = :white)    # only conduit
            scatter!(ax1, x_circ, y_circ, color = :white, markersize=4.0)     # conduit and chamber
            #poly!(ax,points_XpYp)
            
            #stepsize = 10
            #arrows!(ax1, X[1:stepsize:end, 1], Y[1, 1:stepsize:end], U[1:stepsize:end, 1:stepsize:end], V[1:stepsize:end, 1:stepsize:end], arrowsize=7, color = :white)
            display(fig1)

            if save_plt
                save(path * "/$(it).png", fig1)
            end
        end
    end
    return
end
