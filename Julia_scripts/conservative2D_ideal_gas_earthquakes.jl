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

function advection2d(A, Vx, Vy, dx, dy, dt)
    dtadv = dt
    A[1:end-1,:] = A[1:end-1,:] .- (Vx[2:end-1,:] .< 0.0) .* Vx[2:end-1,:] .* diff(A, dims=1) ./ dx .* dtadv
    A[2:end,:] = A[2:end,:] .- (Vx[2:end-1,:] .> 0.0) .* Vx[2:end-1,:] .* diff(A, dims=1) ./ dx .* dtadv
    A[:,1:end-1] = A[:,1:end-1] .- (Vy[:,2:end-1] .< 0.0) .* Vy[:,2:end-1] .* diff(A, dims=2) ./ dy .* dtadv
    A[:,2:end] = A[:,2:end] .- (Vy[:,2:end-1] .> 0.0) .* Vy[:,2:end-1] .* diff(A, dims=2) ./ dy .* dtadv
    return A
end

function update_densityflux(Frhox, Frhoy, Vx, Vy, rho)
    Frhox .= (Vx .> 0.0).*Vx.*rho[1:end-1, :] .+ (Vx .< 0.0).*Vx.*rho[2:end, :]
    Frhoy .= (Vy .> 0.0).*Vy.*rho[:, 1:end-1] .+ (Vy .< 0.0).*Vy.*rho[:, 2:end]
    return Frhox, Frhoy
end

function update_drhodt(drhodt, rho, rho_old, dt)
    drhodt .= (rho .- rho_old)./dt
    return drhodt
end

function update_densityResidual(rhoRes, drhodt, Frhox, Frhoy, dx_non, dy_non, mask)
    rhoRes .= .-drhodt[2:end-1, 2:end-1] .- diff(Frhox[:, 2:end-1], dims=1)./dx_non .- diff(Frhoy[2:end-1, :], dims=2)./dy_non        # updating residual of density
    rhoRes .= rhoRes# .* mask[2:end-1, 2:end-1]                                                                               # applying mask to residual of density
    return rhoRes
end

function update_density(rho, rhoRes, dtrho, CFL_P)
    rho[2:end-1, 2:end-1] .= rho[2:end-1, 2:end-1] .+ rhoRes .* dtrho .* CFL_P                   # updating density
    return rho
end

function update_pressure(P, P_non, ρa_non, ρ)
    P .= P_non ./ ρa_non .* ρ                                                  # equation of state for pressure depending on density
    return P
end

function update_divergenceV(divV, Vx, Vy, dx_non, dy_non)
    divV .= diff(Vx[:,2:end-1], dims=1)./dx_non .+ diff(Vy[2:end-1,:],dims=2)./dy_non                                             # divergence of velocity
    return divV
end

function update_strainrates(εxx, εyy, εxy, Vx, Vy, divV, dx_non, dy_non)
    εxx .= diff(Vx[:,2:end-1],dims=1)./dx_non .- 1.0./3.0.*divV                                                               # strain-rate in x-direction
    εyy .= diff(Vy[2:end-1,:],dims=2)./dy_non .- 1.0./3.0.*divV                                                               # strain-rate in y-direction
    εxy .= 0.5.*(diff(Vy,dims=1)./dx_non .+ diff(Vx,dims=2)./dy_non)                                                             # shear strain-rate in xy-direction
    return εxx, εyy, εxy
end

function update_totalstresses(σxx, σyy, σxy, P, η, εxx, εyy, εxy)
    σxx .= .-P[2:end-1,2:end-1] .+ 2.0 .* η .* εxx                                                                          # total stress (dani class 5 equation)
    σyy .= .-P[2:end-1,2:end-1] .+ 2.0 .* η .* εyy                                                                          # total stress
    σxy .=                         2.0 .* η .* εxy                                                                            # total stress
    return σxx, σyy, σxy
end

function update_advection_momentum(Mx, My, Vx, Vy, ρ)
    Mx .= av_x(ρ).*Vx                                                             # momentum in x-direction
    My .= av_y(ρ).*Vy                                                             # momentum in y-direction
    return Mx, My
end

function update_momentumflux(FMxx, FMyy, FMxy, FMyx, Vx, Vy, Mx, My, σxx, σyy, σxy, dx_non, dy_non)
    FMxx .= (av_x(Vx[ :     ,2:end-1]).> 0.0).*av_x(Vx[ :     ,2:end-1]).*Mx[1:end-1,2:end-1] .+ (av_x(Vx[ :     ,2:end-1]).< 0.0).*av_x(Vx[ :     ,2:end-1]).*Mx[2:end  ,2:end-1]  # upwind advective momentum flux
    FMxy .= (av_x(Vy[2:end-1, :     ]).> 0.0).*av_x(Vy[2:end-1, :     ]).*Mx[2:end-1,1:end-1] .+ (av_x(Vy[2:end-1, :     ]).< 0.0).*av_x(Vy[2:end-1, :     ]).*Mx[2:end-1,2:end  ]  # upwind advective momentum flux
    FMyx .= (av_y(Vx[ :     ,2:end-1]).> 0.0).*av_y(Vx[ :     ,2:end-1]).*My[1:end-1,2:end-1] .+ (av_y(Vx[ :     ,2:end-1]).< 0.0).*av_y(Vx[ :     ,2:end-1]).*My[2:end  ,2:end-1]  # upwind advective momentum flux
    FMyy .= (av_y(Vy[2:end-1, :     ]).> 0.0).*av_y(Vy[2:end-1, :     ]).*My[2:end-1,1:end-1] .+ (av_y(Vy[2:end-1, :     ]).< 0.0).*av_y(Vy[2:end-1, :     ]).*My[2:end-1,2:end  ]  # upwind advective momentum flux
    return FMxx, FMyy, FMxy, FMyx
end

function update_dMdt(dMxdt, dMydt, Mx, My, Mx_old, My_old, dt)
    dMxdt .= (Mx .- Mx_old)./dt
    dMydt .= (My .- My_old)./dt
    return dMxdt, dMydt
end

function update_momentumResidual(MxRes, MyRes, dMxdt, dMydt, FMxx, FMyy, FMxy, FMyx, σxx, σyy, σxy, dx_non, dy_non, g_non, rho, mask1, mask2)
    MxRes .= .-dMxdt[2:end-1,2:end-1] .- diff(FMxx .- σxx,dims=1)./dx_non .- diff(FMxy .- σxy[2:end-1,:],dims=2)./dx_non                                               # updating residual of momentum in y-direction
    MyRes .= .-dMydt[2:end-1,2:end-1] .- diff(FMyy .- σyy,dims=2)./dy_non .- diff(FMyx .- σxy[:,2:end-1],dims=1)./dx_non .- g_non .* av_y(rho[2:end-1,2:end-1])        # updating residual of momentum in x-direction
    MxRes .= MxRes# .* mask1[2:end-1, 2:end-1]                                                                               # applying mask to residual of momentum in x-direction
    MyRes .= MyRes# .* mask2[2:end-1, 2:end-1]                                                                               # applying mask to residual of momentum in y-direction
    return MxRes, MyRes
end

function update_dMdτ(dMxdτ, dMydτ, MxRes, MyRes, ksi)
    dMxdτ .= MxRes .+ dMxdτ .* ksi
    dMydτ .= MyRes .+ dMydτ .* ksi
    return dMxdτ, dMydτ
end

function update_inertia_momentum(Mx, My, dMxdτ, dMydτ, dtPT, CFL_V)
    Mx[2:end-1,2:end-1] .= Mx[2:end-1,2:end-1] .+ dMxdτ.*av_x(dtPT).*CFL_V
    My[2:end-1,2:end-1] .= My[2:end-1,2:end-1] .+ dMydτ.*av_y(dtPT).*CFL_V

    set_momentum_bc!(Mx, My)
    return Mx, My
end

function update_velocities(Vx, Vy, Mx, My, ρ)
    Vx .= Mx ./ av_x(ρ)
    Vy .= My ./ av_y(ρ)

    set_velocity_bc!(Vx)
    return Vx, Vy
end

function set_momentum_bc!(Mx, My)
    Mx[1,:]   .= Mx[2,:]
    Mx[end,:] .= Mx[end-1,:]
    Mx[:,1]   .= Mx[:,2]
    Mx[:,end] .= Mx[:,end-1]

    My[1,:]   .= My[2,:]
    My[end,:] .= My[end-1,:]
    My[:,1]   .= My[:,2]
    My[:,end] .= My[:,end-1]
    return Mx, My
end

function set_velocity_bc!(Vx)
    Vx[1,:] .= 0.0
    Vx[end,:] .= 0.0
    return Vx
end

function conservative2D_earthq()
    # Physics
    Lx      = 1.0e4                             # length of domain in x-direction
    Ly      = Lx                                # length of domain in y-direction
    rho0    = 2.8e3                             # density at rest
    drho    = 3.0e2                             # density difference
    P0      = 1.0e5                             # pressure at rest
    eta_air      = 1.81e-5                           # dynamic viscosity
    eta_solid    = 1.0e21                           # dynamic viscosity
    mu_solid    = 1.0e11                           # shear modulus
    g_y       = 9.81                              # gravitational acceleration
    t = 0.0

    # Numerics
    nx      = 141                               # number of nodes in x-direction
    ny      = 141                               # number of nodes in y-direction
    dx      = Lx/(nx)                           # grid spacing in x-direction
    dy      = Ly/(ny)                           # grid spacing in y-direction
    #nout    = 1
    dt      = 1.0                           # time step size
    CFL_P   = 1.0/16.0                           # Courant number for pressure
    CFL_V   = 1.0/16.0                          # Courant number for velocity
    ksi     = 0.0*0.95                          # relaxation factor for stress
    err_threshold = 1.0e-8                      # error threshold for time integration

    # Nondimensionalization
    CharDim = GEO_units(length=1m, viscosity=1.0e19Pas, stress=100MPa)

    L       = Lx / 2.0                             # length of half the domain
    R       = 8.314 # 287.0                             # gas constant
    T       = 288.0                             # temperature
    M       = 0.028                             # molar mass of air

    dx_geo = dx*m
    dy_geo = dy*m
    L_geo = L*m
    R_geo = R*J/mol/K
    T_geo = T*K
    M_geo = M*kg
    ρs_geo = rho0*kg/m^3
    ρc_geo = drho*kg/m^3
    #βs_geo = beta*1/Pa
    ηa_geo = eta_air*Pa*s
    ηs_geo = eta_solid*Pa*s
    μs_geo = mu_solid*Pa
    g_geo  = g_y*m/s^2
    P_geo  = P0*Pa
    t_geo  = t*s

    dx_non = nondimensionalize(dx_geo, CharDim)
    dy_non = nondimensionalize(dy_geo, CharDim)
    L_non = nondimensionalize(L_geo, CharDim)
    R_non = nondimensionalize(R_geo, CharDim)
    T_non = nondimensionalize(T_geo, CharDim)
    M_non = nondimensionalize(M_geo, CharDim)
    ρs_non = nondimensionalize(ρs_geo, CharDim)
    ρc_non = nondimensionalize(ρc_geo, CharDim)
    #βa_non = nondimensionalize(βa_geo, CharDim)
    ηa_non = nondimensionalize(ηa_geo, CharDim)
    ηs_non = nondimensionalize(ηs_geo, CharDim)
    μs_non = nondimensionalize(μs_geo, CharDim)
    g_non  = nondimensionalize(g_geo, CharDim)
    P_non  = nondimensionalize(P_geo, CharDim)
    t_non  = nondimensionalize(t_geo, CharDim)

    # Reduce allocations
    μ       = zeros(Float64, nx, ny)
    μ_c     = zeros(Float64, nx + 1, ny + 1)
    η       = zeros(Float64, nx, ny)
    η_c     = zeros(Float64, nx + 1, ny + 1)
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
    rho     = ones(Float64, nx + 2, ny + 2) .* ρs_non
    P       = zeros(Float64, nx + 2, ny + 2)
    divV    = zeros(Float64, nx, ny)
    Exx     = zeros(Float64, nx, ny)
    Eyy     = zeros(Float64, nx, ny)
    Ezz     = zeros(Float64, nx, ny)
    Exy     = zeros(Float64, nx + 1, ny + 1)
    Sxx     = zeros(Float64, nx, ny)
    Syy     = zeros(Float64, nx, ny)
    Szz     = zeros(Float64, nx, ny)
    Sxx_old = zeros(Float64, nx, ny)
    Syy_old = zeros(Float64, nx, ny)
    Szz_old = zeros(Float64, nx, ny)
    Sxy     = zeros(Float64, nx + 1, ny + 1)
    Sxy_old = zeros(Float64, nx + 1, ny + 1)
    Mx      = zeros(Float64, nx + 1, ny + 2)
    FMxx    = zeros(Float64, nx, ny)
    FMxy    = zeros(Float64, nx - 1, ny + 1)
    dMxdt   = zeros(Float64, nx + 1, ny + 2)
    MxRes   = zeros(Float64, nx - 1, ny)
    Vx      = zeros(Float64, nx + 1, ny + 2)
    My      = zeros(Float64, nx + 2, ny + 1)
    FMyy    = zeros(Float64, nx, ny)
    FMyx    = zeros(Float64, nx + 1, ny - 1)
    dMydt   = zeros(Float64, nx + 2, ny + 1)
    MyRes   = zeros(Float64, nx, ny - 1)
    Vy      = zeros(Float64, nx + 2, ny + 1)
    rho_old = zeros(Float64, nx + 2, ny + 2)
    Mx_old  = zeros(Float64, nx + 1, ny + 2)
    My_old  = zeros(Float64, nx + 2, ny + 1)
    dMxdtau = zeros(Float64, nx - 1, ny)
    dMydtau = zeros(Float64, nx, ny - 1)
    radrho  = zeros(Float64, nx + 2, ny + 2)
    radVx   = zeros(Float64, nx + 1, ny + 2)
    radVy   = zeros(Float64, nx + 2, ny + 1)
    radV    = zeros(Float64, nx + 1, ny + 1)
    radμc   = zeros(Float64, nx + 1, ny + 1)
    beta_vec= zeros(Float64, nx + 2, ny + 2)
    maskrho_air   = zeros(Float64, nx, ny)
    maskVx_air    = zeros(Float64, nx, ny)
    maskVy_air    = zeros(Float64, nx, ny)
    maskV_solid   = zeros(Float64, nx, ny)
    maskμc_solid  = zeros(Float64, nx + 1, ny + 1)

    maskrho_solid     = zeros(Float64, nx, ny)
    maskVx_solid      = zeros(Float64, nx, ny)
    maskVy_solid      = zeros(Float64, nx, ny)
    maskV_air         = zeros(Float64, nx, ny)
    maskμc_air        = zeros(Float64, nx + 1, ny + 1)

    # Initialization
    Xv      = range(-Lx/2.0, stop=Lx/2.0, step=dx)                      # x-coordinates of velocity nodes
    Xc      = range(-(Lx+dx)/2.0, stop=(Lx+dx)/2.0, step=dx)            # x-coordinates of density nodes
    Yv      = range(-Ly/2.0, stop=Ly/2.0, step=dy)                      # y-coordinates of velocity nodes
    Yc      = range(-(Ly+dy)/2.0, stop=(Ly+dy)/2.0, step=dy)            # y-coordinates of density nodes
    
    # Nondimensional coordinates
    Xv_non  = range(-L_non, stop=L_non, step=dx_non)                    # x-coordinates of velocity nodes
    Xc_non  = range(-(L_non+dx_non), stop=(L_non+dx_non), step=dx_non)  # x-coordinates of density nodes
    Yv_non  = range(-L_non, stop=L_non, step=dy_non)                    # y-coordinates of velocity nodes
    Yc_non  = range(-(L_non+dy_non), stop=(L_non+dy_non), step=dy_non)  # y-coordinates of density nodes

    x2dc, y2dc   = meshgrid(Xc,Yc)                                        # 2d mesh of x- and y-coordinates of density nodes
    x2dVx, y2dVx = meshgrid(Xv,Yc)                                      # 2d mesh of x- and y-coordinates of velocity nodes in x-direction
    x2dVy, y2dVy = meshgrid(Xc,Yv)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction
    x2dV,  y2dV  = meshgrid(Xv,Yv)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction
    x2dμc, y2dμc = meshgrid(Xv,Yv)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction

    x2dc_non, y2dc_non   = meshgrid(Xc_non,Yc_non)                                      # 2d mesh of x- and y-coordinates of density nodes
    x2dVx_non, y2dVx_non = meshgrid(Xv_non,Yc_non)                                      # 2d mesh of x- and y-coordinates of velocity nodes in x-direction
    x2dVy_non, y2dVy_non = meshgrid(Xc_non,Yv_non)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction
    x2dμc_non, y2dμc_non = meshgrid(Xv_non,Yv_non)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction
    
    locX      = 0.0                                                     # x-coordinate of circle
    locY      = -0.35*Ly                                                # y-coordinate of circle
    diam      = 0.2*Lx                                                  # diameter of circle
    radrho    .= sqrt.((x2dc .-locX).^2.0 .+ (y2dc .-locY).^2.0)        # radius of circle
    radVx     .= sqrt.((x2dVx.-locX).^2.0 .+ (y2dVx.-locY).^2.0)        # radius of circle
    radVy     .= sqrt.((x2dVy.-locX).^2.0 .+ (y2dVy.-locY).^2.0)        # radius of circle
    radV      .= sqrt.((x2dV .-locX).^2.0 .+ (y2dV .-locY).^2.0)
    radμc     .= sqrt.((x2dμc.-locX).^2.0 .+ (y2dμc.-locY).^2.0)        # radius of circle
    Xp        = [-1.0/2.0, -1.0/2.0, -0.3, -1.0/8.0, -0.01, -0.01, 0.01, 0.01, 1.0/8.0, 0.3, 1.0/2.0, 1.0/2.0].*Lx          # x-coordinates of polygon
    Yp        = [-1.0/2.0, -1.0/5.0, -1.0/5.0, -0.1, -0.15,  -0.28, -0.28, -0.15,  -0.1, -1.0/5.0, -1.0/5.0, -1.0/2.0].*Ly  # y-coordinates of polygon
    Xp2       = [-0.02, -0.02,  0.02,  0.02].*Lx                                                                            # x-coordinates of polygon
    Yp2       = [-0.15,  -0.35, -0.35, -0.15].*Ly                                                                             # y-coordinates of polygon

    inpolyrho   = inpolygon(x2dc,y2dc,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon
    inpolyVx    = inpolygon(x2dVx,y2dVx,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon
    inpolyVy    = inpolygon(x2dVy,y2dVy,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon
    inpolyV     = inpolygon(x2dV,y2dV,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon
    inpolyμc    = inpolygon(x2dμc,y2dμc,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon

    maskrho_air = 1.0 .- inpolyrho                                             # mask for density
    maskrho_air[radrho .< diam ./ 2.0] .= 1.0                               # mask for density in circle
    #maskrho_air[findall(<(diam ./ 2.0), radrho)] .= 1.0                               # mask for density in circle
    maskVx_air  = 1.0 .- inpolyVx                               # mask for velocity in x-direction
    maskVx_air[radVx .< diam ./ 2.0] .= 1.0                                 # mask for velocity in x-direction in circle
    maskVy_air  = 1.0 .- inpolyVy                               # mask for velocity in y-direction
    maskVy_air[radVy .< diam ./ 2.0] .= 1.0                                 # mask for velocity in y-direction in circle
    #maskV_air   = 1.0 .- inpolyV
    #maskV_air[radV .< diam ./ 2.0] .= 1.0
    maskμc_air = 1.0 .- inpolyμc                                             # mask for density
    maskμc_air[radμc .< diam ./ 2.0] .= 1.0                               # mask for density in circle

    maskrho_solid = 1.0 .- maskrho_air
    maskVx_solid  = 1.0 .- maskVx_air
    maskVy_solid  = 1.0 .- maskVy_air
    maskV_solid   = 1.0 .- maskV_air
    maskμc_solid  = 1.0 .- maskμc_air

    # Initial conditions
    P           .= reverse(cumsum(ρs_non.*g_non.*reverse(maskrho_solid, dims=2)*dy_non,dims=2), dims=2)
    P           .= P_non .+ exp.((.-0.0000005 .* (av_xy(x2dc_non)/ 0.1).^2.0) .+ (.-0.0000005 .* ((av_xy(y2dc_non') .- 5000.0) / 0.1)'.^2.0))          # initial pressure distribution

    for i in 1:ny
        Syy[:, (ny+1)-i] .= ρs_non .* g_non .* (y2dc_non[i] .+ maximum(abs.(y2dc_non)))
    end

    #@. η +=#= eta_air * (maskrho_air[2:end-1, 2:end-1] == 1.0) + =#eta_solid# * (maskrho_solid[2:end-1, 2:end-1] == 1.0)  # initial viscosity distribution 
    #@. η_c += #= eta_air * (maskμc_air[2:end-1, 2:end-1] == 1.0) + =#eta_solid# * (maskμc_solid[2:end-1, 2:end-1] == 1.0)  # initial viscosity distribution
    @. μ += #=eta_air * (maskrho_air[2:end-1, 2:end-1] == 1.0) + =#mu_solid# * (maskrho_solid[2:end-1, 2:end-1] == 1.0)  # initial shear modulus distribution
    @. μ_c += #= eta_air * (maskμc_air[2:end-1, 2:end-1] == 1.0) + =#mu_solid# * (maskμc_solid[2:end-1, 2:end-1] == 1.0)  # initial shear mdoulus distribution

    # Inital plot  
    P_plt = ustrip(dimensionalize(P, Pa, CharDim))
  
    fig = Figure()
    ax = Axis(fig[1, 1], xticks=([-5000, 0, 5000], ["-5000", "0", "5000"]), yticks=([-5000, 0, 5000], ["-5000", "0", "5000"]),
                yticklabelsize=25, xticklabelsize=25, xlabelsize=25, ylabelsize=25)#, title="0")
    p1 = heatmap!(ax, x2dc, y2dc, P_plt, shading=false, colormap=Reverse(:roma))
    #p1 = heatmap!(ax, x2dc, y2dc , data_plt, shading=false, colorrange=(0.0, 350.0), colormap=Reverse(:roma))
    #p1 = heatmap!(ax, x2dc, y2dc , rho, shading=false, colormap=Reverse(:roma))#, colorrange=(P0, P0*2))
    #Colorbar(fig[1, 2], p1, label="Velocity [m/s]", labelsize=25, ticklabelsize=25)
    Colorbar(fig[1, 2], p1, label="Pressure [Pa]", labelsize=25, ticklabelsize=25)
    points_XpYp = Point2f[]
    for i in eachindex(Xp)
        x = Xp[i]
        y = Yp[i]
        point = (x,y)
        push!(points_XpYp, point)
    end
    #lines!(ax, Xp, Yp, color = :white) # only conduit
    display(fig)
    #save("./Plots/Shockwave_1/0.png", fig)

    # Solver
    for it = 1:100                                                   # 60 iterations max for all equilibrated
        @show it
        rho_old .= rho                                              # save old density values for time integration
        Sxx_old .= Sxx                                              # save old stress values in x-direction for time integration
        Syy_old .= Syy
        Sxy_old .= Sxy
        Szz_old .= Szz
        Mx_old .= Mx                                                # save old momentum values in x-direction for time integration
        My_old .= My                                                # save old momentum values in y-direction for time integration
        err = 1.0                                                   # error for time integration
        iter = 0                                                    
        it_counter = 0
        dMxdtau .= 0.0                                              # stress derivative of momentum in x-direction
        dMydtau .= 0.0                                              # stress derivative of momentum in y-direction
        while err > err_threshold
            iter += 1        
            # Compressibility and speed of sound
            beta_vec               .= 1.0 ./ P
            c_loc                  .= 1.0./sqrt.(rho[2:end-1,2:end-1].* beta_vec[2:end-1,2:end-1])                                                    # local speed of sound

            # Time steps
            dt                      = minimum([dx_non./maximum(abs.(Vx)), dy_non./maximum(abs.(Vy)), min(dx_non, dy_non) .* sqrt(maximum(rho .* beta_vec)), min(dx_non, dy_non).^2.0./ηs_non]).*4.0 # time step size
            dtPT                   .= min.(min.(min.(dx_non ./ abs.(av_x(Vx[:, 2:end-1])), dy_non ./ abs.(av_y(Vy[2:end-1, :]))), min(dx_non, dy_non).^2.0 ./ ηs_non), dx_non ./ c_loc) # time step size for pressure and temperature
            dtrho                  .= 1.0 ./ (1.0 ./ dt .+ 1.0 ./ (min(dx_non, dy_non) ./ c_loc ./ 4.1))                                                         # time step size for density
            
            # Conservation of mass
            Frhox, Frhoy            = update_densityflux(Frhox, Frhoy, Vx, Vy, rho)
            drhodt                  = update_drhodt(drhodt, rho, rho_old, dt)
            rhoRes                  = update_densityResidual(rhoRes, drhodt, Frhox, Frhoy, dx_non, dy_non, maskrho_air)
            rho                     = update_density(rho, rhoRes, dtrho, CFL_P) 
            @infiltrate
            # Pressure, strain-rates and stresses
            P                       = update_pressure(P, P_non, ρs_non, rho)
            divV                    = update_divergenceV(divV, Vx, Vy, dx_non, dy_non)
            Exx, Eyy, Exy           = update_strainrates(Exx, Eyy, Exy, Vx, Vy, divV, dx_non, dy_non)
            Sxx, Syy, Sxy           = update_totalstresses(Sxx, Syy, Sxy, P, ηs_non, Exx, Eyy, Exy)
                    
            # Conservation of the x and y-component of momentum
            Mx, My                  = update_advection_momentum(Mx, My, Vx, Vy, rho)
            FMxx, FMyy, FMxy, FMyx  = update_momentumflux(FMxx, FMyy, FMxy, FMyx, Vx, Vy, Mx, My, Sxx, Syy, Sxy, dx_non, dy_non)
            dMxdt, dMydt            = update_dMdt(dMxdt, dMydt, Mx, My, Mx_old, My_old, dt)
            MxRes, MyRes            = update_momentumResidual(MxRes, MyRes, dMxdt, dMydt, FMxx, FMyy, FMxy, FMyx, Sxx, Syy, Sxy, dx_non, dy_non, g_non, rho, maskVx_air, maskVy_air)
            dMxdtau, dMydtau        = update_dMdτ(dMxdtau, dMydtau, MxRes, MyRes, ksi)
            Mx, My                  = update_inertia_momentum(Mx, My, dMxdtau, dMydtau, dtPT, CFL_V)
            Vx, Vy                  = update_velocities(Vx, Vy, Mx, My, rho)
            
            if mod(iter, 1) == 0
                it_counter += 1
                err = maximum(abs.([rhoRes[:]; MxRes[:]; MyRes[:]]))                                                                      # error for time integration concerning density and momentum
                print("PT_iter = $iter, err = $err, rhoRes = $(maximum(abs.(rhoRes_a[:]))), MxRes = $(maximum(abs.(MxRes[:]))), MyRes = $(maximum(abs.(MyRes[:])))\n")
                if isnan.(err) == 1 
                    break
                end
            end

            if err <= err_threshold
                @show it_counter
            end
        end
        t = t + dt

        # Updating plot
        if mod(it-1, 1) == 0
            fig1 = Figure()
            ax1 = Axis(fig1[1,1], xticks=([-5000, 0, 5000], ["-500", "0", "5000"]), yticks=([-5000, 0, 5000], ["-5000", "0", "5000"]),
                    yticklabelsize=25, xticklabelsize=25, xlabelsize=25, ylabelsize=25, title="time = $time")

            X = av_xy(x2dc)
            Y = av_xy(y2dc)
            U = ustrip(dimensionalize(av_y(Vx), m/s, CharDim))
            V = ustrip(dimensionalize(av_x(Vy), m/s, CharDim))
            data_plt = sqrt.(U.^2 .+ V.^2)

            #hm = heatmap!(ax1, x2dc, y2dc, P, shading=false, colorrange=(P0, P0*2))
            hm = heatmap!(ax1, X, Y, data_plt, shading=false, colorrange=(0.0, 350), colormap=Reverse(:roma))
            #Colorbar(fig1[1,2],  hm, label="Pressure [Pa]", labelsize=25, ticklabelsize=25)
            Colorbar(fig1[1,2],  hm, label="Velocity [m/s]", labelsize=25, ticklabelsize=25)
            #lines!(ax1, Xp, Yp, color = :white)    # only conduit
            #poly!(ax,points_XpYp)
            
            #stepsize = 10
            #arrows!(ax1, X[1:stepsize:end, 1], Y[1, 1:stepsize:end], U[1:stepsize:end, 1:stepsize:end], V[1:stepsize:end, 1:stepsize:end], arrowsize=7, color = :white)
            display(fig1)

            #save("./Plots/Shockwave_1/$(it).png", fig1)
        end
    end
    return
end
