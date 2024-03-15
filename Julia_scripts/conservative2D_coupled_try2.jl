#using Pkg
#Pkg.activate(".")
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

function advection2d(A, Vx, Vy, dx_non, dy_non, dt)
    dtadv = dt
    A[1:end-1,:] = A[1:end-1,:] .- (Vx[2:end-1,:] .< 0.0) .* Vx[2:end-1,:] .* diff(A, dims=1) ./ dx_non .* dtadv
    A[2:end,:] = A[2:end,:] .- (Vx[2:end-1,:] .> 0.0) .* Vx[2:end-1,:] .* diff(A, dims=1) ./ dx_non .* dtadv
    A[:,1:end-1] = A[:,1:end-1] .- (Vy[:,2:end-1] .< 0.0) .* Vy[:,2:end-1] .* diff(A, dims=2) ./ dy_non .* d
    tadv
    A[:,2:end] = A[:,2:end] .- (Vy[:,2:end-1] .> 0.0) .* Vy[:,2:end-1] .* diff(A, dims=2) ./ dy_non .* dtadv
    return A
end

function find_boundary_indices(maskrho_air)
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
    return x_circ_ind, y_circ_ind
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
    rhoRes .= rhoRes .* mask[2:end-1, 2:end-1]                                                                               # applying mask to residual of density
    return rhoRes
end

function update_density(rho, rhoRes, dtrho, CFL_P)
    rho[2:end-1, 2:end-1] .= rho[2:end-1, 2:end-1] .+ rhoRes .* dtrho .* CFL_P                   # updating density
    return rho
end

function update_pressure(P, P_non, ρa, ρs, ρ, β, mask1, mask2)
    P .= P_non ./ ρa .* ρ .* mask1                                                 # equation of state for pressure depending on density
    P .+= (P_non .+  (1.0 ./ β) .* log.(ρ ./ ρs)) .* mask2                                                 # equation of state for pressure depending on density
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
    MxRes .= MxRes .* mask1[2:end-1, 2:end-1]                                                                               # applying mask to residual of momentum in x-direction
    MyRes .= MyRes .* mask2[2:end-1, 2:end-1]                                                                               # applying mask to residual of momentum in y-direction
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

function seismic_solver(P, beta_vec, rho, divV, Exx, Eyy, Exy, Vx, Vy, μ, μ_c, η, η_c, dt, maskrho, maskVx, maskVy)
    # Allocations
    dPdt    = zeros(Float64, nx + 2, ny + 2)
    τxx     = zeros(Float64, nx, ny)
    τyy     = zeros(Float64, nx, ny)
    τxy     = zeros(Float64, nx + 1, ny + 1)
    dτxxdt  = zeros(Float64, nx, ny)
    dτyydt  = zeros(Float64, nx, ny)
    dτxydt  = zeros(Float64, nx + 1, ny + 1)
    Pτxx    = zeros(Float64, nx - 1, ny)
    Pτyy    = zeros(Float64, nx, ny - 1)
    dVxdt   = zeros(Float64, nx + 1, ny + 2)
    dVydt   = zeros(Float64, nx + 2, ny + 1)

    dPdt[2:end-1, 2:end-1].= .-(1.0 ./ beta_vec[2:end-1, 2:end-1]) .* divV
    P .+= P .+ dPdt .* dt .* maskrho
    dτxxdt .= 2.0 .* μ .* Exx .- (μ ./ η) .* τxx
    dτyydt .= 2.0 .* μ .* Eyy .- (μ ./ η) .* τyy         
    dτxydt .= 2.0 .* μ_c .* Exy .- (μ_c ./ η_c) .* τxy             
    τxx .= τxx .+ dτxxdt .* dt
    τyy .= τyy .+ dτyydt .* dt
    τxy .= τxy .+ dτxydt .* dt
    Pτxx .=.-diff(P_s[2:end-1,2:end-1], dims=1) ./ dx_non .+ diff(τxx, dims=1) ./ dx_non .+ diff(τxy[2:end-1, :], dims=2) ./ dy_non
    Pτyy .=.-diff(P_s[2:end-1,2:end-1], dims=2) ./ dy_non .+ diff(τyy, dims=2) ./ dy_non .+ diff(τxy[:, 2:end-1], dims=1) ./ dx_non .+ g_non .* (rho[2:end-1, 2:end-2] .* (Vy[2:end-1, 2:end-1] .>= 0.0) .+ rho[2:end-1, 3:end-1] .* (Vy[2:end-1, 2:end-1] .< 0.0))
    dVxdt[2:end-1, 2:end-1] .= (1.0 ./ (rho[2:end-2, 2:end-1] .* (Vx[2:end-1, 2:end-1] .>= 0.0) .+ rho[3:end-1, 2:end-1] .* (Vx[2:end-1, 2:end-1] .< 0.0))) .* (Pτxx) .- Vx[2:end-1, 2:end-1] .* diff(av_x(Vx[:, 2:end-1]), dims=1) ./ dx_non .+ av_x(av_y(Vy[2:end-1, :])) .* diff(av_x(Vx[:, 2:end-1]), dims=1) ./ dy_non
    dVydt[2:end-1, 2:end-1] .= (1.0 ./ (rho[2:end-1, 2:end-2] .* (Vy[2:end-1, 2:end-1] .>= 0.0) .+ rho[2:end-1, 3:end-1] .* (Vy[2:end-1, 2:end-1] .< 0.0))) .* (Pτyy) .- Vy[2:end-1, 2:end-1] .* diff(av_y(Vy[2:end-1, :]), dims=2) ./ dy_non .+ av_y(av_x(Vx[:, 2:end-1])) .* diff(av_y(Vy[2:end-1, :]), dims=2) ./ dx_non
    Vx[2:end-1, :] .= Vx[2:end-1, :] .+ dVxdt[2:end-1, :] .* dt .* maskVx[2:end-1, :] 
    Vy[:, 2:end-1] .= Vy[:, 2:end-1] .+ dVydt[: ,2:end-1] .* dt .* maskVy[: ,2:end-1]
    return P, Vx, Vy
end

function set_pressure_bc!(P, surf_ind_x, surf_ind_y)
    P[surf_ind_x, surf_ind_y] .= P[surf_ind_x, surf_ind_y .+ 1]
    return P
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
    eta      = 1.81e-5                          # dynamic viscosity
    eta_air      = 1.0e-5                          # dynamic viscosity for air
    eta_solid    = 1.0e19                         # dynamic viscosity for solid
    mu_air      = 2.0e-3                           # shear modulus for air
    mu_solid    = 1.0e11                           # shear modulus for solid
    g_y       = 9.81 #1.0e-3                             # gravitational acceleration

    ν  = 0.4                            # Poisson's ratio
    λ  = (2.0 * mu_solid * ν) / (1.0 - 2.0 * ν)     # Lamé parameter 2
    K_s = (2.0 * mu_solid * (1.0 + ν)) / (3.0 * (1.0 - 2.0 * ν))                # model poisson ratio of the solid
    beta_air  = 1.0/141.0e3                         # compressibility
    beta_solid  = 1.0 / K_s                        # compressibility

    # Numerics
    nx      = 301                               # number of nodes in x-direction
    ny      = 301                               # number of nodes in y-direction
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
    CharDim = GEO_units(length=1m, viscosity=1.0e19Pas, stress=100MPa)

    L       = 500.0                             # length of half the domain
    R       = 8.314                             # gas constant
    T       = 288.0                             # temperature
    M       = 0.028                             # molar mass of air

    t_geo = t*s
    t0_geo = t0*s
    td_geo = t_d*s
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

    t_non = nondimensionalize(t_geo, CharDim)
    t0_non = nondimensionalize(t0_geo, CharDim)
    td_non = nondimensionalize(td_geo, CharDim)
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
    beta_vec= zeros(Float64, nx + 2, ny + 2)
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
    rhoRes_a= zeros(Float64, nx, ny)
    rhoRes_s= zeros(Float64, nx, ny)
    rho     = zeros(Float64, nx + 2, ny + 2)
    rho_a   = zeros(Float64, nx + 2, ny + 2)
    rho_s   = zeros(Float64, nx + 2, ny + 2)
    P       = zeros(Float64, nx + 2, ny + 2)
    P_a     = zeros(Float64, nx + 2, ny + 2)
    P_s     = zeros(Float64, nx + 2, ny + 2)
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
    Vx_plt  = zeros(Float64, nx + 1, ny + 2)
    Vx_s    = zeros(Float64, nx + 1, ny + 2)
    Vx      = zeros(Float64, nx + 1, ny + 2)
    My      = zeros(Float64, nx + 2, ny + 1)
    FMyy    = zeros(Float64, nx, ny)
    FMyx    = zeros(Float64, nx + 1, ny - 1)
    dMydt   = zeros(Float64, nx + 2, ny + 1)
    MyRes   = zeros(Float64, nx, ny - 1)
    Vy_plt  = zeros(Float64, nx + 2, ny + 1)
    Vy_s    = zeros(Float64, nx + 2, ny + 1)
    Vy      = zeros(Float64, nx + 2, ny + 1)
    rho_old = zeros(Float64, nx + 2, ny + 2)
    Mx_old  = zeros(Float64, nx + 1, ny + 2)
    My_old  = zeros(Float64, nx + 2, ny + 1)
    dMxdtau = zeros(Float64, nx - 1, ny)
    dMydtau = zeros(Float64, nx, ny - 1)
    radrho  = zeros(Float64, nx + 2, ny + 2)
    radVx   = zeros(Float64, nx + 1, ny + 2)
    radVy   = zeros(Float64, nx + 2, ny + 1)
    radμc   = zeros(Float64, nx + 1, ny + 1)

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
    Xv_non  = range(-L_non, stop=L_non, step=dx_non)                    # x-coordinates of velocity nodes
    Xc_non  = range(-(L_non+dx_non), stop=(L_non+dx_non), step=dx_non)  # x-coordinates of density nodes
    Yv_non  = range(-L_non, stop=L_non, step=dy_non)                    # y-coordinates of velocity nodes
    Yc_non  = range(-(L_non+dy_non), stop=(L_non+dy_non), step=dy_non)  # y-coordinates of density nodes

    x2dc, y2dc   = meshgrid(Xc,Yc)                                        # 2d mesh of x- and y-coordinates of density nodes
    x2dVx, y2dVx = meshgrid(Xv,Yc)                                      # 2d mesh of x- and y-coordinates of velocity nodes in x-direction
    x2dVy, y2dVy = meshgrid(Xc,Yv)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction
    x2dμc, y2dμc = meshgrid(Xv,Yv)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction

    x2dc_non, y2dc_non   = meshgrid(Xc_non,Yc_non)                                      # 2d mesh of x- and y-coordinates of density nodes
    x2dVx_non, y2dVx_non = meshgrid(Xv_non,Yc_non)                                      # 2d mesh of x- and y-coordinates of velocity nodes in x-direction
    x2dVy_non, y2dVy_non = meshgrid(Xc_non,Yv_non)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction
    x2dμc_non, y2dμc_non = meshgrid(Xv_non,Yv_non)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction

    locX       = 0.0                                                     # x-coordinate of circle
    locY       = -0.35*Ly                                                # y-coordinate of circle
    diam       = 0.2*Lx                                                  # diameter of circle
    radrho    .= sqrt.((x2dc .-locX).^2.0 .+ (y2dc .-locY).^2.0)        # radius of circle
    radVx     .= sqrt.((x2dVx.-locX).^2.0 .+ (y2dVx.-locY).^2.0)        # radius of circle
    radVy     .= sqrt.((x2dVy.-locX).^2.0 .+ (y2dVy.-locY).^2.0)        # radius of circle
    radμc     .= sqrt.((x2dμc.-locX).^2.0 .+ (y2dμc.-locY).^2.0)        # radius of circle
    Xp         = [-1.0/2.0, -1.0/2.0, -0.3, -1.0/8.0, -0.01, -0.01, 0.01, 0.01, 1.0/8.0, 0.3, 1.0/2.0, 1.0/2.0].*Lx          # x-coordinates of polygon
    Yp         = [-1.0/2.0, -1.0/5.0, -1.0/5.0, -0.1, -0.15,  -0.28, -0.28, -0.15,  -0.1, -1.0/5.0, -1.0/5.0, -1.0/2.0].*Ly  # y-coordinates of polygon
    Xp2        = [-0.02, -0.02,  0.02,  0.02].*Lx                                                                            # x-coordinates of polygon
    Yp2        = [-0.14,  -0.25, -0.25, -0.14].*Ly                                                                             # y-coordinates of polygon

    inpolyrho   = inpolygon(x2dc,y2dc,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon
    inpolyVx    = inpolygon(x2dVx,y2dVx,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon
    inpolyVy    = inpolygon(x2dVy,y2dVy,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon
    inpolyμc    = inpolygon(x2dμc,y2dμc,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon

    inpolyrho_cond = inpolygon(x2dc,y2dc,Xp2,Yp2)
    inpolyVx_cond  = inpolygon(x2dVx,y2dVx,Xp2,Yp2)
    inpolyVy_cond  = inpolygon(x2dVy,y2dVy,Xp2,Yp2)
    inpolyμc_cond  = inpolygon(x2dμc,y2dμc,Xp2,Yp2)

    maskrho_air = 1.0 .- inpolyrho                                             # mask for density
    maskrho_air[radrho .< diam ./ 2.0] .= 1.0                               # mask for density in circle
    maskμc_air = 1.0 .- inpolyμc                                             # mask for density
    maskμc_air[radμc .< diam ./ 2.0] .= 1.0                               # mask for density in circle
    maskVx_air  = 1.0 .- inpolyVx                               # mask for velocity in x-direction
    maskVx_air[radVx .< diam ./ 2.0] .= 1.0                                 # mask for velocity in x-direction in circle
    maskVy_air  = 1.0 .- inpolyVy                               # mask for velocity in y-direction
    maskVy_air[radVy .< diam ./ 2.0] .= 1.0                                 # mask for velocity in y-direction in circle

    maskrho_solid = 1.0 .- maskrho_air
    maskμc_solid  = 1.0 .- maskμc_air
    maskVx_solid  = 1.0 .- maskVx_air
    maskVy_solid  = 1.0 .- maskVy_air

    x_circ_ind, y_circ_ind = find_boundary_indices(maskrho_air)

    # 90 is the indicie at the surface with 301 nodes
    ind_y = findall(x->x>90, y_circ_ind)
    surf_ind_y = y_circ_ind[ind_y]
    surf_ind_x = x_circ_ind[ind_y]

    # Initial conditions
    P         .= (P_non.*exp.(-g_non.*(y2dc_non[1:end-1, 1:end-1].+1.0./5.0.*L_non).*M_non ./ T_non ./ R_non)) .* maskrho_air .+ reverse(cumsum(ρs_non.*g_non.*reverse(maskrho_solid, dims=2).*dy_non,dims=2), dims=2) .* maskrho_solid                         # barometric setting atmosphere: P = P0*exp(-(g*(h-h0)*M)/(T*R)) M: Mass density of air, T: Temperature, R: Gas constant, h: height, P0: Pressure at sea level

    #rho       .= (rho0 .* P) ./ P0 .* maskrho_air .+ rho_solid .* maskrho_solid                         # equation of state for density depending on pressure
    rho       .= ρa_non ./ P_non .* P .* (maskrho_air .== 1.0) .+ ρs_non .* (maskrho_solid .== 1.0)    # equation of state for density depending on pressure
    rho[radrho .< diam ./ 2.0] .= ρa_non .+ ρc_non                                                          # initial density in the circle
    rho[inpolygon(x2dc,y2dc,Xp2,Yp2) .== 1.0 .&& y2dc.< locY.+diam./2.0] .= ρa_non .+ ρc_non          # initial density of the conduit overlapping with the circular chamber (so in the chamber)
    rho[inpolygon(x2dc,y2dc,Xp2,Yp2) .== 1.0 .&& y2dc.>=locY.+diam./2.0] .= .-0.0 .* (y2dc[inpolygon(x2dc,y2dc,Xp2,Yp2).==1.0 .&& y2dc.>=locY.+diam./2.0].+0.15.*Ly).*ρc_non./(-0.15.-(locY.+diam./2.0)).+ρa_non    # initial density in the conduit
    #P .= reverse(cumsum(ρs_non.*g_non.*reverse(maskrho_solid, dims=2)*dy_non,dims=2), dims=2)
    P        .+= (P_non .* rho) ./ ρa_non .* maskrho_air#.+ (rho[radrho .< diam ./ 2.0] .* g .* depth[radrho .< diam ./ 2.0])#P0 ./ rho0 .* rho
    #P         .+=(P0 .+  (1.0 ./ beta) .* log.(rho ./ rho_solid)) .* maskrho_solid


    # Initial parameter matrices for both phases
    @. η += ηa_non* (maskrho_air[2:end-1, 2:end-1] == 1.0) + ηs_non * (maskrho_solid[2:end-1, 2:end-1] == 1.0)     # initial viscosity distribution
    @. μ += μa_non* (maskrho_air[2:end-1, 2:end-1] == 1.0) + μs_non * (maskrho_solid[2:end-1, 2:end-1] == 1.0)       # initial viscosity distribution
    
    @. η_c += ηa_non* (maskμc_air == 1.0) + ηs_non * (maskμc_solid == 1.0)   # initial viscosity distribution for corner nodes
    @. μ_c += μa_non* (maskμc_air == 1.0) + μs_non * (maskμc_solid == 1.0)     # initial viscosity distribution for corner nodes

    # Inital plot 
    P_dim  = ustrip(dimensionalize(P, Pa, CharDim))
    xc_dim = ustrip(dimensionalize(Xc_non, m, CharDim))
    yc_dim = ustrip(dimensionalize(Yc_non, m, CharDim))
    
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
    
    p1 = heatmap!(ax, xc_dim, yc_dim, P_dim, colormap=Reverse(:roma))
    #p1 = heatmap!(ax, x2dc, y2dc , data_plt, shading=false, colorrange=(0.0, 350.0), colormap=Reverse(:roma))
    #p1 = heatmap!(ax, x2dc, y2dc , rho, shading=false, colormap=Reverse(:roma))#, colorrange=(P0, P0*2))
    #Colorbar(fig[1, 2], p1, label="Velocity", labelsize=25, ticklabelsize=25)
    Colorbar(fig[1, 2], p1, label="Pressure ", labelsize=25, ticklabelsize=25)
    #scatter!(ax, x_circ, y_circ, color = :white, markersize=4.0) # conduit and chamber
    display(fig)
    #if save_plt
    #    mkdir(path)
    save("initial_conditions.png", fig)
    #end

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
        # PT Solver for Atmosphere part
        while err > 1.0e-3 #&& iter < 25*360
            iter += 1
            beta_vec .= 1.0 ./ P
            @infiltrate
            #c_loc .= 1.0./sqrt.(rho[2:end-1,2:end-1].* beta_vec[2:end-1,2:end-1]) .* maskrho_air[2:end-1,2:end-1]                                                    # local speed of sound
            dt = minimum([dx_non./maximum(abs.(Vx)), dy_non./maximum(abs.(Vy)), min(dx_non, dy_non) .* sqrt(maximum(rho .* beta_vec)), min(dx_non, dy_non).^2.0./maximum(η)]).*4.1 # time step size  
            dtPT .= min.(min.(min.(dx_non ./ abs.(av_x(Vx[:, 2:end-1])), dx_non ./ abs.(av_y(Vy[2:end-1, :]))), min(dx_non, dy_non).^2.0 ./ η))#, dx_non ./ c_loc) # time step size for pressure and temperature
            dtrho .= 1.0 ./ (1.0 ./ dt .+ 1.0 ./ (min(dx_non, dy_non) ./ c_loc ./ 4.1))                                                         # time step size for density

            # Conservation of mass
            Frhox, Frhoy            = update_densityflux(Frhox, Frhoy, Vx, Vy, rho)
            drhodt                  = update_drhodt(drhodt, rho, rho_old, dt)
            rhoRes                  = update_densityResidual(rhoRes, drhodt, Frhox, Frhoy, dx_non, dy_non, maskrho_air)
            rho                     = update_density(rho, rhoRes, dtrho, CFL_P) #Frhox .= (Vx .> 0.0).*Vx.*rho[1:end-1, :] .+ (Vx .< 0.0).*Vx.*rho[2:end, :] # mass flux in x-direction

            # Pressure, strain-rates and stresses
            P                       = update_pressure(P, P_non, ρa_non, ρs_non, rho, beta_vec, maskrho_air, maskrho_solid)
            divV                    = update_divergenceV(divV, Vx, Vy, dx_non, dy_non)
            Exx, Eyy, Exy           = update_strainrates(Exx, Eyy, Exy, Vx, Vy, divV, dx_non, dy_non)
            Sxx, Syy, Sxy           = update_totalstresses(Sxx, Syy, Sxy, P, ηa_non, Exx, Eyy, Exy)
                    
            # Conservation of the x and y-component of momentum
            Mx, My                  = update_advection_momentum(Mx, My, Vx, Vy, rho)
            FMxx, FMyy, FMxy, FMyx  = update_momentumflux(FMxx, FMyy, FMxy, FMyx, Vx, Vy, Mx, My, Sxx, Syy, Sxy, dx_non, dy_non)
            dMxdt, dMydt            = update_dMdt(dMxdt, dMydt, Mx, My, Mx_old, My_old, dt)
            MxRes, MyRes            = update_momentumResidual(MxRes, MyRes, dMxdt, dMydt, FMxx, FMyy, FMxy, FMyx, Sxx, Syy, Sxy, dx_non, dy_non, g_non, rho, maskVx_air, maskVy_air)
            dMxdtau, dMydtau        = update_dMdτ(dMxdtau, dMydtau, MxRes, MyRes, ksi)
            Mx, My                  = update_inertia_momentum(Mx, My, dMxdtau, dMydtau, dtPT, CFL_V)
            Vx, Vy                  = update_velocities(Vx, Vy, Mx, My, rho)

            #if it >= 1 && iter ==2290
            #    @infiltrate
            #end

            if mod(iter, 25) == 0
                it_counter += 1
                err = maximum(abs.([rhoRes_a[:]; MxRes[:]; MyRes[:]]))
                print("PT_iter = $iter, err = $err, rhoRes = $(maximum(abs.(rhoRes_a[:]))), MxRes = $(maximum(abs.(MxRes[:]))), MyRes = $(maximum(abs.(MyRes[:])))\n")
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
                #@infiltrate
            end
        end

        # Explicit Solver for Earth part                                                                           # time derivative of density
        #rhoRes_s .= .-drhodt[2:end-1, 2:end-1] .- diff(Frhox[:, 2:end-1], dims=1)./dx_non .- diff(Frhoy[2:end-1, :], dims=2)./dy_non        # updating residual of density
        #rhoRes_s .= rhoRes_s .* maskrho_solid
        #rho_s[2:end-1, 2:end-1] .= rho_s[2:end-1, 2:end-1] .+ rhoRes_s .* dtrho .* CFL_P
        
        #P, Vx, Vy = seismic_solver(P, beta_vec, rho, divV, Exx, Eyy, Exy, Vx, Vy, μ, μ_c, η, η_c, dt, maskrho_solid, maskVx_solid, maskVy_solid)
        
        #=dPdt[2:end-1, 2:end-1].= .-(1.0 ./ beta_vec[2:end-1, 2:end-1]) .* divV
        P .+= P .+ dPdt .* dt .* maskrho_solid
        dτxxdt .= 2.0 .* μ .* Exx .- (μ ./ η) .* τxx
        dτyydt .= 2.0 .* μ .* Eyy .- (μ ./ η) .* τyy         
        dτxydt .= 2.0 .* μ_c .* Exy .- (μ_c ./ η_c) .* τxy             
        τxx .= τxx .+ dτxxdt .* dt
        τyy .= τyy .+ dτyydt .* dt
        τxy .= τxy .+ dτxydt .* dt
        Pτxx .=.-diff(P_s[2:end-1,2:end-1], dims=1) ./ dx_non .+ diff(τxx, dims=1) ./ dx_non .+ diff(τxy[2:end-1, :], dims=2) ./ dy_non
        Pτyy .=.-diff(P_s[2:end-1,2:end-1], dims=2) ./ dy_non .+ diff(τyy, dims=2) ./ dy_non .+ diff(τxy[:, 2:end-1], dims=1) ./ dx_non .+ g_non .* (rho[2:end-1, 2:end-2] .* (Vy_s[2:end-1, 2:end-1] .>= 0.0) .+ rho[2:end-1, 3:end-1] .* (Vy_s[2:end-1, 2:end-1] .< 0.0))
        dVxdt[2:end-1, 2:end-1] .= (1.0 ./ (rho[2:end-2, 2:end-1] .* (Vx_s[2:end-1, 2:end-1] .>= 0.0) .+ rho[3:end-1, 2:end-1] .* (Vx_s[2:end-1, 2:end-1] .< 0.0))) .* (Pτxx) .- Vx_s[2:end-1, 2:end-1] .* diff(av_x(Vx_s[:, 2:end-1]), dims=1) ./ dx_non .+ av_x(av_y(Vy_s[2:end-1, :])) .* diff(av_x(Vx_s[:, 2:end-1]), dims=1) ./ dy_non
        dVydt[2:end-1, 2:end-1] .= (1.0 ./ (rho[2:end-1, 2:end-2] .* (Vy_s[2:end-1, 2:end-1] .>= 0.0) .+ rho[2:end-1, 3:end-1] .* (Vy_s[2:end-1, 2:end-1] .< 0.0))) .* (Pτyy) .- Vy_s[2:end-1, 2:end-1] .* diff(av_y(Vy_s[2:end-1, :]), dims=2) ./ dy_non .+ av_y(av_x(Vx_s[:, 2:end-1])) .* diff(av_y(Vy_s[2:end-1, :]), dims=2) ./ dx_non
        Vx_s[2:end-1, :] .= Vx_s[2:end-1, :] .+ dVxdt[2:end-1, :] .* dt .* maskVx_solid[2:end-1, :] 
        Vy_s[:, 2:end-1] .= Vy_s[:, 2:end-1] .+ dVydt[: ,2:end-1] .* dt .* maskVy_solid[: ,2:end-1]=#
        
        t_non = t_non + dt

        # Bring together both parts of velocity and pressure calculations
        #rho .= rho_a .+ rho_s

        #P .= P_a .+ P_s

        #Vx_plt .= Vx .+ Vx_s
        #Vy_plt .= Vy .+ Vy_s

        # Updating plot
        if mod(it-1, 1) == 0
            Vy_dim = dimensionalize(Vy, m/s, CharDim)
            Vx_dim = dimensionalize(Vx, m/s, CharDim)
            Vx_av = av_x(ustrip(Vx_dim))
            Vy_av = av_y(ustrip(Vy_dim))
            #@infiltrate
            data_plt = sqrt.(Vx_av[:, 2:end-1].^2.0 .+ Vy_av[2:end-1, :].^2.0)
            P_plt = ustrip(dimensionalize(P, Pa, CharDim))
            t_dim = dimensionalize(t_non, s, CharDim)

            fig1 = Figure(size=(800,600))
            ax1 = Axis(fig1[1,1], xticks=([-500, 0, 500], ["-500", "0", "500"]), yticks=([-500, 0, 500], ["-500", "0", "500"]),
                    yticklabelsize=25, xticklabelsize=25, xlabelsize=25, ylabelsize=25, title="time = $t_dim")
            ax2 = Axis(fig1[2,1], xticks=([-500, 0, 500], ["-500", "0", "500"]), yticks=([-500, 0, 500], ["-500", "0", "500"]),
                    yticklabelsize=25, xticklabelsize=25, xlabelsize=25, ylabelsize=25, title="time = $t_dim")

            X = av_xy(x2dc)
            Y = av_xy(y2dc)
            U = av_y(Vx)
            V = av_x(Vy)

            #data_plt = sqrt.(U.^2.0 .+ V.^2.0)

            hm = heatmap!(ax1, xc_dim, yc_dim, P_plt)# colorrange=(P0, P0*2))
            #hm = heatmap!(ax1, X[2:end-1, 2:end-1], Y[2:end-1, 2:end-1], data_plt[2:end-1, 2:end-1], shading=false, colormap=Reverse(:roma))#, colorrange=(0.0, 0.1),)
            hm2 = heatmap!(ax2, xc_dim[2:end-1], yc_dim[2:end-1], data_plt[2:end-1, 2:end-1], colormap=Reverse(:roma), colorrange=(0.0, 1000.0),)
            Colorbar(fig1[1,2],  hm, label="Pressure [Pa]", labelsize=25, ticklabelsize=25)
            #Colorbar(fig1[1,2],  hm, label="Velocity", labelsize=25, ticklabelsize=25)
            Colorbar(fig1[2,2],  hm2, label="Velocity [m/s]", labelsize=25, ticklabelsize=25)
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
end
