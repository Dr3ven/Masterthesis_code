using CairoMakie
using TimerOutputs

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

function conservative2D_coupled()
    # Physics
    Lx      = 1.0e3                             # length of domain in x-direction
    Ly      = Lx                                # length of domain in y-direction
    rho_air    = 1.225                             # density at rest atmosphere
    rho_solid    = 1.225 #2.7e3                             # density solid
    drho    = 0.0 # 3.0e2                             # density difference between chamber and conduit
    Vx0     = 0.0                               # starting velocity in x-direction
    P0      = 1.0e5                             # pressure at rest
    P0_solid        = 1.0e5# 4.02305e6                             # pressure at rest
    beta_air        = 1.0/141.0e3                       # compressibility
    beta_solid      = 1.0/141.0e3                       # compressibility
    eta      = 1.81e-5                           # dynamic viscosity
    eta_air      = 1.81e-5                           # dynamic viscosity
    eta_solid    = 1.81e-5                           # dynamic viscosity
    mu      = 1.0e-5                            # shear modulus
    mu_air          = 1.0e20                            # shear modulus
    mu_solid        = 1.0e20                            # shear modulus
    g       = 9.81                              # gravitational acceleration

    # Numerics
    nx      = 301                               # number of nodes in x-direction
    ny      = 301                               # number of nodes in y-direction
    dx      = Lx/(nx)                           # grid spacing in x-direction
    dy      = Ly/(ny)                           # grid spacing in y-direction
    #nout    = 1
    dt      = 1.0                           # time step size
    CFL_P   = 1.0/16.0                           # Courant number for pressure
    CFL_V   = 1.0/16.0                          # Courant number for velocity
    psi     = 3.0#5.25                               # dampening factor for pseudo transient iteration
    ksi     = 1.0 - (psi / nx) # 0.0*0.95       # relaxation factor for stress when nx=ny

    # Reduce allocations
    β       = zeros(Float64, nx + 2, ny + 2)
    beta_vec= zeros(Float64, nx + 2, ny + 2)
    μ       = zeros(Float64, nx, ny)
    η       = zeros(Float64, nx, ny)
    μ_c       = zeros(Float64, nx + 1, ny + 1)
    η_c       = zeros(Float64, nx + 1, ny + 1)
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
    rho     = ones(Float64, nx + 2, ny + 2)
    P       = zeros(Float64, nx + 2, ny + 2)
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
    radμc   = zeros(Float64, nx + 1, ny + 1)
    maskrho_air   = zeros(Float64, nx, ny)
    maskμc_air    = zeros(Float64, nx + 1, ny + 1)
    maskVx_air    = zeros(Float64, nx, ny)
    maskVy_air    = zeros(Float64, nx, ny)

    maskrho_solid     = zeros(Float64, nx, ny)
    maskμc_air        = zeros(Float64, nx + 1, ny + 1)
    maskVx_solid      = zeros(Float64, nx, ny)
    maskVy_solid      = zeros(Float64, nx, ny)


    # Initialization
    Xv      = range(-Lx/2.0, stop=Lx/2.0, step=dx)                      # x-coordinates of velocity nodes
    Xc      = range(-(Lx+dx)/2.0, stop=(Lx+dx)/2.0, step=dx)            # x-coordinates of density nodes
    Yv      = range(-Ly/2.0, stop=Ly/2.0, step=dy)                      # y-coordinates of velocity nodes
    #Yv      = range(-300.0, stop=700.0, step=dy)
    Yc      = range(-(Ly+dy)/2.0, stop=(Ly+dy)/2.0, step=dy)            # y-coordinates of density nodes
    #Yc      = range(-300.0-(dy/2.0), stop=700.0+(dy/2.0), step=dy)            # y-coordinates of density nodes
    x2dc, y2dc = meshgrid(Xc,Yc)                                        # 2d mesh of x- and y-coordinates of density nodes
    x2dVx, y2dVx = meshgrid(Xv,Yc)                                      # 2d mesh of x- and y-coordinates of velocity nodes in x-direction
    x2dVy, y2dVy = meshgrid(Xc,Yv)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction
    x2dμc, y2dμc = meshgrid(Xv,Yv)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction

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
    Xp2       = [-0.02, -0.02,  0.02,  0.02].*Lx                                                                            # x-coordinates of polygon
    Yp2       = [-0.15,  -0.35, -0.35, -0.15].*Ly                                                                             # y-coordinates of polygon

    inpolyrho   = inpolygon(x2dc,y2dc,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon
    inpolyμc    = inpolygon(x2dμc,y2dμc,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon
    inpolyVx    = inpolygon(x2dVx,y2dVx,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon
    inpolyVy    = inpolygon(x2dVy,y2dVy,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon

    maskrho_air = 1.0 .- inpolyrho                                             # mask for density
    maskrho_air[radrho .< diam ./ 2.0] .= 1.0                               # mask for density in circle
    maskμc_air = 1.0 .- inpolyμc                                             # mask for density
    maskμc_air[radμc .< diam ./ 2.0] .= 1.0                               # mask for density in circle
    maskVx_air  = 1.0 .- inpolyVx                               # mask for velocity in x-direction
    maskVx_air[radVx .< diam ./ 2.0] .= 1.0                                 # mask for velocity in x-direction in circle
    maskVy_air  = 1.0 .- inpolyVy                               # mask for velocity in y-direction
    maskVy_air[radVy .< diam ./ 2.0] .= 1.0                                 # mask for velocity in y-direction in circle

    maskrho_solid = 1.0 .- maskrho_air
    maskμc_solid = 1.0 .- maskμc_air
    maskVx_solid  = 1.0 .- maskVx_air
    maskVy_solid  = 1.0 .- maskVy_air

    x_circ_ind = []
    y_circ_ind = []

    for y= 3:size(maskrho_air)[1] - 2
        for x = 3:size(maskrho_air)[2] - 2
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

    @. β += beta_air * (maskrho_air == 1.0) + beta_solid * (maskrho_solid == 1.0)  # initial viscosity distribution    # So sollte der Conduitverlauf sein: lines!(ax, x, .-(100 .+ x .+ (x.^2.0 ./ 2.0) .+ (x.^3.0 ./ 3.0)))

    P         .= P0.*exp.(-g.*(y2dc.+1.0./5.0.*Ly).*0.028./288.0./8.314) #.*maskrho_air# barometric setting atmosphere: P = P0*exp(-(g*(h-h0)*M)/(T*R)) M: Mass density of air, T: Temperature, R: Gas constant, h: height, P0: Pressure at sea level
    #rho       .= rho_solid .* maskrho_solid
    rho       .+= (rho_air .* P) ./ P0 #.* maskrho_air #rho_air ./ P0 .* P             # equation of state for density depending on pressure
    rho[radrho .< diam ./ 2.0] .= rho_air .+ drho                                       # initial density in the circle/chamber

    #rho[inpolygon(x2dc,y2dc,Xp2,Yp2) .== 1.0 .&& y2dc.< locY.+diam./2.0] .= rho_air.+drho          # initial density of the conduit overlapping with the circular chamber (so in the chamber)
    #rho[inpolygon(x2dc,y2dc,Xp2,Yp2) .== 1.0 .&& y2dc.>=locY.+diam./2.0] .= #=.-0.0 .*=# (y2dc[inpolygon(x2dc,y2dc,Xp2,Yp2).==1.0 .&& y2dc.>=locY.+diam./2.0].+0.15.*Ly).*drho./(-0.15.-(locY.+diam./2.0)).+rho_air    # initial density in the conduit

    #P[radrho .< diam ./ 2.0]         .= (P0 .* rho[radrho .< diam ./ 2.0]) ./ rho_air
    #P         .+= #=((P0 .* rho) ./ rho_air .* maskrho_air) .+ =#((P0_solid .+  (1.0 ./ β) .* log.(rho ./ rho_solid)) .* maskrho_solid) #P0 ./ rho_air .* rho                                                              # equation of state for pressure depending on density

    @. η += eta_air * (maskrho_air[2:end-1, 2:end-1] == 1.0) + eta_solid * (maskrho_solid[2:end-1, 2:end-1] == 1.0)     # initial viscosity distribution
    @. μ += mu_air * (maskrho_air[2:end-1, 2:end-1] == 1.0) + mu_solid * (maskrho_solid[2:end-1, 2:end-1] == 1.0)       # initial viscosity distribution

    @. η_c += eta_air * (maskμc_air == 1.0) + eta_solid * (maskμc_solid == 1.0)   # initial viscosity distribution for corner nodes
    @. μ_c += mu_air * (maskμc_air == 1.0) + mu_solid * (maskμc_solid == 1.0)     # initial viscosity distribution for corner nodes

    time = 0.0

    # Inital plot    
    fig = Figure()
    ax = Axis(fig[1, 1], xticks=([-500, 0, 500], ["-500", "0", "500"]), yticks=([-300, 0, 700], ["-300", "0", "700"]),
                yticklabelsize=25, xticklabelsize=25, xlabelsize=25, ylabelsize=25)#, title="0")

    x_circ = zeros(length(x_circ_ind))
    y_circ = zeros(length(y_circ_ind))
    x_circ = Xc[x_circ_ind]
    y_circ = Yc[y_circ_ind]

    U = av_y(Vx)
    V = av_x(Vy)
    data_plt = sqrt(U.^2.0 .+ V.^2.0)
    
    p1 = heatmap!(ax, x2dc, y2dc, P, shading=false, colormap=Reverse(:roma))
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
    scatter!(ax, x_circ, y_circ, color = :white, markersize=4.0) # conduit and chamber
    display(fig)
    #save("./Plots/Shockwave_1/0.png", fig)

    # Solver
    #Vy        .= 0.0 .* Vx0 ./ 1.0 .* exp.(-1.0e3 * ((x2dVy .- (locX .+ 1.5 .* diam)).^2.0 .+ y2dVy.^2.0)) # initial velocity in y-direction
    #Mx        .= av_x(rho).*Vx  # initial momentum in x-direction
    #My        .= av_y(rho).*Vy  # initial momentum in y-direction
    #rho_ini   .= sum(rho[2:end-1]);
    #Mx_ini    .= sum(Mx);
    #rhoRes    .= zero.(rho[2:end-1,2:end-1])      # residual of density

    # history variables for updating stresses, momentum and density between timesteps without influence of pseudo transient iteration changed in the variables

    nan = false

    reset_timer!()
    for it = 1:350                                                   # 60 iterations max for all equilibrated
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
        while err > 1.0e-3
            iter += 1
            #global Vx = Vx
            #global Vy = Vy
            #global dt = dt
            @timeit "compute beta" beta_vec .= 1.0 ./ P 
            @timeit "compute dt" dt = minimum([dx./maximum(abs.(Vx)), dy./maximum(abs.(Vy)), min(dx, dy) .* sqrt(maximum(rho .* beta_vec)), min(dx, dy).^2.0 ./ maximum(η)]).*4.0 # time step size
            @timeit "compute c_loc" c_loc .= 1.0./sqrt.(rho[2:end-1,2:end-1].* beta_vec[2:end-1,2:end-1])                                                    # local speed of sound
            @timeit "compute dtPT" dtPT .= min.(min.(min.(dx ./ abs.(av_x(Vx[:, 2:end-1])), dx ./ abs.(av_y(Vy[2:end-1, :]))), min(dx, dy).^2.0 ./ η), dx ./ c_loc) # time step size for pressure and temperature
            @timeit "compute dtrho" dtrho .= 1.0 ./ (1.0 ./ dt .+ 1.0 ./ (min(dx, dy) ./ c_loc ./ 4.1))                                                         # time step size for density
            
            # Conservation of mass
            @timeit "compute Frhox" Frhox .= (Vx .> 0.0).*Vx.*rho[1:end-1, :] .+ (Vx .< 0.0).*Vx.*rho[2:end, :] # mass flux in x-direction (upwind scheme)
            @timeit "compute Frhoy" Frhoy .= (Vy .> 0.0).*Vy.*rho[:, 1:end-1] .+ (Vy .< 0.0).*Vy.*rho[:, 2:end] # mass flux in y-direction (upwind scheme)
            @timeit "compute drhodt" drhodt .= (rho .- rho_old)./dt                                                                           # time derivative of density
            @timeit "compute rhoRes1" rhoRes .= .-drhodt[2:end-1, 2:end-1] .- diff(Frhox[:, 2:end-1], dims=1)./dx .- diff(Frhoy[2:end-1, :], dims=2)./dy        # updating residual of density
            @timeit "compute rhoRes2" rhoRes .= rhoRes # .*maskrho_air[2:end-1, 2:end-1]                                                                               # applying mask to residual of density
            @timeit "compute rho1" rho[2:end-1, 2:end-1] .= rho[2:end-1, 2:end-1] .+ rhoRes.*dtrho.*CFL_P                   # updating density
            
            if maximum(isnan.(rho))
                print("rho has NaN\n")
                @show rho
                nan = true
                break
            end
            
            # Boundary conditions Inflow and outflow densities
            #rho[1, :] .= rho[2, :]
            #rho[end, :] .= rho[end-1, :]
            # Boundary conditions impermeable walls
            #rho[:, 1] .= rho[:, 2]
            #rho[:, end] .= rho[:, end-1]

            # Strain-rates and stresses
            P .= zeros(Float64, nx + 2, ny + 2)
            #@timeit "compute P" P              .+= ((P0 ./ rho_air .* rho) .* maskrho_air) .+ ((P0_solid .+  (1.0 ./ beta_vec) .* log.(rho ./ rho_solid)) .* maskrho_solid)         # equation of state for pressure depending on density                                                                   # equation of state for air
            @timeit "compute P" P              .+= (P0 ./ rho_air .* rho)                                                # equation of state for pressure depending on density                                                                   # equation of state for air

            if maximum(isnan.(P))
                print("P has NaN\n")
                @show P
                nan = true
                break
            end

            @timeit "compute divV" divV        .= diff(Vx[:,2:end-1],dims=1)./dx .+ diff(Vy[2:end-1,:],dims=2)./dy                                              # divergence of velocity
            @timeit "compute Exx" Exx          .= diff(Vx[:,2:end-1],dims=1)./dx .- 1.0./3.0.*divV                                                              # strain-rate in x-direction
            @timeit "compute Eyy" Eyy          .= diff(Vy[2:end-1,:],dims=2)./dy .- 1.0./3.0.*divV                                                              # strain-rate in y-direction
            @timeit "compute Ezz" Ezz          .=                                .- 1.0./3.0.*divV                                                              # strain-rate in z-direction
            @timeit "compute Exy" Exy          .= 0.5.*(diff(Vy,dims=1)./dx .+ diff(Vx,dims=2)./dy)                                                             # shear strain-rate in xy-direction
            
            @timeit "compute Sxx" Sxx          .= .-P[2:end-1,2:end-1] .+ 2.0.*η.*Exx                                                                          # total stress (dani class 5 equation)
            @timeit "compute Syy" Syy          .= .-P[2:end-1,2:end-1] .+ 2.0.*η.*Eyy                                                                          # total stress
            @timeit "compute Sxy" Sxy          .=                         2.0.*η_c.*Exy                                                                            # total stress
            @timeit "compute Szz" Szz          .= .-P[2:end-1,2:end-1] .+ 2.0.*η.*Ezz                                                                           # stress in z-direction
            
            #=
            @timeit "compute Sxx" Sxx          .= .-P[2:end-1,2:end-1] .+ 2.0 .* (1.0 ./ ((1.0 ./ η) .+ (1.0 ./ (μ .* dt)))) .* (Exx .+ ((Sxx_old .+ P_old[2:end-1, 2:end-1]) ./ (2.0 .* μ .* dt))) #.*μ.*  Exx .- (μ   ./ η) .*   (Sxx_old .+ P[2:end-1,2:end-1])                                                                          # stress in x-direction
            @timeit "compute Syy" Syy          .= .-P[2:end-1,2:end-1] .+ 2.0 .* (1.0 ./ ((1.0 ./ η) .+ (1.0 ./ (μ .* dt)))) .* (Eyy .+ ((Syy_old .+ P_old[2:end-1, 2:end-1]) ./ (2.0 .* μ .* dt))) #.*μ.*  Eyy .- (μ   ./ η) .*   (Syy_old .+ P[2:end-1,2:end-1])                                             # stress in y-direction
            @timeit "compute Sxy" Sxy          .=                         2.0 .* (1.0 ./ ((1.0 ./ η_c) .+ (1.0 ./ (μ_c .* dt)))) .* (Exy .+ (Sxy_old                          ./ (2.0 .* μ_c .* dt))) #.*μ_c.*Exy .- (μ_c ./ η_c) .* (Sxy_old .+ av_xy(P))                                                  # stress in xy-direction
            @timeit "compute Szz" Szz          .= .-P[2:end-1,2:end-1] .+ 2.0 .* (1.0 ./ ((1.0 ./ η) .+ (1.0 ./ (μ .* dt)))) .* (Ezz .+ ((Szz_old .+ P_old[2:end-1, 2:end-1]) ./ (2.0 .* μ .* dt))) #μ.*  Ezz .- (μ   ./ η) .*   (Szz_old .+ P[2:end-1,2:end-1])                                             # stress in z-direction
            =#
            @timeit "compute dtV" dtV           = 1.0 ./ (1.0 ./ dt .+ 1.0 ./ (min(dx,dy).^2.0 ./ maximum(η) ./ 4.0)) .* CFL_V                                           # time step size for velocity
            
            # Conservation of the x-component of momentum
            @timeit "compute Mx" Mx             .= av_x(rho).*Vx                                                             # momentum in x-direction
            
            if maximum(isnan.(Mx))
                print("Mx has NaN\n")
                nan = true
                break
            end

            @timeit "compute FMxx" FMxx         .= (av_x(Vx[ :     ,2:end-1]).> 0.0).*av_x(Vx[ :     ,2:end-1]).*Mx[1:end-1,2:end-1] .+ (av_x(Vx[ :     ,2:end-1]).< 0.0).*av_x(Vx[ :     ,2:end-1]).*Mx[2:end  ,2:end-1]  # mass flux in x-direction (upwind scheme)
            @timeit "compute FMxy" FMxy         .= (av_x(Vy[2:end-1, :     ]).> 0.0).*av_x(Vy[2:end-1, :     ]).*Mx[2:end-1,1:end-1] .+ (av_x(Vy[2:end-1, :     ]).< 0.0).*av_x(Vy[2:end-1, :     ]).*Mx[2:end-1,2:end  ]  # mass flux in y-direction (upwind scheme)
            @timeit "compute dMxdt" dMxdt       .= (Mx.-Mx_old)./dt                                                                                             # time derivative of momentum in x-direction
            @timeit "compute MxRes1" MxRes      .= .-dMxdt[2:end-1,2:end-1] .- diff((FMxx .- Sxx),dims=1)./dx .- diff(FMxy .- Sxy[2:end-1,:],dims=2)./dy       # updating residual of momentum in x-direction
            @timeit "compute MxRes2" MxRes      .= MxRes # .*maskVx_air[2:end-1,2:end-1]                                                                              # applying mask to residual of momentum in x-direction
            @timeit "compute dMxdtau" dMxdtau   .= MxRes .+ dMxdtau .* ksi                                                                                    # stress derivative of momentum in x-direction
            @timeit "compute Mx2" Mx[2:end-1,2:end-1]  .= Mx[2:end-1,2:end-1] .+ dMxdtau.*av_x(dtPT).*CFL_V                                                     # updating momentum in x-direction
            @timeit "compute Vx" Vx           .= Mx./av_x(rho)                                                              # velocity in x-direction

            # BC fixed walls (normal velocity = 0)
            Mx[:,1]      .= Mx[:,2]
            Mx[:,end]    .= Mx[:,end-1]
            # BC no slip on vertical walls
            #Mx[1,:]      .= Mx[2,:]
            #Mx[end,:]    .= Mx[end-1,:]
            
            Vx[1,:]                           .= Vx0
            Vx[end,:]                         .= Vx0
            
            # Conservation of the y component of momentum
            @timeit "compute My" My             .= av_y(rho).*Vy                                                              # momentum in y-direction
            
            if maximum(isnan.(My))
                print("My has NaN\n")
                nan = true
                break
            end
            
            @timeit "compute FMyy" FMyy         .= (av_y(Vy[2:end-1, :     ]).> 0.0).*av_y(Vy[2:end-1, :     ]).*My[2:end-1,1:end-1] .+ (av_y(Vy[2:end-1, :     ]).< 0.0).*av_y(Vy[2:end-1, :     ]).*My[2:end-1,2:end  ]  # mass flux in y-direction
            @timeit "compute FMyx" FMyx         .= (av_y(Vx[ :     ,2:end-1]).> 0.0).*av_y(Vx[ :     ,2:end-1]).*My[1:end-1,2:end-1] .+ (av_y(Vx[ :     ,2:end-1]).< 0.0).*av_y(Vx[ :     ,2:end-1]).*My[2:end  ,2:end-1]  # mass flux in x-direction
            @timeit "compute dMydt" dMydt       .= (My-My_old)./dt                                                                                              # time derivative of momentum in y-direction
            @timeit "compute MyRes1" MyRes      .= .-dMydt[2:end-1,2:end-1] .- diff(FMyy .- Syy,dims=2)./dy .- diff(FMyx .- Sxy[:,2:end-1],dims=1)./dx - g .* av_y(rho[2:end-1,2:end-1])  # letzter Term war vorher av_y(rho[2:end-1,2:end-1])
            @timeit "compute MyRes2" MyRes      .= MyRes # .*maskVy_air[2:end-1,2:end-1]                                                                              # applying mask to residual of momentum in y-direction
            @timeit "compute dMydtau" dMydtau   .= MyRes .+ dMydtau.*ksi                                                                                      # stress derivative of momentum in y-direction
            @timeit "compute My2" My[2:end-1,2:end-1]  .= My[2:end-1,2:end-1] .+ dMydtau.*av_y(dtPT).*CFL_V                                                     # updating momentum in y-direction
            
            @timeit "compute Vy" Vy             .= My./av_y(rho)                    # updating velocity in y-direction
            
            # BC fixed walls (normal velocity = 0)
            My[1,:]      .= -My[2,:]
            My[end,:]    .= -My[end-1,:]
            # BC no slip on horizontal walls
            My[:,1]      .= My[:,2]
            My[:,end]    .= My[:,end-1]



            if mod(iter,1) == 0
                it_counter += 1
                
                @timeit "compute err" err = maximum(abs.([rhoRes[:]; MxRes[:]; MyRes[:]]))                                                                      # error for time integration concerning density and momentum
                print("iter = $iter, err = $err, rhoRes = $(maximum(abs.(rhoRes[:]))), MxRes = $(maximum(abs.(MxRes[:]))), MyRes = $(maximum(abs.(MyRes[:])))\n")
                if isnan(err) == 1 
                    break
                end
            end

            if err <= 1.0e-3
                @timeit "show it_counter" @show it_counter
            end
        end
        #global time = time
        @timeit "compute time" time = time + dt

        # Updating plot
        if mod(it-1, 1) == 0
            fig1 = Figure()
            ax1 = Axis(fig1[1,1], xticks=([-500, 0, 500], ["-500", "0", "500"]), yticks=([-300, 0, 700], ["-300", "0", "700"]),
                    yticklabelsize=25, xticklabelsize=25, xlabelsize=25, ylabelsize=25, title="time = $time")

            X = av_xy(x2dc)
            Y = av_xy(y2dc)
            U = av_y(Vx)
            V = av_x(Vy)

            data_plt = sqrt.(U.^2.0 .+ V.^2.0)

            #hm = heatmap!(ax1, x2dc, y2dc, P, shading=false, colorrange=(P0, P0*2))
            hm = heatmap!(ax1, X, Y, data_plt, shading=false, colorrange=(0.0, 350), colormap=Reverse(:roma))
            #Colorbar(fig1[1,2],  hm, label="Pressure", labelsize=25, ticklabelsize=25)
            Colorbar(fig1[1,2],  hm, label="Velocity", labelsize=25, ticklabelsize=25)
            #lines!(ax1, Xp, Yp, color = :white)    # only conduit
            scatter!(ax1, x_circ, y_circ, color = :white, markersize=4.0)     # conduit and chamber
            #poly!(ax,points_XpYp)
            
            #stepsize = 10
            #arrows!(ax1, X[1:stepsize:end, 1], Y[1, 1:stepsize:end], U[1:stepsize:end, 1:stepsize:end], V[1:stepsize:end, 1:stepsize:end], arrowsize=7, color = :white)
            display(fig1)

            #save("./Plots/Shockwave_1/$(it).png", fig1)
        end

        if nan
            break
        end
    end
    print_timer()
    return
end
