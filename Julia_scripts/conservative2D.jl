using CairoMakie

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
    A = 0.5 .* (B[2:end,:] .+ B[1:end-1,:])
end

function av_y(B)
    A = 0.5 .* (B[:,2:end] .+ B[:,1:end-1])
end

function av_xy(B)
    A = 0.25 .* (B[2:end,2:end] .+ B[1:end-1,2:end] .+ B[2:end,1:end-1] .+ B[1:end-1,1:end-1])
end

function advection2d(A, Vx, Vy, dx, dy, dt)
    dtadv = dt
    A[1:end-1,:] = A[1:end-1,:] .- (Vx[2:end-1,:] .< 0.0) .* Vx[2:end-1,:] .* diff(A, dims=1) ./ dx .* dtadv
    A[2:end,:] = A[2:end,:] .- (Vx[2:end-1,:] .> 0.0) .* Vx[2:end-1,:] .* diff(A, dims=1) ./ dx .* dtadv
    A[:,1:end-1] = A[:,1:end-1] .- (Vy[:,2:end-1] .< 0.0) .* Vy[:,2:end-1] .* diff(A, dims=2) ./ dy .* dtadv
    A[:,2:end] = A[:,2:end] .- (Vy[:,2:end-1] .> 0.0) .* Vy[:,2:end-1] .* diff(A, dims=2) ./ dy .* dtadv
    return A
end

function conservative2D()
    # Physics
    Lx      = 1.0e4                             # length of domain in x-direction
    Ly      = Lx                                # length of domain in y-direction
    rho0    = 1.225                             # density at rest
    drho    = 3.0e2                             # density difference
    Vx0     = 0                                 # starting velocity in x-direction
    P0      = 1.0e5                             # pressure at rest
    beta    = 1.0/141.0e3                       # compressibility
    mu      = 1.81e-5                           # dynamic viscosity
    g       = 9.81                              # gravitational acceleration

    # Numerics
    nx      = 151                               # number of nodes in x-direction
    ny      = 151                               # number of nodes in y-direction
    dx      = Lx/(nx)                           # grid spacing in x-direction
    dy      = Ly/(ny)                           # grid spacing in y-direction
    dt      = 6.0e-2                            # time step size
    CFL_P   = 1.0/4.0                           # Courant number for pressure
    CFL_V   = 1.0/16.0                          # Courant number for velocity
    ksi     = 0.0*0.95                          # relaxation factor for stress

    # Plots
    st      = ceil(Int, nx / 30)

    # Reduce allocations
    c_loc   = zeros(Float64, nx, ny)
    dtPT    = zeros(Float64, nx, ny)
    dtrho   = zeros(Float64, nx, ny)
    Frhox   = zeros(Float64, nx + 1, ny + 2)
    Frhoy   = zeros(Float64, nx + 2, ny + 1)
    drhodt  = zeros(Float64, nx + 2, ny + 2)
    rhoRes  = zeros(Float64, nx, ny)
    rho     = zeros(Float64, nx + 2, ny + 2)
    P       = zeros(Float64, nx + 2, ny + 2)
    divV    = zeros(Float64, nx, ny)
    Exx     = zeros(Float64, nx, ny)
    Eyy     = zeros(Float64, nx, ny)
    Ezz     = zeros(Float64, nx, ny)
    Exy     = zeros(Float64, nx + 1, ny + 1)
    Sxx     = zeros(Float64, nx, ny)
    Syy     = zeros(Float64, nx, ny)
    Szz     = zeros(Float64, nx, ny)
    Sxy     = zeros(Float64, nx + 1, ny + 1)
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
    maskrho = zeros(Float64, nx, ny)
    maskVx  = zeros(Float64, nx, ny)
    maskVy  = zeros(Float64, nx, ny)

    iters_evo = Float64[]
    errs_evo = Float64[]

    # Initialization
    Xv      = range(-Lx/2.0, stop=Lx/2.0, step=dx)                      # x-coordinates of velocity nodes
    Xc      = range(-(Lx+dx)/2.0, stop=(Lx+dx)/2.0, step=dx)            # x-coordinates of density nodes
    Yv      = range(-Ly/2.0, stop=Ly/2.0, step=dy)                      # y-coordinates of velocity nodes
    Yc      = range(-(Ly+dy)/2.0, stop=(Ly+dy)/2.0, step=dy)            # y-coordinates of density nodes
    x2dc, y2dc = meshgrid(Xc,Yc)                                        # 2d mesh of x- and y-coordinates of density nodes
    x2dVx, y2dVx = meshgrid(Xv,Yc)                                      # 2d mesh of x- and y-coordinates of velocity nodes in x-direction
    x2dVy, y2dVy = meshgrid(Xc,Yv)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction

    Vx        .= Vx0.*ones(Float64, nx + 1,ny + 2)                               # initial velocity in x-direction
    Vx[1,:]   .= Vx0 
    Vx[end,:] .= Vx0
    locX      = 0.0                                                     # x-coordinate of circle
    locY      = -0.35*Ly                                                # y-coordinate of circle
    diam      = 0.2*Lx                                                  # diameter of circle
    radrho    .= sqrt.((x2dc .-locX).^2.0 .+ (y2dc .-locY).^2.0)        # radius of circle
    radVx     .= sqrt.((x2dVx.-locX).^2.0 .+ (y2dVx.-locY).^2.0)        # radius of circle
    radVy     .= sqrt.((x2dVy.-locX).^2.0 .+ (y2dVy.-locY).^2.0)        # radius of circle
    Xp        = [-1.0/2.0, -1.0/2.0, -0.3, -1.0/8.0, -0.02, -0.02, 0.02, 0.02, 1.0/8.0, 0.3, 1.0/2.0, 1.0/2.0].*Lx          # x-coordinates of polygon
    Yp        = [-1.0/2.0, -1.0/5.0, -1.0/5.0, -0.1, -0.15,  -0.35, -0.35, -0.15,  -0.1, -1.0/5.0, -1.0/5.0, -1.0/2.0].*Ly  # y-coordinates of polygon
    Xp2       = [-0.02, -0.02,  0.02,  0.02].*Lx                                                                            # x-coordinates of polygon
    Yp2       = [-0.1,  -0.35, -0.35, -0.1].*Ly                                                                             # y-coordinates of polygon

    inpolyrho = inpolygon(x2dc,y2dc,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon
    inpolyVx = inpolygon(x2dVx,y2dVx,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon
    inpolyVy = inpolygon(x2dVy,y2dVy,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon

    maskrho = 1.0 .- inpolyrho                                             # mask for density
    maskrho[radrho .< diam ./ 2.0] .= 1.0                               # mask for density in circle
    maskVx  = 1.0 .- inpolyVx                               # mask for velocity in x-direction
    maskVx[radVx .< diam ./ 2.0] .= 1.0                                 # mask for velocity in x-direction in circle
    maskVy  = 1.0 .- inpolyVy                               # mask for velocity in y-direction
    maskVy[radVy .< diam ./ 2.0] .= 1.0                                 # mask for velocity in y-direction in circle

    P         = P0.*exp.(-g.*(y2dc.+1.0./5.0.*Ly).*0.028./288.0./8.314) # barometric P setting atmosphere: P = P0*exp(-(g*h*M)/(T*R)) M: Mass density of air, T: Temperature, R: Gas constant, h: height, P0: Pressure at sea level
    rho       = rho0.*exp.(beta.*(P.-P0))             # equation of state for density depending on pressure
    rho[radrho .< diam ./ 2.0] .= rho0 .+ drho        # initial density in the circle

    rho[inpolygon(x2dc,y2dc,Xp2,Yp2) .== 1.0 .&& y2dc.< locY.+diam./2.0] .= rho0.+drho        # initial density below the circle
    rho[inpolygon(x2dc,y2dc,Xp2,Yp2) .== 1.0 .&& y2dc.>=locY.+diam./2.0] .= -0.0 .* (y2dc[inpolygon(x2dc,y2dc,Xp2,Yp2).==1.0 .&& y2dc.>=locY.+diam./2.0].+0.1.*Ly).*drho./(-0.1.-(locY.+diam./2.0)).+rho0                       # initial density above the circle
    P         .= 1.0./beta.*log.(rho./rho0).+P0                         # equation of state for pressure depending on density

    #=global=# time = 0.0

    # Inital plot
    fig = Figure()
    ax  =   (ax_P   = Axis(fig[1, 1][1 ,1]; aspect = DataAspect(), title="P"),
            ax_V    = Axis(fig[1, 2][1 ,1]; aspect = DataAspect(), title="V"),
            ax_open = Axis(fig[2, 1][1 ,1]; aspect = DataAspect(), title="time=$time"),
            ax_err  = Axis(fig[2, 2]; yscale=log10, title="Convergence", xlabel="# iter/nx", ylabel="error"),)

    plt =   (fld = (P = heatmap!(ax.ax_P, x2dc, y2dc, P; shading=false, colorrange=(P0, P0*2)),
                        #lines!(ax.ax_P, Xp, Yp; color = :white),
                        #arrows!(ax, X[1:stepsize:end, 1], Y[1, 1:stepsize:end], U[1:stepsize:end, 1:stepsize:end], V[1:stepsize:end, 1:stepsize:end], arrowsize=7, color = :white),
                    V = heatmap!(ax.ax_V, x2dVy, y2dVy, Vy; shading=false, colorrange=(0.0, 350), colormap=Reverse(:roma)),
                        #lines!(ax.ax_V, Xp, Yp; color = :white),
                        #arrows!(ax, X[1:stepsize:end, 1], Y[1, 1:stepsize:end], U[1:stepsize:end, 1:stepsize:end], V[1:stepsize:end, 1:stepsize:end], arrowsize=7, color = :white),
                    open = heatmap!(ax.ax_open, x2dc, y2dc, P; shading=false, colormap=:turbo, colorrange=(P0, P0*2)),
                        #lines!(ax.ax_open, Xp, Yp; color = :white),
                    err = lines!(ax.ax_err, Point2.(iters_evo, errs_evo), linewidth=4)),)
    Colorbar(fig[1 ,1][1,2], plt.fld.P, label="Pressure [Pa]")
    Colorbar(fig[1 ,2][1,2], plt.fld.V, label="Velocity [m/s]")
    Colorbar(fig[2 ,1][1,2], plt.fld.open, label="open")
    points_XpYp = Point2f[]
    for i in eachindex(Xp)
        x = Xp[i]
        y = Yp[i]
        point = (x,y)
        push!(points_XpYp, point)
    end
    #poly!(ax,points_XpYp)
    lines!(ax.ax_P, Xp, Yp, color = :white)
    lines!(ax.ax_V, Xp, Yp, color = :white)
    lines!(ax.ax_open, Xp, Yp, color = :white)
    display(fig)
    
    # Sovler
    #Vy        .= 0.0 * Vx0 ./ 1.0 .* exp.(-1.0e3 * ((x2dVy .- (locX .+ 1.5 .* diam)).^2.0 .+ y2dVy.^2.0)) # initial velocity in y-direction
    #Mx        .= av_x(rho).*Vx  # initial momentum in x-direction
    #My        .= av_y(rho).*Vy  # initial momentum in y-direction
    #rhoRes    .= zero.(rho[2:end-1,2:end-1])      # residual of density

    @inbounds for it = 1:300
        @show it
        rho_old .= rho                                               # save old density values for time integration
        Mx_old .= Mx                                                # save old momentum values in x-direction for time integration
        My_old .= My                                                # save old momentum values in y-direction for time integration
        err = 1.0                                                   # error for time integration
        iter = 0                                                    
        it_counter = 0
        dMxdtau .= 0.0                                              # stress derivative of momentum in x-direction
        dMydtau .= 0.0                                              # stress derivative of momentum in y-direction        
        while err > 1.0e-3
            iter += 1
            
            dt = minimum([dx/maximum(abs.(Vx)), dy/maximum(abs.(Vy)), min(dx, dy) .* sqrt(maximum(rho).*beta), min(dx, dy).^2.0./mu]).*8.0 # time step size
            c_loc .= 1.0./sqrt.(rho[2:end-1,2:end-1].*beta)                                                    # local speed of sound
            dtPT .= min.(min.(min.(dx ./ abs.(av_x(Vx[:, 2:end-1])), dx ./ abs.(av_y(Vy[2:end-1, :]))), min(dx, dy).^2.0 ./ mu .* ones(Float64, nx, ny)), dx ./ c_loc) # time step size for pressure and temperature
            dtrho .= 1.0 ./ (1.0 ./ dt .+ 1.0 ./ (min(dx, dy) ./ c_loc ./ 4.1))                                                         # time step size for density
            
            # Conservation of mass
            Frhox .= (Vx .> 0.0).*Vx.*rho[1:end-1, :] .+ (Vx .< 0.0).*Vx.*rho[2:end, :] # mass flux in x-direction
            Frhoy .= (Vy .> 0.0).*Vy.*rho[:, 1:end-1] .+ (Vy .< 0.0).*Vy.*rho[:, 2:end] # mass flux in y-direction
            drhodt .= (rho .- rho_old)./dt                                                                           # time derivative of density
            rhoRes .= .-drhodt[2:end-1, 2:end-1] .- diff(Frhox[:, 2:end-1], dims=1)./dx .- diff(Frhoy[2:end-1, :], dims=2)./dy        # updating residual of density
            rhoRes .= rhoRes.*maskrho[2:end-1, 2:end-1]                                                                               # applying mask to residual of density
            rho[2:end-1, 2:end-1] .= rho[2:end-1, 2:end-1] .+ rhoRes.*dtrho.*CFL_P                   # updating density
            # Boundary conditions Inflow and outflow densities
            rho[1, :] .= rho[2, :]
            rho[end, :] .= rho[end-1, :]
            # Boundary conditions impermeable walls
            rho[:, 1] .= rho[:, 2]
            rho[:, end] .= rho[:, end-1]

            # Strain-rates and stresses
            P              .= 1.0 ./ beta .* log.(abs.(rho) ./ rho0) .+ P0                                                  # equation of state for pressure depending on density
            divV        .= diff(Vx[:,2:end-1], dims=1)./dx .+ diff(Vy[2:end-1,:],dims=2)./dy                                             # divergence of velocity
            Exx          .= diff(Vx[:,2:end-1],dims=1)./dx .- 1.0./3.0.*divV                                                               # strain-rate in x-direction
            Eyy          .= diff(Vy[2:end-1,:],dims=2)./dy .- 1.0./3.0.*divV                                                               # strain-rate in y-direction
            Ezz          .=                                .- 1.0./3.0.*divV                                                                      # strain-rate in z-direction
            Exy          .= 0.5.*(diff(Vy,dims=1)./dx .+ diff(Vx,dims=2)./dy)                                                             # shear strain-rate in xy-direction
            Sxx          .= .-P[2:end-1,2:end-1] .+ 2.0.*mu.*Exx                                                                          # stress in x-direction
            Syy          .= .-P[2:end-1,2:end-1] .+ 2.0.*mu.*Eyy                                                                          # stress in y-direction
            Sxy          .=                         2.0.*mu.*Exy                                                                            # stress in xy-direction
            Szz          .= .-P[2:end-1,2:end-1] .+ 2.0.*mu.*Ezz                                                                           # stress in z-direction
            dtV          = 1.0 ./ (1.0 ./ dt .+ 1.0 ./ (min(dx,dy).^2.0 ./ mu ./ 4.0)) .* CFL_V                                           # time step size for velocity
            
            # Conservation of the x-component of momentum
            Mx           .= av_x(rho).*Vx                                                             # momentum in x-direction
            FMxx         .= (av_x(Vx[ :     ,2:end-1]).> 0.0).*av_x(Vx[ :     ,2:end-1]).*Mx[1:end-1,2:end-1] .+ (av_x(Vx[ :     ,2:end-1]).< 0.0).*av_x(Vx[ :     ,2:end-1]).*Mx[2:end  ,2:end-1]  # mass flux in x-direction
            FMxy         .= (av_x(Vy[2:end-1, :     ]).> 0.0).*av_x(Vy[2:end-1, :     ]).*Mx[2:end-1,1:end-1] .+ (av_x(Vy[2:end-1, :     ]).< 0.0).*av_x(Vy[2:end-1, :     ]).*Mx[2:end-1,2:end  ]  # mass flux in y-direction
            dMxdt       .= (Mx.-Mx_old)./dt                                                                                             # time derivative of momentum in x-direction
            MxRes       .= .-dMxdt[2:end-1,2:end-1] .- diff((FMxx .- Sxx),dims=1)./dx .- diff(FMxy .- Sxy[2:end-1,:],dims=2)./dy       # updating residual of momentum in x-direction
            MxRes       .= MxRes.*maskVx[2:end-1,2:end-1]                                                                              # applying mask to residual of momentum in x-direction
            dMxdtau     .= MxRes .+ dMxdtau .* ksi                                                                                    # stress derivative of momentum in x-direction
            Mx[2:end-1,2:end-1]  .= Mx[2:end-1,2:end-1] .+ dMxdtau.*av_x(dtPT).*CFL_V                                                     # updating momentum in x-direction
            # BC fixed walls (normal velocity = 0)
            Mx[:,1]      .= Mx[:,2]
            Mx[:,end]    .= Mx[:,end-1]
            # BC no slip on vertical walls
            # Mx[1,:]      .= Mx[2,:]
            # Mx[end,:]    .= Mx[end-1,:]
            
            Vx           .= Mx./av_x(rho)                                                              # velocity in x-direction
            Vx[1,:]                           .= Vx0
            Vx[end,:]                         .= Vx0
            
            # Conservation of the y component of momentum
            My           .= av_y(rho).*Vy                                                              # momentum in y-direction
            FMyy         .= (av_y(Vy[2:end-1, :     ]).> 0.0).*av_y(Vy[2:end-1, :     ]).*My[2:end-1,1:end-1] .+ (av_y(Vy[2:end-1, :     ]).< 0.0).*av_y(Vy[2:end-1, :     ]).*My[2:end-1,2:end  ]  # mass flux in y-direction
            FMyx         .= (av_y(Vx[ :     ,2:end-1]).> 0.0).*av_y(Vx[ :     ,2:end-1]).*My[1:end-1,2:end-1] .+ (av_y(Vx[ :     ,2:end-1]).< 0.0).*av_y(Vx[ :     ,2:end-1]).*My[2:end  ,2:end-1]  # mass flux in x-direction
            dMydt       .= (My-My_old)./dt                                                                                              # time derivative of momentum in y-direction
            MyRes       .= .-dMydt[2:end-1,2:end-1] .- diff(FMyy .- Syy,dims=2)./dy .- diff(FMyx .- Sxy[:,2:end-1],dims=1)./dx - g .* av_y(rho[2:end-1,2:end-1])  # letzter Term war vorher av_y(rho[2:end-1,2:end-1])
            MyRes       .= MyRes.*maskVy[2:end-1,2:end-1]                                                                              # applying mask to residual of momentum in y-direction
            dMydtau     .= MyRes .+ dMydtau.*ksi                                                                                      # stress derivative of momentum in y-direction
            My[2:end-1,2:end-1]  .= My[2:end-1,2:end-1] .+ dMydtau.*av_y(dtPT).*CFL_V                                                     # updating momentum in y-direction
            Vy           .= My./av_y(rho)                                                                # updating velocity in y-direction
            # BC fixed walls (normal velocity = 0)
            My[1,:]      .= -My[2,:]
            My[end,:]    .= -My[end-1,:]
            # BC no slip on horizontal walls
            My[:,1]      .= My[:,2]
            My[:,end]    .= My[:,end-1]

            if mod(iter, 25) == 0
                it_counter += 1
                print("err = $err\n")
                err = maximum(abs.([rhoRes[:]; MxRes[:]; MyRes[:]]))                                                                      # error for time integration concerning density and momentum
                if isnan(err) == 1 
                    break
                end
            end
            
            push!(iters_evo, iter/nx)
            push!(errs_evo, err)

            if err <= 1.0e-3
                @show it_counter
            end
        end
        time = time + dt

        # Updating plots
        if mod(it-1, 1) == 0
            X = av_xy(x2dc)
            Y = av_xy(y2dc)
            U = av_y(Vx)
            V = av_x(Vy)
            stepsize = 5
            
            plt.fld.P[3] = P
            plt.fld.V[3] = Vy
            plt.fld.open[3] = P
            plt.fld.err[1] = Point2.(iters_evo, errs_evo)

            #=
            fig = Figure()
            ax = (ax_P = Axis(fig[1, 1][1 ,1]; aspect = DataAspect(), title="P"),
                  ax_V = Axis(fig[1, 2][1 ,1]; aspect = DataAspect(), title="V"),
                  ax_open = Axis(fig[2, 1][1 ,1]; aspect = DataAspect(), title="time=$time")
                  ax_err = Axis(fig[2, 2]; aspect = DataAspect(), title="err", xlabel="# iter/nx", ylabel="error"),)

            plt = (hm_P = heatmap!(ax.ax_P, x2dc, y2dc, P, shading=false, colorrange=(P0, P0*2)),
                          lines!(ax.ax_P, Xp, Yp, color = :white),
                          arrows!(ax, X[1:stepsize:end, 1], Y[1, 1:stepsize:end], U[1:stepsize:end, 1:stepsize:end], V[1:stepsize:end, 1:stepsize:end], arrowsize=7, color = :white),
                   hm_V = heatmap!(ax.ax_V, x2dVy, y2dVy, Vy, shading=false, colorrange=(0.0, 350), colormap=:viridis),
                          lines!(ax.ax_V, Xp, Yp, color = :white),
                          arrows!(ax, X[1:stepsize:end, 1], Y[1, 1:stepsize:end], U[1:stepsize:end, 1:stepsize:end], V[1:stepsize:end, 1:stepsize:end], arrowsize=7, color = :white),
                   hm_open = heatmap!(ax.ax_open, x2dc, y2dc, P, shading=false, colormap=:turbo),
                             lines!(ax.ax_open, Xp, Yp, color = :white),
                   err = scatterlines!(ax.ax_err, Point2.(iters_evo, errs_evo), linewidth=4))
            Colorbar(fig[1 ,1][1,2], plt.hm_P, label="Pressure [Pa]")
            Colorbar(fig[1 ,2][1,2], plt.hm_V, label="Velocity [m/s]")
            Colorbar(fig[2 ,1][1,2], plt.hm_P, label="open")
            #poly!(ax,points_XpYp)=#

            display(fig)

            #save("images\\$it_counter.png", fig)
        end
    end
    
end
