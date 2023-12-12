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


function visco_elastic_2D_coupled()
    # Physics
    Lx = 1.0e3                          # domain in x
    Ly = Lx                             # domain in y
    ρ0_air = 1.225                           # density
    ρ0_solid = 3.0e2                           # density
    dρ = 3.0e2 
    β = 1.0/141.0e3                          # compressibility
    β0_air = 1.0/141.0e3   # e-6                          # compressibility
    β0_solid = 1.0                          # compressibility
    η = 1.81e-5
    η0_air = 1.0   # 1.81e-5                    # viscosity air
    η0_solid = 1.0 # 2.5e-6                    # viscosity solid
    μ = 1.0                             # shear modulus
    μ0_air = 1.0                             # shear modulus
    μ0_solid = 1.0                             # shear modulus
    #c_air  = sqrt(1.0/β0/ρ0_air)          # speed of sound / p-wave velocity
    #c_solid = sqrt(1.0/β0/ρ0_solid)
    P0 = 1.0e5                            # initial pressure at all points
    g = 9.81                             # gravity

    # Numerics
    nx = 300                            # number of nodes in x
    ny = 300                            # number of nodes in y
    dx = Lx / nx                        # step size in x
    dy = Ly / ny                        # step size in y
    nt = 1000                           # number of time steps
    CFL_P = 1.0/16.0                           # Courant number for pressure
    CFL_V = 1.0/16.0                          # Courant number for velocity
    ξ = 0.0*0.95                          # relaxation factor for stress
    dt = 1.0 # min(min(min(dx, dy) / c_solid / 4.5, min(dx^2.0, dy^2.0) / ((4.0 / 3.0) * η0_solid / ρ0_solid) / 4.5), min(min(dx, dy) / c_air / 4.5, min(dx^2.0, dy^2.0) / ((4.0 / 3.0) * η0_air / ρ0_air) / 4.5))        # time step size                       

    # Grid definition
    xc = -(Lx - dx) / 2.0:dx:(Lx - dx) / 2.0        # grid nodes in x-direction
    yc = -(Ly - dy) / 2.0:dy:(Ly - dy) / 2.0        # grid nodes in x-direction
    xv =  -Lx       / 2.0:dx: Lx       / 2.0        # grid vertices in x-direction 
    yv =  -Ly       / 2.0:dy: Ly       / 2.0        # grid vertices in x-direction

    # Create copy of this code and try setting the dimensions to the same as in Danis code. Try to run atmosphere then and compare to dani

    # Allocations
    dtρ = zeros(Float64, nx, ny)
    dtPT = zeros(Float64, nx, ny)
    dtV = zeros(Float64, nx, ny)

    ρ = zeros(Float64, nx, ny)
    #ρ_ini = zeros(Float64, nx, ny)
    ρ_old = zeros(Float64, nx, ny)
    Fρx = zeros(Float64, nx + 1, ny)
    Fρy = zeros(Float64, nx, ny + 1)
    β = zeros(Float64, nx, ny)
    #μ = zeros(Float64, nx, ny)
    c = zeros(Float64, nx, ny)
    #η = zeros(Float64, nx, ny)
    P = zeros(Float64, nx, ny)
    divV = zeros(Float64, nx, ny)
    Mx = zeros(Float64, nx + 1, ny)
    My = zeros(Float64, nx, ny + 1)
    Mx_old = zeros(Float64, nx + 1, ny)
    My_old = zeros(Float64, nx, ny + 1)
    FMxx = zeros(Float64, nx, ny)
    FMyy = zeros(Float64, nx, ny)
    FMxy = zeros(Float64, nx - 1, ny + 1)
    FMyx = zeros(Float64, nx + 1, ny - 1)
    Vx = zeros(Float64, nx + 1, ny)
    Vy = zeros(Float64, nx, ny + 1)
    εxx = zeros(Float64, nx, ny)
    εyy = zeros(Float64, nx, ny)
    εxy = zeros(Float64, nx + 1, ny + 1)
    τxx = zeros(Float64, nx, ny)
    τyy = zeros(Float64, nx, ny)
    τxy = zeros(Float64, nx + 1, ny + 1)
    σxx = zeros(Float64, nx, ny)
    σyy = zeros(Float64, nx, ny)
    σxy = zeros(Float64, nx + 1, ny + 1)
    dρdt = zeros(Float64, nx, ny)
    dPdt = zeros(Float64, nx, ny)
    dMxdt = zeros(Float64, nx + 1, ny)
    dMydt = zeros(Float64, nx, ny + 1)
    dMxdξ = zeros(Float64, nx - 1, ny)
    dMydξ = zeros(Float64, nx, ny - 1)
    dVxdt = zeros(Float64, nx + 1, ny)
    dVydt = zeros(Float64, nx, ny + 1)
    dτxxdt = zeros(Float64, nx, ny)
    dτyydt = zeros(Float64, nx, ny)
    dτxydt = zeros(Float64, nx + 1, ny + 1)

    ρRes = zeros(Float64, nx, ny)
    MxRes = zeros(Float64, nx - 1, ny)
    MyRes = zeros(Float64, nx, ny - 1)

    maskrho_air = zeros(Float64, nx, ny)
    maskrho_solid = zeros(Float64, nx, ny)
    maskVx_air = zeros(Float64, nx, ny)
    maskVx_solid = zeros(Float64, nx, ny)
    maskVy_air = zeros(Float64, nx, ny)
    maskVy_solid = zeros(Float64, nx, ny)
    radrho = zeros(Float64, nx, ny)
    radVx = zeros(Float64, nx + 1, ny)
    radVy = zeros(Float64, nx, ny + 1)

    ### Geometrie part ---------------------------
    x2dc, y2dc = meshgrid(xc,yc)                                        # 2d mesh of x- and y-coordinates of density nodes
    x2dVx, y2dVx = meshgrid(xv,yc)                                      # 2d mesh of x- and y-coordinates of velocity nodes in x-direction
    x2dVy, y2dVy = meshgrid(xc,yv)                                      # 2d mesh of x- and y-coordinates of velocity nodes in y-direction

    locX      = 0.0                                                     # x-coordinate of circle
    locY      = -0.35*Ly                                                # y-coordinate of circle
    diam      = 0.2*Lx                                                  # diameter of circle
    radrho    .= sqrt.((x2dc .-locX).^2.0 .+ (y2dc .-locY).^2.0)        # radius of circle
    radVx     .= sqrt.((x2dVx.-locX).^2.0 .+ (y2dVx.-locY).^2.0)        # radius of circle
    radVy     .= sqrt.((x2dVy.-locX).^2.0 .+ (y2dVy.-locY).^2.0)        # radius of circle
    Xp        = [-1.0/2.0, -1.0/2.0, -0.3, -1.0/8.0, -0.01, -0.01, 0.01, 0.01, 1.0/8.0, 0.3, 1.0/2.0, 1.0/2.0].*Lx          # x-coordinates of polygon
    Yp        = [-1.0/2.0, -1.0/5.0, -1.0/5.0, -0.1, -0.15,  -0.28, -0.28, -0.15,  -0.1, -1.0/5.0, -1.0/5.0, -1.0/2.0].*Ly  # y-coordinates of polygon
    Xp2       = [-0.02, -0.02,  0.02,  0.02].*Lx                                                                            # x-coordinates of polygon
    Yp2       = [-0.15,  -0.35, -0.35, -0.15].*Ly                                                                             # y-coordinates of polygon
    Xp3        = [-1.0/2.0, -0.3, -1.0/8.0, -0.01, 0.01, 1.0/8.0, 0.3, 1.0/2.0].*Lx          # x-coordinates of polygon
    Yp3        = [-1.0/5.0, -1.0/5.0, -0.1, -0.15, -0.15,  -0.1, -1.0/5.0, -1.0/5.0].*Ly  # y-coordinates of polygon

    inpolyrho   = inpolygon(x2dc,y2dc,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon
    inpolyVx    = inpolygon(x2dVx,y2dVx,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon
    inpolyVy    = inpolygon(x2dVy,y2dVy,Xp,Yp)                                 # 1: inside polygon, 0: outside polygon

    maskrho_air = 1.0 .- inpolyrho                                             # mask for density
    #maskrho_air[findall(<(diam ./ 2.0), radrho)] .= 1.0                               # mask for density in circle
    maskrho_air[radrho .< diam ./ 2.0] .= 1.0                               # mask for density in circle
    maskVx_air  = 1.0 .- inpolyVx                               # mask for velocity in x-direction
    #maskVx_air[findall(<(diam ./ 2.0), radVx)] .= 1.0                                 # mask for velocity in x-direction in circle
    maskVx_air[radVx .< diam ./ 2.0] .= 1.0                                 # mask for velocity in x-direction in circle
    maskVy_air  = 1.0 .- inpolyVy                               # mask for velocity in y-direction
    #maskVy_air[findall(<(diam ./ 2.0), radVy)] .= 1.0                                 # mask for velocity in y-direction in circle
    maskVy_air[radVy .< diam ./ 2.0] .= 1.0                                 # mask for velocity in y-direction in circle

    maskrho_solid = 1.0 .- maskrho_air
    maskVx_solid  = 1.0 .- maskVx_air
    maskVy_solid  = 1.0 .- maskVy_air

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



    # Initial conditions--------------------------
    #ρ_ini .+= ρ0_air .* (maskrho_air .== 1.0) .+ ρ0_solid .* (maskrho_solid .== 1.0)

    #β .= β0_air .* (maskrho_air .== 1.0) .+ β0_solid .* (maskrho_solid .== 1.0)
    #P .= (P0 .+ exp.((.-5.0e-5 .* (xc).^2.0)./ (2.0 .* 0.5.^2.0) .+ (.-5.0e-5 .* (yc' .+ Ly .* (1.0/2.825)).^2.0) ./ (2.0 .* 0.5.^2.0)))        # initial pressure distribution
    #ρ .+= ρ0_air .* (maskrho_air .== 1.0) .+ ρ0_solid .* (maskrho_solid .== 1.0) 
    
    #---- TEST 
    P         .= P0.*exp.(-g.*(y2dc.+1.0./5.0.*Ly).*0.028./288.0./8.314) # barometric setting atmosphere: P = P0*exp(-(g*(h-h0)*M)/(T*R)) M: Mass density of air, T: Temperature, R: Gas constant, h: height, P0: Pressure at sea level
    ρ       .= (ρ0_air .* P) ./ P0 #ρ0_air ./ P0 .* P             # equation of state for density depending on pressure
    ρ[radrho .< diam ./ 2.0] .= ρ0_air .+ dρ        # initial density in the circle

    ρ[inpolygon(x2dc,y2dc,Xp2,Yp2) .== 1.0 .&& y2dc.< locY.+diam./2.0] .= ρ0_air.+dρ          # initial density of the conduit overlapping with the circular chamber (so in the chamber)
    ρ[inpolygon(x2dc,y2dc,Xp2,Yp2) .== 1.0 .&& y2dc.>=locY.+diam./2.0] .= .-0.0 .* (y2dc[inpolygon(x2dc,y2dc,Xp2,Yp2).==1.0 .&& y2dc.>=locY.+diam./2.0].+0.15.*Ly).*dρ./(-0.15.-(locY.+diam./2.0)).+ρ0_air    # initial density in the conduit
    P         .= (P0 .* ρ) ./ ρ0_air
    #----
    
    # η .+= η0_air .* (maskrho_air .== 1.0) .+ η0_solid .* (maskrho_solid .== 1.0)    # initial viscosity distribution
    #μ .+= μ0_air  .* (maskrho_air .== 1.0) .+ μ0_solid  .* (maskrho_solid .== 1.0)  # initial shear modulus distribution
    #c .+= c_air  .* (maskrho_air .== 1.0) .+ c_solid  .* (maskrho_solid .== 1.0)    # initial sound speed distribution
    t = 0.0
    err = 1.0    
    #---------------------------------------------

    #min_P = minimum(P[P .> 0.0])
    #max_P = maximum(P[P .> 0.0])

    # Initial plotting----------------------------
    x_circ = zeros(length(x_circ_ind))
    y_circ = zeros(length(y_circ_ind))

    x_circ = xc[x_circ_ind]
    y_circ = yc[y_circ_ind]
    U = av_x(Vx)
    V = av_y(Vy)
    data_plt = sqrt.(U.^2.0 .+ V.^2.0)

    #xyticks = ["-400", "-200", "0", "200", "400"]
    #xytick_positions = range(start=-400, stop=400, step=200)

    fig = Figure()
    ax = Axis(fig[1,1][1,1], title="t = $t")#, xticks=(xytick_positions, xyticks), yticks=(xytick_positions, xyticks),
                #yticklabelsize=25, xticklabelsize=25, xlabelsize=25, ylabelsize=25)#, aspect = DataAspect())#, limits=(nothing, nothing, nothing, 1.1))
    #lines!(ax, xc_vec, P[:, Int(nx/2)])
    hm = heatmap!(ax, x2dc, y2dc, P, colormap=Reverse(:roma))#, colorrange=(min_P, max_P))
    #hm = heatmap!(ax, x2dc, y2dc, data_plt, colormap=Reverse(:roma), colorrange=(0.0, 1.0))
    #lines!(ax, Xp, Yp, color=:white)
    scatter!(ax, x_circ, y_circ, color=:white, markersize=4.0)
    Colorbar(fig[1,2][1,1], hm, label="Pressure")#, labelsize=25, ticklabelsize=25)#, vertical=false)
    #Colorbar(fig[1,2][1,1], hm, label="Velocity", labelsize=25, ticklabelsize=25)#, vertical=false)
    #ax3 = Axis3(fig[1,2][1,1], title="time = $t")
    #limits!(ax3, -Lx / 2.0, Lx / 2.0, -Ly / 2.0, Ly / 2.0, -0.1, 1.0)
    #sur2 = surface!(ax3, xc_vec, yc_vec, P)
    #display(fig)
    save("./Plots/Test/0.png", fig)
    #------------------------------------------------

    nan = false

    reset_timer!()
    for i = 1:2#nt
        ρ_old .= ρ                                              # save old density values for time integration
        Mx_old .= Mx                                                # save old momentum values in x-direction for time integration
        My_old .= My                                                # save old momentum values in y-direction for time integration
        err = 1.0
        err2 = zeros(Float64, 3)
        iter = 0
        it_counter = 0
        dMxdξ .= 0.0                                              # stress derivative of momentum in x-direction
        dMydξ .= 0.0                                              # stress derivative of momentum in y-direction
        while err > 1.0e-3
            iter += 1
            # compressibility, local sound of speed and time steps
            @timeit "β" β .= 1.0 ./ P
            @timeit "dt" dt = minimum([dx ./ maximum(abs.(Vx)), dy ./ maximum(abs.(Vy)), min(dx, dy) .* sqrt(maximum(ρ .* β)), min(dx, dy).^2.0 ./ maximum(η)]) .* 4.0 # time step size
            @timeit "c" c .= 1.0 ./ sqrt.(maximum(ρ .* β))                                                    # local speed of sound
            @timeit "dtPT" dtPT .= min.(min.(min.(dx ./ abs.(av_x(Vx)), dx ./ abs.(av_y(Vy))), min(dx, dy).^2.0 ./ μ .* ones(nx, ny)), dx ./ c) # time step size for pressure and temperature
            @timeit "dtρ" dtρ .= 1.0 ./ (1.0 ./ dt .+ 1.0 ./ (min(dx, dy) ./ c ./ 4.1)) 

            # conservation of mass
            @timeit "Fρx" Fρx[2:end-1, :] .= (Vx[2:end-1, :] .> 0.0) .* ρ[1:end-1, :] .* Vx[2:end-1, :] .+ (Vx[2:end-1, :] .< 0.0) .* ρ[2:end, :] .* Vx[2:end-1, :]
            @timeit "Fρy" Fρy[:, 2:end-1] .= (Vy[:, 2:end-1] .> 0.0) .* ρ[:, 1:end-1] .* Vy[: ,2:end-1] .+ (Vy[:, 2:end-1] .< 0.0) .* ρ[:, 2:end] .* Vy[: ,2:end-1]
            @timeit "dρdt" dρdt .= (ρ .- ρ_old) ./ dt
            @timeit "ρRes" ρRes .= ρRes .* maskrho_solid
            @timeit "ρRes" ρRes .= .- dρdt .- diff(Fρx, dims=1) ./ dx .- diff(Fρy, dims=2) ./ dy
            @timeit "ρ" ρ .= ρ .+ ρRes .* dtρ .* CFL_P

            # divergence of velocity, pressure, stresses, strain rates and additional time step
            @timeit "divV" divV .= diff(Vx, dims=1) ./ dx .+ diff(Vy, dims=2) ./ dy
            @timeit "P" P .= P0 ./ ρ0_air .* ρ #.-(1.0 ./ β) .* (ρ_ini ./ ρ).* divV
            #@timeit "dPdt" dPdt .= P0 ./ ρ0_air .* ρ #.-(1.0 ./ β) .* (ρ_ini ./ ρ).* divV
            #@timeit "P" P .= P .+ dPdt .* dt
            @timeit "εxx" εxx .= diff(Vx, dims=1) ./ dx .- (1.0 ./ 3.0) .* divV
            @timeit "εyy" εyy .= diff(Vy, dims=2) ./ dy .- (1.0 ./ 3.0) .* divV
            @timeit "εxy" εxy[2:end-1, 2:end-1] .= 0.5 .* (diff(Vx[2:end-1, :], dims=2) ./ dy .+ diff(Vy[:, 2:end-1], dims=1) ./ dx)
            @timeit "dτxxdt" dτxxdt .= 2.0 .* μ .* εxx #.- (μ ./ η) .* τxx
            @timeit "dτyydt" dτyydt .= 2.0 .* μ .* εyy #.- (μ ./ η) .* τyy
            @timeit "dτxydt" dτxydt[2:end-1, 2:end-1] .= 2.0 .* μ .* εxy[2:end-1, 2:end-1] #.- (av_xy(μ) ./ η) .* τxy[2:end-1, 2:end-1]
            #@timeit "τxx" τxx .= τxx .+ dτxxdt .* dt
            #@timeit "τyy" τyy .= τyy .+ dτyydt .* dt
            #@timeit "τxy" τxy .= τxy .+ dτxydt .* dt
            @timeit "σxx" σxx .= dτxxdt .- P
            @timeit "σyy" σyy .= dτyydt .- P
            #@timeit "σxy" σxy .= τxy .- diff(P, dims=1) ./ dx
            @timeit "dtV" dtV .= 1.0 ./ (1.0 ./ dt .+ 1.0 ./ (min(dx,dy).^2.0 ./ μ ./ 4.0)) .* CFL_V                                           # time step size for velocity

            # conservation of linear momentum
            @timeit "Mx" Mx[2:end-1, :] .= av_x(ρ) .* Vx[2:end-1, :]
            @timeit "My" My[:, 2:end-1] .= av_y(ρ) .* Vy[:, 2:end-1]
            @timeit "FMxx" FMxx .= (av_x(Vx) .> 0.0) .* av_x(Vx) .* Mx[1:end-1, :] .+ (av_x(Vx) .< 0.0) .* av_x(Vx) .* Mx[2:end, :]
            @timeit "FMxy" FMxy[:, 2:end-1] .= (av_x(Vy[:, 2:end-1]) .> 0.0) .* av_x(Vy[:, 2:end-1]) .* Mx[2:end-1, 1:end-1] .+ (av_x(Vy[:, 2:end-1]) .> 0.0) .* av_x(Vy[:, 2:end-1]) .* Mx[2:end-1, 2:end]
            @timeit "FMyy" FMyy .= (av_y(Vy) .> 0.0) .* av_y(Vy) .* My[:, 1:end-1] .+ (av_y(Vy) .< 0.0) .* av_y(Vy) .* My[:, 2:end]
            @timeit "FMyx" FMyx[2:end-1, :] .= (av_y(Vx[2:end-1, :]) .> 0.0) .* av_y(Vx[2:end-1, :]) .* My[1:end-1, 2:end-1] .+ (av_y(Vx[2:end-1, :]) .> 0.0) .* av_y(Vx[2:end-1, :]) .* My[2:end, 2:end-1]
            @timeit "dMxdt" dMxdt .= (Mx .- Mx_old) ./ dt
            @timeit "dMydt" dMydt .= (My .- My_old) ./ dt
            @timeit "MxRes" MxRes .= MxRes .* maskVx_solid[2:end-1, :]
            @timeit "MxRes" MxRes .= .-dMxdt[2:end-1, :] .- diff(FMxx .- σxx, dims=1) ./ dx .- diff(FMxy .- dτxydt[2:end-1, :], dims=2) ./ dy
            @timeit "MyRes" MyRes .= MyRes .* maskVy_solid[:, 2:end-1]
            @timeit "MyRes" MyRes .= .-dMydt[:, 2:end-1] .- diff(FMyy .- σyy, dims=2) ./ dy .- diff(FMyx .- dτxydt[:, 2:end-1], dims=1) ./ dx
            @timeit "dMxdξ" dMxdξ .= MxRes .+ dMxdξ .* ξ                                                                                    # stress derivative of momentum in x-direction
            @timeit "dMydξ" dMydξ .= MyRes .+ dMydξ .* ξ                                                                                    # stress derivative of momentum in y-direction  
            @timeit "update Mx" Mx[2:end-1,:]  .= Mx[2:end-1,:] .+ dMxdξ.*av_x(dtPT).*CFL_V                                                     # updating momentum in x-direction
            @timeit "update My" My[:,2:end-1]  .= My[:,2:end-1] .+ dMydξ.*av_y(dtPT).*CFL_V                                                     # updating momentum in x-direction          
            @timeit "Vx" Vx[2:end-1, :] .= Mx[2:end-1, :] ./ av_x(ρ)
            @timeit "Vy" Vy[:, 2:end-1] .= My[:, 2:end-1] ./ av_y(ρ)
            
            #=
            if maximum(isnan.(P))
                print("P has NaN\n")
                nan = true
                #@show divV
                break
            elseif maximum(isnan.(β))
                print("Compressibility has NaN\n")
                nan = true
                break
            elseif maximum(isnan.(Vx)) || maximum(isnan.(Vy))
                print("Velocity has NaN")
                nan = true
                break
            elseif maximum(isnan.(τxx)) || maximum(isnan.(τyy)) || maximum(isnan.(τxy))
                print("Stress has NaN")
                nan = true
                break
            elseif maximum(isnan.(ρRes)) || maximum(isnan.(MxRes)) || maximum(isnan.(MyRes))
                print("Residual derivations have NaN")
                nan = true
                break
            elseif maximum(isnan.(Mx)) || maximum(isnan.(My))
                print("Momentum have NaN")
                nan = true
                break
            end=#

            if mod(iter, 1) == 0
                it_counter += 1
                @timeit "err" err = maximum(abs.([ρRes[:]; MxRes[:]; MyRes[:]]))    
                err2 .= [abs.(maximum(ρRes[:])); abs.(maximum(MxRes[:])); abs.(maximum(MyRes[:]))]                                                                  # error for time integration concerning density and momentum
                print("it_counter = $it_counter, ρRes = $(err2[1]), MxRes = $(err2[2]), MyRes = $(err2[3])\n")
                if maximum(isnan.(err)) == true 
                    break
                end
            end

            t += dt

            if i % 2 == 0
                X = av_xy(x2dc)
                Y = av_xy(y2dc)
                U = av_x(Vx)
                V = av_y(Vy)

                data_plt = sqrt.(U.^2.0 .+ V.^2.0)

                fig2 = Figure()
                ax2 = Axis(fig2[1,1][1,1], title="t = $t")#, xticks=(xytick_positions, xyticks), yticks=(xytick_positions, xyticks),
                            #yticklabelsize=25, xticklabelsize=25, xlabelsize=25, ylabelsize=25, title="time = $t")#, aspect = DataAspect())#, limits=(nothing, nothing, nothing, 1.1))
                #lines!(ax2, xc_vec, P[:, Int(nx/2)])
                #hm2 = heatmap!(ax2, x2dc, y2dc, P, colormap=Reverse(:roma), colorrange=(P0, maximum(P)))
                hm2 = heatmap!(ax2, x2dc, y2dc, data_plt, colormap=Reverse(:roma))#, colorrange=(min_P, max_P))#, colorrange=(0.0, 1.0))
                scatter!(ax2, x_circ, y_circ, color=:white, markersize=4.0)
                #Colorbar(fig2[1,2][1,1], hm2, label="Pressure", labelsize=25, ticklabelsize=25)#, vertical=false)
                Colorbar(fig2[1,2][1,1], hm2, label="Velocity")#, labelsize=25, ticklabelsize=25)#, vertical=false)
                #ax3 = Axis3(fig2[1,2][1,1], title="time = $t")
                #limits!(ax3, -Lx / 2.0, Lx / 2.0, -Ly / 2.0, Ly / 2.0, -0.1, 1.0)
                #sur2 = surface!(ax3, xc_vec, yc_vec, P)
                #display(fig2)
                save("./Plots/Test/$(i).png", fig2)
            end
        end

        if nan
            break
        end
    end
    print_timer()
    return MxRes, dMxdt, FMxx, FMxy, σxx, dτxydt
end