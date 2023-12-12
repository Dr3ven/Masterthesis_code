using CairoMakie
# In this code we use visco-acoustic wave propagation to get the an incompressible Stokes solution
# Boundary conditions: left and right Vx=1; top and bottom Vy=0; with free slip everywhere
# There is fixed cylinder in the middle of the model (Vx=Vy=0)
# Intital conitions: Vx=1, Vy=0, everywhere. 
# Results: The initial pressure and velocity gradients generate wave propagation.
# Over time these waves dissipate, while they establish a stationary flow field, which is the solution of the Stokes problem. 
function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

# Physics
Lx = 2.0*π                                          # domain length
Ly = 2.0*π                                          # domain length
β  = 1.0                                            # compressibility [1/Pa]
ρ  = 1.0                                            # density [kg/m³]
η  = 0.25             	                            # shear modulus  [Pa]
c  = sqrt(1.0/β/ρ)                                  # p-wave velocity 
P0 = 0.0                                            # inital pressure                                          
Px = Lx/2.0                                         # x coordinate of the centre of the cylinder
Py = Ly/2.0                                         # y coordinate of the centre of the cylinder
R  = Lx/10.0                                        # radius of the cylinder

# Numerics
nx = 150                                            # number of grindpoints
ny = 150                                            # number of grindpoints
dx = Lx / nx                                        # step size
dy = Ly / ny                                        # step size
nt = 25e3                                           # number of time steps
dt = min(min(dx, dy) / (c)/4.5, min(dx^2.0, dy^2.0) / (4.0/3.0*η/ρ)/4.5)                              # Courant-Friedrichs-Lewy limit (CFL-condition)

# Allocation
P = zeros(Float64, nx, ny)                                       # Pressure
divV = zeros(Float64, nx, ny)                                  # velocity in x-direction
Vx = ones(Float64, nx + 1, ny)                                  # velocity in x-direction
Vy = zeros(Float64, nx, ny + 1)                                  # velocity in x-direction
τxx = zeros(Float64, nx, ny)                                     # stress in x-direction
τyy = zeros(Float64, nx, ny)                                     # stress in x-direction
τxy = zeros(Float64, nx + 1, ny + 1)                                     # stress in x-direction
dτxxdt = zeros(Float64, nx, ny)                                  # derivative of stress in on x-plane in x-direction w.r.t time
dτyydt = zeros(Float64, nx, ny)                                  # derivative of stress in on x-plane in x-direction w.r.t time
dτxydt = zeros(Float64, nx + 1, ny + 1)                                  # derivative of stress in on x-plane in x-direction w.r.t time
dVxdt = zeros(Float64, nx + 1, ny)                               # derivative of velocity w.r.t. time
dVydt = zeros(Float64, nx, ny + 1)                               # derivative of velocity w.r.t. time
dPdt = zeros(Float64, nx, ny)                                    # derivative of pressure with respect to time

#Vx_ana = zeros(Float64, nx + 1, ny)                              # velocity in x-direction

# Initilization
xc = dx/2.0:dx:Lx - (dx/2.0)                    # coordinates of cell centers in x (Array)
yc = dy/2.0:dy:Ly - (dy/2.0)                    # coordinates of cell centers in y (Array)
xv = LinRange(0.0, Lx, nx + 1)                  # coordinates of cell vertices
yv = LinRange(0.0, Ly, ny + 1)                  # coordinates of cell vertices
xc2d, yc2d = meshgrid(xc, yc)                   # coordinates of cell ceners in x and y (Matrices)
xVx, yVx = meshgrid(xv, yc)                     # coordinates of cell ceners in x and y (Matrices)
xVy, yVy = meshgrid(xc, yv)                     # coordinates of cell ceners in x and y (Matrices)
mask      = sqrt.((xc2d.-Px).^2.0 .+ (yc2d.-Py).^2.0) .> R
maskVx    = sqrt.((xVx .-Px).^2.0 .+ (yVx .-Py).^2.0) .> R
maskVy    = sqrt.((xVy .-Px).^2.0 .+ (yVy .-Py).^2.0) .> R
#Vx .= sin.(xv)                                  # velocity field (sinus function)
P .= 0.0.*exp.(.-1.0 .* (xc2d.^2.0 .+ yc2d.^2.0))  # pressure field (gaussian distribution)

# Initial plotting
fig = Figure()
ax = Axis(fig[1,1], xlabel ="x [-]", ylabel="y [-]")
hm = heatmap!(ax, xc2d, yc2d, P, colorrange=(0.0, (1.0/3.0)))
Colorbar(fig[1,2], hm, label="pressure [Pa]")
fig

# Solvers
global t = 0.0
for i = 1:nt
    global t += dt

    #Vx_ana .= (1.0 ./ 2.0) .* (sin.(xv .+ c .* t) .+ sin.(xv .- c .* t))

    # leap frog method
    divV            .= diff(Vx, dims=1) ./ dx .+ diff(Vy, dims=2) ./ dy
    dPdt            .= .-(1 ./ β) .* divV
    P              .+= dPdt .* dt
    τxx             .= (2.0 .* η) .* (diff(Vx, dims=1) ./ dx .- divV ./ 3.0)
    τyy             .= (2.0 .* η) .* (diff(Vy, dims=2) ./ dy .- divV ./ 3.0)
    τxy[2:end-1, 2:end-1] .= (2.0 .* η) .* (0.5 .* (diff(Vx[2:end-1, :], dims=2) ./ dy .+ diff(Vy[:, 2:end-1], dims=1) ./ dx))
    dVxdt[2:end-1, :]  .= (1 ./ ρ) .* (.-diff(P, dims=1) ./ dx .+ diff(τxx, dims=1) ./ dx .+ diff(τxy[2:end-1, :], dims=2) ./ dy)
    Vx             .+= dVxdt .* dt
    Vx              .= Vx.*maskVx
    dVydt[:, 2:end-1]  .= (1 ./ ρ) .* (.-diff(P, dims=2) ./ dy .+ diff(τyy, dims=2) ./ dy .+ diff(τxy[:, 2:end-1], dims=1) ./ dx)
    Vy             .+= dVydt .* dt
    Vy              .= Vy.*maskVy

    
    if i % 250 == 0
        fig2 = Figure()
        ax2 = Axis(fig2[1,1], xlabel ="x [-]", ylabel="y [-]")
        hm2 = heatmap!(ax2, xc2d, yc2d, P, label="numerical")
        Colorbar(fig2[1,2], hm2, label="pressure [Pa]")
        display(fig2)
    end
end
Vxc  = (Vx[1:end-1,:] + Vx[2:end,:]) ./ 2.0
Vyc  = (Vy[:,1:end-1] + Vy[:,2:end]) ./ 2.0
st   = 6
X    = xc2d[1:st:end,1:st:end]
Y    = yc2d[1:st:end,1:st:end]
VX   = Vxc[1:st:end,1:st:end]
VY   = Vyc[1:st:end,1:st:end]
Vabs = VX.^2.0 .+ VY.^2.0
fig3 = Figure()
ax     = Axis(fig3[1,1], xlabel ="x [-]", ylabel="y [-]")
hm = heatmap!(ax, xc2d, yc2d, 1.0.-mask, colormap=:grays)
arrows!(X[:], Y[:], VX[:], VY[:], arrowsize = 7, lengthscale = 0.06, linecolor=Vabs[:], arrowcolor=Vabs[:])
display(fig3)