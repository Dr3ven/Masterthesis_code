using CairoMakie

# Physics
Lx = 10.0                           # domain
ρ = 1.0                             # density
β = 1.0                             # compressibility
μ = 1.0                             # shear modulus
c = sqrt(1.0 / (β * ρ))             # speed of sound
P0 = 0.0                          # initial pressure at all points
Vx0 = 0.0                          # initial velocity in x-direction for all points

# Numerics
nx = 100                             # number of nodes in x
dx = Lx / nx                        # step size in x
nt = 500                             # number of time steps
dt = dx / (c * 4.0)                 # time step size

# Grid definition
xc = -(Lx - dx) / 2:dx:(Lx - dx) / 2        # grid nodes in x-direction
xv =   Lx       / 2:dx: Lx       / 2        # grid vertices in x-direction

# Allocations
P = ones(nx)
Vx = zeros(nx + 1)
τxx = zeros(nx)
dPdt = zeros(nx)
dVxdt = zeros(nx - 1)

# Initial conditions
P .= P0 .+ exp.(.- 1.0 .* xc.^2.0)        # initial pressure distribution
global t = 0.0                              # initial time

xc_vec = Vector(xc)

# Initial plotting
fig = Figure()
ax = Axis(fig[1,1], title="t = $t")
lines!(ax, xc_vec, P)
fig

for i = 1:nt
    global t += dt
    dPdt .= .-(1.0 ./ β) .* diff(Vx, dims=1) ./ dx
    P .= P .+ dPdt .* dt
    dVxdt .= .-(1.0 ./ ρ) .* diff(P, dims=1) ./ dx
    Vx[2:end-1] .= Vx[2:end-1] .+ dVxdt .* dt
    
    fig2 = Figure()
    ax2 = Axis(fig2[1,1], title="time = $t")#, limits=(nothing, nothing, -0.25, 0.25))
    ylims!(ax, -1.2, 1.2)
    lines!(ax2, xc_vec, P)
    display(fig2)
end
