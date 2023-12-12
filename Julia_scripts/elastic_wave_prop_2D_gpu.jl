using CairoMakie
using CUDA

# Macros
macro d_xa(A) esc(:($A[ix+1, iy] - $A[ix, iy])) end
macro d_ya(A) esc(:($A[ix, iy+1] - $A[ix, iy])) end
macro av_x(A)  esc(:(0.5 * ($A[ix, iy] + $A[ix+1, iy]))) end
macro av_y(A)  esc(:(0.5 * ($A[ix, iy] + $A[ix, iy+1]))) end
macro av_xy(A) esc(:(0.25 * ($A[ix+1, iy+1] + $A[ix, iy+1] + $A[ix+1, iy] + $A[ix, iy]))) end

# Solvers
function compute_P!(P, dPdt, dt)
    ix = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    iy = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    if  (ix<=size(P, 1) && iy<=size(P, 2))
        P[ix, iy] += dPdt[ix, iy] * dt
    end
    return
end

function compute_dPdt!(dPdt, β, Vx, Vy, dx, dy)
    ix = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    iy = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    if (ix<=size(dPdt, 1) && iy<=size(dPdt, 2))
        dPdt[ix, iy] = -(1.0 / β) * (@d_xa(Vx) ./ dx + @d_ya(Vy) / dy)
    end
    return
end

function compute_dτdt_diag!(dτxxdt, dτyydt, μ, λ, Vx, Vy, dx, dy)
    ix = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    iy = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    if (ix<=size(dτxxdt, 1) && iy<=size(dτxxdt, 2))
        dτxxdt[ix, iy] = (2.0 * μ + λ) * @d_xa(Vx) / dx  + λ * @d_ya(Vy) / dy
    end
    if (ix<=size(dτyydt, 1) && iy<=size(dτyydt, 2))
        dτyydt[ix, iy] = (2.0 * μ + λ) * @d_ya(Vy) / dy  + λ * @d_xa(Vx) / dx
    end
    return
end

function compute_dτdt_minor_diag!(dτxydt, μ, Vx, Vy, dx, dy)
    ix = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    iy = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    if (ix<=size(dτxydt, 1) - 2 && iy<=size(dτxydt, 2) - 2)
        dτxydt[ix + 1, iy + 1] = μ * (@d_xa(@av_xy(Vx)) / dy + @d_ya(@av_xy(Vy)) / dx)
    end
    return
end

function compute_τ_diag!(τxx, τyy, dτxxdt, dτyydt, dt)
    ix = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    iy = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    if (ix<=size(τxx, 1) && iy<=size(τxx, 2))
        τxx[ix, iy] += dτxxdt[ix, iy] * dt
    end
    if (ix<=size(τyy, 1) && iy<=size(τyy, 2))
        τyy[ix, iy] += dτyydt[ix, iy] * dt
    end
    return
end

function compute_τ_minor_diag!(τxy, dτxydt ,dt)
    ix = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    iy = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    if (ix<=size(τxy, 1) - 2 && iy<=size(τxy, 2) - 2)
        τxy[ix + 1, iy + 1] += dτxydt[ix + 1, iy + 1] * dt
    end
    return
end

function compute_dVdt!(dVxdt, dVydt, P, τxx, τyy, τxy, ρ, dx, dy)
    ix = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    iy = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    if (ix<=size(dVxdt, 1) - 2 && iy<=size(dVxdt, 2))
        dVxdt[ix + 1, iy] = -(1.0 / ρ) * (@d_xa(P) / dx + @d_xa(τxx) / dx + @d_ya(τxy[ix + 1, iy]) / dy)
    end
    if (ix<=size(dVydt, 1) && iy<=size(dVydt, 2) - 2)
        dVydt[ix, iy + 1] = -(1.0 / ρ) * (@d_ya(P) / dy + @d_ya(τyy) / dy + @d_xa(τxy[ix, iy + 1]) / dx)
    end
    return
end

function compute_V!(Vx, Vy, dVxdt, dVydt, dt)
    ix = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    iy = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    if (ix<=size(Vx, 1) - 2 && iy<=size(Vx, 2))
        Vx[ix + 1, iy] = Vx[ix + 1, iy] + dVxdt[ix + 1, iy] * dt 
    end
    if (ix<=size(Vy, 1) && iy<=size(Vy, 2) - 2)
        Vy[ix, iy + 1] = Vy[ix, iy + 1] + dVydt[ix, iy + 1] * dt
    end
    return
end

@views function elastic_wave_2D()
    # Select GPU
    CUDA.device!(0)

    # Physics
    Lx = 10.0                           # domain in x
    Ly = Lx                             # domain in y
    ρ = 2800.0                          # density
    β = 3.0e-6                          # compressibility
    μ = 1.0                             # shear modulus
    λ = 1.0                             # Lamé parameter
    #ν = λ / (2 * (λ + μ))              # Poisson ratio
    #Vp = sqrt((λ + 2 * μ) / ρ)
    c = sqrt(1.0 / (β * ρ))             # speed of sound
    P0 = 0.0                            # initial pressure at all points

    # Numerics
    nx = 255                            # number of nodes in x
    ny = 255                            # number of nodes in y
    nthread = (2, 2)                    # number of threads per block
    nblock  = cld.((nx, ny), nthread)   # number of blocks
    dx = Lx / nx                        # step size in x
    dy = Ly / ny                        # step size in y
    nt = 1000                           # number of time steps
    dt = min(dx, dy) / (c * sqrt(8.1))  # time step size                       

    # Grid definition
    xc = -(Lx - dx) / 2:dx:(Lx - dx) / 2        # grid nodes in x-direction
    yc = -(Ly - dy) / 2:dy:(Ly - dy) / 2        # grid nodes in x-direction
    xv =   Lx       / 2:dx: Lx       / 2        # grid vertices in x-direction
    yv =   Ly       / 2:dy: Ly       / 2        # grid vertices in x-direction

    # Allocations
    P = CUDA.zeros(Float64, nx, ny)
    Vx = CUDA.zeros(Float64, nx + 1, ny)
    Vy = CUDA.zeros(Float64, nx, ny + 1)
    τxx = CUDA.zeros(Float64, nx, ny)
    τyy = CUDA.zeros(Float64, nx, ny)
    τxy = CUDA.zeros(Float64, nx + 1, ny + 1)
    dPdt = CUDA.zeros(Float64, nx, ny)
    dVxdt = CUDA.zeros(Float64, nx + 1, ny)
    dVydt = CUDA.zeros(Float64, nx, ny + 1)
    dτxxdt = CUDA.zeros(Float64, nx, ny)
    dτyydt = CUDA.zeros(Float64, nx, ny)
    dτxydt = CUDA.zeros(Float64, nx + 1, ny + 1)
    
    # Initial conditions
    P .= P0 .+ exp.((.-5.0 .* xc.^2.0) .+ (.-5.0 .* yc'.^2.0))              # initial pressure distribution
    global t = 0.0                                                          # initial time

    xc_vec = Vector(xc)
    yc_vec = Vector(yc)

    P_cpu = Array(P)

    # Initial plotting
    fig = Figure()
    ax = Axis(fig[1,1], title="t = $t")#, limits=(nothing, nothing, nothing, 1.1))
    #lines!(ax, xc_vec, P[:, Int(nx/2)])
    hm = heatmap!(ax, xc_vec, yc_vec, P_cpu)
    Colorbar(fig[1,2], hm, label="Pressure")
    fig
    for i = 1:nt
        CUDA.@sync @cuda threads=nthread blocks=nblock compute_dPdt!(dPdt, β, Vx, Vy, dx, dy)
        CUDA.@sync @cuda threads=nthread blocks=nblock compute_P!(P, dPdt, dt)
        CUDA.@sync @cuda threads=nthread blocks=nblock compute_dτdt_diag!(dτxxdt, dτyydt, μ, λ, Vx, Vy, dx, dy)
        CUDA.@sync @cuda threads=nthread blocks=nblock compute_dτdt_minor_diag!(dτxydt, μ, bVx, bVy, dx, dy)
        CUDA.@sync @cuda threads=nthread blocks=nblock compute_τ_diag!(τxx, τyy, dτxxdt, dτyydt, dt)  
        CUDA.@sync @cuda threads=nthread blocks=nblock compute_τ_minor_diag!(τxy, dτxydt, dt)  
        CUDA.@sync @cuda threads=nthread blocks=nblock compute_dVdt!(dVxdt, dVydt, P, τxx, τyy, τxy, ρ, dx, dy)  
        CUDA.@sync @cuda threads=nthread blocks=nblock compute_V!(Vx, Vy, dVxdt, dVydt, dt)  
        
        # Boundary conditions
        Vx[:, 1] .= 0.0 #.-Vx[:, 1]
        Vx[:, end] .= 0.0 #.-Vx[:, end]

        Vy[1, :] .= 0.0 #.-Vy[1, :]
        Vy[end, :] .= 0.0 #.-Vy[end, :]
        
        global t += dt

        # Plotting

        P_cpu = Array(P)

        if i % 5.0 == 0.0
            fig2 = Figure()
            ax2 = Axis(fig2[1,1], title="time = $t")#, limits=(nothing, nothing, nothing, 1.1))
            #ylims!(ax, -1.2, 1.2)
            #lines!(ax2, xc_vec, P[:, Int(nx/2)])
            hm2 = heatmap!(ax2, xc_vec, yc_vec, P_cpu)#, colorrange=(0.0, 1.0))
            Colorbar(fig2[1,2], hm2, label="Pressure")
            display(fig2)
        end
    end
    return
end
