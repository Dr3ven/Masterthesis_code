using BenchmarkTools
using LoopVectorization
using GLMakie

function matrix_mult(ρ, c_loc, β)
    @. c_loc = 1.0 / sqrt(ρ * β)
end

function looping(ρ, c_loc, β)
    @turbo for j in 1:size(ρ)[2], i in 1:size(ρ)[1]
        c_loc[i, j] = 1.0 / sqrt(ρ[i, j] * β[i, j])
    end
    #@show c_loc
end

function speedtest()

    nx = 301
    ny = 301

    rho_solid = 2800.0
    beta_solid = 1.0e-11

    ρ = ones(Float64, nx, ny) .* rho_solid
    β = ones(Float64, nx, ny) .* beta_solid
    c_loc_loop = zeros(Float64, nx, ny)
    c_loc_mult = zeros(Float64, nx, ny)
    err = zeros(Float64, nx, ny)

    @btime looping($ρ, $c_loc_loop, $β)
    @btime matrix_mult($ρ, $c_loc_mult, $β)

    err .= maximum(abs.(c_loc_loop .- c_loc_mult))
end

tauII, epsII = visco_elastic_wave_2D()
tauII2 = abs.(tauII) .* 1.0e-6
epsII2 = abs.(epsII)
fig = Figure()
ax = Axis(fig[1,1], xlabel="Strain", ylabel="Stress (MPa)")
scatter!(ax, epsII2, tauII2, markersize=10.0)
fig