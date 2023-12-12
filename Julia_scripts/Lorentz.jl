using CairoMakie

struct Parameters
    Pr
    r
    b
end

function dydt(Param, y)
    dydt = zeros(3)
    dWdt = Param.Pr .* (y[2] .- y[1])
    dT1dt = .-y[1] .* y[3] .+ Param.r .* y[1] .- y[2]
    dT2dt = .-y[1] .* y[2] .- Param.b .* y[3]
    dydt[1] = dWdt
    dydt[2] = dT1dt
    dydt[3] = dT2dt
    return dydt
end

function rkstep(h, y_1, P)
    k_1 = h .* dydt(P, y_1)
    k_2 = h .* dydt(P, y_1 + k_1./2.0)
    k_3 = h .* dydt(P, y_1 .+ k_2./2.0)
    k_4 = h .* dydt(P, y_1 .+ k_3)
    y_2 = y_1 .+ k_1./6.0 .+ k_2./3.0 .+ k_3./3.0 .+ k_4./6.0
    return y_2
end


Pr = 10.0
y = [10.0, 10.0, 10.0]
ysave = zeros(10000,3)
global time = 0.0
tstop = 50.0
h = 0.005
save_each = 1.0
global nstep = 0.0
global save_step = 0

r = 28.0
b = 8.0 / 3.0

P = Parameters(Pr, r, b)

global fig = Figure()
global ax = Axis3(fig[1,1], xlabel="W", ylabel="T1", zlabel="T2", title="Lorentz problem")

while time < tstop
    if mod(nstep, save_each) == 0
        global save_step += 1
        ysave[save_step, :] .= y
    end
    y .= y .+ rkstep(h, y, P)

    global nstep += 1.0
    global time += h
end

time_p = Vector(0.005:0.005:50)

scatter!(ax, ysave[:,1], ysave[:,2], ysave[:,3])
display(fig)
