using CairoMakie

Pr = 10.0
b = 8.0/3.0
h = 0.005
tstart = 0.0
tend = 50.0
r = 28.0
n = trunc(Int, (tend - tstart) / h)

function f(y)
    ydot = zeros(3)
    ydot[1] = Pr * (y[2] - y[1])
    ydot[2] = y[1] * (Pr - y[3]) - y[2]
    ydot[3] = y[1] * y[2] - b * y[3]
    return ydot
end

function rk4(y)
    k1 = f(y)
    k2 = f(y + 0.5 * h * k1)
    k3 = f(y + 0.5 * h * k2)
    k4 = f(y + h * k3)
    y = y + h * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
    return y
end

x = zeros(n)
y = zeros(n)
z = zeros(n)
t = zeros(n)
x[1] = 0.0
y[1] = 0.5
z[1] = 0.5
t[1] = 0.0

for i in 1:n-1
    y_new = rk4([x[i], y[i], z[i]])
    t[i+1] = t[i] + h
    x[i+1] = y_new[1]
    y[i+1] = y_new[2]
    z[i+1] = y_new[3]
end

fig = Figure(resolution = (800, 600))
ax = Axis3(fig[1, 1]);
lines!(ax, x, y, z, color = :blue, linewidth = 0.5);
fig