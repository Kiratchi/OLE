
using DifferentialEquations
using Plots
using LinearAlgebra
print("Hello")
 
r_1 = 2
r_2 = 2  
K_1 = 1
K_2 = 1

a_1 = 1/2
a_2 = 1/2


# Solution over time
f(dx,x,p,t) = [r_1*x[1]* (1-x[1]/K_1) - a_1*x[1]*x[2],
        r_2*x[2]* (1-x[2]/K_2) - a_2*x[2]*x[1] ]

tspan = 0,10
x0 = [2, 1]
prob = ODEProblem(f,x0,tspan)

sol = solve(prob)
plot(sol)


# Creating vector field
min_x = 0
max_x = 2
step_s = (max_x-min_x)/20

meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))
x1, x2 = meshgrid(min_x:step_s:max_x, min_x:step_s:max_x)

# Define elements of vector field
u1 = similar(x1)
u2 = similar(x2)
for i in eachindex(x1)
    dx = f(0, [x1[i], x2[i]], 0, 0)
    norm_dx = norm(dx) + 1e-8
    norm_dx = 20
    u1[i] = dx[1] / norm_dx
    u2[i] = dx[2] / norm_dx
end

# Plotting vector field
quiver(x1, x2, quiver=(u1, u2), aspect_ratio=1)

# CREATING NULLCLINES
x1_nullcline1 = min_x:0.5:max_x
x2_nullcline1 = @. r_1/a_1 * (1 - x1_nullcline1/K_1)
plot!(x1_nullcline1,x2_nullcline1, label="Nullcline x1")

x2_nullcline2 = min_x:0.5:max_x
x1_nullcline2 = @. r_2/a_2 * (1 - x2_nullcline2/K_2)
plot!(x1_nullcline2,x2_nullcline2, label="Nullcline x2")

plot!(xlims=(min_x, max_x), ylims=(min_x, max_x))