
using DifferentialEquations
using Plots
print("Hello")
 
r_1 = 1.2
r_2 = 1.1
K_1 = 1
K_2 = 1
a_1 = 1
a_2 = 2

f(x,p,t) = [r_1*x[1]* (1-x[1]/K_1) - a_1*x[1]*x[2],
            r_2*x[2]* (1-x[2]/K_1) - a_1*x[2]*x[1] ]
tspan = 0,100
x0 = [2, 1]
prob = ODEProblem(f,x0,tspan)

sol = solve(prob)
plot(sol)