
using DifferentialEquations
using Plots
using LinearAlgebra
using LaTeXStrings

cd("C:/Users/jonat/Documents/Julia/OLE")
pyplot()  # switch to PyPlot backend for fixing plot bug



# Time solution for logistic equation
function plot_log_eq(K,r, x0_vec, filename)
    tspan = (0.0,10.0)
    plot() #Initialize plot

    # Solving ODE for all initial conditions in x0_vec
    for x0 in x0_vec
        f(x, p, t) = r * x * (1 - x / K)
        prob = ODEProblem(f,x0,tspan)
        sol = solve(prob)
        plot!(sol, lw=3)
    end

    plot!(legend = false,
          xlims=tspan,
          ylims= (0,2.1),
          xlabel=L"t",
          ylabel=L"x(t)")
    display(current())      
    savefig(filename)
    println("Figure saved as $filename")
end

plot_log_eq(1.3,1, [1,0.1,0.01,2], "log_eq.png")



# Phase portrait for competative Lotka-Volterra equation
function create_phase_portrait(r_1, r_2, K_1, K_2, a_1, a_2, filename)
    # Solution over time
    f(dx, x, p, t) = [r_1 * x[1] * (1 - x[1] / K_1) - a_1 * x[1] * x[2],
        r_2 * x[2] * (1 - x[2] / K_2) - a_2 * x[2] * x[1]]


    # Creating vector field
    min_x = 0
    max_x = 2
    step_s = (max_x - min_x) / 20

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

    # Creating Nullcline
    x1_nullcline1 = min_x:0.5:max_x
    x2_nullcline1 = @. r_1 / a_1 * (1 - x1_nullcline1 / K_1)
    plot!(x1_nullcline1, x2_nullcline1, label=L"x_1'=0", lw=3)

    x2_nullcline2 = min_x:0.5:max_x
    x1_nullcline2 = @. r_2 / a_2 * (1 - x2_nullcline2 / K_2)
    plot!(x1_nullcline2, x2_nullcline2, label=L"x_2'=0", lw=3)

    plot!(legend = :topright,
          xlims=(min_x, max_x),
          ylims=(min_x, max_x),
          xlabel=L"x_1",
          ylabel=L"x_2")
    
    display(current())
    savefig(filename)
    println("Figure saved as $filename")
end


# Phase portrait for competative lotka-volterra equation for all 9 possible situations
create_phase_portrait(2,2,1,1,0.5,0.5, "phase_portrait_inStable.png")
create_phase_portrait(2,1,2,1,1,1, "phase_portrait_x1Win.png")
create_phase_portrait(1,1,1.7,1.7,1,1, "phase_portrait_nonStable.png")
create_phase_portrait(1,2,1,2,1,1, "phase_portrait_x2Win.png")
create_phase_portrait(1.5,1.5,1.5,1.5,1,1, "phase_portrait_symetric.png")
create_phase_portrait(0.5,1,1,1.5,1,1, "phase_portrait_x2Win_2.png")
create_phase_portrait(1, 0.5,1.5,1,1,1, "phase_portrait_x1Win_2.png")
create_phase_portrait(1.5,1,1,0.5,1,1, "phase_portrait_x1DWin.png")
create_phase_portrait(1, 1.5,0.5,1,1,1, "phase_portrait_x2DWin.png")


# Time solution for competative lotka-volterra equation
function plot_lotkaV(r_1, r_2, K_1, K_2, a_1, a_2, x0_vec, filename)
    tspan = (0.0,10.0)
    plot() #Initialize plot
    color_palette = palette(:tab10)

    # Solving ODE for all initial conditions in x0_vec
    for (i,x0) in enumerate(x0_vec)
        f!(dx, x, p, t) = begin
            dx[1] = r_1 * x[1] * (1 - x[1] / K_1) - a_1 * x[1] * x[2]
            dx[2] = r_2 * x[2] * (1 - x[2] / K_2) - a_2 * x[2] * x[1]
        end
        prob = ODEProblem(f!,x0,tspan)
        sol = solve(prob)
        plot!(sol, idxs=1, color=color_palette[i],label="", linestyle=:dash, lw=3)
        plot!(sol, idxs=2, color=color_palette[i],label="", linestyle=:solid, lw=3)
    end

    # Dummy plots to fix legend
    plot!([0,1], [4,5], color=:black, label=L"x_1", linestyle=:dash, lw=3)
    plot!([0,1], [4,5], color=:black, label=L"x_2", linestyle=:solid, lw=3)


    plot!(legend = :topright,
          xlims=tspan,
          ylims= (0,2.1),
          xlabel=L"t",
          ylabel=L"x(t)")
    display(current())
    savefig(filename)
    println("Figure saved as $filename")
end


# Initial condition
x0_vec = [[0.1, 0.1], [1.5, 0.2], [0.2, 1.5]]

# Time solution for competative lotka-volterra equation for all 9 possible situations
plot_lotkaV(2,2,1,1,0.5,0.5,x0_vec, "lotkaV_inStable.png")
plot_lotkaV(2,1,2,1,1,1,x0_vec, "lotkaV_x1Win.png")
plot_lotkaV(1,1,1.7,1.7,1,1,x0_vec, "lotkaV_nonStable.png")
plot_lotkaV(1,2,1,2,1,1,x0_vec, "lotkaV_x2Win.png")
plot_lotkaV(1.5,1.5,1.5,1.5,1,1,x0_vec, "lotkaV_symetric.png")
plot_lotkaV(0.5,1,1,1.5,1,1,x0_vec, "lotkaV_x2Win_2.png")
plot_lotkaV(1,0.5,1.5,1,1,1,x0_vec, "lotkaV_x1Win_2.png")
plot_lotkaV(1.5,1,1,0.5,1,1,x0_vec, "lotkaV_x1DWin.png")
plot_lotkaV(1,1.5,0.5,1,1,1,x0_vec, "lotkaV_x2DWin.png")



# Extra plot for non-stable equilibrium for more initial conditions
x0_vec = [[0.1, 0.1],[0.1, 0.105],[0.1, 0.11],[0.1, 0.15],]
plot_lotkaV(1,1,1.7,1.7,1,1,x0_vec, "lotkaV_nonStable_2.png")


# Extra plot for symmetric equilibrium with pertubations to parameter
create_phase_portrait(1.5,1.5,1.5,1.5,0.95,1, "phase_portrait_symetric_2.png")
x0_vec = [[0.1, 0.1],[0.1, 0.105],[0.1, 0.11],[0.1, 0.15],]
plot_lotkaV(1.5,1.5,1.5,1.5,0.95,1,x0_vec, "lotkaV_symetric_2.png")