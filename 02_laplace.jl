using Printf
using Plots

function jacobi_method(U_start, Nx, Ny, coef_x, coef_y, denom, eps, max_iter)
    U_old = copy(U_start)
    U_new = copy(U_start)

    iter = 0
    error = Inf

    while error > eps && iter < max_iter
        error = 0.0

        for i in 2:Nx
            for j in 2:Ny
                U_new[i, j] = (
                    coef_x * (U_old[i-1, j] + U_old[i+1, j]) +
                    coef_y * (U_old[i, j-1] + U_old[i, j+1])
                ) / denom

                error = max(error, abs(U_new[i, j] - U_old[i, j]))
            end
        end

        U_old, U_new = U_new, U_old
        iter += 1
    end

    return U_old, iter, error
end

function gauss_seidel_method(U_start, Nx, Ny, coef_x, coef_y, denom, eps, max_iter)
    U = copy(U_start)

    iter = 0
    error = Inf

    while error > eps && iter < max_iter
        error = 0.0

        for i in 2:Nx
            for j in 2:Ny
                old = U[i, j]

                U[i, j] = (
                    coef_x * (U[i-1, j] + U[i+1, j]) +
                    coef_y * (U[i, j-1] + U[i, j+1])
                ) / denom

                error = max(error, abs(U[i, j] - old))
            end
        end

        iter += 1
    end

    return U, iter, error
end

function max_error(U, x, y, Nx, Ny, exact_solution)
    err = 0.0

    for i in 1:Nx+1
        for j in 1:Ny+1
            err = max(err, abs(U[i, j] - exact_solution(x[i], y[j])))
        end
    end

    return err
end

function build_initial_grid(x, y, Nx, Ny, bottom, top, left, right)
    U = zeros(Nx + 1, Ny + 1)

    for i in 1:Nx+1
        U[i, 1] = bottom(x[i])
        U[i, Ny+1] = top(x[i])
    end

    for j in 1:Ny+1
        U[1, j] = left(y[j])
        U[Nx+1, j] = right(y[j])
    end

    return U
end

println("Введіть кількість поділів по x Nx:")
Nx = parse(Int, readline())

println("Введіть кількість поділів по y Ny:")
Ny = parse(Int, readline())

println("Введіть точність eps:")
eps = parse(Float64, readline())

println("Введіть максимальну кількість ітерацій max_iter:")
max_iter = parse(Int, readline())

if Nx < 2 || Ny < 2
    error("Nx і Ny мають бути більші або рівні 2.")
end

if eps <= 0
    error("eps має бути додатним числом.")
end

ax = 0.0
bx = 1.0
ay = 0.0
by = 1.0

hx = (bx - ax) / Nx
hy = (by - ay) / Ny

x = collect(range(ax, bx, length = Nx + 1))
y = collect(range(ay, by, length = Ny + 1))

bottom(x) = 0.0
top(x) = sin(pi * x)
left(y) = 0.0
right(y) = 0.0

exact_solution(x, y) = sinh(pi * y) / sinh(pi) * sin(pi * x)

U = build_initial_grid(x, y, Nx, Ny, bottom, top, left, right)

coef_x = 1.0 / hx^2
coef_y = 1.0 / hy^2
denom = 2.0 * (coef_x + coef_y)

println()
println("Параметри задачі:")
println("Область: [0,1] × [0,1]")
println("Nx = ", Nx, ", Ny = ", Ny)
println("hx = ", hx, ", hy = ", hy)
println("eps = ", eps)
println("max_iter = ", max_iter)

println()
println("Граничні умови:")
println("u(x,0) = 0")
println("u(x,1) = sin(πx)")
println("u(0,y) = 0")
println("u(1,y) = 0")

println()
println("Точний розв'язок:")
println("u(x,y) = sinh(πy)/sinh(π) * sin(πx)")

println()
U_jacobi, iter_jacobi, conv_jacobi = jacobi_method(
    U, Nx, Ny, coef_x, coef_y, denom, eps, max_iter
)

U_gs, iter_gs, conv_gs = gauss_seidel_method(
    U, Nx, Ny, coef_x, coef_y, denom, eps, max_iter
)

error_jacobi = max_error(U_jacobi, x, y, Nx, Ny, exact_solution)
error_gs = max_error(U_gs, x, y, Nx, Ny, exact_solution)

println()
println("Порівняння ітераційних методів:")
println("--------------------------------------------------------------------------")
@printf("%20s %15s %20s %20s\n", "Метод", "Ітерації", "Критерій", "max error")
println("--------------------------------------------------------------------------")
@printf("%20s %15d %20.10e %20.10e\n", "Jacobi", iter_jacobi, conv_jacobi, error_jacobi)
@printf("%20s %15d %20.10e %20.10e\n", "Gauss-Seidel", iter_gs, conv_gs, error_gs)
println("--------------------------------------------------------------------------")

println()
println("Залежність кількості ітерацій методу Гаусса-Зейделя від розміру сітки:")
println("----------------------------------------------------------------")
@printf("%12s %12s %15s %20s\n", "Nx × Ny", "h", "Ітерації", "max error")
println("----------------------------------------------------------------")

grid_values = [10, 50, 100, 200]
gs_iterations = Int[]

for Nk in grid_values
    hxk = (bx - ax) / Nk
    hyk = (by - ay) / Nk

    xk = collect(range(ax, bx, length = Nk + 1))
    yk = collect(range(ay, by, length = Nk + 1))

    Uk = build_initial_grid(xk, yk, Nk, Nk, bottom, top, left, right)

    coef_xk = 1.0 / hxk^2
    coef_yk = 1.0 / hyk^2
    denom_k = 2.0 * (coef_xk + coef_yk)

    U_res, iter_res, _ = gauss_seidel_method(
        Uk, Nk, Nk, coef_xk, coef_yk, denom_k, eps, max_iter
    )

    err_res = max_error(U_res, xk, yk, Nk, Nk, exact_solution)

    push!(gs_iterations, iter_res)

    @printf("%5d × %-5d %12.6f %15d %20.10e\n",
        Nk, Nk, hxk, iter_res, err_res)
end

println("----------------------------------------------------------------")

plot(
    grid_values,
    gs_iterations,
    marker = :circle,
    linewidth = 2,
    xlabel = "N = Nx = Ny",
    ylabel = "Кількість ітерацій",
    title = "Залежність кількості ітерацій від розміру сітки",
    legend = false,
    size = (900, 600),
    dpi = 300
)

savefig("laplace_gs_iterations.png")

heatmap(
    x,
    y,
    U_gs',
    xlabel = "x",
    ylabel = "y",
    title = "Рівняння Лапласа: теплова карта",
    colorbar_title = "u(x,y)",
    aspect_ratio = 1,
    size = (850, 700),
    dpi = 300
)

savefig("laplace_2d_heatmap.png")

surface(
    x,
    y,
    U_gs',
    xlabel = "x",
    ylabel = "y",
    zlabel = "u(x,y)",
    title = "Рівняння Лапласа: поверхня розв'язку",
    size = (900, 700),
    dpi = 300
)

savefig("laplace_2d_surface.png")

contour(
    x,
    y,
    U_gs',
    xlabel = "x",
    ylabel = "y",
    title = "Рівняння Лапласа: лінії рівня",
    fill = true,
    colorbar_title = "u(x,y)",
    aspect_ratio = 1,
    size = (850, 700),
    dpi = 300
)

savefig("laplace_2d_contour.png")

println()