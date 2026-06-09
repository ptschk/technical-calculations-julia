using Printf
using Plots

function thomas_algorithm(lower, main, upper, rhs)
    n = length(rhs)

    c = copy(upper)
    d = copy(rhs)
    b = copy(main)

    for i in 2:n
        m = lower[i-1] / b[i-1]
        b[i] -= m * c[i-1]
        d[i] -= m * d[i-1]
    end

    x = zeros(n)
    x[n] = d[n] / b[n]

    for i in n-1:-1:1
        x[i] = (d[i] - c[i] * x[i+1]) / b[i]
    end

    return x
end

function explicit_scheme(a, h, tau, Nx, Nt, x, t, initial_condition, left_boundary, right_boundary)
    r = a * tau / h^2
    U = zeros(Nx + 1, Nt + 1)

    for i in 1:Nx+1
        U[i, 1] = initial_condition(x[i])
    end

    for n in 1:Nt+1
        U[1, n] = left_boundary(t[n])
        U[Nx+1, n] = right_boundary(t[n])
    end

    for n in 1:Nt
        for i in 2:Nx
            U[i, n+1] = U[i, n] + r * (U[i-1, n] - 2U[i, n] + U[i+1, n])
        end
    end

    return U
end

function implicit_scheme(a, h, tau, Nx, Nt, x, t, initial_condition, left_boundary, right_boundary)
    r = a * tau / h^2
    U = zeros(Nx + 1, Nt + 1)

    for i in 1:Nx+1
        U[i, 1] = initial_condition(x[i])
    end

    for n in 1:Nt+1
        U[1, n] = left_boundary(t[n])
        U[Nx+1, n] = right_boundary(t[n])
    end

    m = Nx - 1

    lower = fill(-r, m - 1)
    main = fill(1 + 2r, m)
    upper = fill(-r, m - 1)

    for n in 1:Nt
        rhs = copy(U[2:Nx, n])

        rhs[1] += r * left_boundary(t[n+1])
        rhs[end] += r * right_boundary(t[n+1])

        U[2:Nx, n+1] = thomas_algorithm(lower, main, upper, rhs)
    end

    return U
end

function compute_exact_solution(x, t, Nx, Nt, exact_solution)
    U_exact = zeros(Nx + 1, Nt + 1)

    for n in 1:Nt+1
        for i in 1:Nx+1
            U_exact[i, n] = exact_solution(x[i], t[n])
        end
    end

    return U_exact
end

function max_error(U_num, U_exact)
    return maximum(abs.(U_num .- U_exact))
end

println("Чисельне розв'язування рівняння теплопровідності")
println("u_t = a * u_xx")

println()
println("Введіть коефіцієнт теплопровідності a:")
a = parse(Float64, readline())

println("Введіть довжину області L:")
L = parse(Float64, readline())

println("Введіть кінцевий час T:")
T = parse(Float64, readline())

println("Введіть кількість поділів по простору Nx:")
Nx = parse(Int, readline())

println("Введіть кількість кроків по часу Nt:")
Nt = parse(Int, readline())

if a <= 0
    error("Коефіцієнт теплопровідності a має бути додатним.")
end

if L <= 0
    error("Довжина області L має бути додатною.")
end

if T <= 0
    error("Кінцевий час T має бути додатним.")
end

if Nx < 2
    error("Nx має бути більшим або рівним 2.")
end

if Nt < 1
    error("Nt має бути більшим або рівним 1.")
end

h = L / Nx
tau = T / Nt

x = collect(range(0.0, L, length = Nx + 1))
t = collect(range(0.0, T, length = Nt + 1))

r = a * tau / h^2

initial_condition(x) = sin(pi * x / L)
left_boundary(t) = 0.0
right_boundary(t) = 0.0
exact_solution(x, t) = exp(-a * (pi / L)^2 * t) * sin(pi * x / L)

println()
println("Параметри задачі:")
println("a = ", a)
println("L = ", L)
println("T = ", T)
println("Nx = ", Nx)
println("Nt = ", Nt)
println("h = ", h)
println("tau = ", tau)
println("r = a*tau/h^2 = ", r)

println()
println("Початкова умова:")
println("u(x,0) = sin(πx/L)")

println()
println("Граничні умови:")
println("u(0,t) = 0")
println("u(L,t) = 0")

println()
println("Точний розв'язок:")
println("u(x,t) = exp(-a(π/L)^2 t) * sin(πx/L)")

explicit_is_stable = r <= 0.5

if explicit_is_stable
    U_explicit = explicit_scheme(
        a, h, tau, Nx, Nt, x, t,
        initial_condition,
        left_boundary,
        right_boundary
    )
else
    U_explicit = nothing
    println()
    println("Явна схема не обчислюється, оскільки r = ", r, " > 0.5")
end

U_implicit = implicit_scheme(
    a, h, tau, Nx, Nt, x, t,
    initial_condition,
    left_boundary,
    right_boundary
)

U_exact = compute_exact_solution(x, t, Nx, Nt, exact_solution)

error_implicit = max_error(U_implicit, U_exact)

println()
println("Порівняння з точним розв'язком:")
println("---------------------------------------------")
@printf("%20s %20s\n", "Метод", "max error")
println("---------------------------------------------")

if explicit_is_stable
    error_explicit = max_error(U_explicit, U_exact)
    @printf("%20s %20.10e\n", "Явна схема", error_explicit)
else
    @printf("%20s %20s\n", "Явна схема", "нестійка")
end

@printf("%20s %20.10e\n", "Неявна схема", error_implicit)
println("---------------------------------------------")

final_time_index = Nt + 1

plot(
    x,
    U_exact[:, final_time_index],
    label = "Точний розв'язок",
    lw = 2,
    linestyle = :dash
)

if explicit_is_stable
    plot!(
        x,
        U_explicit[:, final_time_index],
        label = "Явна схема",
        lw = 2,
        marker = :circle
    )
end

plot!(
    x,
    U_implicit[:, final_time_index],
    label = "Неявна схема",
    lw = 2,
    marker = :diamond
)

title!("Теплопровідність: розв'язок при t = T")
xlabel!("x")
ylabel!("u(x,T)")

savefig("heat_final_time_comparison.png")

if explicit_is_stable
    heatmap(
        x,
        t,
        U_explicit',
        xlabel = "x",
        ylabel = "t",
        title = "Теплопровідність: heatmap явної схеми",
        colorbar_title = "u(x,t)"
    )

    savefig("heat_explicit_heatmap.png")
end

heatmap(
    x,
    t,
    U_implicit',
    xlabel = "x",
    ylabel = "t",
    title = "Теплопровідність: heatmap неявної схеми",
    colorbar_title = "u(x,t)"
)

savefig("heat_implicit_heatmap.png")

println()
println("Таблиця похибок для різних сіток:")
println("--------------------------------------------------------------------------")
@printf("%10s %10s %12s %20s %20s\n", "Nx", "Nt", "r", "explicit error", "implicit error")
println("--------------------------------------------------------------------------")

grid_values = [20, 40, 80]

for Nxk in grid_values
    hk = L / Nxk
    target_r = 0.4
    Ntk = ceil(Int, a * T / (target_r * hk^2))
    tauk = T / Ntk
    rk = a * tauk / hk^2

    xk = collect(range(0.0, L, length = Nxk + 1))
    tk = collect(range(0.0, T, length = Ntk + 1))

    U_exp_k = explicit_scheme(
        a, hk, tauk, Nxk, Ntk, xk, tk,
        initial_condition,
        left_boundary,
        right_boundary
    )

    U_imp_k = implicit_scheme(
        a, hk, tauk, Nxk, Ntk, xk, tk,
        initial_condition,
        left_boundary,
        right_boundary
    )

    U_exact_k = compute_exact_solution(xk, tk, Nxk, Ntk, exact_solution)

    err_exp_k = max_error(U_exp_k, U_exact_k)
    err_imp_k = max_error(U_imp_k, U_exact_k)

    @printf("%10d %10d %12.6f %20.10e %20.10e\n",
        Nxk, Ntk, rk, err_exp_k, err_imp_k)
end

println("--------------------------------------------------------------------------")

selected_times = [0.0, T/4, T/2, 3T/4, T]

p = plot(
    xlabel = "x",
    ylabel = "u(x,t)",
    title = "Теплопровідність: розв'язок у різні моменти часу"
)

for time_value in selected_times
    index = argmin(abs.(t .- time_value))

    plot!(
        p,
        x,
        U_implicit[:, index],
        label = "t = $(round(t[index], digits=4))",
        lw = 2
    )
end

savefig(p, "heat_solution_different_times.png")

println()