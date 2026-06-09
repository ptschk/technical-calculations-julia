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

function solve_bvp_1d(a, b, N, alpha, beta, f, exact)
    h = (b - a) / N
    x = collect(range(a, b, length=N+1))

    n = N - 1

    lower = fill(-1.0 / h^2, n - 1)
    main  = fill( 2.0 / h^2, n)
    upper = fill(-1.0 / h^2, n - 1)

    rhs = [f(x[i]) for i in 2:N]

    rhs[1] += alpha / h^2
    rhs[end] += beta / h^2

    u_inner = thomas_algorithm(lower, main, upper, rhs)

    u_num = vcat(alpha, u_inner, beta)
    u_exact = [exact(xi) for xi in x]

    err = maximum(abs.(u_num .- u_exact))

    return x, u_num, u_exact, err, h
end

println("Введіть ліву межу відрізка a:")
a = parse(Float64, readline())

println("Введіть праву межу відрізка b:")
b = parse(Float64, readline())

println("Введіть кількість поділів N:")
N = parse(Int, readline())

println("Введіть ліву граничну умову u(a) = alpha:")
alpha = parse(Float64, readline())

println("Введіть праву граничну умову u(b) = beta:")
beta = parse(Float64, readline())

if b <= a
    error("Права межа b має бути більшою за ліву межу a.")
end

if N < 2
    error("N має бути більшим або рівним 2.")
end

f(x) = pi^2 * sin(pi * x)
exact(x) = sin(pi * x)

println()
println("Параметри задачі:")
println("a = ", a)
println("b = ", b)
println("N = ", N)
println("u(a) = ", alpha)
println("u(b) = ", beta)

x, u_num, u_exact, err, h = solve_bvp_1d(a, b, N, alpha, beta, f, exact)

println()
println("Побудова сітки:")
println("Крок сітки h = ", h)
println("Кількість вузлів = ", N + 1)

println()
println("Апроксимація другої похідної:")
println("-u''(xᵢ) ≈ (-uᵢ₋₁ + 2uᵢ - uᵢ₊₁) / h²")

println()
println("Порівняння з точним розв'язком:")
println("Максимальна похибка = ", err)

println()
println("Таблиця похибок:")
println("------------------------------------------------")
@printf("%8s %12s %18s\n", "N", "h", "max error")
println("------------------------------------------------")

N_values = [10, 20, 40, 80, 160]

for Nk in N_values
    _, _, _, err_k, h_k = solve_bvp_1d(a, b, Nk, alpha, beta, f, exact)
    @printf("%8d %12.6f %18.10e\n", Nk, h_k, err_k)
end

println("------------------------------------------------")

plot(x, u_num,
    label = "Чисельний розв'язок",
    lw = 2,
    marker = :circle
)

plot!(x, u_exact,
    label = "Точний розв'язок",
    lw = 2,
    linestyle = :dash
)

title!("Одновимірна крайова задача")
xlabel!("x")
ylabel!("u(x)")
grid = true

savefig("bvp_1d_result.png")

println()