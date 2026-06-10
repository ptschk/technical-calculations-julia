using ModelingToolkit
using MethodOfLines
using DomainSets
using OrdinaryDiffEq
using Plots

@parameters t x
@variables u(..)

Dt = Differential(t)
Dxx = Differential(x)^2

a = 1.0
L = 1.0
T = 0.1
Nx = 20
dx = L / Nx

eq = Dt(u(t, x)) ~ a * Dxx(u(t, x))

bcs = [
    u(0, x) ~ sin(pi * x / L),
    u(t, 0) ~ 0.0,
    u(t, L) ~ 0.0
]

domains = [
    t ∈ Interval(0.0, T),
    x ∈ Interval(0.0, L)
]

@named pdesys = PDESystem(eq, bcs, domains, [t, x], [u(t, x)])

discretization = MOLFiniteDifference([x => dx], t)

prob = discretize(pdesys, discretization)

sol = solve(prob, Tsit5(), saveat = 0.01)

println("Задачу теплопровідності розв'язано за допомогою MethodOfLines.jl")
println("Кількість часових точок: ", length(sol.t))

xs = collect(0.0:dx:L)

U = sol[u(t, x)]

if size(U, 1) == length(xs)
    u_final = U[:, end]
else
    u_final = U[end, :]
end

exact_solution(x, t) = exp(-a * (pi / L)^2 * t) * sin(pi * x / L)
u_exact = [exact_solution(xi, T) for xi in xs]

max_err = maximum(abs.(u_final .- u_exact))

println("Максимальна похибка при t = T:")
println(max_err)

plot(
    xs,
    u_exact,
    label = "Точний розв'язок",
    xlabel = "x",
    ylabel = "u(x,T)",
    title = "MethodOfLines.jl: рівняння теплопровідності"
)

plot!(
    xs,
    u_final,
    label = "MethodOfLines.jl",
    linestyle = :dash,
    marker = :circle
)

savefig("heat_methodoflines.png")