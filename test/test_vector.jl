# Vector integrands (NUMFUN ≥ 2): all components share the same subdivision.

@testset "Vector integrand (NUMFUN = 2)" begin
    # ∫₀¹∫₀¹ [x²+y², x·y] = [2/3, 1/4]
    f(u, _) = [u[1]^2 + u[2]^2, u[1] * u[2]]
    prob = IntegralProblem(f, (zeros(2), ones(2)))
    sol = solve(prob, DecuhrAlgorithm(); abstol = 1.0e-9)
    @test isapprox(sol.u[1], 2 / 3; atol = 1.0e-7)
    @test isapprox(sol.u[2], 1 / 4; atol = 1.0e-7)
end

@testset "Singular vector integrand (NUMFUN = 2)" begin
    # A singular and a smooth component side by side:
    # ∫₀¹∫₀¹ [(x·y)^(-1/2), x+y] = [4, 1].
    f(u, _) = [(u[1] * u[2])^(-0.5), u[1] + u[2]]
    prob = IntegralProblem(f, (zeros(2), ones(2)))
    sol = solve(
        prob,
        DecuhrAlgorithm(singul = 2, alpha = -0.5);
        abstol = 1.0e-8, reltol = 1.0e-6, maxiters = 1_000_000
    )
    @test sol.u isa AbstractVector
    @test length(sol.u) == 2
    @test isapprox(sol.u[1], 4.0; rtol = 1.0e-3)
    @test isapprox(sol.u[2], 1.0; rtol = 1.0e-6)
end
