using CorrelatedEquilibria
using Test
using Plots
using JuMP, CDDLib, LazySets, LinearAlgebra

@testset "CorrelatedEquilibria.jl" begin
    println("Test single state game")
    path = joinpath(pwd(), "test/test_game1.yml")
    gamma = 0.9
    
    # test succesful creation of game
    @test_nowarn game = build_game(path)
    game = build_game(path)
    @test CorrelatedEquilibria.maxPayoff(game) == 7.0
    @test CorrelatedEquilibria.minPayoff(game) == 0.0

    # test algorithm
    @test_nowarn (V,X, animArr) = computeCE(game, gamma, 10, true, 0.01)

    println("Test multiple state game")
    path = joinpath(pwd(), "test/test_game2.yml")
    gamma = 0.9
    
    # test succesful creation of game
    @test_nowarn game = build_game(path)
    game = build_game(path)
    @test CorrelatedEquilibria.maxPayoff(game) == 2.0
    @test CorrelatedEquilibria.minPayoff(game) == -2.0

    # test algorithm
    @test_nowarn (V,X, animArr) = computeCE(game, gamma, 10, true, 0.01)
end