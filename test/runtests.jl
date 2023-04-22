using ComplexPaths
using Test

@testset "ComplexPaths.jl" begin
    

    sin_path = Path(t-> t + 2im*sin(t),0,2π)

    @testset "Path Parameterization" begin

        @test sin_path.parameterization(0) ≈ 0 + 0im atol=0.01
        @test sin_path.parameterization(2π) ≈ 2π + 0im atol=0.01
        @test pointsonpath(sin_path,2) ≈ [0.0+0.0im, 2π + 0im] atol=0.01

    end

    @testset "Bad Paths" begin

        @test_throws ErrorException pointsonpath(sin_path, 0)
        @test_throws ErrorException pointsonpath(sin_path, -2)

    end 

    @testset "Path Closure" begin
        
        c_path = CirclePath(1)

        @test isclosed(c_path) == true
        @test isclosed(sin_path) == false

    end

end
