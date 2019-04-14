using ABCDBeamTrace
using Test
import Pkg

# general test setup
f = 100e-3
L = 1000e-3
w0 = 1.0e-3
expander_2x = [ThinLens(f=f), FreeSpace(3f), ThinLens(f=2f)]
system = [expander_2x; FreeSpace(L); reverse(expander_2x)]
beam = GaussianBeam(λ = 1000e-9, w0 = w0)
# gob: Geometric Optics Beam
n_gob = 1.33
k_gob = -0.01
x_gob = 5e-3
gob = GeometricBeam(x = x_gob, n = n_gob, slope = k_gob)

@testset "beamtrace" begin
    @test beamtrace(system, beam) isa Vector{GaussianBeam}
    @test beamtrace(system, gob) isa Vector{GeometricBeam}
end

@testset "discretize" begin
    @test length(discretize(expander_2x,10)) == 12
end

@testset "spotradiusfunc" begin
    @test spotradiusfunc(expander_2x, beam) isa Function
    @test spotradiusfunc(expander_2x, beam)(0) == w0
    @test spotradiusfunc(expander_2x, beam)(3f) ≈ 2w0 rtol=0.01
end

@testset "transform" begin
    @test transform(system, beam) isa GaussianBeam
    @test transform(system, gob) isa GeometricBeam
end

@testset "planes" begin
    # system with no astigmatism: sagittal and parallel plane should
    # be identical (within numerical error)
    @test Matrix(Sag(system)) ≈ Matrix(Tan(system))
end

@testset "accessors" begin
    @test ior(gob) ≈ n_gob
    @test slope(gob) ≈ k_gob
    @test radialpos(gob) ≈ x_gob
    @test location(gob) == 0
    @test ior(beam) ≈ 1.0
    @test slope(beam) == 0
    @test radialpos(beam) == 0
    @test location(beam) == 0
end

@testset "unitful" begin
    if haskey(Pkg.installed(), "Unitful")
        include("unitful.jl")
    else
        @warn string(
            "Skipping Unitful Tests because ",
            "package Unitful is not installed"
        )
    end
end
