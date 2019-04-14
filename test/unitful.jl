using Unitful
using Unitful: nm, mm

function isRayTransferMatrix(m::Matrix)
    if size(m) != (2,2)
        return false
    end
    # assert dimensional compatibility
    if dimension(m[1,1]) != NoDims
        return false
    end
    if dimension(m[2,2]) != NoDims
        return false
    end
    if dimension(m[1,2]*m[2,1]) != NoDims
        return false
    end
    return true
end

function isUnitfulRayTransferMatrix(m::Matrix)
    if size(m) != (2,2)
        return false
    end
    # assert all dimensions
    if dimension(m[1,1]) != NoDims
        return false
    end
    if dimension(m[2,2]) != NoDims
        return false
    end
    if dimension(m[1,2]) != dimension(1Unitful.m)
        @warn "B is not a Length" m[1,2]
        return false
    end
    if dimension(m[2,1]) != dimension(1/Unitful.m)
        @warn "C is not an inverse Length" m[2,1]
        return false
    end
    return true
end

# general test setup, unitful
fu = 100mm
Lu = 1000mm
w0u = 1mm
λu = 1000nm
expander_2x_u = [ThinLens(f=fu), FreeSpace(3fu), ThinLens(f=2fu)]
system_u = [expander_2x_u; FreeSpace(Lu); reverse(expander_2x_u)]
beam_u = Beam(λ = λu, w0 = w0u)

@testset "Ray Transfer Matrices" begin
    @test isUnitfulRayTransferMatrix(RTM(ThinLens(f=100mm)))
    @test isUnitfulRayTransferMatrix(RTM(FreeSpace(500mm)))
    @test isUnitfulRayTransferMatrix(
        RTM(ThinLens(f=100mm)) * RTM(FreeSpace(500mm))
    )
    @test isUnitfulRayTransferMatrix(RTM([
        ThinLens(f=100mm),
        FreeSpace(500mm)
    ]))
end
@testset "transform" begin
    @test dimension(beam_u.z) == dimension(1Unitful.m)
    @test transform(ThinLens(f=100mm), beam_u) isa Beam
    @test transform(FreeSpace(500mm), beam_u) isa Beam
    @test transform(system_u, beam_u) isa Beam
end
@testset "general" begin
    @test RTM(Sag(system_u)) == RTM(Tan(system_u))
    @test isUnitfulRayTransferMatrix(RTM(Sag(system_u)))
    @test isUnitfulRayTransferMatrix(RTM(Tan(system_u)))
    @test beamtrace(system_u, beam_u) isa Vector{Beam}
    @test spotradiusfunc(expander_2x_u, beam_u)(3fu) ≈
        2w0u rtol=0.01
end