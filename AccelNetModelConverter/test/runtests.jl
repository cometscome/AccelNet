using Test
using AccelNetModelConverter

fixture = joinpath(@__DIR__, "fixtures", "n2p2-julia")

@testset "n2p2 ↔ AccelNet round trip" begin
    original = read_n2p2_model(fixture)
    mktempdir() do temporary
        accelnet = joinpath(temporary, "accelnet")
        roundtrip = joinpath(temporary, "n2p2")
        networks = n2p2_to_accelnet(fixture, accelnet)
        @test Set(basename.(networks)) == Set(["O.nn.ascii", "Ti.nn.ascii"])
        ti = read_atomic_network(joinpath(accelnet, "Ti.nn.ascii"))
        @test ti.weights[1:2] == [2.0, 4.0]
        accelnet_to_n2p2(networks, roundtrip)
        converted = read_n2p2_model(roundtrip)
        @test converted.species == original.species
        @test converted.nodes == original.nodes
        @test converted.activations == original.activations
        @test converted.atomic_references == original.atomic_references
        @test converted.mean_energy == original.mean_energy
        @test converted.conv_energy == original.conv_energy
        @test ti.cutoff_type == original.cutoff_type == converted.cutoff_type == 9
        @test ti.cutoff_alpha == original.cutoff_alpha == converted.cutoff_alpha == 0.3
        for symbol in original.species
            @test converted.functions[symbol] == original.functions[symbol]
            @test converted.shifts[symbol] == original.shifts[symbol]
            @test converted.scales[symbol] == original.scales[symbol]
            @test converted.weights[symbol] == original.weights[symbol]
        end
    end
end
