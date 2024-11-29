using RABBIT
using Test 

@testset "RABBIT" begin 

    works = try
        output = RABBIT.read_outputs(@__DIR__)
        @test isapprox(output.powe_data[:,1][1], 208692, rtol=0.10)
        @test isapprox(output.powi_data[:,1][1], 482856, rtol=0.10)
        @test isapprox(output.torqdepo_data[:, :, 1][1], 0.20518352, rtol=0.10)
        @test isapprox(output.nrate_data[:, :, 1][1], 29.10558, rtol=0.10)
        true
    catch
        false
    end
    @test works==true
end