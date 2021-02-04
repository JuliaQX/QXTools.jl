
using Logging

@testset "Test simple amplitude sampling functions" begin
    qubits = 5
    amplitudes = amplitudes_all(qubits)
    @test length(amplitudes) == 2^qubits
    @test all(length.(amplitudes) .== qubits)

    qubits = 5
    N = 10
    amplitudes = amplitudes_uniform(qubits, nothing, N)
    @test length(amplitudes) == N
end

@testset "Test simulation of small circuits" begin
    n = 5
    circ = create_ghz_circuit(n)

    results = with_logger(NullLogger()) do
        run_simulation(circ)
    end

    @test length(results) == 2^n
    @test results["0"^n] ≈ results["1"^n] ≈ 1/sqrt(2)
end