using JLD2
using FileIO
using QXTns

@testset "Test tensor cache" begin
    tc = QXTools.TensorCache()
    a = rand(Float64, 3, 4, 5)
    sym = push!(tc, a)
    @test sym == push!(tc, a)
    @test sym != push!(tc, rand(Float64, 3, 4, 5))
    # test approximate matching
    @test sym == push!(tc, a .+ 0.1*eps(Float64))
    @test sym != push!(tc, a .+ 1.1*eps(Float64))

    @test tc[sym] == a
    @test length(tc) == 3

    # test saving tensor cache
    mktempdir() do path
        save_cache(tc, joinpath(path, "tmp.jld2"))
        loaded_data = load(joinpath(path, "tmp.jld2"))
        @test all([(loaded_data[x] == tc[Symbol(x)]) for x in keys(loaded_data)])
    end
end

@testset "Test dsl header generation" begin

    mktempdir() do path
        tnc = QXTns.TensorNetworkCircuit(1)
        push!(tnc, [1], QXTns.Gates.z())
        push!(tnc, [1], QXTns.Gates.I())
        open(joinpath(path, "test.qx"), "w") do dsl_io
            JLD2.jldopen(joinpath(path, "test.jld2"), "w") do data_io
                QXTools.write_dsl_load_header(tnc, dsl_io, data_io)
            end
        end
        # check the generated DSL is what is expected
        dsl_content = read(open(joinpath(path, "test.qx"), "r"), String)
        @test dsl_content == """# version: $(QXTools.DSL_VERSION)
        outputs 1
        load t1 data_1
        load t2 data_2
        """

        # check the data files match those expected
        data_items = load(joinpath(path, "test.jld2"))
        @test data_items["data_1"] == tensor_data(tnc[:t1], consider_hyperindices=true)
        @test data_items["data_2"] == tensor_data(tnc[:t2], consider_hyperindices=true)
    end
end

@testset "Test write ncon command" begin
    # test case where 2 tensors share one edge which
    # is in the hyper edge group for both tensors
    io = IOBuffer()
    tn = TensorNetwork()
    as = [Index(2), Index(2), Index(2)]
    a_hyper_indices = [[2, 3]]
    bs = [Index(2), as[2]]
    b_hyper_indices = [[1, 2]]
    t1 = push!(tn, QXTensor(as, a_hyper_indices))
    t2 = push!(tn, QXTensor(bs, b_hyper_indices))
    write(io, QXTools.gen_ncon_command(tn, t1, t2, :t3))
    cmd = String(take!(io))
    @test cmd == "ncon t3 1,2 t1 1,2 t2 2\n"

    io = IOBuffer()
    tn = TensorNetwork()
    as = [Index(2), Index(2)]
    a_hyper_indices = [[1, 2]]
    bs = [as[1], as[2]]
    b_hyper_indices = Array{Int64,1}[]
    t1 = push!(tn, QXTensor(as, a_hyper_indices))
    t2 = push!(tn, QXTensor(bs, b_hyper_indices))
    write(io, QXTools.gen_ncon_command(tn, t1, t2, :t3))
    cmd = String(take!(io))
    @test cmd == "ncon t3 0 t1 1 t2 1\n"

    # Test case where A has two sets of hyper indices
    # linked by a hyper index set in B
    tn = TensorNetwork()
    a_i = [Index(2) for _ in 1:9]
    b_i = [a_i[7], Index(2), Index(2), a_i[2]]

    a_hyper = [[3, 4], [5, 6], [7,8]]
    a_hyper_i = [a_i[x] for x in a_hyper]
    b_hyper = [[1, 4], [2, 3]]
    b_hyper_i = [b_i[x] for x in b_hyper]

    push!(tn, QXTensor(a_i, a_hyper))
    push!(tn, QXTensor(b_i, b_hyper))

    io = IOBuffer()
    write(io, QXTools.gen_ncon_command(tn, :t1, :t2, :t3))
    @test String(take!(io)) == "ncon t3 1,2,3,4,5,6 t1 1,2,3,4,5 t2 2,6\n"
end

@testset "Test write view command" begin
    io = IOBuffer()

    tn = TensorNetwork()
    as = [Index(2), Index(2), Index(2)]
    a_hyper_indices = [[2, 3]]
    t1 = push!(tn, QXTensor(as, a_hyper_indices))

    QXTools.write_view_command(io, tn, t1, :t1v1, [as[2]], "v1")

    cmd = String(take!(io))
    @test cmd == "view t1v1 t1 2 v1\n"
end

@testset "Test _permute_and_merge" begin
    # simple case where all independent and merging two sets of two
    cmd = QXTools.ContractCommand("ncon c 3,4,2,1 b 8,9,3,4 c 8,9,2,1")
    cmd.reshape_groups[:] = [1,1,1,1]
    new_index_order = [[1,2], [3,4]]
    QXTools._permute_and_merge!(cmd, new_index_order)
    @test cmd.reshape_groups == [2,2]

    # one set already merged and leaving them as they are
    cmd = QXTools.ContractCommand("ncon c 3,4,2,1 b 8,9,3,4 c 8,9,2,1")
    empty!(cmd.reshape_groups); append!(cmd.reshape_groups, [2,1,1])
    new_index_order = [[1], [2], [3]]
    QXTools._permute_and_merge!(cmd, new_index_order)
    @test cmd.reshape_groups == [2, 1, 1]

    # one set already merged and merging the remaining two
    cmd = QXTools.ContractCommand("ncon c 3,4,2,1 b 8,9,3,4 c 8,9,2,1")
    empty!(cmd.reshape_groups); append!(cmd.reshape_groups, [2,1,1])
    new_index_order = [[1], [3, 2]]
    QXTools._permute_and_merge!(cmd, new_index_order)
    @test cmd.reshape_groups == [2, 2]

    # merging the existing group with one other
    cmd = QXTools.ContractCommand("ncon c 3,4,2,1 b 8,9,3,4 c 8,9,2,1")
    empty!(cmd.reshape_groups); append!(cmd.reshape_groups, [2,1,1])
    new_index_order = [[1, 3], [2]]
    QXTools._permute_and_merge!(cmd, new_index_order)
    @test cmd.reshape_groups == [3, 1]

    # merging the existing group with other two
    cmd = QXTools.ContractCommand("ncon c 3,4,2,1 b 8,9,3,4 c 8,9,2,1")
    empty!(cmd.reshape_groups); append!(cmd.reshape_groups, [2,1,1])
    new_index_order = [[1, 3, 2]]
    QXTools._permute_and_merge!(cmd, new_index_order)
    @test cmd.reshape_groups == [4]
end

@testset "Test find overlaps" begin
    # case with more than one output group
    @test QXTools.find_overlaps([[1,2],[3],[4,5]], [1,2,3,4,5]) == [[1,2], [3], [4,5]]
    @test QXTools.find_overlaps([[1,2],[3],[4,5]], [1,2,4,5]) == [[1,2], [4,5]]
    @test QXTools.find_overlaps([[5,4],[3],[1,2]], [1,2,8,4,5]) == [[1,2],[4],[5]]

    # corner cases with empty arrays
    @test QXTools.find_overlaps([[5,4],[3],[1,2]], Int[]) == Vector{Vector{Int}}()
end

@testset "Test join output idxs function" begin
    cmd = QXTools.ContractCommand("ncon c 3,4,2,1 b 8,9,3,4 c 8,9,2,1")
    empty!(cmd.reshape_groups); append!(cmd.reshape_groups, [2,1,1])
    QXTools.join_output_idxs!(cmd, [[3,4]])
    @test cmd.output_idxs == [3,2,1] && cmd.reshape_groups == [1,1,1]

    cmd = QXTools.ContractCommand("ncon c 3,4,2,1,6,7,8 b 8,9,3,4,2 c 8,9,2,1")
    empty!(cmd.reshape_groups); append!(cmd.reshape_groups, [7])
    QXTools.join_output_idxs!(cmd, [[3,4,2]])
    @test cmd.output_idxs == [3,1,6,7,8] && cmd.reshape_groups == [5]
end

@testset "Test remove repeated in output" begin
    cmd = QXTools.ContractCommand("ncon c 3,4,2,1 b 0 c 0")
    QXTools._remove_repeated_in_output!(cmd, [1,2,3,1])
    @test cmd.output_idxs == [3,4,2]
    @test cmd.reshape_groups == [1,1,1]
end