include("Graph.jl")
include("SimMetric.jl")
include("Exp.jl")

using Printf

function runMain(ags1; resultDir = "../result/")
    G = readLabeledGraph(ags1)
    println(ags1, " ", G.n, " ", G.m)
    # RoleSim
    out1 = open(resultDir*ags1*"RoleSim.txt", "w")
    @time result = topK(G, RoleSim(G))
    for i = 1 : size(result, 1)
        println(out1, i, " ", result[i])
    end
    close(out1)
    # ForestSim-Exact
    out2 = open(resultDir*ags1*"ForestSimExact.txt", "w")
    @time result = topK(G, ForestSim(G))
    for i = 1 : size(result, 1)
        println(out2, i, " ", result[i])
    end
    close(out2)
    # ForestSim-Approx
    out3 = open(resultDir*ags1*"ForestSimApprox.txt", "w")
    @time result = topKofForestSimApprox(G)
    for i = 1 : size(result, 1)
        println(out3, i, " ", result[i])
    end
    close(out3)
    # StructSim
    out4 = open(resultDir*ags1*"StructSim.txt", "w")
    @time result = topK(G, StructSim(G))
    for i = 1 : size(result, 1)
        println(out4, i, " ", result[i])
    end
    close(out4)
end

runMain("karate")
runMain("brazil-airports") 
runMain("europe-airports")
runMain("expressway-starbucks")
runMain("expressway-status")
runMain("actor")