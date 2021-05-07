include("Graph.jl")
include("Selector.jl")

function testRoleSim(G, k)
    T = time()
    f = RoleSimSelector(G)
    foreach(i -> x = f(i, k), 1 : G.n)
    return time()-T
end

function testStructSim(G, k)
    T = time()
    f = StructSimSelector(G)
    foreach(i -> x = f(i, k), 1 : G.n)
    return time()-T
end

function testForestSim(G, k)
    T = time()
    f = ForestSimSelector(G)
    foreach(i -> x = f(i, k), 1 : G.n)
    return time()-T
end

function testForestSimApprox(G, k)
    T = time()
    f = ForestSimApproxSelector(G)
    foreach(i -> x = f(i, k), 1 : G.n)
    return time()-T
end