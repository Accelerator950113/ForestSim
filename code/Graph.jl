struct LabeledGraph
    n :: Int32 # |V|
    m :: Int32 # |E|
    V :: Array{Int32, 1} # V[i] = Real Index of node i
    L :: Array{Int32, 1} # L[i] = Label of node i
    E :: Array{Tuple{Int32, Int32}, 1} # (u, v) in Edge Set
end

function readLabeledGraph(graphName; dataDir = "../data/") # Read a labeled graph
    # Initialized
    n = 0
    V = Vector{Int32}()
    id = Dict{Int32, Int32}()
    edge = Set{Tuple{Int32, Int32}}()

    getid(x) = begin
        if !haskey(id, x)
            n += 1
            push!(V, x)
            id[x] = n
            return n
        end
        return id[x]
    end

    # read graph from edge list file 
    open(dataDir*graphName*".edgelist") do f1
        for line in eachline(f1)
            buf = split(line)
            if size(buf, 1) != 2
                continue
            end
            u = parse(Int32, buf[1])
            v = parse(Int32, buf[2])
            u1 = getid(u)
            v1 = getid(v)
            if u1 == v1
                continue
            end
            if u1 > v1
                u1, v1 = v1, u1
            end
            push!(edge, (u1, v1))
        end
    end

    # read labels
    L = zeros(Int32, n)
    open(dataDir*graphName*".group") do f2
        for line in eachline(f2)
            buf = split(line)
            if size(buf, 1) != 2
                continue
            end
            x = id[parse(Int32, buf[1])]
            y = parse(Int32, buf[2])
            L[x] = y
        end
    end

    m = length(edge)
    E = Array{Tuple{Int32, Int32}, 1}(undef, m)

    i = 0
    for (u, v) in edge
        i = i + 1
        E[i] = (u, v)
    end

    return LabeledGraph(n, m, V, L, E)
end

function readUnLabeledGraph(graphName; dataDir = "../data/") # Read a graph without labels
    # Initialized
    n = 0
    V = Vector{Int32}()
    id = Dict{Int32, Int32}()
    edge = Set{Tuple{Int32, Int32}}()

    getid(x) = begin
        if !haskey(id, x)
            n += 1
            push!(V, x)
            id[x] = n
            return n
        end
        return id[x]
    end

    # read graph from edge list file 
    open(dataDir*graphName*".txt") do f1
        for line in eachline(f1)
            buf = split(line)
            if size(buf, 1) != 2
                continue
            end
            u = parse(Int32, buf[1])
            v = parse(Int32, buf[2])
            u1 = getid(u)
            v1 = getid(v)
            if u1 == v1
                continue
            end
            if u1 > v1
                u1, v1 = v1, u1
            end
            push!(edge, (u1, v1))
        end
    end

    L = zeros(Int32, n)

    m = length(edge)
    E = Array{Tuple{Int32, Int32}, 1}(undef, m)

    i = 0
    for (u, v) in edge
        i = i + 1
        E[i] = (u, v)
    end

    return LabeledGraph(n, m, V, L, E)
end

function adjacentList(G)
    g = Array{Array{Int32, 1}, 1}(undef, G.n)
    foreach(i -> g[i] = [], 1 : G.n)
    for (u, v) in G.E
        push!(g[u], v)
        push!(g[v], u)
    end
    return g
end

function degrees(G)
    d = zeros(Int, G.n)
    for (u, v) in G.E
        d[u] += 1
        d[v] += 1
    end
    return d
end