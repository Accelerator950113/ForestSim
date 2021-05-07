include("Graph.jl")

using SparseArrays
using LinearAlgebra
import Laplacians

function RoleSim(G; delta = 1e-6, beta = 0.15, maxIt = 100)
    n = G.n
    g = adjacentList(G)
    d = degrees(G)
    Spre = zeros(n, n)
    S = zeros(n, n)
    Sp = ones(n, n)

    calc(x, y) = begin
        hx = zeros(Bool, d[x])
        hy = zeros(Bool, d[y])
        A = Vector{Tuple{Int32, Int32, Float64}}()
        for i = 1 : d[x], j = 1 : d[y]
            push!(A, (i, j, Sp[g[x][i], g[y][j]]))
        end
        sort!(A, by=x->-x[3])
        ret = 0.0
        for (u, v, w) in A
            if !hx[u] && !hy[v]
                hx[u] = true
                hy[v] = true
                ret += w
            end
        end
        return ret
    end

    ts = 0
    while true
        ts += 1
        for i = 1 : n, j = 1 : n
            S[i, j] = (1-beta)*calc(i, j)/max(d[i], d[j]) + beta
        end
        if (norm(S - Sp) < delta) || (ts > maxIt)
            break
        end
        for i = 1 : n, j = 1 : n
            Sp[i, j] = S[i, j]
        end
    end
    return S
end

function getW(G)
    IpL = zeros(G.n, G.n)
    for (u, v) in G.E
        IpL[u, u] += 1
        IpL[v, v] += 1
        IpL[u, v] -= 1
        IpL[v, u] -= 1
    end
    foreach(i -> IpL[i, i] += 1, 1 : G.n)
    return inv(IpL)
end

function ForestSim(G)
    n = G.n
    W = getW(G)
    S = zeros(n, n)
    for i = 1 : n, j = 1 : n
        S[i, j] = min(W[i, i], W[j, j]) / max(W[i, i], W[j, j])
    end
    return S
end

function getSparseIpL(G)
    d = ones(G.n)
    for (u, v) in G.E
        d[u] += 1
        d[v] += 1
    end
    Is = zeros(Int32, G.m*2+G.n)
    Js = zeros(Int32, G.m*2+G.n)
    Vs = zeros(G.m*2+G.n)
    ID = 0
    for (u, v) in G.E
        ID += 1
        Is[ID] = u
        Js[ID] = v
        Vs[ID] = -1
        Is[ID + G.m] = v
        Js[ID + G.m] = u
        Vs[ID + G.m] = -1
    end
    for i = 1 : G.n
        Is[G.m + G.m + i] = i
        Js[G.m + G.m + i] = i
        Vs[G.m + G.m + i] = d[i]
    end
    return sparse(Is, Js, Vs, G.n, G.n)
end

function getSparseB(G)
    n = G.n
    m = G.m
    Is = zeros(Int32, m*2)
    Js = zeros(Int32, m*2)
    Vs = zeros(m*2)
    ID = 0
    for (u, v) in G.E
        ID += 1
        Is[ID] = ID
        Js[ID] = u
        Vs[ID] = 1
        Is[ID + m] = ID
        Js[ID + m] = v
        Vs[ID + m] = -1
    end
    return sparse(Is, Js, Vs, m, n)
end

function approxForest(G; eps = 0.1)
    n = G.n
    m = G.m
    IpL = getSparseIpL(G)
    B = getSparseB(G)
    k = round(Int, log2(G.n) / eps^2)
    f = Laplacians.approxchol_sddm(IpL, tol=1e-12)
    ret = zeros(G.n)
    for i = 1 : k
        y1 = B' * randn(m)
        y2 = rand(n)
        z1 = f(y1)
        z2 = f(y2)
        foreach(j -> ret[j] += (z1[j]^2 + z2[j]^2), 1 : G.n)
    end
    ret ./= k
    return ret
end

function calcBins(G, maxk)
    maxbin = ceil(Int32, log2(G.n))
    n = G.n
    g = adjacentList(G)
    d = degrees(G)
    B = zeros(Int32, maxk, n, maxbin)
    dk = zeros(Int32, maxk, n)
    for s = 1 : n
        dist = zeros(Int32, n)
        fill!(dist, 1048576)
        dist[s] = 0
        Q = zeros(Int32, n)
        front = 1
        rear = 1
        Q[1] = s
        while front <= rear
            u = Q[front]
            front += 1
            if dist[u] == maxk
                continue
            end
            for v in g[u]
                if dist[u]+1 < dist[v]
                    dist[v] = dist[u] + 1
                    B[dist[v], s, floor(Int32, log2(d[v]))+1] += 1
                    dk[dist[v], s] += 1 
                    rear += 1
                    Q[rear] = v
                end
            end
        end
    end
    return B, dk
end

function StructSim(G; maxk = 3, wi = 0.5)
    n = G.n
    d = degrees(G)
    B, dk = calcBins(G, maxk)
    maxb = size(B, 3)
    S = zeros(n, n)
    for i = 1 : n, j = 1 : n
        S[i, j] = min(d[i], d[j]) / max(d[i], d[j])
    end
    for k = 1 : maxk
        for i = 1 : n, j = 1 : n
            delta = 0.0
            foreach(z -> delta += min(B[k, i, z], B[k, j, z]), 1 : maxb)
            delta /= max(dk[k, i], dk[k, j])
            S[i, j] = (1.0 - wi) * S[i, j] + wi * delta
        end
    end
    return S
end