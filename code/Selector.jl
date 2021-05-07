include("Graph.jl")
include("SimMetric.jl")

using SparseArrays
using LinearAlgebra
import Laplacians

function getSimSelectorByS(S)
    n = size(S, 1)
    f(u, k) = begin
        ret = Vector{Int32}()
        A = Vector{Tuple{Int32, Float64}}()
        for i = 1 : n
            if i != u
                push!(A, (i, S[u, i]))
            end
        end
        sort!(A, by=x->-x[2])
        for i = 1 : k
            push!(ret, A[i][1])
        end
        return ret
    end
    return f
end

function RoleSimSelector(G)
    S = RoleSim(G)
    return getSimSelectorByS(S)
end

function ForestSimSelector(G)
    S = ForestSim(G)
    return getSimSelectorByS(S)
end

function ForestSimApproxSelector(G)
    S = approxForest(G)
    A = Vector{Tuple{Int32, Float64}}()
    foreach(i -> push!(A, (i, S[i])), 1 : G.n)
    sort!(A, by=x->x[2])
    pos = zeros(Int32, G.n)
    foreach(i -> pos[A[i][1]] = i, 1 : G.n)

    sim(x, y) = min(S[x], S[y]) / max(S[x], S[y])

    f(x, k) = begin
        l = pos[x]-1
        r = pos[x]+1
        ret = Vector{Int32}()
        for rep = 1 : k
            if l == 0
                push!(ret, A[r][1])
                r += 1
            elseif r > G.n
                push!(ret, A[l][1])
                l -= 1
            elseif sim(x, A[l][1]) > sim(x, A[r][1])
                push!(ret, A[l][1])
                l -= 1
            else
                push!(ret, A[r][1])
                r += 1
            end
        end
        return ret
    end

    return f
end

function StructSimSelector(G; maxk = 3, wi = 0.5)
    n = G.n
    d = degrees(G)
    B, dk = calcBins(G, maxk)
    maxb = size(B, 3)

    f(x, k) = begin
        S = zeros(n)
        foreach(i -> S[i] = min(d[i], d[x]) / max(d[i], d[x]), 1 : n)
        for k = 1 : maxk
            for i = 1 : n
                delta = 0.0
                foreach(z -> delta += min(B[k, i, z], B[k, x, z]), 1 : maxb)
                delta /= max(dk[k, i], dk[k, x])
                S[i] = (1.0 - wi) * S[i] + wi * delta
            end
        end
        A = Vector{Tuple{Int32, Float64}}()
        for i = 1 : n
            if i != x
                push!(A, (i, S[i]))
            end
        end
        sort!(A, by=x->-x[2])
        ret = Vector{Int32}()
        for i = 1 : k
            push!(ret, A[i][1])
        end
        return ret
    end

    return f
end