include("Graph.jl")
include("SimMetric.jl")

function topK(G, S)
    dl = Dict{Int32, Int32}()
    foreach(x -> dl[x] = haskey(dl, x) ? dl[x]+1 : 1, G.L)
    maxk = maximum(values(dl))-1
    tot = zeros(maxk)
    totn = zeros(maxk)
    for i = 1 : G.n
        tp = Vector{Tuple{Int32, Float64}}()
        foreach(j -> (i != j) ? push!(tp, (j, S[i, j])) : nothing, 1 : G.n)
        sort!(tp, by=x->-x[2])
        s = 0.0
        for j = 1 : maxk
            s += (G.L[i] == G.L[tp[j][1]]) ? 1 : 0
            totn[j] += 1
            tot[j] += (s / j)
        end
    end
    ret = zeros(maxk)
    foreach(i -> ret[i] = tot[i] / totn[i], 1 : maxk)
    return ret
end

function topKofForestSimApprox(G)
    dl = Dict{Int32, Int32}()
    foreach(x -> dl[x] = haskey(dl, x) ? dl[x]+1 : 1, G.L)
    maxk = maximum(values(dl))-1
    tot = zeros(maxk)
    totn = zeros(maxk)
    S = approxForest(G)
    A = Vector{Tuple{Int32, Float64}}()
    foreach(i -> push!(A, (i, S[i])), 1 : G.n)
    sort!(A, by=x->x[2])
    pos = zeros(Int32, G.n)
    foreach(i -> pos[A[i][1]] = i, 1 : G.n)

    sim(x, y) = min(S[x], S[y]) / max(S[x], S[y])

    predict(x, k) = begin
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

    for i = 1 : G.n
        cntk = maxk
        tp = predict(i, cntk)
        s = 0.0
        for j = 1 : cntk
            s += (G.L[i] == G.L[tp[j]]) ? 1 : 0
            totn[j] += 1
            tot[j] += (s / j)
        end
    end

    ret = zeros(maxk)
    foreach(i -> ret[i] = tot[i] / totn[i], 1 : maxk)
    return ret
end