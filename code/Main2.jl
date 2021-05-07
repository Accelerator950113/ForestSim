include("Graph.jl")
include("Exp2.jl")

using Printf

function runMain(ags1, ags2; resultDir = "../result/")
    lg = open(resultDir*"log.txt", "a")
    G = readUnLabeledGraph(ags1)
    println(lg, ags1, " ", G.n, " ", G.m)
    ds = Dict{String, String}()
    ds["1"] = "RoleSim : "
    ds["2"] = "StructSim : "
    ds["3"] = "ForestSimExact : "
    ds["4"] = "ForestSimApprox : "
    df = Dict{String, Function}()
    df["1"] = testRoleSim
    df["2"] = testStructSim
    df["3"] = testForestSim
    df["4"] = testForestSimApprox
    for x in split(ags2, ',')
        println(lg, ds[x], df[x](G, 10))
    end
    close(lg)
end

runMain("Diseasome", "1,2,3,4")
runMain("EmailUniv", "1,2,3,4")
runMain("Hamster", "1,2,3,4")
runMain("GridWorm", "1,2,3,4")
runMain("GrQc", "1,2,3,4")
runMain("Erdos992", "1,2,3,4")
runMain("Reality", "1,2,3,4")
runMain("Dmela", "1,2,3,4")
runMain("HepPh", "2,3,4")
runMain("PagesCompany", "2,3,4")
runMain("AstroPh", "2,3,4")
runMain("Gplus", "2,3,4")
runMain("Brightkite", "2,4")
runMain("BlogCatalog", "2,4")
runMain("Douban", "2,4")
runMain("MathSciNet", "4")
runMain("Flickr", "4")
runMain("IMDB", "4")
runMain("YoutubeSnap", "4")
runMain("Flixster", "4")
