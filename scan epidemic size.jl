using Random
using Distributions
using CSV
using Dates
using DataFrames
using StatsBase
using Plots
using DimensionalData
# the file can be found in the same repository
tre = CSV.File("https://github.com/wl415/Summer-project-on-time-dependent-infectivity-and-detectability/blob/main/smallLandscape.csv",header=false)|>Tables.matrix
distance = zeros(1111, 1111)
for i in 1:1111
    for j in 1:1111
        distance[i, j] = hypot(tre[i, 1] - tre[j, 1], tre[i, 2] - tre[j, 2])
    end
end

function dispersal(a)
    qq = 1 ./ (1 .+ ((distance ./ a) .^ 2))
    factor = (sum(qq) - 1111)/2
    return qq ./ factor
end


betas = [0, 5.25]
dispersing = dispersal(37)
searchk = [7:20:107;]
searchradius =0:5:140
numrepeats = 1000
simulations2 = DimArray(zeros(length(searchk), length(searchradius), numrepeats), (θ = searchk, rad = searchradius, rep = 1:numrepeats))
mean_simulations2 = DimArray(zeros(length(searchk), length(searchradius)), (θ = searchk, rad = searchradius))
println("start")
for K in searchk
    rnd = Gamma(107/K, K)
    for Rad in searchradius
        println(K," ", Rad)
        for repeats in 1:numrepeats
            println(repeats)
            trees = zeros(1111, 6)
            trees[:, 1:2] = tre[:, 1:2]
            t = 0
            ts = [0]
            Ss = [1101]
            Is = [10]
            Rs = [0]
            initial = sample([1:1111;], 10)
            trees[initial, 3] .= 1
            trees[initial, 4] .= 0
            trees[initial, 6] .= rand(rnd, 10)
            notinf = findall(trees[:, 3] .== 0)
            for i in initial
                trees[notinf, 4] .+= betas[2] * dispersing[i, notinf]
            end
            while Is[end] > 0
                timetosurvey = 90 - (t%90)
                event = cumsum(trees[:, 4])
                rate = event[end]
                if rate > 0
                    timestep = rand(Exponential(1/rate))
                else
                    break
                end
                if timetosurvey < timestep
                    culled = []
                    t += timetosurvey
                    detected = (trees[:, 6] .< t) .&& (trees[:, 3] .== 1)
                    for i in 1:1111
                        if any(distance[detected, i] .<= Rad) && (trees[i, 3] < 2)
                            culled = vcat(culled, i)
                        end 
                    end
                    notinf = findall(trees[:, 3] .== 0)
                    for i in culled
                        if trees[i, 3] == 1
                            trees[notinf, 4] .-= betas[2] .* dispersing[notinf, i] 
                        end
                    end
                    trees[culled, 3] .= 2
                    trees[culled, 4] .= 0
                    ts = vcat(ts, t)
                    Ss = vcat(Ss, sum(trees[:, 3] .== 0))
                    Is = vcat(Is, sum(trees[:, 3] .== 1))
                    Rs = vcat(Rs, sum(trees[:, 3] .== 2))
                else
                    t += timestep
                    ts = vcat(ts, t)
                    whichevent = rate * rand(1)
                    whichosen = findmax(event .> whichevent)[2]
                    trees[whichosen, 3] += 1
                    state = trunc(Int, trees[whichosen, 3])
                    if state == 1
                        trees[whichosen, 5] = t
                        trees[whichosen, 6] = t + rand(rnd)
                    else
                        println("error")
                    end
                    trees[whichosen, 4] = 0
                    db = betas[state + 1] - betas[state]
                    notinf = findall(trees[:, 3] .== 0)
                    trees[notinf, 4] .+= db .* dispersing[whichosen, notinf]
                    Ss = vcat(Ss, sum(trees[:, 3] .== 0))
                    Is = vcat(Is, sum(trees[:, 3] .== 1))
                    Rs = vcat(Rs, sum(trees[:, 3] .== 2))
                end
            end
            simulations2[Dim{:θ}(At(K)), Dim{:rad}(At(Rad)), Dim{:rep}(At(repeats))] = Rs[end]
        end
        mean_simulations2[Dim{:θ}(At(K)), Dim{:rad}(At(Rad))] = mean(simulations2[Dim{:θ}(At(K)), Dim{:rad}(At(Rad))])
    end
end
plot(transpose(mean_simulations))

savefig("fuller optimum radius curve.png")
