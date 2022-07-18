push!(LOAD_PATH, "src")

using UnicodePlots
import JankoUtils

using Test

function plotUnicode(pIn, pOut)

    extremaIn = extrema(pIn; dims=1)
    extremaOut = extrema(pOut; dims=1)

    xMin = min(extremaIn[1][1], extremaOut[1][1])
    xMax = max(extremaIn[1][2], extremaOut[1][2])
    yMin = min(extremaIn[2][1], extremaOut[2][1])
    yMax = max(extremaIn[2][2], extremaOut[2][2])

    xPadding = 0.1 * (xMax - xMin)
    yPadding = 0.1 * (yMax - yMin)
    xMin = xMin - xPadding
    yMin = yMin - yPadding
    xMax = xMax + xPadding
    yMax = yMax + yPadding


    p = UnicodePlots.lineplot([pIn[:,1]; pIn[1,1]], [pIn[:,2]; pIn[1,2]];
                              labels=false, border=:none, grid=false, xlim=(xMin, xMax), ylim=(yMin, yMax))
    UnicodePlots.lineplot!(p, [pOut[:,1]; pOut[1,1]], [pOut[:,2]; pOut[1,2]])
    println(p)
end

function testNoInclusionSimple()
    pIn = [0.0 0;
           1 0;
           1 1;
           0 1]

    pOut = [-1.0 0;
           0.5 0;
           0.5 1.5;
           -1.0 1.5]

    plotUnicode(pIn, pOut)

    @test JankoUtils.intersectIsContained(pIn, pOut) == false
end

function testNoInclusion()
    Nin = 5
    Nout = 10

    Rin = 1
    Rout = 1.1

    pIn = Rin .* [f(t) for t in range(0, 2π, length=Nin+1)[1:end-1], f in [x -> cos(x) - 0.5, sin]]
    pOut = Rout .* [f(t) for t in range(0, 2π, length=Nout+1)[1:end-1], f in [cos, sin]]

    plotUnicode(pIn, pOut)

    @test JankoUtils.intersectIsContained(pIn, pOut) == false
end

function testInclusion()
    Nin = 100
    Nout = 200

    Rin = 1
    Rout = 2

    pIn = Rin .* [f(t) for t in range(0, 2π, length=Nin+1)[1:end-1], f in [cos, sin]]
    pOut = Rout .* [f(t) for t in range(0, 2π, length=Nout+1)[1:end-1], f in [cos, sin]]

    plotUnicode(pIn, pOut)

    @test JankoUtils.intersectIsContained(pIn, pOut) == true
end

function main()
    @info "Test non inclusion simple"
    testNoInclusionSimple()

    @info "Test non inclusion"
    testNoInclusion()

    @info "Test inclusion"
    testInclusion()
end

main()
