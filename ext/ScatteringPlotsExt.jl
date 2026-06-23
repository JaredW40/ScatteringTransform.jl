module ScatteringPlotsExt

using ScatteringTransform
using Plots

include(joinpath(@__DIR__, "..", "src", "scatteringplots.jl"))
export plotOriginalSignal1D, plotZerothLayer1D, plotFirstLayer1DSingleWavelet, gifFirstLayer1D, plotFirstLayer1DAll, plotFirstLayer1D, plotSecondLayer1DOld, plotSecondLayer1DSpecificPath, gifSecondLayer1DSubset, plotSecondLayer1DFixAndVary, plotSecondLayer1D, jointPlot1D
export plotOriginalSignal2D, plotZerothLayer2D, plotFirstLayer2DSingleWavelet, visualizeFirstLayer2D, plotFirstLayer2D, plotFirstLayer2DAll, plotSecondLayer2DSingleWavelet, visualizeSecondLayer2D, plotSecondLayer2D

end