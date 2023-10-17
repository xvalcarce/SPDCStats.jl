module SPDCStats

export param,
	   p_ab,
	   p_a_4b,
	   p_4a_4b,
	   correlator,
	   spdc_correlators

import LinearAlgebra: svdvals

include("./spdc_statistics.jl")

end # module
