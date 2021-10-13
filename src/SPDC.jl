module SPDC

export param,
	   p_ab,
	   p_a_4b,
	   correlator,
	   spdc_correlators

import LinearAlgebra: svd

include("./spdc_statistics.jl")

end # module
