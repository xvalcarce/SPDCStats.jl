module SPDC

export param,
	   p_ab,
	   correlator,
	   spdc_correlators

import LinearAlgebra: svd

include("./spdc_statistics.jl")

end # module
