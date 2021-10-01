using SPDC
using Optim

function chsh(p::Vector{Float64};N=1,ηA=1.0,ηB=1.0)
	X = 2 # Number of Alice's input
	Y = 2 # Number of Bob's input
	p = param(p,X,Y;N=N,ηA=ηA,ηB=ηB)
	AB = spdc_correlators(p,X,Y)
	score = AB[1,1]+AB[1,2]+AB[2,1]-AB[2,2]
	return score
end

function optimize_chsh(;x0=rand(10),N=1,η=1.0)
	r = optimize(x -> -chsh(x,N=N,ηA=η,ηB=η), x0, Optim.Options(iterations=2000))
	return r
end

function chsh_with_loss()
	scores = []
    η = 1.0
    r = optimize_chsh()
    while -r.minimum ≥ 2.0
        push!(scores,[η,-r.minimum])
		@info scores[end]
        η -= .01
        r = optimize_chsh(η=η,x0=r.minimizer)
	end
	return scores
end
