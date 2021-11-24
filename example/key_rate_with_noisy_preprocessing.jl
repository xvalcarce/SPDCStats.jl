using SPDC
using Optim

function h(z::Union{Float64,Complex{Float64}})
	""" Binary entropy """
	return real(-z*log2(z)-(1-z)*log2(1-z))
end

function H(p_ab::Matrix{Float64};noisypp=0.0)
	""" Compute relative entropy H(A|B).

	p_ab : Matrix with element (a,b) = p(a,b|x,y)
	p    : Noisy-preprocessing parameter
	"""
	noisypp = abs(noisypp % 0.5)
	h_b = 0.0 ; h_ab = 0.0
	dim = size(p_ab)
	for b in 1:dim[2]
		p_b = 0
		for a in 1:dim[1]
			p_b += p_ab[a,b]
			ᵃ = a == 1 ? 2 : 1
			p_ᵃb = p_ab[a,b]*(1-noisypp)+p_ab[ᵃ,b]*noisypp
			h_ab -= p_ᵃb*log2(p_ᵃb)
		end
		h_b -= p_b*log2(p_b)
	end
	H_AB = h_ab-h_b
	return H_AB
end

function chsh(p::NamedTuple)
	""" Compute the CHSH value """
	AB = spdc_correlators(p,2,2)
	score = AB[1,1]+AB[1,2]+AB[2,1]-AB[2,2]
	return score
end

function key_rate(S::Float64,H_AB::Float64;noisypp=1e-8)
	""" Key rate in the Ho protocol (noisy-preprocessing)

	S          : CHSH score
	H_AB       : H(A|B)
	noisypp    : noisy preprocessing parameter

	Note
	----
	This function extend the domain of the key rate formula
	to < 2.0 CHSH score. This allows for smoother parameter optimization.
	"""
	if S ≤ 0.0
		return 0.0
	end
	noisypp = abs(noisypp % 0.5)
	if S == 2*√(2)
		return 1.0
	else
		h0 = S > 2.0 ? h((1+ √((S/2)^2-1))/2) : 1 + h((1 +√((S/2)))/2)
		h1 = S > 2.0 ? h((1+ √(Complex(1-noisypp*(1-noisypp)*(8-S^2))))/2) : h((1+ √(Complex(1-noisypp*(1-noisypp)*(S^2))))/2) 
		I_p = h0-h1
		@debug "CHSH: $S | I_p: $I_p | H_B: $H_AB"
		r = 1 - I_p - H_AB
	end
end

function key_rate(p::Vector{Float64};N=1,η=1.0,noisypp=1e-8,v=1.0)
	""" Key rate in the Ho protocol (noisy-preprocessing)
	
	p : parameters
	"""
	X = 2 # Number of Alice inputs
	Y = 3 # Number of Bob inputs
	p = param(p,X,Y,N=N,ηA=η,ηB=η,v=v)
	S = chsh(p) # CHSH is computed using A_1,A_2,B_1,B_2
	p_13 = p_a_4b(p,1,3) # We use A_1 B_3 to generate the key
	if any(p_13 .< 0 ) # Can occur due to numeric precision
		@warn "Found negative probability in $p_13"
		p_13 = abs.(p_13)
	end
	H_AB = H(p_13,noisypp=noisypp)
	r = key_rate(S,H_AB;noisypp=noisypp)
	return r
end


function optim_key_rate_noisypp(;N=1,η=1.0,x0=rand(12),noisypp=1e-8,v=1.0)
	r = optimize(x -> -key_rate(x[1:end-1];N=N,η=η,noisypp=x[end],v=v),
				 [x0...,noisypp], Optim.Options(iterations=3000))
	return r
end

X0 = [0.8637298055783627, 0.7189379837319168, -0.013784943422925924, 0.24108344648959334, 0.5129371064939435, 0.24013915998799606,1.8147086086681588, 0.24010027842186935, 1.2951511547344126, 0.2402467428421382, -0.013653873756858004, 0.24130485427054246, 1e-15]


function key_rate_noisypp_with_loss(;maxstep=1e-2,minstep=1e-5,threshold=1e-6,x0=X0)
	kr = []
    η = 1.0
	if x0 == []
		r = optim_key_rate_noisypp(v=v)
		# Sanity check, since we know the key_rate > 0.25 for η=1.0
		if -r.minimum > 0.25
			for i in 1:20
				r = optim_key_rate_noisypp(x0=r.minimizer[1:end-1],noisypp=r.minimizer[end],v=v)
			end
		else
			return []
		end
	else
		r = optim_key_rate_noisypp(x0=x0[1:end-1],noisypp=x0[end],v=v)
	end
	while -r.minimum ≥ threshold
        push!(kr,[η,-r.minimum,r.minimizer])
		@info kr[end][1:2]
		step_ = 10^(floor(log10(-r.minimum)))
		step = step_ > maxstep ? maxstep : step_
		step = step < minstep ? minstep : step
		η -= step
		r = optim_key_rate_noisypp(;η=η,x0=r.minimizer[1:end-1],noisypp=r.minimizer[end],v=v)
		r = optim_key_rate_noisypp(;η=η,x0=r.minimizer[1:end-1],noisypp=r.minimizer[end],v=v)
	end
	return kr
end

function key_rate_noisypp_with_loss_visibility(;threshold=1e-6,x0=X0,v=1.0)
	kr = []
    η = 1.0
	if x0 == []
		r = optim_key_rate_noisypp(v=v)
		# Sanity check, since we know the key_rate > 0.25 for η=1.0
		if -r.minimum > 0.25
			for i in 1:20
				r = optim_key_rate_noisypp(x0=r.minimizer[1:end-1],noisypp=r.minimizer[end],v=v)
			end
		else
			return []
		end
	else
		r = optim_key_rate_noisypp(x0=x0[1:end-1],noisypp=x0[end],v=v)
	end
	for η in 1.0:-1e-3:0.841
		if -r.minimum ≥ threshold
			push!(kr,[η,v,-r.minimum])
			@info kr[end][1:3]
			η = round(η,digits=3)
			r = optim_key_rate_noisypp(;η=η,x0=r.minimizer[1:end-1],noisypp=r.minimizer[end],v=v)
			r = optim_key_rate_noisypp(;η=η,x0=r.minimizer[1:end-1],noisypp=r.minimizer[end],v=v)
		else
			push!(kr,[η,v,NaN])
			η = round(η,digits=3)
		end
	end
	return kr
end
