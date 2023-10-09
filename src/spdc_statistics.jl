"""
    SPDC source
	-----------
    
    SPDC source emitting a photon-pair. 

                                                Spatial mode number:
										  A_k         (1)
                                          |
			 ↗ a_k , a_k⟂  → ( α_x,ϕ_αx ) ▩ ─ A_k⟂    (2)
    SPDC ⇉⇉⇉⇉
			 ↘ b_k , b_k⟂  → ( β_y,ϕ_βy ) ▩ ─ B_k     (3)
                                          |  
										  B_k⟂        (4)  


    The photon-pair produces photon in coupled modes (a_k,b_k), then send to (Alice,Bob).
    The amount of spatial, frequency and time mode is limited to N. Thus, k∈[1,N].
    
    
    Parameters of the system:
	-------------------------

	X : Alice's number of input (i.e. settings)
	Y : Bob's number of input (i.e. setings)
    
    SPDC source:
	------------
      g  : squeezing parameter for the coupled modes a_k,b_k⟂
      g_ : squeezing parameter for the coupled modes a_k⟂,b_k
	  N  : number of modes (k ∈ [1,N])
    
     Measurement settings:
	 ---------------------
      α_x  : Rotation of the Bloch sphere on Alice's side by α_x
      β_y  : Rotation of the Bloch sphere on Bob's side by β_y
	  ϕ_αx : Phased introduced between Alice's mode
	  ϕ_βy : Phased introduced between Alice's mode
	  
	 Losses and noises:
	 ------------------
	  ηA   : Alice's detectors efficiency
	  ηB   : Bob's detectors efficiency
	  pdcA : Dark count probability on Alice's side
	  pdcB : Dark count probability on Alice's side
      
    Alice can measure her photon in A,A⟂
    Bob can measure his photon in B,B⟂
"""

function M(p::NamedTuple,R::Vector{Float64})
	""" cf. 10.1103/PhysRevA.91.012107 above Eq. 4 """
	RA0,RA1,RB0,RB1 = R
	M11 = √(RA0*RB0)*(p.Tg*cos(p.α)*sin(p.β)*exp(-p.ϕ_β*im)-p.Tg_*sin(p.α)*exp(-p.ϕ_α*im)*cos(p.β))
	M12 = √(RA0*RB1)*(-p.Tg*cos(p.α)*cos(p.β)-p.Tg_*sin(p.α)*exp(-p.ϕ_α*im)*sin(p.β)*exp(p.ϕ_β*im))
	M21 = √(RA1*RB0)*(p.Tg*sin(p.α)*exp(p.ϕ_α*im)*sin(p.β)*exp(-p.ϕ_β*im)+p.Tg_*cos(p.α)*cos(p.β))
	M22 = √(RA1*RB1)*(-p.Tg*sin(p.α)*exp(p.ϕ_α*im)*cos(p.β)+p.Tg_*cos(p.α)*sin(p.β)*exp(p.ϕ_β*im))
	M = [[M11 M12]
		 [M21 M22]]
	return M
end

function p_nc(p::NamedTuple,i::Int64)
	""" Probability of no click in mode i. """
	R = [1.,1.,1.,1.]
	R[i] = i<=2 ? 1-p.ηA : 1-p.ηB
	𝕄 = M(p,R)
	λ1,λ2 = svd(𝕄).S
	p_dc = i<=2 ? (1-p.pdcA) : (1-p.pdcB)
	p_ψ = (((1-p.Tg^2)*(1-p.Tg_^2))/((1-λ1^2)*(1-λ2^2)))^p.N
	c = (p.Tg^2+p.Tg_^2)
	p_ψ1 = (1/(p.Tg^2+p.Tg_^2))*(1-p.Tg^2)*(1-p.Tg_^2)*sum(𝕄.^2)
	p_𝕀 = .25*sum([R[1]R[3],R[1]R[4],R[2]R[3],R[2]R[4]])
	pr = p_dc*(p_ψ + (c*p.v - c)p_ψ1 + c*(1-p.v)*p_𝕀)
	return real(pr)
end

function p_nc(p::NamedTuple,i::Int64,j::Int64)
	""" Probability of no click in modes i,j. """
	R = [1.,1.,1.,1.]
	R[i] = i<=2 ? 1-p.ηA : 1-p.ηB
	R[j] = j<=2 ? 1-p.ηB : 1-p.ηA
	𝕄 = M(p,R)
	λ1,λ2 = svd(𝕄).S
	p_dc = i<=2 ? (1-p.pdcA) : (1-p.pdcB)
	p_dc *= j<=2 ? 1-p.pdcA : 1-p.pdcB
	p_ψ = (((1-p.Tg^2)*(1-p.Tg_^2))/((1-λ1^2)*(1-λ2^2)))^p.N
	c = (p.Tg^2+p.Tg_^2)
	p_ψ1 = (1/(p.Tg^2+p.Tg_^2))*(1-p.Tg^2)*(1-p.Tg_^2)*sum(𝕄.^2)
	p_𝕀 = .25*sum([R[1]R[3],R[1]R[4],R[2]R[3],R[2]R[4]])
	pr = p_dc*(p_ψ + (c*p.v - c)p_ψ1 + c*(1-p.v)*p_𝕀)
	return real(pr)
end

function p_nc(p::NamedTuple,i::Int64,j::Int64,k::Int64)
	""" Probability of no click in modes i,j,k. """
	R = [1.,1.,1.,1.]
	R[i] = 1-p.ηA
	R[j] = j<=2 ? 1-p.ηA : 1-p.ηB
	R[k] = 1-p.ηB
	𝕄 = M(p,R)
	λ1,λ2 = svd(𝕄).S
	p_dc = j<=2 ? (1-p.pdcA)*(1-p.pdcA)*(1-p.pdcB) : (1-p.pdcA)*(1-p.pdcB)*(1-p.pdcB)
	p_ψ = (((1-p.Tg^2)*(1-p.Tg_^2))/((1-λ1^2)*(1-λ2^2)))^p.N
	c = (p.Tg^2+p.Tg_^2)
	p_ψ1 = (1/(p.Tg^2+p.Tg_^2))*(1-p.Tg^2)*(1-p.Tg_^2)*sum(𝕄.^2)
	p_𝕀 = .25*sum([R[1]R[3],R[1]R[4],R[2]R[3],R[2]R[4]])
	pr = p_dc*(p_ψ + (c*p.v - c)p_ψ1 + c*(1-p.v)*p_𝕀)
	return real(pr)
end

function p_nc(p::NamedTuple,i::Int64,j::Int64,k::Int64,l::Int64)
	""" Probability for all detectors to not click """
	R = [1.0-p.ηA,1.0-p.ηA,1.0-p.ηB,1.0-p.ηB]
	𝕄 = M(p,R)
	λ1,λ2 = svd(𝕄).S
	p_dc = (1-p.pdcA)*(1-p.pdcA)*(1-p.pdcB)*(1-p.pdcB)
	p_ψ = (((1-p.Tg^2)*(1-p.Tg_^2))/((1-λ1^2)*(1-λ2^2)))^p.N
	c = (p.Tg^2+p.Tg_^2)
	p_ψ1 = (1/(p.Tg^2+p.Tg_^2))*(1-p.Tg^2)*(1-p.Tg_^2)*sum(𝕄.^2)
	p_𝕀 = .25*sum([R[1]R[3],R[1]R[4],R[2]R[3],R[2]R[4]])
	pr = p_dc*(p_ψ + (c*p.v - c)p_ψ1 + c*(1-p.v)*p_𝕀)
	return real(pr)
end

function param_xy(p::NamedTuple,x::Int64,y::Int64)
	""" Extract parameter for choice of settings x,y """
	p_name = (:Tg,:Tg_,:α,:ϕ_α,:β,:ϕ_β,:N,:pdcA,:pdcB,:ηA,:ηB,:v)
	alice_settings = [getproperty(p,Symbol(field)) for field in ["α$x","ϕ_α$x"]]
	bob_settings = [getproperty(p,Symbol(field)) for field in ["β$y","ϕ_β$y"]]
	param = NamedTuple{p_name}([p.Tg,p.Tg_,alice_settings...,bob_settings...,p.N,p.pdcA,p.pdcB,p.ηA,p.ηB,p.v])
	return param
end

function p_ab(p::NamedTuple,x::Int64,y::Int64)
	""" Return the statistic p(ab|xy) ∀ a,b ∈ [1,2] 

	Probability are binned according to the rule
	(c, nc)   → 0
	otherwise → 1

	"""
	# Get the corresponding measurement settings parameters
	param = param_xy(p,x,y)

	# Compute the relevant no-click probabilities
	p_ = Dict(["$(modes...)" => p_nc(param,modes...) 
			   for modes in [2,4,[1,2],[3,4],[2,4],[1,2,4],[2,3,4],[1,2,3,4]]])

	# Construct p(±1±1|xy)
	p_mm = p_["24"]-p_["124"]-p_["234"]+p_["1234"]
	p_pm = p_["4"]-p_["34"]-p_["24"]+p_["124"]+p_["234"]-p_["1234"]
	p_mp = p_["2"]-p_["12"]-p_["24"]+p_["234"]+p_["124"]-p_["1234"]
	p_pp = 1-(p_mm+p_mp+p_pm)
	p_ab = [p_mm p_mp
			p_pm p_pp]
	return p_ab
end

function p_a_4b(p::NamedTuple,x::Int64,y::Int64)
	""" Return the satistic p(ab|xy) without binning Bob side.
	i.e. ∀ a∈[1,2], b∈[1,2,3,4].

	Alice:
	0 → (c ,nc)
	1 → (nc nc, nc c, c c)

	Bob:
	0 → (nc nc)
	1 → (nc c)
	2 → (c nc)
	3 → (c c)
	"""
	param = param_xy(p,x,y)

	p_ = Dict(["$(modes...)" => p_nc(param,modes...) 
			   for modes in [2,3,4,[1,2],[2,3],[2,4],[3,4],[1,2,3],[1,2,4],[2,3,4],[1,2,3,4]]])

	ps = zeros(2,4)

	ps[1,1] = p_["234"]-p_["1234"]
	ps[1,2] = p_["23"]-p_["123"]-ps[1,1]
	ps[1,3] = p_["24"]-p_["124"]-ps[1,1]
	ps[1,4] = p_["2"]-p_["12"]-sum(ps)

	ps[2,1] = p_["34"]-ps[1,1]
	ps[2,2] = p_["3"]-p_["34"]-ps[1,2]
	ps[2,3] = p_["4"]-p_["34"]-ps[1,3]
	ps[2,4] = 1-p_["3"]-p_["4"]+p_["1234"]-ps[1,4]
	
	return ps
end

function p_4a_4b(p::NamedTuple,x::Int,y::Int)
	""" Return the raw satistic p(ab|xy).
	i.e. ∀ a∈[1,2,3,4], b∈[1,2,3,4].

	Alice:
	0 → (nc nc)
	1 → (nc c)
	2 → (c nc)
	3 → (c c)

	Bob:
	0 → (nc nc)
	1 → (nc c)
	2 → (c nc)
	3 → (c c)
	"""
	param = param_xy(p,x,y)

	p_ = Dict(["$(modes...)" => p_nc(param,modes...) 
	for modes in [1,2,3,4,[1,2],[1,3],[1,4],[2,3],[2,4],[3,4],[1,2,3],[1,2,4],[1,3,4],[2,3,4],[1,2,3,4]]])
	
	ps = zeros(4,4) 
	
	ps[1,1] = p_["1234"]
	ps[1,2] = p_["123"]-p_["1234"]
	ps[1,3] = p_["124"]-p_["1234"]
	ps[1,4] = p_["12"]-ps[1,1]-ps[1,2]-ps[1,3]

	ps[2,1] = p_["134"]-p_["1234"]
	ps[2,2] = p_["13"]-p_["123"]-ps[2,1]
	ps[2,3] = p_["14"]-p_["124"]-ps[2,1]
	ps[2,4] = p_["1"]-p_["12"]-ps[2,1]-ps[2,2]-ps[2,3]

	ps[3,1] = p_["234"]-p_["1234"]
	ps[3,2] = p_["23"]-p_["123"]-ps[3,1]
	ps[3,3] = p_["24"]-p_["124"]-ps[3,1]
	ps[3,4] = p_["2"]-p_["12"]-ps[3,1]-ps[3,2]-ps[3,3]

	ps[4,1] = p_["34"]-ps[1,1]-ps[2,1]-ps[3,1]
	ps[4,2] = p_["3"]-p_["13"]-ps[3,1]-ps[3,2]-ps[4,1]
	ps[4,3] = p_["4"]-p_["14"]-ps[3,1]-ps[3,3]-ps[4,1]
	ps[4,4] = 1-sum(ps)

	return ps
end

function correlator(p::NamedTuple,x::Int64,y::Int64)
	""" Return the correlator <A_x B_y> """
	return correlator(p_ab(p,x,y))
end

function correlator(p::Matrix{Float64})
	""" Construct correlator c = p(a=b) - p(a≢b) """
	corr = p[1,1]+p[2,2]-p[1,2]-p[2,1]
	return corr
end

function spdc_correlators(p::NamedTuple,X::Int64,Y::Int64)
	""" Return all correlators <A_x B_y> for x∈[1,X],y∈[1,Y] """
	corr = [correlator(p_ab(p,x,y)) for x in 1:X, y in 1:Y]
	return corr
end

function param(p::Vector{Float64},X::Int64,Y::Int64;N=1,pdcA=0.,pdcB=0.,ηA=1.0,ηB=1.0,v=1.0)
	""" Construct a NamedTuple containgin all the information of the SPDC system. """
	length_param =  2 + 2*X + 2*Y
	@assert length(p) == length_param "Parameters need to be of size $length_param, but got size $(length(p))"
	param = copy(p)
	param[1] = tanh(param[1])
	param[2] = tanh(param[2])
	alice_settings = [Symbol(a*string(x)) for x in 1:X for a in ["α","ϕ_α"]]
	bob_settings = [Symbol(b*string(y)) for y in 1:Y for b in ["β","ϕ_β"]]
	p_name = (:Tg,:Tg_,alice_settings...,bob_settings...,:N,:pdcA,:pdcB,:ηA,:ηB,:v)
	p = NamedTuple{p_name}([param...,N,pdcA,pdcB,ηA,ηB,v])
	return p
end
