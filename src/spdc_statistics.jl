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

function M(p::NamedTuple,R::Array{Float64,1})
	RA0,RA1,RB0,RB1 = R
	M11 = √(RA0*RB0)*(p.Tg*cos(p.α)*sin(p.β)*exp(-p.ϕ_β*im)-p.Tg_*sin(p.α)*exp(-p.ϕ_α*im)*cos(p.β))
	M12 = √(RA0*RB1)*(-p.Tg*cos(p.α)*cos(p.β)-p.Tg_*sin(p.α)*exp(-p.ϕ_α*im)*sin(p.β)*exp(p.ϕ_β*im))
	M21 = √(RA1*RB0)*(p.Tg*sin(p.α)*exp(p.ϕ_α*im)*sin(p.β)*exp(-p.ϕ_β*im)+p.Tg_*cos(p.α)*cos(p.β))
	M22 = √(RA1*RB1)*(-p.Tg*sin(p.α)*exp(p.ϕ_α*im)*cos(p.β)+p.Tg_*cos(p.α)*sin(p.β)*exp(p.ϕ_β*im))
	M = [[M11 M12]
		 [M21 M22]]
	F = svd(M)
	λ1,λ2 = F.S
	return λ1,λ2
end

function p_nc(p::NamedTuple,i::Int64)
	""" Probability of no click in mode i. """
	R = [1.,1.,1.,1.]
	R[i] = i<=2 ? 1-p.ηA : 1-p.ηB
	λ1,λ2 = M(p,R)
	fact = i<=2 ? (1-p.pdcA) : (1-p.pdcB)
	pr = fact*(((1-p.Tg^2)*(1-p.Tg_^2))/((1-λ1^2)*(1-λ2^2)))^p.N
	return pr
end

function p_nc(p::NamedTuple,i::Int64,j::Int64)
	""" Probability of no click in modes i,j. """
	R = [1.,1.,1.,1.]
	R[i] = i<=2 ? 1-p.ηA : 1-p.ηB
	R[j] = j<=2 ? 1-p.ηB : 1-p.ηA
	λ1,λ2 = M(p,R)
	fact = i<=2 ? 1-p.pdcA : 1-p.pdcB
	fact *= j<=2 ? 1-p.pdcA : 1-p.pdcB
	pr = fact*(((1-p.Tg^2)*(1-p.Tg_^2))/((1-λ1^2)*(1-λ2^2)))^p.N
	return pr
end

function p_nc(p::NamedTuple,i::Int64,j::Int64,k::Int64)
	""" Probability of no click in modes i,j,k. """
	R = [1.,1.,1.,1.]
	R[i] = 1-p.ηA
	R[j] = j<=2 ? 1-p.ηA : 1-p.ηB
	R[k] = 1-p.ηB
	λ1,λ2 = M(p,R)
	fact = j<=2 ? (1-p.pdcA)*(1-p.pdcA)*(1-p.pdcB) : (1-p.pdcA)*(1-p.pdcB)*(1-p.pdcB)
	pr = fact*(((1-p.Tg^2)*(1-p.Tg_^2))/((1-λ1^2)*(1-λ2^2)))^p.N
	return pr
end

function p_nc(p::NamedTuple,i::Int64,j::Int64,k::Int64,l::Int64)
	""" Probability for all detectors to not click """
	R = [1.0-p.ηA,1.0-p.ηA,1.0-p.ηB,1.0-p.ηB]
	λ1,λ2 = M(p,R)
	fact = (1-p.pdcA)*(1-p.pdcA)*(1-p.pdcB)*(1-p.pdcB)
	pr = fact*(((1-p.Tg^2)*(1-p.Tg_^2))/((1-λ1^2)*(1-λ2^2)))^p.N
	return pr
end

function param_xy(p::NamedTuple,x::Int64,y::Int64)
	p_name = (:Tg,:Tg_,:α,:ϕ_α,:β,:ϕ_β,:N,:pdcA,:pdcB,:ηA,:ηB)
	alice_settings = [getproperty(p,Symbol(field)) for field in ["α$x","ϕ_α$x"]]
	bob_settings = [getproperty(p,Symbol(field)) for field in ["β$y","ϕ_β$y"]]
	param = NamedTuple{p_name}([p.Tg,p.Tg_,alice_settings...,bob_settings...,p.N,p.pdcA,p.pdcB,p.ηA,p.ηB])
	return param
end

function p_ab(p::NamedTuple,x::Int64,y::Int64)
	""" Return the statistic p(ab|xy) for a,b in [1,2] 

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

function correlator(p::Matrix{Float64})
	""" Construct correlator c = p(a=b) - p(a≢b) """
	corr = p[1,1]+p[2,2]-p[1,2]-p[2,1]
	return corr
end

function spdc_correlators(p::Vector{Float64},X::Int64,Y::Int64)
	p = param(p,X,Y)
	corr = [correlator(p_ab(p,x,y)) for x in 1:X, y in 1:Y]
	return corr
end

function param(p::Vector{Float64},X::Int64,Y::Int64;N=1,pdcA=0.,pdcB=0.,ηA=1.0,ηB=1.0)
	length_param =  2 + 2*X + 2*Y
	@assert length(p) == length_param "Parameters need to be of size $length_param. Passed size: $(length(p))"
	param = copy(p)
	param[1] = tanh(param[1])
	param[2] = tanh(param[2])
	alice_settings = [Symbol(a*string(x)) for x in 1:X for a in ["α","ϕ_α"]]
	bob_settings = [Symbol(b*string(y)) for y in 1:Y for b in ["β","ϕ_β"]]
	p_name = (:Tg,:Tg_,alice_settings...,bob_settings...,:N,:pdcA,:pdcB,:ηA,:ηB)
	p = NamedTuple{p_name}([param...,N,pdcA,pdcB,ηA,ηB])
	return p
end

function chsh(p::NamedTuple)
	score = correlator(p_ab(p,1,1))+correlator(p_ab(p,1,2))+correlator(p_ab(p,2,1))-correlator(p_ab(p,2,2))
	return score
end

chsh(p::Vector{Float64};N=1,η=1.) =	chsh(param(p,2,2;N=N,ηA=η,ηB=η))
