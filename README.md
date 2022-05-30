# SPDCStats.jl
----------

Julia package to compute statistics arising from a SPDC source and NPRPDs

```
                                            Spatial mode number:
										  A_k         (1)
                                          |
			 ↗ a_k , a_k⟂  → ( α_x,ϕ_αx ) ▩ ─ A_k⟂    (2)
    SPDC ⇉⇉⇉⇉
			 ↘ b_k , b_k⟂  → ( β_y,ϕ_βy ) ▩ ─ B_k     (3)
                                          |  
										  B_k⟂        (4)  

```

System parameters are:
+ `X`,`Y` : Number of measurement for Alice/Bob
+ `N` : Number of modes (k∈[1,N])
+ `ηA`,`ηB` : Detector efficieny of Alice/Bob
+ `pdcA`,`pdcB` : Probability of dark count (thermal noise) in Alices/Bob's detectors

Circuit parameters are:
+ `g`,`g_` : Squeezing parameters
+ `α_x`,`ϕ_αx` : Alice setting for measurement A_x
+ `β_y`,`ϕ_βy` : Alice setting for measurement B_y

## Usage
--------

This package export 4 functions:


`param` : translate a vector of `Float64` to a `NamedTuple`. Contains both circuit paramaters and system parameters (as keyword arguments) with the exception of `X` and `Y`. Default efficiencies are 1.0, default number of mode is 1 and default dark count probability is 0.0. Note that the two first circuit parameters (squeezing) are stored as `Tg = tanh(g)` and `Tg_= tanh(g_)`

```julia
# Construct the NamedTuple from random parametersin a X,Y=(2,3) scenario with detection efficiency 0.8
> p = param(rand(12),2,3;ηA=.8,ηB=.8)
(Tg = 0.2494803207411917, Tg_ = 0.7424573267165877, α1 = 0.010129492977825283, ϕ_α1 = 0.2957634323161866, α2 = 0.8837522471034687, ϕ_α2 = 0.71578185734515, β1 = 0.48774641481304837, ϕ_β1 = 0.1872266067465822, β2 = 0.1703390158609357, ϕ_β2 = 0.8110150871030339, β3 = 0.4627290061607028, ϕ_β3 = 0.5435762981679317, N = 1.0, pdcA = 0.0, pdcB = 0.0, ηA = 0.8, ηB = 0.8)

```

`p_ab` : Compute the set of probabilities `p(ab|xy)` for all `a,b∈[-1,+1]`. Return a matrix with element `[a,b]` corresponding to `p(ab|xy)`.

```julia
# Probabilites for x=1,y=1
> p_11 = p_ab(p,1,1)
2×2 Matrix{Float64}:
 0.00437022  0.021084
 0.31208     0.662466
```


`correlator` : Compute the correlator `<A_x,A_y>` either from a given Matrix of the form arising from `p_ab`, or directly from the `NamedTuple`.

```julia
# Correlator <A_1 B_1> from the probabilities
> A1B1 = correlator(p_11)
0.33367240161689893
# Correlator <A_1 B_1> from the parameters
> A1B1 = correlator(p,1,1)
0.33367240161689893
```
	   
`spdc_correlators` : Return all the `X*Y` correlators. The output is a matrix with element `[x,y]` equal to `<A_x B_y>`.

```julia
> spdc_correlators(p,2,3)
2×3 Matrix{Float64}:
 0.333672  0.058213  0.307882
 0.286395  0.260977  0.263226
```

## Examples and data
--------------------

Examples are available in the example folder.
+ `chsh.jl` : show how to compute and optimize the CHSH value from the SPDC setup.
+ `chsh.py` : same as `chsh.jl` but as an import in Python.
+ `key_rate_with_noisypreprocessing.jl` : show how to compute and optimize the DIQKD key rate, obtained using noisy preprocessing, from the SPDC setup.

We provide some data
+ `key_rate_noisypp_with_loss.csv` : Contains the optimized key rate as well as the circuit parameters (X=2,Y=3) for different value of efficiency.
+ `key_rate_noisypp_with_loss.jld2` : Same as above but in HDF5.


## Installation
---------------

Install `SPDCStats.jl` using Julia `Pkg` package as follow

```julia
> using Pkg
> Pkg.add(url="https://github.com/xvalcarce/SPDCStats.jl")
```
