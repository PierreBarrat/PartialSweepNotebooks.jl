### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ d9b9b640-61de-11ed-0799-1f05d43e1a62
begin
	using Revise
	using Pkg; Pkg.activate("../../../")
	using Chain
	using Measures
	using Plots
	using PlutoUI
	using SIR_partial_sweep
end

# ╔═╡ ce75b1c6-9f8d-479e-9330-9c1c65e4f83c
begin
	include("$(homedir())/.julia/config/startup.jl")
	Plots.default(;pubfig()...)
end

# ╔═╡ b7ad8cb4-8db9-4519-9ce7-6b4a73c33c48
begin
	prefix = let
		f = @__FILE__
		@chain f begin
			replace(r"#==#.*" => "") # replace weird pluto name
			splitdir
			getindex(2)
			replace(r"\.jl" => "")
			*("_")
		end
	end
	savedir = "../../figures/SIR/"
	md"""
	- Save directory: `savedir` = $savedir
	- prefix: `prefix` = $(prefix)
	"""
end

# ╔═╡ a6c9b10c-3406-4d6f-925e-3f88b3bd5e6b
begin
	N = 2 # two viruses
	md"N=$N variants"
end

# ╔═╡ 241f9436-d49e-425f-b3df-58929f909701
begin
	b = 0.6
	f = 0.8
	K = cross_immunity(N, [(b, f)])
	@info "Cross-immunity" K
end

# ╔═╡ 37121b5f-2746-4d5a-89a4-a72b91a65e33
γs = let
	γvals = [.01, .1, .2]
	@bind γ Slider(γvals, default=γvals[1], show_value=true)
end

# ╔═╡ 69321dc4-ab3e-4ab1-9c1b-072e43434e2e
begin
	α = 3
	δ = 1
	# γ = .01
	# md" $\alpha$=$(α) -  γ=$γ, δ=$δ"
	md" α = $(α), γ = $(γ), δ = $(δ)"
end

# ╔═╡ 347d642f-8c2e-4757-9f5f-a522bca883ef
begin
	# Initialize with only the wild-type, and simulate to reach eq.
	X_init = let
		# SIRState with I=0 for the mutant
		params = SIRParameters(; α, γ, δ, N)
		I0 = [1e-6, 0]
		region = SIRRegion(N, I0; K)
		X = SIRState(parameters = params, regions = [region])
		# Simulate it for long enough
		sol = simulate(X, (0, 10/γ))
		sol(sol.tspan[end])
	end
	# Introducing mutant 
	X_init = set_infected(X_init, N, 1e-6)
	X_init.regions[1]
end

# ╔═╡ 4dc11af3-1eb5-48cb-8a48-faf45be484ae
simulation = simulate(X_init, (0, 10/γ));

# ╔═╡ 6feb5d8e-9e89-4ae2-93bb-fb10706f411c
let
	tvals = range(simulation.tspan..., length=100)
	f = frequency(simulation, tvals, N) 
	f_eq = @chain X_init equilibrium frequency(N) # expected eq. freq. 

	p = plot(
		xlabel = "time", 
		ylabel = "Frequency of mutant",
		title = "α=$α, γ=$γ",
	)
	plot!(tvals, f, label="")
	hline!([f_eq], line=(:dash, :black), label="")
	savefig(savedir * prefix * "frequency_vs_time.png")
	p
end

# ╔═╡ 9eee5d0b-8acc-4df5-8155-f7f855a12587
let
	tvals = range(simulation.tspan..., length=100)
	SIR_wt = Dict(g => simulation[tvals, 1, g, 1] for g in (:S, :I, :R))
	SIR_m = Dict(g => simulation[tvals, 1, g, 2] for g in (:S, :I, :R))

	p = plot(
		xlabel = "time", 
		ylabel = "Fraction of host pop.",
		title = "α=$α, γ=$γ",
		yscale=:log10,
		legend = :bottomright,
		ylim = (1e-4, 1.5),
		yticks = [1e-4, 1e-2, 1]
	)
	for (i, g) in enumerate((:S, :I, :R))
		plot!(tvals, SIR_wt[g], label="$g", color=i)
		plot!(tvals, SIR_m[g], label="", color=i, linestyle=:dash)
	end
	savefig(savedir * prefix * "SIR_vs_time.png")
	p
end

# ╔═╡ Cell order:
# ╠═d9b9b640-61de-11ed-0799-1f05d43e1a62
# ╟─b7ad8cb4-8db9-4519-9ce7-6b4a73c33c48
# ╠═ce75b1c6-9f8d-479e-9330-9c1c65e4f83c
# ╟─a6c9b10c-3406-4d6f-925e-3f88b3bd5e6b
# ╠═69321dc4-ab3e-4ab1-9c1b-072e43434e2e
# ╠═241f9436-d49e-425f-b3df-58929f909701
# ╟─347d642f-8c2e-4757-9f5f-a522bca883ef
# ╠═4dc11af3-1eb5-48cb-8a48-faf45be484ae
# ╟─37121b5f-2746-4d5a-89a4-a72b91a65e33
# ╠═6feb5d8e-9e89-4ae2-93bb-fb10706f411c
# ╠═9eee5d0b-8acc-4df5-8155-f7f855a12587
