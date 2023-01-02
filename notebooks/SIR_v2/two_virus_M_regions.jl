### A Pluto.jl notebook ###
# v0.19.19

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

# ╔═╡ 4982e04e-6501-11ed-132e-ad3e395233da
begin
	using Revise
	using Pkg; Pkg.activate()
	using Chain
	using Measures
	using Parameters
	using Plots
	using PlutoUI
	using PartialSweepSIR
end

# ╔═╡ 0bca5e5c-18cf-49f4-a1d2-bfa579b8f81e
begin
	include("$(homedir())/.julia/config/startup.jl")
	Plots.default(;pubfig()...)
end

# ╔═╡ e6e61028-1b62-4369-a7cd-ebd348f30592
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

# ╔═╡ db663438-8523-47ac-a783-88ff59a7ab1d
md"# Setup"

# ╔═╡ 287520e8-97a3-485b-9634-fd2db7a2f4e6
begin
	_Ms = @bind M Slider([2,3,5,10,20], show_value=true, default=2)
	_Cs = @bind c Slider(vcat(0, 10. .^ (-3:.25:-1)), show_value=true, default=.01)
	Ms = md"M = $(_Ms)"
	Cs = md"c = $(_Cs)"
end;

# ╔═╡ f431ba40-9742-4898-92f0-33c3b8e12a1c
params = let
	N = 2
	# M = 10
	α = 3
	γ = .0025
	# c = 0.1
	PSS.Parameters(; N, M, α, γ, c)
end;

# ╔═╡ e38f872f-f9ec-45cf-953c-42a18aa0bab3
@unpack α, γ, δ = params

# ╔═╡ c2a159f4-cd25-4525-a270-2b3065a0dc0e
T = 3/params.γ # simulation time

# ╔═╡ 78b03821-37d4-4a9d-8ff0-349e12f66686
begin
	K0 = let
		b = 1
		f = 1
		[1 b; f 1]
	end
	K_special = let
		b = .5
		f = .4
		[1 b; f 1]
	end
end

# ╔═╡ a62f2142-7d76-4a15-9566-54115aa98e08
regions = let
	S0 = .4
	I0 = 1e-6
	C0 = 1e-9
	v1 = PSS.Virus(; S=S0, I=I0, C=C0) # wt virus
	v2 = PSS.Virus(; S=S0, I=0, C=0.) # mutant: initially not here
	
	R = [PSS.Region(; viruses=[v1,v2], K=K0) for m in 1:params.M]
	r1 = PSS.Region(; viruses=[v1,v2], K=K_special)
	R[1] = r1 # first region has the special cross-immunity
	R
end

# ╔═╡ fd0d5424-074d-4fee-a41c-a80e61591348
state_init ,sol_init = let
	state = PSS.SIRState(; regions, parameters=params)
	sol = PSS.simulate(state, (0, T))
	PSS.set_infected(sol(T), 2, 1e-6), sol
end;

# ╔═╡ c72f0ae1-e9df-4b88-9629-f471ba9f6e24
# ╠═╡ disabled = true
#=╠═╡
state_eq = PSS.equilibrium(state_init);
  ╠═╡ =#

# ╔═╡ 62c7a1db-42d0-463f-9f3d-954d36f1b1d2
sol = PSS.simulate(state_init, (0, T/2));

# ╔═╡ bb4a3499-958f-4b34-b69c-de2874ae8546
md"# Figures"

# ╔═╡ 88b4146d-d8bb-4793-bec2-76a4cbc4e320
Ms

# ╔═╡ d17bc99f-f5cb-49c5-bfe1-e9219b64d19d
Cs

# ╔═╡ 87ffd6f9-3b88-46ba-a2a1-a67c67043aa0
let
	# freq. plot
	tvals = range(sol.tspan..., length=100)
	f_R1 = PSS.frequency(sol, tvals, 1, 2)
	f_R2 = PSS.frequency(sol, tvals, 2, 2)
	f_R = PSS.frequency(sol, tvals, 2)

	lw = 4

	p = plot(legend = :bottomright)
	plot!(tvals, f_R1, label="R1 (special)", linewidth=lw)
	plot!(tvals, f_R2, label="Other regions", linewidth=lw)
	plot!(tvals, f_R, label="Overall", linewidth=lw)

	savefig(p, savedir * prefix * "traj_frequency_mutant.png")
	p
end

# ╔═╡ add14358-6d60-47a5-876f-8e0651781c57
pS = let
	g = :S
	tvals = range(sol.tspan..., length=400)
	lw = 4
	
	X_R1 = sol[tvals, 1, 2, g]
	X_R2 = sol[tvals, 2, 2, g]
	# X_R = sol[tvals, 2, g]
	
	p = plot(
		ylabel="Susceptibles", 
		xlabel="time",
		legend=:topright,
		yticks = [.2, .3, .4, .5, .6],
	)
	plot!(tvals, X_R1, label="Special region", linewidth = lw)
	plot!(tvals, X_R2, label="Other regions", linewidth = lw)
	hline!([δ/α], line=(lw+2, :black, 0.3), label="")

	p
end

# ╔═╡ 2e5f1f19-d5f2-4e54-8337-49db178b13c9
pI = let
	g = :I
	tvals = range(sol.tspan..., length=400)
	lw = 4
	
	X_R1 = sol[tvals, 1, 2, g]
	X_R2 = sol[tvals, 2, 2, g]
	# X_R = sol[tvals, 2, g]
	
	p = plot(
		ylabel="Infected", 
		xlabel="time",
		legend=:bottomright,
		yscale=:log10,
	)
	plot!(tvals, X_R1, label="Special region", linewidth = lw)
	plot!(tvals, X_R2, label="Other regions", linewidth = lw)
	# hline!([δ/α], line=(lw+2, :black, 0.3), label="")

	p
end

# ╔═╡ 48be9f2c-97dc-43ad-bf4e-4a0db37acf59
let
	p = plot(pS, pI, layout=(1,2))
	savefig(p, savedir * prefix * "traj_S_I.png")
	p
end

# ╔═╡ Cell order:
# ╠═4982e04e-6501-11ed-132e-ad3e395233da
# ╠═e6e61028-1b62-4369-a7cd-ebd348f30592
# ╠═0bca5e5c-18cf-49f4-a1d2-bfa579b8f81e
# ╟─db663438-8523-47ac-a783-88ff59a7ab1d
# ╠═287520e8-97a3-485b-9634-fd2db7a2f4e6
# ╠═f431ba40-9742-4898-92f0-33c3b8e12a1c
# ╠═e38f872f-f9ec-45cf-953c-42a18aa0bab3
# ╠═c2a159f4-cd25-4525-a270-2b3065a0dc0e
# ╠═78b03821-37d4-4a9d-8ff0-349e12f66686
# ╠═a62f2142-7d76-4a15-9566-54115aa98e08
# ╠═fd0d5424-074d-4fee-a41c-a80e61591348
# ╠═c72f0ae1-e9df-4b88-9629-f471ba9f6e24
# ╠═62c7a1db-42d0-463f-9f3d-954d36f1b1d2
# ╟─bb4a3499-958f-4b34-b69c-de2874ae8546
# ╠═88b4146d-d8bb-4793-bec2-76a4cbc4e320
# ╠═d17bc99f-f5cb-49c5-bfe1-e9219b64d19d
# ╠═87ffd6f9-3b88-46ba-a2a1-a67c67043aa0
# ╠═add14358-6d60-47a5-876f-8e0651781c57
# ╠═2e5f1f19-d5f2-4e54-8337-49db178b13c9
# ╠═48be9f2c-97dc-43ad-bf4e-4a0db37acf59
