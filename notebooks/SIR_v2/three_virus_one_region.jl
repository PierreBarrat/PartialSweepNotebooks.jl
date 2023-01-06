### A Pluto.jl notebook ###
# v0.19.19

using Markdown
using InteractiveUtils

# ╔═╡ 4982e04e-6501-11ed-132e-ad3e395233da
begin
	using Revise
	using Pkg; Pkg.activate()
	using Chain
	using Measures
	using Parameters
	using Plots
	using PartialSweepSIR
end

# ╔═╡ 3b91ddde-982e-4570-bba8-39d416f2f4ac
begin
	include("$(homedir())/.julia/config/startup.jl")
	Plots.default(;pubfig()...)
end

# ╔═╡ 403e857b-6532-42fa-bbcc-caecf3508ae9
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

# ╔═╡ 50e61003-b631-484c-a63b-486ef81695e9
md"# Parameters, setup"

# ╔═╡ f431ba40-9742-4898-92f0-33c3b8e12a1c
params = let
	N = 3
	M = 1
	α = 3
	γ = .01
	PSS.Parameters(; N, M, α, γ)
end;

# ╔═╡ ef61baf5-b7a4-4997-9e54-fec7d895317a
@unpack α, γ, δ = params;

# ╔═╡ 78b03821-37d4-4a9d-8ff0-349e12f66686
K = let
	b2 = .6
	f2 = .8

	b3 = .8
	f3 = .6
	[
		1 b2 b3; 
		f2 1 b3;
		f3 f3 1
	]
end

# ╔═╡ a62f2142-7d76-4a15-9566-54115aa98e08
region = let
	S0 = .4
	I0 = 1e-6
	C0 = 1e-9
	R0 = 1 - S0 - I0 - C0
	v1 = PSS.Virus(; S=S0, I=I0, C=C0, R=R0)
	v2 = PSS.Virus(; S=1, I=0, C=0., R=0.)
	v3 = PSS.Virus(; S=1, I=0, C=0., R=0.)
	PSS.Region(; viruses=[v1, v2, v3], K)
end;

# ╔═╡ ccc3dffa-04c1-443f-9af6-725d9efbbe75
md"# Simulation"

# ╔═╡ 272fb236-08b4-4628-843c-61d3147de5ff
begin
	teq = 5/γ # equilibration time
	T2 = 1/γ # time to introduce v2
	T3 = T2 + teq # time to introduce v3
end;

# ╔═╡ fd0d5424-074d-4fee-a41c-a80e61591348
state_1 = let
	state = PSS.SIRState(; regions=[region], parameters=params)
	sol = PSS.simulate(state, (0, teq))
	sol(teq)
end;

# ╔═╡ 2cb64cc2-1c94-4e54-af12-6d799c6dddd7
state_2, sol_1 = let
	sol = PSS.simulate(state_1, (0, T2))
	state = PSS.set_infected(sol(T2), 2, 1e-6)
	state, sol
end;

# ╔═╡ db2abde6-6d06-497e-a3b0-3df626e90bfe
state_3, sol_2 = let
	sol = PSS.simulate(state_2, (T2, T3))
	state = PSS.set_infected(sol(T3), 3, 1e-6)
	state, sol
end;

# ╔═╡ c2439ce7-226c-4a11-8589-6e67b0d03f40
sol_3 = let
	PSS.simulate(state_3, (T3, T3+teq))
end;

# ╔═╡ 4060bfc5-d9e4-4bf3-a26d-3b8dcb62490f
sols = [sol_1, sol_2, sol_3];

# ╔═╡ 04845384-01c0-4ef2-8d60-2a76d7e70731
md"# Plots"

# ╔═╡ 4a0f195e-7bf9-4f98-a98b-f5843a1628ef
md"## S and I"

# ╔═╡ 6ac7a34a-1104-47c6-a42e-b072d256b708
let
	x = 2
	tvals = range(sols[x].tspan..., length=50)

	a = 2
	X = sols[x][tvals, 1, a, :I]
	plot(tvals, X)
end


# ╔═╡ 2101d558-f51b-4dbb-970b-529f00f3d3c3
v(a) = if a == 1 
	"w.t."
elseif a == 2
	"Mutant 1"
elseif a == 3
	"Mutant 2"
end

# ╔═╡ bd013ae4-900a-46d4-a40e-fd3c8a0d4cbf
pI = let
	ts = [range(s.tspan..., length=100) for s in sols]

	g = :I
	p = plot(
		# yscale=:log10, 
		ylim = (1e-5, 1e-2), 
		xlabel="time", 
		ylabel="Infected",
		legend=:topright,
	)
	for a in 1:3
		for x in 1:3
			X = sols[x][ts[x], 1, a, g]
			plot!(ts[x], X, color = a, label= (x == 1 ? "$(v(a))" : ""))
		end
	end
	p
end

# ╔═╡ cba5b364-4620-4341-bc06-370e8a38f12a
md"## Frequency"

# ╔═╡ 030863e6-324f-4b30-8f9f-a9b527d83569
md"I assume that mutant 2 comes from the background of mutant 1"

# ╔═╡ c5ef295c-0eab-492b-bfd2-59fa1c0f95c8
f1, f2, f3 = let
	ts = [range(s.tspan..., length=100) for s in sols]
	I = map(1:3) do a 
		I = []
		for x in 1:3
			X = sols[x][ts[x], 1, a, :I]
			append!(I, X)
		end
		I
	end
	I
	f = [I[a] ./ sum(I) for a in 1:3]
	f
end

# ╔═╡ eecf1a0b-6503-4a73-a85b-2168357a0eda
let
	ts = vcat([range(s.tspan..., length=100) for s in sols]...)
	p = plot(
		xlabel = "time",
		ylabel = "Frequency",
		legend = nothing,
		size = (1500,1200),
		xlim = extrema(ts),
		ylim=(0,1),
	)
	plot!(ts, f2 + f3, label="w.t.", fillrange=1)
	plot!(ts, f3, label="Mutant 1", fillrange=f2+f3)
	plot!(ts, zeros(length(f3)), fillrange=f3, label="Mutant 2")

	savefig(p, savedir * prefix * "frequency_mutants_area.png")
	p
end

# ╔═╡ Cell order:
# ╠═4982e04e-6501-11ed-132e-ad3e395233da
# ╠═3b91ddde-982e-4570-bba8-39d416f2f4ac
# ╟─403e857b-6532-42fa-bbcc-caecf3508ae9
# ╟─50e61003-b631-484c-a63b-486ef81695e9
# ╠═f431ba40-9742-4898-92f0-33c3b8e12a1c
# ╠═ef61baf5-b7a4-4997-9e54-fec7d895317a
# ╠═78b03821-37d4-4a9d-8ff0-349e12f66686
# ╠═a62f2142-7d76-4a15-9566-54115aa98e08
# ╟─ccc3dffa-04c1-443f-9af6-725d9efbbe75
# ╠═272fb236-08b4-4628-843c-61d3147de5ff
# ╠═fd0d5424-074d-4fee-a41c-a80e61591348
# ╠═2cb64cc2-1c94-4e54-af12-6d799c6dddd7
# ╠═db2abde6-6d06-497e-a3b0-3df626e90bfe
# ╠═c2439ce7-226c-4a11-8589-6e67b0d03f40
# ╠═4060bfc5-d9e4-4bf3-a26d-3b8dcb62490f
# ╟─04845384-01c0-4ef2-8d60-2a76d7e70731
# ╟─4a0f195e-7bf9-4f98-a98b-f5843a1628ef
# ╠═6ac7a34a-1104-47c6-a42e-b072d256b708
# ╠═2101d558-f51b-4dbb-970b-529f00f3d3c3
# ╠═bd013ae4-900a-46d4-a40e-fd3c8a0d4cbf
# ╟─cba5b364-4620-4341-bc06-370e8a38f12a
# ╠═030863e6-324f-4b30-8f9f-a9b527d83569
# ╠═c5ef295c-0eab-492b-bfd2-59fa1c0f95c8
# ╠═eecf1a0b-6503-4a73-a85b-2168357a0eda
