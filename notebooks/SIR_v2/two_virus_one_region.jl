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

# ╔═╡ bac23cd9-57af-441f-ac7d-eaa4c95af133
sf = false

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
	N = 2
	M = 1
	α = 3
	γ = .005
	PSS.Parameters(; N, M, α, γ)
end

# ╔═╡ ef61baf5-b7a4-4997-9e54-fec7d895317a
@unpack α, γ, δ = params

# ╔═╡ 78b03821-37d4-4a9d-8ff0-349e12f66686
K = begin
	b = 0
	f = .5
	[1 b; f 1]
end

# ╔═╡ 272fb236-08b4-4628-843c-61d3147de5ff
begin
	teq = 5/γ
	tsim = 5/γ
	t0 = tsim/5
end;

# ╔═╡ a62f2142-7d76-4a15-9566-54115aa98e08
region = let
	S0 = .4
	I0 = 1e-6
	C0 = 1e-9
	R0 = 1 - S0 - I0 - C0
	v1 = PSS.Virus(; S=S0, I=I0, C=C0, R=R0)
	v2 = PSS.Virus(; S=1, I=0, C=0., R=0.)
	PSS.Region(; viruses=[v1,v2], K)
end

# ╔═╡ fd0d5424-074d-4fee-a41c-a80e61591348
state_init, sol_init = let
	state = PSS.SIRState(; regions=[region], parameters=params)
	sol = PSS.simulate(state, (0, teq))
	sol_init = PSS.simulate(sol(teq), (-t0, 0)) # simulating without mut for show
	state_init = PSS.set_infected(sol(teq), 2, 1e-6) # infect with mutant
	(state_init, sol_init)
end;

# ╔═╡ ccc3dffa-04c1-443f-9af6-725d9efbbe75
md"# Simulation"

# ╔═╡ c72f0ae1-e9df-4b88-9629-f471ba9f6e24
state_eq = PSS.equilibrium(state_init);

# ╔═╡ ccc623a7-7a31-4bcf-9f09-9b7a97f8bf02
sol = PSS.simulate(state_init, (0, tsim));

# ╔═╡ 04845384-01c0-4ef2-8d60-2a76d7e70731
md"# Plots"

# ╔═╡ 4a0f195e-7bf9-4f98-a98b-f5843a1628ef
md"## S and I"

# ╔═╡ 489692aa-9449-43f1-8809-d8b979d11ed0
pS, pI = let
	tpast = range(sol_init.tspan..., length=50)
	tvals = range(sol.tspan..., length=200)
	lw = 6
		v(a) = a == 1 ? "w.t." : "Mutant"
	
	# I 
	pI = plot(
		# yscale=:log10, 
		ylim = (0, 5e-3), 
		xlabel="time", 
		ylabel="",
		title="Infected",
		legend=:topright,
		yticks = ([2e-3, 4e-3], ["2e-3","4e-3"]),
	)
	g = :I
	for a in 1:2
		X = sol_init[tpast, 1, a, g]
		plot!(tpast, X, color=a, label="$(v(a))", line=(lw))
		X = sol[tvals, 1, a, g]
		plot!(tvals, X, color=a, label="", line=(lw))
	end
	vline!(pI, [0.], line=(2, :dashdot, :black), label="")
	
	# S 
	pS = plot(
		ylabel="",
		title="Susceptibles",
		xlabel="time",
		legend=:topright,
		yticks = [0.3,0.33,.34,.35,.38,.4]
	)
	g = :S
	for a in 1:2
		X = sol_init[tpast, 1, a, g]
		plot!(tpast, X, color=a, label="$(v(a))", line=(:dash, lw))
		X = sol[tvals, 1, a, g]
		plot!(tvals, X, color=a, label="", line=(:dash, lw))
	end
	vline!(pS, [0.], line=(2, :dashdot, :black), label="")
	# hline!([δ/α], line=(1.5*lw, :black), alpha=0.4, label="")
	
	(pS, pI)
end;

# ╔═╡ bd381c78-74ff-439b-87c0-b46f0ea8467d
1/(1+f*(1-α/δ))

# ╔═╡ 0006792d-5632-46d2-9fd4-49719f27aad6
let
	p = plot(pS, pI, layout=(1,2), size = (1800,900), margin=5mm, bottom_margin=12mm)
	sf && savefig(p, savedir * prefix * "traj_S_I.png")
	p
end

# ╔═╡ 75677613-21f4-4165-9257-dacd7febf42d
savedir

# ╔═╡ e2252084-7070-42a2-a5dd-695e6d973543
K

# ╔═╡ cba5b364-4620-4341-bc06-370e8a38f12a
md"## Frequency"

# ╔═╡ 5c3787e1-a188-410b-8a35-389c350ea9f6
sol_init[[sol_init.tspan...], 1, :wt, :I]

# ╔═╡ 3187cc3c-7b47-4391-ac00-7b289ed0cd51
savedir

# ╔═╡ ee10c49a-f996-4bee-8dc3-84e74ef6853a
md"## Checking equilibrium values"

# ╔═╡ 5ef29b01-7313-4558-8380-63e7f8101265
let
	# I
	tvals = range(sol.tspan..., length=100)

	g = :I
	p = plot(yscale=:log10, title="$g")
	for a in 1:2
		X = sol[tvals, 1, a, g]
		plot!(tvals, X, color=a, label="")
		hline!([state_eq[1, a, g]], label="", color=a, linestyle=:dash)
	end
	hline!([0.1 / 50], line=(:black, :dash), label="flu prevalence")
	p
end

# ╔═╡ 2599199d-0542-4556-a2a3-ad8a155d8b6f
inv(K)

# ╔═╡ 9bd09e9b-3c09-4728-94b1-508716e1a55a
let
	# C
	tvals = range(sol.tspan..., length=100)

	g = :C
	p = plot(yscale=:log10, title="$g", ylim=(1e-9, 1))
	for a in 1:2
		X = sol[tvals, 1, a, g]
		plot!(tvals, X, color=a, label="")
		hline!([state_eq[1, a, g]], label="", color=a, linestyle=:dash)
	end
	p
end

# ╔═╡ ec1076d9-3cfc-420e-a985-bc0575c99be1
md"### Testing initial growth rate"

# ╔═╡ a51806b9-392c-47f7-9b54-46dfaf515338
let
	tpast = range(sol_init.tspan..., length=50)
	tvals = range(sol.tspan..., length=200)
	lw = 6
		v(a) = a == 1 ? "w.t." : "Mutant"
	
	# I 
	pI = plot(
		# yscale=:log10, 
		ylim = (1e-5, 3e-3), 
		xlabel="time", 
		ylabel="",
		title="Infected",
		legend=:bottomright,
		yticks = ([1e-4, 1e-3], ["1e-4","1e-3"]),
		yscale=:log10,
	)

	X = sol_init[tpast, 1, :m, :I]
	plot!(tpast, X, color=2, label="Mutant", line=(lw))
	X = sol[tvals, 1, :m, :I]
	plot!(tvals, X, color=2, label="", line=(lw))

	R0 = α/δ
	Rm = δ * (1-f)*(R0-1) / (1 + f*(R0-1))
	Im_0 = sol[[0],1,:m,:I][1]
	plot!(tvals, Im_0 * exp.(Rm * tvals), label="", line=(:black, :dashdot))
end

# ╔═╡ 321548a2-8223-4826-9475-cca8a6568d92
function logistic_growth(s, I0, tvals)
	return exp.(s*tvals) ./ (1/I0 - 1 .+ exp.(s*tvals))
end

# ╔═╡ 87ffd6f9-3b88-46ba-a2a1-a67c67043aa0
p_freq = let
	lw = 4

	tpast = range(sol_init.tspan..., length=50)
	tvals = range(sol.tspan..., length=100)
	I1 = sol[tvals, 1, 1, :I]
	I2 = sol[tvals, 1, 2, :I]
	freq = I2 ./ (I1 + I2)

	p = plot(xlabel="Time", title="Frequency of mutant")
	plot!(tpast, 0*similar(tpast), color=1, label="")
	plot!(tvals, freq, label="", line=(lw), color=1)
	# hline!([(1-f)/(2-b-f)], label="", line=(:black, lw+2, 0.3)) # Eq. hline
	# vline!(p, [0.], line=(2, :dashdot, :black), label="") # t=0 vline

	# Comparison with logistic growth
	R0 = α/δ
	Rm = δ * (1-f)*(R0-1) / (1 + f*(R0-1))
	Im_0 = sol[[0],1,:m,:I][1]
	Iwt_0 = sol_init[tpast[end:end], 1, :wt, :I][1]
	plot!(
		tvals, logistic_growth(Rm, Im_0 / (Im_0 + Iwt_0), tvals),
		label="", line=(:black, :dash, 2, 0.7)
	)
	
	sf && savefig(p, savedir * prefix * "traj_frequency_mutant.png")

	p
end

# ╔═╡ f831c365-1fdb-4286-98f1-214379b904d6
let
	p = plot(pS, pI, p_freq, layout=(1,3), size = (2500,900), margin=5mm, bottom_margin=20mm)
	sf && savefig(p, savedir * prefix * "traj_S_I_freq.png")
	sf && savefig(p, "/home/pierrebc/Documents/BaleLabo/Notes/ExpiringFitness/figures/two_strains_partial_sweep.png")
	p
end

# ╔═╡ 5c9fe291-962e-4cd4-bc8d-75ff7c632cb4
md"# Tests"

# ╔═╡ 47723643-4c69-4f68-9b03-c6464fa65f4b
let 
	tvals = range(sol.tspan..., length=50)
	p = plot()
	plot!(tvals, logistic_growth(3e-2, 1e-6, tvals), yscale=:log10)
end

# ╔═╡ Cell order:
# ╠═4982e04e-6501-11ed-132e-ad3e395233da
# ╠═bac23cd9-57af-441f-ac7d-eaa4c95af133
# ╠═3b91ddde-982e-4570-bba8-39d416f2f4ac
# ╠═403e857b-6532-42fa-bbcc-caecf3508ae9
# ╟─50e61003-b631-484c-a63b-486ef81695e9
# ╠═f431ba40-9742-4898-92f0-33c3b8e12a1c
# ╠═ef61baf5-b7a4-4997-9e54-fec7d895317a
# ╠═78b03821-37d4-4a9d-8ff0-349e12f66686
# ╠═272fb236-08b4-4628-843c-61d3147de5ff
# ╠═a62f2142-7d76-4a15-9566-54115aa98e08
# ╠═fd0d5424-074d-4fee-a41c-a80e61591348
# ╟─ccc3dffa-04c1-443f-9af6-725d9efbbe75
# ╠═c72f0ae1-e9df-4b88-9629-f471ba9f6e24
# ╠═ccc623a7-7a31-4bcf-9f09-9b7a97f8bf02
# ╟─04845384-01c0-4ef2-8d60-2a76d7e70731
# ╟─4a0f195e-7bf9-4f98-a98b-f5843a1628ef
# ╠═489692aa-9449-43f1-8809-d8b979d11ed0
# ╠═bd381c78-74ff-439b-87c0-b46f0ea8467d
# ╠═0006792d-5632-46d2-9fd4-49719f27aad6
# ╠═f831c365-1fdb-4286-98f1-214379b904d6
# ╠═75677613-21f4-4165-9257-dacd7febf42d
# ╠═e2252084-7070-42a2-a5dd-695e6d973543
# ╟─cba5b364-4620-4341-bc06-370e8a38f12a
# ╠═87ffd6f9-3b88-46ba-a2a1-a67c67043aa0
# ╠═5c3787e1-a188-410b-8a35-389c350ea9f6
# ╠═3187cc3c-7b47-4391-ac00-7b289ed0cd51
# ╟─ee10c49a-f996-4bee-8dc3-84e74ef6853a
# ╠═5ef29b01-7313-4558-8380-63e7f8101265
# ╠═2599199d-0542-4556-a2a3-ad8a155d8b6f
# ╠═9bd09e9b-3c09-4728-94b1-508716e1a55a
# ╠═ec1076d9-3cfc-420e-a985-bc0575c99be1
# ╠═a51806b9-392c-47f7-9b54-46dfaf515338
# ╠═321548a2-8223-4826-9475-cca8a6568d92
# ╟─5c9fe291-962e-4cd4-bc8d-75ff7c632cb4
# ╠═47723643-4c69-4f68-9b03-c6464fa65f4b
