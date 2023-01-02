### A Pluto.jl notebook ###
# v0.19.15

using Markdown
using InteractiveUtils

# ╔═╡ 13ddd21e-6998-11ed-0f20-235f1cee6e83
begin
	using Pkg; Pkg.activate()
	using BackwardCoalescent
	using TreeTools
end

# ╔═╡ c56ab812-b14b-4bee-b05f-31db10116ccd
BackwardCoalescent.genealogy

# ╔═╡ 72c88c94-ab48-4869-aeef-f0f967164408
begin
	N = 200
	β1 = .4
	β2 = .05
	ρ1 = N*β1^2
	ρ2 = N*β2^2
end

# ╔═╡ 5f34e7d0-b634-440c-8872-3149c37c1be1
ρ2

# ╔═╡ f0ab6f49-0592-4390-a1aa-e436bae39bdc
n = 30

# ╔═╡ ac4cd45d-bcfd-4f62-a1df-4da603a94cb3
ps1 = EFCoalescent(n, β1, ρ1)

# ╔═╡ bf73169a-a580-4f5a-937b-341d9b6a5808
ps2 = EFCoalescent(n, β2, ρ2)

# ╔═╡ 3380dc47-4fa9-4bc8-943b-8ae93c5440b0
kg = KingmanCoalescent(n, N)

# ╔═╡ efdfc65a-6aea-4d59-80ab-db9b87d0e0f3
let
	t_kg = genealogy(kg)
	TreeTools.write_newick("/tmp/tree_kingman.nwk", t_kg)
end

# ╔═╡ 0f41e1d8-7062-46ae-a1cb-3c4b4f14633d
let
	t_ps = genealogy(ps1)
	TreeTools.write_newick("/tmp/tree_partialsweep_1.nwk", t_ps)
end

# ╔═╡ eb95d3a7-fdd3-4526-9092-00674db1115c
let
	t_ps = genealogy(ps2)
	TreeTools.write_newick("/tmp/tree_partialsweep_2.nwk", t_ps)
end

# ╔═╡ 282f9e1a-d0dd-43b1-b6d0-479e7e7ec312


# ╔═╡ Cell order:
# ╠═13ddd21e-6998-11ed-0f20-235f1cee6e83
# ╠═c56ab812-b14b-4bee-b05f-31db10116ccd
# ╠═72c88c94-ab48-4869-aeef-f0f967164408
# ╠═5f34e7d0-b634-440c-8872-3149c37c1be1
# ╠═f0ab6f49-0592-4390-a1aa-e436bae39bdc
# ╠═ac4cd45d-bcfd-4f62-a1df-4da603a94cb3
# ╠═bf73169a-a580-4f5a-937b-341d9b6a5808
# ╠═3380dc47-4fa9-4bc8-943b-8ae93c5440b0
# ╠═efdfc65a-6aea-4d59-80ab-db9b87d0e0f3
# ╠═0f41e1d8-7062-46ae-a1cb-3c4b4f14633d
# ╠═eb95d3a7-fdd3-4526-9092-00674db1115c
# ╠═282f9e1a-d0dd-43b1-b6d0-479e7e7ec312
