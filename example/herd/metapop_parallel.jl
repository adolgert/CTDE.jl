# You would run it like this. The -p argument is the number of processes.
# julia -p 4 sirs_parallel.jl 50 20 4.0 1.0 0.5 34

# Uncomment for HDF5. Pkg.add("HDF5")
# using HDF5, JLD

require("metapop_herd.jl")

individual_cnt=int(ARGS[1])
block_cnt=int(ARGS[2])
block_length=int(ARGS[3])
run_cnt=int(ARGS[4])
fenceline=float(ARGS[5])
mixed=float(ARGS[6])
seed=int(ARGS[7])

for init_idx in 2:nprocs()
    remotecall(init_idx, set_rng, seed)
end
@everywhere save_file="metapop$(myid()).h5"

extra_params={"individual_cnt"=>individual_cnt, "block_cnt"=>block_cnt,
        "block_length"=>block_length, "fenceline"=>fenceline,
        "mixed"=>mixed}
work=Array(Any,run_cnt)
for i in 1:run_cnt
    work[i]=extra_params
end
println("run count ", run_cnt)
r=pmap(work) do package
  herd_model_remote(package)
end
