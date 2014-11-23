# You would run it like this. The -p argument is the number of processes.
# julia -p 4 sirs_parallel.jl 50 20 4.0 1.0 0.5 34

# Uncomment for HDF5. Pkg.add("HDF5")
# using HDF5, JLD

if length(ARGS)<6
    println("Need some args. Try 50 20 4.0 1.0 0.5 34")
    exit(1)
end

individual_cnt=int(ARGS[1])
run_cnt=int(ARGS[2])
beta=float(ARGS[3])
gamma=float(ARGS[4])
wane=float(ARGS[5])
seed=int(ARGS[6])

require("sirs.jl")

disease_exponential={
    'i'=>beta/individual_cnt, # rate of infection of neighbor
    'r'=>gamma, # infectious to removed
    'w'=>wane, # removed to susceptible
}
observation_times=Time[5.0, 10.0, 15.0]

for init_idx in 2:nprocs()
    println("Setting seed for $init_idx")
    remotecall(init_idx, set_rng, seed)
end

work=Array(Any,run_cnt)
for i in 1:run_cnt
    work[i]=(disease_exponential, individual_cnt, observation_times)
end

r=pmap(work) do package
  apply(herd_single, package)
end

results=zeros(Int, run_cnt, 3, length(observation_times))
# Each entry is [(s, i, r, t)*] where t=0 for no entry.
for (run_idx, entry) in enumerate(r)
    for (obs_idx, obs) in enumerate(entry)
        if obs[4]>0.0001
            results[run_idx, 1, obs_idx]=obs[1]
            results[run_idx, 2, obs_idx]=obs[2]
            results[run_idx, 3, obs_idx]=obs[3]
        else
            results[run_idx, 1, obs_idx]=-1
            results[run_idx, 2, obs_idx]=-1
            results[run_idx, 3, obs_idx]=-1
        end
    end
end
writedlm("z.txt", results, ',')
# h5write("z.h5", "sir", results)
# Use "h5dump z.h5" to see what's in there.
