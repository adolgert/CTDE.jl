using HDF5
using Dates

function directory_git(directory)
    (repo, version)=cd(directory) do
        repo=readall(`git remote show origin`)
        version=strip(readall(`git log --pretty=%H HEAD^..HEAD`))
        (repo, version)
    end
    repository=""
    m=match(r"URL:\s*(\S+)\n", repo)
    if m!=nothing
        repository=m.captures[1]
    end
    Dict{String,String}({"repository"=>repository, "version"=>version})
end

directory_git()=directory_git(".")

function system_info()
    return Dict{String,String}({"OS"=>string(OS_NAME), "CORES"=>string(CPU_CORES),
            "JuliaVersion"=>string(VERSION)})
end

function invocation()
    return Dict{String,String}({"Arguments"=>join(ARGS, ":"),
            "Location"=>string(functionloc(directory_git)),
            "DateTime"=>string(Dates.now())})
end

function metadata()
    md=Dict{String,String}()
    merge!(md, directory_git())
    merge!(md, system_info())
    merge!(md, invocation())
    md
end

function save_trajectory(state::Array{Int64,2}, times::Array{Float64,1},
        filename)
    mode="w"
    if isfile(filename)
        mode="r+"
    end
    h5open(filename, mode) do file
        last=0
        dirnames=names(file["/"])
        if length(dirnames)>0
            last=maximum(Int[int(x) for x in dirnames])
        end
        #println("last entry ", last)
        g=g_create(file, string(last+1))
        g["state"]=state
        g["time"]=times
        for (key, value) in metadata()
            attrs(g)[key]=convert(ASCIIString, value)
        end
    end
end


function test_save_trajectory()
    state=Int[0 1; 2 3; 4 5]
    times=Float64[0.1, 9.2, 16.3]
    save_trajectory(state, times, "z.h5", 3)
end

