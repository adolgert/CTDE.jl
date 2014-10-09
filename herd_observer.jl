include("herd_trajectory.jl")

type IndividualDiseaseObserver
    t::Array{Float64,2}
    latent::Float64 # When latent was started.
    previous_time::Float64
    idx
    IndividualDiseaseObserver(run_cnt)=new(zeros(Float64, run_cnt, 4), 0., 0., 1)
end

function observe(eo::IndividualDiseaseObserver, state)
    last_fired=state.last_fired
    delta=state.current_time-eo.previous_time
    const ind_fired={ 'f'=>1, 'a'=>2, 'd'=>3, 'e'=>4, 'z'=>5 }
    idx=ind_fired[state.last_fired]
    # nonzero=Char[]
    # for p in Char['l', 'i', 'r', 'n', 'c']
    #   if length(state.marking, p)>0
    #       push!(nonzero, p)
    #   end
    # end
    # println(join(nonzero), " ", eo.idx, " ", idx, " ", state.current_time, " ", eo.latent)
    if idx<5
        eo.t[eo.idx, idx]=state.current_time-eo.latent
    elseif idx==5
        # Back to the beginning.
        eo.latent=state.current_time
        eo.idx+=1
    end

    eo.previous_time=state.current_time
    eo.idx<=size(eo.t)[1]
end


function plot_observer(eo::IndividualDiseaseObserver, title)
    names=["time to infectious", "time to clinical", "time to recovered"]
    cumulative=smoothed(100, eo.t[:,1:3], names)
    df=DataFrame(Times=times, Infect=cumulative[:,1], Clinical=cumulative[:,2],
            Removed=cumulative[:,3])
    plot_density(df, title, names)
end


########################################

typealias Time Float64

type TrajectoryStruct
  s::Int64
  l::Int64
  i::Int64
  r::Int64
  n::Int64
  c::Int64
  t::Time
  TrajectoryStruct()=new(0,0,0,0,0,0,0.)
  TrajectoryStruct(s_,l_,i_,r_,n_,c_,t_)=new(s_,l_,i_,r_,n_,c_,t_)
end

type HerdDiseaseObserver
    state::Array{Int64,2}
    time::Array{Time,1}
    cnt::Int
    sir::TrajectoryStruct
    previous_time::Float64
    HerdDiseaseObserver(cnt)=new(zeros(Int64, 10_000,6),
            fill(-1., 10_000),
            1, TrajectoryStruct(int(cnt-1), 0, 1, 0, cnt, 0, 0.), 0.)
end

function observe(eo::HerdDiseaseObserver, state)
    last_fired=state.last_fired
    delta=state.current_time-eo.previous_time
    shiftstore= xo->begin
        xo.sir.t=state.current_time
        xo.state[xo.cnt,:]=Int[xo.sir.s, xo.sir.l, xo.sir.i, xo.sir.r, xo.sir.n,
                xo.sir.c]
        xo.time[xo.cnt]=xo.sir.t
        xo.cnt+=1
    end
    key=last_fired[3]
    if key=='d'
        eo.sir.i-=1
        eo.sir.r+=1
        shiftstore(eo)
    elseif key=='g'
        eo.sir.s-=1
        eo.sir.l+=1
        shiftstore(eo)
    elseif key=='f'
        eo.sir.l-=1
        eo.sir.i+=1
        shiftstore(eo)
    elseif key=='a'
        eo.sir.n-=1
        eo.sir.c+=1
        shiftstore(eo)
    elseif key=='e'
        eo.sir.c-=1
        eo.sir.n+=1
        shiftstore(eo)
    end
    if eo.cnt>length(eo.state)
        new_len=2*eo.cnt
        new_t=zeros(Int64, new_len, 6)
        new_t[1:length(eo.state)]=eo.state
        eo.state=new_t
        new_time=zeros(Time, new_len)
        new_time[1:length(eo.time)]=eo.time
        eo.time=new_time
    end
    eo.previous_time=state.current_time
    true # run until you can't run any more.
end

function show(eo::HerdDiseaseObserver)
  for i in 1:(eo.cnt-1)
      println(eo.state[i,1], " ", eo.state[i,2], " ", eo.state[i,3], " ",
        eo.state[i,4], " ", eo.state[i,5], " ", eo.time[i])
  end
end


function plot_observer(eo::HerdDiseaseObserver, title)
    names=["susceptible", "latent", "infectious", "removed",
        "clinical", "subclinical"]
    state=eo.state[1:(eo.cnt-1),:]
    df=DataFrame(Times=eo.time[1:(eo.cnt-1)], Susceptible=state[:,1],
        Latent=state[:,2], Infectious=state[:,3],
        Removed=state[:,4],Subclinical=state[:,5],
        Clinical=state[:,6])
    plot_trajectory(df, title, names)
end

function save_observer(eo::HerdDiseaseObserver, filename)
    save_trajectory(eo.state[1:(eo.cnt-1),:], eo.time[1:(eo.cnt-1)], filename)
end
