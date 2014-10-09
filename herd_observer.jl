
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
    df=smoothed(100, eo.t[:,1:3], names)
    plot_density(df, title, names)
end


########################################

typealias Time Float64
typealias TrajectoryEntry (Int64,Int64,Int64,Float64)
type TrajectoryStruct
  s::Int64
  i::Int64
  r::Int64
  t::Time
  TrajectoryStruct()=new(0,0,0,0.)
  TrajectoryStruct(s_,i_,r_,t_)=new(s_,i_,r_,t_)
end

function convert(::Type{TrajectoryEntry}, x::TrajectoryStruct)
    (x.s, x.i, x.r, x.t)
end

type HerdDiseaseObserver
    t::Array{TrajectoryEntry,1}
    cnt::Int
    sir::TrajectoryStruct
    previous_time::Float64
    HerdDiseaseObserver(cnt)=new(Array(TrajectoryEntry, 10_000),
            1, TrajectoryStruct(int(cnt-1), 1, 0, 0.), 0.)
end

function observe(eo::HerdDiseaseObserver, state)
    last_fired=state.last_fired
    delta=state.current_time-eo.previous_time
    const ind_fired={ 'f'=>1, 'a'=>2, 'd'=>3, 'e'=>4, 'z'=>5, 'g'=>6 }
    shiftstore= xo->begin
        xo.sir.t=state.current_time
        xo.t[xo.cnt]=convert(TrajectoryEntry, xo.sir)
        xo.cnt+=1
    end
    if length(last_fired)==2
        if last_fired[2]=='d'
            eo.sir.i-=1
            eo.sir.r+=1
            shiftstore(eo)
        end
    else
        eo.sir.s-=1
        eo.sir.i+=1
        shiftstore(eo)
    end
    if eo.cnt>length(eo.t)
        new_len=2*eo.cnt
        new_t=Array(TrajectoryEntry, new_len)
        new_t[1:length(eo.t)]=eo.t
        eo.t=new_t
    end
    eo.previous_time=state.current_time
    true # run until you can't run any more.
end

function show(eo::HerdDiseaseObserver)
  for i in 1:(eo.cnt-1)
      println(eo.t[i][1], " ", eo.t[i][2], " ", eo.t[i][3], " ", eo.t[i][4])
  end
end