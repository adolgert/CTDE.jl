# For individual model.

function plot_density(df, title, names)
    #myplot=plot(df, x="Times", y="Original", Geom.line)
    myplot=plot(df, layer(x="Times", y="Infect", Geom.line,
            Theme(default_color=color("blue"))),
        layer(x="Times", y="Clinical", Geom.line,
            Theme(default_color=color("green"))),
        layer(x="Times", y="Removed", Geom.line,
            Theme(default_color=color("orange"))),
        Guide.xlabel("time interval since infection [days]"),
        Guide.ylabel("probability distribution of firing"),
        Guide.title(title))
    filename=join(matchall(r"[A-Za-z]", title))
    draw(PDF("$(filename).pdf", 20cm, 15cm), myplot)
end

function plot_one_observer(d, name)
    sort!(d)
    plot_density(d, name)
end

function smoothed(plot_cnt, data::Array{Float64,2}, names)
    col_cnt=length(names)
    row_cnt=size(data)[1]
    cumulative=zeros(Float64, plot_cnt, col_cnt)
    times=zeros(Float64, plot_cnt)
    bandwidth=0.2
    lambda=1/bandwidth
    vmax=maximum(data)
    kernel=SmoothingKernels.kernels[:epanechnikov]
    for i in 1:plot_cnt
        dx=(i-1)*1.1*vmax/plot_cnt
        for col_idx in 1:col_cnt
            cumulative[i, col_idx]=sum(lambda*kernel(lambda*(data[:,col_idx]-dx)))
        end
        times[i]=dx
    end
    cumulative/=row_cnt
end


function plot_trajectory(df, title, names)
    myplot=plot(df, layer(x="Times", y="Susceptible", Geom.line,
            Theme(default_color=color("blue"))),
        layer(x="Times", y="Latent", Geom.line,
            Theme(default_color=color("orange"))),
        layer(x="Times", y="Infectious", Geom.line,
            Theme(default_color=color("red"))),
        layer(x="Times", y="Removed", Geom.line,
            Theme(default_color=color("black"))),
        layer(x="Times", y="Clinical", Geom.line,
            Theme(default_color=color("violet"))),
        layer(x="Times", y="Subclinical", Geom.line,
            Theme(default_color=color("green")))
    )
    filename=join(matchall(r"[A-Za-z]", title))
    draw(PDF("$(filename).pdf", 20cm, 15cm), myplot)
end
