require("Rif")
using Rif
initr()

r_base = Rif.importr("base")
r_stats = Rif.importr("stats")
r_graphics = Rif.importr("graphics")

m = r_base.matrix(r_stats.rnorm(100); nrow=20)

# A Julia matrix mj of type (Array{Float64, 2}) could
# be used with
# m = RArray{Float64,2}(mj)

d = r_stats.dist(m)
hc = r_stats.hclust(d)
r_graphics.plot(hc; 
                sub=cR(""),
                xlab=cR(""))