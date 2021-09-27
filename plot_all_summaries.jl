include("common.jl")

config = Dict{String, Any}()
for d in readdir()
    if isdir(d) && isfile(joinpath(d, "config.jl"))
        include(joinpath(d, "config.jl"))
    end
end

# Setup environment for making automated plots
ENV["GKS_ENCODING"] = "utf8"
ENV["GKSwstype"]    = "100"
ENV["PLOTS_TEST"]   = "true"

for (key, cfg) in pairs(config)
    plot_convergence_summary(cfg)
end
