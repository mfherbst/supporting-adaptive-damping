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

# Special case for AlVac:
# For AlVac_Kerker no calculation is converged, so use the converged energy from AlVac
# to get a meaningful plot
cfg_AlVac        = config[joinpath(@__DIR__, "AlVac")]
cfg_AlVac_Kerker = config[joinpath(@__DIR__, "AlVac_Kerker")]
plot_convergence_summary(cfg_AlVac_Kerker, mintotals=collect_minimal_total_energy(cfg_AlVac))
