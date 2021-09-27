include("common.jl")

if DFTK.mpi_nprocs() > 1
    disable_threading()
else
    setup_threading()
end

function plots_noniteractive()
    # Setup environment for making automated plots
    ENV["GKS_ENCODING"] = "utf8"
    ENV["GKSwstype"]    = "100"
    ENV["PLOTS_TEST"]   = "true"
end


config = Dict{String, Any}()
for d in readdir()
    if isdir(d) && isfile(joinpath(d, "config.jl"))
        include(joinpath(d, "config.jl"))
    end
end

if !isinteractive()
    plots_noniteractive()
    run_step_on_cases(dump_guess, config)
    run_step_on_cases(scf, config)
end
