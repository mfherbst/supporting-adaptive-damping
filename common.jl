using DFTK
using MPI
using JSON
using JLD2
using Printf
using Plots
using Statistics
using LinearAlgebra
using Dates
using Logging
include("quantum_espresso.jl")

function build_magnetic_moments(atoms::AbstractArray; magmoms...)
    magmoms = Dict(magmoms)
    map(atoms) do (element, positions)
        element => fill(magmoms[element.symbol], length(positions))
    end
end

function attach_pseudos(atoms::AbstractArray; pseudomap...)
    pseudomap = Dict(pseudomap)

    map(atoms) do (element, position)
        pspfile = get(pseudomap, element.symbol, nothing)
        ElementPsp(element.symbol, psp=load_psp(pspfile)) => position
    end
end

function standard_cases()
    cases = Dict{String, Any}(
        "adaptive"    => (damping=DFTK.AdaptiveDamping(), ),  # Parameter-free adaptive damping
        "adaptive_05" => (damping=DFTK.AdaptiveDamping(0.5), ),
    )
    for α in 1:9
        cases["standard_0$α"] = (damping=DFTK.FixedDamping(α / 10), )
    end
    cases["standard_10"] = (damping=DFTK.FixedDamping(1.0), )
    cases
end

function random_logfile(prefix)
    logfile = "$(prefix)_$(@sprintf "%05d" rand(UInt16)).log"
    while isfile(logfile)
        logfile = "$(prefix)_$(@sprintf "%05d" rand(UInt16)).log"
    end
    logfile
end


function run_step_on_cases(payload, cases)
    spay = string(payload)
    header_printed = false

    anymissing = false
    for dir in sort(collect(keys(cases)))
        config = cases[dir]

        MPI.Barrier(MPI.COMM_WORLD)
        if isfile(joinpath(config.directory, string(payload) * ".done"))
            continue
        end

        for i in 1:10  # To avoid infinite loop
            MPI.Barrier(MPI.COMM_WORLD)
            res = nothing
            logfile = nothing
            if mpi_master()
                if !header_printed
                    println("#####" * "#"^length(spay) * "#####")
                    println("#--  $(string(payload))  --#")
                    println("#####" * "#"^length(spay) * "#####")
                    header_printed = true
                end

                logfile = random_logfile(joinpath(config.directory, string(payload)))
                println(basename(dir), "  ->  ", basename(dir) * "/" * basename(logfile))
                res = open(logfile, "w") do fp
                    redirect_stdout(fp) do
                        with_logger(SimpleLogger(stdout, Logging.Debug)) do
                            payload(config)
                        end
                    end
                end
                println("    ... $res\n")
            else
                res = payload(config)
            end

            if res == :exhausted
                if mpi_master()
                    rm(logfile)
                    anymissing || touch(joinpath(config.directory, string(payload) * ".done"))
                end
                break
            elseif res == :missing
                anymissing = true
                mpi_master() && rm(logfile)
            end
        end
    end

    MPI.Barrier(MPI.COMM_WORLD)
    nothing
end


function has_results(fileprefix)
    resultsfile = fileprefix * ".scfres.json"
    isfile(resultsfile) && open(JSON.parse, resultsfile, "r")["done"]
end


function basis_from_config(config)
    model = model_PBE(config.lattice, config.atoms;
                      magnetic_moments=config.magnetic_moments, config.kwargs_model...)
    fft_size = compute_fft_size(model, config.Ecut, supersampling=config.supersampling)
    if !haskey(config, :kgrid)
        kgrid = kgrid_size_from_minimal_spacing(config.lattice, 2π * 0.022)
    else
        kgrid=config.kgrid
    end
    PlaneWaveBasis(model, config.Ecut, kgrid=kgrid, kshift=[0, 0, 0], fft_size=fft_size)
end


function scf(config)
    basis = basis_from_config(config)

    scfkey = nothing
    if mpi_master()
        for scf in sort(collect(keys(config.cases)))
            lockfile = joinpath(config.directory, scf * ".running")
            fileprefix = joinpath(config.directory, scf)
            (isfile(lockfile) || has_results(fileprefix))  && continue
            write(lockfile, read(`hostname`, String))
            scfkey = scf
            break
        end
    end
    scfkey = MPI.bcast(scfkey, 0, MPI.COMM_WORLD)

    scfkey === nothing && return :exhausted
    fileprefix  = joinpath(config.directory, scfkey)
    lockfile    = joinpath(config.directory, scfkey * ".running")
    resultsfile = joinpath(config.directory, scfkey * ".scfres.json")
    guessfile   = joinpath(config.directory, "init.jld2")

    !isfile(guessfile) && return :missing
    ρ, ψ0, V = load_guess(basis, guessfile)

    try
        if mpi_master()
            println("#\n#-- SCF $(basename(config.directory)) -- case=$scfkey\n#")
            println("    keywords: $((; config.cases[scfkey]..., config.kwargs_scf...))")
        end

        MPI.Barrier(MPI.COMM_WORLD)
        mpi_master() && DFTK.reset_timer!(DFTK.timer)
        if haskey(config.cases[scfkey], :type) && config.cases[scfkey].type == :qe
            caseargs = Dict(k => v for (k, v) in pairs(config.cases[scfkey]) if k != :type)
            run_qe(basis, fileprefix, resultsfile;
                   magnetic_moments=config.magnetic_moments, caseargs...)
        else
            run_dftk(basis, resultsfile; ρ=ρ, ψ0=ψ0, V=V,
                     mixing=config.mixing, config.cases[scfkey]..., config.kwargs_scf...)
        end
        MPI.Barrier(MPI.COMM_WORLD)
        mpi_master() && println(DFTK.timer)
    catch
        mpi_master() && rm(joinpath(config.directory, "scf.done"), force=true)
        rethrow()
    finally
        mpi_master() && rm(lockfile)
    end
    mpi_master() && plot_convergence_summary(config)

    :done
end

function dump_guess(config)
    basis     = basis_from_config(config)
    guessfile = joinpath(config.directory, "init.jld2")
    isfile(guessfile) && return :exhausted

    ρ = guess_density(basis, config.magnetic_moments)
    _, ham = energy_hamiltonian(basis, nothing, nothing; ρ=ρ)

    V = nothing
    isnothing(config.potential) || (V = config.potential(basis))
    if !isnothing(V)
        ham = DFTK.hamiltonian_with_total_potential(ham, V)
    end

    n_bands = (haskey(config, :n_bands) ? (n_bands = config.n_bands)
                                        : DFTK.default_n_bands(basis.model))
    res_ψ0 = DFTK.next_density(ham; n_bands=DFTK.default_n_bands(basis.model),
                               n_ep_extra=3, miniter=2, tol=0.0001)
    diagiter = DFTK.mpi_mean(mean(res_ψ0.diagonalization.iterations), basis.comm_kpts)
    mpi_master() && println("guess:    diag = $diagiter  diagtol = 0.0001")
    isnothing(V) && (V = DFTK.total_local_potential(ham))

    # Gather ψ array on master
    ψ = DFTK.gather_kpts(res_ψ0.ψ, basis)
    if mpi_master()
        JLD2.jldopen(guessfile, "w") do jld
            jld["ρ"] = ρ
            jld["V"] = V
            jld["ψ"] = ψ
        end
    end

    :done
end


function load_guess(basis, guessfile)
    JLD2.jldopen(guessfile, "r") do jld
        ψ = jld["ψ"][basis.krange_thisproc]
        jld["ρ"], ψ, jld["V"]
    end
end


function run_dftk(basis, resultsfile;
                  ρ, ψ0, V, mixing=SimpleMixing(),
                  damping, tol=1e-11, type=:potmix, diagtol_max=1e-3, kwargs...)

    n_bands = (haskey(config, :n_bands) ? (n_bands = config.n_bands)
                                        : DFTK.default_n_bands(basis.model))
    if mpi_master()
        println()
        println("hostname        = $(read(`hostname`, String))")
        println("started on      = $(Dates.now())")
        println("julia threads   = $(Threads.nthreads())")
        println("BLAS threads    = $(BLAS.get_num_threads())")
        println("MPI procs       = $(mpi_nprocs())")
        println()
        println("temperature     = $(basis.model.temperature)")
        println("smearing        = $(basis.model.smearing)")
        println("lattice         = $(round.(basis.model.lattice, sigdigits=4))")
        println("Ecut            = $(basis.Ecut)")
        println("fft_size        = $(basis.fft_size)")
        println("kgrid           = $(basis.kgrid)")
        println("kshift          = $(basis.kshift)")
        println("n_irreducible_k = $(sum(length, basis.krange_allprocs))")
        println("n_bands         = $(n_bands)")
        println("n_electrons     = $(basis.model.n_electrons)")
        println("type            = $type")
        println("mixing          = $mixing")
        println("damping         = $damping")
        flush(stdout)
    end

    # Collection of data to dump later
    save_dict = Dict(
        "temperature"  => basis.model.temperature,
        "smearing"     => string(basis.model.smearing),
        "lattice"      => basis.model.lattice,
        "Ecut"         => basis.Ecut,
        "fft_size"     => basis.fft_size,
        "kgrid"        => basis.kgrid,
        "n_kirred"     => sum(length, basis.krange_allprocs),
        "n_bands"      => n_bands,
        "n_electr"     => basis.model.n_electrons,
        "mixing"       => string(mixing),
        "damping"      => string(damping),
        "type"         => string(type),
        "kwargs"       => string(kwargs),
        "done"         => false,
        #
        "energies"     => Float64[],
        "residuals"    => Float64[],
        "ndiag"        => Int[],
        "α"            => Float64[],
    )

    # Empty existing results
    if mpi_master()
        open(fp -> JSON.print(fp, save_dict), resultsfile, "w")
    end

    function ExtractCallback()
        function callback(info)
            if info.stage != :finalize
                push!(save_dict["energies"],  info.energies.total)
                push!(save_dict["ndiag"],     length(info.diagonalization))
                push!(save_dict["α"], isnan(info.α) ? DFTK.trial_damping(damping) : info.α)

                if type == :potmix
                    push!(save_dict["residuals"], norm(info.Vout - info.Vin))
                elseif type == :densmix
                    push!(save_dict["residuals"], norm(info.ρout - info.ρin))
                end
            end

            if mpi_master()
                open(fp -> JSON.print(fp, save_dict), resultsfile, "w")
            end
        end
    end

    determine_diagtol = DFTK.ScfDiagtol(; ratio_ρdiff=0.01, diagtol_max, diagtol_min=1e-12)
    callback = ExtractCallback() ∘ DFTK.ScfDefaultCallback()

    if type == :potmix
        run_scf = ((damping isa DFTK.AdaptiveDamping) ? DFTK.scf_potential_mixing_adaptive
                                                      : DFTK.scf_potential_mixing)
        scfres = run_scf(basis; tol=tol, ψ=ψ0, determine_diagtol=determine_diagtol,
                         mixing=mixing, callback=callback, damping=damping, ρ=ρ, V=V,
                         diag_miniter=2, kwargs...)
    elseif type == :densmix
        @assert damping isa Number
        scfres = self_consistent_field(basis; tol=tol, ψ=ψ0,
                                       determine_diagtol=determine_diagtol,
                                       mixing=mixing, callback=callback,
                                       damping=damping, ρ=ρ, kwargs...)
    else
        error("Unsupported SCF type: $type")
    end

    if mpi_master()
        println()
        println(scfres.energies)
        flush(stdout)
        save_dict["done"] = true
        open(fp -> JSON.print(fp, save_dict), resultsfile, "w")
    end

    scfres
end


function run_qe(basis, qe_directory, resultsfile; magnetic_moments,
                mixing_beta=0.7, kwargs...)
    if !mpi_master()
        MPI.Barrier(MPI.COMM_WORLD)
        return nothing
    end
    # This only runs on MPI master

    println()
    println("hostname        = $(read(`hostname`, String))")
    println("started on      = $(Dates.now())")
    println()
    flush(stdout)

    model = basis.model
    outputfile = joinpath(qe_directory, "espresso.pwo")
    if !isfile(outputfile)
        qeres = scf_quantum_espresso(model.lattice, model.atoms;
                                     Ecut=basis.Ecut,
                                     fileprefix=joinpath(qe_directory, "espresso"),
                                     temperature=model.temperature,
                                     mixing_beta=mixing_beta,
                                     kgrid=basis.kgrid,
                                     magnetic_moments=magnetic_moments,
                                     kwargs...)
    @assert qeres.outputfile == outputfile
    end

    # Copy output file to stdout
    println.(readlines(outputfile))
    flush(stdout)

    # Extract results
    save_dict = Dict(
        "done"         => true,
        "energies"     => Float64[],
        "residuals"    => Float64[],
        "ndiag"        => Int[],
        "α"            => Float64[],
    )
    for iteration in parse_quantum_espresso_pwo(outputfile)
        push!(save_dict["energies"],  iteration.total_energy)
        push!(save_dict["residuals"], iteration.accuracy)
        push!(save_dict["ndiag"],     1)
        push!(save_dict["α"],         mixing_beta)
    end
    open(fp -> JSON.print(fp, save_dict), resultsfile, "w")

    MPI.Barrier(MPI.COMM_WORLD)
    nothing
end


function collect_minimal_total_energy(config::NamedTuple; kwargs...)
    collect_minimal_total_energy(config.directory, config.cases; kwargs...)
end
function collect_minimal_total_energy(directory::AbstractString, cases::AbstractDict;
                                      kwargs...)
    collect_minimal_total_energy(directory, collect(keys(cases)); kwargs...)
end
function collect_minimal_total_energy(directory::AbstractString, cases::AbstractVector;
                                      tol_resid=0.05, tol_energy=1e-4)
    mintotals = Float64[]
    for scfkey in sort(cases)
        startswith(scfkey, "espresso_") && continue
        resultfile = joinpath(directory, scfkey * ".scfres.json")
        isfile(resultfile) || continue

        residuals = open(JSON.parse, resultfile, "r")["residuals"]
        energies  = open(JSON.parse, resultfile, "r")["energies"]
        energies === nothing && continue
        isempty(energies) && continue

        minidx = argmin(energies[max(end-4, 1):end]) + max(length(energies) - 4, 1) - 1
        residuals[minidx] > tol_resid && continue

        min_matching = findfirst(e -> abs(e - energies[minidx]) < tol_energy, mintotals)
        if min_matching === nothing
            push!(mintotals, energies[minidx])
        else
            mintotals[min_matching] = min(mintotals[min_matching], energies[minidx])
        end
    end

    mintotals
end


function parse_scfkey_label(scfkey)
    if scfkey == "adaptive"
        ("adaptive", 0.2)
    elseif startswith(scfkey, "adaptive_") && all(isnumeric, scfkey[10:end])
        digits = length(scfkey[10:end]) - 1
        alpha  = parse(Int, scfkey[10:end]) / 10^digits
        ("adaptive", alpha)
    elseif startswith(scfkey, "adaptive_")
        ("adaptive", scfkey[10:end])
    elseif startswith(scfkey, "standard_") && all(isnumeric, scfkey[10:end])
        digits = length(scfkey[10:end]) - 1
        alpha  = parse(Int, scfkey[10:end]) / 10^digits
        ("fixed", alpha)
    elseif startswith(scfkey, "espresso_") && all(isnumeric, scfkey[10:end])
        digits = length(scfkey[10:end]) - 1
        alpha  = parse(Int, scfkey[10:end]) / 10^digits
        ("espresso", alpha)
    elseif startswith(scfkey, "highprec_") && all(isnumeric, scfkey[10:end])
        digits = length(scfkey[10:end]) - 1
        alpha  = parse(Int, scfkey[10:end]) / 10^digits
        ("highprec", alpha)
    else
        (scfkey, nothing)
    end
end


function plot_convergence_summary(config; iteration=Inf, tol=1e-10, damping=true,
                                  outputfile=joinpath(config.directory, "summary.svg"),
                                  scfkeys=nothing, kwargs...)
    function deduce_label(scfkey)
        method, alpha = parse_scfkey_label(scfkey)
        method * (isnothing(alpha) ? "" : " (α = $alpha)")
    end

    !mpi_master() && error("Plotting only on MPI master supported.")
    mintotals = collect_minimal_total_energy(config)
    isempty(mintotals) && return plot()  # Return on no data

    if scfkeys === nothing
        scfkeys = collect(keys(config.cases))
    end

    errors   = Dict{String, Vector{Float64}}()
    cumdiags = Dict{String, Vector{Int}}()
    dampings = Dict{String, Vector{Float64}}()
    whichmin = Dict{String, Int}()
    for scfkey in sort(scfkeys)
        method, _ = parse_scfkey_label(scfkey)
        resultfile = joinpath(config.directory, scfkey * ".scfres.json")
        isfile(resultfile) || continue

        energies  = open(JSON.parse, resultfile, "r")["energies"]
        ndiags    = open(JSON.parse, resultfile, "r")["ndiag"]
        damps     = open(JSON.parse, resultfile, "r")["α"]

        energies === nothing && continue
        if method == "espresso"
            values = open(JSON.parse, resultfile, "r")["residuals"]
            mintrial  = minimum(energies[max(end-4, 1):end])
            bestmatch = sortperm(abs.(mintrial .- mintotals))[1]
        else
            isempty(energies) && continue
            # Find best matching minimum
            mintrial  = minimum(energies[max(end-4, 1):end])
            bestmatch = sortperm(abs.(mintrial .- mintotals))[1]
            values = abs.(energies .- mintotals[bestmatch])
            idx = findfirst(values .< 1e-16)
            if !isnothing(idx)
                values = values[1:idx-1]
                ndiags = ndiags[1:idx-1]
                damps  = damps[1:idx-1]
            end
        end
        isempty(values) && continue

        whichmin[scfkey] = bestmatch
        errors[scfkey]   = values
        dampings[scfkey] = damps
        cumdiags[scfkey] = cumsum(ndiags)
    end

    # If no data return
    isempty(errors) && return plot()

    any_adaptive = false
    ymax = 0
    maxiter = 0
    colors = Dict{Any, Int}()
    p = plot(; yaxis=:log, legend=:topright, title=basename(config.directory), kwargs...)
    for (i, scfkey) in enumerate(sort(collect(keys(errors))))
        label = deduce_label(scfkey)
        maxiter = max(maxiter, length(errors[scfkey]))
        ymax = max(ymax, maximum(errors[scfkey]))

        if length(mintotals) > 1
            label *= " (min $(whichmin[scfkey]))"
        end

        method, alpha = parse_scfkey_label(scfkey)
        if isnothing(alpha)
            color = i
        elseif alpha in keys(colors)
            color = colors[alpha]
        else
            colors[alpha] = i
            color = colors[alpha]
        end

        mark = :x
        if method == "adaptive"
            mark = :o
            any_adaptive = true
        elseif method == "espresso"
            mark = :rtriangle
        elseif method == "highprec"
            mark = :ltriangle
        end
        plot!(p, cumdiags[scfkey], errors[scfkey];
              label=label, color=color, mark=mark, lw=1.3, markersize=4)
    end
    ylims!(p, (tol, 2 * ymax))
    ylabel!(p, "Total energy absolute error")

    if damping && any_adaptive
        q = twinx()
        for scfkey in sort(collect(keys(errors)))
            method, alpha = parse_scfkey_label(scfkey)
            method == "adaptive" || continue
            (alpha in keys(colors)) || continue
            plot!(q, cumdiags[scfkey], abs.(dampings[scfkey]), mark=:+,
                  color=colors[alpha], label="", linestyle=:dot)
        end
        ylims!(q, (0.0, 1.0))
        ylabel!(q, "Damping parameter")
    end
    xlims!(p, (0, min(iteration, maxiter)))
    xlabel!(p, "Iteration / Hamiltonian diagonalization")

    if isnothing(outputfile)
        p
    else
        savefig(p, outputfile)
    end
end
