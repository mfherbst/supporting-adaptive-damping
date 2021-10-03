include("common.jl")

function setup_plots()
    # Setup environment for making automated plots
    ENV["GKS_ENCODING"] = "utf8"
    ENV["GKSwstype"]    = "100"
    ENV["PLOTS_TEST"]   = "true"

    gr()
    default(size=tuple(Int.(ceil.(0.75 .* [600, 400]))...),
            guidefontsize=10, grid=false)
end

function dump_table_convergence(fn::AbstractString; kwargs...)
    open(fn, "w") do io
        dump_table_convergence(io; kwargs...)
    end
end
function dump_table_convergence(; kwargs...)
    dump_table_convergence(stdout; kwargs...)
end
function dump_table_convergence(io::IO;
                                selected=["Fe2CrGa"], allscfkeys=nothing,
                                tol=1e-10, ncthresh=100, labels=Dict{String,String}(),
                                precons=nothing,
                               )
    separator_after = String[]
    for (i, case) in enumerate(selected)
        if i+1 <= length(selected) && selected[i+1] == :separator
            push!(separator_after, case)
        end
    end
    selected = filter!(a -> a != :separator, selected)

    if isnothing(allscfkeys)
        allscfkeys = union(Set.(default_scfkeys(case) for case in selected)...)
    end
    allscfkeys = sort([allscfkeys...], by=parse_scfkey_label)
    adaptkeys  = filter(startswith("adaptive"), allscfkeys)
    fixedkeys  = filter(startswith("standard"), allscfkeys)
    allkeys    = vcat(fixedkeys, adaptkeys)
    n_adaptive = length(adaptkeys)
    n_standard = length(fixedkeys)

    if length(adaptkeys) == 1
        adaptheader = "adaptive"
    else
        adaptheader = "adaptive damping"
    end


    println(io, raw"\begin{tabular}{ll|*{" * string(n_standard) * raw"}{c}|*{"
            * string(n_adaptive) *  raw"}{c}}")
    println(io, raw"System&Precond." *
                raw"&\multicolumn{" * string(n_standard) * raw"}{c|}{fixed damping $\alpha$}" *
                raw"&\multicolumn{" * string(n_adaptive) * raw"}{c}{" * adaptheader * raw"}" *
                raw"\\\\")

    print(io, "&")
    for scfkey in fixedkeys
      _, alpha = parse_scfkey_label(scfkey)
      print(io, "& $alpha ")
    end
    if length(adaptkeys) == 1
        print(io, "&damping")
    else
        for scfkey in adaptkeys
          _, alpha = parse_scfkey_label(scfkey)
          print(io, raw"& $\atmin = " * string(alpha) * raw"$ ")
        end
    end
    println(io, "\\\\")
    println(io, raw"\hline")

    for case in selected
        directory = joinpath(@__DIR__, case)
        collectkeys = default_scfkeys(case)
        mintotals = collect_minimal_total_energy(directory, collectkeys)

        min_n_diags_fixed = typemax(Int)
        n_diags = Vector{Union{Missing, Int}}()
        for scfkey in allkeys
            if isempty(mintotals)
                push!(n_diags, typemax(Int))
                continue
            end

            method, _ = parse_scfkey_label(scfkey)
            jsonfile = joinpath(directory, "$scfkey.scfres.json")
            if !isfile(jsonfile)
                push!(n_diags, missing)
                continue
            end

            data = open(JSON.parse, jsonfile)
            energies = data["energies"]
            mintrial  = minimum(energies[max(end-4, 1):end])
            bestmatch = sortperm(abs.(mintrial .- mintotals))[1]
            values = abs.(energies .- mintotals[bestmatch])
            n_conv = findfirst(values .< tol)
            if isnothing(n_conv) && values[end] < 3tol
                n_conv = length(values)
            end
            if case == "Fe2CrGa" && scfkey == "standard_09"
                # Dirty hack, because I manually verified this to be (barely) not converged
                n_conv = nothing
            end
            if isnothing(n_conv)
                push!(n_diags, typemax(Int))
            else
                # unit_cell_volume = abs(det(hcat(data["lattice"]...)))
                # dvol = unit_cell_volume / prod(data["fft_size"])
                # println(case, "   ", scfkey, "   ", data["residuals"][1] * sqrt(dvol),
                #               "   ", data["residuals"][n_conv] * sqrt(dvol),
                #               "   ", data["residuals"][n_conv]/data["residuals"][1])

                n_diag = sum(data["ndiag"][1:n_conv])
                push!(n_diags, n_diag)
                if method == "fixed"
                    min_n_diags_fixed = min(min_n_diags_fixed, n_diag)
                end
            end
        end

        if case in keys(labels)
            print(io, labels[case], raw"  &  ")
        else
            print(io, raw"\ce{", case, raw"}  &  ")
        end
        if case in keys(precons)
            print(io, precons[case], raw"  &  ")
        else
            print(io, raw"--  &  ")
        end

        for (scfkey, n_diag) in zip(allkeys, n_diags)
            method, _ = parse_scfkey_label(scfkey)
            endchar = scfkey == allkeys[end] ? raw" \\\\" : " & "

            if ismissing(n_diag)
                print(io, "--", endchar)
            elseif n_diag > ncthresh
                print(io, raw"$\times$", endchar)
            elseif method == "fixed" && n_diag == min_n_diags_fixed
                print(io, raw"\textbf{", n_diag, "}", endchar)
            else
                print(io, n_diag, endchar)
            end
        end
        println(io)
        if case in separator_after
            println(io, raw"\hline")
        end
    end
    println(io, raw"\end{tabular}")
end

function dump_latex_tables()
    file = "convergence_table.tex"
    allscfkeys = ["adaptive"]
    for i in 1:9
        push!(allscfkeys, "standard_0$i")
    end
    push!(allscfkeys, "standard_10")
    selected = ["Al_nodiis_Kerker", "Al_nodiis", "Al_Kerker", "Al", "AlVac_Kerker", "AlVac",
                :separator,
                "GaAs",
                :separator,
                "CoFeMnGa", "Fe2CrGa", "Fe2MnAl", "FeNiF6", "Mn2RuGa", "Mn3Si", "Mn3Si_AFM",
                :separator,
                "Cr19", "WFe"
               ]
    labels = Dict(
        "Al_nodiis"        => (raw"\ce{Al8} supercell"),
        "Al_nodiis_Kerker" => (raw"\ce{Al8} supercell"),
        "Al"               => (raw"\ce{Al40} supercell"),
        "Al_Kerker"        => (raw"\ce{Al40} supercell"),
        "GaAs"             => (raw"\ce{Ga20As20} supercell"),
        "AlVac_Kerker"     => (raw"\ce{Al40} surface"),
        "AlVac"            => (raw"\ce{Al40} surface"),
        "Mn3Si_AFM"        => (raw"\ce{Mn3Si}$^\text{AFM}$"),
        "Cr19"             => (raw"\ce{Cr19} defect"),
        "WFe"              => (raw"\ce{Fe28W8} multilayer\cite{Marks2021}"),
    )
    precons = Dict(
        "Al_nodiis"        => raw"None$^\dagger$",
        "Al_nodiis_Kerker" => raw"Kerker$^\dagger$",
        "Al_Kerker"        => "Kerker",
        "Al"               => "None",
        "GaAs"             => "None",
        "AlVac"            => "None",
        "AlVac_Kerker"     => "Kerker",
        "CoFeMnGa"         => "Kerker",
        "Fe2CrGa"          => "Kerker",
        "Fe2MnAl"          => "Kerker",
        "FeNiF6"           => "Kerker",
        "Mn2RuGa"          => "Kerker",
        "Mn3Si"            => "Kerker",
        "Mn3Si_AFM"        => "Kerker",
        "Cr19"             => "Kerker",
        "WFe"              => "Kerker",
    )
    dump_table_convergence(file; selected, labels, precons, allscfkeys)
end

function purify_maxcond(jsonfile, mintotals)
    energies = open(JSON.parse, jsonfile, "r")["energies"]
    ndiags   = open(JSON.parse, jsonfile, "r")["ndiag"]

    @assert length(mintotals) == 3
    minenes = [abs.(energies .- mintotals[1]),
               abs.(energies .- mintotals[2]),
               abs.(energies .- mintotals[3])]
    ibest = sortperm(minimum.(minenes))[1]
    values = minenes[ibest]
    idx = findfirst(values .< 1e-16)
    if !isnothing(idx)
        values = values[1:idx-1]
        ndiags = ndiags[1:idx-1]
    end

    (; values, cumdiags=cumsum(ndiags))
end


function plot_maxcond(; iteration=80, tol=1e-10, legend=:bottomleft, label=nothing)
    mintotals = [
        -333.3052730117624,
        -333.30261995802,
        -333.30262027079226,
    ]

    errors   = Dict{String, Vector{Float64}}()
    cumdiags = Dict{String, Vector{Int}}()
    for maxcond in (2, 4, 6, 8, 10)
        key = @sprintf "%02i" maxcond
        jsonfile = joinpath(@__DIR__, "Fe2CrGa_maxcond$key/standard_02.scfres.json")
        dlabel = "maxcond = 1e$maxcond"
        pf = purify_maxcond(jsonfile, mintotals)
        errors[dlabel]   = pf.values
        cumdiags[dlabel] = pf.cumdiags
    end

    for key in ("1e2", "1e4", "1e6", "1e8", "1e10")
        jsonfile = joinpath(@__DIR__, "Fe2CrGa_40guess/maxcond_$key.scfres.json")
        dlabel = "maxcond = $key (restart)"
        pf = purify_maxcond(jsonfile, mintotals)
        errors[dlabel]   = pf.values
        cumdiags[dlabel] = 40 .+ pf.cumdiags
    end

    titleargs = ()
    if !isnothing(label)
        titleargs = (title=label, titlelocation=:left, titlefontsize=10)
    end

    ymax = 0
    maxiter = 0
    p = plot(; yaxis=:log, legend, titleargs...)
    for (i, maxcond) in enumerate([2, 4, 6, 8, 10])
        dlabel = "maxcond = 1e$maxcond"
        maxiter = max(maxiter, length(errors[dlabel]))
        ymax    = max(ymax,    maximum(errors[dlabel]))

        plot!(p, cumdiags[dlabel], errors[dlabel];
              label=dlabel, color=i, mark=:x, linestyle=:solid, markerstrokecolor=i,
              lw=1.5, markersize=4)

        rlabel = dlabel * " (restart)"
        if rlabel in keys(errors)
            plot!(p, cumdiags[rlabel], errors[rlabel];
                  label=rlabel, color=i, mark=:+, linestyle=:dot, markerstrokecolor=i,
                  lw=1.5, markersize=4)
        end
    end
    ylims!(p, (tol, 2 * ymax))
    ylabel!(p, "Total energy absolute error")

    xlims!(p, (0, min(iteration, maxiter)))
    xlabel!(p, "Number of Hamiltonian diagonalizations")

    p
end

function default_scfkeys(case)
    directory  = joinpath(@__DIR__, case)
    [key[1:end-12] for key in readdir(directory)
     if (    endswith(key, ".scfres.json")
         && !occursin("adaptive_slowstart", key)
         && !occursin("espresso", key)
         && !occursin("density", key)
         && !occursin("nlsolve", key)
         && !occursin("highprec", key))
    ]
end

function plot_convergence(case; iteration=60, tol=1e-10,
                          scfkeys=default_scfkeys(case),
                          colormap=Dict{Float64,Int}(),
                          legend=:topright, label=nothing,
                          mintotals=nothing, showmins=true, show_αmin=true,
                          show_damping=false)
    # First pass: For DFTK runs collect mimimal total energy
    directory = joinpath(@__DIR__, case)
    collectkeys = default_scfkeys(case)
    if isnothing(mintotals)
        mintotals = collect_minimal_total_energy(directory, collectkeys)
    end
    isempty(mintotals) && return plot()

    errors   = Dict{String, Vector{Float64}}()
    cumdiags = Dict{String, Vector{Int}}()
    dampings = Dict{String, Vector{Float64}}()
    whichmin = Dict{String, Int}()
    for scfkey in sort(scfkeys)
        method, _ = parse_scfkey_label(scfkey)
        method in ("espresso", "highprec") && error(
            "espresso or highprec not supported."
        )

        resultfile = joinpath(directory, scfkey * ".scfres.json")
        isfile(resultfile) || error("Invalid scfkey $scfkey")

        energies  = open(JSON.parse, resultfile, "r")["energies"]
        ndiags    = open(JSON.parse, resultfile, "r")["ndiag"]
        damps     = open(JSON.parse, resultfile, "r")["α"]

        energies === nothing && continue
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
        isempty(values) && continue

        whichmin[scfkey] = bestmatch
        errors[scfkey]   = values
        dampings[scfkey] = damps
        cumdiags[scfkey] = cumsum(ndiags)
    end
    isempty(errors) && return plot()

    titleargs = ()
    if !isnothing(label)
        titleargs = (title=label, titlelocation=:left, titlefontsize=10)
    end

    ymax = 0
    maxiter = 0
    p = plot(; yaxis=:log, legend, titleargs...)
    for (i, scfkey) in enumerate(sort(collect(keys(errors)), by=parse_scfkey_label))
        method, alpha = parse_scfkey_label(scfkey)
        if method == "adaptive" && show_αmin
            label = method * " (α̃ₘᵢₙ = $alpha)"
        elseif method == "adaptive"
            label = "adaptive"
        elseif method == "fixed"
            label = method * " (α = $alpha)"
        else
            error("Unknown method: $method")
        end
        if length(mintotals) > 1 && showmins
            label *= " (min $(whichmin[scfkey]))"
        end
        maxiter = max(maxiter, length(errors[scfkey]))
        ymax    = max(ymax,    maximum(errors[scfkey]))

        if alpha in keys(colormap)
            color = colormap[alpha]
        else
            color = isempty(colormap) ? 1 : 1 + maximum(values(colormap))
            colormap[alpha] = color
        end
        if method == "adaptive"
            mark = :o
        elseif method == "fixed"
            mark = :x
        else
            error("Unknown method")
        end
        linestyle = :solid

        plot!(p, cumdiags[scfkey], errors[scfkey];
              label, color, mark, linestyle, markerstrokecolor=color,
              lw=1.5, markersize=4)
    end
    ylims!(p, (tol, 2 * ymax))
    ylabel!(p, "Total energy absolute error")

    if show_damping
        q = twinx()
        for (i, scfkey) in enumerate(sort(collect(keys(errors)), by=parse_scfkey_label))
            method, alpha = parse_scfkey_label(scfkey)
            method == "adaptive" || continue
            (alpha in keys(colormap)) || continue
            plot!(q, cumdiags[scfkey], abs.(dampings[scfkey]), mark=:+,
                  color=colormap[alpha], label="", linestyle=:dot)
        end
        ylims!(q, (0.0, 1.0))
        ylabel!(q, "Damping parameter")
    end

    xlims!(p, (0, min(iteration, maxiter)))
    xlabel!(p, "Number of Hamiltonian diagonalizations")

    p
end

function dump_plots()
    setup_plots()

    let kwargs = ()
        case  = "Al_nodiis_Kerker"
        label = "(b) Al₈ supercell (Kerker, no Anderson)"
        scfkeys = ["standard_02", "standard_05", "standard_09",
                   "adaptive",    "adaptive_05", "adaptive_09"]
        colormap = Dict{Float64,Int}()
        p = plot_convergence(case; colormap, scfkeys, label, kwargs...)
        ylims!(p, (1e-10, 10))
        savefig(p, case * ".pdf")

        case  = "Al_nodiis"
        label = "(a) Al₈ supercell (no preconditioner, no Anderson)"
        scfkeys = ["standard_01", "standard_02", "standard_03",
                   "adaptive_01", "adaptive_03", "adaptive_05"]
        p = plot_convergence(case; colormap, scfkeys, label, kwargs...)
        ylims!(p, (1e-10, 10))
        savefig(p, case * ".pdf")
    end

    let kwargs = ()
        mintotals = [-333.3052730117624, -333.30261995802, -333.30262027079226]
        case  = "Fe2CrGa"
        label = "(b) Fe₂CrGa (Kerker, Anderson)"
        colormap = Dict{Float64,Int}()
        scfkeys = ["standard_02", "standard_03", "standard_04", "standard_05",
                   "adaptive",    "adaptive_03", "adaptive_04", "adaptive_05"]
        p = plot_convergence(case; colormap, scfkeys, label, showmins=false, mintotals, kwargs...)
        savefig(p, case * ".pdf")

        case  = "Fe2CrGa_nodiis"
        label = "(a) Fe₂CrGa (Kerker, no Anderson)"
        scfkeys = ["standard_005", "standard_010", "standard_015",
                   "adaptive_01", "adaptive", "adaptive_03"]
        p = plot_convergence(case; colormap, scfkeys, label,
                             legend=:bottomleft, mintotals, showmins=false, kwargs...)
        savefig(p, case * ".pdf")
    end

    let kwargs = ()
        colormap = Dict{Float64,Int}()
        case = [("GaAs", "Ga₂₀As₂₀ supercell (no preconditioner, Anderson)"),
                ("AlVac", "Al₄₀ surface (no preconditioner, Anderson)"),
                ("Mn3Si", "Mn₃Si (Kerker, Anderson)")]

        for (case, label) in case
            scfkeys = ["standard_01", "standard_05", "standard_08",
                       "adaptive",    "adaptive_05"]
            case == "AlVac" && filter!(!isequal("adaptive_05"), scfkeys)
            p = plot_convergence(case; colormap, scfkeys, label, showmins=false, kwargs...)
            savefig(p, case * ".pdf")
        end
    end
end

function main()
    dump_plots()
    dump_latex_tables()
end

(abspath(PROGRAM_FILE) == @__FILE__) && main()
