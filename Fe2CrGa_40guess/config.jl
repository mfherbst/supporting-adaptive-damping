lattice = load_lattice(joinpath(@__DIR__, "CrFe2Ga_qe.in"))
atoms   = load_atoms(  joinpath(@__DIR__, "CrFe2Ga_qe.in"))
atoms   = attach_pseudos(atoms, Cr="hgh/pbe/cr-q14.hgh",
                                Fe="hgh/pbe/fe-q16.hgh",
                                Ga="hgh/pbe/ga-q3.hgh")
magnetic_moments = build_magnetic_moments(atoms; Ga=0.0, Fe=5.0, Cr=5.0)

cases = Dict(
    "maxcond_1e2" => (damping=DFTK.FixedDamping(0.3), diagtol_max=1e-8,
                      acceleration=DFTK.AndersonAcceleration(; maxcond=1e2)),
    "maxcond_1e4" => (damping=DFTK.FixedDamping(0.3), diagtol_max=1e-8,
                      acceleration=DFTK.AndersonAcceleration(; maxcond=1e4)),
    "maxcond_1e6" => (damping=DFTK.FixedDamping(0.3), diagtol_max=1e-8,
                      acceleration=DFTK.AndersonAcceleration(; maxcond=1e6)),
    "maxcond_1e8" => (damping=DFTK.FixedDamping(0.3), diagtol_max=1e-8,
                      acceleration=DFTK.AndersonAcceleration(; maxcond=1e8)),
    "maxcond_1e10" => (damping=DFTK.FixedDamping(0.3), diagtol_max=1e-8,
                       acceleration=DFTK.AndersonAcceleration(; maxcond=1e10)),
)

function load_Fe2CrGa_potential_40(basis)
    guessfile = joinpath(@__DIR__, "Fe2CrGa40_initV.jld2")
    if !isfile(guessfile)
        ρ, ψ0, V = load_guess(basis, joinpath(@__DIR__, "../Fe2CrGa/init.jld2"))
        resultsfile = joinpath(@__DIR__, "init.scfres.json")

        scfres = run_dftk(basis, resultsfile; ρ=ρ, ψ0=ψ0, V=V,
                          mixing=KerkerMixing(), damping=DFTK.FixedDamping(0.3),
                          acceleration=DFTK.AndersonAcceleration(; maxcond=1e6), maxiter=40)
        V = DFTK.total_local_potential(scfres.ham)
        jldopen(guessfile, "w") do file
            file["V"] = V
        end
    end
    jldopen(file -> file["V"], guessfile, "r")
end


config[@__DIR__] = (
    directory        =  @__DIR__,
    lattice          = lattice,
    atoms            = atoms,
    magnetic_moments = magnetic_moments,
    #
    kwargs_model     = (smearing=Smearing.Gaussian(), temperature=0.01),
    Ecut             = 45,
    kgrid            = [13, 13, 13],
    supersampling    = 2.0,
    potential        = load_Fe2CrGa_potential_40,
    kwargs_scf       = (; ),
    mixing           = KerkerMixing(),
    cases            = cases,
)
