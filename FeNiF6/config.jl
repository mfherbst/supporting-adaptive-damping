lattice = load_lattice(joinpath(@__DIR__, "FeNiF6.in"))
atoms   = load_atoms(  joinpath(@__DIR__, "FeNiF6.in"))
atoms   = attach_pseudos(atoms, Ni="hgh/pbe/ni-q18.hgh",
                                Fe="hgh/pbe/fe-q16.hgh",
                                 F="hgh/pbe/f-q7.hgh")
magnetic_moments = build_magnetic_moments(atoms; F=0.0, Fe=5.0, Ni=5.0)

cases = standard_cases()
cases["adaptive_slowstart"] = (damping=DFTK.AdaptiveDamping(;α_trial_init=0.5), )

config[@__DIR__] = (
    directory        =  @__DIR__,
    lattice          = lattice,
    atoms            = atoms,
    magnetic_moments = magnetic_moments,
    #
    kwargs_model     = (smearing=Smearing.Gaussian(), temperature=0.01),
    Ecut             = 45,
    kgrid            = [11, 11, 11],
    supersampling    = 2.0,
    potential        = nothing,
    kwargs_scf       = (; ),
    mixing           = KerkerMixing(),
    cases            = cases,
)
