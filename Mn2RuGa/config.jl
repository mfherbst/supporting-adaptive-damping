lattice = load_lattice(joinpath(@__DIR__, "./Mn2RuGa.in"))
atoms   = load_atoms(  joinpath(@__DIR__, "./Mn2RuGa.in"))
atoms   = attach_pseudos(atoms, Mn="hgh/pbe/mn-q15.hgh",
                                Ru="hgh/pbe/ru-q16.hgh",
                                Ga="hgh/pbe/ga-q3.hgh")
magnetic_moments = build_magnetic_moments(atoms; Ga=0.0, Mn=5.0, Ru=5.0)


config[@__DIR__] = (
    directory        =  @__DIR__,
    lattice          = lattice,
    atoms            = atoms,
    magnetic_moments = magnetic_moments,
    #
    kwargs_model     = (smearing=Smearing.Gaussian(), temperature=0.01),
    Ecut             = 35,
    kgrid            = [13, 13, 13],
    supersampling    = 2.0,
    potential        = nothing,
    kwargs_scf       = (; ),
    mixing           = KerkerMixing(),
    cases            = standard_cases(),
)
