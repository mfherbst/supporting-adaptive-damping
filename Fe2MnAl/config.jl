lattice = load_lattice(joinpath(@__DIR__, "./AlFe2Mn_qe.in"))
atoms   = load_atoms(  joinpath(@__DIR__, "./AlFe2Mn_qe.in"))
atoms   = attach_pseudos(atoms, Mn="hgh/pbe/mn-q15.hgh",
                                Fe="hgh/pbe/fe-q16.hgh",
                                Al="hgh/pbe/al-q3.hgh")
magnetic_moments = build_magnetic_moments(atoms; Al=0.0, Fe=5.0, Mn=5.0)


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
    potential        = nothing,
    kwargs_scf       = (; ),
    mixing           = KerkerMixing(),
    cases            = standard_cases(),
)
