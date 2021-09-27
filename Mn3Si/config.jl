lattice = load_lattice(joinpath(@__DIR__, "./Mn3Si.in"))
atoms   = load_atoms(  joinpath(@__DIR__, "./Mn3Si.in"))
atoms   = attach_pseudos(atoms, Mn="hgh/pbe/mn-q15.hgh",
                                Si="hgh/pbe/si-q4.hgh")
magnetic_moments = build_magnetic_moments(atoms; Si=0.0, Mn=5.0)

cases = standard_cases()
for i in [3, 4, 6]
    cases["adaptive_0$i"] = (damping=DFTK.AdaptiveDamping(i / 10), )
end

config[@__DIR__] = (
    directory        = @__DIR__,
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
    cases            = cases,
)
