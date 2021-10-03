lattice = load_lattice(joinpath(@__DIR__, "./Mn3Si.in"))
atoms   = load_atoms(  joinpath(@__DIR__, "./Mn3Si.in"))
atoms   = attach_pseudos(atoms, Mn="hgh/pbe/mn-q15.hgh",
                                Si="hgh/pbe/si-q4.hgh")

# AFM ordering, configuration 5(c) in Jiang and Yang J. Alloys Compounds (2021)
@assert atoms[1][1].symbol == :Mn
@assert atoms[2][1].symbol == :Si
magnetic_moments = [atoms[1][1] => [5.0, -5.0, 5.0], atoms[2][1] => [0.0]]

cases = standard_cases()

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
