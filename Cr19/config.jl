using Unitful
using DFTK

lattice = load_lattice(joinpath(@__DIR__, "Cr19.json"))
atoms   = load_atoms(  joinpath(@__DIR__, "Cr19.json"))
atoms   = attach_pseudos(atoms, Cr="hgh/pbe/cr-q14.hgh")
magnetic_moments = build_magnetic_moments(atoms; Cr=5.0)

cases = standard_cases()
cases["standard_02"] = (damping=DFTK.FixedDamping(0.2), tol=1e-12, )
cases["standard_03"] = (damping=DFTK.FixedDamping(0.3), tol=1e-12, )
cases["standard_04"] = (damping=DFTK.FixedDamping(0.4), tol=1e-12, )

config[@__DIR__] = (
    directory        =  @__DIR__,
    lattice          = lattice,
    atoms            = atoms,
    magnetic_moments = magnetic_moments,
    #
    kwargs_model     = (smearing=Smearing.Gaussian(), temperature=0.01),
    Ecut             = 30,
    kgrid            = kgrid_size_from_minimal_spacing(lattice, 2π * 0.04 / u"Å"),
    supersampling    = 1.5,
    potential        = nothing,
    kwargs_scf       = (; ),
    mixing           = KerkerMixing(),
    cases            = cases,
)
