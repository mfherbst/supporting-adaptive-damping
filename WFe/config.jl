using Unitful
using DFTK

lattice = load_lattice(joinpath(@__DIR__, "WFe.json"))
atoms   = load_atoms(  joinpath(@__DIR__, "WFe.json"))
atoms   = attach_pseudos(atoms, Fe="hgh/pbe/fe-q16.hgh", W="hgh/pbe/w-q14.hgh")
magnetic_moments = build_magnetic_moments(atoms; Fe=3.0, W=2.0)

cases = standard_cases()
cases["standard_06"] = (damping=DFTK.FixedDamping(0.6), tol=1e-12, )

config[@__DIR__] = (
    directory        =  @__DIR__,
    lattice          = lattice,
    atoms            = atoms,
    magnetic_moments = magnetic_moments,
    #
    kwargs_model     = (smearing=Smearing.Gaussian(), temperature=0.01),
    Ecut             = 25,
    kgrid            = (7, 7, 1),
    supersampling    = 1.5,
    potential        = nothing,
    kwargs_scf       = (; ),
    mixing           = KerkerMixing(),
    cases            = cases,
)
