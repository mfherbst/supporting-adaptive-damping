lattice = load_lattice(joinpath(@__DIR__, "CrFe2Ga_qe.in"))
atoms   = load_atoms(  joinpath(@__DIR__, "CrFe2Ga_qe.in"))
atoms   = attach_pseudos(atoms, Cr="hgh/pbe/cr-q14.hgh",
                                Fe="hgh/pbe/fe-q16.hgh",
                                Ga="hgh/pbe/ga-q3.hgh")
magnetic_moments = build_magnetic_moments(atoms; Ga=0.0, Fe=5.0, Cr=5.0)

cases = standard_cases()
cases["espresso_02"] = (mixing_beta=0.2, type=:qe, )
cases["espresso_03"] = (mixing_beta=0.3, type=:qe, )

for i in 1:8
    cases["adaptive_0$i"] = (damping=DFTK.AdaptiveDamping(i / 10), )
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
    potential        = nothing,
    kwargs_scf       = (; ),
    mixing           = KerkerMixing(),
    cases            = cases,
)
