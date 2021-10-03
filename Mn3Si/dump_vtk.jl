using DFTK
using WriteVTK

function build_magnetic_moments(atoms::AbstractArray; magmoms...)
    magmoms = Dict(magmoms)
    map(atoms) do (element, positions)
        element => fill(magmoms[element.symbol], length(positions))
    end
end

function attach_pseudos(atoms::AbstractArray; pseudomap...)
    pseudomap = Dict(pseudomap)

    map(atoms) do (element, position)
        pspfile = get(pseudomap, element.symbol, nothing)
        ElementPsp(element.symbol, psp=load_psp(pspfile)) => position
    end
end


lattice = load_lattice(joinpath(@__DIR__, "./Mn3Si.in"))
atoms   = load_atoms(  joinpath(@__DIR__, "./Mn3Si.in"))
atoms   = attach_pseudos(atoms, Mn="hgh/pbe/mn-q15.hgh",
                                Si="hgh/pbe/si-q4.hgh")
magnetic_moments = build_magnetic_moments(atoms; Si=0.0, Mn=5.0)

model = model_PBE(lattice, atoms; magnetic_moments=magnetic_moments,
                  smearing=Smearing.Gaussian(), temperature=0.01)
basis = PlaneWaveBasis(model; Ecut=45, kgrid=[13, 13, 13])

determine_diagtol = DFTK.ScfDiagtol(; ratio_ρdiff=0.01, diagtol_max=1e-3, diagtol_min=1e-12)
ρ = guess_density(basis, magnetic_moments)
scfres = DFTK.scf_potential_mixing(basis; tol=1e-11, determine_diagtol=determine_diagtol,
                                   mixing=KerkerMixing(), damping=DFTK.FixedDamping(0.7),
                                   diag_miniter=2, ρ=ρ)

save_scfres("Mn3Si.vts", scfres; save_ψ=false);
