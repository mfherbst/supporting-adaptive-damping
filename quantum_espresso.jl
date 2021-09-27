using MPI
using MPICH_jll
using PyCall
using QuantumEspresso_jll
using Unitful
using UnitfulAtomic

const QE_PSEUDODIR = joinpath(@__DIR__, "qe_pseudos")
const QE_PSEUDOMAP = Dict(
   "Si" => "Si.pbe-hgh.UPF",
   "O"  => "O.pbe-hgh.UPF",
   "Al" => "Al.pbe-hgh.UPF",
   "Cu" => "Cu.pbe-d-hgh.UPF",
   "Ga" => "Ga.pbe-d-hgh.UPF",
   "As" => "As.pbe-hgh.UPF",
   "H"  => "H.pbe-hgh.UPF",
   "Cr" => "Cr.pbe-sp-hgh.UPF",
   "Ni" => "Ni.pbe-sp-hgh.UPF",
   "Co" => "Co.pbe-sp-hgh.UPF",
   "Fe" => "Fe.pbe-sp-hgh.UPF",
   "Mn" => "Mn.pbe-sp-hgh.UPF",
   "Ge" => "Ge.pbe-hgh.UPF",
   "Pt" => "Pt.pbe-hgh.UPF",
   "B"  => "B.pbe-hgh.UPF",
   "Hf" => "Hf.pbe-sp-hgh.UPF",
)


# TODO This should be transformed into a proper function inside DFTK
#      ... a better interface would be based on a `basis` or a `model`
"""Ecut is in Hartree, tol is energy tolerance in Hartree, kspacing is in inverse Bohrs"""
function scf_quantum_espresso(lattice, atoms;
                              Ecut=20, tol=1e-11, mixing_mode="plain",
                              fileprefix="espresso", smearing="gaussian",
                              temperature=0, mixing_beta=0.7, mixing_ndim=10,
                              maxiter=100, kgrid=nothing, magnetic_moments=nothing,
                              input_dft="pbe", pseudomap=QE_PSEUDOMAP,
                              pseudodir=QE_PSEUDODIR, kwargs...)
    if !mpi_master()
        error("Running QE does not support MPI parallelisation")
    end


    # TODO Pseudos are just hardcoded for now
    #      ... better look at `atoms` and map to the HGH QE files
    if isnothing(magnetic_moments) || isempty(magnetic_moments)
        system = ase_atoms(lattice, atoms)
    else
        system = ase_atoms(lattice, atoms, magnetic_moments)
    end

    # Units used by ASE and QE
    Ry = 0.5
    eV = austrip(u"1eV")
    system.calc = pyimport("ase.calculators.espresso").Espresso(;
        label=fileprefix,
        input_dft=input_dft,
        pseudopotentials=pseudomap,
        kpts=kgrid,
        ecutwfc=Ecut / Ry,  # QE uses Rydbergs
        tstress=true,       # Compute stresses
        tprnfor=true,       # Compute forces
        mixing_mode=mixing_mode,
        mixing_beta=mixing_beta,
        conv_thr=tol / Ry,
        occupations=(temperature > 0 ? "smearing" : "fixed"),
        smearing=smearing,
        degauss=temperature / Ry,
        electron_maxstep=maxiter,
        mixing_ndim=mixing_ndim,
        kwargs...
    )

    results = Dict{Symbol, Any}(
        :inputfile  => fileprefix * ".pwi",
        :outputfile => fileprefix * ".pwo",
    )
    MPICH_jll.mpiexec() do mpirun
        QuantumEspresso_jll.pwscf() do pwscf
            qe_command = "$mpirun -np $(mpi_nprocs()) $pwscf -in PREFIX.pwi > PREFIX.pwo"
            @info "Running QuantumEspresso: $qe_command"
            system.calc.command = qe_command
            withenv("ESPRESSO_PSEUDO" => pseudodir,
                    "OMP_NUM_THREADS" => BLAS.get_num_threads()) do
                results[:energies_total] = system.get_potential_energy() * eV
            end
        end
    end

    (; results...)
end


function parse_quantum_espresso_pwo(file)
    Ry = 0.5  # Ry -> Hartree

    iterations = []
    current_iter = nothing

    for line in readlines(file)
        miter = match(r"^ *iteration # *([0-9]+)", line)
        if miter !== nothing  # New iteration starts
            if !isnothing(current_iter)
                push!(iterations, (; current_iter...))
            end
            current_iter = Dict{Symbol, Any}(:n_iter => parse(Int, miter[1]))
        end

        menergy = match(r"^ *total energy *= *([0-9.eE+-]+) *Ry", line)
        if menergy !== nothing
            current_iter[:total_energy] = Ry * parse(Float64, menergy[1])
        end

        maccuracy = match(r"^ *estimated scf accuracy *< *([0-9.eE+-]+) *Ry", line)
        if maccuracy !== nothing
            current_iter[:accuracy] = Ry * parse(Float64, maccuracy[1])
        end
    end

    iterations
end
