lattice = load_lattice(joinpath(@__DIR__, "GaAs20.json"))
atoms   = load_atoms(  joinpath(@__DIR__, "GaAs20.json"))
atoms   = attach_pseudos(atoms, Ga="hgh/pbe/ga-q3.hgh", As="hgh/pbe/as-q5.hgh")

cases = standard_cases()
cases["adaptive_03"] = (damping=DFTK.AdaptiveDamping(;α_trial_min=0.3, α_min=0.05), tol=1e-12)
cases["adaptive_04"] = (damping=DFTK.AdaptiveDamping(;α_trial_min=0.4, α_min=0.05), tol=1e-12)
cases["adaptive_05"] = (damping=DFTK.AdaptiveDamping(;α_trial_min=0.5, α_min=0.05), tol=1e-12)

function load_GaAs_20_potential(basis)
    guessfile = joinpath(@__DIR__, "GaAs20_initV.jld2")
    if !isfile(guessfile)
        ρ0 = guess_density(basis)
        _, ham = DFTK.energy_hamiltonian(basis, nothing, nothing; ρ=ρ0)
        V = DFTK.total_local_potential(ham)
        V .+= 3.0rand(Float64, size(V)...)

        jldopen(guessfile, "w") do file
            file["V"] = V
        end
    end

    jldopen(guessfile, "r") do file
        file["V"]
    end
end

# Load guess from file
config[@__DIR__] = (
    directory        = @__DIR__,
    lattice          = lattice,
    atoms            = atoms,
    magnetic_moments = [],
    #
    kwargs_model     = (smearing=Smearing.None(), temperature=0.0),
    Ecut             = 20,
    kgrid            = [4, 4, 1],
    supersampling    = 1.5,
    potential        = load_GaAs_20_potential,
    kwargs_scf       = (; ),
    mixing           = SimpleMixing(),
    cases            = cases,
)
