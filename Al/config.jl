lattice = load_lattice(joinpath(@__DIR__, "Al40.json"))
atoms   = load_atoms(  joinpath(@__DIR__, "Al40.json"))
atoms   = attach_pseudos(atoms, Al="hgh/pbe/al-q3.hgh")

cases = standard_cases()
cases["adaptive_05"] = (damping=DFTK.AdaptiveDamping(;α_trial_min=0.5, α_min=0.05), )

function load_Al40_potential(basis)
    guessfile = joinpath(@__DIR__, "Al40_initV.jld2")
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


config[@__DIR__] = (
    directory        =  @__DIR__,
    lattice          = lattice,
    atoms            = atoms,
    magnetic_moments = [],
    #
    kwargs_model     = (smearing=Smearing.Gaussian(), temperature=0.001),
    Ecut             = 20,
    kgrid            = [3, 3, 1],
    supersampling    = 1.5,
    potential        = load_Al40_potential,
    kwargs_scf       = (; ),
    mixing           = SimpleMixing(),
    cases            = cases,
)

