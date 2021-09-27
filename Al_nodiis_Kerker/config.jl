lattice = load_lattice(joinpath(@__DIR__, "Al8.json"))
atoms   = load_atoms(  joinpath(@__DIR__, "Al8.json"))
atoms   = attach_pseudos(atoms, Al="hgh/pbe/al-q3.hgh")

function load_Al_nodiis_potential(basis)
    guessfile = joinpath(@__DIR__, "Al_nodiis_initV.jld2")
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


cases = Dict{String, Any}(
    "adaptive" => (damping=DFTK.AdaptiveDamping(), ),  # Parameter-free adaptive damping
)
for α in 1:9
    cases["standard_0$α"] = (damping=DFTK.FixedDamping(α / 10), )
end
cases["standard_10"] = (damping=DFTK.FixedDamping(1.0), )
for i in 3:9
    cases["adaptive_0$i"] = (damping=DFTK.AdaptiveDamping(i / 10), )
end

config[@__DIR__] = (
    directory        =  @__DIR__,
    lattice          = lattice,
    atoms            = atoms,
    magnetic_moments = [],
    #
    kwargs_model     = (smearing=Smearing.Gaussian(), temperature=0.001),
    Ecut             = 20,
    kgrid            = kgrid_size_from_minimal_spacing(lattice, 2π * 0.04 / u"Å"),
    supersampling    = 1.5,
    potential        = load_Al_nodiis_potential,
    kwargs_scf       = (acceleration=DFTK.AndersonAcceleration(;m=0), ),
    mixing           = KerkerMixing(),
    cases            = cases,
)

