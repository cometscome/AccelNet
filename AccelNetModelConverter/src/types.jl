struct ConversionError <: Exception
    message::String
end
Base.showerror(io::IO, error::ConversionError) = print(io, error.message)

const PERIODIC_TABLE = split("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og")
const ATOMIC_NUMBER = Dict(symbol => index for (index, symbol) in enumerate(PERIODIC_TABLE))
const ACTIVATION_TO_CODE = Dict('l'=>0, 't'=>1, 's'=>2, 'p'=>3, 'r'=>5,
                                'g'=>6, 'c'=>7, 'S'=>8, 'e'=>9, 'h'=>10)
const CODE_TO_ACTIVATION = Dict(value => key for (key, value) in ACTIVATION_TO_CODE)

struct SymmetryFunction
    central::String
    kind::Int
    neighbor1::String
    neighbor2::Union{Nothing,String}
    eta::Float64
    shift::Float64
    lambda::Float64
    zeta::Float64
    cutoff::Float64
end

identity(function_::SymmetryFunction) =
    (function_.kind, function_.neighbor1, function_.neighbor2, function_.eta,
     function_.shift, function_.lambda, function_.zeta, function_.cutoff)

n2p2_sort_key(function_::SymmetryFunction) =
    (function_.kind, function_.cutoff, function_.eta, function_.shift,
     function_.zeta, function_.lambda, ATOMIC_NUMBER[function_.neighbor1],
     isnothing(function_.neighbor2) ? 0 : ATOMIC_NUMBER[function_.neighbor2])

function aenet_sort_key(function_::SymmetryFunction)
    first = ATOMIC_NUMBER[function_.neighbor1]
    second = isnothing(function_.neighbor2) ? 0 : ATOMIC_NUMBER[function_.neighbor2]
    function_.kind == 2 && return (first, 0, 0, 0, n2p2_sort_key(function_))
    return (min(first, second), 1, max(first, second), function_.kind == 3 ? 0 : 1,
            n2p2_sort_key(function_))
end

function sorted_unique(functions::Vector{SymmetryFunction}, by)
    result = sort(copy(functions); by)
    identities = identity.(result)
    length(unique(identities)) == length(identities) ||
        throw(ConversionError("duplicate symmetry functions cannot be converted safely"))
    return result
end

function reorder_values(values::Vector{Float64}, old::Vector{SymmetryFunction},
                        new::Vector{SymmetryFunction})
    length(values) == length(old) || throw(ConversionError("descriptor array length mismatch"))
    positions = Dict(identity(function_) => index for (index, function_) in enumerate(old))
    try
        return [values[positions[identity(function_)]] for function_ in new]
    catch error
        error isa KeyError || rethrow()
        throw(ConversionError("symmetry-function sets differ during reordering"))
    end
end

function reorder_first_layer(weights::Vector{Float64}, nodes::Vector{Int},
                             old::Vector{SymmetryFunction}, new::Vector{SymmetryFunction})
    nodes[1] == length(old) == length(new) ||
        throw(ConversionError("first-layer size does not match symmetry functions"))
    nout = nodes[2]
    first_layer_size = (nodes[1] + 1) * nout
    length(weights) >= first_layer_size || throw(ConversionError("network weights are too short"))
    positions = Dict(identity(function_) => index for (index, function_) in enumerate(old))
    result = Float64[]
    for function_ in new
        old_index = get(positions, identity(function_), 0)
        old_index > 0 || throw(ConversionError("symmetry-function sets differ during weight reordering"))
        start = (old_index - 1) * nout + 1
        append!(result, @view weights[start:start+nout-1])
    end
    bias_start = nodes[1] * nout + 1
    append!(result, @view weights[bias_start:first_layer_size])
    append!(result, @view weights[first_layer_size+1:end])
    return result
end

struct AtomicNetwork
    atomtype::String
    description::String
    species::Vector{String}
    atomic_references::Vector{Float64}
    nodes::Vector{Int}
    activations::Vector{Int}
    weights::Vector{Float64}
    functions::Vector{SymmetryFunction}
    descriptor_minimum::Vector{Float64}
    descriptor_maximum::Vector{Float64}
    descriptor_shift::Vector{Float64}
    descriptor_scale::Vector{Float64}
    minimum_radius::Float64
    maximum_radius::Float64
    energy_scale::Float64
    energy_shift::Float64
end

struct N2P2Model
    species::Vector{String}
    atomic_references::Dict{String,Float64}
    nodes::Vector{Int}
    activations::Vector{Char}
    functions::Dict{String,Vector{SymmetryFunction}}
    minima::Dict{String,Vector{Float64}}
    maxima::Dict{String,Vector{Float64}}
    shifts::Dict{String,Vector{Float64}}
    scales::Dict{String,Vector{Float64}}
    weights::Dict{String,Vector{Float64}}
    mean_energy::Float64
    conv_energy::Float64
end
