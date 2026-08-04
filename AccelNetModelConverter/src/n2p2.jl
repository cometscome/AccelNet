function read_settings(path::AbstractString)
    entries = Dict{String,Vector{Vector{String}}}()
    for raw in eachline(path)
        line = strip(first(split(raw, '#'; limit=2)))
        isempty(line) && continue
        fields = split(line)
        push!(get!(entries, fields[1], Vector{String}[]), fields[2:end])
    end
    return entries
end

values_for(entries, key) = get(entries, key, Vector{String}[])

function single_value(entries, key, default=nothing)
    values = values_for(entries, key)
    if isempty(values)
        isnothing(default) && throw(ConversionError("input.nn is missing keyword '$key'"))
        return string(default)
    end
    length(values) == 1 && length(values[1]) == 1 ||
        throw(ConversionError("input.nn keyword '$key' must have one value"))
    return values[1][1]
end

function parse_symmetry_function(fields::Vector{String}, species::Vector{String})
    length(fields) >= 6 || throw(ConversionError("invalid symfunction_short entry"))
    central = fields[1]
    central in species || throw(ConversionError("unknown central species '$central'"))
    kind = parse(Int, fields[2])
    if kind == 2
        length(fields) == 6 || throw(ConversionError("type-2 symmetry function needs 6 arguments"))
        neighbor = fields[3]
        neighbor in species || throw(ConversionError("unknown neighbor species '$neighbor'"))
        eta, shift, cutoff = parse.(Float64, fields[4:6])
        return SymmetryFunction(central, 2, neighbor, nothing, eta, shift, 0.0, 0.0, cutoff)
    elseif kind in (3, 9)
        length(fields) in (8, 9) ||
            throw(ConversionError("type-$kind symmetry function needs 8 or 9 arguments"))
        neighbor1, neighbor2 = fields[3:4]
        neighbor1 in species && neighbor2 in species ||
            throw(ConversionError("unknown angular neighbor species"))
        eta, lambda, zeta, cutoff = parse.(Float64, fields[5:8])
        shift = length(fields) == 9 ? parse(Float64, fields[9]) : 0.0
        shift == 0.0 || throw(ConversionError("shifted angular functions are not supported"))
        if ATOMIC_NUMBER[neighbor1] > ATOMIC_NUMBER[neighbor2]
            neighbor1, neighbor2 = neighbor2, neighbor1
        end
        return SymmetryFunction(central, kind, neighbor1, neighbor2, eta, shift,
                                lambda, zeta, cutoff)
    end
    throw(ConversionError("only n2p2 symmetry-function types 2, 3, and 9 are convertible"))
end

function read_n2p2_model(directory::AbstractString)
    directory = abspath(directory)
    entries = read_settings(joinpath(directory, "input.nn"))
    declared = parse(Int, single_value(entries, "number_of_elements"))
    elements = values_for(entries, "elements")
    length(elements) == 1 && length(elements[1]) == declared ||
        throw(ConversionError("number_of_elements and elements are inconsistent"))
    all(haskey(ATOMIC_NUMBER, symbol) for symbol in elements[1]) ||
        throw(ConversionError("unknown chemical element in input.nn"))
    species = sort(copy(elements[1]); by=symbol -> ATOMIC_NUMBER[symbol])
    lowercase(single_value(entries, "nnp_type", "2G")) == "2g" ||
        throw(ConversionError("only n2p2 2G models are convertible"))
    cutoff_entries = values_for(entries, "cutoff_type")
    length(cutoff_entries) == 1 && length(only(cutoff_entries)) in (1, 2) ||
        throw(ConversionError("cutoff_type requires a type and optional alpha"))
    cutoff_type = parse(Int, only(cutoff_entries)[1])
    cutoff_alpha = length(only(cutoff_entries)) == 2 ? parse(Float64, only(cutoff_entries)[2]) :
        parse(Float64, single_value(entries, "cutoff_alpha", "0"))
    0 <= cutoff_type <= 9 || throw(ConversionError("cutoff_type must be between 0 and 9"))
    0.0 <= cutoff_alpha < 1.0 || throw(ConversionError("cutoff alpha must satisfy 0 <= alpha < 1"))
    cutoff_type == 9 && cutoff_alpha <= 0.0 &&
        throw(ConversionError("fractional cutoff alpha=h/Rc must be positive"))
    isempty(values_for(entries, "normalize_nodes")) ||
        throw(ConversionError("normalize_nodes is not yet supported by the Julia converter; use the Fortran converter"))
    for key in keys(entries)
        startswith(key, "element_") && endswith(key, "_short") &&
            throw(ConversionError("per-element n2p2 topologies are not yet supported by the Julia converter; use the Fortran converter"))
    end

    hidden_layers = parse(Int, single_value(entries, "global_hidden_layers_short"))
    node_entries = values_for(entries, "global_nodes_short")
    hidden_nodes = isempty(node_entries) ? Int[] : parse.(Int, only(node_entries))
    activation_entries = values_for(entries, "global_activation_short")
    length(activation_entries) == 1 || throw(ConversionError("global_activation_short is required"))
    activations = only.(only(activation_entries))
    length(hidden_nodes) == hidden_layers && length(activations) == hidden_layers + 1 ||
        throw(ConversionError("inconsistent n2p2 network topology"))
    all(haskey(ACTIVATION_TO_CODE, value) for value in activations) ||
        throw(ConversionError("unsupported n2p2 activation function"))

    functions = Dict(symbol => SymmetryFunction[] for symbol in species)
    for fields in values_for(entries, "symfunction_short")
        function_ = parse_symmetry_function(fields, species)
        push!(functions[function_.central], function_)
    end
    for symbol in species
        functions[symbol] = sorted_unique(functions[symbol], n2p2_sort_key)
        isempty(functions[symbol]) && throw(ConversionError("no symmetry functions for $symbol"))
    end
    counts = unique(length(functions[symbol]) for symbol in species)
    length(counts) == 1 || throw(ConversionError("global topology requires equal descriptor counts"))
    nodes = [only(counts); hidden_nodes; 1]

    references = Dict(symbol => 0.0 for symbol in species)
    for values in values_for(entries, "atom_energy")
        length(values) == 2 && haskey(references, values[1]) ||
            throw(ConversionError("invalid atom_energy entry"))
        references[values[1]] = parse(Float64, values[2])
    end

    use_scale = !isempty(values_for(entries, "scale_symmetry_functions"))
    use_center = !isempty(values_for(entries, "center_symmetry_functions"))
    sigma_scale = !isempty(values_for(entries, "scale_symmetry_functions_sigma"))
    sigma_scale && (use_scale || use_center) && throw(ConversionError("invalid scaling combination"))
    scale_min = parse(Float64, single_value(entries, "scale_min_short", "0"))
    scale_max = parse(Float64, single_value(entries, "scale_max_short", "1"))
    minima = Dict(symbol => zeros(length(functions[symbol])) for symbol in species)
    maxima = Dict(symbol => zeros(length(functions[symbol])) for symbol in species)
    shifts = Dict(symbol => zeros(length(functions[symbol])) for symbol in species)
    scales = Dict(symbol => ones(length(functions[symbol])) for symbol in species)
    if use_scale || use_center || sigma_scale
        path = joinpath(directory, "scaling.data")
        isfile(path) || throw(ConversionError("scaling.data is required"))
        seen = Dict(symbol => Set{Int}() for symbol in species)
        for raw in eachline(path)
            line = strip(first(split(raw, '#'; limit=2)))
            isempty(line) && continue
            fields = split(line)
            length(fields) >= 5 || throw(ConversionError("invalid scaling.data row"))
            element_index, sf_index = parse.(Int, fields[1:2])
            1 <= element_index <= length(species) || throw(ConversionError("invalid element index"))
            symbol = species[element_index]
            1 <= sf_index <= length(functions[symbol]) && !(sf_index in seen[symbol]) ||
                throw(ConversionError("invalid or duplicate scaling index"))
            minimum, maximum, mean = parse.(Float64, fields[3:5])
            sigma = length(fields) > 5 ? parse(Float64, fields[6]) : 0.0
            push!(seen[symbol], sf_index)
            minima[symbol][sf_index] = minimum
            maxima[symbol][sf_index] = maximum
            if sigma_scale
                sigma != 0.0 || throw(ConversionError("zero sigma in scaling.data"))
                factor = (scale_max - scale_min) / sigma
                shifts[symbol][sf_index] = mean - scale_min / factor
                scales[symbol][sf_index] = factor
            elseif use_scale
                maximum != minimum || throw(ConversionError("zero range in scaling.data"))
                factor = (scale_max - scale_min) / (maximum - minimum)
                origin = use_center ? mean : minimum
                shifts[symbol][sf_index] = origin - scale_min / factor
                scales[symbol][sf_index] = factor
            else
                shifts[symbol][sf_index] = mean
            end
        end
        all(length(seen[symbol]) == length(functions[symbol]) for symbol in species) ||
            throw(ConversionError("scaling.data is incomplete"))
    end

    weights = Dict{String,Vector{Float64}}()
    expected = sum((nodes[layer] + 1) * nodes[layer + 1] for layer in 1:length(nodes)-1)
    for symbol in species
        path = joinpath(directory, "weights.$(lpad(ATOMIC_NUMBER[symbol], 3, '0')).data")
        isfile(path) || throw(ConversionError("missing $(basename(path))"))
        values = Float64[]
        for raw in eachline(path)
            line = strip(first(split(raw, '#'; limit=2)))
            isempty(line) || push!(values, parse(Float64, first(split(line))))
        end
        length(values) == expected || throw(ConversionError("wrong weight count in $(basename(path))"))
        weights[symbol] = values
    end
    mean_values, conv_values = values_for(entries, "mean_energy"), values_for(entries, "conv_energy")
    isempty(mean_values) == isempty(conv_values) ||
        throw(ConversionError("mean_energy and conv_energy must be used together"))
    mean_energy = isempty(mean_values) ? 0.0 : parse(Float64, only(only(mean_values)))
    conv_energy = isempty(conv_values) ? 1.0 : parse(Float64, only(only(conv_values)))
    conv_energy != 0.0 || throw(ConversionError("conv_energy must be nonzero"))
    return N2P2Model(species, references, nodes, activations, functions, minima,
                     maxima, shifts, scales, weights, mean_energy, conv_energy,
                     cutoff_type, cutoff_alpha)
end
