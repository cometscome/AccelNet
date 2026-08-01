function n2p2_to_accelnet(input_directory::AbstractString, output_directory::AbstractString)
    model = read_n2p2_model(input_directory)
    mkpath(output_directory)
    written = String[]
    for symbol in model.species
        source = model.functions[symbol]
        target = sorted_unique(source, aenet_sort_key)
        network = AtomicNetwork(symbol, "Converted n2p2 2G-HDNNP model", model.species,
            [model.atomic_references[item] for item in model.species], model.nodes,
            [ACTIVATION_TO_CODE[value] for value in model.activations],
            reorder_first_layer(model.weights[symbol], model.nodes, source, target), source,
            reorder_values(model.minima[symbol], source, target),
            reorder_values(model.maxima[symbol], source, target),
            reorder_values(model.shifts[symbol], source, target),
            reorder_values(model.scales[symbol], source, target), 0.1,
            maximum(function_.cutoff for function_ in source), model.conv_energy, model.mean_energy)
        path = joinpath(output_directory, "$symbol.nn.ascii")
        write_atomic_network(path, network)
        push!(written, path)
    end
    write(joinpath(output_directory, "networks.list"),
          join(("$symbol $symbol.nn.ascii" for symbol in model.species), "\n") * "\n")
    return written
end

function validate_network_set(networks::Vector{AtomicNetwork})
    isempty(networks) && throw(ConversionError("at least one AccelNet network is required"))
    first_network = networks[1]
    Set(network.atomtype for network in networks) == Set(first_network.species) &&
        length(networks) == length(first_network.species) ||
        throw(ConversionError("one network is required for every embedded species"))
    for network in networks[2:end]
        network.species == first_network.species &&
        network.atomic_references == first_network.atomic_references &&
        network.nodes == first_network.nodes && network.activations == first_network.activations &&
        network.energy_scale == first_network.energy_scale &&
        network.energy_shift == first_network.energy_shift ||
            throw(ConversionError("AccelNet networks have inconsistent global metadata"))
    end
    return sort(copy(first_network.species); by=symbol -> ATOMIC_NUMBER[symbol])
end

function write_n2p2_weights(path, weights, nodes)
    lines = ["# Neural network connections converted from AccelNet."]
    index = 1; connection = 1
    for layer in 1:length(nodes)-1
        for source in 1:nodes[layer], target in 1:nodes[layer+1]
            push!(lines, @sprintf("%.17e a %d %d %d %d %d", weights[index], connection,
                                  layer-1, source, layer, target))
            index += 1; connection += 1
        end
        for target in 1:nodes[layer+1]
            push!(lines, @sprintf("%.17e b %d %d %d", weights[index], connection, layer, target))
            index += 1; connection += 1
        end
    end
    write(path, join(lines, "\n") * "\n")
end

function accelnet_to_n2p2(network_paths::Vector{<:AbstractString}, output_directory::AbstractString)
    networks = read_atomic_network.(network_paths)
    species = validate_network_set(networks)
    by_atom = Dict(network.atomtype => network for network in networks)
    first_network = networks[1]
    mkpath(output_directory)
    settings = ["# n2p2 2G-HDNNP model converted from AccelNet/ænet networks",
                "number_of_elements $(length(species))", "elements $(join(species, ' '))"]
    reference = Dict(zip(first_network.species, first_network.atomic_references))
    append!(settings, [@sprintf("atom_energy %s %.17e", symbol, reference[symbol]) for symbol in species])
    if first_network.energy_scale != 1.0 || first_network.energy_shift != 0.0
        append!(settings, [@sprintf("mean_energy %.17e", first_network.energy_shift),
                           @sprintf("conv_energy %.17e", first_network.energy_scale), "conv_length 1.0"])
    end
    append!(settings, ["cutoff_type 1", "cutoff_alpha 0.0", "scale_symmetry_functions_sigma",
        "scale_min_short 0.0", "scale_max_short 1.0",
        "global_hidden_layers_short $(length(first_network.nodes)-2)",
        "global_nodes_short $(join(first_network.nodes[2:end-1], ' '))",
        "global_activation_short $(join((CODE_TO_ACTIVATION[value] for value in first_network.activations), ' '))", ""])
    scaling = ["# e_index sf_index sf_min sf_max sf_mean sf_sigma"]
    written = String[]
    for (element_index, symbol) in enumerate(species)
        network = by_atom[symbol]
        source = sorted_unique(network.functions, aenet_sort_key)
        target = sorted_unique(network.functions, n2p2_sort_key)
        target_minima = reorder_values(network.descriptor_minimum, source, target)
        target_maxima = reorder_values(network.descriptor_maximum, source, target)
        target_shifts = reorder_values(network.descriptor_shift, source, target)
        target_scales = reorder_values(network.descriptor_scale, source, target)
        target_weights = reorder_first_layer(network.weights, network.nodes, source, target)
        for function_ in target
            if function_.kind == 2
                push!(settings, @sprintf("symfunction_short %s 2 %s %.17e %.17e %.17e",
                    symbol, function_.neighbor1, function_.eta, function_.shift, function_.cutoff))
            else
                push!(settings, @sprintf("symfunction_short %s %d %s %s %.17e %.17e %.17e %.17e 0.0",
                    symbol, function_.kind, function_.neighbor1, function_.neighbor2,
                    function_.eta, function_.lambda, function_.zeta, function_.cutoff))
            end
        end
        for index in eachindex(target)
            scale = target_scales[index]
            scale != 0.0 || throw(ConversionError("zero descriptor scale cannot be represented"))
            sigma = inv(scale); minimum = target_minima[index]; maximum = target_maxima[index]
            if minimum == maximum
                minimum = target_shifts[index] - 1.0e6 * abs(sigma)
                maximum = target_shifts[index] + 1.0e6 * abs(sigma)
            end
            push!(scaling, @sprintf("%d %d %.17e %.17e %.17e %.17e", element_index,
                index, minimum, maximum, target_shifts[index], sigma))
        end
        path = joinpath(output_directory, "weights.$(lpad(ATOMIC_NUMBER[symbol], 3, '0')).data")
        write_n2p2_weights(path, target_weights, network.nodes)
        push!(written, path)
    end
    input_path = joinpath(output_directory, "input.nn")
    scaling_path = joinpath(output_directory, "scaling.data")
    write(input_path, join(settings, "\n") * "\n")
    write(scaling_path, join(scaling, "\n") * "\n")
    return [input_path, scaling_path, written...]
end
