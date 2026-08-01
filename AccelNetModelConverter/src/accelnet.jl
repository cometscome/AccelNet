mutable struct TextRecordReader
    path::String
    lines::Vector{String}
    position::Int
end
TextRecordReader(path) = TextRecordReader(path, readlines(path), 1)

function nextline!(reader::TextRecordReader)
    reader.position <= length(reader.lines) || throw(ConversionError("unexpected end of $(basename(reader.path))"))
    value = reader.lines[reader.position]
    reader.position += 1
    return value
end

function tokens!(reader::TextRecordReader, count::Int)
    result = String[]
    while length(result) < count
        append!(result, split(nextline!(reader)))
    end
    length(result) == count || throw(ConversionError("unexpected record length in $(basename(reader.path))"))
    return result
end

mutable struct BinaryRecordReader
    path::String
    data::Vector{UInt8}
    position::Int
    swap::Bool
end

native_int32(bytes) = reinterpret(Int32, bytes)[1]

function BinaryRecordReader(path)
    data = read(path)
    length(data) >= 12 || throw(ConversionError("unsupported Fortran file $(basename(path))"))
    for swap in (false, true)
        raw = native_int32(data[1:4])
        size = Int(swap ? bswap(raw) : raw)
        if 0 < size <= length(data) - 8
            trailing = native_int32(data[5+size:8+size])
            trailing = swap ? bswap(trailing) : trailing
            trailing == size && return BinaryRecordReader(path, data, 1, swap)
        end
    end
    throw(ConversionError("$(basename(path)) is not a supported Fortran sequential file"))
end

function record!(reader::BinaryRecordReader)
    reader.position + 7 <= length(reader.data) || throw(ConversionError("unexpected binary EOF"))
    raw = native_int32(reader.data[reader.position:reader.position+3])
    size = Int(reader.swap ? bswap(raw) : raw)
    start = reader.position + 4
    stop = start + size - 1
    stop + 4 <= length(reader.data) || throw(ConversionError("truncated Fortran record"))
    trailing = native_int32(reader.data[stop+1:stop+4])
    trailing = reader.swap ? bswap(trailing) : trailing
    trailing == size || throw(ConversionError("Fortran record marker mismatch"))
    reader.position = stop + 5
    return reader.data[start:stop]
end

function integers!(reader::BinaryRecordReader, count)
    bytes = record!(reader)
    length(bytes) == 4count || throw(ConversionError("unexpected binary integer record"))
    values = collect(reinterpret(Int32, bytes))
    reader.swap && (values = bswap.(values))
    return Int.(values)
end
integer!(reader) = only(integers!(reader, 1))

function reals!(reader::BinaryRecordReader, count)
    bytes = record!(reader)
    length(bytes) == 8count || throw(ConversionError("unexpected binary real record"))
    if reader.swap
        words = bswap.(collect(reinterpret(UInt64, bytes)))
        return reinterpret(Float64, words) |> collect
    end
    return collect(reinterpret(Float64, bytes))
end
real!(reader) = only(reals!(reader, 1))

function characters!(reader::BinaryRecordReader, count)
    bytes = record!(reader)
    length(bytes) % count == 0 || throw(ConversionError("unexpected binary character record"))
    width = length(bytes) ÷ count
    return [rstrip(String(bytes[(index-1)*width+1:index*width]), [' ', '\0']) for index in 1:count]
end
character!(reader) = only(characters!(reader, 1))

function atomic_network_from_components(path, nlayers, maxnodes, nweights, nodes,
        activations, weights, description, atomtype, environment_names,
        minimum_radius, maximum_radius, descriptor_name, nsf, nparam, kinds,
        raw_parameters, raw_environments, minima, maxima, averages, moments,
        energy_scale, energy_shift, species, atomic_references)
    lowercase(strip(descriptor_name)) == "behler2011" ||
        throw(ConversionError("$(basename(path)) uses $descriptor_name; n2p2 requires Behler2011"))
    maxnodes == maximum(nodes) && nodes[1] == nsf || throw(ConversionError("inconsistent topology"))
    expected = sum((nodes[layer] + 1) * nodes[layer+1] for layer in 1:nlayers-1)
    expected == nweights || throw(ConversionError("inconsistent weight count"))
    all(haskey(CODE_TO_ACTIVATION, code) for code in activations) ||
        throw(ConversionError("activation has no n2p2 equivalent"))
    functions = SymmetryFunction[]
    for index in 1:nsf
        parameter = raw_parameters[(index-1)*nparam+1:index*nparam]
        environment = raw_environments[(index-1)*2+1:index*2]
        1 <= environment[1] <= length(environment_names) || throw(ConversionError("invalid environment"))
        neighbor1 = environment_names[environment[1]]
        if kinds[index] == 2
            push!(functions, SymmetryFunction(atomtype, 2, neighbor1, nothing,
                                              parameter[3], parameter[2], 0.0, 0.0, parameter[1]))
        elseif kinds[index] in (4, 5)
            1 <= environment[2] <= length(environment_names) || throw(ConversionError("invalid angular environment"))
            neighbor2 = environment_names[environment[2]]
            if ATOMIC_NUMBER[neighbor1] > ATOMIC_NUMBER[neighbor2]
                neighbor1, neighbor2 = neighbor2, neighbor1
            end
            push!(functions, SymmetryFunction(atomtype, kinds[index] == 4 ? 3 : 9,
                                              neighbor1, neighbor2, parameter[4], 0.0,
                                              parameter[2], parameter[3], parameter[1]))
        else
            throw(ConversionError("descriptor kind $(kinds[index]) has no n2p2 equivalent"))
        end
    end
    scales = [begin
        variance = max(moment - average^2, 0.0)
        variance > 0 ? inv(sqrt(variance)) : 1.0
    end for (average, moment) in zip(averages, moments)]
    return AtomicNetwork(atomtype, description, species, atomic_references, nodes,
                         activations, weights, functions, minima, maxima, averages,
                         scales, minimum_radius, maximum_radius, energy_scale, energy_shift)
end

function read_ascii_atomic_network(path)
    reader = TextRecordReader(path)
    nlayers = parse(Int, nextline!(reader)); maxnodes = parse(Int, nextline!(reader))
    nweights = parse(Int, nextline!(reader)); nextline!(reader)
    nodes = parse.(Int, tokens!(reader, nlayers))
    activations = parse.(Int, tokens!(reader, nlayers-1))
    tokens!(reader, nlayers); tokens!(reader, nlayers)
    weights = parse.(Float64, tokens!(reader, nweights))
    description = strip(nextline!(reader)); atomtype = strip(nextline!(reader))
    nenv = parse(Int, nextline!(reader)); environment_names = tokens!(reader, nenv)
    minimum_radius = parse(Float64, nextline!(reader)); maximum_radius = parse(Float64, nextline!(reader))
    descriptor_name = strip(nextline!(reader)); nsf = parse(Int, nextline!(reader)); nparam = parse(Int, nextline!(reader))
    kinds = parse.(Int, tokens!(reader, nsf))
    parameters = parse.(Float64, tokens!(reader, nparam*nsf))
    environments = parse.(Int, tokens!(reader, 2nsf)); nextline!(reader)
    minima = parse.(Float64, tokens!(reader, nsf)); maxima = parse.(Float64, tokens!(reader, nsf))
    averages = parse.(Float64, tokens!(reader, nsf)); moments = parse.(Float64, tokens!(reader, nsf))
    nextline!(reader); nextline!(reader)
    energy_scale = parse(Float64, nextline!(reader)); energy_shift = parse(Float64, nextline!(reader))
    ntypes = parse(Int, nextline!(reader)); species = tokens!(reader, ntypes)
    references = parse.(Float64, tokens!(reader, ntypes))
    nextline!(reader); nextline!(reader); tokens!(reader, 3)
    return atomic_network_from_components(path, nlayers, maxnodes, nweights, nodes,
        activations, weights, description, atomtype, environment_names, minimum_radius,
        maximum_radius, descriptor_name, nsf, nparam, kinds, parameters, environments,
        minima, maxima, averages, moments, energy_scale, energy_shift, species, references)
end

function read_binary_atomic_network(path)
    reader = BinaryRecordReader(path)
    nlayers = integer!(reader); maxnodes = integer!(reader); nweights = integer!(reader); integer!(reader)
    nodes = integers!(reader, nlayers); activations = integers!(reader, nlayers-1)
    integers!(reader, nlayers); integers!(reader, nlayers); weights = reals!(reader, nweights)
    description = character!(reader); atomtype = character!(reader); nenv = integer!(reader)
    environment_names = characters!(reader, nenv); minimum_radius = real!(reader); maximum_radius = real!(reader)
    descriptor_name = character!(reader); nsf = integer!(reader); nparam = integer!(reader)
    kinds = integers!(reader, nsf); parameters = reals!(reader, nparam*nsf)
    environments = integers!(reader, 2nsf); integer!(reader)
    minima = reals!(reader, nsf); maxima = reals!(reader, nsf)
    averages = reals!(reader, nsf); moments = reals!(reader, nsf)
    character!(reader); record!(reader)
    energy_scale = real!(reader); energy_shift = real!(reader); ntypes = integer!(reader)
    species = characters!(reader, ntypes); references = reals!(reader, ntypes)
    integer!(reader); integer!(reader); reals!(reader, 3)
    return atomic_network_from_components(path, nlayers, maxnodes, nweights, nodes,
        activations, weights, description, atomtype, environment_names, minimum_radius,
        maximum_radius, descriptor_name, nsf, nparam, kinds, parameters, environments,
        minima, maxima, averages, moments, energy_scale, energy_shift, species, references)
end

function read_atomic_network(path::AbstractString)
    bytes = read(path)
    text = String(copy(bytes))
    isvalid(text) || return read_binary_atomic_network(path)
    return read_ascii_atomic_network(path)
end

format_floats(values) = join((@sprintf("%.17e", value) for value in values), " ")
format_ints(values) = join(values, " ")

function write_atomic_network(path, network::AtomicNetwork)
    nlayers = length(network.nodes); maxnodes = maximum(network.nodes)
    offsets = Int[]; offset = 0
    for layer in 1:nlayers-1
        push!(offsets, offset); offset += (network.nodes[layer] + 1) * network.nodes[layer+1]
    end
    push!(offsets, offset)
    species_index = Dict(symbol => index for (index, symbol) in enumerate(network.species))
    kinds, environments = Int[], Int[]; parameters = Float64[]
    for function_ in network.functions
        if function_.kind == 2
            push!(kinds, 2); append!(parameters, [function_.cutoff, function_.shift, function_.eta, 0.0])
            append!(environments, [species_index[function_.neighbor1], 0])
        else
            push!(kinds, function_.kind == 3 ? 4 : 5)
            append!(parameters, [function_.cutoff, function_.lambda, function_.zeta, function_.eta])
            append!(environments, [species_index[function_.neighbor1], species_index[function_.neighbor2]])
        end
    end
    moments = network.descriptor_shift.^2 .+ inv.(network.descriptor_scale).^2
    records = [string(nlayers), string(maxnodes), string(length(network.weights)), string(maxnodes*nlayers),
        format_ints(network.nodes), format_ints(network.activations), format_ints(offsets), format_ints(offsets),
        format_floats(network.weights), network.description, network.atomtype, string(length(network.species)),
        join(network.species, " "), @sprintf("%.17e", network.minimum_radius), @sprintf("%.17e", network.maximum_radius),
        "Behler2011", string(length(network.functions)), "4", format_ints(kinds), format_floats(parameters),
        format_ints(environments), "0", format_floats(network.descriptor_minimum),
        format_floats(network.descriptor_maximum), format_floats(network.descriptor_shift), format_floats(moments),
        "converted from n2p2", "T", @sprintf("%.17e", network.energy_scale), @sprintf("%.17e", network.energy_shift),
        string(length(network.species)), join(network.species, " "), format_floats(network.atomic_references),
        "0", "0", "0.0 0.0 0.0"]
    write(path, join(records, "\n") * "\n")
end
