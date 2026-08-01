function usage(io=stdout)
    println(io, "Usage:")
    println(io, "  accelnet-model-converter.jl n2p2-to-accelnet INPUT_DIR OUTPUT_DIR")
    println(io, "  accelnet-model-converter.jl accelnet-to-n2p2 OUTPUT_DIR NETWORK...")
end

function main(arguments=ARGS)
    if isempty(arguments) || arguments[1] in ("-h", "--help")
        usage(); return 0
    end
    try
        if arguments[1] == "n2p2-to-accelnet" && length(arguments) == 3
            written = n2p2_to_accelnet(arguments[2], arguments[3])
        elseif arguments[1] == "accelnet-to-n2p2" && length(arguments) >= 3
            written = accelnet_to_n2p2(arguments[3:end], arguments[2])
        else
            usage(stderr); return 2
        end
        foreach(println, written)
        return 0
    catch error
        if error isa ConversionError || error isa SystemError || error isa ArgumentError
            print(stderr, "error: "); showerror(stderr, error); println(stderr); return 2
        end
        rethrow()
    end
end
