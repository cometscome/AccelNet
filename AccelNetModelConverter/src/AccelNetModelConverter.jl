module AccelNetModelConverter

using Printf

export ConversionError, read_n2p2_model, read_atomic_network
export n2p2_to_accelnet, accelnet_to_n2p2, main

include("types.jl")
include("n2p2.jl")
include("accelnet.jl")
include("convert.jl")
include("cli.jl")

end
