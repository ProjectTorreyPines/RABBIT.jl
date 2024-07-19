module RABBIT

include("inputs.jl")

include("outputs.jl")

include("RABBIT_utils.jl")

const document = Dict()
document[Symbol(@__MODULE__)] = [name for name in Base.names(@__MODULE__, all=false, imported=false) if name != Symbol(@__MODULE__)]

end
