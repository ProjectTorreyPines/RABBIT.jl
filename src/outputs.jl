Base.@kwdef mutable struct RABBIToutputs 
    powe_data::Union{Vector{<:Real},Vector{Vector{<:Real}},Matrix{<:Real},Missing} = missing
    powi_data::Union{Vector{<:Real},Vector{Vector{<:Real}},Matrix{<:Real},Missing} = missing
    rho_data::Union{Vector{<:Real},Vector{Vector{<:Real}},Matrix{<:Real},Missing} = missing
    time_data::Union{Vector{<:Real},Vector{Vector{<:Real}},Matrix{<:Real},Missing} = missing
end

function read_outputs()
    struct_len = 4 
    result_path = abspath(joinpath(dirname(@__DIR__), "rtfi_result_oav.bin"))
    
    open(string(result_path), "r") do f
        seekend(f)
        file_size = position(f)
        seekstart(f) 
    
        while !eof(f)
            ntime = Int(struct_unpack(read(f, struct_len)))
            nrho = Int(struct_unpack(read(f, struct_len)))
            nv = Int(struct_unpack(read(f, struct_len)))
    
            time_data = Float32[struct_unpack(read(f, struct_len)) for _ in 1:ntime]
            rho_data = Float32[struct_unpack(read(f, struct_len)) for _ in 1:nrho]
    
            function read_and_reshape(f, element_count; dims)
                data = Float32[struct_unpack(read(f, struct_len)) for _ in 1:element_count]
                return reshape(data, dims)
            end
    
            bdens_data = read_and_reshape(f, nrho * ntime, dims=(ntime, nrho))
            press_data = read_and_reshape(f, nrho * ntime, dims=(ntime, nrho))
            powe_data = read_and_reshape(f, nrho * ntime, dims=(ntime, nrho))            
            powi_data = read_and_reshape(f, nrho * ntime, dims=(ntime, nrho))
            
            return powe_data, powi_data, rho_data, time_data
        
            break
        end
    end  
    
end

