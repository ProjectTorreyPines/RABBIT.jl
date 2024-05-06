Base.@kwdef mutable struct RABBIToutputs 
    powe_data::Union{Vector{<:Real},Vector{Vector{<:Real}},Matrix{<:Real},Missing} = missing
    powi_data::Union{Vector{<:Real},Vector{Vector{<:Real}},Matrix{<:Real},Missing} = missing
    jnbcd_data::Union{Vector{<:Real}, Vector{Vector{<:Real}},Matrix{<:Real},Missing} = missing
    rho_data::Union{Vector{<:Real},Vector{Vector{<:Real}},Matrix{<:Real},Missing} = missing
    time_data::Union{Vector{<:Real},Vector{Vector{<:Real}},Matrix{<:Real},Missing} = missing
end

function read_outputs(path::String; filename::String ="run")
    struct_len = 4 
    result_path = abspath(joinpath(path, "$filename/beam1/rtfi_result_oav.bin"))
    
    open(string(result_path), "r") do f
        seekend(f)
        file_size = position(f)
        seekstart(f) 
    
        while !eof(f)
            ntime = Int(struct_unpack(read(f, struct_len)))
            nrho = Int(struct_unpack(read(f, struct_len)))
            nv = Int(struct_unpack(read(f, struct_len)))
    
            time_data = Float64[struct_unpack(read(f, struct_len)) for _ in 1:ntime]
            rho_data = Float64[struct_unpack(read(f, struct_len)) for _ in 1:nrho]
    
            function read_and_reshape(f, element_count; dims)
                data = Float64[struct_unpack(read(f, struct_len)) for _ in 1:element_count]
                return reshape(data, dims)
            end
    
            bdens_data = read_and_reshape(f, nrho * ntime, dims=(ntime, nrho))
            press_data = read_and_reshape(f, nrho * ntime, dims=(ntime, nrho))
            powe_data = read_and_reshape(f, nrho * ntime, dims=(ntime, nrho))            
            powi_data = read_and_reshape(f, nrho * ntime, dims=(ntime, nrho))

            jfi_data = read_and_reshape(f, nrho * ntime, dims=(ntime, nrho))
            jnbcd_data = read_and_reshape(f, nrho * ntime, dims=(ntime, nrho))
    
            next_line = read(f, struct_len)
            if isempty(next_line)
                break
            end
            dV_data = read_and_reshape(f, nrho * ntime, dims=(ntime, nrho))
    
            bdep_data = read_and_reshape(f, nrho * nv * ntime, dims=(ntime, nv, nrho))
            bdep_k1_data = read_and_reshape(f, nrho * nv * ntime, dims=(ntime, nv, nrho))
    
            next_line = read(f, struct_len)
            if isempty(next_line)
                break
            end
            pheatI_data = Float64[struct_unpack(read(f, struct_len)) for _ in 1:ntime]
    
            pheatE_data = Float64[struct_unpack(read(f, struct_len)) for _ in 1:ntime]
            pheat_data = Float64[struct_unpack(read(f, struct_len)) for _ in 1:ntime]
            pshine_data = Float64[struct_unpack(read(f, struct_len)) for _ in 1:ntime]
    
            next_line = read(f, struct_len)
            if isempty(next_line)
                break
            end
    
            prot_data = read_and_reshape(f, nrho * ntime, dims=(ntime,nrho))
            ploss_data = read_and_reshape(f, nrho * ntime, dims=(ntime,nrho))
            
            return powe_data, powi_data, jnbcd_data, rho_data, time_data
        
            break
        end
    end  
    
end