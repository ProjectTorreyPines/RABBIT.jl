Base.@kwdef mutable struct RABBIToutput 
    powe_data::Union{Matrix{<:Real},Missing} = missing
    powi_data::Union{Matrix{<:Real},Missing} = missing
    jnbcd_data::Union{Matrix{<:Real},Missing} = missing
    bdep_data::Union{Matrix{<:Real},Array{<:Real},Missing} = missing
    torque_data::Union{Matrix{<:Real},Array{<:Real},Missing} = missing
    rho_data::Union{Vector{<:Real},Missing} = missing
    time_data::Union{Vector{<:Real},Missing} = missing 

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

            dV_data = zeros(Float64, (ntime, nrho))
            dV_data[1] = struct_unpack(next_line)
            
            for i in 2:(ntime * nrho)
                dV_data[i] = struct_unpack(read(f, struct_len))
            end
            bdep_data = read_and_reshape(f, nrho * nv * ntime, dims=(ntime, nv, nrho))
            bdep_k1_data = read_and_reshape(f, nrho * nv * ntime, dims=(ntime, nv, nrho))

            next_line = read(f, struct_len)
            if isempty(next_line)
                break
            end
            pheatI_data = zeros(Float64, ntime)
            for i in 2:ntime
                pheatI_data[i] = struct_unpack(read(f, struct_len))
            end

            pheatE_data = Float64[struct_unpack(read(f, struct_len)) for _ in 1:ntime]
            pheat_data = Float64[struct_unpack(read(f, struct_len)) for _ in 1:ntime]
            pshine_data = Float64[struct_unpack(read(f, struct_len)) for _ in 1:ntime]
    
            # next_line = read(f, struct_len)
            if isempty(next_line)
                break
            end
    
            prot_data = zeros(Float64, ntime)
            prot_data[1] = struct_unpack(read(f, struct_len))
            for i in 2:ntime
                prot_data[i] = struct_unpack(read(f, struct_len))[1]
            end
            ploss_data = read_and_reshape(f, ntime, dims=(ntime))

            next_line = read(f, struct_len)
            if isempty(next_line)
                break
            end
            pcx_data = ones(Float64, ntime)
            pcx_data[1] = struct_unpack(next_line)
            
            for i in 2:ntime
                pcx_data[i] = struct_unpack(read(f, struct_len))
            end

            if file_size - position(f) > nrho * ntime * struct_len
                # rabbit_version_strlen_bytes = read(f, Int32)
                rabbit_version_strlen = reinterpret(Int32, read(f, sizeof(Int32)))
                @show rabbit_version_strlen
                rabbit_version_strlen = 10
                if rabbit_version_strlen > 0
            
                # rabbit_version_bytes = read(f, rabbit_version_strlen)
                # rabbit_version = reinterpret(UInt8, read(f, rabbit_version_strlen))[1]
                    rabbit_version_bytes = read(f, rabbit_version_strlen)
                    rabbit_version = String(rabbit_version_bytes)
                    @show rabbit_version_bytes
                    @show rabbit_version

                    dArea = zeros(Float64, nrho * ntime)
                    for i in 1:length(dArea)
                        dArea[i] = struct_unpack(read(f, struct_len))
                    end
                    @show dArea

                    torqdepo_data = read_and_reshape(f, nrho * nv * ntime, dims=(ntime, nv, nrho))
                    torqjxb_data = read_and_reshape(f, nrho * nv * ntime, dims=(ntime, nv, nrho))

                    @show torqjxb_data
                end
            
            return powe_data, powi_data, jnbcd_data, bdep_data, torqdepo_data, rho_data, time_data
        
            break
            end
        end
    end   
end