Base.@kwdef mutable struct RABBIToutputs 
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

            @show pheatI_data
    
            pheatE_data = Float64[struct_unpack(read(f, struct_len)) for _ in 1:ntime]
            pheat_data = Float64[struct_unpack(read(f, struct_len)) for _ in 1:ntime]
            pshine_data = Float64[struct_unpack(read(f, struct_len)) for _ in 1:ntime]
    
            next_line = read(f, struct_len)
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
            pcx_data = zeros(Float64, ntime)
            pcx_data[1] = struct_unpack(next_line)
            
            for i in 2:ntime
                pcx_data[i] = struct_unpack(read(f, struct_len))
            end

            if file_size - position(f) > nrho * ntime * struct_len
                rabbit_version_strlen_bytes = read(f, struct_len)
                rabbit_version_strlen = 10
            
                rabbit_version_bytes = read(f, rabbit_version_strlen)
        
                dArea = zeros(Float64, nrho * ntime)
                for i in 1:length(dArea)
                    dArea[i] = struct_unpack(read(f, struct_len))
                end
                dArea_data = reshape(dArea, (ntime, nrho))
        
                torqdepo = zeros(Float64, nrho * nv * ntime)
                for i in 1:length(torqdepo)
                    torqdepo[i] = struct_unpack(read(f, struct_len))
                end
                torqdepo_data = reshape(torqdepo, (ntime, nv, nrho))
        
                torqjxb = zeros(Float64, nrho * nv * ntime)
                for i in 1:length(torqjxb)
                    torqjxb[i] = struct_unpack(read(f, struct_len))
                end
                torqjxb_data = reshape(torqjxb, (ntime, nv, nrho))
        
                torqe = zeros(Float64, nrho * ntime)
                for i in 1:length(torqe)
                    torqe[i] = struct_unpack(read(f, struct_len))
                end
                torque_data = reshape(torqe, (ntime, nrho))
        
                torqi = zeros(Float64, nrho * ntime)
                for i in 1:length(torqi)
                    torqi[i] = struct_unpack(read(f, struct_len))
                end
                torqi_data = reshape(torqi, (ntime, nrho))
            
                next_line = read(f, struct_len)
                if isempty(next_line)
                    break
                end
                        
                if file_size - position(f) >= nrho * ntime * struct_len
                    wfi_par = zeros(Float64, nrho * ntime)
                    for i in 1:length(wfi_par)
                        wfi_par[i] = struct_unpack(read(f, struct_len))
                    end
                    wfi_par_data = reshape(wfi_par, (ntime, nrho))
                    wfi_perp_data = 1.5f0 .* press_data - wfi_par_data
                end
            
                if position(f) < file_size
                    println("Warning: $(file_size - position(f)) bytes remain in file!")
                end
            
            return powe_data, powi_data, jnbcd_data, bdep_data, torqdepo_data, rho_data, time_data
        
            break
            end
        end
    end   
end