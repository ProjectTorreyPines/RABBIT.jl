Base.@kwdef mutable struct RABBIToutput 
    powe_data::Union{Array{<:Real},Missing} = missing
    powi_data::Union{Array{<:Real},Missing} = missing
    jnbcd_data::Union{Array{<:Real},Missing} = missing
    bdep_data::Union{Array{<:Real},Missing} = missing
    torqdepo_data::Union{Array{<:Real},Missing} = missing
    rho_data::Union{Vector{<:Real},Missing} = missing
    time_data::Union{Vector{<:Real},Missing} = missing 
    nrate_data::Union{Array{<:Real}, Missing} = missing

end

function read_outputs(path::String; filename::String ="run")
    struct_len = 4 
    result_path = abspath(joinpath(path, "$filename/beam1/rtfi_result_oav.bin"))
    nrate_path = abspath(joinpath(path, "$filename/rtfi_nrate_oav.bin"))

    output = RABBIToutput()

    function read_and_reshape(f, element_count; dims)
        data = Float64[struct_unpack(read(f, struct_len)) for _ in 1:element_count]
        return reshape(data, dims)
    end
    
    open(string(result_path), "r") do f
        seekend(f)
        file_size = position(f)
        seekstart(f) 
    
        while !eof(f)
            ntime = Int(struct_unpack(read(f, struct_len)))
            nrho = Int(struct_unpack(read(f, struct_len)))
            nv = Int(struct_unpack(read(f, struct_len)))
    
            output.time_data = Float64[struct_unpack(read(f, struct_len)) for _ in 1:ntime]
            output.rho_data = Float64[struct_unpack(read(f, struct_len)) for _ in 1:nrho]
    
            bdens_data = read_and_reshape(f, nrho * ntime, dims=(nrho, ntime))
            press_data = read_and_reshape(f, nrho * ntime, dims=(nrho, ntime))
            output.powe_data = read_and_reshape(f, nrho * ntime, dims=(nrho, ntime))            
            output.powi_data = read_and_reshape(f, nrho * ntime, dims=(nrho, ntime))

            jfi_data = read_and_reshape(f, nrho * ntime, dims=(nrho, ntime))
            output.jnbcd_data = read_and_reshape(f, nrho * ntime, dims=(nrho, ntime))
    
            dV_data = read_and_reshape(f, ntime * nrho, dims = (nrho, ntime))
            output.bdep_data = read_and_reshape(f, nrho * nv * ntime, dims=(nrho, ntime, nv))
            bdep_k1_data = read_and_reshape(f, nrho * nv * ntime, dims=(ntime, nv, nrho))

            pheatI_data = read_and_reshape(f, ntime, dims=(ntime))
            pheatE_data = read_and_reshape(f, ntime, dims=(ntime))
            pheat_data = read_and_reshape(f, ntime, dims=(ntime))

            pshine_data = read_and_reshape(f, ntime, dims=(ntime))
            prot_data = read_and_reshape(f, ntime, dims=(ntime))
            ploss_data = read_and_reshape(f, ntime, dims=(ntime))
            pcx_data = read_and_reshape(f, ntime, dims=(ntime))

            if file_size - position(f) > nrho * ntime * struct_len
                rabbit_version_strlen = reinterpret(Int32, read(f, sizeof(Int32)))
                rabbit_version_bytes = read(f, rabbit_version_strlen[1])

                dArea = read_and_reshape(f, nrho * ntime, dims=(nrho, ntime))
                output.torqdepo_data = read_and_reshape(f, nrho * nv * ntime, dims=(nrho, ntime, nv))
                torqjxb_data = read_and_reshape(f, nrho * nv * ntime, dims=(nrho, ntime, nv))
                    
            break
            end
        end
    end 
    
    open(string(nrate_path), "r") do f
        seekstart(f)

        nrho = length(output.rho_data)
        ntime = length(output.time_data)
        while !eof(f)
            output.nrate_data = read_and_reshape(f, nrho * ntime, dims = (nrho, ntime))
        end
    end

    return output

end