macro time_and_print(expr)
    return quote
        local elapsed_time = @elapsed begin
            $(esc(expr))
        end
        @info "Time elapsed: $(elapsed_time) seconds per file."
    end
end

macro setup_logging(logfile)
    return quote
        # Open log file
        log_file = open($logfile, "w")

        # Create a Logger
        file_logger = SimpleLogger(log_file)
        tee_logger = TeeLogger(file_logger, global_logger())
        global_logger(tee_logger)

        # Make sure closing the logger at the end of the program.
        atexit(() -> close(log_file))
    end
end

function get_analysis_info(directory::String, filename::String)
    if !(endswith(directory,"/")) || !(endswith(directory,"\\"))
        if Sys.iswindows()
            directory *= "\\"
        elseif Sys.isapple() || Sys.islinux()
            directory *= "/"
        else
            error("SystemError: Unknown System Kernal.")
        end
    end
    log_info = Dict{String,Any}()
    log_info["Analysis Date"] = today()
    log_info["System Kernal"] = Sys.KERNEL
    log_info["File name"] = filename
    log_info["File full path"] = directory*filename
    log_info["Filesize (MB)"] = (filesize(directory*filename)/(1024*1024))
    return log_info
end

function First_logging()
    @info "Start Logging...\nPHANTOMJulia analysis Module\n  Version: 0.0.1\n    Make by Wei-Shan Su, 2024\n"
end

function initial_logging(info::Dict{String,Any})
    log_message = join(["$key: $value" for (key, value) in info], "\n")
    @info "----------------Information of analysis----------------\n" log_message
    @info "\nStart analysis..."
end
function last_logging()
    @info "End analysis."
end

for name in names(PHANTOMJulia; all=true)
    if name ∉ (:include, :eval) && isdefined(PHANTOMJulia, name)
        @eval export $name
    end
end

