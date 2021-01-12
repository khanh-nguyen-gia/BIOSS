using Logging

"""
    formatMetaProg(title,code)

    Print the code produced in metaprogramming formated,
    with auto comments removed

"""
function formatMetaProg(title,code)

    tempStr = string(code)
    tempStr = split(tempStr,"\n")

    out = title * "\n \n "

    for i in 1:length(tempStr)
        if !occursin("#",tempStr[i])
            out *=replace(tempStr[i], "\n" => "")*"\n "
        end
    end

    return out
end



function initDebug()
    logger = SimpleLogger(stdout, Logging.Debug)
    LogLevel(-1000)
    old_logger = global_logger(logger)
end
