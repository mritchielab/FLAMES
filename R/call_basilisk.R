callBasilisk <- function(env_name, FUN, ...) {
    proc <- basiliskStart(env_name)
    on.exit(basiliskStop(proc))

    return_value <- basiliskRun(proc, FUN, ...)

    return_value
}