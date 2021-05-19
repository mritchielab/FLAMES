#' Internal utility function for simplifying
#' calls to basiliskRun using a given basilisk environment
#' @param env_name the name of the basilisk env (made through BasiliskEnvironment) to execute code within
#' @param FUN the function to execute from with the basilisk environment
#' @param ... extra parameters required by FUN
#' @return the result of `FUN`
#' @importFrom basilisk basiliskStart basiliskStop basiliskRun
callBasilisk <- function(env_name, FUN, ...) {
    proc <- basiliskStart(env_name)
    on.exit(basiliskStop(proc))

    return_value <- basiliskRun(proc, FUN, ...)

    return_value
}
