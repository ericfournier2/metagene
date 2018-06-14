parameter_manager <- R6Class("parameter_manager",
    public=list(
        initialize = function(param_values, param_validations=NULL,
                              overall_validation=NULL, locked=TRUE) {
            private$parameter_values = param_values
            if(!is.null(param_validations)) {
                private$parameter_validations = param_validations
            }
            
            if(!is.null(overall_validation)) {
                private$overall_validation = overall_validation
            }
            private$locked = locked
        },
        get = function(param_name) {
            private$validate_exists(param_name)
            return(private$parameter_values[[param_name]])
        },
        get_all = function() {
            return(private$parameter_values)
        },
        set = function(param_name, param_value) {
            private$validate_exists(param_name)
        
            # Don't do anything for NA params.
            if(!private$test_for_na(param_value)) {
                if(is.null(param_value)) {
                    # Cannot directly set NULL values inside a list outside of construction.
                    # See https://stackoverflow.com/questions/7944809/assigning-null-to-a-list-element-in-r
                    private$parameter_values[param_name] <- list(NULL)
                } else {
                    if(!is.null(private$parameter_validations[[param_name]])) {
                        private$parameter_validations[[param_name]](param_value)
                        #if(!private$parameter_validations[[param_name]](param_value)) {
                        #    warning("Parameter validation failed for ", param_name)
                        #    return(NULL)
                        #}
                    }
                    private$parameter_values[[param_name]] <- param_value
                }
            }
        },
        have_params_changed = function(...) {
            # This prologue makes it possible to infer parameter names from the
            # name of the variable it is passed in. This allows us to avoid
            # design=design, bin_count=bin_count repetitive code.
            #
            # It cannot be factorized into a function, since in any further call,
            # the argument list will deparse as "list(...)".
            param_names_alt = sapply( substitute(list(...)), deparse)[-1]
            arg_list = list(...)
            if(is.null(names(arg_list))) {
                names(arg_list) = param_names_alt
            } else {
                names(arg_list) = ifelse(names(arg_list)=="", param_names_alt, names(arg_list))
            }
            return(private$have_params_changed_internal(arg_list))
        },
        update_params = function(...) {
            # This prologue makes it possible to infer parameter names from the
            # name of the variable it is passed in. This allows us to avoid
            # design=design, bin_count=bin_count repetitive code.
            #
            # It cannot be factorized into a function, since in any further call,
            # the argument list will deparse as "list(...)".
            param_names_alt = sapply( substitute(list(...)), deparse)[-1]
            arg_list = list(...)
            if(is.null(names(arg_list))) {
                names(arg_list) = param_names_alt
            } else {
                names(arg_list) = ifelse(names(arg_list)=="", param_names_alt, names(arg_list))
            }

            private$update_params_internal(arg_list)
        }
    ),
    private=list(
        locked = TRUE,
        parameter_values=list(),
        parameter_validations=list(),
        overall_validation=NULL,
        have_params_changed_internal = function(param_list) {
            ret_val = FALSE
            for(i in names(param_list)) {
                param = param_list[[i]]
                # NA value means "keep what we had", so obviously that did not change.
                if(!private$test_for_na(param)) {
                    ret_val = ret_val || !identical(self$get(i), param)
                }
            }
        
            return(ret_val)        
        },
        update_params_internal = function(param_list) {
            if(private$have_params_changed_internal(param_list)) {
                for(i in names(param_list)) {
                    self$set(i, param_list[[i]])
                }
                return(TRUE)
            } else {
                return(FALSE)
            }
        },
        validate_exists = function(param_name) {
            if(private$locked && !(param_name %in% names(private$parameter_values))) {
                stop("Trying to access unknown parameter ", param_name)
            }
        },
        # Complex parameters such as data-frames will return
        # multiple values from is.na, which can't be used within if statement
        # without a warning. 
        # To test for a single NA value, we first eliminate the obvious problem
        # cases (NULL, lists) before testing for NA directly.
        test_for_na = function(value) {
            return(!(is.null(value) || is.list(value) || !is.na(value)))
        }
    )
)