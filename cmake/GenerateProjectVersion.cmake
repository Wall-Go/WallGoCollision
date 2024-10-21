# Sets INOUT_VERSION in parent scope
function(GenerateProjectVersion INOUT_VERSION)

    if (DEFINED SKBUILD_PROJECT_VERSION)
        # Get version number from scikit-build
        set(LOCAL_VERSION ${SKBUILD_PROJECT_VERSION})
    else()
    
        # Attempt to read git tag
        find_package(Git QUIET)
        if (GIT_FOUND)
            # Will fail if ran outside a git repo, GIT_DESCRIBE_RESULT for success check
            execute_process(
                # get most recent tag
                COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0
                OUTPUT_VARIABLE GIT_TAG
                RESULT_VARIABLE GIT_DESCRIBE_RESULT
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
            )

            if (GIT_DESCRIBE_RESULT EQUAL 0)
                # tag is vX.Y.Z, strip the v
                string(REGEX REPLACE "^v" "" GIT_TAG ${GIT_TAG})
            endif()

            # Strip all extra from semantic version. So leaves only X.Y.Z
            string(REGEX REPLACE "([0-9]+\\.[0-9]+\\.[0-9]+).*" "\\1" LOCAL_VERSION ${GIT_TAG})
         
        endif()
    endif()

    if (NOT DEFINED LOCAL_VERSION)
        # Fallback in case everything else failed
        set(LOCAL_VERSION "0.0.0")
        message(WARNING "Failed to determine project version from Git tags. Defaulting to \"${LOCAL_VERSION}\".")
    endif()

    set(${INOUT_VERSION} ${LOCAL_VERSION} PARENT_SCOPE)
endfunction()