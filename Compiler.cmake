# public settings

if ( WIN32 )
    
    if ( NOT CYGWIN )
        if ( MSVC )
            
            if ( MSVC80 )
                
                message(" The current complier is Visual Studio 8 ...")
            
                SET(EXECUTABLE_OUTPUT_PATH ${GADGETRON_SOURCE_DIR}/prod/bin/vc8 CACHE STRING "Where to put the executables")
                SET(LIBRARY_OUTPUT_PATH ${GADGETRON_SOURCE_DIR}/prod/bin/vc8 CACHE STRING "Where to put the libraries")

                if ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )
                    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
                else ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )
                    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
                endif ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )
                
                if ( HAS_OPENMP )
                    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
                endif ( HAS_OPENMP )

            else ( MSVC80 )
    
                if ( MSVC_VERSION GREATER 1500 ) # vc10
                
                    message(" The current complier is Visual Studio 10 ...")
                    SET(EXECUTABLE_OUTPUT_PATH ${GADGETRON_SOURCE_DIR}/prod/bin/vc10 CACHE STRING "Where to put the executables")
                    SET(LIBRARY_OUTPUT_PATH ${GADGETRON_SOURCE_DIR}/prod/bin/vc10 CACHE STRING "Where to put the libraries")
                                    
                    if ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )
                        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
                    else ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )
                        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
                    endif ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )

                    if ( HAS_OPENMP )
                        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
                    endif ( HAS_OPENMP )
                
                else ( MSVC_VERSION GREATER 1500 )
                
                    if ( MSVC_VERSION GREATER 1400 ) # vc9
                    
                        message(" The current complier is Visual Studio 9 ...")
                        
                        SET(EXECUTABLE_OUTPUT_PATH ${GADGETRON_SOURCE_DIR}/prod/bin/vc9 CACHE STRING "Where to put the executables")
                        SET(LIBRARY_OUTPUT_PATH ${GADGETRON_SOURCE_DIR}/prod/bin/vc9 CACHE STRING "Where to put the libraries")
                                        
                        if ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )
                            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
                        else ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )
                            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
                        endif ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )

                        if ( HAS_OPENMP )
                            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
                        endif ( HAS_OPENMP )
                        
                    else ( MSVC_VERSION GREATER 1400 )
                        
                        if ( MSVC60 )
                        
                            message(" The current complier is Visual Studio 6 ...")
                        
                            SET(EXECUTABLE_OUTPUT_PATH ${GADGETRON_SOURCE_DIR}/prod/bin/vc6 CACHE STRING "Where to put the executables")
                            SET(LIBRARY_OUTPUT_PATH ${GADGETRON_SOURCE_DIR}/prod/bin/vc6 CACHE STRING "Where to put the libraries")

                            if ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )
                                SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
                            else ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )
                                SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
                            endif ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )
                            
                            if ( HAS_OPENMP )
                                SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
                            endif ( HAS_OPENMP )

                        else ( MSVC60 )
            
                            if ( MSVC70 ) # visual studio 2003.net
                            
                                message(" The current complier is Visual Studio 7 .NET ...")
                            
                                SET(EXECUTABLE_OUTPUT_PATH ${GADGETRON_SOURCE_DIR}/prod/bin/vc2003 CACHE STRING "Where to put the executables")
                                SET(LIBRARY_OUTPUT_PATH ${GADGETRON_SOURCE_DIR}/prod/bin/vc2003 CACHE STRING "Where to put the libraries")    

                                if ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )
                                    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
                                else ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )
                                    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W0")
                                endif ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )
                                    
                                if ( HAS_OPENMP )
                                    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
                                endif ( HAS_OPENMP )
                                
                            else ( MSVC70 )
            
                                message(" The current complier is Visual Studio ...")
            
                                SET(EXECUTABLE_OUTPUT_PATH ${GADGETRON_SOURCE_DIR}/prod/bin/vc CACHE STRING "Where to put the executables")
                                SET(LIBRARY_OUTPUT_PATH ${GADGETRON_SOURCE_DIR}/prod/bin/vc CACHE STRING "Where to put the libraries")
                                SET(CMAKE_INSTALL_PREFIX ${GADGETRON_SOURCE_DIR}/prod/Install/vc)
        
                            endif ( MSVC70 )
                        endif ( MSVC60 )
                    endif ( MSVC_VERSION GREATER 1400 )
                endif ( MSVC_VERSION GREATER 1500 )
            endif ( MSVC80 )
        endif ( MSVC )

    else ( NOT CYGWIN )
    
        message(" The current complier is Cygwin ...")
    
        SET(EXECUTABLE_OUTPUT_PATH ${GADGETRON_SOURCE_DIR}/prod/bin/cygwin CACHE STRING "Where to put the executables")
        SET(LIBRARY_OUTPUT_PATH ${GADGETRON_SOURCE_DIR}/prod/bin/cygwin CACHE STRING "Where to put the libraries")

        if ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS ) 
          SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wextra -Wall -Wno-deprecated -Wno-uninitialized -Wundef -Wpointer-arith -Wwrite-strings -Wconversion -Wsign-compare -Wmissing-noreturn -Wmissing-format-attribute -Wredundant-decls")
        else ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )
          SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -w")
        endif ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )

        if ( HAS_OPENMP )
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        endif ( HAS_OPENMP )
    
    endif ( NOT CYGWIN )
    
else ( WIN32 )

    if ( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCC )

        message(" The current complier is GCC ...")

        ADD_DEFINITIONS(-DBUILD_PLATFORM_LINUX)

        if ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )
          SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wextra -Wall -Wno-deprecated -Wno-uninitialized -Wundef -Wpointer-arith -Wwrite-strings -Wconversion -Wsign-compare -Wmissing-noreturn -Wmissing-format-attribute -Wredundant-decls")
        else ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )
          SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -w")
        endif ( NOT FLAG_FTK_SUPPRESS_COMPILING_WARNINGS )

        if ( HAS_OPENMP )
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        endif ( HAS_OPENMP )    

    else ( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCC )

        message(" The current complier is unknown ...")

        SET(EXECUTABLE_OUTPUT_PATH ${GADGETRON_SOURCE_DIR}/prod/bin/others CACHE STRING "Where to put the executables")
        SET(LIBRARY_OUTPUT_PATH ${GADGETRON_SOURCE_DIR}/prod/bin/others CACHE STRING "Where to put the libraries")

    endif ( CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCC )

endif ( WIN32 )

