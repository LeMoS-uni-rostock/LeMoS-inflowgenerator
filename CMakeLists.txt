if (INSIGHT_BUILD_PYTHONBINDINGS) # required for refdata

project(inflowgen)
set(INSIGHT_INSTALL_COMPONENT inflowgen)

if (INSIGHT_BUILD_TOOLKIT AND INSIGHT_BUILD_OPENFOAM)

    set(inflowgen_SOURCES 
        inflowflatplatebl.cpp inflowflatplatebl.h
        inflowchannel.cpp inflowchannel.h
#        inflowstatisticstest.cpp
    )

    add_library(inflowgen SHARED ${inflowgen_SOURCES})
    target_link_libraries(inflowgen testcases)
    target_include_directories(inflowgen 
    PUBLIC ${CMAKE_CURRENT_BINARY_DIR} 
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    )
    add_PDL(inflowgen "${inflowgen_SOURCES}")

    install(TARGETS inflowgen LIBRARY DESTINATION lib COMPONENT ${INSIGHT_INSTALL_COMPONENT})

    install_shared_file(inflowgen_module_cfg inflow-generator.module modules.d)

    add_subdirectory(extensions/openfoam/inflowGeneratorBC)
    add_subdirectory(extensions/openfoam/checkInflowGeneratorStatistics)
    add_subdirectory(extensions/openfoam/randomInlet)
#     #add_subdirectory(extensions/openfoam/initInflowGenerator)
#     add_subdirectory(extensions/openfoam/insertTurbulentStructure)
#     add_subdirectory(extensions/openfoam/tuneInflowGenerator)
#    add_subdirectory(extensions/openfoam/superposeIsotropicTurbulence)
   
endif()

add_subdirectory(test)

endif()
