# Add subdirectories

add_library(Methods "")
target_sources(Methods PUBLIC
        InputTreatment.h
        InputTreatment.cpp
        MeshInit.h
        MeshInit.cpp
        Output.cpp
        Output.h
        Solver.cpp
        Solver.h
	ScriptTools.cpp
	ScriptTools.h
        )

target_link_libraries(Methods PUBLIC Tools)
target_include_directories(Methods PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(Simulations)

