#This file creates a library with tools
target_sources(Tools PUBLIC
    DataStructure.h
    Tools.h
    Tools.cpp
    computStatist.h
    computStatist.cpp
)

target_link_libraries(Tools PUBLIC NeoPZ::pz)
target_include_directories(Tools PUBLIC ${PZ_INCLUDE_DIRS})
target_include_directories(Tools PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
