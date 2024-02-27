target_sources(Tools PUBLIC
        TPZCreateHybridizedMixedSpaces.cpp
        TPZCreateHybridizedMixedSpaces.h
        TPZCreateHybridizedH1Spaces.cpp
        TPZCreateHybridizedH1Spaces.h
        )

target_include_directories(Tools PUBLIC ${CMAKE_CURRENT_LIST_DIR})