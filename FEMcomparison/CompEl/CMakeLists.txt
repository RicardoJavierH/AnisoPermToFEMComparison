target_sources(Tools PUBLIC
        LCC_MultiphysicsInterfaceEl.cpp
        LCC_MultiphysicsInterfaceEl.h
        )

target_include_directories(Tools PUBLIC ${CMAKE_CURRENT_LIST_DIR})