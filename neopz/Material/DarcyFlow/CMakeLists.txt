# @file neopz/Material/DarcyFlow/CMakeLists.txt

set(public_headers
    TPZAnisotropicPermeability.h 
    TPZIsotropicPermeability.h
    TPZDarcyFlow.h
    TPZDarcyFlowAniso.h
    TPZHybridDarcyFlow.h
    TPZHybridDarcyFlowAniso.h
    TPZMixedDarcyFlowAniso.h
    TPZMixedDarcyFlow.h
    #TPZMixedDarcyFractureFlowAniso.h
    )

set(headers
    TPZAnisotropicPermeability.h 
    TPZIsotropicPermeability.h
    TPZDarcyFlow.h
    TPZDarcyFlowAniso.h
    TPZHybridDarcyFlow.h
    TPZHybridDarcyFlowAniso.h
    TPZMixedDarcyFlowAniso.h
    TPZMixedDarcyFlow.h
    #TPZMixedDarcyFractureFlowAniso.h
    )
set(sources
    TPZAnisotropicPermeability.cpp 
    TPZIsotropicPermeability.cpp
    TPZDarcyFlow.cpp
    TPZDarcyFlowAniso.cpp
    TPZHybridDarcyFlow.cpp
    TPZHybridDarcyFlowAniso.cpp
    TPZMixedDarcyFlowAniso.cpp
    TPZMixedDarcyFlow.cpp
    #TPZMixedDarcyFractureFlowAniso.cpp
    )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Material/DarcyFlow)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()

set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)
