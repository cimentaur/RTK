WRAP_GGO(rtkpolyquant_GGO_C rtkpolyquant.ggo ../rtkinputprojections_section.ggo ../rtk3Doutputimage_section.ggo ${RTK_BINARY_DIR}/rtkVersion.ggo)
add_executable(rtkpolyquant rtkpolyquant.cxx alg_polyquant.cxx ${rtkpolyquant_GGO_C})
target_link_libraries(rtkpolyquant RTK)

if (RTK_USE_CUDA)
  target_link_libraries(rtkpolyquant rtkcuda)
endif ()

# Installation code
if(NOT RTK_INSTALL_NO_EXECUTABLES)
  foreach(EXE_NAME rtkpolyquant) 
    install(TARGETS ${EXE_NAME}
      RUNTIME DESTINATION ${RTK_INSTALL_RUNTIME_DIR} COMPONENT Runtime
      LIBRARY DESTINATION ${RTK_INSTALL_LIB_DIR} COMPONENT RuntimeLibraries
      ARCHIVE DESTINATION ${RTK_INSTALL_ARCHIVE_DIR} COMPONENT Development)
  endforeach() 
endif()

