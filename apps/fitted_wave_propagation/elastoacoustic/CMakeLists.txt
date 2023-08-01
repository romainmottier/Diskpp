add_executable(elastoacoustic elastoacoustic.cpp ${fitted_waves_headers} ${fitted_waves_sources})
target_link_libraries(elastoacoustic PRIVATE ${LINK_LIBS} Spectra::Spectra)

file(COPY ${meshes} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes/)
