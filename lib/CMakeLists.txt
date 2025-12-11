# -------- sources --------
set(LIB_SOURCES       
    ${CMAKE_CURRENT_SOURCE_DIR}/lodepng/lodepng.cpp  
)

set(LIB_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/simd/simd.h
    ${CMAKE_CURRENT_SOURCE_DIR}/lodepng/lodepng.h
)

add_library(structty_utils STATIC
    ${LIB_SOURCES}
    ${LIB_HEADERS}
)

set_target_properties(structty_utils PROPERTIES POSITION_INDEPENDENT_CODE ON)

# include paths
target_include_directories(structty_utils
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/simd
        ${CMAKE_CURRENT_SOURCE_DIR}/lodepng
        ${CMAKE_CURRENT_SOURCE_DIR}/gemmi/include
)

# external deps
target_link_libraries(structty_utils
    PUBLIC
        gemmi::gemmi_cpp
)
