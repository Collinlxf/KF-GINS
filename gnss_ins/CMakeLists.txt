PROJECT(GnssInsLcLib)
CMAKE_MINIMUM_REQUIRED(VERSION 3.0)
CMAKE_POLICY(SET CMP0015 NEW)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fPIC")

add_definitions(-DKALMANX=16)
add_definitions(-DCMAKE_POSITION_INDEPENDENT_CODE=ON)

# include_directories(
#     ${PROJECT_SOURCE_DIR}/include
#     ${PROJECT_SOURCE_DIR}/src
# )

add_subdirectory(src/PublicAlgo)
add_subdirectory(src/LooselyCoupled)
add_subdirectory(src/EgoMotion)
add_subdirectory(src/Calibrator)

file(GLOB SRC_FILES "src/*.cpp" "src/*.h")
add_library(GnssInsLcLib SHARED ${SRC_FILES})
target_link_libraries(GnssInsLcLib yaml-cpp PublicAlgoLib LooselyCoupledLib EgoMotionLib EgoMotionLib ${LOC_PUBLIC_LIB})
# target_link_libraries(CmnLib libGeographiccc ${YAML_CPP_LIBRARIES})

install(TARGETS GnssInsLcLib 
        RUNTIME DESTINATION ./localization/fusion_localization/bin
        ARCHIVE DESTINATION ./localization/fusion_localization/lib
        LIBRARY DESTINATION ./localization/fusion_localization/lib
        INCLUDES DESTINATION ./localization/fusion_localization/include
        )

install(FILES src/DeadReckoning.h DESTINATION ./localization/fusion_localization/include)