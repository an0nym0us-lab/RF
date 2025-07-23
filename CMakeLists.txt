cmake_minimum_required(VERSION 3.10)

# set the project name
project(CostSketch VERSION 0.1)

# add sub dir as lib
add_subdirectory(util)

# 添加 telescoping-filter 目录中的源文件
set(TELESCOPING_FILTER_SOURCES
    ${CMAKE_SOURCE_DIR}/nonlearnedfilter/telescoping-filter/src/bit_util.c
    ${CMAKE_SOURCE_DIR}/nonlearnedfilter/telescoping-filter/src/murmur3.c
    ${CMAKE_SOURCE_DIR}/nonlearnedfilter/telescoping-filter/src/rsqf.c
    ${CMAKE_SOURCE_DIR}/nonlearnedfilter/telescoping-filter/src/set.c
    ${CMAKE_SOURCE_DIR}/nonlearnedfilter/telescoping-filter/src/taf.c
    ${CMAKE_SOURCE_DIR}/nonlearnedfilter/telescoping-filter/src/utaf.c
    ${CMAKE_SOURCE_DIR}/nonlearnedfilter/telescoping-filter/src/arcd.c
)
# 添加库
add_library(telescoping_filter ${TELESCOPING_FILTER_SOURCES})
# add the executable


add_executable(experiment_semi_adaptive experiment_semi_adaptive.cpp)

add_executable(experiment_full_adaptive experiment_full_adaptive.cpp)
# add_executable(experiment_ModuHashNum experiment_ModuHashNum.cpp)
# add_executable(experiment_InitHashNum experiment_InitHashNum.cpp)
# add_executable(experiment_ratioDelete experiment_ratioDelete.cpp)
# add_executable(experiment_delete experiment_delete.cpp)
# add_executable(experiment_aver experiment_aver.cpp)
# add_executable(experiment_latency experiment_latency.cpp)
# add_executable(experiment_oldIndex experiment_oldIndex.cpp)
# add_executable(experiment_dynamicHashSize experiment_dynamicHashSize.cpp)
# add_executable(experiment_counterVersion experiment_counterVersion.cpp)
# add_executable(experiment_latency_counterV_posneg experiment_latency_counterV_posneg.cpp)

# link executable to library
target_link_libraries(experiment_semi_adaptive PUBLIC util)
target_link_libraries(experiment_full_adaptive  util telescoping_filter)

# target_link_libraries(experiment_numOfDelete PUBLIC util)
# target_link_libraries(experiment_ModuHashNum PUBLIC util)
# target_link_libraries(experiment_InitHashNum PUBLIC util)
# target_link_libraries(experiment_ratioDelete PUBLIC util)
# target_link_libraries(experiment_delete PUBLIC util)
# target_link_libraries(experiment_aver PUBLIC util)
# target_link_libraries(experiment_latency PUBLIC util)
# target_link_libraries(experiment_oldIndex PUBLIC util)
# target_link_libraries(experiment_dynamicHashSize PUBLIC util)
# target_link_libraries(experiment_counterVersion PUBLIC util)
# target_link_libraries(experiment_latency_counterV_posneg PUBLIC util)

# add the binary tree to the search path for include files
# target_include_directories(experiment_numOfDelete PUBLIC
#                           "${PROJECT_BINARY_DIR}"
#                           "${PROJECT_SOURCE_DIR}/util"
#                           ) 
target_include_directories(experiment_semi_adaptive PUBLIC
                            "${PROJECT_BINARY_DIR}"
                            "${PROJECT_SOURCE_DIR}/util"
                            )                    
target_include_directories(experiment_full_adaptive PUBLIC
                            "${PROJECT_BINARY_DIR}"
                            "${PROJECT_SOURCE_DIR}/util"
                            )                       
# 包含头文件目录
# target_include_directories(experiment_counterVersion PUBLIC
#                           "${PROJECT_BINARY_DIR}"
#                           "${PROJECT_SOURCE_DIR}/util"
#                           )
# target_include_directories(experiment_numOfDelete_temp PUBLIC
#                           "${PROJECT_BINARY_DIR}"
#                           "${PROJECT_SOURCE_DIR}/util"
#                           )
# target_include_directories(experiment_numOfDelete_temp1 PUBLIC
#                           "${PROJECT_BINARY_DIR}"
#                           "${PROJECT_SOURCE_DIR}/util"
#                           )
# target_include_directories(experiment_numOfDelete_uoset PUBLIC
#                           "${PROJECT_BINARY_DIR}"
#                           "${PROJECT_SOURCE_DIR}/util"
#                           )
# target_include_directories(experiment_ModuHashNum PUBLIC
#                           "${PROJECT_BINARY_DIR}"
#                           "${PROJECT_SOURCE_DIR}/util"
#                           )
# target_include_directories(experiment_InitHashNum PUBLIC
#                           "${PROJECT_BINARY_DIR}"
#                           "${PROJECT_SOURCE_DIR}/util"
#                           )
# target_include_directories(experiment_ratioDelete PUBLIC
#                           "${PROJECT_BINARY_DIR}"
#                           "${PROJECT_SOURCE_DIR}/util"
#                           )

# target_include_directories(experiment_delete PUBLIC
#                           "${PROJECT_BINARY_DIR}"
#                           "${PROJECT_SOURCE_DIR}/util"
#                           )
# target_include_directories(experiment_aver PUBLIC
#                           "${PROJECT_BINARY_DIR}"
#                           "${PROJECT_SOURCE_DIR}/util"
#                           )
# target_include_directories(experiment_latency PUBLIC
#                           "${PROJECT_BINARY_DIR}"
#                           "${PROJECT_SOURCE_DIR}/util"
#                           )
# target_include_directories(experiment_oldIndex PUBLIC
#                           "${PROJECT_BINARY_DIR}"
#                           "${PROJECT_SOURCE_DIR}/util"
#                           )
# target_include_directories(experiment_dynamicHashSize PUBLIC
#                           "${PROJECT_BINARY_DIR}"
#                           "${PROJECT_SOURCE_DIR}/util"
#                           )

# target_include_directories(experiment_latency_counterV_posneg PUBLIC
#                           "${PROJECT_BINARY_DIR}"
#                           "${PROJECT_SOURCE_DIR}/util"
#                           )
# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-O3")