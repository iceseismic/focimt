cmake_minimum_required(VERSION 2.8)

file(GLOB foci_SRC
    "*.h"
    "*.cpp"
    )

add_executable(focimt ${foci_SRC})
target_link_libraries(focimt cairo)
