cmake_minimum_required(VERSION 2.8)

project(MALT C)

set(standard "-std=c99 -D_XOPEN_SOURCE=600")

# Clang doesn't know -Wno-maybe-uninitialized, there aren't any
# reports of such anyway.
#set(warnings "-Wall -Wextra -Wwrite-strings -Wvla -Wno-maybe-uninitialized")
set(warnings "-Wall -Wextra -Wwrite-strings -Wvla")

set(CMAKE_C_FLAGS "${standard} ${warnings} ${CMAKE_C_FLAGS} ")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

add_executable(malt
               call_spice.c
               config.c
               corners.c
               define.c
               gplot.c
               list.c
               malt.c
               margins.c
               marg_opt_yield.c
               numerical.c
               optimize.c
               space.c
               stat_math.c
               toml.c)

target_link_libraries(malt m pthread)
