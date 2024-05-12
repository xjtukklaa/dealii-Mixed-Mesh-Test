cmake_minimum_required(VERSION 3.27)
find_package(deal.II 9.5.0 REQUIRED HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})

set(TARGET Mixed_Test)
set(TARGET_SRC ${TARGET}.cpp)

deal_ii_initialize_cached_variables()
project(${TARGET})
deal_ii_invoke_autopilot()