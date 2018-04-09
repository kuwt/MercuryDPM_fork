add_executable(periodicChute PeriodicChute.cc PeriodicChuteDriver.cpp)
target_link_libraries(periodicChute  Chute)

add_executable(MaserParallel MaserParallel.cpp )
target_link_libraries(MaserParallel  Chute)