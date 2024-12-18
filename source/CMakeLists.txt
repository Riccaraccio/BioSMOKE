add_executable(BioSMOKEpp ${CMAKE_CURRENT_LIST_DIR}/BioSMOKE.cpp)

target_include_directories(BioSMOKEpp PUBLIC ${PROJECT_SOURCE_DIR} ${OPENSMOKE_INCLUDE_DIR} ${COMPULSORY_INCLUDE})

target_link_libraries(BioSMOKEpp
  ${COMPULSORY_LIBRARIES}
  Eigen3::Eigen
  Boost::date_time
  Boost::filesystem
  Boost::program_options
  Boost::system
  Boost::regex
  Boost::timer
  Boost::chrono
  # Linear Algebra
  ${LINEAR_ALGEBRA_LIB}
)

target_compile_definitions(
  BioSMOKEpp
  PUBLIC OPENSMOKE_USE_MKL=0
         OPENSMOKE_USE_OPENBLAS=0
         OPENSMOKE_USE_SUNDIALS=0
         OPENSMOKE_USE_UMFPACK=0
         OPENSMOKE_USE_BZZMATH=0
         OPENSMOKE_USE_SUPERLU_SERIAL=0
         OPENSMOKE_USE_RADAU=0
         OPENSMOKE_USE_ODEPACK=0
         OPENSMOKE_USE_DVODE=0
         OPENSMOKE_USE_DASPK=0
         OPENSMOKE_USE_MEBDF=0)
