find_program(cppcheck_exe NAMES cppcheck REQUIRED)
if(cppcheck_exe)
  set(cppcheck_opts --enable=all --inline-suppr --quiet
                    --suppressions-list=${PROJECT_SOURCE_DIR}/.cppcheck)
  set(CMAKE_CXX_CPPCHECK ${cppcheck_exe} --std=c++20 ${cppcheck_opts})
else()
  error(
    "cppcheck static analyzer cannot be find please download it from: https://cppcheck.sourceforge.io/"
  )
endif()
