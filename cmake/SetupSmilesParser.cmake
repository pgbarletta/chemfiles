find_package(BISON)
find_package(FLEX)

if(MSVC)
add_definitions("/D YY_NO_UNISTD_H")
endif()

if(FLEX_EXECUTABLE)
  flex_target(SmilesL ${PROJECT_SOURCE_DIR}/scripts/smiles.ll
              ${CMAKE_CURRENT_BINARY_DIR}/include/lex.yysmiles.cpp
             COMPILE_FLAGS "-Pyysmiles_" )
  set(FLEX_OUTPUT_FILES ${FLEX_SmilesL_OUTPUTS})
else()
  configure_file(${PROJECT_SOURCE_DIR}/scripts/lex.yysmiles.cpp
                 ${CMAKE_CURRENT_BINARY_DIR}/include/lex.yysmiles.cpp COPYONLY)
  FILE(GLOB FLEX_OUTPUT_FILES "${CMAKE_CURRENT_BINARY_DIR}/include/lex.*.cpp")
endif()

if(BISON_EXECUTABLE)
  bison_target(SmilesY ${PROJECT_SOURCE_DIR}/scripts/smiles.yy
               ${CMAKE_CURRENT_BINARY_DIR}/include/smiles.tab.cpp
               COMPILE_FLAGS "-pyysmiles_" )
  set(BISON_OUTPUT_FILES ${BISON_SmilesY_OUTPUTS})
  if(FLEX_EXECUTABLE)
    ADD_FLEX_BISON_DEPENDENCY(SmilesL SmilesY)
  endif(FLEX_EXECUTABLE)
else()
  configure_file(${PROJECT_SOURCE_DIR}/scripts/smiles.tab.hpp
                 ${CMAKE_CURRENT_BINARY_DIR}/include/smiles.tab.hpp  COPYONLY)
  configure_file(${PROJECT_SOURCE_DIR}/scripts/smiles.tab.cpp
                 ${CMAKE_CURRENT_BINARY_DIR}/include/smiles.tab.cpp COPYONLY)
  FILE(GLOB BISON_OUTPUT_FILES "${CMAKE_CURRENT_BINARY_DIR}/include/*.tab.?pp")
endif()

add_library(lexer_parser_objs OBJECT ${FLEX_OUTPUT_FILES} ${BISON_OUTPUT_FILES})

target_include_directories(lexer_parser_objs PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
)

# We are dealing with generated C code here!
set_source_files_properties( ${FLEX_OUTPUT_FILES} ${BISON_OUTPUT_FILES} PROPERTIES COMPILE_FLAGS
  "-Wno-old-style-cast -Wno-sign-conversion -Wno-conversion -Wno-sign-compare -Wno-unused-function" )

target_include_directories(lexer_parser_objs SYSTEM BEFORE PRIVATE ${EXTERNAL_INCLUDES})
