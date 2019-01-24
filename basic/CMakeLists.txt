cmake_minimum_required(VERSION 3.1)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(basic_HEADERS
    assertions.h
    attribute_adapter.h
    attribute_copier.h
    attribute_life_cycle.h
    attribute_manager.h
    attribute_serializer.h
    attribute_store.h
    attribute.h
    basic_common.h
    basic_types.h
    canvas.h
    color.h
    counted.h
    dlist.h
    file_utils.h
    generic_attributes_io.h
    line_stream.h
    logger.h
    pointer_iterator.h
    progress.h
    rat.h
    raw_attribute_store.h
    record_id.h
    smart_pointer.h
    stop_watch.h
    )

set(basic_SOURCES
    assertions.cpp
    attribute_adapter.cpp
    attribute_life_cycle.cpp
    attribute_manager.cpp
    attribute_serializer.cpp
    attribute_store.cpp
    basic_types.cpp
    counted.cpp
    file_utils.cpp
    logger.cpp
    progress.cpp
    rat.cpp
    raw_attribute_store.cpp
    stop_watch.cpp
    )


add_library(${PROJECT_NAME} SHARED ${basic_SOURCES} ${basic_HEADERS})
set_target_properties(${PROJECT_NAME} PROPERTIES
                      FOLDER "PolyFit")


target_compile_definitions(${PROJECT_NAME} PRIVATE BASIC_EXPORTS)


if (MSVC)
    target_compile_definitions(${PROJECT_NAME} PRIVATE
        _CRT_SECURE_NO_WARNINGS
        _CRT_SECURE_NO_DEPRECATE
        )
endif()
