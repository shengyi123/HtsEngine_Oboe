# For more information about using CMake with Android Studio, read the
# documentation: https://d.android.com/studio/projects/add-native-code.html

# Sets the minimum version of CMake required to build the native library.

cmake_minimum_required(VERSION 3.4.1)

# Creates and names a library, sets it as either STATIC
# or SHARED, and provides the relative paths to its source code.
# You can define multiple libraries, and CMake builds them for you.
# Gradle automatically packages shared libraries with your APK.

### INCLUDE OBOE LIBRARY ###

# Set the path to the Oboe library directory
set (OBOE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../..)


# Add the Oboe library as a subproject. Since Oboe is an out-of-tree source library we must also
# specify a binary directory
add_subdirectory(${OBOE_DIR} ./oboe-bin)

# Include the Oboe headers
include_directories(${OBOE_DIR}/include)

### END OBOE INCLUDE SECTION ###

# Debug utilities
set (DEBUG_UTILS_PATH "../debug-utils")
set (DEBUG_UTILS_SOURCES ${DEBUG_UTILS_PATH}/trace.cpp)
include_directories(${DEBUG_UTILS_PATH})

# App specific sources
set (APP_DIR src/main/cpp)
file (GLOB_RECURSE APP_SOURCES
        ${APP_DIR}/hts_Engine.c
        ${APP_DIR}/jni_bridge.cpp
        ${APP_DIR}/PlayAudioEngine.cpp
        ${APP_DIR}/SineGenerator.cpp
        ${APP_DIR}/HTS_audio.c
        ${APP_DIR}/HTS_engine1.c
        ${APP_DIR}/HTS_gstream.c
        ${APP_DIR}/HTS_label.c
        ${APP_DIR}/HTS_misc.c
        ${APP_DIR}/HTS_model.c
        ${APP_DIR}/HTS_orthogonal_expansion.c
        ${APP_DIR}/HTS_pstream.c
        ${APP_DIR}/HTS_sstream.c
        ${APP_DIR}/HTS_vocoder.c
        ${APP_DIR}/HTS_Xlabel.c
        ${APP_DIR}/oec.c
        )

# Build the libhello-oboe library


# App specific sources

add_library( # Sets the name of the library.
        native-lib

        # Sets the library as a shared library.
        SHARED

        # Provides a relative path to your source file(s).
        ${DEBUG_UTILS_SOURCES}
        ${APP_SOURCES})

# Searches for a specified prebuilt library and stores the path as a
# variable. Because CMake includes system libraries in the search path by
# default, you only need to specify the name of the public NDK library
# you want to add. CMake verifies that the library exists before
# completing its build.

find_library( # Sets the name of the path variable.
        log-lib

        # Specifies the name of the NDK library that
        # you want CMake to locate.
        log)

# Specifies libraries CMake should link to your target library. You
# can link multiple libraries, such as libraries you define in this
# build script, prebuilt third-party libraries, or system libraries.

target_link_libraries( # Specifies the target library.
        native-lib
        android
        # Links the target library to the log library
        # included in the NDK.
        ${log-lib}
        oboe )

# Enable optimization flags: if having problems with source level debugging,
# disable -Ofast ( and debug ), re-enable after done debugging.
target_compile_options(native-lib PRIVATE -Wall -Werror "$<$<CONFIG:RELEASE>:-Ofast>")