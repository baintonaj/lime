# ==============================================================================
#  Lime — CMake build
#
#  Targets:
#    lime_dsp        static, JUCE-free double-precision DSP library
#    Lime            JUCE audio plugin (AU / VST3 / Standalone)   [LIME_BUILD_PLUGIN]
#    <Name>Tests     six standalone console test suites            [LIME_BUILD_TESTS]
#
#  The tests and lime_dsp never touch JUCE: configuring with
#  -DLIME_BUILD_PLUGIN=OFF neither fetches nor requires JUCE at all.
#
#  JUCE is fetched from GitHub (tag 9.0.0) by default. To use a local copy
#  instead, pass -DLIME_JUCE_SOURCE_DIR=/path/to/JUCE.
#
#  CMAKE_OSX_ARCHITECTURES is deliberately not set here; pass e.g.
#  -DCMAKE_OSX_ARCHITECTURES="arm64;x86_64" on the command line for
#  universal release builds.
# ==============================================================================

cmake_minimum_required(VERSION 3.22)

set(CMAKE_OSX_DEPLOYMENT_TARGET "10.13" CACHE STRING "Minimum macOS version")

project(Lime VERSION 1.0.2 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(LIME_BUILD_PLUGIN "Build the JUCE plugin (fetches JUCE)" ON)
option(LIME_BUILD_TESTS  "Build the standalone DSP test suites" ON)

# ------------------------------------------------------------------------------
# lime_dsp — JUCE-free double-precision DSP library
# ------------------------------------------------------------------------------

add_library(lime_dsp STATIC
    Lime/Source/dsp/AtypeEngine.cpp
    Lime/Source/dsp/Bark.cpp
    Lime/Source/dsp/Calibration.cpp
    Lime/Source/dsp/Fft64.cpp
    Lime/Source/dsp/FilterDesign.cpp
    Lime/Source/dsp/FixedBand.cpp
    Lime/Source/dsp/ModulationControl.cpp
    Lime/Source/dsp/PhaseRotator.cpp
    Lime/Source/dsp/SlidingLayer.cpp
    Lime/Source/dsp/SrEngine.cpp
    Lime/Source/dsp/SrSideChain.cpp
    Lime/Source/dsp/TapeChannel.cpp)

target_include_directories(lime_dsp PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/Lime/Source")

if(APPLE)
    # Fft64 uses vDSP when LIME_HAS_ACCELERATE (defaults on via __APPLE__).
    target_link_libraries(lime_dsp PUBLIC "-framework Accelerate")
endif()

# ------------------------------------------------------------------------------
# Tests — plain console executables, one per suite, no JUCE dependency
# ------------------------------------------------------------------------------

if(LIME_BUILD_TESTS)
    enable_testing()

    set(lime_test_suites
        DspTests
        CalibrationTests
        TapeTests
        BarkTests
        AtypeTests
        PhaseTests)

    foreach(suite IN LISTS lime_test_suites)
        add_executable(${suite} Lime/Source/tests/${suite}.cpp)
        target_link_libraries(${suite} PRIVATE lime_dsp)
        add_test(NAME ${suite} COMMAND ${suite})
    endforeach()
endif()

# ------------------------------------------------------------------------------
# Plugin — JUCE is only acquired when this is ON
# ------------------------------------------------------------------------------

if(LIME_BUILD_PLUGIN)
    set(LIME_JUCE_SOURCE_DIR "" CACHE PATH "Use a local JUCE checkout instead of FetchContent")

    if(LIME_JUCE_SOURCE_DIR)
        add_subdirectory("${LIME_JUCE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}/JUCE" EXCLUDE_FROM_ALL)
    else()
        include(FetchContent)
        FetchContent_Declare(JUCE
            GIT_REPOSITORY https://github.com/juce-framework/JUCE.git
            GIT_TAG        9.0.0
            GIT_SHALLOW    TRUE)
        FetchContent_MakeAvailable(JUCE)
    endif()

    juce_add_plugin(Lime
        PRODUCT_NAME "Lime"
        DESCRIPTION "Spectral companding dynamics"
        COMPANY_NAME "Aptitude Audio"
        COMPANY_COPYRIGHT "Andy Bainton (Aptitude Audio)"
        COMPANY_WEBSITE "https://github.com/baintonaj/lime"
        BUNDLE_ID com.AptitudeAudio.Lime
        PLUGIN_MANUFACTURER_CODE APTI
        PLUGIN_CODE Lim1
        FORMATS AU VST3 Standalone
        AU_MAIN_TYPE kAudioUnitType_Effect
        VST3_CATEGORIES Fx Dynamics
        COPY_PLUGIN_AFTER_BUILD FALSE
        IS_SYNTH FALSE
        NEEDS_MIDI_INPUT FALSE
        NEEDS_MIDI_OUTPUT FALSE)

    juce_generate_juce_header(Lime)

    target_sources(Lime PRIVATE
        Lime/Source/PluginProcessor.cpp
        Lime/Source/PluginEditor.cpp)

    target_compile_definitions(Lime PRIVATE
        JUCE_STRICT_REFCOUNTEDPOINTER=1
        JUCE_VST3_CAN_REPLACE_VST2=0
        JUCE_WEB_BROWSER=0
        JUCE_USE_CURL=0)

    target_link_libraries(Lime
        PRIVATE
            lime_dsp
            juce::juce_audio_utils
            juce::juce_dsp
            juce::juce_gui_extra
        PUBLIC
            juce::juce_recommended_config_flags
            juce::juce_recommended_warning_flags)
endif()
