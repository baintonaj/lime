# Lime — a thin convenience wrapper around the CMake build.
#
# The two-command CMake flow in the README remains the primary route; this file
# only saves typing it. Installing is deliberately its own target rather than a
# side effect of building: `auval` reports on whatever is installed, so the copy
# must be a step you chose to run, never one the build ran for you.
#
#   make            configure and build (Release, native architecture)
#   make universal  the same, built for both arm64 and x86_64
#   make test       build, then run all six DSP test suites through ctest
#   make install    copy the AU and VST3 into ~/Library/Audio/Plug-Ins,
#                   then validate the AU with auval
#   make uninstall  remove both installed plug-ins
#   make clean      delete the build directory
#
# BUILD_DIR overrides the build directory (default: build). ARCHS passes an
# architecture list straight through to CMake — `make universal` is shorthand
# for ARCHS="arm64;x86_64". CMake caches the architecture choice, so a later
# plain `make` or `make install` in the same build directory keeps it.

BUILD_DIR ?= build
ARCHS ?=

CMAKE_FLAGS := -DCMAKE_BUILD_TYPE=Release
ifneq ($(ARCHS),)
CMAKE_FLAGS += -DCMAKE_OSX_ARCHITECTURES="$(ARCHS)"
endif

ARTEFACTS := $(BUILD_DIR)/Lime_artefacts/Release
AU_DEST   := $(HOME)/Library/Audio/Plug-Ins/Components
VST3_DEST := $(HOME)/Library/Audio/Plug-Ins/VST3

.PHONY: all universal test install uninstall clean

all:
	cmake -B "$(BUILD_DIR)" $(CMAKE_FLAGS)
	cmake --build "$(BUILD_DIR)" --parallel

universal:
	$(MAKE) ARCHS="arm64;x86_64" all

test: all
	ctest --test-dir "$(BUILD_DIR)" --output-on-failure

# The rm before each cp matters: cp -R into an existing bundle nests the new
# bundle inside the old one instead of replacing it.
install: all
	mkdir -p "$(AU_DEST)" "$(VST3_DEST)"
	rm -rf "$(AU_DEST)/Lime.component"
	cp -R "$(ARTEFACTS)/AU/Lime.component" "$(AU_DEST)/"
	rm -rf "$(VST3_DEST)/Lime.vst3"
	cp -R "$(ARTEFACTS)/VST3/Lime.vst3" "$(VST3_DEST)/"
	auval -v aufx Lim1 APTI
	@echo "Installed. Restart your host so it rescans its plug-in list."

uninstall:
	rm -rf "$(AU_DEST)/Lime.component"
	rm -rf "$(VST3_DEST)/Lime.vst3"

clean:
	rm -rf "$(BUILD_DIR)"
