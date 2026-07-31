# Lime — a thin convenience wrapper around the CMake build.
#
# The two-command CMake flow in the README remains the primary route; this file
# only saves typing it. Installing is deliberately its own target rather than a
# side effect of building: `auval` reports on whatever is installed, so the copy
# must be a step you chose to run, never one the build ran for you.
#
#   make all        the whole flow in one word: build, run all six test
#                   suites, install, validate — stopping at the first failure
#   make            configure and build (Release, native architecture)
#   make universal  the same, built for both arm64 and x86_64
#   make test       build, then run all six DSP test suites through ctest
#   make install    copy the AU into ~/Library/Audio/Plug-Ins/Components and
#                   the VST3 into /Library/Audio/Plug-Ins/VST3 (sudo asks for
#                   your password for that one), then validate the AU with
#                   auval
#   make uninstall  remove both installed plug-ins
#   make clean      delete the build directory
#   make help       print this summary in the terminal
#
# BUILD_DIR overrides the build directory (default: build). ARCHS passes an
# architecture list straight through to CMake — `make universal` is shorthand
# for ARCHS="arm64;x86_64". CMake caches the architecture choice, so a later
# plain `make` or `make install` in the same build directory keeps it.
# AU_DEST and VST3_DEST override the install folders.

BUILD_DIR ?= build
ARCHS ?=

CMAKE_FLAGS := -DCMAKE_BUILD_TYPE=Release
ifneq ($(ARCHS),)
CMAKE_FLAGS += -DCMAKE_OSX_ARCHITECTURES="$(ARCHS)"
endif

ARTEFACTS := $(BUILD_DIR)/Lime_artefacts/Release

# The AU registrar scans the per-user folder, so the component installs
# without privileges. The VST3 goes to the system-wide folder instead: it is
# where hosts reliably look — a machine's existing VST3s live there — and it
# is root-owned, so that copy runs under sudo, which asks for the password in
# the terminal. Override either destination if your setup differs.
AU_DEST   ?= $(HOME)/Library/Audio/Plug-Ins/Components
VST3_DEST ?= /Library/Audio/Plug-Ins/VST3

USER_VST3 := $(HOME)/Library/Audio/Plug-Ins/VST3

.PHONY: all build universal test install uninstall clean help tools

# Plain `make` builds and nothing more — installing is never a side effect of
# a target that doesn't say so.
.DEFAULT_GOAL := build

# A fresh Mac ships with neither a compiler nor CMake. Naming the missing
# tool and the command that installs it beats the page of configure errors
# CMake prints when it finds no compiler.
tools:
	@xcode-select -p >/dev/null 2>&1 || { \
	    echo "No C++ compiler: the Xcode Command Line Tools are not installed."; \
	    echo "Run   xcode-select --install   accept the dialog, then re-run make."; \
	    exit 1; }
	@command -v cmake >/dev/null 2>&1 || { \
	    echo "CMake is not installed. Get it with Homebrew — brew install cmake —"; \
	    echo "or from https://cmake.org/download/, then re-run make."; \
	    exit 1; }

build: tools
	cmake -B "$(BUILD_DIR)" $(CMAKE_FLAGS)
	cmake --build "$(BUILD_DIR)" --parallel

universal:
	$(MAKE) ARCHS="arm64;x86_64" build

# Build, prove it, then install it: make stops at the first failing step, so a
# plug-in that fails its own test suite never reaches the Library folders.
all: test install

test: build
	ctest --test-dir "$(BUILD_DIR)" --output-on-failure

# The rm before each cp matters: cp -R into an existing bundle nests the new
# bundle inside the old one instead of replacing it. The VST3 steps escalate
# to sudo only when the destination refuses the plain command, so an
# overridden writable destination never sees a password prompt. A leftover
# per-user copy from an earlier release is removed too — two Lime.vst3s on
# the search path would shadow each other.
install: build
	mkdir -p "$(AU_DEST)"
	rm -rf "$(AU_DEST)/Lime.component"
	cp -R "$(ARTEFACTS)/AU/Lime.component" "$(AU_DEST)/"
	@if [ -w "$(VST3_DEST)" ]; then \
	    rm -rf "$(VST3_DEST)/Lime.vst3" && \
	    cp -R "$(ARTEFACTS)/VST3/Lime.vst3" "$(VST3_DEST)/"; \
	else \
	    echo "$(VST3_DEST) needs administrator rights — sudo will ask for your password."; \
	    sudo mkdir -p "$(VST3_DEST)" && \
	    sudo rm -rf "$(VST3_DEST)/Lime.vst3" && \
	    sudo cp -R "$(ARTEFACTS)/VST3/Lime.vst3" "$(VST3_DEST)/"; \
	fi
	@if [ "$(VST3_DEST)" != "$(USER_VST3)" ]; then rm -rf "$(USER_VST3)/Lime.vst3"; fi
	auval -v aufx Lim1 APTI
	@echo "Installed. Restart your host so it rescans its plug-in list."

uninstall:
	rm -rf "$(AU_DEST)/Lime.component"
	rm -rf "$(USER_VST3)/Lime.vst3"
	@if [ -w "$(VST3_DEST)" ] || [ ! -e "$(VST3_DEST)/Lime.vst3" ]; then \
	    rm -rf "$(VST3_DEST)/Lime.vst3"; \
	else \
	    echo "$(VST3_DEST) needs administrator rights — sudo will ask for your password."; \
	    sudo rm -rf "$(VST3_DEST)/Lime.vst3"; \
	fi

clean:
	rm -rf "$(BUILD_DIR)"

# The same summary the header comment carries, reachable without opening the
# file.
help:
	@echo "Lime — build and install from source. Targets:"
	@echo ""
	@echo "  make all        the whole flow in one word: build, run all six test"
	@echo "                  suites, install, validate — stopping at the first failure"
	@echo "  make            configure and build (Release, native architecture)"
	@echo "  make universal  the same, built for both arm64 and x86_64"
	@echo "  make test       build, then run all six DSP test suites through ctest"
	@echo "  make install    copy the AU into ~/Library/Audio/Plug-Ins/Components and"
	@echo "                  the VST3 into /Library/Audio/Plug-Ins/VST3 (sudo asks for"
	@echo "                  your password for that one), then validate the AU with auval"
	@echo "  make uninstall  remove both installed plug-ins"
	@echo "  make clean      delete the build directory"
	@echo "  make help       this summary"
	@echo ""
	@echo "Variables, overridable per invocation (make install VST3_DEST=...):"
	@echo ""
	@echo "  BUILD_DIR   build directory                 (default: build)"
	@echo "  ARCHS       CMake architecture list         (default: native)"
	@echo "  AU_DEST     Audio Unit install folder       (default: ~/Library/Audio/Plug-Ins/Components)"
	@echo "  VST3_DEST   VST3 install folder             (default: /Library/Audio/Plug-Ins/VST3)"
	@echo ""
	@echo "Prerequisites: the Xcode Command Line Tools (xcode-select --install) and"
	@echo "CMake (brew install cmake). Every build target checks for both and says"
	@echo "which is missing."
	@echo ""
	@echo "The Makefile is a thin wrapper over CMake; see README.md, section Building."
