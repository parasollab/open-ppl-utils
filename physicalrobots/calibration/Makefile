# Externals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Set up your externals here.

# OpenCV
OPENCV_LIBS := -lopencv_calib3d -lopencv_core \
               -lopencv_features2d -lopencv_flann -lopencv_highgui \
	       -lopencv_imgproc -lopencv_ml -lopencv_objdetect \
	       -lopencv_photo -lopencv_stitching -lopencv_superres \
	       -lopencv_video -lopencv_videostab -lopencv_videoio

# Aggregated externals.
EXTERNAL_INCL :=
EXTERNAL_DEFS :=
EXTERNAL_LIBS := $(OPENCV_LIBS)


# Compiler setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Configure compiler, flags, and options.

CXX   := g++
FLAGS := -O3 -g3 -std=c++11
WARN  := -Wall -Werror -Wno-error=unused-variable
DEFS  := $(EXTERNAL_DEFS)
INCL  := -I. $(EXTERNAL_INCL)

COMPILE := $(CXX) $(FLAGS) $(WARN) $(DEFS) $(INCL)


# Linker setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Configure linker and objects.

LIBS := $(EXTERNAL_LIBS)

LINK := $(CXX) $(WARN)


# Directory layout ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Configure the directory tree.
# - Object files will be placed in BUILD_DIR.
# - All other directories will be assumed to contain source code.

BUILD_DIR := build
SRC_DIRS  := $(shell find . -type d ! \( -name '.' -o -name $(BUILD_DIR) \) )

$(BUILD_DIR):
	@mkdir $@

vpath %.h $(SRC_DIRS)
vpath %.cpp $(SRC_DIRS)


# Source configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Scan the local directory tree to look for source files with suffix '.cpp'.

SRCS := $(shell find . -type f \( -name '*.cpp' ! -name 'main.cpp' \) \
                       -printf '%f ')


# Object configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Create an object name for each source file and define object compilation
# recipes.

OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o, $(SRCS))

$(BUILD_DIR)/%.o: %.cpp %.h | $(BUILD_DIR)
	@echo Compiling $<...
	@$(COMPILE) -c $< -o $@

$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	@echo Compiling $< \(no associated header found\)...
	@$(COMPILE) -c $< -o $@


# Library ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Define the executable name and the recipe for linking it.

LIBRARY := libparasol_calibration.a

$(LIBRARY): $(OBJS)
	@echo Linking static library...
	@ar rcs $@ $^


# Executable ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Define the executable name and the recipe for linking it.

EXECUTABLE := calibrate

.DEFAULT_GOAL := $(EXECUTABLE)

$(EXECUTABLE): $(BUILD_DIR)/main.o $(LIBRARY)
	@echo Linking $@...
	@$(LINK) $^ $(LIBS) -o $@


# Cleanup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Define clean-up recipes.

.PHONY: clean
clean:
	@rm -rf $(LIBRARY) $(EXECUTABLE) $(BUILD_DIR)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
