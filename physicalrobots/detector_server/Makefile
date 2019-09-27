# Externals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Set up your externals here.

# Aruco
ARUCO_DIR   := ../aruco
ARUCO_INCL  := -I$(ARUCO_DIR)/src
ARUCO_LIBS  := -L$(ARUCO_DIR)/build -laruco

ARUCO := $(ARUCO_DIR)/build/libaruco.a
$(ARUCO):
	@echo Building aruco library...
	@cd $(ARUCO_DIR) && $(MAKE)


# pthread
THREAD_LIBS := -lpthread

# open cv
OPENCV_LIBS := -lopencv_calib3d -lopencv_core \
               -lopencv_features2d -lopencv_flann -lopencv_highgui \
	       -lopencv_imgproc -lopencv_ml -lopencv_objdetect \
	       -lopencv_photo -lopencv_stitching -lopencv_superres \
	       -lopencv_video -lopencv_videostab -lopencv_videoio

# calibration
CALIBRATION_DIR  := ../calibration
CALIBRATION_INCL := -I$(CALIBRATION_DIR)
CALIBRATION_LIBS := -L$(CALIBRATION_DIR) -lparasol_calibration

CALIBRATION := $(CALIBRATION_DIR)/libparasol_calibration.a
$(CALIBRATION):
	@echo Building calibration library...
	@cd $(CALIBRATION_DIR) && $(MAKE)

# Aggregated externals.
EXTERNAL_INCL := $(ARUCO_INCL) $(CALIBRATION_INCL)
EXTERNAL_LIBS := $(ARUCO_LIBS) $(OPENCV_LIBS) $(THREAD_LIBS) $(CALIBRATION_LIBS)
EXTERNAL_DEFS :=


# Directory layout ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Configure the directory tree.
# - Object files will be placed in BUILD_DIR.
# - All other directories will be assumed to contain source code.

BUILD_DIR := build
SRC_DIRS  := $(shell find . -type d ! \( -name '.' -o -name $(BUILD_DIR) \) )

$(BUILD_DIR):
	@mkdir $@

vpath %.h   $(SRC_DIRS)
vpath %.cpp $(SRC_DIRS)


# Dependency Tracking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Set up automatic dependency tracking.

DEP_FILE   := $(BUILD_DIR)/Dependencies
DEP_GLOB    = $(BUILD_DIR)/$*.d
TRACK_DEPS  = cat $(DEP_GLOB) >> $(DEP_FILE)

-include $(DEP_FILE)


# Compiler setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Configure compiler, flags, and options.

CXX   := g++
FLAGS := -O2 -g3 -std=c++11
WARN  := -Wall -Wextra -Werror -Wno-deprecated
DEFS  := $(EXTERNAL_DEFS)
INCL  := -I. $(EXTERNAL_INCL)
DEPS   = -MMD -MF $(DEP_GLOB)

COMPILE = $(CXX) $(FLAGS) $(WARN) $(DEFS) $(INCL) $(DEPS)


# Linker setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Configure linker and objects.

LINKER    := $(CXX)
LINK_OPTS :=
LIBS      := $(EXTERNAL_LIBS)

LINK := $(LINKER) $(LINK_OPTS)


# Source configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Scan the local directory tree to look for source files with suffix '.cpp' but
# no definition for the main function.

SRCS := $(shell grep -rLIz -m 1 --include='*.cpp' 'int\s*main' .)


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


# Executables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Define the executable names and the recipes for linking them.

.DEFAULT_GOAL = detector

detector: $(BUILD_DIR)/main_server.o $(OBJS) | $(ARUCO) $(CALIBRATION)
	@echo Linking $@...
	@$(CXX) $(FLAGS) $(WARN) $(DEFS) $^ $(LIBS) -o $@

client: $(BUILD_DIR)/main_client.o $(OBJS) | $(ARUCO) $(CALIBRATION)
	@echo Linking $@...
	@$(CXX) $(FLAGS) $(WARN) $(DEFS) $^ $(LIBS) -o $@


# Cleanup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Define clean-up recipes.

.PHONY: clean
clean:
	@rm -rf $(BUILD_DIR) detector client


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
