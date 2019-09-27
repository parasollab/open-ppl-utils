# Compiler setup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

CXX   := g++
FLAGS := -O2 -g3 -std=c++11
WARN  := -Wno-deprecated-declarations -Wno-ignored-attributes -Wno-deprecated

DEFS  :=

EIGEN3_INCL := -I./3rdparty/eigen3   # Include the local eigen3 library.
INCL  := $(EIGEN3_INCL)

COMPILE = $(CXX) $(FLAGS) $(WARN) $(DEFS) $(INCL)

# Directory layout ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Final outputs will land in BUILD_DIR.
# Intermediate objects will land in a corresponding subdirectory, e.g.,
# 	BUILD_DIR/LIB_DIR for library intermediates.
# We aren't using anything in the tests or utils_gl directories, so they have
# 	been removed from the codebase.

BUILD_DIR       := build

LIB_DIR         := src
UTILS_DIR       := utils
CALIBRATION_DIR := utils_calibration
MARKERMAP_DIR   := utils_markermap

SRC_DIRS        := $(shell find $(LIB_DIR) -type d) \
		   $(shell find $(UTILS_DIR) -type d) \
		   $(shell find $(CALIBRATION_DIR) -type d) \
		   $(shell find $(MARKERMAP_DIR) -type d)

$(BUILD_DIR):
	@mkdir $@ $(patsubst %,$@/%, $(SRC_DIRS))

INCL += $(patsubst %,-I%, $(SRC_DIRS))

# Source configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

LIB_SRCS         := $(shell find $(LIB_DIR) -type f -name '*.cpp')
UTILS_SRCS       := $(shell grep -rLI --include='*.cpp' 'int main' $(UTILS_DIR))
CALIBRATION_SRCS := $(shell grep -rLI --include='*.cpp' 'int main' $(CALIBRATION_DIR))
MARKERMAP_SRCS   := $(shell grep -rLI --include='*.cpp' 'int main' $(MARKERMAP_DIR))

# Object configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

LIB_OBJS         := $(patsubst %.cpp,$(BUILD_DIR)/%.o, $(LIB_SRCS))
UTILS_OBJS       := $(patsubst %.cpp,$(BUILD_DIR)/%.o, $(UTILS_SRCS))
CALIBRATION_OBJS := $(patsubst %.cpp,$(BUILD_DIR)/%.o, $(CALIBRATION_SRCS))
MARKERMAP_OBJS   := $(patsubst %.cpp,$(BUILD_DIR)/%.o, $(MARKERMAP_SRCS))

$(BUILD_DIR)/%.o: %.cpp %.h | $(BUILD_DIR)
	@echo Compiling $<...
	@$(COMPILE) -c $< -o $@

$(BUILD_DIR)/%.o: %.cpp | $(BUILD_DIR)
	@echo Compiling $< \(no associated header found\)...
	@$(COMPILE) -c $< -o $@

# Libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

LIBRARY_NAME := aruco
LIBRARY := $(BUILD_DIR)/lib$(LIBRARY_NAME).a

$(LIBRARY): $(LIB_OBJS)
	@echo Linking static library $@...
	@ar rcs $@ $^

# Executables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

OPENCV_LIBS := -lopencv_calib3d -lopencv_contrib -lopencv_core \
               -lopencv_features2d -lopencv_flann -lopencv_gpu -lopencv_highgui \
	       -lopencv_imgproc -lopencv_legacy -lopencv_ml -lopencv_objdetect \
	       -lopencv_ocl -lopencv_photo -lopencv_stitching -lopencv_superres \
	       -lopencv_ts -lopencv_video -lopencv_videostab

LIBS := -L$(BUILD_DIR) -l$(LIBRARY_NAME) $(OPENCV_LIBS)


EXECUTABLES := $(shell grep -rlI --include='aruco*' -e 'int main' . \
                              | sed 's|^|$(BUILD_DIR)/|' | sed 's|\.cpp||' )

$(BUILD_DIR)/$(UTILS_DIR)/%: $(UTILS_DIR)/%.cpp $(UTILS_OBJS) $(LIBRARY)
	@echo Compiling and linking executable $@...
	@$(COMPILE) $< $(UTILS_OBJS) $(LIBS) -o $@

$(BUILD_DIR)/$(CALIBRATION_DIR)/%: $(CALIBRATION_DIR)/%.cpp $(CALIBRATION_OBJS) $(LIBRARY)
	@echo Compiling and linking executable $@...
	@$(COMPILE) $< $(CALIBRATION_OBJS) $(LIBS) -o $@

$(BUILD_DIR)/$(MARKERMAP_DIR)/%: $(MARKERMAP_DIR)/%.cpp $(MARKERMAP_OBJS) $(LIBRARY)
	@echo Compiling and linking executable $@...
	@$(COMPILE) $< $(MARKERMAP_OBJS) $(LIBS) -o $@

# Targets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

.DEFAULT_GOAL = $(LIBRARY)

all: $(LIBRARY) $(EXECUTABLES)

# Cleanup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

.PHONY: clean
clean:
	@rm -rf $(BUILD_DIR)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
