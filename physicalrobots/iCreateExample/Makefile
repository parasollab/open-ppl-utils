# Externals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Aruco
ARUCO_INCL  := -I../aruco/src
ARUCO_LIBS  := -L../aruco/build/ -laruco

# Player
PLAYER_DIR  := ../player/build/install
PLAYER_INCL := -I$(PLAYER_DIR)/include/player-3.0
PLAYER_LIBS := -L$(PLAYER_DIR)/lib64 -lplayerc++ -lplayerc -lplayercommon \
               -lplayercore -lplayerinterface -lplayerjpeg -lplayertcp \
	       -lplayerudp -lplayerwkb

# pthreads
THREAD_LIBS := -lpthread

# OpenCV
OPENCV_LIBS := -lopencv_calib3d -lopencv_core \
               -lopencv_features2d -lopencv_flann -lopencv_highgui \
	       -lopencv_imgproc -lopencv_ml -lopencv_objdetect \
	       -lopencv_photo -lopencv_stitching -lopencv_superres \
	       -lopencv_ts -lopencv_video -lopencv_videostab

# Boost
BOOST_LIBS  := -lboost_system -lboost_thread -lboost_signals
BOOST_DEFS  := -DBOOST_SIGNALS_NO_DEPRECATION_WARNING


# Aggregated externals
EXTERNAL_INCL := $(ARUCO_INCL) $(PLAYER_INCL)
EXTERNAL_DEFS := $(BOOST_DEFS)
EXTERNAL_LIBS := $(ARUCO_LIBS) $(PLAYER_LIBS) $(THREAD_LIBS) $(OPENCV_LIBS) \
                 $(BOOST_LIBS)


# Directory Layout ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

OBJ_DIR  := build
SRC_DIRS := $(patsubst ./%,%, \
	      $(shell find . -type d ! \( -name $(OBJ_DIR) -o -name '.' \) ))

$(OBJ_DIR):
	@mkdir $@

vpath %.h $(SRC_DIRS)
vpath %.cpp $(SRC_DIRS)


# Compiler Config ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

CXX   := g++
FLAGS := -std=c++11 -g3 -O2
WARN  := -Wall -Wno-deprecated-declarations -Wno-deprecated -Werror
DEFS  := $(EXTERNAL_DEFS)
INCL  := -I. $(EXTERNAL_INCL)

COMPILE := $(CXX) $(FLAGS) $(WARN) $(DEFS) $(INCL)


# Linker Config ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

LINKER    := $(CXX)
LINK_OPTS := -Wl --start-group
LIBS      := $(EXTERNAL_LIBS)

LINK := $(LINKER) $(LINKER_OPTS)


# Source Config ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

SRCS := $(shell find . -type f -name '*.cpp' ! -name 'main.cpp' -printf "%f ")


# Object File Config ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

OBJS  := $(patsubst %.cpp,$(OBJ_DIR)/%.o, $(SRCS))

$(OBJ_DIR)/%.o: %.cpp %.h | $(OBJ_DIR)
	@echo Compiling $<...
	@$(COMPILE) -c $< -o $@


# Executables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

EXECUTABLE := example
.DEFAULT_GOAL := $(EXECUTABLE)

$(EXECUTABLE): main.cpp $(OBJS)
	@echo Building executable $@...
	@$(COMPILE) $(LINKER_OPTS) $^ $(LIBS) -o $@


# Cleanup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

.PHONY: clean
clean:
	@rm -rf $(LIBRARY) $(EXECUTABLE) $(OBJ_DIR)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
