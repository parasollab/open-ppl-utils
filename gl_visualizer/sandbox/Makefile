# Externals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Set up externals here.

# Qt
BASE         := /usr/include
QT_INCL      := -I$(BASE)/Qt -I$(BASE)/QtCore -I$(BASE)/QtGui -I$(BASE)/QtOpenGL
QT_LIBS      := -lQtCore -lQtGui -lQtOpenGL

# GL
GL_LIBS      := -lGLU -lGL

# nonstd
NONSTD_DIR   := ../nonstd
NONSTD_INCL  := -I$(NONSTD_DIR)
NONSTD_LIBS  := -L$(NONSTD_DIR) -lnonstd

# glutils
GLUTILS_DIR  := ../glutils
GLUTILS_INCL := -I$(GLUTILS_DIR)
GLUTILS_LIBS := -L$(GLUTILS_DIR) -lglutils

# Aggregated externals.
EXTERNAL_INCL := -I. $(QT_INCL) $(NONSTD_INCL) $(GLUTILS_INCL)
EXTERNAL_DEFS :=
EXTERNAL_LIBS := $(QT_LIBS) $(GLUTILS_LIBS) $(GL_LIBS) $(NONSTD_LIBS)


# Directory Layout ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Configure the directory tree.
# - Object files will be placed in OBJ_DIR.
# - Moc'd sources will go into MOC_DIR.
# - Source code is in SRC_DIR.
# - Documentation will be built in DOC_DIR.

OBJ_DIR  := build
MOC_DIR  := $(OBJ_DIR)/moc
SRC_DIRS := $(patsubst ./%,%, $(shell find . -type d \
		    ! \( -name $(OBJ_DIR) -o -path ./$(MOC_DIR) -o -name '.' \) ))
DOC_DIR  := docs

$(OBJ_DIR):
	@mkdir $@
$(MOC_DIR): | $(OBJ_DIR)
	@mkdir $@

vpath %.h       $(SRC_DIRS)
vpath %.cpp     $(SRC_DIRS)
vpath moc_%.cpp $(MOC_DIR)


# Dependency Tracking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Set up automatic dependency tracking.

DEP_FILE   := $(OBJ_DIR)/Dependencies
DEP_GLOB    = $(OBJ_DIR)/$*.d
TRACK_DEPS  = cat $(DEP_GLOB) >> $(DEP_FILE) ; rm -f $(DEP_GLOB)

-include $(DEP_FILE)


# Compiler configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Configure compiler, flags, and options.

CXX   := g++
FLAGS := -std=c++11 -fPIC -g3 -O3
WARN  := -Wall -Werror -Wextra
DEFS  := $(EXTERNAL_DEFS)
INCL  := -I. $(EXTERNAL_INCL)
DEPS   = -MMD -MF $(DEP_GLOB)

COMPILE = $(CXX) $(FLAGS) $(WARN) $(DEFS) $(INCL) $(DEPS)


# Linker configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Configure linker and objects.

LINKER    := $(CXX)
LINK_OPTS :=
LIBS      := $(EXTERNAL_LIBS)

LINK := $(LINKER) $(LINK_OPTS)


# Qt MOC configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Set the Qt4 MOC command.

MOC := moc


# Source file configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Scan the source and test directories to look for files with suffix '.cpp'.

SRCS     := $(patsubst ./%,%, $(shell find . -type f -name '*.cpp' \
		    ! -name 'moc_*.cpp' ! -name 'main.cpp' -printf "%f "))
MOC_SRCS := $(patsubst ./%.h,moc_%.cpp, $(shell find . -type f \
		    -name '*.h' -execdir grep -rl Q_OBJECT \{} + ))

$(MOC_DIR)/moc_%.cpp: %.h | $(MOC_DIR)
	@echo Moc\'ing $(shell echo $< | sed 's|.*/||' )...
	@$(MOC) $< -o $@


# Object file configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Create an object name for each source file and define object compilation
# recipes.

OBJS     := $(patsubst %.cpp,$(OBJ_DIR)/%.o, $(SRCS))
MOC_OBJS := $(patsubst %.cpp,$(OBJ_DIR)/%.o, $(MOC_SRCS))

$(OBJ_DIR)/moc_%.o: $(MOC_DIR)/moc_%.cpp %.h | $(MOC_DIR)
	@echo Compiling $(shell echo $< | sed 's|.*/||' )...
	@$(CXX) $(FLAGS) $(WARN) $(INCL) -c $< -o $@

$(OBJ_DIR)/%.o: %.cpp %.h | $(OBJ_DIR)
	@echo Compiling $(shell echo $< | sed 's|.*/||' )...
	@$(CXX) $(FLAGS) $(WARN) $(INCL) -c $< -o $@

$(OBJ_DIR)/%.o: %.cpp | $(OBJ_DIR)
	@echo Compiling $(shell echo $< | sed 's|.*/||' ) \
	      \(no associated header\)...
	@$(CXX) $(FLAGS) $(WARN) $(INCL) -c $< -o $@


# Libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Define the output libraries and recipes.

NAME   := glsandbox
SHARED := lib$(NAME).o
STATIC := lib$(NAME).a

.PHONY: shared
shared: $(SHARED)
$(SHARED): $(OBJS) $(MOC_OBJS) $(patsubst %,$(MOC_DIR)/%,$(MOC_SRCS))
	@echo Linking shared library $(SHARED)...
	@$(LINK) -shared $^ -o $@

.PHONY: static
static: $(STATIC)
$(STATIC): $(OBJS) $(MOC_OBJS) $(patsubst %,$(MOC_DIR)/%,$(MOC_SRCS))
	@echo Linking static library $(STATIC)...
	@ar rcs $@ $^


# Executable ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Define the example executable.

EXECUTABLE := sandbox
.DEFAULT_GOAL := $(EXECUTABLE)

$(EXECUTABLE): $(OBJ_DIR)/main.o $(STATIC)
	@+cd $(NONSTD_DIR) && $(MAKE)
	@+cd $(GLUTILS_DIR) && $(MAKE)
	@echo Linking $@...
	@$(LINK) $< -L. -l$(NAME) $(LIBS) -o $@


# Docs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Define recipe for building the documentation.

.PHONY: docs
docs: cleandocs
	@echo Building docs...
	@doxygen $(DOC_DIR)/Doxyfile


# Clean ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# Define cleanup recipes.

.PHONY: clean
clean: cleandocs
	@rm -rf $(EXECUTABLE) $(OBJ_DIR) $(STATIC) $(SHARED)
	@cd $(GLUTILS_DIR) && $(MAKE) clean
	@cd $(NONSTD_DIR) && $(MAKE) clean

.PHONY: cleandocs
cleandocs:
	@rm -rf $(DOC_DIR)/html


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
