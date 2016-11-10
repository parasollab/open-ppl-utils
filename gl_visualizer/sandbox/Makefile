# Externals ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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

# Compiler configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

CXX   := g++
FLAGS := -std=c++11 -g3 -O2
WARN  := -Wall -Werror

INCL := -I. $(QT_INCL) $(NONSTD_INCL) $(GLUTILS_INCL)

# Linker configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

LIBS := $(QT_LIBS) $(GLUTILS_LIBS) $(GL_LIBS) $(NONSTD_LIBS)

# Qt4 MOC configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

MOC := moc

# Directory Layout ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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

# Source file configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

SRCS     := $(patsubst ./%,%, $(shell find . -type f -name '*.cpp' \
		    ! -name 'moc_*.cpp' ! -name 'main.cpp' -printf "%f "))
MOC_SRCS := $(patsubst ./%.h,moc_%.cpp, $(shell find . -type f \
		    -name '*.h' -execdir grep -rl Q_OBJECT \{} + ))

$(MOC_DIR)/moc_%.cpp: %.h | $(MOC_DIR)
	@echo Moc\'ing $(shell echo $< | sed 's|.*/||' )...
	@$(MOC) $< -o $@

# Object file configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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

NAME   := glsandbox
SHARED := lib$(NAME).o
STATIC := lib$(NAME).a

.PHONY: shared
shared: $(SHARED)
$(SHARED): $(OBJS) $(MOC_OBJS) $(patsubst %,$(MOC_DIR)/%,$(MOC_SRCS))
	@echo Linking shared library $(SHARED)...
	@$(CXX) -shared $^ -o $@

.PHONY: static
static: $(STATIC)
$(STATIC): $(OBJS) $(MOC_OBJS) $(patsubst %,$(MOC_DIR)/%,$(MOC_SRCS))
	@echo Linking static library $(STATIC)...
	@ar rcs $@ $^

# Executable ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

EXECUTABLE := sandbox

.DEFAULT_GOAL := $(EXECUTABLE)

# Use the static library to build the executable.
$(EXECUTABLE): $(OBJ_DIR)/main.o $(STATIC)
	@cd $(NONSTD_DIR) ; make
	@cd $(GLUTILS_DIR) ; make
	@echo Linking $@...
	@$(CXX) $< -L. -l$(NAME) $(LIBS) -o $@

# Clean ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

.PHONY: clean
clean: cleandocs
	@rm -rf $(EXECUTABLE) $(OBJ_DIR) $(STATIC) $(SHARED)
	@cd $(GLUTILS_DIR) ; make clean
	@cd $(NONSTD_DIR) ; make clean

.PHONY: cleandocs
cleandocs:
	@rm -rf $(DOC_DIR)/html

# Docs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

.PHONY: docs
docs:
	@echo Building docs...
	@doxygen $(DOC_DIR)/Doxyfile

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
