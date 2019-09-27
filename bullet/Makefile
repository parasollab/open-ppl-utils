# This is a weird makefile that runs a bash script to set up bullet.
# Some day we should unify these.

platform ?= LINUX_gcc

# Script Launcher ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

BUILD_DIR := build
LIB_DIR   := $(BUILD_DIR)/install/lib
$(LIB_DIR):
	@echo Building bullet library...
	@sh parasol_setup.sh $(platform)

# Cleanup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

.PHONY: clean
clean:
	@echo Cleaning bullet library...
	@rm -rf $(BUILD_DIR)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
