# This makefile is a dispatcher for the four sub-units in this utility.

.PHONY: all
all:
	@cd player && $(MAKE)


.PHONY: clean
clean:
	@cd aruco && $(MAKE) clean
	@cd player && $(MAKE) clean
	@cd iCreateExample && $(MAKE) clean
	@cd detector_server && $(MAKE) clean
	@cd calibration && $(MAKE) clean
