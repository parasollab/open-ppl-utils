all: rapid pqp vclip

rapid:
	@echo "compiling rapid..."
	@cd RAPID; gmake library
	@echo "rapid done"

pqp:
	@echo "compiling pqp..."
	@cd PQP_v1.2; gmake library
	@echo "pqp done"

qhull:
	@echo "compiling qhull..."
	@cd qhull3.0/src; gmake all
	@echo "qhull done"

vclip: qhull
	@echo "compiling vclip..."
	@cd VCLIP; gmake all
	@echo "vcliip done"

clean: clean_rapid clean_pqp clean_qhull clean_vclip

clean_rapid:
	@echo "cleaning rapid..."
	@cd RAPID; gmake clean
	@echo "rapid done"

clean_pqp:
	@echo "cleaning pqp..."
	@cd PQP_v1.2; gmake clean
	@echo "pqp done"

clean_qhull:
	@echo "cleaning qhull..."
	@cd qhull3.0/src; gmake cleanall
	@echo "qhull done"

clean_vclip:
	@echo "cleaning vclip..."
	@cd VCLIP; gmake clean
	@echo "vclip done"
