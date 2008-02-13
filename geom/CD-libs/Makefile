all: rapid pqp vclip

rapid:
	@echo "compiling rapid..."
	@cd RAPID; ${MAKE} 
	@echo "rapid done"

pqp:
	@echo "compiling pqp..."
	@cd PQP_v1.2; ${MAKE}
	@echo "pqp done"

qhull:
	@echo "compiling qhull..."
	@cd qhull3.0/src; ${MAKE}
	@echo "qhull done"

vclip: qhull
	@echo "compiling vclip..."
	@cd VCLIP; ${MAKE}
	@echo "vclip done"

clean: clean_rapid clean_pqp clean_qhull clean_vclip

clean_rapid:
	@echo "cleaning rapid..."
	@cd RAPID; ${MAKE} clean
	@echo "rapid done"

clean_pqp:
	@echo "cleaning pqp..."
	@cd PQP_v1.2; ${MAKE} clean
	@echo "pqp done"

clean_qhull:
	@echo "cleaning qhull..."
	@cd qhull3.0/src; ${MAKE} clean
	@echo "qhull done"

clean_vclip:
	@echo "cleaning vclip..."
	@cd VCLIP; ${MAKE} clean
	@echo "vclip done"
