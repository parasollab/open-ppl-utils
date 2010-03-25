all: rapid pqp vclip solid

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

SOLIDLIB := ${PWD}/solid-3.5.6

solid: 
	@echo "comiling solid..."
	@cd solid-3.5.6; ./configure --prefix=${SOLIDLIB} --exec-prefix=${SOLIDLIB}; ${MAKE}; ${MAKE} -f MakefileSpecial cleanlibrary; ${MAKE} -f MakefileSpecial library
	@echo "solid done"

clean: clean_rapid clean_pqp clean_qhull clean_vclip clean_solid

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

clean_solid:
	@echo "cleaning solid..."
	@cd solid-3.5.6; ${MAKE} clean; ${MAKE} -f MakefileSpecial cleanlibrary
	@echo "solid done"

