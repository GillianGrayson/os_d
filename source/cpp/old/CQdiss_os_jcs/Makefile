default: Makefile.inc build
	

# Is Makefile.inc available ?
Makefile.inc:
	@echo "######################################################################"
	@echo "# BEFORE COMPILING TEST, YOU SHOULD HAVE AN APPROPRIATE FILE"
	@echo "# Makefile.inc AVALAIBLE. "
	@echo "######################################################################"
	@exit 1

include Makefile.inc
	
build: 
	(cd ./read_config; make)
	(cd ./SparceM_op; make)
	(cd ./CQdiss; make)
	
clean:
	(cd ./read_config; make clean)
	(cd ./SparceM_op;  make clean)
	(cd ./CQdiss;      make clean)
