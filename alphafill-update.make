# Makefile to create and/or update a local alphafill databank
#
# This makefile assumes alphafill is installed and can be found
# in the path. If not, change the ALPHAFILL variable below or
# specify it on the command line

SHELL = /bin/bash
ALPHAFILL = /home/maarten/.local/bin/alphafill

firstTarget: update

CIF_DIR     = /DATA/afdb-v4/
DEST_DIR    = /DATA/af-filled-ng/

SUB_CIF_DIRS    = $(wildcard $(CIF_DIR)/??)

empty=
space=$(empty) $(empty)
VPATH += $(subst $(space),:,$(SUB_CIF_DIRS))

CIF_IDS     = $(wildcard $(CIF_DIR)/*/*.cif.gz)
CIF_IDS     := $(subst $(CIF_DIR),$(DEST_DIR),$(CIF_IDS))
CIF_IDS     := $(sort $(CIF_IDS))

TARGETS = $(CIF_IDS:%-model_v4.cif.gz=%-filled_v4.cif.gz)

.PHONY: test
test:
	echo $(TARGETS)

.PHONY: update
update: $(TARGETS)

$(DEST_DIR)/%-filled_v4.cif.gz: $(CIF_DIR)/%-model_v4.cif.gz
	$(ALPHAFILL) process --quiet --threads=1 $? $@ |& tee $(@:%.cif.gz=%.log)

