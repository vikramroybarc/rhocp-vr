###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
#
###############################################################################
# Use the MOOSE submodule if it exists and MOOSE_DIR is not set
print-%:
	@echo '$* = $($*)'

MOOSE_SUBMODULE    := $(CURDIR)/moose
ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
  MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
else
  MOOSE_DIR        ?= $(shell dirname `pwd`)/moose
endif

# framework
FRAMEWORK_DIR      := $(MOOSE_DIR)/framework
include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

################################## MODULES ####################################
# To use certain physics included with MOOSE, set variables below to
# yes as needed.  Or set ALL_MODULES to yes to turn on everything (overrides
# other set variables).

ALL_MODULES                := no

CHEMICAL_REACTIONS         := no
CONTACT                    := no
ELECTROMAGNETICS           := no
EXTERNAL_PETSC_SOLVER      := no
FLUID_PROPERTIES           := no
FSI                        := no
FUNCTIONAL_EXPANSION_TOOLS := no
GEOCHEMISTRY               := no
HEAT_CONDUCTION            := no
LEVEL_SET                  := no
MISC                       := no
NAVIER_STOKES              := no
PERIDYNAMICS               := no
PHASE_FIELD                := yes
POROUS_FLOW                := no
RAY_TRACING                := no
REACTOR                    := no
RDG                        := no
RICHARDS                   := no
STOCHASTIC_TOOLS           := no
THERMAL_HYDRAULICS         := no
HEAT_TRANSFER              := yes
SOLID_MECHANICS            := yes
XFEM                       := no

include $(MOOSE_DIR)/modules/modules.mk
###############################################################################

# ----------------------- Couple with Raccoon (external app) ------------------
# Detect a raccoon submodule under this repo; otherwise assume a sibling checkout.
RACCOON_SUBMODULE  := $(CURDIR)/raccoon-vr
ifneq ($(wildcard $(RACCOON_SUBMODULE)/Makefile),)
  RACCOON_DIR      ?= $(RACCOON_SUBMODULE)
else
  # Meaning: set RACCOON_DIR to "<parent of current>/raccoon" if not provided
  RACCOON_DIR      ?= $(shell dirname `pwd`)/raccoon-vr
endif

# Raccoon headers and library
RAC_INCLUDE_DIRS   := $(RACCOON_DIR)/include
RAC_LIB_DIR        := $(RACCOON_DIR)/lib
RAC_LIB_FILE       := $(RAC_LIB_DIR)/libraccoon-$(METHOD).la
# If your build produces a static archive instead, use:
# RAC_LIB_FILE     := $(RAC_LIB_DIR)/libraccoon-$(METHOD).a

ADDITIONAL_INCLUDES += $(addprefix -I, $(shell find $(RACCOON_DIR)/include -type d))
ADDITIONAL_LIBS     += $(RAC_LIB_FILE)
EXTERNAL_LIBS       += $(RAC_LIB_FILE)

# Add runtime path so you don’t need LD_LIBRARY_PATH
ADDITIONAL_LDFLAGS += -Wl,-rpath,$(RACCOON_DIR)/lib

# (Optional) auto-build the Raccoon lib if missing
$(RAC_LIB_FILE):
	$(MAKE) -C $(RACCOON_DIR) BUILD_LIB:=yes METHOD=$(METHOD) lib
# ---------------------------------------------------------------------------



# dep apps (this application)
APPLICATION_DIR    := $(CURDIR)
APPLICATION_NAME   := rhocp
BUILD_EXEC         := yes
# If you also want to produce a reusable library for RhoCP:
# BUILD_LIB        := yes
GEN_REVISION       := no

include            $(FRAMEWORK_DIR)/app.mk

###############################################################################
# Additional special case targets should be added here
