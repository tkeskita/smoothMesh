sinclude $(GENERAL_RULES)/module-path-user

/* Failsafe - user location */
ifeq (,$(strip $(FOAM_MODULE_APPBIN)))
    FOAM_MODULE_APPBIN = $(FOAM_USER_APPBIN)
endif
ifeq (,$(strip $(FOAM_MODULE_LIBBIN)))
    FOAM_MODULE_LIBBIN = $(FOAM_USER_LIBBIN)
endif

EXE_INC = \
    $(VERSION_SPECIFIC_INC) \
    -I$(LIB_SRC)/finiteVolume/lnInclude \
    -I$(LIB_SRC)/meshTools/lnInclude \
    -I$(LIB_SRC)/surfMesh/lnInclude \

EXE_LIBS = \
    -lfiniteVolume \
    -lmeshTools \
    -lsurfMesh
