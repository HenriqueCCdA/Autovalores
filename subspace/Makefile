#!/bin/sh
# $URL: svn+ssh://henrique@146.164.26.19/mnt/shared/repositorio/mefpar/branches/geoacustica/mef/Makefile $
# $Id: Makefile 2 2009-04-03 17:46:50Z helano $

#make
#------------------Variaveis Globais-----------------------
PATH_MOD=""
PATH_INCLUDE=""    
PRENAME=subsp
FC=gfortran
OPENMP=yes
OT=O3
#------------------gerando o nome do excutavel-------------
ifeq ($(FC),ifort)
  COMPILER_NAME=intel_${OT}
endif

ifeq ($(FC),gfortran)
  COMPILER_NAME=gnu_${OT}
endif
  
NAME+=$(PRENAME)_$(COMPILER_NAME)
#-------------------Fontes---------------------------------
src_mef = src/Dsolv.for\
src/Subsp.for\
src/mmio.for\
src/Coo.for\
src/Vmop.for


#mod_mef=../src/Malloc.f
#-------------------Flags necessarios--------------------------------
NFLAGS=    
LDFLAGS=
#--------------------compiladores------------------------------------
# intel ifort
ifeq ($(FC),ifort)
  LDFLAGS += 
  OFLAGS  +=   
  ifeq ($(OPENMP),yes)     
    OFLAGS  += -openmp
  endif
endif
# gnu gcc
ifeq ($(FC),gfortran)
  LDFLAGS += 
  OFLAGS  += -Wall 
  ifeq ($(OPENMP),yes)
    OFLAGS  += -fopenmp
  endif
endif
#--------------------------------------------------------------------
#---------------------------Debug------------------------------------
ifeq ($(DEBUG),yes)
  OFLAGS += -g  	
else
  OFLAGS += -${OT} 
endif
#--------------------------------------------------------------------
FFLAGS= $(NFLAGS) $(OFLAGS) 

.SUFFIXES: 
.SUFFIXES: .for .fi .o

objs_src= $(patsubst %.for,%.o,$(src_mef))

build:	$(objs_src)	
	ls bin || mkdir -p bin
	$(FC) $(FFLAGS) $(objs_src) -o bin/$(NAME) $(LDFLAGS)   

$(objs_src): %.o:%.for
	$(FC) $(FFLAGS) -c $< -o $@

tags:
	ctags -R src/*.f include/*.fi

.PHONY: cleantags
cleantags:
	@rm -fv tags
	
.PHONY: clean
clean:  
	@rm -fv src/*.o
	@rm -fv ../src/*.o
	@rm -fv ../src/*.mod
	@rm -fv bin/$(NAME)

.PHONY: cleanall
cleanall:  
	@rm -fv tags
	@rm -fv src/*.o
	@rm -fv ../src/*.o
	@rm -fv ../src/*.mod
	@rm -fv bin/$(NAME)

.PHONY: cleanmod
cleanmod:  
	@rm -fv ../src/*.mod

.PHONY: help
help:
	@echo "Autor :$(AUTHOR)                              "
	@echo "Makefile para prepar para sitemas linux.      "
	@echo -e "\E[7;32mOpcoes:\E[1;0m                      "
	@echo "build         - compila o prepar              "
	@echo "build_modules - gera os modulos               "
	@echo "tags          - gera os tags                  "
	@echo "cleantags     - limpa os tags                 "
	@echo "clean         - limpa os obj, bin e mod       "
	@echo "cleaall       - limpa tudo obj,bin,mod e tags "
	@echo "cleanmod      - limpa os modulos              "

# DO NOT DELETE

# DO NOT DELETE
