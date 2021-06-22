# Compile an external processes, and amend the process tables
DFM_ROOT=/home/rusty/src/dfm/r68819-dbg

LT=$(DFM_ROOT)/bin/libtool

F90=$(DFM_ROOT)/bin/mpif90
FFLAGS=-DWITH_DELFTONLINE -g -O0 -ffixed-line-length-none

all: liblocalprocess.so proc_def.def

%.lo: %.f
	/bin/bash $(LT) --tag=F77 --mode=compile $(F90) $(FFLAGS) -c $< -o $@

LOCAL_INSTALL=$(realpath .)

# Figure out how to get libtool to churn out a shared library.
liblocalprocess.so: ageconc.lo
	/bin/bash $(LT) --tag=F77 --mode=link $(F90) $(FFLAGS) $< -o liblocalprocess.la -ldl -lrt -rpath $(LOCAL_INSTALL)
	/bin/bash $(LT) --tag=F77 --mode=install install -c liblocalprocess.la $(LOCAL_INSTALL)/liblocalprocess.la

test_compile:
	gfortran $(FFLAGS) -c ageconc.f

test_compile2:
	gfortran $(FFLAGS) -c nitrif.f

ORIG_CSV_DIR=../waq_tables

CSVS=inputs.csv disps.csv grpsub.csv fortran.csv con_sub.csv con_pro.csv velocs.csv stochi.csv outputs.csv outpflx.csv items.csv \
 config.csv

%.csv: loc-%.csv
	cp $(ORIG_CSV_DIR)/$@ .
	-tail -n+2 $< | grep -v -e '#' >> $@

proc_def.def: $(CSVS)
	$(DFM_ROOT)/bin/waqpbexport


.PHONY: clean

clean:
	-rm $(CSVS)
	-rm liblocalprocess.so* liblocalprocess.a liblocalprocess.la
	-rm proc_def.def proc_def.dat
	-rm ageconc.lo ageconc.o
	-rm -r .libs
	-rm *.pptex coefedit.dat
	-rm *.log
	-rm procesm.asc

