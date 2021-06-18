# Compile an external processes, and amend the process tables
DFM_ROOT=/home/rusty/src/dfm/r68819-dbg

LT=$(DFM_ROOT)/libtool

F90=$(DFM_ROOT)/bin/mpif90
FFLAGS=-DWITH_DELFTONLINE -g -O0 -ffixed-line-length-none

all: local_proc.so proc_def.def


%.lo: %.f90
	$(LT) --tag=F77 --mode= $(F90) $(FFLAGS) -c -o $< $@

# /bin/bash ../../../../../../libtool  --tag=F77   --mode=compile /home/rusty/src/dfm/r68819-dbg/bin/mpif90  -I../../../../../../engines_gpl/waq/include -I../../../../../../engines_gpl/waq/packages/waq_utils_f/src/mod_timings -I../../../../../../engines_gpl/waq/packages/waq_utils_f/src/waq_utils -I../../../../../../engines_gpl/waq/packages/waq_process/src/bloom -I../../../../../../utils_lgpl/deltares_common/packages/deltares_common/src -DWITH_DELFTONLINE -g -O0 -ffixed-line-length-none -c -o nitrif.lo nitrif.f

test_compile:
	gfortran $(FFLAGS) -c ageconc.f

test_compile2:
	gfortran $(FFLAGS) -c nitrif.f

ORIG_CSV_DIR=../waq_tables

csvs: inputs.csv disps.csv grpsub.csv fortran.csv con_sub.csv con_pro.csv velocs.csv stochi.csv outputs.csv outpflx.csv items.csv \
 config.csv

%.csv: loc-%.csv
	cp $(ORIG_CSV_DIR)/$@ .
	-tail -n+2 $< | grep -v -e '#' >> $@

proc_def.def: csvs
	$(DFM_ROOT)/bin/waqpbexport

# Failing around the velocities I think.
# veloit(ivelo) 'VxUpAgeC00'
# subsid(isubs) $7 = 'AgeC00    '
# velosu(ivelo) $8 = 'AgeC00    '
# Could be lack of ^M endings?
#  Probably that the velocities need to be
#  in items, maybe in outputs, too?
