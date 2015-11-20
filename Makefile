# Make the HSE Integrator

SUNINCDIR = /home/eugene/local/sundials/instdir/include
SUNLIBDIR = /home/eugene/local/sundials/instdir/lib
LSEOSDIR  = /home/eugene/codes/astro/lseos_v2.7_f90

F90 	= gfortran
FFLAGS 	= -O3 #-ffpe-trap=invalid,zero,overflow
SUNINCLUDE = -I ${SUNINCDIR}
LSEINCLUDE = -I ${LSEOSDIR}
INCLUDE = 
F90LINK	= gfortran
CC	= gfortran

LLIBSS	= -L${SUNLIBDIR} -lsundials_fcvode -lsundials_cvode -lsundials_fnvecserial -lsundials_nvecserial
LLIBSM	= -L/usr/lib64 -llapack -lblas

SRCS 	= cvodehse.f90 cvode_indices.f90 integrator.f90 physical_constants.f90 polytrope_eos.f90 data_wrangler.f90 parameters.f90 ls_wrap_eos.f90 eos.f90 cvode_parameters.f90 eos_types.f90
OBJS	= cvodehse.o cvode_indices.o cvode_parameters.o data_wrangler.o el_eos_data.o eos_m4c_data.o eos.o eos_types.o force_data.o integrator.o lseos_v2.7.o ls_wrap_eos.o maxwel_data.o parameters.o physical_constants.o polytrope_eos.o

MAIN	= inthse

all:	${MAIN}
	@echo	Compilation complete!

%.mod : %.o
	@echo Compiled module.

cvodehse.o : cvodehse.f90 cvode_indices.mod eos.mod physical_constants.mod
	${F90} ${FFLAGS} -c $< -o $@

integrator.o : integrator.f90 polytrope_eos.mod cvodehse.mod cvode_indices.mod data_wrangler.mod parameters.mod cvode_parameters.mod
	${F90} ${FFLAGS} -c $< -o $@

cvode_parameters.o : cvode_parameters.f90
	${F90} ${FFLAGS} -c $< -o $@

eos_types.o : eos_types.f90
	${F90} ${FFLAGS} -c $< -o $@

eos.o : eos.f90 eos_types.mod ls_wrap_eos.mod polytrope_eos.mod
	${F90} ${FFLAGS} -c $< -o $@

ls_wrap_eos.o : ls_wrap_eos.f90 eos_types.mod physical_constants.mod
	${F90} ${FFLAGS} ${LSEINCLUDE} -c $< -o $@

data_wrangler.o : data_wrangler.f90 polytrope_eos.mod cvode_indices.mod
	${F90} ${FFLAGS} -c $< -o $@

parameters.o : parameters.f90 eos.mod cvode_parameters.mod
	${F90} ${FFLAGS} -c $< -o $@

cvode_indices.o : cvode_indices.f90
	${F90} ${FFLAGS} -c $< -o $@

physical_constants.o : physical_constants.f90
	${F90} ${FFLAGS} -c $< -o $@

polytrope_eos.o : polytrope_eos.f90 eos_types.mod
	${F90} ${FFLAGS} -c $< -o $@

$(MAIN): ${OBJS}
	@echo '${CC} ${INCLUDE} -o ${MAIN} ${OBJS} ${LLIBSS} ${LLIBSM}'
	${CC} ${INCLUDE} -o ${MAIN} ${OBJS} ${LLIBSS} ${LLIBSM}

clean:
	rm cvodehse.o  cvode_indices.o  cvode_parameters.o  data_wrangler.o  eos.o  eos_types.o  integrator.o  ls_wrap_eos.o  parameters.o  physical_constants.o  polytrope_eos.o  *.mod *~ ${MAIN}

depend: ${SRCS}
	makedepend ${INCLUDE} $^
