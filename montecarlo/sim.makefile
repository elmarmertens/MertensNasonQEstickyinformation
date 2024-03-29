THIS ?= cambridgesimSIthetaTVPlambdaTVP

T ?= 200
Ngrid ?= 320
doInflationNoise ?= 1
Nparticles ?= 10000
DATALABEL ?= cambridge2018GDPD
FCmode ?= prod
doSmoother ?= 0
Nsmoother ?= 1000
smootherNparticles ?= 1000
doSecondResamplingStep ?= 0
doStoreDGP ?= 0
dofAPFfattail ?= 0

toolboxes=vslbox.o embox.o blaspackbox.o timerbox.o  gibbsbox.o statespacebox.o  densitybox.o cambridgebox.o

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
  # mac
  # for ifort 2020: add -warn noexternals
  FCdebugseq=ifort -mkl -warn all -WB -check all -check noarg_temp_created -static-intel -Wl,-stack_size,0x20000000 -Wl,-rpath,$(MKLROOT)/../compiler/lib/
  FCdebug=ifort -mkl -warn all -WB -check all -check noarg_temp_created -static-intel -Wl,-stack_size,0x20000000  -qopenmp -Wl,-rpath,$(MKLROOT)/../compiler/lib/
  FCprod=ifort -O3 -mkl -nocheck -qopenmp -static-intel -xHost -Wl,-stack_size,0x80000000 -Wl,-rpath,$(MKLROOT)/../compiler/lib/
else
  # linux
  FCdebugseq=ifort -mkl -warn all -WB -check all -check noarg_temp_created  -shared-intel
  FCdebug=ifort -mkl -warn all -WB -check all -check noarg_temp_created  -shared-intel -qopenmp
  FCprod=ifort -O3 -mkl -nocheck -qopenmp  -shared-intel  -xHost
endif


ifeq ($(FCmode),debug)
  FC=$(FCdebug)
else ifeq ($(FCmode),debugseq)
  FC=$(FCdebugseq)
else
  FC=$(FCprod)
endif


toolboxdir=../fortranbox/
vslbox=INTELvslbox
timerbox=OMPtimerbox

main  : $(THIS)
toolboxes : $(toolboxes)

$(THIS)  : $(THIS).f90 $(toolboxes)
	$(FC) $(THIS).f90 -o $(THIS) $(toolboxes)




cambridgebox.o  : ../main/cambridgebox.f90 embox.o vslbox.o
	$(FC) -c ../main/cambridgebox.f90

statespacebox.o :: $(toolboxdir)statespacebox.f90 blaspackbox.o embox.o vslbox.o
	$(FC) -c $(toolboxdir)statespacebox.f90

embox.o  : $(toolboxdir)embox.f90 vslbox.o
	$(FC) -c $(toolboxdir)embox.f90

densitybox.o  : $(toolboxdir)densitybox.f90 embox.o blaspackbox.o vslbox.o statespacebox.o
	$(FC) -c $(toolboxdir)densitybox.f90

blaspackbox.o  : $(toolboxdir)blaspackbox.f90  embox.o
	$(FC) -c $(toolboxdir)blaspackbox.f90

timerbox.o : $(toolboxdir)$(timerbox).f90 embox.o
	$(FC) -c $(toolboxdir)$(timerbox).f90 -o timerbox.o

gibbsbox.o  : $(toolboxdir)gibbsbox.f90 timerbox.o blaspackbox.o statespacebox.o
	$(FC) -c $(toolboxdir)gibbsbox.f90

vslbox.o  : $(toolboxdir)$(vslbox).f90
	$(FC) -c $(toolboxdir)$(vslbox).f90 -o vslbox.o

compile	: $(THIS)

run	: $(THIS)
	rm -f *.debug
	time -p ./$(THIS) $(Ngrid) $(T) $(doInflationNoise) $(Nparticles) $(doStoreDGP) $(doSmoother) $(Nsmoother) $(smootherNparticles) $(doSecondResamplingStep) $(dofAPFfattail)


edit :
	aquamacs $(THIS).f90

clean	:
	rm -f $(THIS)
	rm -f *.debug
	rm -f *_genmod.f90
	rm -f *.log
	rm -f *.mod
	rm -f *.o
	rm -f *~

cleanall :
	rm -f *.dat
	$(MAKE) clean
