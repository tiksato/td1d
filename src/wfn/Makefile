include ../../config/makeconfig

all: general.lib x2e.lib hf.lib interface.lib

general.lib:
	cd $(@:.lib=); $(MAKE); cd ..
x2e.lib:
	cd $(@:.lib=); $(MAKE); cd ..
hf.lib:
	cd $(@:.lib=); $(MAKE); cd ..
interface.lib:
	cd $(@:.lib=); $(MAKE); cd ..

clean:
	rm -f *.o
fullclean:
	rm -f *.o *__genmod.* *~ */*.o */*__genmod.* */*~

.SUFFIXES: .f90 .o

.f90.o:
	$(F90) $(FLAG) -I$(TD1D_SRC)/global -c $<
	@echo
