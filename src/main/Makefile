include ../../config/makeconfig

OBJ1 = td1d.o
EXE1 = td1d.exe
all: exe1
exe1: $(OBJ1)
	$(F90) $(FLAG) \
	../global/*.o \
	../file/*.o \
	../field/*.o \
	../grid/*.o \
	../guess/*.o \
	../init/*.o \
	../input/*.o \
	../prop/*.o \
	../tdse/*.o \
	../util/*.o \
	../wfn/general/*.o \
	../wfn/x2e/*.o \
	../wfn/hf/*.o \
	../wfn/interface/*.o \
	td1d.o \
	$(LIB) -o $(EXE1)
	@echo
clean:
	rm -f $(OBJ)

fullclean:
	rm -f $(OBJ) *__genmod* *~

.SUFFIXES: .F90 .o

.F90.o:
	$(F90) $(FLAG) -I$(TD1D_SRC)/global -c $<
	@echo

