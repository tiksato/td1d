include ../config/makeconfig

MAKE = make

all: exe

exe: mod lib
	cd ./main; $(MAKE); cd ..

mod:
	cd ./global; $(MAKE); cd ..

lib: file.lib field.lib grid.lib guess.lib init.lib input.lib prop.lib tdse.lib util.lib wfn.lib

file.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
field.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
grid.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
guess.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
init.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
input.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
prop.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
tdse.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
util.lib: 
	cd $(@:.lib=); $(MAKE); cd ..
wfn.lib: 
	cd $(@:.lib=); $(MAKE); cd ..

clean:
	rm -f */*.o */*/*.o */*.mod
fullclean:
	rm -f */*.o */*/*.o */*.mod */*__genmod.* */*/*__genmod.* *~ */*~ */*/*~
