include ../make.cfg/flags.mak

PROGNAME=frb_search
# OK : -lcfitsio  -L/home/msok/askap/craft//lib/ -lfrblib  -L/opt/caastro/bighorns//lib -lfitslib -lnova
CRAFT=/home/msok/askap/craft/
CMN_LIBS=-L$(CRAFT)/lib/ -lfrblib $(COMMON_LIBS)

#  -I ../craft/tests/ - Sigproc (filterbank file)
INCLUDES=-D_UNIX -I$(SRCDIR)/cmn/baselib/ -I$(SRCDIR)/cmn/mathlib -I../fitslib/ -I ../craft/tests/dada2fil/

all : create_dynaspec dynaspec_search calcstat inject_frb calcstat read_data_cube find_sources_fits
	g++ mwa_fits.cpp -c $(INCLUDES)   $(OPT)        
	g++ sources.cpp -c $(INCLUDES)   $(OPT)
	g++ dedisp_search.cpp -c $(INCLUDES)   $(OPT)
	g++ main.cpp mwa_fits.o dedisp_search.o sources.o -o  $(PROGNAME) $(OPT) -I /opt/caastro/bighorns/software/common -I../fitslib/ $(INCLUDES) -lcfitsio -L$(BIGHORNS)/lib -lfitslib $(CMN_LIBS)
	cp $(PROGNAME) $(BIGHORNS)/bin/
	cp calcstat $(BIGHORNS)/bin/
	cp inject_frb  $(BIGHORNS)/bin/
	cp dynaspec_search $(BIGHORNS)/bin/
	cp create_dynaspec $(BIGHORNS)/bin/
	chmod +x *.sh
	cp *.sh $(BIGHORNS)/bin/
	cp filter*.py $(BIGHORNS)/bin/
	cp read_data_cube $(BIGHORNS)/bin/

test :
	g++ main_test.cpp mwa_fits.cpp dedisp_search.cpp sources.cpp -o  dynspec_test $(OPT) -I /opt/caastro/bighorns/software/common -I../fitslib/ $(INCLUDES) -lcfitsio -L$(BIGHORNS)/lib -lfitslib $(CMN_LIBS)
	cp dynspec_test $(BIGHORNS)/bin/

show_dm_pixels : main_show_dm_pixels.cpp mwa_fits.cpp dedisp_search.cpp sources.cpp 
	g++ main_show_dm_pixels.cpp mwa_fits.cpp dedisp_search.cpp sources.cpp -o  show_dm_pixels $(OPT) -I /opt/caastro/bighorns/software/common -I../fitslib/ $(INCLUDES) -lcfitsio -L$(BIGHORNS)/lib -lfitslib $(CMN_LIBS)
	cp show_dm_pixels $(BIGHORNS)/bin/

# program searching for signal in dynamic spectrum fits files :
dynaspec_search : main_test_dynaspec.cpp dedisp_search.cpp mwa_fits.cpp
	g++ main_test_dynaspec.cpp mwa_fits.cpp dedisp_search.cpp sources.cpp -o  dynaspec_search $(OPT) -I /opt/caastro/bighorns/software/common -I../fitslib/ $(INCLUDES) -lcfitsio  $(CMN_LIBS) -L$(BIGHORNS)/lib -lfitslib
	cp dynaspec_search $(BIGHORNS)/bin/


create_dynaspec : mwa_fits.cpp sources.cpp dedisp_search.cpp main_dynaspec.cpp myprogress.cpp
	g++ mwa_fits.cpp -c $(INCLUDES)   $(OPT)        
	g++ sources.cpp -c $(INCLUDES)   $(OPT)
	g++ dedisp_search.cpp -c $(INCLUDES)   $(OPT)
	g++ myprogress.cpp -c $(INCLUDES)   $(OPT)
	g++ main_dynaspec.cpp mwa_fits.o dedisp_search.o sources.o myprogress.o -o create_dynaspec $(OPT) -I /opt/caastro/bighorns/software/common -I../fitslib/ $(INCLUDES) -lcfitsio -L$(BIGHORNS)/lib -lfitslib $(CMN_LIBS)
	cp create_dynaspec $(BIGHORNS)/bin/


calcstat : calcstat.cpp mwa_fits.cpp
	g++ mwa_fits.cpp -c $(INCLUDES)   $(OPT)   
	g++ sources.cpp -c $(INCLUDES)   $(OPT)
	g++ calcstat.cpp mwa_fits.o sources.o -o calcstat  -I /opt/caastro/bighorns/software/common -I../fitslib/ $(INCLUDES) -lcfitsio -L$(BIGHORNS)/lib -lfitslib $(CMN_LIBS)
	

inject_frb : main_inject_frb.cpp dedisp_search.cpp
	g++ main_inject_frb.cpp mwa_fits.cpp dedisp_search.cpp sources.cpp -o  inject_frb $(OPT) -I /opt/caastro/bighorns/software/common -I../fitslib/ $(INCLUDES) -lcfitsio -L$(BIGHORNS)/lib -lfitslib $(CMN_LIBS)
	cp inject_frb $(BIGHORNS)/bin/
	
find_sources_fits : main_find_sources.cpp dedisp_search.cpp mwa_fits.cpp
	g++ main_find_sources.cpp mwa_fits.cpp dedisp_search.cpp sources.cpp -o find_sources_fits $(OPT) -I /opt/caastro/bighorns/software/common -I../fitslib/ $(INCLUDES) -lcfitsio -L$(BIGHORNS)/lib -lfitslib $(CMN_LIBS)
	cp find_sources_fits $(BIGHORNS)/bin/
	
clean :
	rm -f $(PROGNAME) *.o inject_frb calcstat create_dynaspec dynaspec_search test read_data_cube test_cubes
	
galaxy :
	rsync -avP *cpp galaxy:~/mwa_soft/frb_search/

read_data_cube : read_data_cube.cpp mwa_fits.cpp sources.cpp
	g++ read_data_cube.cpp mwa_fits.cpp sources.cpp -o read_data_cube $(OPT) -I /opt/caastro/bighorns/software/common -I../fitslib/ $(INCLUDES) -lcfitsio $(CMN_LIBS) -L$(BIGHORNS)/lib -lfitslib -lnova

test_cubes : read_data_cube_test.cpp
	g++ read_data_cube_test.cpp mwa_fits.cpp sources.cpp -o read_data_cube_test $(OPT) -I /opt/caastro/bighorns/software/common -I../fitslib/ $(INCLUDES) -lcfitsio -L$(BIGHORNS)/lib -lfitslib $(CMN_LIBS) 


	