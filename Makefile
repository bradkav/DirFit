NESTLIBDIR = /Users/bradkav/Code/MultiNest_v3.10/
DECODMDIR = /Users/bradkav/Code/DeCoDM/

LIBS =  -L/opt/local/lib -I/opt/local/include -L/opt/local/lib -L$(NESTLIBDIR) -L$(DECODMDIR)  -lm -lgsl -lgslcblas -lnest3 -llapack -lstdc++ -ldecodm 

FC = gfortran

OBJFILES = DirFit.o

all: DirFit 

%.o: %.cc
	$(CXX) $(CFLAGS) -I/Users/bradkav/Code/DeCoDM/include $(LIBS) -c $*.cc 
 
DirFit: $(OBJFILES)
	$(FC) $(FFLAGS) -o DirFit  $(OBJFILES) $(LIBS)

clean:
	rm -f *.o *.mod DirFit
