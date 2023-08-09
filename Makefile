SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

LEMONDIR	  = /home/tool/lemon/lemon-1.3.1/include/ 
CPLEXDIR      = /home/tool/ILOG/CPLEX_Studio1263/cplex/include
CONCERTDIR    = /home/tool/ILOG/CPLEX_Studio1263/concert/include
LEMONLIBDIR	  = /home/tool/lemon/
CPLEXLIBDIR   = $(CPLEXDIR)/../lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/../lib/$(SYSTEM)/$(LIBFORMAT)



CCPATH   = g++
CXXOPTS  = -m64 -O3 -std=c++14 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -Wno-ctor-dtor-privacy -fopenmp
CCFLAG = $(CXXOPTS)  -I$(CPLEXDIR) -I$(CONCERTDIR) -I$(LEMONDIR)

SYSLIBS  = -ldl
DEBUG = -g -gstrict-dwarf -gdwarf-2

CCLNFLAG = $(SYSLIBS) -L${CPLEXLIBDIR} -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -L$(LEMONLIBDIR) -lm -pthread -fopenmp -std=c++0x -static-libstdc++
 


# Link the executable
mbff: mbff.cpp 
	$(CCPATH) $(DEBUG) $(CXXOPTS) $(CCFLAG) mbff.cpp -o mbff \
        $(CCLNFLAG) 


rng: RNG.cpp
	$(CCPATH) $(CXXOPTS) RNG.cpp -o rng


clean: 
	@/bin/rm  -rf rng
	@/bin/rm  -rf mbff
