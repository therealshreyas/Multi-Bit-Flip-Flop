SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

LEMONDIR	  = /home/tool/lemon/lemon-1.3.1/include/ 
LEMONLIBDIR	  = /home/tool/lemon/
ORDIR = /home/tool/ortools/install/CentOS7-gcc9/include
ORLIBDIR = /home/tool/ortools/install/CentOS7-gcc9/lib64

CCPATH   = g++
CXXOPTS  = -m64 -std=c++17 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD -Wno-ctor-dtor-privacy -fopenmp -O3
CCFLAG = $(CXXOPTS)  -I$(LEMONDIR) -I$(ORDIR)

SYSLIBS  = -ldl
DEBUG = -g -gstrict-dwarf -gdwarf-2

CCLNFLAG = $(SYSLIBS) -L$(LEMONLIBDIR) -lm -L$(ORLIBDIR) -lortools -fopenmp -std=c++17 -static-libstdc++
 

# Link the executable
mbff: main.cpp mbff.cpp 
	$(CCPATH) $(DEBUG) $(CXXOPTS) $(CCFLAG) main.cpp mbff.cpp -o mbff  \
        $(CCLNFLAG) 

clean: 
	@/bin/rm  -rf mbff
