## PARAMETERS

# Include for Matrix and GSL.
UINC = -I./include -I$(HOME)/Code/Matrix -I$(HOME)/Code/RNG

GINC = -I$(HOME)/Code/include
GLIB = -L$(HOME)/Code/lib

RINC = $(shell R CMD config --cppflags)
RLNK = $(shell R CMD config --ldflags)

# make USE_R=-DUSE_R whatever
# USE_R = -DUSE_R
# USE_R = 

# ## For BLAS and LAPACK -- I need to figure out how to find these libraries for any system.
# MKLLNK = -undefined dynamic_lookup -single_module -multiply_defined suppress \
#        -L/Library/Frameworks/R.framework/Resources/lib/x86_64 -lRblas -lgfortran \
#        -L/Library/Frameworks/R.framework/Resources/lib/x86_64 -lRlapack

VANLNK = -lblas -llapack

MATLNK = $(VANLNK)

ifndef USE_R
INC = $(UINC) $(GINC)
LNK = $(GLIB) $(MATLNK) -lgsl
DEP = GRNG.o 
else
INC = $(UINC) $(RINC)
LNK = $(MATLNK) $(RLNK) 
DEP = RRNG.o 
endif

OPTS = -O2 $(USE_R) -g -DNDEBUG
OPTS = -g -Wall -std=c++98 -pedantic $(USE_R) #-Wshadow
## OPTS = -g -Wall -pedantic $(USE_R)
## Note: pedantic kills eigen.

sobj = BayesLogit.so
obj  = BridgeWrapper.o BridgeRegression.o retstable.o Matrix.o MatrixFrame.o RNG.o $(DEP)
obj2 = retstable.o Matrix.o MatrixFrame.o RNG.o $(DEP)

test_bridge : test_bridge.cpp $(obj)
	g++ -o test_bridge test_bridge.cpp $(OPTS) $(INC) $(obj) $(LNK) $(GLIB) -lgsl

## Link gsl last so that we will use R libraries when compiling for R.
BayesBridge.so : $(obj) 
	g++ -o BayesBridge.so $(obj) $(OPTS) -shared -fPIC $(LNK) $(GLIB) -lgsl

BridgeWrapper.o : BridgeWrapper.h BridgeWrapper.cpp BridgeRegression.o retstable.o Matrix.o
	g++ $(OPTS) $(INC) -fPIC -o BridgeWrapper.o -c BridgeWrapper.cpp 

BridgeRegression.o : BridgeRegression.h BridgeRegression.cpp
	g++ $(OPTS) $(INC) -fPIC -o BridgeRegression.o -c BridgeRegression.cpp

retstable.o : retstable.cpp retstable.h RNG.o
	g++ $(OPTS) $(INC) -fPIC -o retstable.o -c retstable.cpp

Matrix.o : MatrixFrame.o 
	g++ $(OPTS) $(INC) -fPIC -c $(HOME)/Code/Matrix/Matrix.cpp -o Matrix.o

MatrixFrame.o : 
	g++ $(OPTS) $(INC) -fPIC -c $(HOME)/Code/Matrix/MatrixFrame.cpp -o MatrixFrame.o

RNG.o : $(DEP)
	g++ $(OPTS) $(INC) -fPIC -c $(HOME)/Code/RNG/RNG.cpp -o RNG.o 

GRNG.o: 
	g++ $(INC) $(OPTS) -fPIC -c $(HOME)/Code/RNG/GRNG.cpp -o GRNG.o

RRNG.o: 
	g++ $(INC) $(OPTS) -fPIC -c $(HOME)/Code/RNG/RRNG.cpp -o RRNG.o

clean :
	rm *.o

clean.so :
	rm *.so

################################################################################