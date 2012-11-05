# Include for Matrix and GSL.
# INC = -I$(HOME)/RV-Project/Code/C_Examples/MyLib/Matrix/ \
# 	-I$(HOME)/RV-Project/Code/C_Examples/MyLib/RNG/ \
# 	-I$(HOME)/Code/include/ \
# 	-I$\./

INC = -I$(HOME)/Code/Matrix/ \
	-I$(HOME)/Code/RNG/ \
	-I$(HOME)/Code/include/ \
	-I$\./

GLIB = -L$(HOME)/Code/lib

RINC = $(shell R CMD config --cppflags)
RLNK = $(shell R CMD config --ldflags)

# OPT = -g -Wall
OPT = -g -Wall -O2 -DNDEBUG -Wunused -Wuninitialized
# OPT = -03

CC = g++

## -DYA_BLAS -DYA_LAPACK -DYA_BLASMULT -lblas -llapack -framework vecLib

# rbridge : Bridge.h Bridge.cpp
all :
	$(CC) Bridge.cpp $(INC) $(RINC) -DUSE_R $(OPT) -shared -fPIC -lblas -llapack -o libbridge.so $(RLNK)

test1 : test.c libHmcSampler.so
	$(CC) test.c $(INC) $(OPT) $(GLIB) -L. -o test1 -lgsl -lHmcSampler

test2 : test2.cpp BridgeRegression.h libBridge.so
	$(CC) test2.cpp $(INC) $(OPT) $(GLIB) -L. -lBridge \
	-lHmcSampler -lgsl -lblas -llapack  -o test2

# libBridge.so : BridgeWrapper.h BridgeWrapper.cpp libHmcSampler.so
# 	$(CC) BridgeWrapper.cpp $(INC) $(GLIB) $(OPT) -L. -shared -fPIC \
# 	-lblas -llapack -lHmcSampler -lgsl -o libBridge.so

# libBridge.so : BridgeWrapper.h BridgeWrapper.cpp libHmcSampler.so
# 	$(CC) BridgeWrapper.cpp $(INC) $(GLIB) $(OPT) -L. -shared -fPIC \
# 	-undefined dynamic_lookup -single_module -multiply_defined suppress \
# 	-DYA_BLAS -DYA_LAPACK -DYA_BLASMULT -latlas -lblas -llapack -framework vecLib	\
# 	-lblas -llapack -lHmcSampler -lgsl -o libBridge.so

libBridge.so : BridgeWrapper.h BridgeWrapper.cpp BridgeRegression.h libHmcSampler.so
	$(CC) BridgeWrapper.cpp $(INC) $(GLIB) $(OPT) -L. -shared -fPIC \
	-undefined dynamic_lookup -single_module -multiply_defined suppress \
	-L/usr/local/lib -framework R -Wl,-framework -Wl,CoreFoundation \
	-L/Library/Frameworks/R.framework/Resources/lib/x86_64 -lRblas -lgfortran \
	-L/Library/Frameworks/R.framework/Resources/lib/x86_64 -lRlapack -F/Library/Frameworks/R.framework/.. \
	-lgsl -lHmcSampler -o libBridge.so 
	# -lblas -llapack

GBridge2 : BridgeWrapper.h BridgeWrapper.cpp libHmcSampler.so
	$(CC) BridgeWrapper.cpp $(INC) $(GLIB) $(OPT) -L. -shared -fPIC -lgsl -lblas -llapack -lHmcSampler -o Bridge.so

GBridge : BridgeWrapper.h BridgeWrapper.cpp libHmcSampler.so
	$(CC) BridgeWrapper.cpp $(INC) $(GLIB) $(OPT) -L. -shared -fPIC \
	-undefined dynamic_lookup -single_module -multiply_defined suppress -L/usr/local/lib -L/Library/Frameworks/R.framework/Resources/lib/x86_64 -lRblas -lgfortran -L/Library/Frameworks/R.framework/Resources/lib/x86_64 -lRlapack -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation \
	-lgsl  -lHmcSampler -o Bridge.so 
	# -lblas -llapack	

RBridge : BridgeWrapper.h BridgeWrapper.cpp libHmcSampler.so
	$(CC) BridgeWrapper.cpp $(INC) $(RINC) $(RLNK) $(OPT) -L. -shared -fPIC -DUSE_R \
	-undefined dynamic_lookup -single_module -multiply_defined suppress -L/usr/local/lib -L/Library/Frameworks/R.framework/Resources/lib/x86_64 -lRblas -lgfortran -L/Library/Frameworks/R.framework/Resources/lib/x86_64 -lRlapack -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation \
	 -lHmcSampler -o Bridge.so
	#-lblas -llapack

RBridge2 : BridgeWrapper.h BridgeWrapper.cpp libHmcSampler.so
	$(CC) BridgeWrapper.cpp $(INC) $(RINC) $(RLNK) $(OPT) -L. -shared -fPIC -DUSE_R \
	 -lHmcSampler -o Bridge.so -lblas -llapack -latlas -Bdirect,--hash-style=both,-Wl,-O1	

libbridge.dylib : Bridge.h Bridge.cpp
	$(CC) Bridge.cpp $(INC) $(GLIB) $(OPT) -dynamiclib -fPIC -lgsl -lblas -llapack -o libbridge.dylib

test : test.cpp BridgeRegression.h
	$(CC) test.cpp $(INC) $(GLIB) $(OPT) -lgsl -lblas -llapack -o test

libHmcSampler.so : HmcSampler.cpp HmcSampler.h
	$(CC) HmcSampler.cpp $(INC) $(OPT) -shared -fPIC -o libHmcSampler.so

clean :
	rm *.so

################################################################################