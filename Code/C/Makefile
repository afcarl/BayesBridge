# Include for Matrix and GSL.
INC = -I$(HOME)/RV-Project/Code/C_Examples/MyLib/Matrix/ \
	-I$(HOME)/RV-Project/Code/C_Examples/MyLib/RNG/ \
	-I$(HOME)/Code/include/ \
	-I$\./

GLIB = -L$(HOME)/Code/lib

RINC = $(shell R CMD config --cppflags)
RLNK = $(shell R CMD config --ldflags)

# OPT = -g -Wall
OPT = -g -Wall -O2 -DNDEBUG
# OPT = -03

CC = g++

# rbridge : Bridge.h Bridge.cpp
all :
	$(CC) Bridge.cpp $(INC) $(RINC) -DUSE_R $(OPT) -shared -fPIC -lblas -llapack -o libbridge.so $(RLNK)

test1 : test.c libHmcSampler.so
	$(CC) test.c $(INC) $(OPT) $(GLIB) -L. -o test1 -lgsl -lHmcSampler

test2 : test2.cpp BridgeRegression.h gbridge2
	$(CC) test2.cpp $(INC) $(OPT) $(GLIB) -L. -lBridge -lHmcSampler -lgsl -lblas -llapack  -o test2

gbridge2 : BridgeWrapper.h BridgeWrapper.cpp libHmcSampler.so
	$(CC) BridgeWrapper.cpp $(INC) $(GLIB) $(OPT) -L. -shared -fPIC -lgsl -lblas -llapack -lHmcSampler -o libBridge.so

GBridge : BridgeWrapper.h BridgeWrapper.cpp
	$(CC) BridgeWrapper.cpp $(INC) $(GLIB) $(OPT) -shared -fPIC -lgsl -lblas -llapack -o Bridge.so

RBridge : BridgeWrapper.h BridgeWrapper.cpp libHmcSampler.so
	$(CC) BridgeWrapper.cpp $(INC) $(RINC) $(RLNK) $(OPT) -L. -shared -fPIC -DUSE_R \
	-lblas -llapack -lHmcSampler -o Bridge.so

libbridge.dylib : Bridge.h Bridge.cpp
	$(CC) Bridge.cpp $(INC) $(GLIB) $(OPT) -dynamiclib -fPIC -lgsl -lblas -llapack -o libbridge.dylib

test : test.cpp BridgeRegression.h
	$(CC) test.cpp $(INC) $(GLIB) $(OPT) -lgsl -lblas -llapack -o test

libHmcSampler.so : HmcSampler.cpp HmcSampler.h
	$(CC) HmcSampler.cpp $(INC) $(OPT) -shared -fPIC -o libHmcSampler.so

clean :
	rm *.so

################################################################################