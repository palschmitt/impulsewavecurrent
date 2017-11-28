TARGET = alglib_test

CFLAGS = -O3 -I ~/proj/alglib/src -Wall

LFLAGS =

CC = g++

SRC = test_i.cpp
SRC += ~/proj/alglib/src/ap.cpp
SRC += ~/proj/alglib/src/alglibinternal.cpp
SRC += ~/proj/alglib/src/alglibmisc.cpp
SRC += ~/proj/alglib/src/dataanalysis.cpp
SRC += ~/proj/alglib/src/diffequations.cpp
SRC += ~/proj/alglib/src/fasttransforms.cpp
SRC += ~/proj/alglib/src/integration.cpp
SRC += ~/proj/alglib/src/interpolation.cpp
SRC += ~/proj/alglib/src/linalg.cpp
SRC += ~/proj/alglib/src/optimization.cpp
SRC += ~/proj/alglib/src/solvers.cpp
SRC += ~/proj/alglib/src/specialfunctions.cpp
SRC += ~/proj/alglib/src/statistics.cpp

TMP = $(SRC:.cpp=.o)
OBJ = $(TMP:.c=.o)

$(TARGET): $(OBJ)
@echo Linking $@
@$(CC) $(LFLAGS) $(OBJ) -o $@
@echo Build Complete

.cpp.o: $<
@echo Compiling $<
@$(CC) -c $(CFLAGS) $< -o $@

clean:
@rm -f $(OBJ) $(TARGET)
@echo All object files and binaries removed
