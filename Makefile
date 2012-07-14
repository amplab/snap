ifndef CXXFLAGS
  CXXFLAGS = -O3
endif

CXXFLAGS += -MMD -ISNAPLib

LDFLAGS += -pthread

UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
  LIBS += -lrt
endif

CXX = g++

LIB_SRC = $(wildcard SNAPLib/*.cpp)
LIB_OBJ = $(patsubst %.cpp, %.o, $(LIB_SRC))

SNAP_SRC = $(wildcard apps/snap/*.cpp)
TEST_SRC = $(wildcard tests/*.cpp)

SNAP_OBJ = $(patsubst %.cpp, %.o, $(SNAP_SRC))
TEST_OBJ = $(patsubst %.cpp, %.o, $(TEST_SRC))

ALL_OBJ = $(LIB_OBJ) $(SNAP_OBJ) $(TEST_OBJ)

DEPS = $(pathsubst %.o, %.d, $(ALL_OBJ))

EXES = snap unit_tests

default: $(EXES)

-include $(pathsubst %.o, %.d, $(ALL_OBJ))

$(OBJS): %.o : %.cpp
	$(CXX) -o $@ $(CXXFLAGS) -c $< 

snap: $(LIB_OBJ) $(SNAP_OBJ)
	$(CXX) -o $@ $(CXXFLAGS) -Itests $(LDFLAGS) $^ $(LIBS)

unit_tests: $(LIB_OBJ) $(TEST_OBJ)
	$(CXX) -o $@ $(CXXFLAGS) -Itests $(LDFLAGS) $^ $(LIBS)

clean:
	rm -f $(ALL_OBJ) $(DEPS) $(EXES)

.phony: clean default
