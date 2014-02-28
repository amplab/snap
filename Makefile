ifndef CXXFLAGS
  CXXFLAGS = -O3 -Wno-format
endif

CXXFLAGS += -MMD -ISNAPLib -msse

LDFLAGS += -pthread

#LIBHDFS_HOME = ../hadoop-2.2.0-src/hadoop-hdfs-project/hadoop-hdfs/src/main/native/libhdfs
#JAVA_HOME = /usr/lib/jvm/java-7-oracle

ifdef LIBHDFS_HOME
  CXXFLAGS += -DSNAP_HDFS -I$(LIBHDFS_HOME)
  LDFLAGS += -L$(LIBHDFS_HOME) -L$(JAVA_HOME)/jre/lib/amd64/server -L$(JAVA_HOME)/jre/lib/amd64
  LIBS +=  -lhdfs -ljvm
endif

UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
  LIBS += -lrt -lz
endif

ifeq ($(UNAME), Darwin)
  LIBS += -lz
endif

CXX = g++

LIB_SRC = $(wildcard SNAPLib/*.cpp)
LIB_OBJ = $(patsubst %.cpp, %.o, $(LIB_SRC))

SNAP_SRC = $(wildcard apps/snap/*.cpp)
TEST_SRC = $(wildcard tests/*.cpp)
ROC_SRC = $(wildcard apps/ComputeROC/*.cpp)

SNAP_OBJ = $(patsubst %.cpp, %.o, $(SNAP_SRC))
TEST_OBJ = $(patsubst %.cpp, %.o, $(TEST_SRC))
ROC_OBJ = $(patsubst %.cpp, %.o, $(ROC_SRC))

ALL_OBJ = $(LIB_OBJ) $(SNAP_OBJ) $(TEST_OBJ)

DEPS = $(pathsubst %.o, %.d, $(ALL_OBJ))

EXES = snap unit_tests

default: $(EXES)

-include $(pathsubst %.o, %.d, $(ALL_OBJ))

$(OBJS): %.o : %.cpp
	$(CXX) -o $@ $(CXXFLAGS) -c $< 

snap: $(LIB_OBJ) $(SNAP_OBJ)
	$(CXX) -o $@ $(CXXFLAGS) -Itests $(LDFLAGS) $^ $(LIBS)

roc: $(LIB_OBJ) $(ROC_OBJ)
	$(CXX) -o $@ $(CXXFLAGS) -Itests $(LDFLAGS) $^ $(LIBS)

unit_tests: $(LIB_OBJ) $(TEST_OBJ)
	$(CXX) -o $@ $(CXXFLAGS) -Itests $(LDFLAGS) $^ $(LIBS)

clean:
	rm -f $(ALL_OBJ) $(DEPS) $(EXES)

.phony: clean default
