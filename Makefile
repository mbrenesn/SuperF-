default : RLevel.x

CC_FILES := $(wildcard src/*/*.cc)
OBJ_FILES := $(addprefix obj/,$(notdir $(CC_FILES:.cc=.o)))

CXXLINKER = icpc
CXX = icpc
CXXFLAGS = -g -O2 -DMKL_ILP64

RLevel.x : $(OBJ_FILES)
	$(CXXLINKER) $^ -o $@ -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

obj/%.o : src/*/%.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $< -I${MKLROOT}/include

clean :
	rm obj/*.o *.x
