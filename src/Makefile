CXX=g++
LINKER=g++
ROOTCONFIG := root-config
UCESB_DIR := /u/land/fake_cvmfs/10/cvmfs_fairroot_v18.8.0_fs_nov22p1_extra/ucesb
CFLAGS := $(shell $(ROOTCONFIG) --cflags)
CFLAGS += -I$(UCESB_DIR)/hbook
CFLAGS += --std=c++17 -g -O0 -fexceptions
LDFLAGS := $(shell $(ROOTCONFIG) --libs)
LDFLAGS += -lRHTTP
LDFLAGS += -L$(UCESB_DIR)/hbook -lext_data_clnt
LDFLAGS += -g
LDFLAGS += -L$(ROOT_LIBRARY_DIR) 

EXEC=event_display
OBJ=event_display.o

all: $(EXEC)

$(EXEC): $(OBJ)
	${LINKER} -o $@ $^ ${LDFLAGS}
	echo " COMP $@"

%.o : %.cpp
	$(MAKEDEPEND)
	${CXX} ${CFLAGS} $(DIR_INC) -c $< -o $@
	echo "  CXX $@"

clean:
	rm -f *.o $(EXEC)
