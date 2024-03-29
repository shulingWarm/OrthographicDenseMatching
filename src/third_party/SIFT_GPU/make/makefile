#################################################################
#    SiftGPU congiruation: CG, CUDA, SSE, TIMING  
#################################################################
# enable cg-based SiftGPU?
siftgpu_enable_cg = 0
#enable siftgpu serve
siftgpu_enable_server = 1

#------------------------------------------------------------------------------------------------
# enable CUDA-based SiftGPU?
siftgpu_enable_cuda = 1
CUDA_INSTALL_PATH = /usr/local/cuda-8.0
#change  additional  settings, like SM version here if it is not 1.0 (eg. -arch sm_13 for GTX280)
#siftgpu_cuda_options = -Xopencc -OPT:unroll_size=200000
siftgpu_cuda_options = -arch sm_52
#--------------------------------------------------------------------------------------------------
# enable SSE optimization for GL-based implementations
siftgpu_enable_sse = 1
siftgpu_sse_options = -march=native -mfpmath=sse
#-------------------------------------------------------------------------------------------------
# timing function, 1 for clock, 0 for gettimeofday
siftgpu_timing_use_clock = 0
#--------------------------------------------------------------------------------------------------
# openGL context creation.  1 for glut, 0 for xlib
siftgpu_prefer_glut = 0
#------------------------------------------------------------------------------------------------
#where SimpleSIFT is using runtime loading of libsiftgpu.so or static linking of libsiftgpu.a
simplesift_runtime_load = 1
#################################################################




# detect OS
OSUPPER = $(shell uname -s 2>/dev/null | tr [:lower:] [:upper:])
OSLOWER = $(shell uname -s 2>/dev/null | tr [:upper:] [:lower:])
DARWIN = $(strip $(findstring DARWIN, $(OSUPPER)))


SHELL = /bin/sh
INC_DIR = ../include
BIN_DIR = ./bin
SRC_SIFTGPU = ../SiftGPU
SRC_DRIVER = ../TestWin
SRC_SERVER = ../ServerSiftGPU
CC = g++
CFLAGS =  -I$(INC_DIR) -fPIC -L/usr/lib -L./bin -Wall -Wno-deprecated -pthread  

ifneq ($(siftgpu_enable_sse), 0)
	CFLAGS += $(siftgpu_sse_options)
endif

ifneq ($(siftgpu_timing_use_clock), 0)
	CFLAGS += -DTIMING_BY_CLOCK
endif

ifneq ($(siftgpu_prefer_glut), 0)
	CFLAGS += -DWINDOW_PREFER_GLUT
endif


ODIR_SIFTGPU = siftgpu


# external header files
_HEADER_EXTERNAL = GL/glew.h GL/glut.h IL/il.h  
# siftgpu header files
_HEADER_SIFTGPU = FrameBufferObject.h GlobalUtil.h GLTexImage.h ProgramGPU.h ShaderMan.h ProgramGLSL.h SiftGPU.h SiftPyramid.h SiftMatch.h PyramidGL.h LiteWindow.h
# siftgpu library header files for drivers
_HEADER_SIFTGPU_LIB = SiftGPU.h  

ifneq ($(DARWIN),) 
#librarys for SiftGPU
LIBS_SIFTGPU = -lGLEW -framework GLUT -lIL -framework OpenGL 
CFLAGS += -fPIC -L/opt/local/lib -L/usr/local/lib 
else
#librarys for SiftGPU
LIBS_SIFTGPU = -lGLEW -lglut -lIL -lGL
endif
 
#Obj files for SiftGPU
_OBJ_SIFTGPU = FrameBufferObject.o GlobalUtil.o GLTexImage.o ProgramGLSL.o ProgramGPU.o ShaderMan.o SiftGPU.o SiftPyramid.o PyramidGL.o SiftMatch.o

#add cg options 
ifneq ($(siftgpu_enable_cg), 0)
	CFLAGS += -DCG_SIFTGPU_ENABLED
	_OBJ_SIFTGPU += ProgramCG.o 
	_HEADER_SIFTGPU += ProgramCG.h 
	ifneq ($(DARWIN),) 
		LIBS_SIFTGPU += -framework CG
	else
		LIBS_SIFTGPU += -lCg -lCgGL 
	endif
endif

#add cuda options
ifneq ($(siftgpu_enable_cuda), 0)
	ifdef CUDA_BIN_PATH
		NVCC = $(CUDA_BIN_PATH)/nvcc
	else
		NVCC = $(CUDA_INSTALL_PATH)/bin/nvcc
	endif

	ifndef CUDA_INC_PATH
		CUDA_INC_PATH = $(CUDA_INSTALL_PATH)/include
	endif

	ifndef CUDA_LIB_PATH
		CUDA_LIB_PATH = $(CUDA_INSTALL_PATH)/lib64
	endif

	CFLAGS += -DCUDA_SIFTGPU_ENABLED -I$(CUDA_INC_PATH) -L$(CUDA_LIB_PATH)
	LIBS_SIFTGPU += -lcudart -lX11
	_OBJ_SIFTGPU += CuTexImage.o PyramidCU.o SiftMatchCU.o
	_HEADER_SIFTGPU += CuTexImage.h ProgramCU.h PyramidCU.h
endif
 
all: makepath siftgpu server  driver 
 

#the dependencies of SiftGPU library 
DEPS_SIFTGPU = $(patsubst %, $(SRC_SIFTGPU)/%, $(_HEADER_SIFTGPU))


#rules for the rest of the object files
$(ODIR_SIFTGPU)/%.o: $(SRC_SIFTGPU)/%.cpp $(DEPS_SIFTGPU) 
	$(CC) -o $@ $< $(CFLAGS) -c 


ifneq ($(siftgpu_enable_cuda), 0)
NVCC_FLAGS = -I$(INC_DIR) -I$(CUDA_INC_PATH) -DCUDA_SIFTGPU_ENABLED -O2
ifdef siftgpu_cuda_options
	NVCC_FLAGS += $(siftgpu_cuda_options)
endif
#build rule for CUDA 
$(ODIR_SIFTGPU)/ProgramCU.o: $(SRC_SIFTGPU)/ProgramCU.cu $(DEPS_SIFTGPU)
	$(NVCC) $(NVCC_FLAGS) -o $@ $< --compiler-options "-fPIC" -c
_OBJ_SIFTGPU += ProgramCU.o
endif


ifneq ($(siftgpu_enable_server), 0)
$(ODIR_SIFTGPU)/ServerSiftGPU.o: $(SRC_SERVER)/ServerSiftGPU.cpp $(DEPS_SIFTGPU)
	$(CC) -o $@ $< $(CFLAGS) -DSERVER_SIFTGPU_ENABLED -c
_OBJ_SIFTGPU += ServerSiftGPU.o
endif

OBJ_SIFTGPU = $(patsubst %,$(ODIR_SIFTGPU)/%,$(_OBJ_SIFTGPU))
LIBS_DRIVER = $(BIN_DIR)/libsiftgpu.a $(LIBS_SIFTGPU) 
SRC_TESTWIN = $(SRC_DRIVER)/TestWinGlut.cpp $(SRC_DRIVER)/BasicTestWin.cpp  
DEP_TESTWIN = $(SRC_DRIVER)/TestWinGlut.h $(SRC_DRIVER)/BasicTestwin.h $(SRC_DRIVER)/GLTransform.h 



ifneq ($(simplesift_runtime_load), 0)
LIBS_SIMPLESIFT = -ldl -DSIFTGPU_DLL_RUNTIME
else
LIBS_SIMPLESIFT = $(LIBS_DRIVER) -DSIFTGPU_STATIC
endif

siftgpu: makepath $(OBJ_SIFTGPU)
	ar rcs $(BIN_DIR)/libsiftgpu.a $(OBJ_SIFTGPU) 
	$(CC) -shared -o $(BIN_DIR)/libsiftgpu.so $(OBJ_SIFTGPU) $(LIBS_SIFTGPU) $(CFLAGS)

 
driver: makepath 
	$(CC) -o $(BIN_DIR)/TestWinGlut $(SRC_TESTWIN) $(LIBS_DRIVER) $(CFLAGS)
	$(CC) -o $(BIN_DIR)/SimpleSIFT $(SRC_DRIVER)/SimpleSIFT.cpp $(LIBS_SIMPLESIFT) $(CFLAGS) 
	$(CC) -o $(BIN_DIR)/speed $(SRC_DRIVER)/speed.cpp $(LIBS_DRIVER) $(CFLAGS) 
	$(CC) -o $(BIN_DIR)/MultiThreadSIFT $(SRC_DRIVER)/MultiThreadSIFT.cpp $(LIBS_DRIVER) $(CFLAGS)  -pthread
	
ifneq ($(siftgpu_enable_server), 0)
server: makepath
	$(CC) -o $(BIN_DIR)/server_siftgpu $(SRC_SERVER)/server.cpp $(LIBS_DRIVER) $(CFLAGS)
else
server: 

endif	
	
makepath:
	mkdir -p $(ODIR_SIFTGPU)
	mkdir -p $(BIN_DIR) 
	ln ../data/ -d -s -f
 
clean:
	rm -f $(ODIR_SIFTGPU)/*.o
	rm -f $(BIN_DIR)/libsiftgpu.a
	rm -f $(BIN_DIR)/libsiftgpu.so
	rm -f $(BIN_DIR)/TestWinGlut
	rm -f $(BIN_DIR)/SimpleSIFT
	rm -f $(BIN_DIR)/speed
	rm -f $(BIN_DIR)/server_siftgpu
	rm -f $(BIN_DIR)/MultiThreadSIFT
	rm -f ProgramCU.linkinfo

