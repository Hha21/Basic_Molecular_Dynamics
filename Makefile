# COMPILER, FLAGS -Wall -Wextra
CXX = g++
NVCC = nvcc
CXXFLAGS = -Ofast -march=native -funroll-loops -std=c++17
LDFLAGS = -lboost_program_options
NVCCFLAGS = -O3 -std=c++17

# DIRECTORIES
SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin
DOC_DIR = doc

SRCS_COMMON = $(wildcard $(SRC_DIR)/*.cpp)
SRCS_PAR = $(filter-out $(SRC_DIR)/SolverSerial.cpp $(SRC_DIR)/SolverSerial.h $(SRC_DIR)/Unit_Tests.cpp $(SRC_DIR)/Unit_Tests_Par.cpp $(SRC_DIR)/main.cpp $(SRC_DIR)/SolverSerialCUDA.cpp $(SRC_DIR)/computeLJ.cu, $(SRCS_COMMON))
SRCS_SERIAL = $(filter-out $(SRC_DIR)/Solver.cpp $(SRC_DIR)/Solver.h $(SRC_DIR)/Unit_Tests.cpp $(SRC_DIR)/Unit_Tests_Par.cpp $(SRC_DIR)/main.cpp $(SRC_DIR)/SolverSerialCUDA.cpp $(SRC_DIR)/computeLJ.cu, $(SRCS_COMMON))
SRCS_CUDA = $(SRC_DIR)/SolverSerialCUDA.cu 

OBJS_PAR = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS_PAR))
OBJS_SERIAL = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS_SERIAL))
OBJS_CUDA = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(filter %.cpp, $(SRCS_CUDA))) \
            $(patsubst $(SRC_DIR)/%.cu, $(BUILD_DIR)/%.o, $(filter %.cu, $(SRCS_CUDA)))


#UNIT TESTS
TEST_SRC_PAR = $(SRC_DIR)/Unit_Tests_Par.cpp
TEST_SRC_SERIAL = $(SRC_DIR)/Unit_Tests.cpp
TEST_PAR = $(BIN_DIR)/unittests_par
TEST_SERIAL = $(BIN_DIR)/unittests_serial

#TARGET EXECUTABLES
TARGET_PAR = $(BIN_DIR)/mdpar
TARGET_SERIAL = $(BIN_DIR)/md
TARGET_CUDA = $(BIN_DIR)/mdcuda

# DEFAULT
default: $(TARGET_PAR) $(TARGET_SERIAL) $(TARGET_CUDA)

mdpar: $(TARGET_PAR)
md: $(TARGET_SERIAL)
mdcuda: $(TARGET_CUDA)

unittests: $(TEST_PAR) $(TEST_SERIAL)
	$(TEST_PAR)
	$(TEST_SERIAL)

$(TARGET_PAR): $(OBJS_PAR) $(BUILD_DIR)/main.o
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -fopenmp $^ -o $@ $(LDFLAGS)

$(TARGET_SERIAL): $(OBJS_SERIAL) $(BUILD_DIR)/main.o
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(TARGET_CUDA): $(OBJS_CUDA) $(BUILD_DIR)/main.o
	mkdir -p $(BIN_DIR)
	$(NVCC) $(NVCCFLAGS) $^ -o $@ $(LDFLAGS) 

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cu
	mkdir -p $(BUILD_DIR)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(TEST_PAR): $(TEST_SRC_PAR) $(OBJS_PAR)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) -fopenmp $^ -o $@ $(LDFLAGS)

$(TEST_SERIAL): $(TEST_SRC_SERIAL) $(OBJS_SERIAL)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

doc:
	doxygen Doxyfile

clean-doc:
	rm -rf docs/html docs/latex

run-par: $(TARGET_PAR)
	./$(TARGET_PAR)

run-serial: $(TARGET_SERIAL)
	./$(TARGET_SERIAL)

run-cuda: $(TARGET_CUDA)
	./$(TARGET_CUDA)