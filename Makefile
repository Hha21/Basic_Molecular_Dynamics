# COMPILER, FLAGS -Wall -Wextra
CXX = g++
CXXFLAGS = -O3 -std=c++17 
LDFLAGS = -lboost_program_options

# DIRECTORIES
BUILD_DIR = build
SRC_DIR = src
BIN_DIR = bin
DOC_DIR = doc

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))	 # .cpp -> .o in ./build
TARGET = $(BIN_DIR)/md

all: $(TARGET)

$(TARGET): $(OBJS)
	mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@ -lboost_program_options

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

doc:
	doxygen Doxyfile

clean-doc:
	rm -rf docs/html docs/latex

run: $(TARGET)
	./$(TARGET)