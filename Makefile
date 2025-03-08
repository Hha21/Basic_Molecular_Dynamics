# COMPILER, FLAGS -Wall -Wextra
CXX = g++
CXXFLAGS = -std=c++17 -O2
LDFLAGS = -lboost_program_options

# DIRECTORIES
BUILD_DIR = build
SRC_DIR = .

SRCS = main.cpp ICScenario.cpp Solver.cpp Domain.cpp Particle.cpp Logger.cpp
OBJS = $(SRCS:%.cpp=$(BUILD_DIR)/%.o)	 # .cpp -> .o in ./build
TARGET = md  # md.out

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(TARGET)

run: $(TARGET)
	./$(TARGET)