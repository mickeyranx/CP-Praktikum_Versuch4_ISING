# Compiler
CXX = g++
CXXFLAGS = -std=c++20 -O2 -Wall -Wextra

# Executable name
TARGET = run_ising

# Source files
SRCS = main.cpp Simulation.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Default target
all: $(TARGET)

# Link
$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET)

# Compile
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(OBJS) $(TARGET)

rebuild: clean all