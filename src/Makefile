# Compiler and flags
CXX = g++
# CXXFLAGS = -std=c++17
#LDFLAGS = -L/path/to/libraries -lfftw3  # Adjust path and library flags 
LDFLAGS = -L/C:/fftw-3.3.5-dll64 -lfftw3-3
TARGET = run

# Define the paths to FFTW
FFTW_DIR = C:/fftw-3.3.5-dll64

# Define the compilation flags
CXXFLAGS = -std=c++17 -fopenmp -o -I$(FFTW_DIR)

# CXXFLAGS = -std=c++17 -I$(FFTW_DIR)

# Define the linker flags
LDFLAGS = -L$(FFTW_DIR) -lfftw3-3 -lfftw3f-3 -lfftw3l-3 -mconsole

# Source and object files
# SRCS = main.cpp electrostatics.cpp io.cpp
SRCS = main.cpp electrostatics.cpp io.cpp
OBJS = $(SRCS:.cpp=.o)

# Default target
all: $(TARGET)

# Build target
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

# Compile source files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Debug build
debug: CXXFLAGS += -g
debug: $(TARGET)
	gdb $(TARGET)

# Clean up build files
clean:
	del $(TARGET) $(OBJS)
	rmdir /s /q output energy rho relative deposition

# Create directories
prepare:
	mkdir output
	mkdir energy
	mkdir rho
	mkdir relative
	mkdir deposition