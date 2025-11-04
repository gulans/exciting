# This Makefile is for compiling the C++ parts of the Faddeeva library.
# It creates a static library archive file named "libfaddeeva.a".

# Get compiler definitions from the main exciting configuration
include ../../build/make.inc

# List of the C++ source files
CXX_SOURCES = Faddeeva.cc faddeeva_wrapper.cpp

# Pattern to turn source files (.cc, .cpp) into object files (.o)
CXX_OBJECTS = $(CXX_SOURCES:.cc=.o)
CXX_OBJECTS := $(CXX_OBJECTS:.cpp=.o)

# The final library name
TARGET_LIB = libfaddeeva.a

# Default command: build the library
all: $(TARGET_LIB)

# Rule to build the library archive
$(TARGET_LIB): $(CXX_OBJECTS)
	ar rcs $(TARGET_LIB) $(CXX_OBJECTS)
	@echo "Faddeeva C++ library compiled successfully: $(TARGET_LIB)"

# Rule to compile a .cc file into a .o file
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Rule to compile a .cpp file into a .o file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Command to clean up compiled files
clean:
	rm -f *.o $(TARGET_LIB)