# file Makefile
CXX ?= g++

# path #
SRC_PATH = src
BUILD_PATH = build
BIN_PATH = $(BUILD_PATH)/bin

# extensions #
SRC_EXT = cpp

COMPILE_FLAGS = -std=c++17 -Wall -Wextra -g -O0
INCLUDES = -I include/ 

SRCFILES= $(wildcard $(SRC_PATH)/*.$(SRC_EXT))
OBJFILES= $(patsubst %.$(SRC_EXT), %.o, $(SRCFILES))
PROGFILES= $(patsubst %.$(SRC_EXT), %, $(SRCFILES))

%.exe : $(SRCFILES)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MP -MMD -c $< -o $@

## eof Makefile

