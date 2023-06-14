CXX = g++-13

OMPFLAG = -fopenmp

OPTFLAG = -O3

CXXFLAGS = $(OPTFLAG) $(OMPFLAG) -std=c++17 -DSTB_IMAGE_WRITE_IMPLEMENTATION -g

TARGET = tutte

SRC = tutte_embedding.cpp

all: $(TARGET) tutte 

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	$(RM) $(TARGET)
	$(RM) flattened.obj

