CXX = g++-13

OMPFLAG = -fopenmp

OPTFLAG = -O3

CXXFLAGS = $(OPTFLAG) $(OMPFLAG) -std=c++17 -DSTB_IMAGE_WRITE_IMPLEMENTATION -g

TARGET = seam

SRC = seam_carving.cpp

all: $(TARGET) seam 

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	$(RM) $(TARGET)
	$(RM) energy_map.jpg
	$(RM) output.jpg

