CXX = g++-13

OMPFLAG = -fopenmp

OPTFLAG = -O3

CXXFLAGS = -I./src -I./include -I./lib $(OPTFLAG) $(OMPFLAG) -std=c++17 -DSTB_IMAGE_WRITE_IMPLEMENTATION

TARGET = main

SRC = main.cpp lib/lbfgs.c

all: $(TARGET) main 

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

create_video: $(TARGET)
	ffmpeg -r 10 -f image2 -i ./animation_frames/animation%d.png -vcodec libx264 -crf 25  -pix_fmt yuv420p animation.mp4

clean:
	$(RM) $(TARGET)
	$(RM) animation*.png
	$(RM) animation.mp4
	$(RM) voronoi.svg
