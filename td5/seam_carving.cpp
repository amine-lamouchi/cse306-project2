#define STB_IMAGE_IMPLEMENTATION
#include "../include/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../include/stb_image_write.h"

#include <iostream>
#include <vector>
#include <algorithm>

typedef struct {
    unsigned char r, g, b;
} RGB;

int main() {
    int width, height, bpp;
    RGB* rgb_image = (RGB*)stbi_load("image.jpg", &width, &height, &bpp, 3);
    
    if(rgb_image == nullptr) {
        std::cout << "Couldn't open image.";
        return 1;
    }

    int nun_seams_to_remove = 20;
    for(int iter = 0; iter < nun_seams_to_remove; ++iter) {
        std::vector<std::vector<int>> energy_map(height, std::vector<int>(width, 0));
        std::vector<std::vector<int>> cumulative_map(height, std::vector<int>(width, 0));
        std::vector<int> seam(height);

        for(int y = 1; y < height - 1; ++y) {
            for(int x = 1; x < width - 1; ++x) {
                RGB l = rgb_image[y * width + x - 1];
                RGB r = rgb_image[y * width + x + 1];
                RGB u = rgb_image[(y - 1) * width + x];
                RGB d = rgb_image[(y + 1) * width + x];
                energy_map[y][x] = abs(int(r.r) - int(l.r)) + abs(int(r.g) - int(l.g)) + abs(int(r.b) - int(l.b))
                                 + abs(int(d.r) - int(u.r)) + abs(int(d.g) - int(u.g)) + abs(int(d.b) - int(u.b));
            }
        }

        int max_energy = 0;
        for(const auto& row : energy_map) {
            max_energy = std::max(max_energy, *std::max_element(row.begin(), row.end()));
        }

        if (iter == 0) {
            std::vector<unsigned char> energy_greyscale(height * width);
            for(int y = 0; y < height; ++y) {
                for(int x = 0; x < width; ++x) {
                    unsigned char intensity = static_cast<unsigned char>(255 * energy_map[y][x] / max_energy);
                    energy_greyscale[y * width + x] = intensity;
                }
            }
            stbi_write_jpg("energy_map.jpg", width, height, 1, energy_greyscale.data(), 100);
        }

        for(int y = 1; y < height; ++y) {
            for(int x = 1; x < width - 1; ++x) {
                cumulative_map[y][x] = energy_map[y][x] + std::min({cumulative_map[y - 1][x - 1], cumulative_map[y - 1][x], cumulative_map[y - 1][x + 1]});
            }
        }

        seam[height - 1] = std::min_element(cumulative_map[height - 1].begin() + 1, cumulative_map[height - 1].end() - 1) - cumulative_map[height - 1].begin();
        for(int y = height - 2; y >= 0; --y) {
            int prev_x = seam[y + 1];
            seam[y] = std::min_element(cumulative_map[y].begin() + std::max(prev_x - 1, 1), cumulative_map[y].begin() + std::min(prev_x + 2, width - 1)) - cumulative_map[y].begin();
        }

        RGB* new_rgb_image = new RGB[height * (width - 1)];
        for(int y = 0; y < height; ++y) {
            for(int x = 0; x < seam[y]; ++x) {
                new_rgb_image[y * (width - 1) + x] = rgb_image[y * width + x];
            }
            for(int x = seam[y]; x < width - 1; ++x) {
                new_rgb_image[y * (width - 1) + x] = rgb_image[y * width + x + 1];
            }
        }
        delete[] rgb_image;
        rgb_image = new_rgb_image;
        --width;
    }

    stbi_write_jpg("output.jpg", width, height, 3, rgb_image, 100);
    delete[] rgb_image;
    return 0;
}
