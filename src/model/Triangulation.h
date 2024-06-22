//
// Created by faliszewskii on 17.06.24.
//

#ifndef MATERIAL_PAINTER_TRIANGULATION_H
#define MATERIAL_PAINTER_TRIANGULATION_H

#include <vector>
#include <eigen3/Eigen/Core>

namespace TAGen {
    struct Triangulation {
        struct Triangle {
            std::array<Eigen::Vector2f, 3> vertices;
            std::array<unsigned int, 3> indices;
        };

        std::vector<Triangle> triangles;
        int verticesCount;
    };
}

#endif //MATERIAL_PAINTER_TRIANGULATION_H
