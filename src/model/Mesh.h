//
// Created by faliszewskii on 17.06.24.
//

#ifndef MATERIAL_PAINTER_MESH_H
#define MATERIAL_PAINTER_MESH_H

#include <vector>
#include <eigen3/Eigen/Core>

namespace TAGen {
    struct Mesh {
        std::vector<Eigen::Vector3f> vertices;
        std::vector<unsigned int> indices;
    };
}

#endif //MATERIAL_PAINTER_MESH_H
