//
// Created by faliszewskii on 17.06.24.
//

#ifndef MATERIAL_PAINTER_UVMAP_H
#define MATERIAL_PAINTER_UVMAP_H

#include <vector>
#include <Eigen/Core>

namespace TAGen {
    struct UvMap {
        std::vector<Eigen::Vector2f> uv;
    };
}

#endif //MATERIAL_PAINTER_UVMAP_H
