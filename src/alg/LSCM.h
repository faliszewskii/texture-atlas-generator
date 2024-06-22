//
// Created by faliszewskii on 17.06.24.
//

#ifndef MATERIAL_PAINTER_LSCM_H
#define MATERIAL_PAINTER_LSCM_H

#include <expected>
#include "../model/UvMap.h"
#include "../model/Triangulation.h"
#include "../model/Mesh.h"
#include "../../../../src/interface/logger/Logger.h"

namespace TAGen {
    class LSCM {
    public:
        Triangulation getFullTriangulation(Mesh mesh, /*TODO DEBUG*/Logger &logger);
        std::expected<UvMap, std::string> parametrizeChart(Triangulation triangulation, Logger &logger);

    };
}

#endif //MATERIAL_PAINTER_LSCM_H
