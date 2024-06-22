//
// Created by faliszewskii on 17.06.24.
//

#include "LSCM.h"
#include <ranges>
#include <algorithm>
#include <Eigen/Core>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

TAGen::Triangulation TAGen::LSCM::getFullTriangulation(TAGen::Mesh mesh, /*TODO DEBUG*/Logger &logger) {
    using namespace Eigen;
    auto view = mesh.indices | std::views::chunk(3);

    Triangulation triangulation;
    std::ranges::transform(view, std::back_inserter(triangulation.triangles), [&](const auto& chunk) {
        auto i0 = chunk[0];
        auto i1 = chunk[1];
        auto i2 = chunk[2];
        Vector3f p0 = mesh.vertices[i0];
        Vector3f p1 = mesh.vertices[i1];
        Vector3f p2 = mesh.vertices[i2];

        Vector3f normal = (p1 - p0).cross(p2 - p0).normalized();
        Vector3f tangent = (p1 - p0).normalized();
        Vector3f biTangent = tangent.cross(normal);

        Matrix3f basisAxis;
        basisAxis <<
                  tangent.x(), tangent.y(), tangent.z(),
                normal.x(), normal.y(), normal.z(),
                biTangent.x(), biTangent.y(), biTangent.z();
//            tangent.x(), biTangent.x(), normal.x(),
//            tangent.y(), biTangent.y(), normal.y(),
//            tangent.z(), biTangent.z(), normal.z();

        Vector3f translation = -basisAxis * p0;

        Matrix4f basisChange;
        basisChange.setIdentity();
        basisChange.block<3,3>(0, 0) = basisAxis;
        basisChange.block<3, 1>(0, 3) = translation;

        Vector4f v0, v1, v2;
        v0 << p0, 1;
        v1 << p1, 1;
        v2 << p2, 1;

        auto v0Prim = basisChange * v0;
        auto v1Prim = basisChange * v1;
        auto v2Prim = basisChange * v2;

        Triangulation::Triangle triangle;

        triangle.indices[0] = i0;
        triangle.indices[1] = i1;
        triangle.indices[2] = i2;
        triangle.vertices[0] = Vector2f(v0Prim.x(), v0Prim.z());
        triangle.vertices[1] = Vector2f(v1Prim.x(), v1Prim.z());
        triangle.vertices[2] = Vector2f(v2Prim.x(), v2Prim.z());

        return triangle;
    });
    triangulation.verticesCount = mesh.vertices.size();
    return triangulation;
}

std::expected<TAGen::UvMap, std::string> TAGen::LSCM::parametrizeChart(Triangulation triangulation, Logger &logger) {
    using namespace Eigen;

    int n = triangulation.verticesCount;
    int nPrim = triangulation.triangles.size();

    MatrixXd MReal;
    MatrixXd MImag;
    MReal.resize(nPrim, n);
    MImag.resize(nPrim, n);
    MReal.setZero();
    MReal.setZero();

    for(int i = 0; i < nPrim; i++) {
        auto &vertices = triangulation.triangles[i].vertices;
        auto &indices = triangulation.triangles[i].indices;

        float d = (vertices[0].x() * vertices[1].y() - vertices[0].y() * vertices[1].x()) +
                  (vertices[1].x() * vertices[2].y() - vertices[1].y() * vertices[2].x()) +
                  (vertices[2].x() * vertices[0].y() - vertices[2].y() * vertices[0].x());
        float dSqrt = sqrt(abs(d));

        for(int j = 0; j < 3; j++) {
            float WReal, WImag;
            int prev = (j + 2) % 3;
            int next = (j + 1) % 3;

            WReal = vertices[prev].x() - vertices[next].x();
            WImag = vertices[prev].y() - vertices[next].y();
            auto t = WReal / dSqrt;
            auto t1 = WImag / dSqrt;
            if(dSqrt == 0) return std::unexpected("Area of a triangle is 0!");
            MReal(i, indices[j]) = t;
            MImag(i, indices[j]) = t1;
        }
    }

    int p  = 2;

    MatrixXd MRealP = MReal.block(0, 0, nPrim, p);
    MatrixXd MRealF = MReal.block(0, p, nPrim, n - p);
    MatrixXd MImagP = MImag.block(0, 0, nPrim, p);
    MatrixXd MImagF = MImag.block(0, p, nPrim, n - p);

    MatrixXd A;
    A.resize(2 * nPrim, 2 * (n - p));
    A.block(0, 0, nPrim, n - p) = MRealF;
    A.block(nPrim, 0, nPrim, n - p) = MImagF;
    A.block(0, n - p, nPrim, n - p) = -MImagF;
    A.block(nPrim, n - p, nPrim, n - p) = MRealF;

    MatrixXd tempM;
    tempM.resize(2 * nPrim, 2 * p);
    tempM.block(0, 0, nPrim, p) = MRealP;
    tempM.block(nPrim, 0, nPrim, p) = MImagP;
    tempM.block(0, p, nPrim, p) = -MImagP;
    tempM.block(nPrim, p, nPrim, p) = MRealP;

    VectorXd uValues;
    uValues.resize(p);
    uValues(0) = {0.0};
    uValues(1) = {0.0};
    VectorXd vValues;
    vValues.resize(p);
    vValues(0) = {0.0};
    vValues(1) = {1.0};

    VectorXd U;
    U.resize(2 * p);
    U << uValues, vValues;

    VectorXd b = -tempM * U;
    assert(b.rows() == 2 * nPrim && b.cols() == 1);

//    VectorXd x = (A.transpose() * A).inverse() * A.transpose() * b;

    SparseMatrix<double> ASparce(A.rows(), A.cols());
    for(int row = 0; row < A.rows(); row++)
        for(int col = 0; col < A.cols(); col++)
            if(A(row, col) != 0)
                ASparce.insert(row, col) = A(row, col);

    LeastSquaresConjugateGradient<SparseMatrix<double>> lscg;
    lscg.compute(ASparce);
    VectorXd x = lscg.solve(b);


//    VectorXd x = A.completeOrthogonalDecomposition().solve(b);

    VectorXd xU = x.head(nPrim);
    VectorXd xV = x.tail(nPrim);

    std::vector<Vector2f> uv;
    for(int i = 0; i < p; i++)
        uv.emplace_back(uValues(i), vValues(i));
    for(int i = 0; i < nPrim; i++)
        uv.emplace_back(xU(i), xV(i));

    return UvMap{uv};
}

