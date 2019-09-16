#include "Stripes.h"
#include <stripes/Mesh.h>
#include <stripes/MeshIO.h>

namespace DDG {

namespace {

class MeshEigen : public Mesh {
public:
    int init(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
    void setDirectionField(const Eigen::MatrixXd &D);
};

// Return 0 for success
int MeshEigen::init(const Eigen::MatrixXd &V,
                const Eigen::MatrixXi &F)
{
    MeshData meshData;
    meshData.positions.reserve(V.rows());
    for (int i = 0; i < V.rows(); ++i) {
        meshData.positions.emplace_back(V(i, 0), V(i, 1), V(i, 2));
    }
    meshData.indices.reserve(F.rows());
    for (int i = 0; i < F.rows(); ++i) {
        meshData.indices.emplace_back(F.cols());
        for (int lv = 0; lv < F.cols(); ++lv) {
            meshData.indices.back()[lv] = Index(F(i, lv), -1, -1);
        }
    }
    if (auto ret = MeshIO::buildMesh(meshData, *this)) {
        return ret;
    }
    indexElements();
    initializeSmoothStructure();
    buildMassMatrix();

    return 0;
}

void MeshEigen::setDirectionField(const Eigen::MatrixXd &D) {
    // Ignore input direction field for now
    computeCurvatureAlignedSection();
}

}  // anonymous namespace

int computeStripePatterns(const Eigen::MatrixXd &V,
                          Eigen::MatrixXi &F,
                          const Eigen::MatrixXd &directionField,
                          double frequency,
                          int numCoords,
                          Eigen::VectorXi &branchIndex,
                          Eigen::MatrixXd &parameterization,
                          Eigen::MatrixXi &zeroIndex,
                          std::vector<bool> &isBorder)
{
    assert(V.cols() == 3);
    assert(F.cols() == 3);
    if (!(numCoords == 1 || numCoords == 2)) {
        throw std::runtime_error("Invalid number of coordinate functions.");
    }

    // Compute field-aligned parametrization
    MeshEigen mesh;
    if (auto ret = mesh.init(V, F)) {
        return ret;
    }
    mesh.setDirectionField(directionField);
    mesh.nCoordinateFunctions = 2;
    mesh.lambda = frequency;
    mesh.parameterize();
    if (numCoords == 2) {
        mesh.glueParameterization();
    }

    // Convert result
    branchIndex.resize(F.rows());
    branchIndex.setZero();
    parameterization.resize(3 * F.rows(), numCoords);
    parameterization.setZero();
    zeroIndex.resize(F.rows(), numCoords);
    zeroIndex.setZero();
    isBorder.resize(F.rows());

    int faceIndex = 0;
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++, ++faceIndex) {
        isBorder[faceIndex] = f->isBoundary();
        if (f->isBoundary()) continue;

        double k = f->fieldIndex(2.);
        branchIndex(faceIndex) = k;
        for (int p = 0; p < numCoords; ++p) {
            zeroIndex(faceIndex, p) = f->paramIndex[p];
        }

        // We don't do anything special for faces with singularities for now
        int i = 0;
        HalfEdgeCIter he = f->he;
        do {
            Complex g = he->texcoord;
            F(faceIndex, i) = he->vertex->index;
            parameterization(3 * faceIndex + i, 0) = g.re;
            if (numCoords == 2) {
                parameterization(3 * faceIndex + i, 1) = g.im;
            }
            i++;
            he = he->next;
        } while (he != f->he);
    }

    return 0;
}

}  // namespace DDG
