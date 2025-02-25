#pragma once

#include <vector>

class ProcessInfo {
public:
    double DeltaTime() const { return 0.01; }
};

class CSRMatrix {
    std::vector<int> mIndex1;
    std::vector<int> mIndex2;
    std::vector<double> mValues;
public:
    int* index1_data() { return mIndex1.data(); }
    int* index2_data() { return mIndex2.data(); }
    double* values_data() { return mValues.data(); }
};

using Vector = std::vector<double>;

class Formulation {
public:
    virtual void Build(CSRMatrix& rLHSMatrix, Vector& rRHSVector, ProcessInfo& rProcessInfo) = 0;
    virtual void SetEquationId(ProcessInfo& rProcessInfo) = 0;
    virtual void SolveExplicit(ProcessInfo& rProcessInfo) = 0;
    virtual void LHSMult(const Vector& rX, Vector& rResult) = 0;
};;

// The ThermalVoxelBlockFormulation class is a specialization of the Formulation class for thermal problems
template <class TVoxelBlockMeshType>
class ThermalVoxelBlockFormulation : public Formulation {
    TVoxelBlockMeshType& mMesh;
public:
    using VoxelBlockType = typename TVoxelBlockMeshType::VoxelBlockType;

    ThermalVoxelBlockFormulation(TVoxelBlockMeshType& TheMesh) : mMesh(TheMesh) {}

    void Build(CSRMatrix& rLHSMatrix, Vector& rRHSVector, ProcessInfo& rProcessInfo) override {
        // // initializing the lhs matrix to have 27 nonzeros per row
        // // fill the index1 array with row index to start at 0 and increment by 27
        // rLHSMatrix.index1_data()[0] = 0;
        // for (int i = 1; i < mMesh.Size(); i++) {
        //     rLHSMatrix.index1_data()[i] = rLHSMatrix.index1_data()[i-1] + 27;
        // }

        // // Calculating the equation id for each node
        // // Here we can also perform the elimination of the fixed and slave dofs
        // int equation_id = 0;
        // for (auto& block : mMesh) {
        //     for (int i = 0; i < block.BlockSize() + 1; i++) { // Todo: change to 1 continuous loop
        //         for (int j = 0; j < block.BlockSize() + 1; j++) {
        //             for (int k = 0; k < block.BlockSize() + 1; k++) {
        //                 block.GetNodalData(i, j, k).SetEquationId(equation_id++);
        //             }
        //         }
        //     }
        // }

        // constexpr char node_i_offset[8]={0, 1, 1, 0, 0, 1, 1, 0};
        // constexpr char node_j_offset[8]={0, 0, 1, 1, 0, 0, 1, 1};
        // constexpr char node_k_offset[8]={0, 0, 0, 0, 1, 1, 1, 1};
        // constexpr char assemble_row_positions[64]={13,14,17,16,22,23,26,25, 12,13,16,15,21,22,25,24, 9,10,13,12,18,19,22,21, 10,11,14,13,19,20,23,22, 4,5,8,7,13,14,17,16, 3,4,7,6,12,13,16,15, 0,1,4,3,9,10,13,12, 1,2,5,4,10,11,14,13};
        // constexpr char k_coeficients[64]={4, 0, -1, 0, 0, -1, -1, -1, 0, 4, 0, -1, -1, 0, -1, -1, -1, 0, 4, 0, -1, -1, 0, -1, 0, -1, 0, 4, -1, -1, -1, 0, 0, -1, -1, -1, 4, 0, -1, 0, -1, 0, -1, -1, 0, 4, 0, -1, -1, -1, 0, -1, -1, 0, 4, 0, -1, -1, -1, 0, 0, -1, 0, 4};
        // constexpr int number_of_elements_in_block = VoxelBlockType::size();

        // // loop over all blocks in the mesh
        // std::for_each(std::execution::par_unseq, mMesh.begin(), mMesh.end(), [&](VoxelBlockType& block) {
        //     // loop over elements in the block
        //     for(int i = 0 ; i < number_of_elements_in_block ; i++) {
        //         // get the elemental data
        //         auto& elemental_data = block.GetElementalData(i);
        //         double conductivity = elemental_data;
        //         // loop over nodes in the element
        //         for(int l = 0; l < 8; l++) {
        //             // get the nodal data
        //             auto& nodal_data = block.GetElementNodalData(i, l);
        //             double temperature = nodal_data.GetValue(Kratos::TEMPERATURE);
        //             int equation_id = nodal_data.GetEquationId();
        //             int row_index = rLHSMatrix.index1_data()[equation_id];
        //             for(int m = 0; m < 8; m++) {
        //                 rLHSMatrix.value_data()[row_index + assemble_row_positions[l * 8 + m]] += conductivity * static_cast<double>(k_coeficients[l * 8 + m]);
        //             }

        //         }
        //     }
        // });
    }

    void SetEquationId(ProcessInfo& rProcessInfo){
        // Calculating the equation id for each node
        // Here we can also perform the elimination of the fixed and slave dofs
        int equation_id = 0;
        for (auto& block : mMesh) {
            block.SetStartEquationId(equation_id);
            equation_id += block.NumberOfNodes();

            // for (int i = 0; i < block.BlockSize() + 1; i++) { // Todo: change to 1 continuous loop
            //     for (int j = 0; j < block.BlockSize() + 1; j++) {
            //         for (int k = 0; k < block.BlockSize() + 1; k++) {
            //             block.GetNodalData(i, j, k).SetEquationId(equation_id++);
            //         }
            //     }
            // }
        }
    }    

    void SolveExplicit(ProcessInfo& rProcessInfo) override {}
    void LHSMult(const Vector& rX, Vector& rResult) override{
        constexpr int number_of_elements_in_block = VoxelBlockType::size();
        constexpr int block_size = VoxelBlockType::BlockSize();

        const double l = mMesh.GetElementSize();


        std::for_each(std::execution::par_unseq, mMesh.begin(), mMesh.end(), [&](VoxelBlockType& block) {
            const auto& conductivities = block.GetElementalData().Get(Kratos::CONDUCTIVITY);
            const auto& densities = block.GetElementalData().Get(Kratos::DENSITY);

            for(int k = 0 ; k < block_size ; k++) {
                for(int j = 0 ; j < block_size ; j++) {
                    const int elemental_index = j * block_size + k * block_size * block_size;

                    const int node_0_index(j * (block_size + 1) + k * (block_size + 1) * (block_size + 1));
                    const int node_1_index=node_0_index+1;
                    const int node_2_index=node_0_index+block_size + 2;
                    const int node_3_index=node_0_index+block_size + 1;
                    const int node_4_index=node_0_index+(block_size + 1) * (block_size + 1);
                    const int node_5_index=node_0_index+1+(block_size + 1) * (block_size + 1);
                    const int node_6_index=node_0_index+block_size + 2+(block_size + 1) * (block_size + 1);
                    const int node_7_index=node_0_index+block_size + 1+(block_size + 1) * (block_size + 1);
                    const int start_id = block.GetStartEquationId();
                    const int equation_id_0 = start_id + node_0_index;
                    const int equation_id_1 = start_id + node_1_index;
                    const int equation_id_2 = start_id + node_2_index;
                    const int equation_id_3 = start_id + node_3_index;
                    const int equation_id_4 = start_id + node_4_index;
                    const int equation_id_5 = start_id + node_5_index;
                    const int equation_id_6 = start_id + node_6_index;
                    const int equation_id_7 = start_id + node_7_index;

                    for(int i = 0 ; i < block_size ; i++) {
                        // get the elemental data
                        // auto& elemental_data = block.GetElementalData(i, j, k);
                        // double conductivity = elemental_data;
                        // loop over nodes in the element

                        double d = conductivities[elemental_index + i];
                        double rho = densities[elemental_index + i];

                        const double& X_0 = rX[equation_id_0+i];
                        const double& X_1 = rX[equation_id_1+i];
                        const double& X_2 = rX[equation_id_2+i];
                        const double& X_3 = rX[equation_id_3+i];
                        const double& X_4 = rX[equation_id_4+i];
                        const double& X_5 = rX[equation_id_5+i];
                        const double& X_6 = rX[equation_id_6+i];
                        const double& X_7 = rX[equation_id_7+i];

                        const double x0 = d*l;
                        const double l3 = l*l*l;
                        const double x2 = rho*l3;
                        const double x3 = (1.0/3.0)*x0 + (1.0/27.0)*x2;
                        const double x4 = (1.0/12.0)*x0;
                        const double x5 = (1.0/216.0)*rho*l3 - x4;
                        const double x6 = (1.0/108.0)*rho*l3 - x4;
                        const double x7 = X_5*x6;
                        const double x8 = X_7*x6;
                        const double x9 = (1.0/54.0)*x2;
                        const double x10 = X_1*x9;
                        const double x11 = X_3*x9;
                        const double x12 = x10 + x11 + x7 + x8;
                        const double x13 = X_2*x6 + X_4*x9;
                        const double x14 = X_4*x6;
                        const double x15 = X_6*x6;
                        const double x16 = X_0*x9;
                        const double x17 = X_2*x9;
                        const double x18 = x14 + x15 + x16 + x17;
                        const double x19 = X_3*x6 + X_5*x9;
                        const double x20 = X_0*x6 + X_6*x9;
                        const double x21 = X_1*x6 + X_7*x9;
                        const double x22 = x19 + x21;
                        const double x23 = x13 + x20;
                        rResult[equation_id_0+i]+= X_0*x3 + X_6*x5 + x12 + x13 ;
                        rResult[equation_id_1+i]+= X_1*x3 + X_7*x5 + x18 + x19 ;
                        rResult[equation_id_2+i]+= X_2*x3 + X_4*x5 + x12 + x20 ;
                        rResult[equation_id_3+i]+= X_3*x3 + X_5*x5 + x18 + x21 ;
                        rResult[equation_id_4+i]+= X_2*x5 + X_4*x3 + x15 + x16 + x22 ;
                        rResult[equation_id_5+i]+= X_3*x5 + X_5*x3 + x10 + x23 + x8 ;
                        rResult[equation_id_6+i]+= X_0*x5 + X_6*x3 + x14 + x17 + x22 ;
                        rResult[equation_id_7+i]+= X_1*x5 + X_7*x3 + x11 + x23 + x7 ;
                    }
                }
            }
        });
    }

        std::array<double,27> CalculateNodeContribution(int I, int J, int K, VoxelBlockType& Block) {
        std::array<double,27> lhs_row;
        int index = 0;
        constexpr char k_coeficients[64]={4, 0, -1, 0, 0, -1, -1, -1, 0, 4, 0, -1, -1, 0, -1, -1, -1, 0, 4, 0, -1, -1, 0, -1, 0, -1, 0, 4, -1, -1, -1, 0, 0, -1, -1, -1, 4, 0, -1, 0, -1, 0, -1, -1, 0, 4, 0, -1, -1, -1, 0, -1, -1, 0, 4, 0, -1, -1, -1, 0, 0, -1, 0, 4};
        constexpr char reference_node_local_indices[8]={6, 7, 5, 4, 2, 3, 1, 0};
        constexpr char elements_node_positions[64]={0,1,4,3,9,10,13,12, 1,2,5,4,10,11,14,13, 3,4,7,6,12,13,16,15, 4,5,8,8,13,14,17,16, 9,10,13,12,18,19,22,21, 10,11,14,13,19,20,23,22, 12,13,16,15,21,22,25,24, 13,14,17,16,22,23,26,25};

        // Loop over elements in relative positions -1,0 in k,j,i
        for(int k = -1; k < 1; k++) {
            for(int j = -1; j < 1; j++) {
                for(int i = -1; i < 1; i++) {
                    if(Block.HasElementInRelativePosition(I, J, K, i, j, k)) {
                        double conductivity = Block.GetElementalData(I + i, J + j, K + k);
                        int reference_node_index = reference_node_local_indices[index];
                        int reference_node_equation_id = Block.GetNodalData(I + i, J + j, K + k).GetEquationId();
                        int reference_node_k_index = reference_node_index * 8;
                        int element_node_index = index * 8;
                        for(int l = 0; l < 8; l++) {
                            lhs_row[elements_node_positions[element_node_index+l]] += conductivity * static_cast<double>(k_coeficients[reference_node_k_index + l]);
                        }
                    }
                    index++;
                }
            }
        }
        return lhs_row;
    }

    std::array<int,27> CalculateNodeContributionColumn(int I, int J, int K, VoxelBlockType& Block) {
        std::array<int,27> lhs_col;
        int index = 0;
        int current_equation_id = Block.GetNodalData(I, J, K).GetEquationId();

        // Loop over relative positions -1,0,1 in i, j, k
        for(int k = -1; k < 1; k++) {
            for(int j = -1; j < 1; j++) {
                for(int i = -1; i < 1; i++) {
                    if (Block.HasNodeInRelativePosition(I, J, K, i, j, k)) {
                        lhs_col[index] = Block.GetNodalData(I + i, J + j, K + k).GetEquationId();
                    }
                    else { // If the node does not exist we will have 0 in this position
                        lhs_col[index] = current_equation_id; // Adding the diagonal element for 0 values to have cache efficiency
                    }
                    index++;
                }
            }
        }

        return lhs_col;
    }


    double CalculateNodeRhs(int i, int j, int k, VoxelBlockType& block) {
        return 0.0;
    }

};