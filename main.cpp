#include <iostream>
#include <array>
#include <vector>
#include <algorithm>
#include <execution>
#include <chrono>
#include <string>
#include <unordered_map>
#include <iomanip>

#include "includes/thermal_voxel_block_nodal_data.h"
#include "includes/thermal_voxel_block_elemental_data.h"
#include "includes/voxel_block.h"
#include "includes/voxel_block_mesh.h" 
#include "includes/formulation.h"
class Timer {
public:
    Timer() : mStart(std::chrono::high_resolution_clock::now()) {}
    void Reset() { mStart = std::chrono::high_resolution_clock::now(); }
    double Elapsed() const {
        return std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - mStart).count();
    }
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> mStart;
};

template <class T>
double norm_2(const std::vector<T>& rVector) {
    double result = 0.0;
    for (const auto& value : rVector) {
        result += value * value;
    }
    return std::sqrt(result);
}

template<int N>
static double MultiplyBenchmark() {
    std::cout << N << '\t';

    Timer build_timer;

    const int size = 256 / N;
    double element_size = 0.216; // 0.1 is the size of the domain
    const std::array<double, 3> min_point = {-0.1, -0.1, -0.1};
    VoxelBlockMesh<VoxelBlock<ThermalVoxelBlockNodalData, ThermalVoxelBlockElementalData, N>> mesh(min_point,element_size,size, size, size);
    ThermalVoxelBlockFormulation<decltype(mesh)> formulation(mesh);
    // seting the conductivity
    for (auto& block : mesh) {
        auto& conductivities = block.GetElementalData().Get(Kratos::CONDUCTIVITY);
        auto& densities = block.GetElementalData().Get(Kratos::DENSITY);
        for(int i = 0 ; i < block.size() ; i++) {
            conductivities[i] = 100; 
            densities[i] = 1000; // 2700 is the density of aluminum
        }
    }
    const int number_of_nodes = mesh.Size() * (N + 1) * (N + 1) * (N + 1);
    std::cout << int(number_of_nodes/1e6) << "M\t\t";
    std::vector<double> x_vector(number_of_nodes, 1.0);
    // add perturmbation to x_vector
    // for(int i = 0 ; i < number_of_nodes ; i+=2) {
    //     x_vector[i] =-1.0;
    // }

    ProcessInfo process_info;
    std::vector<double> rhs_vector(number_of_nodes, 0.0);
    formulation.SetEquationId(process_info);
    std::cout << std::fixed << std::setprecision(4) << build_timer.Elapsed() << '\t' << '\t'; // Set precision for this output

    Timer multiply_timer;
    for (int i = 0; i < 100; i++)
        formulation.LHSMult(x_vector, rhs_vector);
    std::cout << std::fixed << std::setprecision(4) << multiply_timer.Elapsed() << '\t'; // Set precision for this output

    double norm = norm_2(rhs_vector);
    std::cout << "\t" << norm << std::endl;
    return norm;
}

int main() {
    std::cout << "N\tN. Nodes\tBuild\t\tMultiply\tNorm2" << std::endl;
    MultiplyBenchmark<4>();
    MultiplyBenchmark<8>();
    // MultiplyBenchmark<16>();
    // MultiplyBenchmark<32>();
    // MultiplyBenchmark<64>();
    return 0;
}
