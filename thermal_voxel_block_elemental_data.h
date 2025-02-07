#pragma once

#include "kratos.h"

template <int BlockSize>
class ThermalVoxelBlockElementalData {
    std::array<double, BlockSize * BlockSize * BlockSize> mConductivities;
    std::array<double, BlockSize * BlockSize * BlockSize> mDensities;
public:
    ThermalVoxelBlockElementalData() {
        std::fill(mConductivities.begin(), mConductivities.end(), 0.0);
        std::fill(mDensities.begin(), mDensities.end(), 0.0);
    }

    template <class TVariableType>
    auto& Get(TVariableType const& rVariable) {
        if (rVariable == Kratos::CONDUCTIVITY) {
            return mConductivities;
        }
        else if (rVariable == Kratos::DENSITY) {
            return mDensities;
        }
        throw std::runtime_error("Variable not found");
    }
};
