#pragma once

#include "kratos.h"

template <int BlockSize>
class ThermalVoxelBlockNodalData {
    std::array<double, BlockSize * BlockSize * BlockSize> mTemperatures;
    std::array<int, BlockSize * BlockSize * BlockSize> mEquationIds;
    std::array<bool, BlockSize * BlockSize * BlockSize> mIsFixed;
public:
    ThermalVoxelBlockNodalData() {
        std::fill(mTemperatures.begin(), mTemperatures.end(), 0.0);
        std::fill(mEquationIds.begin(), mEquationIds.end(), 0);
    }

    template <class TVariableType>
    auto& Get(TVariableType const& rVariable) {
        if (rVariable == Kratos::TEMPERATURE) {
            return mTemperatures;
        }
    }

    int& GetEquationId(int i) {
        return mEquationIds[i];
    }

    const int& GetEquationId(int i) const {
        return mEquationIds[i];
    }


    bool IsFixed(int i) const {
        return mIsFixed[i];
    }

    void Fix(int i) {
        mIsFixed[i] = true;
    }

    void Free(int i) {
        mIsFixed[i] = false;
    }
};
