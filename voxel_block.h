#pragma once

// VoxelBlock class which stores a 3D array of voxels
template <template<int TBlockSize> class TNodalDataType, template<int TBlockSize> class TElemntalDataType, int TBlockSize>
class VoxelBlock {
    int mStartEquationId;
    double mElementSize;
    TElemntalDataType<TBlockSize> mElementalData;
    TNodalDataType<TBlockSize> mNodalData;

public:
    using ElementalDataType = TElemntalDataType<TBlockSize>;
    using NodalDataType = TNodalDataType<TBlockSize>;

    VoxelBlock(double ElementSize) 
    : mElementSize(ElementSize), mElementalData(), mNodalData() {}

    ElementalDataType& GetElementalData() {
        return mElementalData;
    }

    const ElementalDataType& GetElementalData() const {
        return mElementalData;
    }

    NodalDataType& GetNodalData() {
        return mNodalData;
    }

    const NodalDataType& GetNodalData() const {
        return mNodalData;
    }

    int GetNodeIndex(int I, int J, int K) const {
        return I + J * (TBlockSize + 1) + K * (TBlockSize + 1) * (TBlockSize + 1);
    }

    int GetElementNodeIndex(int ElementIndex, int NodeLocalIndex) {
        constexpr int node_i_offset[8]={0, 1, 1, 0, 0, 1, 1, 0};
        constexpr int node_j_offset[8]={0, 0, 1, 1, 0, 0, 1, 1};
        constexpr int node_k_offset[8]={0, 0, 0, 0, 1, 1, 1, 1};

        const int k = ElementIndex / (TBlockSize * TBlockSize);
        const int j = (ElementIndex - k * TBlockSize * TBlockSize) / TBlockSize;
        const int i = ElementIndex - k * TBlockSize * TBlockSize - j * TBlockSize;

        return GetNodeIndex(i + node_i_offset[NodeLocalIndex], j + node_j_offset[NodeLocalIndex], k + node_k_offset[NodeLocalIndex]);
       
        // This is the code with shadow elements
        // constexpr char node_offset[8]={0, 1, (TBlockSize + 2), (TBlockSize + 1), (TBlockSize + 1) * (TBlockSize + 1), (TBlockSize + 1) * (TBlockSize + 1) + 1, (TBlockSize + 1) * (TBlockSize + 1) + (TBlockSize + 2), (TBlockSize + 1) * (TBlockSize + 1) + (TBlockSize + 1)};
        // return ElementIndex + node_offset[NodeLocalIndex];
    }

    constexpr static int size() {
        return  TBlockSize * TBlockSize * TBlockSize;
    }

    constexpr static int BlockSize() {
        return TBlockSize;
    }

    bool HasNodeInRelativePosition(int I, int J, int K, int RelativeI, int RelativeJ, int RelativeK) const {
        return I + RelativeI >= 0 && I + RelativeI < (TBlockSize + 1) &&
               J + RelativeJ >= 0 && J + RelativeJ < (TBlockSize + 1) &&
               K + RelativeK >= 0 && K + RelativeK < (TBlockSize + 1);
    }

    bool HasElementInRelativePosition(int I, int J, int K, int RelativeI, int RelativeJ, int RelativeK) const {
        return I + RelativeI >= 0 && I + RelativeI < TBlockSize &&
               J + RelativeJ >= 0 && J + RelativeJ < TBlockSize &&
               K + RelativeK >= 0 && K + RelativeK < TBlockSize;  
    }

    constexpr static int NumberOfNodes(){
        return (TBlockSize + 1) * (TBlockSize + 1) * (TBlockSize + 1);
    }

    constexpr static int NumberOfElements(){
        return TBlockSize * TBlockSize * TBlockSize;
    }

    constexpr static int NumberOfNodesInElement(){
        return 8;
    }

    int GetStartEquationId() const {
        return mStartEquationId;
    }

    void SetStartEquationId(int StartEquationId) {
        mStartEquationId = StartEquationId;
    }
};
