#pragma once

// VoxelMesh class which stores a 3D array of VoxelBlocks
template <class TVoxelBlockType>
class VoxelMesh {
    int mSizeX, mSizeY, mSizeZ;
    double mElementSize;
    std::vector<TVoxelBlockType> mBlocks;
    std::array<double, 3> mMinPoint;

public:
    using VoxelBlockType = TVoxelBlockType;

    template<typename TPointType>
    VoxelMesh(TPointType const& MinPoint, double ElementSize, int SizeX, int SizeY, int SizeZ)
    : mSizeX(SizeX), mSizeY(SizeY), mSizeZ(SizeZ)
    , mElementSize(ElementSize)
    , mMinPoint({MinPoint[0], MinPoint[1], MinPoint[2]}),
    mBlocks(SizeX * SizeY * SizeZ, TVoxelBlockType(ElementSize)) {}


    TVoxelBlockType& GetBlock(int x, int y, int z) {
        return mBlocks[x + y * mSizeX + z * mSizeX * mSizeY];
    }

    const TVoxelBlockType& GetBlock(int x, int y, int z) const {
        return mBlocks[x + y * mSizeX + z * mSizeX * mSizeY];
    }

    int SizeX() const {
        return mSizeX;
    }

    int SizeY() const {
        return mSizeY;
    }

    int SizeZ() const {
        return mSizeZ;
    }

    int Size() const {
        return mBlocks.size();
    }

    TVoxelBlockType& operator[](int i) {
        return mBlocks[i];
    }

    const TVoxelBlockType& operator[](int i) const {
        return mBlocks[i];
    }

    auto begin() {
        return mBlocks.begin();
    }

    auto end() {
        return mBlocks.end();
    }

    auto begin() const {
        return mBlocks.begin();
    }

    auto end() const {
        return mBlocks.end();
    }

    double GetElementSize() const {
        return mElementSize;
    }
};
