#ifndef CLUSTER_ALGORITHMS_H
#define CLUSTER_ALGORITHMS_H

#include <vector>

// unsigned 128-bit integer type
typedef unsigned __int128 uint128_t;

// A simple 2D array class
template <typename T>
class Array2D {
    public:
        std::vector<T> data;
        int ncols;
        Array2D(int nrow, int ncol) : data(nrow * ncol), ncols(ncol) {}
        T* operator[](int row) {
            return &data[row * ncols];
        }
};

// A simple upper triangular matrix class
template <typename T>
class UpperTriangularMatrix {
    public:
        std::vector<T> data;
        int size;
        UpperTriangularMatrix(int n) : data((n * (n - 1)) / 2), size(n) {}
        T* operator[](int row) {
            if (row < size - 1) {
                // row r starts after rows 0..r-1, of lengths size-1, ..., size-r
                return &data[(row * (size - 1)) - ((row * (row - 1)) / 2)];
            } else {
                return nullptr; // The last row has no i<j pairs
            }
        }
};

#endif
