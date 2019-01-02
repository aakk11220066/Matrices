#ifndef EX3_MTMMATTRIAG_H
#define EX3_MTMMATTRIAG_H


#include <vector>
#include "MtmExceptions.h"
#include "Auxilaries.h"
#include "MtmMatSq.h"

using std::size_t;
using namespace MtmMath;

namespace MtmMath {


    template <typename T>
    class MtmMatTriag : public MtmMatSq {
        //checks if given MtmMatSq is also upper triangular
        bool isUpper(const MtmMatSq<T>& mat) const;

        //zeroes out top/bottom triangle of matrix
        void triangulate(MtmMatSq<T>& mat, bool makeUpper) const;
    public:

        /*
         * Triangular Matrix constructor, m is the number of rows and columns in the matrix,
         * val is the initial value for the matrix elements and isUpper_ is whether it is upper
         * Rectangular matrix (true means it is)
         */
        MtmMatTriag<T> (size_t m, const T& val=T(), bool isUpper_t=true);

        //copy constructor
        MtmMatTriag<T> (const MtmMatTriag& original) = default;

        //destructor
        ~MtmMatTriag<T> () = default;

        //operator=
        MtmMatTriag<T>& operator=(const MtmMatTriag&) = default;

        //constructor for normal/square matrices to triangular
        MtmMatTriag<T>& (const MtmMatSq& original);
    };

}


//implementation begins here
template <typename T>
MtmMatTriag<T>::MtmMatTriag(size_t m, const T& val, bool isUpper_t) :
    MtmMatSq(m, val){

    triangulate (*this, isUpper_t);
}

template <typename T>
MtmMatTriag<T>::MtmMatTriag(const MtmMatSq& original) :
    MtmMatSq(original) {

    //if upper - triangulate to upper
    //if lower - leave, otherwise triangulate to upper TODO: MAKE SURE and then finish function
    if isUpper(original) return;
    ...
}

template <typename T>
void MtmMatTriag<T>::triangulate(MtmMatSq<T> &mat, bool makeUpper) const {
    int m = mat.dimensions.getRow();
    for (int row = 0; row < m; ++row) { //zero out top or bottom of matrix
        for (int col = 0; col < m; ++col) {
            if ((row < col && !makeUpper) || (row > col && makeUpper)) {
                mat[row][col] = 0;
            }
        }
    }
}

template<typename T>
bool MtmMatTriag<T>::isUpper(const MtmMatSq<T> &mat) const {
    int m = mat.dimensions.getRow();
    for (int row = 0; row < m; ++row) {
        for (int col = 0; col < m; ++col) {
            if (row>col && mat[row][col]!=0) { //found nonzero beneath diagonal
                return false;
            }
        }
    }
    return true;
}


#endif //EX3_MTMMATTRIAG_H
