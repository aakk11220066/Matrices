#ifndef EX3_MTMMAT_H
#define EX3_MTMMAT_H


#include <vector>
#include "MtmExceptions.h"
#include "Auxilaries.h"
#include "MtmVec.h"

using std::size_t;

namespace MtmMath {

    template <typename T>
    class MtmMat : MtmVec<MtmVec<T>> {
    protected:
        const Dimensions dimensions; //Roi
    public:
        /*
         * Matrix constructor, dim_t is the dimension of the matrix and val is the initial value for the matrix elements
         */
        MtmMat<T>(Dimensions dim_t, const T& val=T()); //implemented

        //copy constructor
        MtmMat<T>(const MtmMat<T>& original) = default;

        //assignment
        MtmMat<T>& operator=(const MtmMat& original) = default;

        //destructor (Note: MtmVec matrix will destroy its own sub-lists)
        virtual ~MtmMat<T>() = default;

        //Roi :
        //operator[] inherited from MtmVec

        /* Matrix multiplication
         * creates new matrix
         * technique: when multiplying matrices A*B
         * row i of new matrix equals sum of scalar products of jth element of
         * row i of A * row j of B
         */
        MtmMat<T> operator*(MtmMat<T> other);

        /*
         * Function that get function object f and uses it's () operator on each element in the matrix columns.
         * It outputs a vector in the size of the matrix columns where each element is the final output
         * by the function object's * operator
         */
        template <typename Func>
        MtmVec<T> matFunc(Func& f) const;

        /*
         * resizes a matrix to dimension dim, new elements gets the value val.
         */
        virtual void resize(Dimensions dim, const T& val=T());

        /*
         * reshapes matrix so linear elements value are the same without changing num of elements.
         */
        virtual void reshape(Dimensions newDim);

        /*
         * Performs transpose operation on matrix
         */
        virtual void transpose();

    };
}

//implementations begin here
using namespace MtmMath;

template<typename T>
MtmMat<T>::MtmMat(Dimensions dim_t, const T &val) : dimensions(dim_t),
    MtmVec<MtmVec<T>>(dim_t.getRow(), MtmVec<T>(dim_t.getCol(), val)) {
        //TODO: make sure vectors are horizontal!  Perhaps it would be prudent to transpose them here, for instance...
    }

template <typename T>
MtmMat<T> MtmMat<T>::operator*(MtmMat<T> other){
    const size_t n = dimensions.getCol();
    if (n != other.dimensions.getRows()){
        throw MtmExceptions::DimensionMismatch(dimensions, other.dimensions);
    }

    const size_t numRows = dimensions.getRow(),numCols = other.dimensions.getCol();
    MtmMat<T> answer(Dimensions(numRows, numCols), 0);
    for (int answerRow=0; answerRow < numRows; ++answerRow){
        for (int runner=0; runner<n; ++runner){
            answer[answerRow] = answer[answerRow]
                    + ((*this)[answerRow][runner] * other[runner]);
        }
    }

    return answer;
}


#endif //EX3_MTMMAT_H
