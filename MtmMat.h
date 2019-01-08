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
        const Dimensions dimensions;
    private:
        //generates a square matrix full of zeroes except for the secondary
        //  diagonal, which has 1s
        MtmMat<T> altMatrix(size_t dimension) const;

        T& linearIndexToReference(unsigned int linearIndex){
            const unsigned int numRows = dimensions.getRow();
            const unsigned int row=linearIndex%numRows, col=linearIndex/numRows;
            if (row>=numRows || col>=dimensions.getCol()){
                throw MtmExceptions::AccessIllegalElement();
            }
            return (*this)[row][col];
        }
        unsigned int coordinatesToLinearIndex(size_t row, size_t col){
            return dimensions.getRow()*col + row;
        }
    public:
        /*
         * Matrix constructor, dim_t is the dimension of the matrix and val is
         * the initial value for the matrix elements
         */
        MtmMat<T>(Dimensions dim_t, const T& val=T()); //implemented

        //copy constructor
        MtmMat<T>(const MtmMat<T>& original) = default;

        //assignment
        virtual MtmMat<T>& operator=(const MtmMat& original) = default;

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
        friend MtmMat<T> operator*(const MtmMat<T> matrix1, //implemented
                const MtmMat<T> matrix2);

        //Matrix addition
        friend MtmMat<T> operator+(const MtmMat<T> matrix1, //implemented
                const MtmMat<T> matrix2);


    /*
     * Function that get function object f and uses it's () operator on each
     * element in the matrix columns.
     * It outputs a vector in the size of the matrix columns where each element
     * is the final output
     * by the function object's * operator
     */
        template <typename Func>
        MtmVec<T> matFunc(Func& f) const;

        /*
         * resizes a matrix to dimension dim, new elements gets the value val.
         */
        virtual void resize(Dimensions dim, const T& val=T());

        /*
         * reshapes matrix so linear elements value are the same without
         * changing num of elements.
         */
        virtual void reshape(Dimensions newDim);

        /*
         * Performs transpose operation on matrix
         */
        virtual void transpose();

        class iterator{
            //Roi : POTENTIAL ERROR: code replication?
        protected:
            MtmMat* self;
        private:
            unsigned int linearIndex;

            friend iterator begin();
            friend iterator end();
            explicit iterator(MtmMat* self, unsigned int startIndex) :
                self(self), linearIndex(startIndex){}
        public:
            T& operator*(){
                return self->linearIndexToReference(linearIndex);
            }
            friend bool operator==(iterator& me, iterator& other);
            friend bool operator!=(iterator& me, iterator& other);
            virtual iterator& operator++(){
                if (*this != self->end()) ++linearIndex;
            }
        };
        iterator begin(){
            return iterator(this, 0);
        }
        iterator end(){
            return iterator(this,coordinatesToLinearIndex(dimensions.getRow()-1,
                    dimensions.getCol()-1) + 1);
        }

        class nonzero_iterator : iterator{
        public:
            nonzero_iterator& operator++() override{
                if (*this != self->end()) iterator::operator++(*this);
                while(*this != self->end() && this->operator*() == 0){
                    iterator::operator++(*this);
                }
            }
        private:
            friend nonzero_iterator nzbegin();
            friend nonzero_iterator nzend();
            explicit nonzero_iterator(MtmMat* self, unsigned int startIndex) :
                iterator(self,startIndex){
                    if (this->operator*() == 0) ++(*this);
            }
        };
        nonzero_iterator nzbegin();
        nonzero_iterator nzend();

    };
}

//implementations begin here
using namespace MtmMath;

template <typename T>
bool MtmMat<T>::iterator::operator!=(iterator& me, iterator& other){
    return !me==other;
}

template <typename T>
bool MtmMat<T>::iterator::operator==(iterator& me, iterator& other){
    return me.linearIndex==other.linearIndex;
}

template<typename T>
MtmMat<T>::MtmMat(Dimensions dim_t, const T &val) : dimensions(dim_t),
    MtmVec<MtmVec<T>>(dim_t.getRow(), MtmVec<T>(dim_t.getCol(), val)) {
        //make sub-vectors horizontal
        for (int row=0; row<dim_t.getRow(); ++row){
            (*this)[row].transpose(); //MtmVec transpose
        }
}

template <typename T>
MtmMat<T> operator*(const MtmMat<T> matrix1, const MtmMat<T> matrix2){
    const size_t n = matrix1.dimensions.getCol();
    if (n != matrix2.dimensions.getRows()){
        throw MtmExceptions::DimensionMismatch(matrix1.dimensions,
                matrix2.dimensions);
    }

    const size_t numRows = matrix1.dimensions.getRow();
    const size_t numCols = matrix2.dimensions.getCol();
    MtmMat<T> answer(Dimensions(numRows, numCols), 0);
    for (int answerRow=0; answerRow < numRows; ++answerRow){
        for (int runner=0; runner<n; ++runner){
            answer[answerRow] = answer[answerRow]
                    + (matrix1[answerRow][runner] * matrix2[runner]);
        }
    }

    return answer;
}

template <typename T>
MtmMat<T> operator+(const MtmMat<T> matrix1, const MtmMat<T> matrix2){
    const size_t numRows = matrix1.dimensions.getRows();
    const size_t numCols = matrix1.dimensions.getCols();
    if (matrix2.dimensions.getRows() != numRows
        || matrix2.dimensions.getCols() != numCols){

        throw MtmExceptions::DimensionMismatch(matrix1.dimensions,
                matrix2.dimensions);
    }

    MtmMat<T> answer(Dimensions(numRows, numCols), 0);
    for (int row=0; row<numRows; ++row){
        answer[row] = matrix1[row] + matrix2[row]; //vector addition
    }
    return answer;
}

template <typename T>
MtmMat<T> MtmMat<T>::altMatrix(size_t dimension) const {
    MtmMat<T> answer(Dimensions(dimension, dimension), 0);
    for (int i=0; i<dimension; i++) {
        answer[dimension-i-1][i] = 1;
    }
    return answer;
}

template <typename T>
void MtmMat<T>::transpose(){
    //POTENTINAL ERROR: replacing matrix (rather than modifying) may reset
    // additional properties of subclasses

    //mathematically proven to be equivalent to transposition
    *this = altMatrix(dimensions.getRow())
            * (*this)
            * altMatrix(dimensions.getCol());
}

template <typename T>
void MtmMat<T>::reshape(Dimensions newDim){
    //POTENTINAL ERROR: replacing matrix (rather than modifying) may reset
    // additional properties of subclasses

    const size_t numRows = newDim.getRow(), numCols = newDim.getCol();
    const unsigned int maxIndex =
            MtmMat<T>::coordinatesToLinearIndex(numRows-1, numCols-1);
    //sanitize inputs
    if ((newDim.getRow()<=0 || newDim.getCol()<=0)
        || (dimensions.getCol()*dimensions.getRow() != numRows*numCols)) {

        throw MtmExceptions::ChangeMatFail(dimensions, newDim);
    }

    MtmMat<T> replacement(newDim, 0);
    for (int index=0; index<maxIndex; ++index){
        replacement.linearIndexToReference(index)=linearIndexToReference(index);
    }

    *this = replacement;
}

template <typename T>
void MtmMat<T>::resize(Dimensions dim, const T& val) {
    //POTENTINAL ERROR: replacing matrix (rather than modifying) may reset
    // additional properties of subclasses
    //TODO:Add input sanitization (for dimension 0, for example)

    MtmMat<T> replacement(dim, val);
    for (int row=0; row<dimensions.getRow() && row<dim.getRow(); ++row){
        for (int col=0; col<dimensions.getCol() && col<dim.getCol(); ++col){
            replacement[row][col] = (*this)[row][col];
        }
    }
    (*this) = replacement;
}


#endif //EX3_MTMMAT_H
