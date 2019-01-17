#ifndef EX3_MTMMAT_H
#define EX3_MTMMAT_H


#include <vector>
#include "MtmExceptions.h"
#include "Auxilaries.h"
#include "MtmVec.h"
#include "RootVector.h"
#include <cassert>

using std::size_t;

namespace MtmMath {
    const size_t defaultElement = 0, firstIndex = 0, errorValue = 8998;

    //forward declarations
    template <typename T>
    class MtmMat;

    template<typename T>
    class MtmMat : public RootVector<T> {
    protected:
        Dimensions dimensions;

        //generates vector copy of some column of matrix
        const MtmVec<T> getColAsVector(size_t col) const;
    private:
        //generates a square matrix full of zeroes except for the secondary
        //  diagonal, which has 1s
        MtmMat<T> exchangeMatrix(size_t dimension) const;

        //generates reference to element based on lniear index
        T &linearIndexToReference(size_t linearIndex) {
            const size_t numRows = dimensions.getRow();
            const size_t row = linearIndex % numRows, col = linearIndex / numRows;
            return (*this)[row][col];
        }

        //calculates linear index of matrix coordinates
        size_t coordinatesToLinearIndex(size_t row, size_t col) {
            return dimensions.getRow() * dimensions.getCol();
        }

    public:
        /*
         * Matrix constructor, dim_t is the dimension of the matrix and val is
         * the initial value for the matrix elements
         */
        explicit MtmMat<T>(Dimensions dim_t, const T &val = T());

        //copy constructor
        MtmMat<T>(const MtmMat<T> &original) = default;

        //vector-to-matrix constructor
        MtmMat<T>(const MtmVec<T> &original);

        //assignment
        virtual MtmMat<T> &operator=(const MtmMat &original){
            if (this == &original) return *this;
            dimensions = original.dimensions;
            this->RootVector<T>::operator=(original);
            return *this;
        }

        //destructor (Note: MtmVec matrix will destroy its own sub-lists)
        virtual ~MtmMat<T>() = default;

        //operator[] inherited from MtmVec

        /* Matrix multiplication
         * creates new matrix
         * technique: when multiplying matrices A*B
         * row i of new matrix equals sum of scalar products of jth element of
         * row i of A * row j of B
         */
        friend MtmMat<T> operator*(const MtmMat<T> &matrix1,
            const MtmMat<T> &matrix2){

            const size_t n = matrix1.dimensions.getCol();
            if (n != matrix2.dimensions.getRow()) {
                throw MtmExceptions::DimensionMismatch(matrix1.dimensions,
                                                       matrix2.dimensions);
            }

            const size_t numRows = matrix1.dimensions.getRow();
            const size_t numCols = matrix2.dimensions.getCol();
            MtmMat<T> answer(Dimensions(numRows, numCols), defaultElement);
            for (int i = firstIndex; i < numRows; ++i) {
                for (int j = firstIndex; j < numCols; ++j) {
                    for (int runner = firstIndex; runner < n; ++runner) {
                        answer[i][j] += matrix1[i][runner] *
                                        matrix2[runner][j];
                    }
                }
            }
            return answer;
        }

        //scalar multiplication
        MtmMat<T> operator*(const T& scalar) const {
            MtmMat<T> answer(*this);
            for (size_t row = firstIndex; row<dimensions.getRow(); ++row){
                answer[row] = (*this)[row] * scalar;
            }
            return answer;
        }

        //Matrix addition
        friend MtmMat<T> operator+(const MtmMat<T> &matrix1,
                                   const MtmMat<T> &matrix2){

            const size_t numRows = matrix1.dimensions.getRow();
            const size_t numCols = matrix1.dimensions.getCol();
            if (matrix2.dimensions != matrix1.dimensions) {
                throw MtmExceptions::DimensionMismatch(matrix1.dimensions,
                                                       matrix2.dimensions);
            }

            MtmMat<T> answer(matrix1);
            for (int row = firstIndex; row < numRows; ++row) {
                answer[row] = matrix1[row] + matrix2[row]; //vector addition
            }
            return answer;
        }

        //scalar addition
        MtmMat<T> operator+(const T& scalar) const {
            MtmMat<T> answer(*this);
            for (size_t row = firstIndex; row<dimensions.getRow(); ++row){
                answer[row] = (*this)[row] + scalar;
            }
            return answer;
        }

        /*
         * Function that get function object f and uses it's () operator on each
         * element in the matrix columns.
         * It outputs a vector in the size of the matrix columns where each element
         * is the final output
         * by the function object's * operator
         */
        template<typename Func>
        MtmVec<T> matFunc(Func &f) const;

        /*
         * resizes a matrix to dimension dim, new elements gets the value val.
         */
        virtual void resize(Dimensions dim, const T &val = T());

        /*
         * reshapes matrix so linear elements value are the same without
         * changing num of elements.
         */
        virtual void reshape(Dimensions newDim);

        /*
         * Performs transpose operation on matrix
         */
        virtual void transpose();

        class iterator {
        protected:
            MtmMat *self;
            size_t linearIndex;

        private:
            friend class MtmMat<T>;

            template<T>
            friend iterator begin();

            template<T>
            friend iterator end();

            explicit iterator(MtmMat* self, size_t startIndex) :
                    self(self), linearIndex(startIndex) {}

        public:
            T& operator*() {
                return self->linearIndexToReference(linearIndex);
            }

            bool operator==(const iterator &other){
                return linearIndex == other.linearIndex;
            }

            bool operator!=(const iterator &other){
                return !((*this) == other);
            }

            virtual iterator &operator++() {
                if (*this != self->end()) ++linearIndex;
                return *this;
            }
        };

        iterator begin() {
            return iterator(this, firstIndex);
        }

        iterator end() {
            return iterator(this,
                            coordinatesToLinearIndex(dimensions.getRow() - 1,
                                                     dimensions.getCol() - 1) + 1);
        }

        class nonzero_iterator : public iterator {
        public:
            nonzero_iterator &operator++() override {
                if (*this != this->self->end()) iterator::operator++(*this);
                while (*this != this->self->end() && this->operator*() == 0) {
                    iterator::operator++(*this);
                }
            }

        private:
            nonzero_iterator nzbegin();

            nonzero_iterator nzend();

            explicit nonzero_iterator(MtmMat *self, size_t startIndex) :
                    iterator(self, startIndex) {
                if (this->operator*() == 0) ++(*this);
            }
        };

        nonzero_iterator nzbegin() {
            return nonzero_iterator(this, firstIndex);
        }

        nonzero_iterator nzend() {
            iterator temp = end();
            return nonzero_iterator(this, temp.linearIndex);
        }

        Dimensions getDimensions() const{
            return dimensions;
        }

        void printMatrix(){
            cout << endl;
            const MtmMat<T>& mat = (*this);
            for (int i = 0; i <mat.getDimensions().getRow(); i++){
                for (int j = 0; j<mat.getDimensions().getCol(); j++){
                    cout << mat[i][j] << " ";
                }
                cout << endl;
            }
        } //DEBUG
    };

    //matrix-vector multiplication (promotion on hold because class is generic)
    template<typename T>
    MtmMat<T> operator*(const MtmVec<T> &vector1, const MtmVec<T> &vector2);

    template<typename T>
    MtmMat<T> operator*(const MtmVec<T> &vector1, const MtmMat<T> &matrix2)  {
        return MtmMat<T>(vector1) * matrix2;
    }

    template<typename T> //TODO: flip coordinates thrown in case of error (in all commutativeness functions like this)
    MtmMat<T> operator*(const MtmMat<T> &matrix1, const MtmVec<T> &vector2)  {
        return MtmMat<T>(vector2) * matrix1;
    }

    template <typename T>
    MtmMat<T> operator*(const T &scalar, const MtmMat<T> &matrix)  {
        return matrix*scalar;
    }

    //matrix-vector addition (promotion on hold because class is generic)
    template <typename T>
    MtmMat<T> operator+(const MtmVec<T>& vector1, const MtmMat<T>& matrix2){
        return MtmMat<T>(vector1)+matrix2;
    }

    template <typename T>
    MtmMat<T> operator+(const MtmMat<T>& matrix1, const MtmVec<T>& vector2) {
        return vector2+matrix1;
    }

    template <typename T>
    MtmMat<T> operator+(const T &scalar, const MtmMat<T> &matrix)  {
        return matrix+scalar;
    }



    //implementations begin here
    template<typename T>
    MtmMat<T>::MtmMat(Dimensions dim_t, const T &val) : dimensions(dim_t),
                                                        RootVector<T>(dim_t.getRow(), MtmVec<T>(dim_t.getCol(), val)) {

        //make sub-vectors horizontal
        for (int row = firstIndex; row < dim_t.getRow(); ++row) {
            (*this)[row].transpose(); //MtmVec transpose
        }
    }

    template <typename T>
    void swap(T& a, T& b){
        T temp = a;
        a = b;
        b = temp;
    }

    template<typename T>
    MtmMat<T>::MtmMat(const MtmVec<T> &original) : MtmMat<T>(Dimensions(
            (original.getIsColumn()? 1 : original.getSize()),
            (original.getIsColumn()? original.getSize() : 1)),0) {
        const size_t rows = dimensions.getRow(), cols = dimensions.getCol();
        for (int i = firstIndex; i < ((original.getIsColumn()) ? cols :
                                      rows); ++i) {
            if (original.getIsColumn()) {
                (*this)[firstIndex][i] = original[i];
            } else {
                (*this)[i][firstIndex] = original[i];
            }
        }
    }

    template<typename T>
    MtmMat<T> operator*(const MtmVec<T> &vector1, const MtmVec<T> &vector2) {
        return MtmMat<T>(vector1) * MtmMat<T>(vector2);
    }

    template<typename T>
    MtmMat<T> MtmMat<T>::exchangeMatrix(size_t dimension) const {
        MtmMat<T> answer(Dimensions(dimension, dimension), defaultElement);
        for (int i = firstIndex; i < dimension; i++) {
            answer[dimension - i - 1][i] = 1;
        }
        return answer;
    }

    template<typename T>
    void MtmMat<T>::transpose() {
        //POTENTIAL ERROR: replacing matrix (rather than modifying) may reset
        // additional properties of subclasses

        const size_t rows = dimensions.getRow(), cols = dimensions.getCol(),
                maxDim = (rows > cols) ? rows : cols;
        resize(Dimensions(maxDim, maxDim), defaultElement); //make matrix square

        for (int i=firstIndex; i<maxDim; ++i){
            for (int j=i; j<maxDim; j++){
                swap((*this)[i][j], (*this)[j][i]);
            }
        }

        resize(Dimensions(cols, rows), errorValue);
    }

    template<typename T>
    void MtmMat<T>::reshape(Dimensions newDim) {
        //POTENTIAL ERROR: replacing matrix (rather than modifying) may reset
        // additional properties of subclasses

        const size_t numRows = newDim.getRow(), numCols = newDim.getCol();
        const size_t maxIndex =
                MtmMat<T>::coordinatesToLinearIndex(numRows - 1, numCols - 1);
        //sanitize inputs
        cout << numRows << " " << numCols << " " << maxIndex << endl;
        if ((newDim.getRow() <= firstIndex || newDim.getCol() <= firstIndex)
            ||(dimensions.getCol() * dimensions.getRow() != numRows * numCols)){

            throw MtmExceptions::ChangeMatFail(dimensions, newDim);
        }

        MtmMat<T> replacement(newDim, defaultElement);
        for (size_t index = firstIndex; index < maxIndex; ++index) {
            replacement.linearIndexToReference(index) =
                    linearIndexToReference(index);
        }

        *this = replacement;
    }

    template<typename T>
    void MtmMat<T>::resize(Dimensions dim, const T &val) {
        //POTENTIAL ERROR: replacing matrix (rather than modifying) may reset
        // additional properties of subclasses
        if (dimensions.getRow() <= firstIndex || dimensions.getCol() <= firstIndex) {
            throw MtmExceptions::ChangeMatFail(dimensions, dim);
        }

        MtmMat<T> replacement(dim, val);
        for (int row = firstIndex; row < dimensions.getRow() && row < dim.getRow(); ++row) {
            for (int col = firstIndex; col < dimensions.getCol() && col < dim.getCol(); ++col) {
                replacement[row][col] = (*this)[row][col];
                assert(replacement[row][col] != errorValue); //DEBUG
            }
        }
        (*this) = replacement;
    }

    template<typename T>
    template<typename Func>
    MtmVec<T> MtmMat<T>::matFunc(Func &f) const {
        //create vector to hold result
        MtmVec<T> answer(dimensions.getCol(), defaultElement);

        //create a vector of each column and call f on it
        const size_t numCols = dimensions.getCol();
        for (size_t col = firstIndex; col < numCols; ++col) {
            answer[col] = getColAsVector(col).vecFunc(f);
        };

        answer.transpose();
        return answer;
    }

    template<typename T>
    const MtmVec<T> MtmMat<T>::getColAsVector(size_t col) const {
        MtmVec<T> answer(dimensions.getRow(), defaultElement);
        for (int row = firstIndex; row < dimensions.getRow(); ++row) {
            const T& cell = (*this)[row][col];
            answer[row] = cell;
        }
        return answer;
    }
}


#endif //EX3_MTMMAT_H
