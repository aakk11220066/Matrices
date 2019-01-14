#ifndef EX3_MTMVEC_H
#define EX3_MTMVEC_H

#include <vector>
#include "MtmExceptions.h"
#include "Auxilaries.h"
#include "complex.h"

//DEBUG
#include <iostream>
#include <string>

using std::size_t;
#define firstIndex 0 //DEBUG
#define defaultElement 0 //DEBUG

namespace MtmMath {

    template<typename T>
    class MtmVec {
    public:
        T *data; //FIXME
        bool locked=false;
        size_t lockStartIndex = 0;
        size_t lockEndIndex = size-1; //FIXME
    protected:
        bool is_column = true;
        size_t size;
    public:
        /*
         * Vector constructor, m is the number of elements in it and val is the initial value for the matrix elements
         */
        explicit MtmVec(size_t m = 1, const T &val = T());

        //destructor
        virtual ~MtmVec();

        //copy constructor
        MtmVec<T>(const MtmVec<T> &original) : size(original.size) {
            try {
                data = new T[size];
            } catch (std::bad_alloc) {
                throw MtmExceptions::OutOfMemory();
            }

            is_column = original.is_column;
            for (int i = firstIndex; i < original.size; ++i) {
                (*this)[i] = original[i];
            }

            locked = original.locked;
            lockStartIndex = original.lockStartIndex;
            lockEndIndex = original.lockEndIndex;
        }

        //assignment operator
        MtmVec<T> &operator=(const MtmVec<T> original) {
            if (this == &original) return *this;
            delete[] data;
            size = original.size;
            try {
                data = new T[size];
            } catch (std::bad_alloc) {
                throw MtmExceptions::OutOfMemory();
            }
            is_column = original.is_column;
            for (int i = firstIndex; i < original.size; ++i) {
                (*this)[i] = original[i];
            }
            return *this;
        }

        //Scalar multiplication (returns new vector)
        virtual MtmVec<T> operator*(const T &scalar) const;

        //Scalar addition
        virtual MtmVec<T> operator+(const T &scalar) const {
            MtmVec<T> answer(size, defaultElement);
            for (int i = firstIndex; i < size; ++i) {
                answer[i] = (*this)[i] + scalar;
            }
            return answer;
        }

        //Vector addition
        virtual MtmVec<T> operator+(const MtmVec<T> other) const;

        //negation
        MtmVec<T> operator-() {
            return (-1) * (*this);
        }

        /*
         * Function that get function object f and uses it's () operator on each element in the vectors.
         * It outputs the function object's * operator after iterating on all the vector's elements
         */
        template<typename Func>
        T vecFunc(Func &f) const {
            for (int i = firstIndex; i < size; i++) {
                f((*this)[i]);
            }
            return *f;
        }

        /*
         * Resizes a vector to dimension dim, new elements gets the value val.
         * Notice vector cannot transpose through this method.
         */
        virtual void resize(Dimensions dim, const T &val = T()) {
            int cols = 1, rows = 1;
            if ((is_column && dim.getRow() != 1)
                || (!is_column && dim.getCol() != 1)
                || dim.getCol() <= 0
                || dim.getRow() <= 0) {
                if (is_column) {
                    rows = size;
                } else {
                    cols = size;
                }
                throw MtmExceptions::ChangeMatFail(Dimensions(cols, rows), dim);
            }

            int newSize = (dim.getRow() > dim.getCol()) ?
                          dim.getRow() : dim.getCol();
            T *newData = NULL;
            try {
                newData = new T[newSize];
            } catch (std::bad_alloc) {
                throw MtmExceptions::OutOfMemory();
            }
            for (int i = 0; i < newSize && i < size; ++i) {
                newData = data;
            }
            delete[] data;
            data = newData;
        }

        /*
         * Performs transpose operation on matrix
         */
        virtual void transpose() {
            is_column = !is_column;
        }

        virtual T &operator[](int index) {
            if (index < 0 || index >= size || (locked &&
            lockStartIndex>index && lockEndIndex<index)) {
                throw MtmExceptions::AccessIllegalElement();
            }
            return data[index];
        }

        virtual const T &operator[](int index) const {//code duplicated!
            if (index < 0 || index >= size) {
                throw MtmExceptions::AccessIllegalElement();
            }
            return data[index];
        }

        bool setLock(bool status) {
            return locked = status;
        }

        int setLockStartIndex(int index) {
            return lockStartIndex = index;
        }

        int setLockEndIndex(int index) {
            return lockEndIndex = index;
        }

        class iterator{
            int index=0;
            MtmVec* self = NULL;

            iterator(MtmVec<T>* self, size_t startIndex=0) : self(self),
                                                             index(startIndex){};

            template <T>
            friend iterator begin();

            template <T>
            friend iterator end();
        public:
            T& operator*(){
                return (*this)[index];
            }
            friend bool operator==(iterator& iterator1, iterator& iterator2){
                return iterator1.index == iterator2.index;
            }
            friend bool operator!=(iterator& iterator1, iterator& iterator2){
                return !(iterator1==iterator2);
            }
            iterator& operator++() {
                if (*this != this->end()) {
                    return (*this)[index];
                }
            }
        };

        iterator begin(){
            return iterator(this);
        }

        iterator end(){
            return iterator(this, size);
        }

        //DEBUG
        friend std::ostream& operator<<(std::ostream& outstream, MtmVec<T>& me){
            std::string outstr("[");
            for (const T& elem : me) outstr+=to_string(elem)+=", ";
            return outstream << (outstr+=std::string("]"));
        }
    };

    //implementations begin here
    using namespace MtmMath;

    template<typename T>
    MtmVec<T> MtmVec<T>::operator+(const MtmVec<T> other) const {
        if (size != other.size || is_column != other.is_column) {
            size_t cols = 1, rows = 1;
            if (is_column) {
                rows = size;
            } else {
                cols = size;
            }
            size_t otherRows = 1, otherCols = 1;
            if (other.is_column) {
                otherRows = other.size;
            } else {
                otherCols = other.size;
            }
            throw MtmExceptions::DimensionMismatch(Dimensions(rows, cols),
                                                   Dimensions(otherRows, otherCols));
        }
        MtmVec<T> answer(size, defaultElement);
        for (int i = firstIndex; i < size; ++i) answer[i] = (*this)[i] + other[i];
        return answer;
    }

    template<typename T>
    MtmVec<T> operator+(const T &scalar, const MtmVec<T> vector) {
        return vector + scalar;
    }

    template<typename T>
    MtmVec<T> MtmVec<T>::operator*(const T &scalar) const {
        MtmVec<T> answer(*this);
        for (int i = firstIndex; i < size; ++i) {
            answer[i] *= scalar;
        }
        return answer;
    }

    template<typename T>
    MtmVec<T> operator*(const T &scalar, const MtmVec<T> vector) {
        return vector * scalar;
    }

    template<typename T>
    MtmVec<T>::MtmVec(size_t m, const T &val) : size(m) {
        if (m <= 0) throw MtmExceptions::IllegalInitialization();
        try {
            data = new T[m];
        } catch (std::bad_alloc) {
            throw MtmExceptions::OutOfMemory();
        }

        for (int i = firstIndex; i < size; i++) {
            data[i] = val;
        }
    }

    template<typename T>
    MtmVec<T>::~MtmVec() {
        delete[] data;
    }
}

#undef firstIndex //DEBUG
#undef defaultElement //DEBUG

#endif //EX3_MTMVEC_H
