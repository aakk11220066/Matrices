



#ifndef EX3_MTMVEC_H
#define EX3_MTMVEC_H

#include <vector>
#include "MtmExceptions.h"
#include "Auxilaries.h"
#include "complex.h"
#include "MtmVec.h"

using std::size_t;
#define firstIndex 0 //DEBUG
#define defaultElement 0 //DEBUG

namespace MtmMath {

    template <typename T>
    class MtmVec {
        const long long int DISABLED = -1;
        const long long int ENABLED = 0;
    private:
        T* data;
    protected:
        bool is_column = true;
        long long int* locked;
        size_t size;
    public:
        /*
         * Vector constructor, m is the number of elements in it and val is the initial value for the matrix elements
         */
        explicit MtmVec(size_t m, const T& val=T());

        //destructor
        virtual ~MtmVec();

        //copy constructor
        MtmVec<T>(const MtmVec<T>& original){
            MtmVec<T> answer(original.size, defaultElement);
            answer.is_column = original.is_column;
            for (int i=firstIndex; i<original.size; ++i) {
                answer[i] = original[i];
                answer.locked[i] = original.locked[i];
            }
            return answer;
        }

        //assignment operator
        MtmVec<T>& operator=(const MtmVec<T> original){
            if (this == &original) return *this;
            delete[] data;
            delete[] locked;
            size=original.size;
            data = new T[size];
            locked = new long long int[size];
            is_column = original.is_column;
            for (int i=firstIndex; i<original.size; ++i){
                (*this)[i] = original[i];
                locked[i] = original.locked[i];
            }
        }

        //Mat Locker
        virtual MtmVec<T>& lockVector(int begin_index, int end_index){
            for (int i =firstIndex; i< size; i++){
                if ((i >= begin_index)&&(i<=end_index)){
                    locked[i] = ENABLED;
                }
            }
        }

        //Scalar multiplication (returns new vector)
        virtual MtmVec<T> operator*(const T& scalar) const;

        //Scalar addition
        virtual MtmVec<T> operator+(const T& scalar) const {
            MtmVec<T> answer(size, defaultElement);
            for (int i=firstIndex; i<size; ++i){
                answer = (*this) + scalar;
            }
            return answer;
        }

        //Vector addition
        virtual MtmVec<T> operator+(const MtmVec<T> other) const;

        //negation
        MtmVec<T> operator-(){
            return (-1)*(*this);
        }

        /*
         * Function that get function object f and uses it's () operator on each element in the vectors.
         * It outputs the function object's * operator after iterating on all the vector's elements
         */
        template <typename Func>
        T vecFunc(Func& f) const{
            for (int i=firstIndex; i<size; i++){
                f((*this)[i]);
            }
            return *f;
        }

        /*
         * Resizes a vector to dimension dim, new elements gets the value val.
         * Notice vector cannot transpose through this method.
         */
        virtual void resize(Dimensions dim, const T& val=T()){
            int cols=1, rows=1;
            if ((is_column && dim.getRow()!=1)
                || (!is_column && dim.getCol()!=1)
                || dim.getCol() <= 0
                || dim.getRow() <= 0){
                if (is_column){
                    rows = size;
                }
                else{
                    cols = size;
                }
                throw MtmExceptions::ChangeMatFail(Dimensions(cols, rows), dim);
            }
            int newSize = (dim.getRow()>dim.getCol())?
                          dim.getRow() : dim.getCol();
            T* newData = new T[newSize];
            for (int i=0; i<newSize && i<size; ++i){
                newData = data;
            }
            delete[] data;
            data = newData;
        }

        /*
         * Performs transpose operation on matrix
         */
        virtual void transpose(){
            is_column = !is_column;
        }

        virtual T& operator[](int index){
            if (index<0 || index>size) {
                throw MtmExceptions::AccessIllegalElement();
            }
            return data[index];
        }
//don't we need to add throw also here??
        virtual const T& operator[](int index) const{
            return (*this)[index];
        }
    };

    //complementary symmetry (necessary because class is generic, so no
    // promotion)
    template <typename T>
    MtmVec<T> operator*(const T& scalar, const MtmVec<T> vector){
        return vector*scalar;
    }

    template <typename T>
    MtmVec<T> operator+(const T& scalar, const MtmVec<T> vector){
        return vector+scalar;
    }
}

//implementations begin here
using namespace MtmMath;

template <typename T>
MtmVec<T> MtmVec<T>::operator+(const MtmVec<T> other) const {
    if (size!=other.size || is_column!=other.is_column){
        int cols=1, rows=1;
        if (is_column) {
            rows=size;
        }
        else{
            cols = size;
        }
        size_t otherRows=1, otherCols=1;
        if (other.is_column){
            otherRows = other.size;
        }
        else{
            otherCols = other.size;
        }
        throw MtmExceptions::DimensionMismatch(Dimensions(rows, cols),
                                               Dimensions(otherRows,
                                                       otherCols));
    }
    MtmVec<T> answer(size,defaultElement);
    for (int i=firstIndex; i<size; ++i) answer[i] = (*this)[i]+other[i];
    return answer;
}
//added error in case of illegal size. we need to verify this.
template <typename T>
MtmVec<T>::MtmVec(size_t m, const T& val) : size(m) {
    try {
        if (m < 0) throw MtmExceptions::IllegalInitialization();
        locked = new long long int[m];
        data = new T[m];
    } catch (std::bad_alloc) {
        throw MtmExceptions::OutOfMemory();
    }
    for (int i=firstIndex; i<size; i++){
        locked[i] = DISABLED;
        data[i] = val;
    }
}

template <typename T>
MtmVec<T>::~MtmVec(){
    delete[] data;
}

template <typename T>
MtmVec<T> MtmVec<T>::operator*(const T& scalar) const {
    MtmVec<T> answer(*this);
    for (int i=firstIndex; i<size; ++i){
        answer[i] *= scalar;
    }
    return answer;
}

#undef firstIndex //DEBUG
#undef defaultElement //DEBUG

#endif //EX3_MTMVEC_H
