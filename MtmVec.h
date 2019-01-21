#ifndef EX3_MTMVEC_H
#define EX3_MTMVEC_H

#include <vector>
#include "MtmExceptions.h"
#include "Auxilaries.h"
#include "complex.h"

using std::size_t;

namespace MtmMath {
    const size_t defaultElement = 0, firstIndex = 0, errorValue = 8998;

    template<typename T>
    class MtmVec {
    protected:
        size_t size;
    private:
        MtmVec *self;
        T *data;
        bool locked=false;
    public:
        bool isLocked() const;

    private:
        size_t lockStartIndex = 0;
        size_t lockEndIndex = size-1;
    protected:
        bool is_column = true;
    public:
        /*
         * Vector constructor, m is the number of elements in it and val is the
         * initial value for the matrix elements
         */
        explicit MtmVec(size_t m = 1, const T &val = T());

        //destructor
        virtual ~MtmVec() noexcept;

        //copy constructor
        MtmVec<T>(const MtmVec<T> &original) : size(original.size) {
            try {
                data = new T[size];
            } catch (std::bad_alloc) {
                throw MtmExceptions::OutOfMemory();
            }

            is_column = original.is_column;
            for (size_t i = firstIndex; i < original.size; ++i) {
                (*this)[i] = original[i];
            }

            locked = original.locked;
            lockStartIndex = original.lockStartIndex;
            lockEndIndex = original.lockEndIndex;
        }

        //assignment operator
        MtmVec<T> &operator=(const MtmVec<T>& original) {
            if (this == &original) return *this;
            delete[] data;
            size = original.size;
            try {
                data = new T[size];
            } catch (std::bad_alloc) {
                throw MtmExceptions::OutOfMemory();
            }
            is_column = original.is_column;
            locked = false;
            for (size_t i = firstIndex; i < original.size; ++i) {
                (*this)[i] = original[i];
            }

            locked = original.locked;
            lockStartIndex = original.lockStartIndex;
            lockEndIndex = original.lockEndIndex;

            return *this;
        }


        //Scalar multiplication (returns new vector)
        virtual MtmVec<T> operator*(const T &scalar) const;

        //Scalar addition
        virtual MtmVec<T> operator+(const T &scalar) const {
            MtmVec<T> answer(*this);
            answer.setLock(false);
            for (size_t i = firstIndex; i < size; ++i) {
                answer[i] = (*this)[i] + scalar;
            }
            return answer;
        }

        //Vector addition
        virtual MtmVec<T> operator+(const MtmVec<T> other) const;

        //negation
        MtmVec<T> operator-() const {
            return (-1) * (*this);
        }

        /*
         * Function that get function object f and uses it's () operator on each
         * element in the vectors.
         * It outputs the function object's * operator after iterating on all
         * the vector's elements
         */
        template<typename Func>
        T vecFunc(Func &f) const {
            Func g = f;
            for (size_t i = firstIndex; i < size; i++) {
                g((*this)[i]);
            }
            return *g;
        }

        /*
         * Resizes a vector to dimension dim, new elements gets the value val.
         * Notice vector cannot transpose through this method.
         */
        virtual void resize(Dimensions dim, const T &val = T()) {
            size_t cols = 1, rows = 1;
            if ((!is_column && dim.getRow() != 1)||(is_column&&dim.getCol()!=1)
                ||dim.getCol()<=0||dim.getRow()<=0){
                if (is_column) {
                    rows = size;
                } else {
                    cols = size;
                }
                throw MtmExceptions::ChangeMatFail(Dimensions(cols, rows), dim);
            }

            size_t newSize = (dim.getRow() > dim.getCol()) ?
                          dim.getRow() : dim.getCol();
            T *newData = NULL;
            try {
                newData = new T[newSize];
            } catch (std::bad_alloc) {
                throw MtmExceptions::OutOfMemory();
            }
            size_t i;
            for (i = 0; i < newSize && i < size; ++i) {
                newData[i] = data[i];
            }
            for (; i< newSize; i++){
                newData[i] = val;
            }
            delete[] data;
            data = newData;
            size = newSize;
        }

        /*
         * Performs transpose operation on matrix
         */
        virtual void transpose() {
            is_column = !is_column;
        }

        virtual T &operator[](size_t index) {
            if (index < 0 || index >= size
                || (locked && (lockStartIndex>index || lockEndIndex<index))) {

                throw MtmExceptions::AccessIllegalElement();
            }
            return data[index];
        }

        virtual const T &operator[](size_t index) const {
            if (index < 0 || index >= size) {
                throw MtmExceptions::AccessIllegalElement();
            }
            return data[index];
        }

        bool setLock(bool status) {
            return locked = status;
        }

        size_t setLockStartIndex(size_t index) {
            return lockStartIndex = index;
        }

        size_t setLockEndIndex(size_t index) {
            return lockEndIndex = index;
        }

        class iterator;

        iterator begin();

        iterator end();

        class nonzero_iterator : public iterator {
        public:
            nonzero_iterator &operator++() override {
                if (*this != this->self->end()) this->iterator::operator++();
                while (*this != this->self->end() && this->operator*() == 0) {
                    this->iterator::operator++();
                }
                return *this;
            }

            friend nonzero_iterator MtmVec<T>::nzbegin();

            friend nonzero_iterator MtmVec<T>::nzend();

        private:
            explicit nonzero_iterator(MtmVec *self, size_t startIndex) :
                    iterator(self, startIndex) {
                if (startIndex!=self->end().getIndex() && this->operator*()==0){
                    ++(*this);
                }
            }

            nonzero_iterator(const MtmVec<T>::iterator& original) :
                nonzero_iterator(original.getSelf(), original.getIndex()){}
        };

        nonzero_iterator nzbegin() {
            return nonzero_iterator(begin());
        }

        nonzero_iterator nzend() {
            return nonzero_iterator(end());
        }

        bool getIsColumn() const {
            return is_column;
        }

        size_t getSize() const {
            return size;
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
        MtmVec<T> answer(*this);
        answer.setLock(false);
        for (size_t i = firstIndex; i<size;++i) answer[i]=(*this)[i] + other[i];
        return answer;
    }

    template<typename T>  //REVERSAL - treated
    MtmVec<T> operator+(const T &scalar, const MtmVec<T> vector) {
        try {
            return vector + scalar;
        } catch (MtmExceptions::MtmExceptions& e){
            e.reverseDescription();
            throwError(e);
        }
        return vector; //THIS LINE SHOULD NEVER BE REACHED!!!
    }

    template<typename T>
    MtmVec<T> MtmVec<T>::operator*(const T &scalar) const {
        MtmVec<T> answer(*this);
        bool wasLocked = locked;
        answer.setLock(false);
        for (size_t i = firstIndex; i < size; ++i) {
            answer[i] = answer[i] * scalar;
        }
        answer.setLock(wasLocked);
        return answer;
    }

    template<typename T> //REVERSAL - treated
    MtmVec<T> operator*(const T &scalar, const MtmVec<T> vector) {
        try {
            return vector * scalar;
        } catch (MtmExceptions::MtmExceptions& e){
            e.reverseDescription();
            throwError(e);
        }
        return vector; //THIS LINE SHOULD NEVER BE REACHED!!!
    }

    template <typename T>
    MtmVec<T> operator-(const T& scalar, const MtmVec<T> vector) {
        return scalar + -vector;
    }

    template <typename T>
    MtmVec<T> operator-(const MtmVec<T> vector, const T& scalar){
        return vector + -scalar;
    }

    template <typename T>
    MtmVec<T> operator-(const MtmVec<T> vector1, const MtmVec<T> vector2){
        return vector1 + -vector2;
    }

    template<typename T>
    MtmVec<T>::MtmVec(size_t m, const T &val) : size(m) {
        if (m <= 0) throw MtmExceptions::IllegalInitialization();
        try {
            data = new T[m];
        } catch (std::bad_alloc) {
            throw MtmExceptions::OutOfMemory();
        }

        for (size_t i = firstIndex; i < size; i++) {
            data[i] = val;
        }
    }

    template<typename T>
    MtmVec<T>::~MtmVec() noexcept {
        delete[] data;
    }

    template <typename T>
    class MtmVec<T>::iterator{
    protected:
        size_t index=0;
        MtmVec* self = nullptr;

        iterator(MtmVec<T>* self, size_t startIndex=0) :
                index(startIndex), self(self) {};
        friend iterator MtmVec::begin();

        friend iterator MtmVec::end();

    public:
        size_t getIndex() const{
            return index;
        }

        MtmVec* getSelf() const {
            return self;
        }

        T& operator*(){
            return (*self)[index];
        }
        const T& operator*() const{
            return (*self)[index];
        }
        bool operator==(iterator iterator2){
            return index == iterator2.index;
        }
        bool operator!=(iterator iterator2){
            return !((*this)==iterator2);
        }
        virtual iterator& operator++() {
            if (*this != self->end()) {
                ++index;
            }
            return *this;
        }
    };

    template <typename T>
    typename MtmVec<T>::iterator MtmVec<T>::begin(){
        return iterator(this);
    }

    template <typename T>
    typename MtmVec<T>::iterator MtmVec<T>::end(){
        return iterator(this, size);
    }

    template<typename T>
    bool MtmVec<T>::isLocked() const {
        return locked;
    }
}

#endif //EX3_MTMVEC_H
