#ifndef MATRICES_ROOTVECTOR_H
#define MATRICES_ROOTVECTOR_H

#include <vector>
#include "MtmExceptions.h"
#include "Auxilaries.h"
#include "complex.h"


#include "MtmVec.h"
#include <vector>
using std::size_t;

namespace MtmMath {
    template <typename T> //describes a vector of vectors
    class RootVector {
        size_t size;
        std::vector<MtmVec<T>> vectors;
    public:
        //ctor
        RootVector(size_t columns, const MtmVec<T>& val);

        //copy ctor
        RootVector(const RootVector<T>& original);

        //dtor
        virtual ~RootVector() = default;

        //operator=
        virtual RootVector<T>& operator=(const RootVector<T>& original)=default;

        //operator[]
        virtual MtmVec<T>& operator[](size_t index);

        //const operator[]
        virtual const MtmVec<T>& operator[](size_t index) const;

        bool isLocked(){
            if (size == 0) return false;
            return (*this)[0].isLocked();
        }

        void setLock(bool newStatus){
            for (size_t i=0; i<size; i++){
                (*this)[i].MtmVec<T>::setLock(newStatus);
            }
        }
    };

    template <typename T>
    RootVector<T>::RootVector(size_t columns, const MtmVec<T>& val) try:
        size(columns), vectors(size, val) {
            if (size<=0) throw MtmExceptions::IllegalInitialization();
        }
        catch(std::bad_alloc){
            throw MtmExceptions::OutOfMemory();
        }

    template <typename T>
    RootVector<T>::RootVector(const MtmMath::RootVector<T> &original) :
        RootVector(original.size, MtmVec<T>()){

        for (size_t i=0; i<size; ++i){
            (*this)[i] = original[i];
        }
    }

    template <typename T>
    const MtmVec<T>& RootVector<T>::operator[](size_t index) const {
        if (index >= size) throw MtmExceptions::AccessIllegalElement();
        return vectors.at(index);
    }

    template <typename T>
    MtmVec<T>& RootVector<T>::operator[](size_t index) {
        if (index >= size) throw MtmExceptions::AccessIllegalElement();
        return vectors.at(index);
    }
}
#endif //MATRICES_ROOTVECTOR_H
