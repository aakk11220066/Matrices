#ifndef EX3_MTMVEC_H
#define EX3_MTMVEC_H

#include <vector>
#include "MtmExceptions.h"
#include "Auxilaries.h"
#include "complex.h"

using std::size_t;

namespace MtmMath {
    template <typename T>
    class MtmVec {
    public:
        /*
         * Vector constructor, m is the number of elements in it and val is the initial value for the matrix elements
         */
        MtmVec(size_t m, const T& val=T());


        //Roi : please implement (or else I will).  Should be 1 line, use vecFunc
        //Scalar multiplication (returns new vector)
        MtmVec<T> operator*(const T& scalar);

        //Roi : don't forget to implement unary operator (they said we need this)
        MtmVec<T> operator-();

        /*
         * Function that get function object f and uses it's () operator on each element in the vectors.
         * It outputs the function object's * operator after iterating on all the vector's elements
         */
        template <typename Func>
        T vecFunc(Func& f) const;

        /*
         * Resizes a vector to dimension dim, new elements gets the value val.
         * Notice vector cannot transpose through this method.
         */
        virtual void resize(Dimensions dim, const T& val=T());

        /*
         * Performs transpose operation on matrix
         */
        virtual void transpose();

        //Roi
        virtual T& operator[](int index);
    };
}

#endif //EX3_MTMVEC_H
