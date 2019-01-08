#ifndef EX3_MTMMATREC_H
#define EX3_MTMMATREC_H

#include <vector>
#include "MtmExceptions.h"
#include "Auxilaries.h"
#include "MtmMat.h"

using std::size_t;

namespace MtmMath {

    template <typename T>
    class MtmMatSq : public MtmMat<T> {
    public:

        /*
         * Rectangular Matrix constructor, m is the number of rows and columns in the matrix
         * and val is the initial value for the matrix elements
         */
        MtmMatSq (size_t m, const T& val=T()):MtmMat<T>(Dimensions(m,m), val){}

        //copy constructor
        MtmMatSq (const MtmMatSq& original) = default;

        //constructor for building square matrix based on matrix that happens to have square dimensions
        MtmMatSq (const MtmMat<T>& original) : MtmMat<T>(original){
            if (original.dimensions().getRow()
                != original.dimensions().getCol()) {
                throw MtmExceptions::IllegalInitialization();
            }
        }

        //destructor
        virtual ~MtmMatSq() = default;

        //operator=
        virtual MtmMatSq& operator=(const MtmMatSq& original) = default;

        //override: prohibit resizing to non-square
        void resize(Dimensions newDim, const T& val=T()) override {
            if (newDim.getRow() != newDim.getCol()){
                throw MtmExceptions::ChangeMatFail(this->dimensions, newDim);
            }
            MtmMat<T>::resize(newDim, val);
        }

            //override: prohibit calling reshape on MtmMatSq
        void reshape(Dimensions newDim) override {
            throw MtmExceptions::ChangeMatFail(this->dimensions, newDim);
        }
    };

}

#endif //EX3_MTMMATREC_H
