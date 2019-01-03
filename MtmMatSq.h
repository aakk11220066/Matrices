#ifndef EX3_MTMMATREC_H
#define EX3_MTMMATREC_H

#include <vector>
#include "MtmExceptions.h"
#include "Auxilaries.h"
#include "MtmMat.h"

using std::size_t;
using namespace MtmMath;

namespace MtmMath {

    template <typename T>
    class MtmMatSq : public MtmMat<T> {
    public:

        /*
         * Rectangular Matrix constructor, m is the number of rows and columns in the matrix
         * and val is the initial value for the matrix elements
         */
        MtmMatSq (size_t m, const T& val=T()) : MtmMat<T>(Dimensions(m,m), val) {}

        //constructor for building square matrix based on matrix that happens to have square dimensions
        MtmMatSq (const MtmMat<T>& original) : MtmMat<T>(original){
            if (original.dimensions().getRow()
                != original.dimensions().getCol()) {
                throw MtmExceptions::IllegalInitialization(); //FIXME: return error info
            }
        }

        //copy constructor
        MtmMatSq (const MtmMatSq& original) = default;

        //destructor
        ~MtmMatSq() = default;

        //operator=
        MtmMatSq& operator=(const MtmMatSq& original) = default;

        //throw ChangeMatFail if called reshape on MtmMatSq
        void reshape(Dimensions) override {
            throw MtmExceptions::ChangeMatFail(); //FIXME: return error info
        }
    };

}

#endif //EX3_MTMMATREC_H
