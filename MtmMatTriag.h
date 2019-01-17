#ifndef EX3_MTMMATTRIAG_H
#define EX3_MTMMATTRIAG_H


#include <vector>
#include "MtmExceptions.h"
#include "Auxilaries.h"
#include "MtmMatSq.h"

using std::size_t;

namespace MtmMath {
    enum TriangleType {
        NEITHER, UPPER, LOWER
    };

    template<typename T>
    class MtmMatTriag : public MtmMatSq<T> {
        /* checks if given MtmMatSq is also upper/lower triangular or neither.
         * By default selects upper
         */
        TriangleType isUpperOrLower(const MtmMatSq<T> &mat) const;

        //takes square matrix and zeroes out top/bottom triangle
        //also handles locking
        void triangulate(MtmMatSq<T>& target, bool makeUpper);

        bool upper = true;

    public:

        /*
         * Triangular Matrix constructor, m is the number of rows and columns in the matrix,
         * val is the initial value for the matrix elements and isUpper_ is whether it is upper
         * Rectangular matrix (true means it is)
         */
        MtmMatTriag<T>(size_t m, const T &val = T(), bool isUpper_t = true);

        //copy constructor
        MtmMatTriag<T>(const MtmMatTriag &original) = default;

        //constructor for normal/square matrices to triangular
        MtmMatTriag<T>(const MtmMatSq<T> &original);

        //destructor
        virtual ~MtmMatTriag<T>() = default;

        //operator=
        virtual MtmMatTriag<T> &operator=(const MtmMatTriag &) = default;

        //override: note that triangle is now opposite kind of triangle
        void transpose() override {
            upper = !upper;
            MtmMatSq<T>::transpose();
        }

        //override: triangulate after resizing
        void resize(Dimensions dim, const T &val = T()) override {
            MtmMatSq<T>::resize(dim, val);
            triangulate(*this, upper);
        }
    };

    //implementation begins here
    using namespace MtmMath;

    template<typename T>
    MtmMatTriag<T>::MtmMatTriag(size_t m, const T &val, bool isUpper_t) :
            MtmMatSq<T>(m, val) {

        triangulate(*this, isUpper_t);
        upper = isUpper_t;
    }

    template<typename T>
    MtmMatTriag<T>::MtmMatTriag(const MtmMatSq<T> &original) :
            MtmMatSq<T>(original) {

        TriangleType triangleType = isUpperOrLower(original);
        if (triangleType == NEITHER) throw MtmExceptions::IllegalInitialization();
        upper = (triangleType == UPPER);

        triangulate(*this, upper);
    }

    template<typename T>
    TriangleType MtmMatTriag<T>::isUpperOrLower(const MtmMatSq<T> &mat) const {
        size_t m = mat.dimensions.getRow();
        bool isUpper = true, isLower = true;
        for (int row = 0; row < m; ++row) {
            for (int col = 0; col < m; ++col) {
                //found nonzero beneath (next line: above) diagonal
                if (row > col && mat[row][col] != 0) isUpper = false;
                if (row < col && mat[row][col] != 0) isLower = false;
            }
        }
        TriangleType pre_answer = (isLower) ? LOWER : NEITHER;
        return (isUpper) ? UPPER : pre_answer;
    }

    template<typename T>
    void MtmMatTriag<T>::triangulate(MtmMatSq<T>& target, bool makeUpper) {
        const size_t m = target.getDimensions().getRow();
        for (size_t row = 0; row < m; ++row) {
            target[row].setLock(false); //unlock row
            for (int col = 0; col < m; ++col) {
                if ((row > col && makeUpper) || (row < col && !makeUpper)) {
                    target[row][col] = 0;
                }
            }
            //lock
            if (makeUpper) {
                target[row].setLockStartIndex(row);
            }
            else {
                target[row].setLockEndIndex(row);
            }
            target[row].setLock(true); //lock row
        }
    }
}
#endif //EX3_MTMMATTRIAG_H
