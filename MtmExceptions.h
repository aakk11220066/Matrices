#ifndef EX3_MTMEXCEPTIONS_H
#define EX3_MTMEXCEPTIONS_H

#include <exception>
#include <string>
#include <iostream>
#include "Auxilaries.h"
#include <string.h>

using std::to_string;

namespace MtmMath {
    namespace MtmExceptions {
        enum errorType {MISMATCH, MAT_FAIL, ILLEGAL_INIT, OUT_MEM, ACCESS};

        class MtmExceptions : public std::exception {
        protected:
            std::string description;
        public:
            errorType error;

            MtmExceptions(std::string description, errorType error)
                : description(description), error(error){}
            const char* what() const noexcept {
                return description.c_str();
            }
            virtual ~MtmExceptions() throw() {}

            void reverseDescription();
        };

        /*
         * Exception for illegal initialization of an object, needs to output
         * "MtmError: Illegal initialization values" in what() class function
         */
        class IllegalInitialization : public MtmExceptions {
        public:
            IllegalInitialization()
                : MtmExceptions("MtmError: Illegal initialization values",
                        ILLEGAL_INIT){}
        };

        /*
         * Exception for Memory allocation failure for an object, needs to output
         * "MtmError: Out of memory" in what() class function
         */
        class OutOfMemory : public MtmExceptions {
        public:
            OutOfMemory()
                    : MtmExceptions("MtmError: Out of memory", OUT_MEM){}
        };

        /*
         * Exception for dimension mismatch during a mathematical function, needs to output
         * "MtmError: Dimension mismatch: (<mat 1 row>,<mat 1 col>) (<mat 2 row>,<mat 2 col>)"
         * in what() class function
         */
        class DimensionMismatch : public MtmExceptions {
        public:
            DimensionMismatch(Dimensions mat1Dimensions,
                    Dimensions mat2Dimensions) : MtmExceptions( std::string(
                            "MtmError: Dimension mismatch: (")
                            += (to_string(mat1Dimensions.getRow()) += ",")
                            += (to_string(mat1Dimensions.getCol()) += ") (")
                            += (to_string(mat2Dimensions.getRow()) += ",")
                            += to_string(mat2Dimensions.getCol()) += ")"
                            , MISMATCH) {}
        };

        /*
         * Exception for error for changing matrix/vector shape in reshape and resize, needs to output
         * "MtmError: Change matrix shape failed from: (<mat row>,<mat col>) (<new mat row>,<new mat col>)"
         * in what() class function
         */
        class ChangeMatFail : public MtmExceptions {
        public:
            ChangeMatFail(Dimensions oldMatDimensions,
                    Dimensions newMatDimensions) : MtmExceptions( std::string(
                    "MtmError: Change matrix shape failed from (")
                    += (to_string(oldMatDimensions.getRow()) += ",")
                    += (to_string(oldMatDimensions.getCol()) += ") to (")
                    += (to_string(newMatDimensions.getRow()) += ",")
                    += (to_string(newMatDimensions.getCol()) += ")")
                    , MAT_FAIL) {}
        };

        /*
         * Exception for accessing an illegal element in matrix or vector, needs to output
         * "MtmError: Attempt access to illegal element" in what() class function
         */
        class AccessIllegalElement : public MtmExceptions {
        public:
            AccessIllegalElement() : MtmExceptions(
                    "MtmError: Attempt access to illegal element", ACCESS) {}
        };

        void throwError(MtmExceptions& e){
            switch (e.error){
                case ILLEGAL_INIT:
                    throw dynamic_cast<IllegalInitialization&>(e);
                case OUT_MEM:
                    throw dynamic_cast<OutOfMemory&>(e);
                case MISMATCH:
                    throw dynamic_cast<DimensionMismatch&>(e);
                case MAT_FAIL:
                    throw dynamic_cast<ChangeMatFail&>(e);
                case ACCESS:
                    throw dynamic_cast<AccessIllegalElement&>(e);
            }
        }

        void MtmExceptions::reverseDescription() {
            std::string pieces[5];
            if (error!=MISMATCH && error!=MAT_FAIL) return;

            //split string up into pieces
            std::string toFlip = description;
            pieces[0] = toFlip.substr(0, toFlip.find("(")+1);
            toFlip = toFlip.substr(pieces[0].length());

            pieces[1] = toFlip.substr(0, toFlip.find(")"));
            toFlip = toFlip.substr(pieces[1].length());

            pieces[2] = toFlip.substr(0, toFlip.find("(")+1);
            toFlip = toFlip.substr(pieces[2].length());

            pieces[3] = toFlip.substr(0, toFlip.find(")"));

            pieces[4] = toFlip.substr(pieces[2].length());

            //reverse order of parentheses
            description = pieces[0] += pieces[3] += pieces[2]
                += pieces[1] += pieces[4];

        }
    }
}


#endif //EX3_MTMEXCEPTIONS_H
