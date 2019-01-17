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
        class MtmExceptions : public std::exception {
        protected:
            std::string description;
            enum errorType {MISMATCH, MAT_FAIL, OTHER};
            errorType error = OTHER;
        public:
            MtmExceptions(std::string description) : description(description){}
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
                : MtmExceptions("MtmError: Illegal initialization values"){}
        };

        /*
         * Exception for Memory allocation failure for an object, needs to output
         * "MtmError: Out of memory" in what() class function
         */
        class OutOfMemory : public MtmExceptions {
        public:
            OutOfMemory()
                    : MtmExceptions("MtmError: Out of memory"){}
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
                            ) {
                error = MISMATCH;
            }
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
                    "MtmError: Change matrix shape failed from: (")
                    += (to_string(oldMatDimensions.getRow()) += ",")
                    += (to_string(oldMatDimensions.getCol()) += ") to (")
                    += (to_string(newMatDimensions.getRow()) += ",")
                    += (to_string(newMatDimensions.getCol()) += ")")
                    ) {
                error = MAT_FAIL;
            }
        };

        /*
         * Exception for accessing an illegal element in matrix or vector, needs to output
         * "MtmError: Attempt access to illegal element" in what() class function
         */
        class AccessIllegalElement : public MtmExceptions {
        public:
            AccessIllegalElement() : MtmExceptions(
                    "MtmError: Attempt access to illegal element") {}
        };

        void MtmExceptions::reverseDescription() {
            std::string pieces[5];
            if (error==OTHER) return;

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
