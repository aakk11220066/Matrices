#ifndef EX3_MTMEXCEPTIONS_H
#define EX3_MTMEXCEPTIONS_H

#include <exception>
#include <string>
#include <iostream>
#include "Auxilaries.h"

using std::to_string;

namespace MtmMath {
    namespace MtmExceptions {
        class MtmExceptions : public std::exception {
        protected:
            std::string description;
        public:
            MtmExceptions(std::string description) : description(description){}
            const char* what() const noexcept {
                return description.c_str();
            }
            virtual ~MtmExceptions() throw() {}
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
                            += (to_string(mat1Dimensions.getRow()) += ",")
                            += to_string(mat2Dimensions.getCol()) += ")"
                            ) {}
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
                    += (to_string(oldMatDimensions.getCol()) += ") (")
                    += (to_string(newMatDimensions.getRow()) += ",")
                    += (to_string(newMatDimensions.getCol()) += ")")
                    ){}
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
    }
}


#endif //EX3_MTMEXCEPTIONS_H
