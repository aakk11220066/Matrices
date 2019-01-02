#ifndef EX3_MTMEXCEPTIONS_H
#define EX3_MTMEXCEPTIONS_H

#include <exception>
#include <string>
#include <iostream>
#include "Auxilaries.h"

namespace MtmMath {
    namespace MtmExceptions {
        class MtmExceptions : public std::exception {
        protected:
            std::string description("MY IDIOT PROGRAMMER IS OUT TO LUNCH.");
        public:
            virtual ~MtmExceptions() throw() {}
        };

        /*
         * Exception for illegal initialization of an object, needs to output
         * "MtmError: Illegal initialization values" in what() class function
         */
        class IllegalInitialization : public MtmExceptions {
                //TODO: add constructor which sets description
        };

        /*
         * Exception for Memory allocation failure for an object, needs to output
         * "MtmError: Out of memory" in what() class function
         */
        class OutOfMemory : public MtmExceptions {
            //TODO: add constructor which sets description
        };

        /*
         * Exception for dimension mismatch during a mathematical function, needs to output
         * "MtmError: Dimension mismatch: (<mat 1 row>,<mat 1 col>) (<mat 2 row>,<mat 2 col>)"
         * in what() class function
         */
        class DimensionMismatch : public MtmExceptions {
            //TODO: add constructor which sets description

        };

        /*
         * Exception for error for changing matrix/vector shape in reshape and resize, needs to output
         * "MtmError: Dimension mismatch: (<mat row>,<mat col>) (<new mat row>,<new mat col>)"
         * in what() class function
         */
        class ChangeMatFail : public MtmExceptions {
            //TODO: add constructor which sets description

        };

        /*
         * Exception for accessing an illegal element in matrix or vector, needs to output
         * "MtmError: Illegal initialization values" in what() class function
         */
        class AccessIllegalElement : public MtmExceptions {
            //TODO: add constructor which sets description
        };
    }
}


#endif //EX3_MTMEXCEPTIONS_H
