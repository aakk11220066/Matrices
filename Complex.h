#ifndef EX3_COMPLEX_H
#define EX3_COMPLEX_H
namespace MtmMath {
    class Complex {
        double re,im;
    public:
        Complex(double r=0, double i=0);
        Complex(const Complex&) ;
        Complex& operator=(const Complex&);

    };
}
#endif //EX3_COMPLEX_H
