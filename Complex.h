#ifndef EX3_COMPLEX_H
#define EX3_COMPLEX_H

//DEBUG
#include <string>

namespace MtmMath {
    class Complex {
        double re, im;
    public:
        Complex(double r = 0, double i = 0) :
                re(r), im(i) {}

        virtual ~Complex() = default;

        Complex(const Complex &) = default;

        Complex &operator=(const Complex &) = default;

        Complex &operator+=(const Complex &c);

        Complex &operator*=(const Complex &c);

        Complex &operator-=(const Complex &c);

        Complex operator-() const;

        bool operator==(const Complex &c) const;

        bool operator!=(const Complex& c) const{
            return !((*this) == c);
        } //Roi: added this

        std::string to_string(){ //DEBUG
            const int accuracy = 2;
            return std::to_string(re).substr(0,accuracy)+=std::string("+")+=std::to_string(im).substr(0,accuracy)+=std::string("i");
        }
    };

    std::ostream& operator<<(std::ostream& stream, Complex& c){
        return stream<<c.to_string();
    }//DEBUG

    MtmMath::Complex operator+(const Complex &a, const Complex &b);

    MtmMath::Complex operator*(const Complex &a, const Complex &b);

    MtmMath::Complex operator-(const Complex &a, const Complex &b);


    MtmMath::Complex &Complex::operator+=(const Complex &c) {
        re += c.re;
        im += c.im;
        return *this;
    }

    MtmMath::Complex &Complex::operator*=(const Complex &c) {
        auto oldReal = re;
        re = (re*c.re) - (im*c.im);
        im = (oldReal*c.im) + (im*c.re);
        return *this;
    }

    MtmMath::Complex &Complex::operator-=(const Complex &c) {
        return this->operator+=(-c); // or *this += -c
    }

    MtmMath::Complex Complex::operator-() const {
        return Complex(-re, -im);
    }

    MtmMath::Complex operator+(const Complex &a, const Complex &b) {
        Complex c = a;
        return c += b;
    }

    MtmMath::Complex operator*(const Complex &a, const Complex &b) {
        Complex c = a;
        return c *= b;
    }

    bool MtmMath::Complex::operator==(const Complex &c) const {
        return c.re == re && c.im == im;
    }

}
#endif //EX3_COMPLEX_H
