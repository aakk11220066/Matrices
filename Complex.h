#ifndef EX3_COMPLEX_H
#define EX3_COMPLEX_H


namespace MtmMath {
    class Complex {
        double re, im;
    public:
        Complex(double r = 0, double i = 0) : //Roi : I removed "virtual" keyword, constructor should accept automatic promotion (e.g. putting 0 into matrix of complex numbers)
                re(r), im(i) {}

        virtual ~Complex() = default; //Roi : all destructors should be virtual (prevents bugs when making polymorphisms)

        Complex(const Complex &) = default;

        Complex &operator=(const Complex &) = default;

        Complex &operator+=(const Complex &c);

        Complex &operator*=(const Complex &c);

        Complex &operator-=(const Complex &c);

        Complex operator-() const;

        bool operator==(const Complex &c) const;

    };

    MtmMath::Complex operator+(const Complex &a, const Complex &b);

    MtmMath::Complex operator*(const Complex &a, const Complex &b);

    MtmMath::Complex operator-(const Complex &a, const Complex &b);


    MtmMath::Complex &Complex::operator+=(const Complex &c) {
        re += c.re;
        im += c.im;
        return *this;
    }

    MtmMath::Complex &Complex::operator*=(const Complex &c) {
        re *= c.re;
        im *= c.im;
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
