#include <iostream>

#include "MtmVec.h"
#include "MtmMat.h"
//#include "MtmMatSq.h"
//#include "MtmMatTriag.h"
#include "complex.h"

#include <assert.h>
using namespace MtmMath;
using std::cout;
using std::endl;
using std::size_t;
/*
void exceptionsTest() {
    try {
        MtmVec<float> v(0,5);
        assert(false);
    }
    catch (MtmExceptions::IllegalInitialization& e){
        cout<< e.what() <<endl;
    }

    try {
        //MtmMat<int> m1(Dimensions(10,100000000000),5);
        //assert(false);
    }
    catch (MtmExceptions::OutOfMemory& e){
        cout<< e.what() <<endl;
    }

    try {
        MtmVec<int> v1(3,5);
        //MtmMat<int> m1(Dimensions(3,3),5);
        //MtmMat<int> m3=m1+v1;
        //assert(false);
    }
    catch (MtmExceptions::DimensionMismatch& e){
        cout<< e.what() <<endl;
    }
    try {
        //MtmMat<int> m1(Dimensions(3,3),5);
        //m1.reshape(Dimensions(2,3));
        //assert(false);
    }
    catch (MtmExceptions::ChangeMatFail& e){
        cout<< e.what() <<endl;
    }
    /*try {
        MtmMat<int> m1(Dimensions(3,3),5);
        m1[4][3]=5;
        assert(false);
    }
    catch (MtmExceptions::AccessIllegalElement& e){
        cout<< e.what() <<endl;
    }*/
//}

void constructors() {
    MtmVec<int> v1(5,3);
    MtmMat<int> m1(Dimensions(3,3),0);
    //MtmMatSq<int> m2(3,1);
    //MtmMatTriag<int> m3(3,1,true);
}
/*
void dataTypes() {
    MtmVec<int> v1(5,3);
    MtmVec<double > v2(5,3);
    MtmVec<float> v3(5,3);
    //MtmVec<Complex> v4(5,Complex(3,4));
}
*/

//MtmVec<int> v(5,0);
//v[1]=3;v[2]=-7;v[3]=-1;//v[4]=2;
/*MtmMat<int> m(Dimensions(2,3),0);
m[0][1]=1;m[1][0]=2;
m[0][1]=3;m[1][1]=2;
m[0][2]=5;m[1][2]=-6;
maxAbsolute f;
assert (v.vecFunc(f)==7);
MtmVec<int> res(m.matFunc(f));
assert(res[0]==2 and res[1]==3 and res[2]==6);*/


/*void iterators() {
    MtmMatSq<int> m(2,0);
    m[1][0]=1;m[1][1]=2;
    int res_it[]={0,1,0,2};
    int res_nz_it[]={1,2};
    int i=0;
    for (MtmMatSq<int>::iterator it=m.begin();it!=m.end();++it) {
        assert (res_it[i]==(*it));
        ++i;
    }
    i=0;
    for (MtmMatSq<int>::nonzero_iterator it=m.nzbegin();it!=m.nzend();++it) {
        assert (res_nz_it[i]==(*it));
        ++i;
    }
}*/

int main() {
    try {
        /*cout<<"vector1 = "<<endl;
        MtmVec<int> v1(4, 8);
        v1.print_vec();
        cout<<"v1 + 2 = "<<endl;
        (v1 + 3).print_vec();*/
        //MtmVec<int> v2(4,7);
        //MtmVec<int> v4 (2,6);

        //MtmMat<Complex> m1(Dimensions(4,9), Complex(14,16));
        MtmMat<int> m1(Dimensions(3,4), 2);
        m1[0][0] = m1[1][1] = m1[2][2] = 8;
        cout << "Matrix 1:";
        m1.printMatrix();
        cout << endl;
        //m1 = Complex(3,-2) * m1;

        MtmMat<int> m2 = 9*(m1*3 + -5);
        m2[0][0] = m2[1][1] = m2[2][2] = 40;
        m2[2][3] = 13;
        m2[1][0] = 20;
        m2[0][1] = 15;
        m2.transpose();
        cout << "Matrix 2:";
        m2.printMatrix();
        cout << endl;
       //MtmMat<int> m3 = m1*m2;

       MtmMat<int> m4 = m1+m2;

       cout << "Matrix1 + 9(3Matrix1 - 5) =";
       m4.printMatrix();

        //v2[2] = 3;
        //v2[3] = 2;
        //v1 = v2;
        //v1.operator=(v2);
        //v1.print_vec();
        //MtmVec<int> v3 = v1 +v2+4;
        //int g = *v2.data;
        //int h = int(v2.lockEndIndex);
        //v2.transpose();
        //v2.print_vec();
        //v2.resize(Dimensions(9,1), 6);
        //v2.print_vec();
        //v3.print_vec();

        //cout << g << endl;
        //cout << h << endl;
    }
    catch (MtmExceptions::IllegalInitialization& e){
        cout<< e.what() <<endl;
    }
    catch (MtmExceptions::DimensionMismatch& e){
        cout<< e.what() <<endl;
    }
}
