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
        MtmVec<int> v1(4, 8);
        v1[1]=7;
        //v1.print_vec();
        int num;
        int i = 0;
        /*
        MtmVec<int>::iterator it =v1.begin();
        for (;i<= 3;++it) {
            num = *it;
            i++;
            cout << num << " " << i << endl;
        }
        if (it == v1.end()) cout << "alright" << endl;
         */
        //MtmVec<int> v2(4,7);
        //MtmVec<int> v4 (2,6);
        //MtmVec<int> v3 = v1 + v4;

        int a =3, b=3, count = 0;
        MtmMat<int> m1(Dimensions(a,b),5);

        for (int j = 0; j <b; j++){
            for (int i = 0; i<a; i++){
                m1[i][j] = count;
                count++;
            }
            cout << endl;
        }

        //m1.reshape(Dimensions (2,8));
        MtmMat<int>::iterator it = m1.begin();
        for (int i = 0; i < 20; i++, ++it) {
            //num = *it;
            //cout << num << "hello" << endl;
        }
        if (it==m1.end()) cout <<"works" << endl;
        if (++(m1.end())==(m1.end())) cout << "works2" << endl;
        m1.printMatrix();
        cout << endl;
        cout << endl;
        cout << endl;
        //MtmVec<int> v2 = 3* v1;

        //MtmMat<int> m2 = v1*m1;
        //m2.printMatrix();
        //MtmMat<int> m2 = m1;
        // MtmMat<int> m3 = m1*m2;
/*
        cout << endl;
        for (int i = 0; i <b; i++){
            for (int j = 0; j<a; j++){
                cout << m2[i][j] << " ";
            }
            cout << endl;
        }
*/
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
