#include <iostream>

#include "MtmVec.h"
#include "MtmMat.h"
#include "MtmMatSq.h"
#include "MtmMatTriag.h"
#include "Complex.h"

#include <assert.h>
using namespace MtmMath;
using std::cout;
using std::endl;



bool operator ==(MtmMat<int> &m1, MtmMat<int>& m2){

    if(m1.getDimensions().getRow()!=m2.getDimensions().getRow() || m1.getDimensions().getCol()!=m2.getDimensions().getCol())
        return false;

    for(size_t i=0; i<m1.getDimensions().getRow(); i++) {
        for (size_t j = 0; j < m2.getDimensions().getCol(); j++) {

            if(m1[i][j]!=m2[i][j])
                return false;

        }
    }

    return true;

}

void exceptionsTest() {
    try {
        MtmVec<float> v(0,5);
        assert(false);
    }
    catch (MtmExceptions::IllegalInitialization& e){
        cout<< e.what() <<endl;
    }

    try {
        MtmMat<int> m1(Dimensions(10,100000000000),5);
        assert(false);
    }
    catch (MtmExceptions::OutOfMemory& e){
        cout<< e.what() <<endl;
    }

    try {
        MtmVec<int> v1(3,5);
        MtmMat<int> m1(Dimensions(3,3),5);
        MtmMat<int> m3=m1+v1;
        assert(false);
    }
    catch (MtmExceptions::DimensionMismatch& e){
        cout<< e.what() <<endl;
    }
    try {
        MtmMat<int> m1(Dimensions(3,3),5);
        m1.reshape(Dimensions(2,3));
        assert(false);
    }
    catch (MtmExceptions::ChangeMatFail& e){
        cout<< e.what() <<endl;
    }
    try {
        MtmMat<int> m1(Dimensions(3,3),5);
        m1[4][3]=5;
        assert(false);
    }
    catch (MtmExceptions::AccessIllegalElement& e){
        cout<< e.what() <<endl;
    }
    try{

        MtmMatTriag<int>m(5,0);
        m[1][0]=1;

    }
    catch (MtmExceptions::AccessIllegalElement& e){
        cout<< e.what() <<endl;
    }

    try{

        MtmMatTriag<int>m(5,0,false);
        m[0][1]=1;

    }
    catch (MtmExceptions::AccessIllegalElement& e){
        cout<< e.what() <<endl;
    }

    try{

        MtmMatSq<int>m(5,0);
        m[0][1]=1;
        m[1][0]=1;

        MtmMatTriag<int>m2(m);

    }
    catch (MtmExceptions::IllegalInitialization& e){
        cout<< e.what() <<endl;
    }

    try{

        MtmMatTriag<int>m(5,0);
        m[1][0]=1;

    }
    catch (MtmExceptions::AccessIllegalElement& e){
        cout<< e.what() <<endl;
    }

    try{

        MtmMatSq<int>m(0,0);


    }
    catch (MtmExceptions::IllegalInitialization& e){
        cout<< e.what() <<endl;
    }

    try {
        MtmMat<int> m1(Dimensions(0,3),5);
        assert(false);
    }
    catch (MtmExceptions::IllegalInitialization& e){
        cout<< e.what() <<endl;
    }

    try {
        MtmMat<int> m1(Dimensions(3,0),5);
        assert(false);
    }
    catch (MtmExceptions::IllegalInitialization& e){
        cout<< e.what() <<endl;
    }
}


void constructors() {
    MtmVec<int> v1(5,3);
    MtmMat<int> m1(Dimensions(3,3),0);
    MtmMatSq<int> m2(3,1);
    MtmMatTriag<int> m3(3,1,true);
}

void dataTypes() {
    MtmVec<int> v1(5,3);
    MtmVec<double > v2(5,3);
    MtmVec<float> v3(5,3);
    MtmVec<Complex> v4(5,Complex(3,4));
}

void testOperator(){

    MtmMat<int> m1(Dimensions(2, 5), 0);
    MtmMat<int> m2(Dimensions(2, 5), 1);
    m2[1][3] = 4;
    m2[0][2] = 9;
    m1[0][0] = 1;
    m1 = m1 + m2;
    m2[0][0] += 1;
    assert(m1 == m2);
    m1 = m1 + 2;
    assert(m1[0][0] == 4);
    assert(m1[1][3] == 6);
    m1 = m1 - 2;
    assert(m1 == m2);
    m1 = 2 + m1;
    assert(m1[0][0] == 4);
    assert(m1[1][3] == 6);
    m1 = -2 + m1;
    assert(m1 == m2);
    MtmMat<int> m3(Dimensions(4, 1), 1);
    m3[1][0] = 2;
    m3[2][0] = 3;
    m3[3][0] = 4;
    //cout<<"m3 before reshape";
    //m3.printMatrix();
    MtmMat<int> m4(m3);
    m3.reshape(Dimensions(1, 4));
    //cout<<"m3 after reshape, before transpose";
    //m3.printMatrix();
    m3.transpose();
    assert(m3 == m4);
    m3.reshape(Dimensions(4, 1));
    //cout<<endl<<"m3 after second reshape";
    //m3.printMatrix();
    assert(m3 == m4);
    MtmMat<int> m5(Dimensions(6, 2), 1);
    MtmMat<int> m6(m5);
    //cout<<endl<<"m6 = ";
    //m6.printMatrix();
    m6.transpose();
    //cout<<endl<<"m6 after transpose = ";
    //m6.printMatrix();
    //cout<<endl<<"m5 = ";
    //m5.printMatrix();
    m5[5][1] = 3;
    m5[5][0] = 2;
    m5[4][0] = 7;
    //cout<<endl<<"m5 after changes= ";
    //m5.printMatrix();
    m5.reshape(Dimensions(2,6));
    //cout<<endl<<"m5 after reshape = ";
    //m5.printMatrix();
    m6[0][2]=7;m6[1][2]=2; m6[1][5]=3;
    //printm(m6);
    //printm(m5);
    assert(m5 == m6);
    m6.transpose();
    m5.reshape(Dimensions(6, 2));
    MtmVec<int> v1(3, 1);
    MtmMat<int> v2(Dimensions(1, 3),1);
    MtmMat<int> m7(Dimensions(3, 3), 1);
    cout<<endl<<"m7 = ";
    m7.printMatrix();
    cout<<endl<<"v1 = "<<endl;
    v1.print_vec();
    cout<<endl<<"v2 = "<<endl;
    v2.printMatrix();
    m6 = v1 * v2;
    assert(m6 == m7);
    m6 = v2 * v1;
    MtmMat<int> m8(Dimensions(1, 1), 3);
    assert(m6 == m8);
    MtmMat<int> m9(Dimensions(2, 3), 1);
    m9[0][0] = 1;
    m9[1][0] = 1;
    m9[0][1] = 2;
    m9[1][1] = 2;
    m9[0][2] = 3;
    m9[1][2] = 3;
    v1[1] = 2;
    v1[2] = 3;
    MtmMat<int> m10(Dimensions(2, 1), 14);
    m4 = m9 * v1;
    assert(m4 == m10);
    MtmMat<int> m11(m9);
    m11.transpose();
    m6 = m9 * m11;
    MtmMat<int> res(Dimensions(2, 2), 14);
    assert(m6 == res);
    MtmMat<int> m12(Dimensions(3, 6), 1);
    m8 = m7 * m12;
    MtmMat<int> res2(Dimensions(3, 6), 3);

}

void FuncExample() {
    class maxAbsolute {
        int currMax;
    public:
        maxAbsolute() : currMax(0) {}
        void operator()(int x) {
            int absX = x>(-x) ? x : -x;
            if (currMax<absX) {currMax=absX;}
        }
        int operator*() { return  currMax;}
    };

    MtmVec<int> v(5,0);
    v[1]=3;v[2]=-7;v[3]=-1;v[4]=2;
    MtmMat<int> m(Dimensions(2,3),0);
    m[0][1]=1;m[1][0]=2;
    m[0][1]=3;m[1][1]=2;
    m[0][2]=5;m[1][2]=-6;
    maxAbsolute f;
    assert (v.vecFunc(f)==7);
    MtmVec<int> res(m.matFunc(f));
    assert(res[0]==2 and res[1]==3 and res[2]==6);
}

void iterators() {
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

        int x = *it;
        assert (res_nz_it[i]==(x));
        ++i;
    }

    i=0;

    MtmMatTriag<int> m1(5,1,false);

    int res[]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    for(MtmMatTriag<int>::nonzero_iterator it=m1.nzbegin();it!=m1.end();++it) {

        int t1 = res[i];
        int t2 = *it;
        *it=15;
        assert(t1 == t2 && 15>i);
        i++;
        //5
    }

    MtmMatTriag<int> m2(5,0,false);
    m2[4][4]=2;
    MtmMatTriag<int>::nonzero_iterator it=m2.nzbegin();
    int t = *it;
    assert(t==2);
    ++it;
    assert(it==m2.nzend());


}



int main() {
    constructors();
    dataTypes();
    FuncExample();
    iterators();
    testOperator();
    exceptionsTest();


    return 0;

}
