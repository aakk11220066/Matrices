#include <iostream>

#include "MtmVec.h"
#include "MtmMat.h"
#include "MtmMatSq.h"
#include "MtmMatTriag.h"
#include "Complex.h"
#include <vector>


#include <assert.h>
using namespace MtmMath;
using std::cout;
using std::endl;

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
        assert (res_nz_it[i]==(*it));
        ++i;
    }
}

using std::ostream;

ostream& operator<<(ostream& outStream, const Complex& num){
    outStream<<num.re<<"+"<<num.im<<"i";
    return outStream;
}

void warmUp(){
    exceptionsTest();
    iterators();
    constructors();
    dataTypes();
    FuncExample();
    cout<<"warmUp Done! \n";
}
void triangle_test(){
    try {
        MtmMatTriag<Complex> m(4,1, false);
        m[4][3]=5;
        assert(false);
    }
    catch (MtmExceptions::AccessIllegalElement& e){
        cout<< e.what() <<endl;
    }
    try {
        MtmMatTriag<Complex> m(4,1, false);
        m[0][3]=5;
        assert(false);
    }
    catch (MtmExceptions::AccessIllegalElement& e){
        cout<< e.what() <<endl;
    }
    MtmMatTriag<Complex> m(4,1, false);
    m[2][1]=5;
    assert(m[2][1]==5 and m[2][0]==1);
    m.resize(Dimensions(5,5),6);
    assert(m[4][0]==6 and m[4][1]==6 and m[4][2]==6 and m[3][0]==1);
    try {
        m[0][4]=5;
        assert(false);
    }
    catch (MtmExceptions::AccessIllegalElement& e){
        cout<< e.what() <<endl;
    }
    m.transpose();
    try {
        m[4][0]=5;
        assert(false);
    }
    catch (MtmExceptions::AccessIllegalElement& e){
        cout<< e.what() <<endl;
    }
    try {
        m[-1][0]=5;
        assert(false);
    }
    catch (MtmExceptions::AccessIllegalElement& e){
        cout<< e.what() <<endl;
    }
    assert(m[0][4]==6 and m[1][4]==6 and m[2][4]==6 and m[0][3]==1);
    m.transpose();
    MtmMatTriag<Complex> m2(5,1);
    MtmMat<Complex> res(m+m2);
    assert(res[0][4]==1 and res[4][1]==6 and res[4][2]==6 and res[0][3]==1);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j <5 ; ++j) {
            res[i][j];
        }
    }
    MtmMatTriag<Complex> mt(1,3);
    MtmVec<Complex> v(1,3);
    MtmMat<Complex> res2(mt+v);
    assert(res2[0][0]==6);
    MtmMat<Complex> res3(m+0);
    for (int k = 0; k < 5; ++k) {
        for (int i = 0; i < 5; ++i) {
            res3[k][i];
        }
    }
    MtmMat<Complex> res4=mt*v;
    assert(res4[0][0]==9);
    MtmMat<Complex> res5=m*5;
    //assert(res5[0][0]==5 and res5[2][1]==25 and res5[4][3]==30 and res5[0][2]==0); //TODO: Roi : I believe they were supposed to access these elements with const access, here they're not using const access
    const MtmMat<Complex> const_res5 = res5;
    assert(const_res5[0][0]==5 and const_res5[2][1]==25 and const_res5[4][3]==30 and const_res5[0][2]==0);
    MtmMatTriag<Complex> newM(res5);
    assert(res5[0][0]==5 and res5[2][1]==25 and res5[4][3]==30);
    try {
        newM[0][2]=5;
        assert(false);
    }
    catch (MtmExceptions::AccessIllegalElement& e){
        cout<< e.what() <<endl;
    }
    cout<<"res5 = "<<res5<<endl;
    try {
        res5[0][1]=1;
        MtmMatTriag<Complex> new2(res5);
        assert(false);
    }
    catch (MtmExceptions::IllegalInitialization& e){
        cout<< e.what() <<endl;
    }
    res5.transpose();
    res5[1][0]=0;
    MtmMatTriag<Complex> new2(res5);
    assert(new2[0][2]!=0);
    try {
        new2[2][0]=5;
        assert(false);
    }
    catch (MtmExceptions::AccessIllegalElement& e){
        cout<< e.what() <<endl;
    }
    const MtmMatTriag<Complex> const_tria(new2);
    for (int l = 0; l < 5; ++l) {
        for (int i = 0; i < 5; ++i) {
            const_tria[l][i];
        }
    }
    try {
        MtmMatTriag<Complex> matTriag(3,1);
        matTriag.reshape(Dimensions(3,3));
    }
    catch (MtmExceptions::ChangeMatFail& e){
        cout<< e.what() <<endl;
    }
    MtmMatTriag<int> mat1(3,1, false);
    MtmMatTriag<int> mat2(3,1);
    mat1.resize(Dimensions(4,4),0);
    mat2.resize(Dimensions(4,4),2);
    MtmMat<int> res7(mat1*mat2);
    res7.transpose();
    res7=res7+1;
    res7=1+res7;
    res7=res7*2;
    res7=-res7;
    res7=1-res7;
    res7=2*res7;
    res7=res7-mat2;
    MtmMatTriag<int> m9(4,8,false);
    m9.transpose();
    m9[3][3]=0;
    res7=res7-m9;
    std::vector<int> vec={5,14,14,18,5,9,18,26,5,9,13,34,0,0,0,8};
    std::vector<int> nz_vec={5,14,14,18,5,9,18,26,5,9,13,34,8};
    int k=0;
    for(int i :res7){
        assert(i==vec[k]);
        k++;
    }
    assert(k==vec.size());
    k=0;
    for (MtmMat<int>::nonzero_iterator n = res7.nzbegin(); n != res7.nzend(); ++n) {
        assert(*n==nz_vec[k]);
        k++;
    }
    assert(k==nz_vec.size());


    cout<<"triangle test Done!! \n";

}

void vector_test(){
    MtmVec<int> v(5,1);
    try {
        v.resize(Dimensions(1, 5), 2);
        assert(false);
    }
    catch (MtmExceptions::ChangeMatFail& e){
        cout<< e.what() <<endl;
    }
    MtmVec<int> v1(5,2);
    MtmVec<int> res=v1-v;
    res.resize(Dimensions(8,1),0);
    res.resize(Dimensions(10,1),1);
    for (int i = 0; i < 5; ++i) {
        assert(res[i]==1);
    }
    int i=0;
    for(MtmVec<int>::nonzero_iterator nz_it=res.nzbegin();nz_it!=res.nzend();++nz_it,i++){
        assert(*nz_it==1);
    }
    assert(i==7);
    try {
        MtmVec<int> v10(1,0);
        MtmVec<int> v11(1,0);
        v11.transpose();
        v10+v11;
        assert(false);
    }
    catch (MtmExceptions::DimensionMismatch& e){
        cout<< e.what() <<endl;
    }
    MtmVec<int> v2(4,2);
    MtmVec<int> v3(6,4);
    v2.resize(Dimensions(7,1),4);
    v3.transpose();
    v3.resize(Dimensions(1,7),4);
    for (int j = 0; j < 6; ++j) {
        v2[j]=v3[j];
    }
    try {
        v2+v3;
        assert(false);
    }
    catch (MtmExceptions::DimensionMismatch& e){
        cout<< e.what() <<endl;
    }
    MtmVec<int>::iterator it1=v2.begin();
    int k = 0;
    for (MtmVec<int>::iterator it2=v3.begin(); it2!=v3.end(); ++it1,++it2,++k) {
        assert(*it1==*it2);
        v3=1+v3;
        assert(v3[0]==5);
        v3=v3-1;
        assert(v3[0]==4);
        v3=v3*(-1);
        assert(v3[0]==-4);
        v3=(-1)*v3;
        assert(v3[0]==4);
        v3=1-v3;
        assert(v3[0]==-3);
        v3=v3-1;
        assert(v3[0]==-4);
        v3=-v3;
    }
    assert(k==7);
    v3.resize(Dimensions(1,10),0);
    v3.resize(Dimensions(1,20),4);
    v3.resize(Dimensions(1,22),0);
    k=0;
    for (MtmVec<int>::nonzero_iterator it2=v3.nzbegin(); it2!=v3.nzend(); ++it2,++k) {
        assert(*it2==4);
    }
    assert(k==17);
    cout<<"vector test done!\n";
}

int main() {
    warmUp();
    triangle_test();
    vector_test();
    cout<<"tests done!\n";
    return 0;

}

