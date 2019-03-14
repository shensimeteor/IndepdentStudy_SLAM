#include <iostream>
using namespace std;



class A{
public: 
    A() {} 
//    A(const A& a) {
//        cout<<"Copy Constructor of A"<<endl;
//    }
    int x=5;
};

class B: public A{
public:
    B() {}
    B(const A& a): A(a) {
        cout << "Copy Constructor of B"<<endl;
    }
    double x=4.4444;
};



int main(){
    B b;
    cout<<b.x<<endl;
    A a;
    A a2(a);
    cout<<a2.x<<endl;
    A a3=a;
    cout<<a3.x<<endl;
    B b2(a);
}
