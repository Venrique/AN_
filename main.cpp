#include <math.h>
#include <iostream>
#include <complex>
#include <iomanip>


using std::cout;	using std::endl;
using std::complex;	using std::pow;
using std::sqrt;	using std::abs;
using std::exp;		using std::cos;
using std::real;	using std::setprecision;
using std::imag;

//Colocar aquÃ­ la funciÃ³n
inline complex<double> f(complex<double> x){
    return pow(x,3)+(3.0*pow(x,2))-1.0;
}

void muller(){
    //------------------ Valores iniciales ------------------
    complex<double> p0(1.0, 0);
    complex<double> p1(2.0, 0);
    complex<double> p2(3.0, 0);
    //------------------ Valores iniciales ------------------

    complex<double> h1 = p1 - p0;
    complex<double> h2 = p2 - p1;
    complex<double> delta1 = (f(p1) - f(p0))/h1;
    complex<double> delta2 = (f(p2) - f(p1))/h2;
    complex<double> d = (delta2 - delta1)/(h2 + h1);
    complex<double> b(0, 0);
    complex<double> D(0, 0);
    complex<double> E(0, 0);
    complex<double> h(0, 0);
    complex<double> p(0, 0);

    for(int i = 1; i <= 1000; i++){
        b = delta2 + h2*d;
        D = sqrt( pow(b, 2.0) -4.0*d*f(p2) );

        if( abs(b - D) < abs(b + D) ) E = b + D;
        else E = b - D;

        h = (-2.0*f(p2))/E;
        p = p2 + h;

        cout << setprecision(8) << i << ": " << "p0: " << p0 << " p1: " << p1 << " p2: " <<  p2 << " p: " << p << " Eabs: " << abs(h) << endl;

        if(abs(h) < abs(pow(10, -5.0))){
            cout << "Luego de " << i <<" iteraciones." << endl;
            break;
        }

        p0 = p1;
        p1 = p2;
        p2 = p;
        h1 = p1 - p0;
        h2 = p2 - p1;
        delta1 = (f(p1) - f(p0))/h1;
        delta2 = (f(p2) - f(p1))/h2;
        d = (delta2 - delta1)/(h2 + h1);
    }

    complex<double> mycomplex(1.0, 5.0);
    cout << f(mycomplex) << endl;
}

void lagrange(){
    int n,i,j;
    float mult,sum=0,x[10],f[10],a;

    x[0] = 1950;
    x[1] = 1960;
    x[2] = 1970;
    x[3] = 1980;
    x[4] = 1990;
    x[5] = 2000;

    f[0] = 151326;
    f[1] = 179323;
    f[2] = 203302;
    f[3] = 226542;
    f[4] = 249633;
    f[5] = 281422;

    //IMPORTANTESSSSSSSSSS
    n = 6;

    a = 2020;

    for(i=0;i<=n-1;i++)
    {
        mult=1;
        for(j=0;j<=n-1;j++)
        {
            if(j!=i)
                mult*=(a-x[j])/(x[i]-x[j]);
        }
        sum+=mult*f[i];
    }
    cout<<"\nThe estimated value of f(x) = "<<sum;
}

void diferencias(){

    double datos[4][5] = {{0.6,-0.17694460},{0.7,0.01375227},{0.8,0.22363362},{1.0,0.65809197}};
    int filas = (sizeof(datos)/ sizeof(datos[0]));
    int columnas = (sizeof(datos[0])/ sizeof(datos[0][0]));

    for (int i = 0; i < filas; ++i) {
        for (int j = 2; j < columnas; ++j) {
            if(j<=i+1){
                datos[i][j] = (datos[i][j-1]-datos[i-1][j-1])/(datos[i][0]-datos[i-j+1][0]);
            }
        }
    }

    for (int k = 0; k < filas; ++k) {
        for (int i = 0; i < columnas; ++i) {
            cout<<datos[k][i]<<"    ";
        }
        cout<<endl;
    }

}


int main(){
    //muller();

    //lagrange();

    diferencias();

    return 0;
}