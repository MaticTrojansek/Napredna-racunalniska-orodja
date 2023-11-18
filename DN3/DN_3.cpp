#include <iostream>
#include <cmath>

// Taylorjev razvoj funkcije arctan(x/2)
double atan_taylor_series(double x) {
    
    double result = 0.0;
    double term = x;
    double x_squared = (x * x);
    double divisor = 1.0;

    for (int i = 1; i <= 100; ++i) {  // Uporabim prvih 100 èlenov Taylorjeve vrste
        
        // Èe je naslednji èlen deljiv z 2 (vsak drugi èlen), potem odšteje èlen
        // V nasprotnem primeru ga prišteje
        if (i % 2 == 0) {
            result -= term / divisor; 
        }
        else {
            result += term / divisor;
        }
        term *= x_squared;
        divisor += 2.0;
    }

    return result;
}

// Integralska funkcija
double f(double x) {
    return exp(3 * x) * atan_taylor_series(x / 2);
}

// Funkcija trapezne metode
double trapezoidalMethod(double a, double b, int n) {
    double h = (b - a) / n;  // Korak h
    double integral = (f(a) + f(b)) / 2; // Integral

    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        integral += f(x);
    }

    return integral * h;
}

int main() {
    // Meji obmoèja
    double a = 0.0;
    double b = 3.1415926535 / 4.0;

    // Število intervalof
    int n;

    std::cout << "Vpisi stevilo intervalov: ";
    std::cin >> n; // Vhodno vrednost vpišemo(št.intervalof)

    // Izraèun integrala z uporabo trapezne funkcije
    double result = trapezoidalMethod(a, b, n);

    // Izpiši rezutat
    std::cout << "Aproksimirana vrednost integrala: " << result;

    return 0;
}