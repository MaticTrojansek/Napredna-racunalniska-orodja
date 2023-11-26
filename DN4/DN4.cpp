#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

//Funkcija diferenčne sheme naprej
double DiferenčnaShemaNaprej(std::vector<double> x, std::vector<double> fx, int i) {
    double h = x[i + 1] - x[i];
    double odvod1 = (-3 * fx[i] + 4 * fx[i + 1] - fx[i + 2]) / (2 * h);

    return odvod1;
}
//Funkcija centralne diferenčne sheme
double CentralnaDiferenčnaShema(std::vector<double> x, std::vector<double> fx, int i) {
    double h = x[i] - x[i - 1];
    double odvod2 = (fx[i + 1] - fx[i - 1]) / (2 * h);

    return odvod2;
}
//Funkcija diferenčne sheme nazaj
double DiferenčnaShemaNazaj(std::vector<double> x, std::vector<double> fx, int i) {
    double h = x[i] - x[i - 1];
    double odvod3 = (3 * fx[i] - 4 * fx[i - 1] + fx[i - 2]) / (2 * h);

    return odvod3;
}

int main() {

    //Definiram vektorje ID, x, f(x)
    std::vector<double> xv;
    std::vector<double> fxv;
    std::vector<int> idv;

    //Začnem brati datoteko poly.txt
    std::ifstream file("poly.txt");
    std::string v1;

    if (file.is_open()) {

        //Preberem število točk iz prve vrstice in shranim kot integer
        std::getline(file, v1);
        int st_tock = std::stoi(v1);

        //Berem vrstico po vrstico in si shranjujem ID,x in f(x)
        for (int i = 0; i < st_tock; i++) {
            std::string vrstica;
            std::getline(file, vrstica);

            //Uporabim stringstream za razčlenitev podatkov
            std::istringstream s1(vrstica);

            std::string id_str;
            std::getline(s1, id_str, ':');//Berem do :
            int ID = std::stoi(id_str);
            idv.push_back(ID);

            //Preostanek vrstice (od : dalje)
            std::string preostanek_vrstice;
            std::getline(s1, preostanek_vrstice);

            //Uporabim stringstream za branje preostanka vrstice
            std::istringstream preostanek_ss(preostanek_vrstice);

            //Preberem vrednosti x in f(x)
            double x, fx;
            preostanek_ss >> x >> fx;

            //Shranim vrednosti v vektorja x in f(x)
            xv.push_back(x);
            fxv.push_back(fx);

            //Izpišem na konzolo vrednosti ID,x in f(x), da preverim, če je datoteka ustrezno prebrana
            //std::cout << "ID: " << idv[i] << " x: " << xv[i] << " f(x): " << fxv[i] << std::endl;

        }
        file.close();
    }
    else {
        std::cerr << "Datoteke ni bilo mogoče odpreti: poly.txt" << std::endl;
    }

    //Ustvarim datoteko diff_poly.txt
    std::ofstream outputfile("diff_poly.txt");
    if (outputfile.is_open()) {
        //Prva vrstica...
        outputfile << "Odvodi:" << std::endl;
        //Odvod 1. točke
        outputfile << "f'(x=" << xv[0] << ") = " << DiferenčnaShemaNaprej(xv, fxv, 0) << std::endl;
        //Vmesni odvodi
        for (int i = 1; i < xv.size() - 1; i++) {
            outputfile << "f'(x=" << xv[i] << ") = " << CentralnaDiferenčnaShema(xv, fxv, i) << std::endl;
        }
        //Odvod zadnje točke
        outputfile << "f'(x=" << xv[xv.size()-1] << ") = " << DiferenčnaShemaNazaj(xv, fxv, xv.size() - 1);

        //Na konzolo izpišem, uspešen zapis datoteke
        std::cout << "Datoteka diff_poly.txt je bila uspesno ustvarjena in se nahaja v vasi mapi projekta!" << std::endl;
        std::cout << "V njej so shranjeni vsi odvodi.";
        outputfile.close();
    }
    else {
        std::cerr << "Datoteke ni bilo mogoče ustvariti: diff_poly.txt" << std::endl;
    }
    return 0;
}