#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <matplot/matplot.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <complex>

// domyślny namespace
using namespace std;

// aliasy namespace'ów do użytku
namespace py = pybind11;
namespace mp = matplot;

// zdefiniowanie liczby PI
const double PI = acos(-1.0);

// DEKLARACJE FUNKCJI
vector<double> gen_sig(int, double, double, double, double, int, double);
void rys(vector<double>, vector<double>, string);
void test_gen_sig(int, double, double, double, double, int, double);
void low_f_filter_real(vector<double>&, const vector<double>&, double);
void low_f_filter_im(vector<complex<double>>&, int);
vector<complex<double>> dft_vector(vector<double>);
vector<double> amplitude(vector<complex<double>>, int);
vector<double> frequency(vector<complex<double>>, int);

vector<double> generator(double, double, double, double, double, double, double, double, double, int);
void dft_test1(double, double, double, double, double, double, double, double, double, int);
void dft_test3(int, double, double, double, double, double, double, double, double, double, double);
void dft_test5(int, int, double, double, double, double, double, double, double, double, double);

vector<double> gen_sig(int rodzaj_funkcji, double amplituda, double czestotliwosc, double przesuniecie_faz, double ruch_y, int liczba_probek, double dlugosc) {
    vector<double> x = mp::linspace(0, dlugosc * PI, liczba_probek);
    vector<double> y(liczba_probek);

    switch (rodzaj_funkcji) {
        case 1: {
            y = mp::transform(x, [&](auto x) { return (ruch_y + amplituda * sin(czestotliwosc * x - przesuniecie_faz * PI)); });
            break;
        }
        case 2: {
            y = mp::transform(x, [&](auto x) { return (ruch_y + amplituda * cos(czestotliwosc * x - przesuniecie_faz * PI)); });
            break;
        }
        case 3: {
            y = mp::transform(x, [&](auto x) { return (sin(czestotliwosc * x - przesuniecie_faz * PI)); });
            for (int i = 0; i < liczba_probek; ++i) {
                if (y[i] >= 0) y[i] = amplituda + ruch_y;
                else y[i] = -amplituda + ruch_y;
            }
            break;
        }
        case 4: {
            double petla = 0;
            bool check = false;
            double param = 0;
            bool rot = false;

            if (amplituda < 0) {
                amplituda *= -1;
                ruch_y *= -1;
                rot = true;
            }

            for (int i = 0; i < liczba_probek; ++i) {
                if (i > 0 && y[i - 1] < -amplituda + ruch_y && check) {
                    petla = i;
                    przesuniecie_faz = 0;
                    param = 0;
                }

                y[i] = 2 * amplituda * param + ruch_y + amplituda - (amplituda * czestotliwosc / PI) * (x[i] - x[petla] - (przesuniecie_faz * PI));

                if (!check) {
                    while (y[i] + param * 2 * amplituda > amplituda) param -= 1;
                    while (y[i] + param * 2 * amplituda < -amplituda) param += 1;
                    check = true;
                }
            }

            if (rot) {
                for (int i = 0; i < liczba_probek; ++i) y[i] *= -1;
            }
            break;
        }
    }

    return y;
}

void rys(vector<double> x, vector<double> y, string name) {
    mp::plot(x, y, "-");
    mp::title("Wykres funkcji " + name);
    mp::xlabel("x");
    mp::ylabel(name);
    mp::grid(true);
    mp::show();
}

void test_gen_sig(int rodzaj_funkcji, double amplituda, double czestotliwosc, double przesuniecie_faz, double ruch_y, int liczba_probek, double dlugosc) {
    vector<double> x = mp::linspace(0, dlugosc * PI, liczba_probek);
    vector<double> y;

    string name;
    y = gen_sig(rodzaj_funkcji, amplituda, czestotliwosc, przesuniecie_faz, ruch_y, liczba_probek, dlugosc);

    switch (rodzaj_funkcji) {
        case 1: name = "sin"; break;
        case 2: name = "cos"; break;
        case 3: name = "square"; break;
        case 4: name = "saw"; break;
    }

    rys(x, y, name);
}

void low_f_filter_real(vector<double>& amplituda, const vector<double>& czestotliwosc, double filter) {
    int i = 0;
    while (czestotliwosc[i] < filter) {
        amplituda[i] = 0;
        ++i;
    }
}

void low_f_filter_im(vector<complex<double>>& spectrum, int filtracja) {
    int k = spectrum.size();
    for (int i = 0; i <= filtracja; ++i) {
        spectrum[i] = 0.0;
    }
    for (int i = k - 1; i >= k - filtracja; --i) {
        spectrum[i] = 0.0;
    }
}

vector<complex<double>> dft_vector(vector<double> zlozenie) {
    int probki = zlozenie.size();
    vector<complex<double>> spektrum(probki);

    for (int m = 0; m < probki; ++m) {
        complex<double> suma = 0.0;
        for (int n = 0; n < probki; ++n) {
            double kat = 2 * PI * m * n / probki;
            suma += zlozenie[n] * polar(1.0, -kat);
        }
        spektrum[m] = suma;
    }

    return spektrum;
}

vector<double> amplitude(vector<complex<double>> spektrum, int podzialka) {
    int liczba_probek = spektrum.size();
    vector<double> amplituda(liczba_probek / podzialka);
    for (int i = 0; i < liczba_probek / podzialka; ++i) {
        amplituda[i] = abs(spektrum[i]) * 2 / liczba_probek;
    }

    return amplituda;
}

vector<double> frequency(vector<complex<double>> spektrum, int podzialka) {
    int liczba_probek = spektrum.size();
    vector<double> czestotliwosc(liczba_probek / podzialka);
    for (int i = 0; i < liczba_probek / podzialka; ++i) {
        czestotliwosc[i] = i;
    }

    return czestotliwosc;
}

void rdft_vector(const vector<complex<double>>& spektrum, vector<double>& zlozenie) {
    int dlugosc = zlozenie.size();
    for (int m = 0; m < dlugosc; ++m) {
        complex<double> suma = 0.0;
        for (int n = 0; n < dlugosc; ++n) {
            double kat = 2 * PI * m * n / dlugosc;
            suma += spektrum[n] * polar(1.0, kat);
        }
        zlozenie[m] = suma.real() / dlugosc;
    }
}

void dft_test1(double amplituda1, double czestotliwosc1, double przesuniecie_faz1, double amplituda2, double czestotliwosc2, double przesuniecie_faz2, double amplituda3, double czestotliwosc3,
    double przesuniecie_faz3, int liczba_probek) {

    vector<double> x = mp::linspace(0, 2 * PI, liczba_probek);

    vector<double> y1 = mp::transform(x, [&](auto x) { return (amplituda1 * sin(czestotliwosc1 * x - przesuniecie_faz1 * PI)); });
    vector<double> y2 = mp::transform(x, [&](auto x) { return (amplituda2 * sin(czestotliwosc2 * x - przesuniecie_faz2 * PI)); });
    vector<double> y3 = mp::transform(x, [&](auto x) { return (amplituda3 * sin(czestotliwosc3 * x - przesuniecie_faz3 * PI)); });
    vector<double> y = mp::transform(x, [&](auto x) { return (amplituda1 * sin(czestotliwosc1 * x - przesuniecie_faz1 * PI) + amplituda2 * sin(czestotliwosc2 * x - przesuniecie_faz2 * PI) + amplituda3 * sin(czestotliwosc3 * x - przesuniecie_faz3 * PI)); });

    rys(x, y, "sin1 + sin2 + sin3");

    vector<complex<double>> spektrum = dft_vector(y);
    vector<double> amplituda = amplitude(spektrum, 2);
    vector<double> czestotliwosc = frequency(spektrum, 2);

    rys(czestotliwosc, amplituda, "spectrum");
}

void dft_test3(int rodzaj_funkcji, double amplituda1, double czestotliwosc1, double przesuniecie_faz1, double amplituda2, double czestotliwosc2, double przesuniecie_faz2, double amplituda3, double czestotliwosc3, double przesuniecie_faz3, int liczba_probek) {
    vector<double> x = mp::linspace(0, 2 * PI, liczba_probek);

    vector<double> y1 = gen_sig(rodzaj_funkcji, amplituda1, czestotliwosc1, przesuniecie_faz1, 0, liczba_probek, 2);
    vector<double> y2 = gen_sig(rodzaj_funkcji, amplituda2, czestotliwosc2, przesuniecie_faz2, 0, liczba_probek, 2);
    vector<double> y3 = gen_sig(rodzaj_funkcji, amplituda3, czestotliwosc3, przesuniecie_faz3, 0, liczba_probek, 2);
    vector<double> y = mp::transform(x, [&](auto x, int i) { return y1[i] + y2[i] + y3[i]; });

    rys(x, y, "signal1 + signal2 + signal3");

    vector<complex<double>> spektrum = dft_vector(y);
    vector<double> amplituda = amplitude(spektrum, 2);
    vector<double> czestotliwosc = frequency(spektrum, 2);

    rys(czestotliwosc, amplituda, "spectrum");
}

void dft_test5(int rodzaj_funkcji, int ilosc_sygnalow, double amplituda1, double czestotliwosc1, double przesuniecie_faz1, double amplituda2, double czestotliwosc2, double przesuniecie_faz2, double amplituda3, double czestotliwosc3, double przesuniecie_faz3, int liczba_probek) {
    vector<double> x = mp::linspace(0, 2 * PI, liczba_probek);

    vector<double> y(liczba_probek);
    if (ilosc_sygnalow == 2) {
        vector<double> y1 = gen_sig(rodzaj_funkcji, amplituda1, czestotliwosc1, przesuniecie_faz1, 0, liczba_probek, 2);
        vector<double> y2 = gen_sig(rodzaj_funkcji, amplituda2, czestotliwosc2, przesuniecie_faz2, 0, liczba_probek, 2);
        y = mp::transform(x, [&](auto x, int i) { return y1[i] + y2[i]; });
    } else if (ilosc_sygnalow == 3) {
        vector<double> y1 = gen_sig(rodzaj_funkcji, amplituda1, czestotliwosc1, przesuniecie_faz1, 0, liczba_probek, 2);
        vector<double> y2 = gen_sig(rodzaj_funkcji, amplituda2, czestotliwosc2, przesuniecie_faz2, 0, liczba_probek, 2);
        vector<double> y3 = gen_sig(rodzaj_funkcji, amplituda3, czestotliwosc3, przesuniecie_faz3, 0, liczba_probek, 2);
        y = mp::transform(x, [&](auto x, int i) { return y1[i] + y2[i] + y3[i]; });
    }

    rys(x, y, "combined signal");

    vector<complex<double>> spektrum = dft_vector(y);
    vector<double> amplituda = amplitude(spektrum, 2);
    vector<double> czestotliwosc = frequency(spektrum, 2);

    rys(czestotliwosc, amplituda, "spectrum");
}

PYBIND11_MODULE(signal_analysis, m) {
    m.def("gen_sig", &gen_sig, "Generate signal");
    m.def("rys", &rys, "Draw plot");
    m.def("test_gen_sig", &test_gen_sig, "Test signal generation");
    m.def("low_f_filter_real", &low_f_filter_real, "Low-pass filter for real signals");
    m.def("low_f_filter_im", &low_f_filter_im, "Low-pass filter for complex signals");
    m.def("dft_vector", &dft_vector, "Discrete Fourier Transform");
    m.def("amplitude", &amplitude, "Get amplitude spectrum");
    m.def("frequency", &frequency, "Get frequency spectrum");
    m.def("dft_test1", &dft_test1, "Test DFT with 1 combined signal");
    m.def("dft_test3", &dft_test3, "Test DFT with 3 combined signals");
    m.def("dft_test5", &dft_test5, "Test DFT with multiple signals");
}
