#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <matplot/matplot.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <complex>
#include <AudioFile.h>

//domyślny namespace
using namespace std;

//aliasy namespace'ów do użytku
namespace py = pybind11;
namespace mp = matplot;

// zdefiniowanie liczby PI
const double PI = acos(-1.0);

//DEKLARACJE FUNKCJI

vector<double> gen_sig(int, double, double, double, double, int, double);
void rys(vector<double>, vector<double>, string);
void test_gen_sig(int, double, double, double, double, int, double);
void low_f_filter_real(vector<double>&, const vector<double>&, double);
void low_f_filter_im(vector<complex<double>>&, int);
vector<complex<double>> dft_vector(vector<double>);
vector<double> amplitude(vector<complex<double>>, int);
vector<double> frequency(vector<complex<double>>, int);
vector<double> rdft_vector(vector<complex<double>>);
void dft_test1(double, double, double, int);
void dft_test3(double, double, double, int, double);
void dft_test5(double, double, double, int, int);




vector<double> gen_sig(int rodzaj_funkcji, double amplituda, double czestotliwosc, double przesuniecie_faz, double ruch_y, int liczba_probek, double dlugosc)
{
    vector<double> x = mp::linspace(0, dlugosc * PI, liczba_probek);
    vector<double> y;

    switch (rodzaj_funkcji)
    {
        case 1:
        {
            y = mp::transform(x, [&](auto x) { return (ruch_y + amplituda * sin(czestotliwosc * x - przesuniecie_faz * PI)); });
            break;
        }
        case 2:
        {
            y = mp::transform(x, [&](auto x) { return (ruch_y + amplituda * cos(czestotliwosc * x - przesuniecie_faz * PI)); });
            break;
        }
        case 3:
        {
            y = mp::transform(x, [&](auto x) { return (sin(czestotliwosc * x - przesuniecie_faz * PI)); });

            for (int i = 0; i < liczba_probek; ++i)
            {
                if (y[i] >= 0)
                    y[i] = amplituda + ruch_y;
                else
                    y[i] = -amplituda + ruch_y;
            }
            break;
        }
        case 4:
        {
            double petla = 0;
            bool check = false;
            double param = 0;
            bool rot = false;

            if (amplituda < 0)
            {
                amplituda *= -1;
                ruch_y *= -1;

                rot = true;
            }


            for (int i = 0; i < liczba_probek; ++i)
            {

                if (y[i - 1] < -amplituda && check == true)
                {
                    petla = i;
                    przesuniecie_faz = 0;
                    param = 0;
                }

                y[i] = 2 * amplituda * param + ruch_y + amplituda - (amplituda * czestotliwosc / PI) * (x[i] - x[petla] - (przesuniecie_faz * PI));


                if (check == false)
                {
                    while (y[i] + param * 2 * amplituda > amplituda)
                    {
                        param -= 1;
                    }

                    while (y[i] + param * 2 * amplituda < -amplituda)
                    {
                        param += 1;
                    }

                    check = true;
                }

            }

            if (rot == true)
            {
                for (int i = 0; i < liczba_probek; ++i)
                {
                    y[i] *= -1;
                }
            }
            break;
        }

    }

    return y;
}


void rys(vector<double> x, vector<double> y, string name)
{
	mp::plot(x, y, "-");
	mp::title("Wykres funkcji " + name);
	mp::xlabel("x");
	mp::ylabel(name);
	mp::grid(true);
	mp::show();

}

void test_gen_sig(int rodzaj_funkcji, double amplituda, double czestotliwosc, double przesuniecie_faz, double ruch_y, int liczba_probek, double dlugosc)
{
    vector<double> x = mp::linspace(0, dlugosc * PI, liczba_probek);
    vector<double> y;

    y = gen_sig(rodzaj_funkcji, amplituda, czestotliwosc, przesuniecie_faz, ruch_y, liczba_probek, dlugosc);

    rys(x, y, to_string(rodzaj_funkcji));
}

void low_f_filter_real(vector<double>& amplituda, const vector<double>& czestotliwosc, double filter)
{
    int i = 0;

    while (czestotliwosc[i] < filter)
    {
        amplituda[i] = 0;

        ++i;
    }
}

void low_f_filter_im(vector<complex<double>>& spectrum, int filtracja) 
{
    int k = spectrum.size();

    for (int i = 0; i <= filtracja; ++i)
    {
        spectrum[i] = 0.0;
    }

    for (int i = k; i >= k - filtracja; --i)
    {
        spectrum[i] = 0.0;
    }
}




vector<complex<double>> dft_vector(vector<double> zlozenie)
{
    int probki = zlozenie.size();

    vector<complex<double>> spektrum(probki);

    for (int m = 0; m < probki; ++m)
    {
        complex<double> suma = 0.0;

        for (int n = 0; n < probki; ++n)
        {
            double kat = 2 * PI * m * n / probki;
            suma += zlozenie[n] * polar(1.0, -kat);
        }
        spektrum[m] = suma;
    }

    return spektrum;
}

vector<double> amplitude(vector<complex<double>> spektrum, int podzialka)
{

    int liczba_probek = spektrum.size();
    vector<double> amplituda(liczba_probek / podzialka);
    for (int i = 0; i < liczba_probek / podzialka; ++i) 
    {
        amplituda[i] = abs(spektrum[i]) * 2 / liczba_probek;
    }

    return amplituda;
}

vector<double> frequency(vector<complex<double>> spektrum, int podzialka)
{

    int liczba_probek = spektrum.size();
    vector<double> czestotliwosc(liczba_probek / podzialka);
    for (int i = 0; i < liczba_probek / podzialka; ++i) 
    {
        czestotliwosc[i] = i;
    }

    return czestotliwosc;
}

vector<double> rdft_vector(vector<complex<double>> spektrum)
{
    int dlugosc = spektrum.size();

    vector<double> zlozenie2(dlugosc);

    for (int m = 0; m < dlugosc; ++m)
    {
        complex<double> suma = 0.0;
        for (int n = 0; n < dlugosc; ++n)
        {
            double kat = 2 * PI * m * n / dlugosc;
            suma += spektrum[n] * polar(1.0, kat);
        }
        zlozenie2[m] = suma.real() / dlugosc;
    }

    return zlozenie2;
}

void dft_test1(double czest1, double czest2, double czest3, int liczba_probek)
{

    vector<double> x = mp::linspace(0, 2 * PI, liczba_probek);

    vector<double> zlozenie = mp::transform(x, [&](auto x) { return (sin(czest1 * x) + sin(czest2 * x) + sin(czest3 * x)); });


    rys(x, zlozenie, "zlozenie funkcji");
}

void dft_test3(double czest1, double czest2, double czest3, int a, double filtr)
{
    int liczba_probek = 10 * a;
    int zakres = a;

    vector<double> x = mp::linspace(0, 2 * PI, liczba_probek);

    vector<double> zlozenie = mp::transform(x, [&](auto t) { return sin(czest1 * t) + sin(czest2 * t) + sin(czest3 * t); });

    vector<complex<double>> spektrum(liczba_probek);

    spektrum = dft_vector(zlozenie);

    vector<double> amplituda(liczba_probek / 10);
    
    amplituda = amplitude(spektrum, 10);

    vector<double> czestotliwosc(liczba_probek / 10);
    
    czestotliwosc = frequency(spektrum, 10);
  
    low_f_filter_real(amplituda, czestotliwosc, filtr);

    rys(czestotliwosc, amplituda, "dft");

}

void dft_test5(double czest1, double czest2, double czest3, int a, int ft)
{
    int liczba_probek = 10 * a;
    int zakres = a;

    vector<double> x = mp::linspace(0, 2 * PI, liczba_probek);

    vector<double> zlozenie = mp::transform(x, [&](auto t) { return sin(czest1 * t) + sin(czest2 * t) + sin(czest3 * t); });

    vector<complex<double>> spektrum(liczba_probek);

    spektrum = dft_vector(zlozenie);

    low_f_filter_im(spektrum, ft);

    vector<double> zlozenie2(liczba_probek);

    zlozenie2 = rdft_vector(spektrum);

    rys(x, zlozenie2, "dft");

}
/////////////////////////////////////////////////////TESTY
PYBIND11_MODULE(projekt, m) {

	m.def("signal", &gen_sig);
	m.def("show", &rys);
    m.def("test_sig", &test_gen_sig);
    m.def("filter_real", &low_f_filter_real);
    m.def("dft", &dft_vector);
    m.def("amplitude", &amplitude);
    m.def("frequency", &frequency);
    m.def("rdft", &rdft_vector);
    m.def("splot_test", &dft_test1);
    m.def("dft_test", &dft_test3);
    m.def("rdft_test", &dft_test5);
}