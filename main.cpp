#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <matplot/matplot.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <complex>

// zdefiniowanie liczby PI
#define PI 3.14159265

//aliasy namespace'ów do użytku
namespace py = pybind11;
namespace mp = matplot;

//domyślny namespace
using namespace std;

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
    vector<double> y(liczba_probek);



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


    rys(x, y, to_string(rodzaj_funkcji));
}

void low_f_filter(vector<double>& amplituda, const vector<double>& czestotliwosc, double filter)
{
    int i = 0;

    while (czestotliwosc[i] < filter)
    {
        amplituda[i] = 0;

        ++i;
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
            suma += complex<double>(zlozenie[n] * cos(kat), -zlozenie[n] * sin(kat));
        }
        spektrum[m] = suma;
    }

    return spektrum;
}

vector<double> amplitude(vector<complex<double>> spektrum)
{

    int liczba_probek = spektrum.size();
    vector<double> amplituda(liczba_probek / 2);
    for (int i = 0; i < liczba_probek / 2; ++i) 
    {
        amplituda[i] = abs(spektrum[i]) * 2 / liczba_probek;
    }

    return amplituda;
}

vector<double> frequency(vector<complex<double>> spektrum)
{

    int liczba_probek = spektrum.size();
    vector<double> czestotliwosc(liczba_probek / 2);
    for (int i = 0; i < liczba_probek / 2; ++i) 
    {
        czestotliwosc[i] = i;
    }

    return czestotliwosc;
}

void dft_test1(double czest1, double czest2, double czest3)
{
    int liczba_probek = 1000;

    vector<double> x = mp::linspace(0, 4 * PI, liczba_probek);

    vector<double> f1 = mp::transform(x, [&](auto x) { return sin(czest1 * x); });
    vector<double> f2 = mp::transform(x, [&](auto x) { return sin(czest2 * x); });
    vector<double> f3 = mp::transform(x, [&](auto x) { return sin(czest3 * x); });

    vector<double> zlozenie(liczba_probek);

    for (int i = 0; i < liczba_probek; ++i)
    {
        zlozenie[i] = f1[i] + f2[i] + f3[i];
    }

    rys(x, zlozenie, "zlozenie funkcji");
}

void dft_test3(double czest1, double czest2, double czest3, int a, double filtr)
{
    int liczba_probek = 2 * a;
    int czest_probk = 2 * a;
    vector<double> x = mp::linspace(0, 2 * PI, liczba_probek);

    vector<double> f1 = mp::transform(x, [&](auto x) { return sin(czest1 * x); });
    vector<double> f2 = mp::transform(x, [&](auto x) { return sin(czest2 * x); });
    vector<double> f3 = mp::transform(x, [&](auto x) { return sin(czest3 * x); });

    vector<double> zlozenie(liczba_probek);

    for (int i = 0; i < liczba_probek; ++i)
    {
        zlozenie[i] = f1[i] + f2[i] + f3[i];
    }

    vector<complex<double>> spektrum(liczba_probek);

    for (int m = 0; m < liczba_probek; ++m)
    {
        complex<double> suma = 0.0;

        for (int n = 0; n < liczba_probek; ++n)
        {
            double kat = 2 * PI * m * n / liczba_probek;
            suma += complex<double>(zlozenie[n] * cos(kat), -zlozenie[n] * sin(kat));
        }
        spektrum[m] = suma;
    }

    vector<double> amplituda(liczba_probek / 2);
    for (int i = 0; i < liczba_probek / 2; ++i) 
    {
        amplituda[i] = abs(spektrum[i]) * 2 / liczba_probek;
    }

    vector<double> czestotliwosc(liczba_probek / 2);
    for (int i = 0; i < liczba_probek / 2; ++i) 
    {
        czestotliwosc[i] = i;
    }

    low_f_filter(amplituda, czestotliwosc, filtr);

    rys(czestotliwosc, amplituda, "dft");

}

void dft_test4(double czest1, double czest2, double czest3, int a)
{
    int liczba_probek = 2 * a;
    int czest_probk = 2 * a;
    vector<double> x = mp::linspace(0, 2 * PI, liczba_probek);

    vector<double> f1 = mp::transform(x, [&](auto x) { return sin(czest1 * x); });
    vector<double> f2 = mp::transform(x, [&](auto x) { return sin(czest2 * x); });
    vector<double> f3 = mp::transform(x, [&](auto x) { return sin(czest3 * x); });

    vector<double> zlozenie(liczba_probek);

    for (int i = 0; i < liczba_probek; ++i)
    {
        zlozenie[i] = f1[i] + f2[i] + f3[i];
    }

    vector<complex<double>> spektrum(liczba_probek);

    for (int m = 0; m < liczba_probek; ++m)
    {
        complex<double> suma = 0.0;

        for (int n = 0; n < liczba_probek; ++n)
        {
            double kat = 2 * PI * m * n / liczba_probek;
            suma += complex<double>(zlozenie[n] * cos(kat), -zlozenie[n] * sin(kat));
        }
        spektrum[m] = suma;
    }

    vector<double> amplituda(liczba_probek / 2);
    for (int i = 0; i < liczba_probek / 2; ++i) 
    {
        amplituda[i] = abs(spektrum[i]) * 2 / liczba_probek;
    }

    vector<double> czestotliwosc(liczba_probek / 2);
    for (int i = 0; i < liczba_probek / 2; ++i) 
    {
        czestotliwosc[i] = i;
    }

    vector<double> zlozenie2(liczba_probek);

    for (int n = 0; n < liczba_probek; ++n)
    {
        complex<double> suma = 0.0;
        for (int k = 0; k < liczba_probek; ++k) 
        {
            double kat = 2 * PI * n * k / liczba_probek;
            suma += spektrum[k] * polar(1.0, kat);
        }
        zlozenie2[n] = suma.real() / liczba_probek;
    }

    rys(x, zlozenie2, "dft");

}

PYBIND11_MODULE(projekt, m) {

	m.def("signal", &gen_sig);
	m.def("show", &rys);
    m.def("test_sig", &test_gen_sig);
    m.def("splot_test", &dft_test1);
    m.def("dft_test", &dft_test3);
    m.def("filtr", &low_f_filter);
    m.def("rdft_test", &dft_test4);
}