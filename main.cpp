#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <matplot/matplot.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <string>

// zdefiniowanie liczby PI
#define PI 3.14159265

//aliasy namespace'ów do użytku
namespace py = pybind11;
namespace mp = matplot;

//domyślny namespace
using namespace std;

//testowa funkcja do usunięcia po prostu sprawdzałem czy się kompiluje jak należy
float square(float x) { return x * x; }
//////////////////////////////////////

////mój test DO POTĘŻNEJ PRZEBUDOWY
//int test_generowania_sin()
//{
//	int n = 0;
//	double amp = 1 /* NIE MOŻE BYĆ RÓWNE 0 */, wsp_cz = 1 /* NIE MOŻE BYĆ RÓWNE 0 */, faza = 0;
//	cin >> n;
//	double* list = new double[n] {0};
//	for (int i = 0; i < n; i++)
//	{
//		//Generowanie dyskretnej funkcji sinusa (ciągła funkcja jest dla komputera nieosiągalna)
//		list[i] = amp * sin(wsp_cz * (i * 2 * PI / (n - 1) + faza));
//		//Korekcja dokładności dla 0
//		if (list[i] < 1e-6 && list[i] > -1e-6) list[i] = 0;
//	}
//	//------DEBUG------
//	/*for (int i = 0; i < n; i++)
//	{
//		cout << list[i] << ' ';
//	}
//	cout << endl << "#probek/okres: " << n << " , skok probki: " << 360. / (n - 1) << " stopni" << endl;
//	for (int i = 0; i < n; i++)
//	{
//		list[i] = NULL;
//	}*/
//	//-----------------
//	delete[] list;
//	list = nullptr;
//	return 0;
//}



//int wykres() {
//	using namespace matplot;
//	std::vector<double> x = linspace(0, 2 * pi);
//	std::vector<double> y = transform(x, [](auto x) { return sin(x); });
//
//	plot(x, y, "-o");
//	hold(on);
//	plot(x, transform(y, [](auto y) { return -y; }), "--xr");
//	plot(x, transform(x, [](auto x) { return x / pi - 1.; }), "-:gs");
//	plot({ 1.0, 0.7, 0.4, 0.0, -0.4, -0.7, -1 }, "k");
//
//	show();
//	return 0;
//}

//void rys_sin(int f, double dlugosc, int liczba_punktow)
//{
//	vector<double> x = mp::linspace(0, dlugosc * PI, liczba_punktow);
//	vector<double> y = mp::transform(x, [&](auto x) { return sin(f*x); });
//
	//mp::plot(x, y, "-");
	//mp::title("Wykres funkcji sin(x) o częstotliwości " + to_string(f));
	//mp::xlabel("x");
	//mp::ylabel("sin(x)");
	//mp::grid(true);
	//mp::show();
//
//}

vector<double> gen_sig(int rodzaj_funkcji, double amplituda, double czestotliwosc, double przesuniecie_faz, double ruch_y, int liczba_probek, double dlugosc)
{
    vector<double> x = mp::linspace(0, dlugosc * PI, liczba_probek);;
    vector<double> y;



    switch (rodzaj_funkcji)
    {
        case 1:
        {
            y = mp::transform(x, [&](auto x) { return (ruch_y + amplituda * sin(czestotliwosc * x + przesuniecie_faz * PI)); });
            break;
        }
        case 2:
        {
            y = mp::transform(x, [&](auto x) { return (ruch_y + amplituda * cos(czestotliwosc * x + przesuniecie_faz * PI)); });
            break;
        }
        case 3:
        {
            y = mp::transform(x, [&](auto x) { return (sin(czestotliwosc * x + przesuniecie_faz * PI)); });

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

            if (amplituda > 0)
            {
                for (int i = 0; i < liczba_probek; ++i)
                {
                    y[i] = amplituda - (amplituda * czestotliwosc / (2 * PI)) * ((i - przesuniecie_faz - petla));

                    if (y[i] <= -amplituda) petla = i;
                }
            }
            else
            {
                for (int i = 0; i < liczba_probek; ++i)
                {
                    y[i] = amplituda + (amplituda * czestotliwosc / (2 * PI)) * ((i - przesuniecie_faz - petla));

                    if (y[i] <= amplituda) petla = i;
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
    vector<double> x = mp::linspace(0, dlugosc * PI, liczba_probek);;
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
        case 4: // Piłokształtny nie ma przesunięcia fazowego
        {
            double petla = 0;
            double temp = 1;
            bool check = false;
            double param = 0;
            bool rot = false;

            if (amplituda < 0)
            {
                amplituda *= -1;

                rot = true;
            }


            for (int i = 0; i < liczba_probek; ++i)
            {
                
                    /*y[i] = amplituda - (amplituda * czestotliwosc / (2 * PI)) * ((i - przesuniecie_faz - petla));
                    if (y[i] <= -amplituda) petla = i;*/

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
//////////////////////////////////////////////////////////////////////////przekazywanie bez kopiowania

PYBIND11_MODULE(projekt, m) {
	//również do usunięcia póki co zostawiam jako template
	m.def("square", &square);
	///////////////////////////////////
	/*m.def("wykres", &wykres);*/

	/*m.def("rys_sin", &rys_sin);*/

	m.def("signal", &gen_sig);
	m.def("show", &rys);
    m.def("test_sig", &test_gen_sig);
}