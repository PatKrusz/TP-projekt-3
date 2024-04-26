#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <iostream>
#include <cmath>
#include <vector>

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

//mój test DO POTĘŻNEJ PRZEBUDOWY
int test_generowania_sin()
{
	int n = 0;
	double amp = 1 /* NIE MOŻE BYĆ RÓWNE 0 */, wsp_cz = 1 /* NIE MOŻE BYĆ RÓWNE 0 */, faza = 0;
	cin >> n;
	double* list = new double[n] {0};
	for (int i = 0; i < n; i++)
	{
		//Generowanie dyskretnej funkcji sinusa (ciągła funkcja jest dla komputera nieosiągalna)
		list[i] = amp * sin(wsp_cz * (i * 2 * PI / (n - 1) + faza));
		//Korekcja dokładności dla 0
		if (list[i] < 1e-6 && list[i] > -1e-6) list[i] = 0;
	}
	//------DEBUG------
	/*for (int i = 0; i < n; i++)
	{
		cout << list[i] << ' ';
	}
	cout << endl << "#probek/okres: " << n << " , skok probki: " << 360. / (n - 1) << " stopni" << endl;
	for (int i = 0; i < n; i++)
	{
		list[i] = NULL;
	}*/
	//-----------------
	delete[] list;
	list = nullptr;
	return 0;
}



int wykres() {
	using namespace matplot;
	std::vector<double> x = linspace(0, 2 * pi);
	std::vector<double> y = transform(x, [](auto x) { return sin(x); });

	plot(x, y, "-o");
	hold(on);
	plot(x, transform(y, [](auto y) { return -y; }), "--xr");
	plot(x, transform(x, [](auto x) { return x / pi - 1.; }), "-:gs");
	plot({ 1.0, 0.7, 0.4, 0.0, -0.4, -0.7, -1 }, "k");

	show();
	return 0;
}



PYBIND11_MODULE(projekt, m) {
	//również do usunięcia póki co zostawiam jako template
	m.def("square", &square);
	///////////////////////////////////
	m.def("wykres", &wykres);
}