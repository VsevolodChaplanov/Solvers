#ifndef __VECTOROPERATIONS__
#define __VECTOROPERATIONS__

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

// Процедуры основных векторных и матричных операций

// Скалярное произведение векторов
double DotProduct(const std::vector<double> &a, const std::vector<double> &b)
{
    double Result = 0;

    for (size_t i = 0; i < a.size(); i++)
    {
        Result += a[i] * b[i];
    } // i

    return Result;
}

// Умножение вектора на число
std::vector<double> Mult_N(const std::vector<double> &a, const double &b)
{
    std::vector<double> Result(a.size(), 0.);
    for (size_t i = 0; i < a.size(); i++)
    {
        Result[i] = a[i] * b;
    } // i
    return Result;
}

// Сложение векторов
std::vector<double> VSum(const std::vector <double> &a, const std::vector <double> &b)
{
    std::vector<double> Result(a.size(), 0.);
    for (size_t i = 0; i < a.size(); i++)
    {
        Result[i] = a[i] + b[i];
    }
    return Result;
}

// Вычитание векторов
std::vector<double> VDiff(const std::vector<double> &a, const std::vector<double> &b)
{
    std::vector<double> Result(a.size(), 0.);
    for (size_t i = 0; i < a.size(); i++)
    {
        Result[i] = a[i] - b[i];
    } // i
    return Result;
}

// Запись вектора a в файл с названием Filename.csv
void WriteInFile(const std::vector<double> &a, const std::string &Filename)
{
	std::ofstream file;
    file.open(Filename + ".csv");
    if (file.is_open())
    {
		for (auto elem : a)
		{
			file << elem << std::endl;
		} // elem
    } //endif
    file.close();
}

// Норма по максимальному отклонению || r ||_{max}
double NormMax(const std::vector<double> &a, const std::vector<double> &b)
{
	double MaxElem = 0;
	for (size_t i = 0; i < a.size(); i++)
	{
		double diff = fabs(a[i] - b[i]);
		diff > MaxElem ? MaxElem = diff : MaxElem = MaxElem;
	}
	return MaxElem;
}

// Вторая норма || r ||_{2}
double SecondNorm(const std::vector<double> &a)
{
	double Result = DotProduct(a,a);
	return sqrt(Result);
}

bool CheckMatSym(CMatrix &Matrix)
{
    bool result = true;
    for (size_t i = 0; i < Matrix.size(); i++)
    {
        for (auto elem : Matrix[i])
        {
            if (elem.second != Matrix.GetValue(elem.first, i))
            {
                result = false;
            }
        }
    }
    return result;
}

#endif