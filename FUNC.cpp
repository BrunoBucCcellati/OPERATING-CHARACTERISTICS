#include "pch.h"
#include <random>
#define M_PI 3.14159265358979323846
template <typename T>
short Sign(T Value)
{
    if (Value == 0.)  return 0;
    if (Value > 0.)  return 1;
    else return -1;
}
extern "C" __declspec(dllexport) double HillFunc(double x)
{
    const unsigned short HH_COUNT = 14;
    double a = 0.0;
    double b = 1.0;
    //
    srand((unsigned short)time(NULL));
    std::default_random_engine generator(rand());
    std::uniform_real_distribution<double> distribution(a, b);
    //
    double res = -1.1 + distribution(generator) * 2.0;
    //
    unsigned short i = 1;
    //
    while (i <= HH_COUNT)
    {
        res += (-1.1 + distribution(generator) * 2.0) * sin(M_PI * double(i) * x * 2.0) + (-1.1 + distribution(generator) * 2.0) * cos(M_PI * double(i) * x * 2.0);
        i++;
    }
    return res;
}
extern "C" __declspec(dllexport) double ShekelFunc(double x)
{
    const unsigned short SH_COUNT = 10;
    double a = 0.0;
    double b = 1.0;
    //
    srand((unsigned short)time(NULL));
    std::default_random_engine generator(rand());
    std::uniform_real_distribution<double> distribution(a, b);
    //
    double res = 0.0;
    //
    unsigned short i = 1;
    //
    while (i <= SH_COUNT)
    {
        res -= 1.0 / ((5.0 + 20.0 * distribution(generator)) * pow((x - 10.0 * distribution(generator)), 2.0) + 1.0 + 0.2 * distribution(generator) * 2.0);
        i++;
    }
    //
    return res;
}
extern "C" __declspec(dllexport) double GrishaginFunc(double x1, double x2)
{
    const unsigned short GR_COUNT = 7;
    double a = -1.0;
    double b = 1.0;
    //
    srand((unsigned short)time(NULL));
    std::default_random_engine generator(rand());
    std::uniform_real_distribution<double> distribution(a, b);
    //
    double part1 = 0.0, part2 = 0.0;
    //
    unsigned short i = 1, j = 1;
    //
    while (i <= GR_COUNT)
    {
        while (j <= GR_COUNT)
        {
            part1 += distribution(generator) * sin(M_PI * double(i) * x1) * sin(M_PI * double(j) * x2) + distribution(generator) * cos(M_PI * double(i) * x1) * cos(M_PI * double(j) * x2);
            part2 += distribution(generator) * sin(M_PI * double(i) * x1) * sin(M_PI * double(j) * x2) - distribution(generator) * cos(M_PI * double(i) * x1) * cos(M_PI * double(j) * x2);
            j++;
        }
        i++;
    }
    return -sqrt(pow(part1, 2.0) + pow(part2, 2.0));
}
extern "C" __declspec(dllexport) double Characteristic(double _m, double x1, double x2, double y1, double y2, unsigned short _N)
{
    return ((pow(abs(x2 - x1), (1 / double(_N)))) + pow((y2 - y1), 2) / (pow(_m, 2) * (pow(abs(x2 - x1), (1 / double(_N))))) - 2 * (y2 + y1) / _m);
}
extern "C" __declspec(dllexport) double Shag(double _m, double x1, double x2, double y1, double y2, unsigned short _N)
{
    return ((x1 + x2) / 2) - Sign(y2 - y1) * pow((abs(y2 - y1) / (2 * _m)), _N);
}