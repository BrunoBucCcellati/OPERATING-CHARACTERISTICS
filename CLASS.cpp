#include "pch.h"
#include "framework.h"
#include "Windows.h"
#include<algorithm>
#include <random>
using namespace System::Collections::Generic;
namespace OC
{
    public ref class Otrezok
    {
    private:
        KeyValuePair<double, double> start;
        KeyValuePair<double, double> end;
        double M;
        double R;
        Otrezok^ Next;
        Otrezok^ Previous;
        double m;
    public:
        Otrezok(KeyValuePair<double, double> _start, KeyValuePair<double, double> _end, unsigned short _N)
        {
            start = _start;
            end = _end;
            M = abs((end.Value - start.Value) / pow(abs(end.Key - start.Key), (1.0 / double(_N))));
            R = 0.0;
            Next = this;
            Previous = this;
        }
        void ChangeCharacteristic(double _m, unsigned short _N)
        {
            HINSTANCE load_function = LoadLibrary(L"FUNC.dll");
            typedef double (*characteristic) (double, double, double, double, double, unsigned short);
            characteristic Characteristic = (characteristic)GetProcAddress(load_function, "Characteristic");
            R = Characteristic(_m, start.Key, end.Key, start.Value, end.Value, _N);
            FreeLibrary(load_function);
        }
        double GetCharacteristic()
        {
            return R;
        }
        void SetEnd(KeyValuePair<double, double> _end, unsigned short _N)
        {
            end = _end;
            M = abs((end.Value - start.Value) / pow(abs(end.Key - start.Key), (1.0 / double(_N))));
        }
        KeyValuePair<double, double> GetEnd()
        {
            return end;
        }
        KeyValuePair<double, double> GetStart()
        {
            return start;
        }
        Otrezok^ GetNext()
        {
            return Next;
        }
        void SetNext(Otrezok^ _Next)
        {
            Next = _Next;
        }
        Otrezok^ GetPrevious()
        {
            return Previous;
        }
        void SetPrevious(Otrezok^ _Previous)
        {
            Previous = _Previous;
        }
        double GetM()
        {
            return M;
        }
        void SetM(double _M)
        {
            M = _M;
        }
    };
    //
    public ref class Otrezki
    {
    private:
        Otrezok^ Head;
        double Mmax;
        double m;
        double r;
        KeyValuePair<double, double> x_Rmax;
        KeyValuePair<double, double> y_Rmax;
        double dmax;
    public:
        Otrezki(Otrezok^ _Head, double _r, unsigned short _N)
        {
            Head = _Head;
            r = _r;
            if (Head->GetM() != 0.0)
            {
                Mmax = Head->GetM();
                m = r * Mmax;
            }
            else
            {
                Mmax = 0.0;
                m = 1.0;
            }
            Head->ChangeCharacteristic(m, _N);
            x_Rmax = KeyValuePair<double, double>((Head->GetStart()).Key, (Head->GetEnd()).Key);
            y_Rmax = KeyValuePair<double, double>((Head->GetStart()).Value, (Head->GetEnd()).Value);
        }
        void Add(KeyValuePair<double, double> tmp, unsigned short _N, double _eta_0, unsigned short _global_local_iterations, unsigned short _K, double _a, double _b, double _epsilon_granichnoe)
        {
            Otrezok^ curr = Head;
            while ((curr->GetEnd()).Key < tmp.Key)
            {
                curr = curr->GetNext();
            }
            Otrezok^ curr1 = gcnew Otrezok(tmp, curr->GetEnd(), _N);
            curr->SetEnd(tmp, _N);
            curr1->SetPrevious(curr);
            curr->GetNext()->SetPrevious(curr1);
            curr1->SetNext(curr->GetNext());
            curr->SetNext(curr1);
            //
            if (_K >= _global_local_iterations && _K % 2 == 0)
            {
                if (curr->GetM() >= curr1->GetM())
                {
                    double eta_shtrih = max(curr->GetM(), max(curr->GetPrevious()->GetM(), curr->GetNext()->GetM()));
                    double eta_2shtrih = Mmax * pow(((curr->GetEnd()).Key - (curr->GetStart()).Key), double(1.0 / _N)) / dmax;
                    m = max((1.0 / r) * eta_shtrih + ((r - 1.0) / r) * eta_2shtrih, _eta_0) * r;
                    x_Rmax = KeyValuePair<double, double>((curr->GetStart()).Key, (curr->GetEnd()).Key);
                    y_Rmax = KeyValuePair<double, double>((curr->GetStart()).Value, (curr->GetEnd()).Value);
                }
                else
                {
                    double eta_shtrih = max(curr1->GetM(), max(curr1->GetPrevious()->GetM(), curr1->GetNext()->GetM()));
                    double eta_2shtrih = Mmax * pow(((curr1->GetEnd()).Key - (curr1->GetStart()).Key), double(1.0 / _N)) / dmax;
                    m = max((1.0 / r) * eta_shtrih + ((r - 1.0) / r) * eta_2shtrih, _eta_0) * r;
                    x_Rmax = KeyValuePair<double, double>((curr1->GetStart()).Key, (curr1->GetEnd()).Key);
                    y_Rmax = KeyValuePair<double, double>((curr1->GetStart()).Value, (curr1->GetEnd()).Value);
                }
            }
            else
            {
                curr = Head;
                Otrezok^ Otrezok_Rmax = curr;
                Mmax = curr->GetM();
                curr = curr->GetNext();
                while (curr != Head)
                {
                    if (curr->GetM() > Mmax)
                    {
                        Mmax = curr->GetM();
                        Otrezok_Rmax = curr;
                    }
                    curr = curr->GetNext();
                }
                Otrezok_Rmax->SetM(0.0);
                if (Mmax != 0.0)
                {
                    m = r * Mmax;
                }
                else
                {
                    m = 1.0;
                }
                //
                curr = Head;
                while (curr->GetNext() != Head)
                {
                    curr = curr->GetNext();
                    curr->ChangeCharacteristic(m, _N);
                }
                curr = curr->GetNext();
                curr->ChangeCharacteristic(m, _N);
                //
                double Rmax = -DBL_MAX;
                HINSTANCE load_function = LoadLibrary(L"FUNC.dll");
                typedef double (*shag) (double, double, double, double, double, unsigned short);
                shag Shag = (shag)GetProcAddress(load_function, "Shag");
                if (_a + _epsilon_granichnoe <= Shag(m, (curr->GetStart()).Key, (curr->GetEnd()).Key, (curr->GetStart()).Value, (curr->GetEnd()).Value, _N) <= _b - _epsilon_granichnoe)
                {
                    Rmax = curr->GetCharacteristic();
                    x_Rmax = KeyValuePair<double, double>((curr->GetStart()).Key, (curr->GetEnd()).Key);
                    y_Rmax = KeyValuePair<double, double>((curr->GetStart()).Value, (curr->GetEnd()).Value);
                }
                dmax = pow((curr->GetEnd()).Key - (curr->GetStart()).Key, double(1.0 / _N));
                curr = curr->GetNext();
                while (curr != Head)
                {
                    if (curr->GetCharacteristic() > Rmax && _a + _epsilon_granichnoe <= Shag(m, (curr->GetStart()).Key, (curr->GetEnd()).Key, (curr->GetStart()).Value, (curr->GetEnd()).Value, _N) <= _b - _epsilon_granichnoe)
                    {
                        Rmax = curr->GetCharacteristic();
                        x_Rmax = KeyValuePair<double, double>((curr->GetStart()).Key, (curr->GetEnd()).Key);
                        y_Rmax = KeyValuePair<double, double>((curr->GetStart()).Value, (curr->GetEnd()).Value);
                    }
                    if (pow((curr->GetEnd()).Key - (curr->GetStart()).Key, double(1.0 / _N)) > dmax)
                    {
                        dmax = pow((curr->GetEnd()).Key - (curr->GetStart()).Key, double(1.0 / _N));
                    }
                    curr = curr->GetNext();
                }
                FreeLibrary(load_function);
            }
        }
        double Getm()
        {
            return m;
        }
        KeyValuePair<double, double> GetX_Rmax()
        {
            return x_Rmax;
        }
        KeyValuePair<double, double> GetY_Rmax()
        {
            return y_Rmax;
        }
    };
    //
    enum List { Top, Dawn, Left, Right };
    public ref class CenterSquare
    {
    private:
        KeyValuePair<double, double> x1x2;
        CenterSquare^ DawnLeft;
        CenterSquare^ TopLeft;
        CenterSquare^ TopRight;
        CenterSquare^ DawnRight;
        List Type;
        CenterSquare^ Start;
        CenterSquare^ End;
        CenterSquare^ Next;
        CenterSquare^ Previous;
    public:
        CenterSquare(double _x1, double _x2, List _Type, CenterSquare^ _Next, CenterSquare^ _Previous, unsigned short i, double h1, double h2)
        {
            x1x2 = KeyValuePair<double, double>(_x1, _x2);
            Type = _Type;
            if (_Next == nullptr)
            {
                Next = this;
            }
            else
            {
                Next = _Next;
            }
            if (_Previous == nullptr)
            {
                Previous = this;
            }
            else
            {
                Previous = _Previous;
            }
            if (i != 0)
            {
                i--;
                if (Type == Right)
                {
                    TopLeft = gcnew CenterSquare(_x1 - h1, _x2 + h2, List::Dawn, nullptr, nullptr, i, h1 / 2.0, h2 / 2.0);
                    TopRight = gcnew CenterSquare(_x1 + h1, _x2 + h2, List::Right, TopLeft, nullptr, i, h1 / 2.0, h2 / 2.0);
                    DawnRight = gcnew CenterSquare(_x1 + h1, _x2 - h2, List::Right, TopRight, nullptr, i, h1 / 2.0, h2 / 2.0);
                    DawnLeft = gcnew CenterSquare(_x1 - h1, _x2 - h2, List::Top, DawnRight, nullptr, i, h1 / 2.0, h2 / 2.0);
                    TopLeft->SetPrevious(TopRight);
                    TopRight->SetPrevious(DawnRight);
                    DawnRight->SetPrevious(DawnLeft);
                    Start = DawnLeft;
                    End = TopLeft;
                }
                if (Type == Top)
                {
                    DawnRight = gcnew CenterSquare(_x1 + h1, _x2 - h2, List::Left, nullptr, nullptr, i, h1 / 2.0, h2 / 2.0);
                    TopRight = gcnew CenterSquare(_x1 + h1, _x2 + h2, List::Top, DawnRight, nullptr, i, h1 / 2.0, h2 / 2.0);
                    TopLeft = gcnew CenterSquare(_x1 - h1, _x2 + h2, List::Top, TopRight, nullptr, i, h1 / 2.0, h2 / 2.0);
                    DawnLeft = gcnew CenterSquare(_x1 - h1, _x2 - h2, List::Right, TopLeft, nullptr, i, h1 / 2.0, h2 / 2.0);
                    DawnRight->SetPrevious(TopRight);
                    TopRight->SetPrevious(TopLeft);
                    TopLeft->SetPrevious(DawnLeft);
                    Start = DawnLeft;
                    End = DawnRight;
                }
                if (Type == Left)
                {
                    DawnRight = gcnew CenterSquare(_x1 + h1, _x2 - h2, List::Top, nullptr, nullptr, i, h1 / 2.0, h2 / 2.0);
                    DawnLeft = gcnew CenterSquare(_x1 - h1, _x2 - h2, List::Left, DawnRight, nullptr, i, h1 / 2.0, h2 / 2.0);
                    TopLeft = gcnew CenterSquare(_x1 - h1, _x2 + h2, List::Left, DawnLeft, nullptr, i, h1 / 2.0, h2 / 2.0);
                    TopRight = gcnew CenterSquare(_x1 + h1, _x2 + h2, List::Dawn, TopLeft, nullptr, i, h1 / 2.0, h2 / 2.0);
                    DawnRight->SetPrevious(DawnLeft);
                    DawnLeft->SetPrevious(TopLeft);
                    TopLeft->SetPrevious(TopRight);
                    Start = TopRight;
                    End = DawnRight;
                }
                if (Type == Dawn)
                {
                    TopLeft = gcnew CenterSquare(_x1 - h1, _x2 + h2, List::Right, nullptr, nullptr, i, h1 / 2.0, h2 / 2.0);
                    DawnLeft = gcnew CenterSquare(_x1 - h1, _x2 - h2, List::Dawn, TopLeft, nullptr, i, h1 / 2.0, h2 / 2.0);
                    DawnRight = gcnew CenterSquare(_x1 + h1, _x2 - h2, List::Dawn, DawnLeft, nullptr, i, h1 / 2.0, h2 / 2.0);
                    TopRight = gcnew CenterSquare(_x1 + h1, _x2 + h2, List::Left, DawnRight, nullptr, i, h1 / 2.0, h2 / 2.0);
                    TopLeft->SetPrevious(DawnLeft);
                    DawnLeft->SetPrevious(DawnRight);
                    DawnRight->SetPrevious(TopRight);
                    Start = TopRight;
                    End = TopLeft;
                }
            }
        }
        double Get_x1()
        {
            return x1x2.Key;
        }
        double Get_x2()
        {
            return x1x2.Value;
        }
        List GetType()
        {
            return Type;
        }
        CenterSquare^ GetStart()
        {
            return Start;
        }
        CenterSquare^ GetEnd()
        {
            return End;
        }
        CenterSquare^ GetNext()
        {
            return Next;
        }
        CenterSquare^ GetPrevious()
        {
            return Previous;
        }
        void SetNext(CenterSquare^ _Next)
        {
            Next = _Next;
        }
        void SetPrevious(CenterSquare^ _Previous)
        {
            Previous = _Previous;
        }
    };
    //
    public ref class Hilbert_Curve_2D
    {
    private:
        CenterSquare^ Head;
        CenterSquare^ End;
    public:
        Hilbert_Curve_2D(double a, double b, double c, double d, unsigned short m, List Type)
        {
            Head = gcnew CenterSquare((b - a) / 2.0, (d - c) / 2.0, Type, nullptr, nullptr, 1, (b - a) / 4.0, (d - c) / 4.0);
            End = Head->GetEnd();
            Head = Head->GetStart();
            m--;
            double h1 = (b - a) / 4.0;
            double h2 = (d - c) / 4.0;
            while (m != 0)
            {
                m--;
                h1 = h1 / 2.0;
                h2 = h2 / 2.0;
                CenterSquare^ Curr = Head;
                CenterSquare^ Current = gcnew CenterSquare(Curr->Get_x1(), Curr->Get_x2(), Curr->GetType(), Curr->GetNext(), Curr->GetPrevious(), 1, h1, h2);
                Curr = Curr->GetNext();
                Head = Current->GetStart();
                while (Curr->GetNext() != Curr)
                {
                    CenterSquare^ Current1 = gcnew CenterSquare(Curr->Get_x1(), Curr->Get_x2(), Curr->GetType(), Curr->GetNext(), Curr->GetPrevious(), 1, h1, h2);
                    Current1->GetStart()->SetPrevious(Current->GetEnd());
                    Current->GetEnd()->SetNext(Current1->GetStart());
                    Current = Current1;
                    Curr = Curr->GetNext();
                }
                CenterSquare^ Current1 = gcnew CenterSquare(Curr->Get_x1(), Curr->Get_x2(), Curr->GetType(), Curr->GetNext(), Curr->GetPrevious(), 1, h1, h2);
                Current1->GetStart()->SetPrevious(Current->GetEnd());
                Current->GetEnd()->SetNext(Current1->GetStart());
                Current = Current1;
                End = Current->GetEnd();
            }
        }
        CenterSquare^ GetHead()
        {
            return Head;
        }
        CenterSquare^ GetEnd()
        {
            return End;
        }
    };
    //
    public ref class Extr
    {
    private:
        KeyValuePair<double, double> min_xy;
        Extr^ _Extr;
    public:
        Extr(){}
        void Base_LNA_1_2_Mer_AGP(double a, double b, double epsilon, double r, unsigned short N, unsigned short m, double epsilon_granichnoe, double eta_0, unsigned short global_local_iterations, unsigned short max_iterations, bool mode)
        {
            unsigned short schetchick = 0;
            if (N == 1)
            {
                HINSTANCE load_function = LoadLibrary(L"FUNC.dll");
                typedef double (*sh) (double);
                sh ShekelFunc = (sh)GetProcAddress(load_function, "ShekelFunc");
                typedef double (*shag) (double, double, double, double, double, unsigned short);
                shag Shag = (shag)GetProcAddress(load_function, "Shag");
                //
                KeyValuePair<double, double> start = KeyValuePair<double, double>(a, ShekelFunc(a));
                KeyValuePair<double, double> end = KeyValuePair<double, double>(b, ShekelFunc(b));
                //
                Otrezok^ otrezok = gcnew Otrezok(start, end, N);
                //
                Otrezki^ interval = gcnew Otrezki(otrezok, r, N);
                //
                std::pair<double, double> pred_i_sled_shag = std::pair<double, double>(a, b);
                //
                while ((abs(pred_i_sled_shag.second - pred_i_sled_shag.first) > epsilon && max_iterations == 0) || (schetchick != max_iterations))
                {
                    min_xy = KeyValuePair<double, double>(pred_i_sled_shag.second, ShekelFunc(pred_i_sled_shag.second));
                    //
                    pred_i_sled_shag.first = pred_i_sled_shag.second;
                    //
                    pred_i_sled_shag.second = Shag(interval->Getm(), (interval->GetX_Rmax()).Key, (interval->GetX_Rmax()).Value, (interval->GetY_Rmax()).Key, (interval->GetY_Rmax()).Value, N);
                    if (pred_i_sled_shag.second >= b)
                    {
                        srand((unsigned short)time(NULL));
                        pred_i_sled_shag.second = (b - epsilon_granichnoe) + (double)rand() / RAND_MAX * ((b - epsilon_granichnoe / 2.0) - (b - epsilon_granichnoe));
                    }
                    if (pred_i_sled_shag.second <= a)
                    {
                        srand((unsigned short)time(NULL));
                        pred_i_sled_shag.second = (a + epsilon_granichnoe / 2.0) + (double)rand() / RAND_MAX * ((a + epsilon_granichnoe) - (a + epsilon_granichnoe / 2.0));
                    }
                    //
                    KeyValuePair<double, double> promejutochnaya_tochka = KeyValuePair<double, double>(pred_i_sled_shag.second, ShekelFunc(pred_i_sled_shag.second));
                    if (mode == true)
                    {
                        interval->Add(promejutochnaya_tochka, N, eta_0, global_local_iterations, schetchick, a, b, epsilon_granichnoe);
                    }
                    else
                    {
                        interval->Add(promejutochnaya_tochka, N, -1.0, 0, -1, a, b, epsilon_granichnoe);
                    }
                    schetchick++;
                }
                FreeLibrary(load_function);
            }
            else
            {
                HINSTANCE load_function = LoadLibrary(L"FUNC.dll");
                typedef double (*grsh) (double, double);
                grsh GrishaginFunc = (grsh)GetProcAddress(load_function, "GrishaginFunc");
                typedef double (*shag) (double, double, double, double, double, unsigned short);
                shag Shag = (shag)GetProcAddress(load_function, "Shag");
                //
                Hilbert_Curve_2D^ Curve_2D = gcnew Hilbert_Curve_2D(a, b, a, b, m, OC::List::Top);
                Hilbert_Curve_2D^ Curve_2D_PI_Na_Dva = gcnew Hilbert_Curve_2D(a, b, a, b, m, OC::List::Left);
                Hilbert_Curve_2D^ Curve_2D_Minus_PI_Na_Dva = gcnew Hilbert_Curve_2D(a, b, a, b, m, OC::List::Right);
                //
                KeyValuePair<double, double> start = KeyValuePair<double, double>(a, GrishaginFunc(Curve_2D->GetHead()->Get_x1(), Curve_2D->GetHead()->Get_x2()));
                KeyValuePair<double, double> end = KeyValuePair<double, double>(b, GrishaginFunc(Curve_2D->GetEnd()->Get_x1(), Curve_2D->GetEnd()->Get_x2()));
                KeyValuePair<double, double> start_PI_Na_Dva = KeyValuePair<double, double>(a, GrishaginFunc(Curve_2D_PI_Na_Dva->GetHead()->Get_x1(), Curve_2D_PI_Na_Dva->GetHead()->Get_x2()));
                KeyValuePair<double, double> end_PI_Na_Dva = KeyValuePair<double, double>(b, GrishaginFunc(Curve_2D_PI_Na_Dva->GetEnd()->Get_x1(), Curve_2D_PI_Na_Dva->GetEnd()->Get_x2()));
                KeyValuePair<double, double> start_Minus_PI_Na_Dva = KeyValuePair<double, double>(a, GrishaginFunc(Curve_2D_Minus_PI_Na_Dva->GetHead()->Get_x1(), Curve_2D_Minus_PI_Na_Dva->GetHead()->Get_x2()));
                KeyValuePair<double, double> end_Minus_PI_Na_Dva = KeyValuePair<double, double>(b, GrishaginFunc(Curve_2D_Minus_PI_Na_Dva->GetEnd()->Get_x1(), Curve_2D_Minus_PI_Na_Dva->GetEnd()->Get_x2()));
                //
                Otrezok^ otrezok = gcnew Otrezok(start, end, N);
                Otrezok^ otrezok_PI_Na_Dva = gcnew Otrezok(start_PI_Na_Dva, end_PI_Na_Dva, N);
                Otrezok^ otrezok_Minus_PI_Na_Dva = gcnew Otrezok(start_Minus_PI_Na_Dva, end_Minus_PI_Na_Dva, N);
                //
                Otrezki^ interval = gcnew Otrezki(otrezok, r, N);
                Otrezki^ interval_PI_Na_Dva = gcnew Otrezki(otrezok_PI_Na_Dva, r, N);
                Otrezki^ interval_Minus_PI_Na_Dva = gcnew Otrezki(otrezok_Minus_PI_Na_Dva, r, N);
                //
                std::pair<double, double> pred_i_sled_shag = std::pair<double, double>(a, b);
                std::pair<double, double> pred_i_sled_shag_PI_Na_Dva = std::pair<double, double>(a, b);
                std::pair<double, double> pred_i_sled_shag_Minus_PI_Na_Dva = std::pair<double, double>(a, b);
                //
                while ((max(abs(pred_i_sled_shag.second - pred_i_sled_shag.first), max(abs(pred_i_sled_shag_PI_Na_Dva.second - pred_i_sled_shag_PI_Na_Dva.first), abs(pred_i_sled_shag_Minus_PI_Na_Dva.second - pred_i_sled_shag_Minus_PI_Na_Dva.first))) > epsilon && max_iterations == 0) || (schetchick != max_iterations))
                {
                    if (pred_i_sled_shag.second == b)
                    {
                        min_xy = KeyValuePair<double, double>(pred_i_sled_shag.second, min(end.Value, min(end_PI_Na_Dva.Value, end_Minus_PI_Na_Dva.Value)));
                    }
                    else
                    {
                        unsigned int number = pred_i_sled_shag.second / ((b - a) / pow(2.0, m * N));
                        if (number <= pow(2.0, m * N) / 2.0)
                        {
                            CenterSquare^ Curr = Curve_2D->GetHead();
                            while (number != 0)
                            {
                                number--;
                                Curr = Curr->GetNext();
                            }
                            min_xy = KeyValuePair<double, double>(pred_i_sled_shag.second, GrishaginFunc(Curr->Get_x1(), Curr->Get_x2()));
                        }
                        else
                        {
                            CenterSquare^ Curr = Curve_2D->GetEnd();
                            while (pow(2.0, m * N) - number - 1 != 0)
                            {
                                number++;
                                Curr = Curr->GetPrevious();
                            }
                            min_xy = KeyValuePair<double, double>(pred_i_sled_shag.second, GrishaginFunc(Curr->Get_x1(), Curr->Get_x2()));
                        }
                        number = pred_i_sled_shag_PI_Na_Dva.second / ((b - a) / pow(2.0, m * N));
                        if (number <= pow(2.0, m * N) / 2.0)
                        {
                            CenterSquare^ Curr = Curve_2D_PI_Na_Dva->GetHead();
                            while (number != 0)
                            {
                                number--;
                                Curr = Curr->GetNext();
                            }
                            if (GrishaginFunc(Curr->Get_x1(), Curr->Get_x2()) < min_xy.Value)
                            {
                                min_xy = KeyValuePair<double, double>(pred_i_sled_shag_PI_Na_Dva.second, GrishaginFunc(Curr->Get_x1(), Curr->Get_x2()));
                            }
                        }
                        else
                        {
                            CenterSquare^ Curr = Curve_2D_PI_Na_Dva->GetEnd();
                            while (pow(2.0, m * N) - number - 1 != 0)
                            {
                                number++;
                                Curr = Curr->GetPrevious();
                            }
                            if (GrishaginFunc(Curr->Get_x1(), Curr->Get_x2()) < min_xy.Value)
                            {
                                min_xy = KeyValuePair<double, double>(pred_i_sled_shag_PI_Na_Dva.second, GrishaginFunc(Curr->Get_x1(), Curr->Get_x2()));
                            }
                        }
                        number = pred_i_sled_shag_Minus_PI_Na_Dva.second / ((b - a) / pow(2.0, m * N));
                        if (number <= pow(2.0, m * N) / 2.0)
                        {
                            CenterSquare^ Curr = Curve_2D_Minus_PI_Na_Dva->GetHead();
                            while (number != 0)
                            {
                                number--;
                                Curr = Curr->GetNext();
                            }
                            if (GrishaginFunc(Curr->Get_x1(), Curr->Get_x2()) < min_xy.Value)
                            {
                                min_xy = KeyValuePair<double, double>(pred_i_sled_shag_Minus_PI_Na_Dva.second, GrishaginFunc(Curr->Get_x1(), Curr->Get_x2()));
                            }
                        }
                        else
                        {
                            CenterSquare^ Curr = Curve_2D_Minus_PI_Na_Dva->GetEnd();
                            while (pow(2.0, m * N) - number - 1 != 0)
                            {
                                number++;
                                Curr = Curr->GetPrevious();
                            }
                            if (GrishaginFunc(Curr->Get_x1(), Curr->Get_x2()) < min_xy.Value)
                            {
                                min_xy = KeyValuePair<double, double>(pred_i_sled_shag_Minus_PI_Na_Dva.second, GrishaginFunc(Curr->Get_x1(), Curr->Get_x2()));
                            }
                        }
                    }
                    //
                    pred_i_sled_shag.first = pred_i_sled_shag.second;
                    pred_i_sled_shag_PI_Na_Dva.first = pred_i_sled_shag_PI_Na_Dva.second;
                    pred_i_sled_shag_Minus_PI_Na_Dva.first = pred_i_sled_shag_Minus_PI_Na_Dva.second;
                    //
                    pred_i_sled_shag.second = Shag(interval->Getm(), (interval->GetX_Rmax()).Key, (interval->GetX_Rmax()).Value, (interval->GetY_Rmax()).Key, (interval->GetY_Rmax()).Value, N);
                    pred_i_sled_shag_PI_Na_Dva.second = Shag(interval_PI_Na_Dva->Getm(), (interval_PI_Na_Dva->GetX_Rmax()).Key, (interval_PI_Na_Dva->GetX_Rmax()).Value, (interval_PI_Na_Dva->GetY_Rmax()).Key, (interval_PI_Na_Dva->GetY_Rmax()).Value, N);
                    pred_i_sled_shag_Minus_PI_Na_Dva.second = Shag(interval_Minus_PI_Na_Dva->Getm(), (interval_Minus_PI_Na_Dva->GetX_Rmax()).Key, (interval_Minus_PI_Na_Dva->GetX_Rmax()).Value, (interval_Minus_PI_Na_Dva->GetY_Rmax()).Key, (interval_Minus_PI_Na_Dva->GetY_Rmax()).Value, N);
                    if (pred_i_sled_shag.second >= b)
                    {
                        srand((unsigned short)time(NULL));
                        pred_i_sled_shag.second = (b - (b - a) / pow(2.0, m * N) + epsilon_granichnoe) + (double)rand() / RAND_MAX * ((b - epsilon_granichnoe) - (b - (b - a) / pow(2.0, m * N) + epsilon_granichnoe));
                    }
                    if (pred_i_sled_shag.second <= a)
                    {
                        srand((unsigned short)time(NULL));
                        pred_i_sled_shag.second = (a + epsilon_granichnoe) + (double)rand() / RAND_MAX * ((a + (b - a) / pow(2.0, m * N) - epsilon_granichnoe) - (a + epsilon_granichnoe));
                    }
                    if (pred_i_sled_shag_PI_Na_Dva.second >= b)
                    {
                        srand((unsigned short)time(NULL));
                        pred_i_sled_shag_PI_Na_Dva.second = (b - (b - a) / pow(2.0, m * N) + epsilon_granichnoe) + (double)rand() / RAND_MAX * ((b - epsilon_granichnoe) - (b - (b - a) / pow(2.0, m * N) + epsilon_granichnoe));
                    }
                    if (pred_i_sled_shag_PI_Na_Dva.second <= a)
                    {
                        srand((unsigned short)time(NULL));
                        pred_i_sled_shag_PI_Na_Dva.second = (a + epsilon_granichnoe) + (double)rand() / RAND_MAX * ((a + (b - a) / pow(2.0, m * N) - epsilon_granichnoe) - (a + epsilon_granichnoe));
                    }
                    if (pred_i_sled_shag_Minus_PI_Na_Dva.second >= b)
                    {
                        srand((unsigned short)time(NULL));
                        pred_i_sled_shag_Minus_PI_Na_Dva.second = (b - (b - a) / pow(2.0, m * N) + epsilon_granichnoe) + (double)rand() / RAND_MAX * ((b - epsilon_granichnoe) - (b - (b - a) / pow(2.0, m * N) + epsilon_granichnoe));
                    }
                    if (pred_i_sled_shag_Minus_PI_Na_Dva.second <= a)
                    {
                        srand((unsigned short)time(NULL));
                        pred_i_sled_shag_Minus_PI_Na_Dva.second = (a + epsilon_granichnoe) + (double)rand() / RAND_MAX * ((a + (b - a) / pow(2.0, m * N) - epsilon_granichnoe) - (a + epsilon_granichnoe));
                    }
                    //
                    std::pair<double, double>* promejutochnaya_tochka = new std::pair<double, double>(pred_i_sled_shag.second, double());
                    std::pair<double, double>* promejutochnaya_tochka_PI_Na_Dva = new std::pair<double, double>(pred_i_sled_shag_PI_Na_Dva.second, double());
                    std::pair<double, double>* promejutochnaya_tochka_Minus_PI_Na_Dva = new std::pair<double, double>(pred_i_sled_shag_Minus_PI_Na_Dva.second, double());
                    unsigned short flag = 0;
                    CenterSquare^ Curr;
                    CenterSquare^ Curr1;
                    CenterSquare^ Curr2;
                    unsigned int number = (pred_i_sled_shag.second) / ((b - a) / pow(2.0, m * N));
                    pred_i_sled_shag.second = number * ((b - a) / pow(2.0, m * N));
                    if (number <= pow(2.0, m * N) / 2.0)
                    {
                        Curr = Curve_2D->GetHead();
                        while (number != 0)
                        {
                            number--;
                            Curr = Curr->GetNext();
                        }
                        if (schetchick % 10 == 0 && schetchick != 0)
                        {
                            promejutochnaya_tochka->second = promejutochnaya_tochka_PI_Na_Dva->second = promejutochnaya_tochka_Minus_PI_Na_Dva->second = GrishaginFunc(Curr->Get_x1(), Curr->Get_x2());
                        }
                        else
                        {
                            promejutochnaya_tochka->second = GrishaginFunc(Curr->Get_x1(), Curr->Get_x2());
                        }
                    }
                    else
                    {
                        Curr = Curve_2D->GetEnd();
                        while (pow(2.0, m * N) - number - 1 != 0)
                        {
                            number++;
                            Curr = Curr->GetPrevious();
                        }
                        if (schetchick % 10 == 0 && schetchick != 0)
                        {
                            promejutochnaya_tochka->second = promejutochnaya_tochka_PI_Na_Dva->second = promejutochnaya_tochka_Minus_PI_Na_Dva->second = GrishaginFunc(Curr->Get_x1(), Curr->Get_x2());
                        }
                        else
                        {
                            promejutochnaya_tochka->second = GrishaginFunc(Curr->Get_x1(), Curr->Get_x2());
                        }
                    }
                    number = (pred_i_sled_shag_PI_Na_Dva.second) / ((b - a) / pow(2.0, m * N));
                    pred_i_sled_shag_PI_Na_Dva.second = number * ((b - a) / pow(2.0, m * N));
                    if (number <= pow(2.0, m * N) / 2.0)
                    {
                        Curr1 = Curve_2D_PI_Na_Dva->GetHead();
                        while (number != 0)
                        {
                            number--;
                            Curr1 = Curr1->GetNext();
                        }
                        if (schetchick % 10 == 0 && schetchick != 0)
                        {
                            if (GrishaginFunc(Curr1->Get_x1(), Curr1->Get_x2()) < promejutochnaya_tochka->second)
                            {
                                promejutochnaya_tochka->second = promejutochnaya_tochka_PI_Na_Dva->second = promejutochnaya_tochka_Minus_PI_Na_Dva->second = GrishaginFunc(Curr1->Get_x1(), Curr1->Get_x2());
                                flag = 1;
                            }
                        }
                        else
                        {
                            promejutochnaya_tochka_PI_Na_Dva->second = GrishaginFunc(Curr1->Get_x1(), Curr1->Get_x2());
                        }
                    }
                    else
                    {
                        Curr1 = Curve_2D_PI_Na_Dva->GetEnd();
                        while (pow(2.0, m * N) - number - 1 != 0)
                        {
                            number++;
                            Curr1 = Curr1->GetPrevious();
                        }
                        if (schetchick % 10 == 0 && schetchick != 0)
                        {
                            if (GrishaginFunc(Curr1->Get_x1(), Curr1->Get_x2()) < promejutochnaya_tochka->second)
                            {
                                promejutochnaya_tochka->second = promejutochnaya_tochka_PI_Na_Dva->second = promejutochnaya_tochka_Minus_PI_Na_Dva->second = GrishaginFunc(Curr1->Get_x1(), Curr1->Get_x2());
                                flag = 1;
                            }
                        }
                        else
                        {
                            promejutochnaya_tochka_PI_Na_Dva->second = GrishaginFunc(Curr1->Get_x1(), Curr1->Get_x2());
                        }
                    }
                    number = (pred_i_sled_shag_Minus_PI_Na_Dva.second) / ((b - a) / pow(2.0, m * N));
                    pred_i_sled_shag_Minus_PI_Na_Dva.second = number * ((b - a) / pow(2.0, m * N));
                    if (number <= pow(2.0, m * N) / 2.0)
                    {
                        Curr2 = Curve_2D_Minus_PI_Na_Dva->GetHead();
                        while (number != 0)
                        {
                            number--;
                            Curr2 = Curr2->GetNext();
                        }
                        if (schetchick % 10 == 0 && schetchick != 0)
                        {
                            if (GrishaginFunc(Curr2->Get_x1(), Curr2->Get_x2()) < promejutochnaya_tochka->second)
                            {
                                promejutochnaya_tochka->second = promejutochnaya_tochka_PI_Na_Dva->second = promejutochnaya_tochka_Minus_PI_Na_Dva->second = GrishaginFunc(Curr2->Get_x1(), Curr2->Get_x2());
                                flag = 2;
                            }
                        }
                        else
                        {
                            promejutochnaya_tochka_Minus_PI_Na_Dva->second = GrishaginFunc(Curr2->Get_x1(), Curr2->Get_x2());
                        }
                    }
                    else
                    {
                        Curr2 = Curve_2D_Minus_PI_Na_Dva->GetEnd();
                        while (pow(2.0, m * N) - number - 1 != 0)
                        {
                            number++;
                            Curr2 = Curr2->GetPrevious();
                        }
                        if (schetchick % 10 == 0 && schetchick != 0)
                        {
                            if (GrishaginFunc(Curr2->Get_x1(), Curr2->Get_x2()) < promejutochnaya_tochka->second)
                            {
                                promejutochnaya_tochka->second = promejutochnaya_tochka_PI_Na_Dva->second = promejutochnaya_tochka_Minus_PI_Na_Dva->second = GrishaginFunc(Curr2->Get_x1(), Curr2->Get_x2());
                                flag = 2;
                            }
                        }
                        else
                        {
                            promejutochnaya_tochka_Minus_PI_Na_Dva->second = GrishaginFunc(Curr2->Get_x1(), Curr2->Get_x2());
                        }
                    }
                    if (schetchick % 10 == 0 && schetchick != 0)
                    {
                        if (flag == 0)
                        {
                            promejutochnaya_tochka->first = pred_i_sled_shag.second;
                            number = 0;
                            Curr1 = Curve_2D_PI_Na_Dva->GetHead();
                            while (Curr->Get_x1() != Curr1->Get_x1() && Curr->Get_x2() != Curr1->Get_x2())
                            {
                                number++;
                                Curr1 = Curr1->GetNext();
                            }
                            promejutochnaya_tochka_PI_Na_Dva->first = pred_i_sled_shag_PI_Na_Dva.second = number * ((b - a) / pow(2.0, m * N));
                            number = 0;
                            Curr2 = Curve_2D_Minus_PI_Na_Dva->GetHead();
                            while (Curr->Get_x1() != Curr2->Get_x1() && Curr->Get_x2() != Curr2->Get_x2())
                            {
                                number++;
                                Curr2 = Curr2->GetNext();
                            }
                            promejutochnaya_tochka_Minus_PI_Na_Dva->first = pred_i_sled_shag_Minus_PI_Na_Dva.second = number * ((b - a) / pow(2.0, m * N));
                        }
                        if (flag == 1)
                        {
                            promejutochnaya_tochka_PI_Na_Dva->first = pred_i_sled_shag_PI_Na_Dva.second;
                            number = 0;
                            Curr = Curve_2D->GetHead();
                            while (Curr1->Get_x1() != Curr->Get_x1() && Curr1->Get_x2() != Curr->Get_x2())
                            {
                                number++;
                                Curr = Curr->GetNext();
                            }
                            promejutochnaya_tochka->first = pred_i_sled_shag.second = number * ((b - a) / pow(2.0, m * N));
                            number = 0;
                            Curr2 = Curve_2D_Minus_PI_Na_Dva->GetHead();
                            while (Curr1->Get_x1() != Curr2->Get_x1() && Curr1->Get_x2() != Curr2->Get_x2())
                            {
                                number++;
                                Curr2 = Curr2->GetNext();
                            }
                            promejutochnaya_tochka_Minus_PI_Na_Dva->first = pred_i_sled_shag_Minus_PI_Na_Dva.second = number * ((b - a) / pow(2.0, m * N));
                        }
                        else
                        {
                            promejutochnaya_tochka_Minus_PI_Na_Dva->first = pred_i_sled_shag_Minus_PI_Na_Dva.second;
                            number = 0;
                            Curr = Curve_2D->GetHead();
                            while (Curr2->Get_x1() != Curr->Get_x1() && Curr2->Get_x2() != Curr->Get_x2())
                            {
                                number++;
                                Curr = Curr->GetNext();
                            }
                            promejutochnaya_tochka->first = pred_i_sled_shag.second = number * ((b - a) / pow(2.0, m * N));
                            number = 0;
                            Curr1 = Curve_2D_PI_Na_Dva->GetHead();
                            while (Curr2->Get_x1() != Curr1->Get_x1() && Curr2->Get_x2() != Curr1->Get_x2())
                            {
                                number++;
                                Curr1 = Curr1->GetNext();
                            }
                            promejutochnaya_tochka_PI_Na_Dva->first = pred_i_sled_shag_PI_Na_Dva.second = number * ((b - a) / pow(2.0, m * N));
                        }
                    }
                    if (mode == true)
                    {

                        interval->Add(KeyValuePair<double, double>(promejutochnaya_tochka->first, promejutochnaya_tochka->second), N, eta_0, global_local_iterations, schetchick, a, b, epsilon_granichnoe);
                        interval_PI_Na_Dva->Add(KeyValuePair<double, double>(promejutochnaya_tochka_PI_Na_Dva->first, promejutochnaya_tochka_PI_Na_Dva->second), N, eta_0, global_local_iterations, schetchick, a, b, epsilon_granichnoe);
                        interval_Minus_PI_Na_Dva->Add(KeyValuePair<double, double>(promejutochnaya_tochka_Minus_PI_Na_Dva->first, promejutochnaya_tochka_Minus_PI_Na_Dva->second), N, eta_0, global_local_iterations, schetchick, a, b, epsilon_granichnoe);
                    }
                    else
                    {
                        interval->Add(KeyValuePair<double, double>(promejutochnaya_tochka->first, promejutochnaya_tochka->second), N, -1.0, 0, -1, a, b, epsilon_granichnoe);
                        interval_PI_Na_Dva->Add(KeyValuePair<double, double>(promejutochnaya_tochka_PI_Na_Dva->first, promejutochnaya_tochka_PI_Na_Dva->second), N, -1.0, 0, -1, a, b, epsilon_granichnoe);
                        interval_Minus_PI_Na_Dva->Add(KeyValuePair<double, double>(promejutochnaya_tochka_Minus_PI_Na_Dva->first, promejutochnaya_tochka_Minus_PI_Na_Dva->second), N, -1.0, 0, -1, a, b, epsilon_granichnoe);
                    }
                    schetchick++;
                }
                FreeLibrary(load_function);
            }
        }
        KeyValuePair<double, double> GetMin()
        {
            return min_xy;
        }
    };
};