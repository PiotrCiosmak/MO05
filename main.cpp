#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>

using namespace std;

int selectMatrix(const vector<vector<double>> &matrix, int n);

void decompositionLU(vector<vector<double>> &U, vector<vector<double>> &L, vector<double> &b);

void fun1(vector<double> &Y, const vector<vector<double>> &L);

void fun2(vector<double> &X, const vector<double> &Y, const vector<vector<double>> &U);

void showVector(const vector<double> &vector, string label);

void showMatrix(const vector<vector<double>> &matrix, string label);

int main()
{
    vector<vector<double>> A{{1.0,   -20.0,  30.0,   -4.0},
                             {2.0,   -40.0,  -6.0,   50.0},
                             {9.0,   -180.0, 11.0,   -12.0},
                             {-16.0, 15.0,   -140.0, 13.0}};

    std::vector<std::vector<double>> L(A.size(), std::vector<double>(A.size(), 0.0));
    vector<vector<double>> U{A};
    vector<double> b{35.0, 104.0, -366.0, -354.0};
    vector<double> X(4);

    showMatrix(A, "MATRIX A");
    showVector(b, "WEKTOR b");

    decompositionLU(U, L, b);

    showMatrix(L, "MATRIX L");
    showMatrix(U, "MATRIX U");
    showVector(b, "KONCOWY WEKTOR b");

    vector<double> Y{b};
    fun1(Y, L);
    showVector(Y, "KONCOWY WEKTOR Y");
    fun2(X, Y, U);
    showVector(X, "WEKTOR X ROZWIAZANIE");
}

int selectMatrix(const vector<vector<double>> &matrix, int n)
{
    int tmp;
    for (int i = n + 1; i < matrix.size(); ++i)
    {
        if (fabs(matrix[n][i]) < fabs(matrix[n][i + 1]))
            tmp = i + 1;
        else
            tmp = i;
    }
    return tmp;
}

void decompositionLU(vector<vector<double>> &U, vector<vector<double>> &L, vector<double> &b)
{
    unsigned long long size{U.size()};
    double tmpElement;
    int maxElementIndex;
    double tmp;

    for (int i = 0; i < size; ++i)
    {
        if (U[i][i] == 0.0)
        {
            maxElementIndex = selectMatrix(U, i);
            for (int j = 0; j < size; ++j)
            {
                tmpElement = U[i][j];
                U[i][j] = U[maxElementIndex][j];
                U[maxElementIndex][j] = tmpElement;
                tmpElement = L[i][j];
                L[i][j] = L[maxElementIndex][j];
                L[maxElementIndex][j] = tmpElement;
            }
            tmpElement = b[i];
            b[i] = b[maxElementIndex];
            b[maxElementIndex] = tmpElement;
        }
        for (int j = i; j < size - 1; ++j)
        {
            tmp = U[j + 1][i] / U[i][i];
            for (int k = i; k < size; ++k)
                U[j + 1][k] = U[j + 1][k] - U[i][k] * tmp;
            L[j + 1][i] = tmp;
        }
        L[i][i] = 1;
    }
}

void fun1(vector<double> &Y, const vector<vector<double>> &L)
{
    unsigned long long size{L.size()};
    for (int i = 1; i < size; ++i)
        for (int j = i; j < size; ++j)
            Y[j] -= Y[i - 1] * L[j][i - 1];
}

void fun2(vector<double> &X, const vector<double> &Y, const vector<vector<double>> &U)
{
    unsigned long long size{U.size() - 1};
    X[size] = Y[size] / U[size][size];
    double el_Y = 0;
    for (int i = size - 1; i >= 0; --i)
    {
        for (int j = size; j > i; --j)
            el_Y = el_Y + X[j] * U[i][j];
        X[i] = (Y[i] - el_Y) / U[i][i];
        el_Y = 0;
    }
}

void showVector(const vector<double> &vector, string label)
{
    cout << label << endl;
    for (auto element: vector)
        cout << "|" << setw(10) << element << "|" << endl;
    cout << endl;
    cout << endl;
}

void showMatrix(const vector<vector<double>> &matrix, string label)
{
    cout << label << endl;
    for (auto row: matrix)
    {
        for (auto element: row)
            cout << "|" << setw(10) << element;
        cout << "|" << endl;
    }
    cout << endl;
    cout << endl;
}