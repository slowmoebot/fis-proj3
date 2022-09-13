#pragma once
#include <iostream>
#include <string>
#include <cassert>
#include <fstream>
#include <vector>
#include <math.h>

typedef std::vector<double> vec;

struct mat
{
    std::string mat_type;
    int N;
    int M;
    int nnz;
    std::vector<int> I;
    std::vector<int> J;
    std::vector<double> val;
};

struct mat_coo
{
    int N;
    int M;
    std::vector<int> I;
    std::vector<int> J;
    std::vector<double> val;
};

struct entry
{
    int row;
    int col;
    double val;
};

void disp_mat(mat M)
{
    printf("N=%d, M=%d, nnz=%d\n", M.N, M.M, M.nnz);
    printf("Idx ");
    for (auto &&i : M.I)
        printf("%d, ", i);

    printf("\nCol Vec:\n");
    for (auto &&j : M.J)
        printf("%d, ", j);

    printf("\nVal Vec:\n");
    for (auto &&e : M.val)
        printf("% .3e, ", e);

    printf("\n\n\n");
}

void disp_vec(vec v)
{
    int sz = v.size();
    if (sz > 20)
    {
        printf("vec: [%e, %e, %e, ... ,%e, %e, %e]\n", v[0], v[1], v[2], v[sz - 3], v[sz - 2], v[sz - 1]);
    }
    else
    {
        printf("vec: [");
        for (size_t i = 0; i < sz - 1; i++)
        {
            printf("%e, ", v[i]);
        }
        printf("%e]\n", v[sz - 1]);
    }
}

mat read_mat(std::string path)
{
    mat matrix;
    std::ifstream fin(path);
    fin >> matrix.mat_type;
    fin >> matrix.N >> matrix.nnz;

    matrix.M = matrix.N;

    matrix.I.resize(matrix.N + 1, 0);
    for (size_t i = 0; i < matrix.N + 1; i++)
    {
        fin >> matrix.I[i];
        matrix.I[i]--;
    }

    matrix.J.resize(matrix.nnz, 0);
    matrix.val.resize(matrix.nnz, 0);
    for (size_t i = 0; i < matrix.nnz; i++)
    {
        fin >> matrix.J[i] >> matrix.val[i];
        matrix.J[i]--;
    }

    fin.close();
    return matrix;
}

double dot(vec a, vec b)
{
    int sz = a.size();
    assert(b.size() == sz);
    double res = 0;
    for (size_t i = 0; i < sz; i++)
    {
        res += a[i] * b[i];
    }
    return res;
}

vec scale(vec v, double a)
{
    for (size_t i = 0; i < v.size(); i++)
    {
        v[i] *= a;
    }
    return v;
}

vec vec_add(vec v, vec w)
{
    int sz = v.size();
    assert(sz == w.size());

    for (size_t i = 0; i < sz; i++)
    {
        v[i] += w[i];
    }
    return v;
}

double norm(vec x)
{
    int sz = x.size();
    double res = 0;
    for (size_t i = 0; i < sz; i++)
    {
        res += pow(x[i], 2.0);
    }
    return sqrt(res);
}

vec normalize(vec x)
{
    int sz = x.size();
    double x_abs = norm(x);
    for (size_t i = 0; i < sz; i++)
    {
        x[i] /= x_abs;
    }
    return x;
}

vec vecprod_csr(mat matrix, vec vec)
{
    int sz = vec.size();
    assert(matrix.N == sz);

    std::vector<double> res(matrix.M, 0);

    for (size_t i = 0; i < matrix.M; i++)
    {
        int i1 = matrix.I[i];
        int i2 = matrix.I[i + 1];
        int temp_len = i2 - i1;
        std::vector<double> temp1(temp_len);
        std::vector<double> temp2(temp_len);
        for (size_t j = i1; j < i2; j++)
        {
            temp1[j - i1] = matrix.val[j];
            temp2[j - i1] = vec[matrix.J[j]];
        }
        double temp_res = dot(temp1, temp2);
        res[i] = temp_res;
    }

    return res;
}

std::vector<double> vecprod_csc_t(mat matrix, std::vector<double> vec)
{
    size_t sz = vec.size();

    assert(sz == matrix.N);

    std::vector<double> res(matrix.M, 0);

    for (size_t i = 0; i < matrix.M; i++)
    {
        int i1 = matrix.I[i];
        int i2 = matrix.I[i + 1];
        for (size_t j = i1; j < i2; j++)
        {
            if (i != matrix.J[j])
                res[matrix.J[j]] += matrix.val[j] * vec[i];
        }
    }
    return res;
}

std::vector<double> vecprod(mat matrix, std::vector<double> vec)
{
    // Get result from vector product in compressed sparse row format

    std::vector<double> res = vecprod_csr(matrix, vec);

    // If matrix is symmetric the transpose of the matrix needs to be multiplied as well
    if (matrix.mat_type == "s")
    {

        // Get result from transposed vector product without diagonal
        std::vector<double> temp = vecprod_csc_t(matrix, vec);

        // Add results
        for (size_t i = 0; i < matrix.M; i++)
            res[i] += temp[i];
    }

    // Return final result
    return res;
}