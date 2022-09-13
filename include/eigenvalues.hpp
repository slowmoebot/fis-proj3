#pragma once
#include <iostream>
#include <fstream>
#include <tuple>
#include <chrono>
#include <iomanip>
#include "matrix.hpp"

std::tuple<mat, std::vector<vec>> getKrylov(mat A, vec v, int m)
{
    int sz = A.N;
    assert(A.M == sz);
    assert(v.size() == sz);
    assert(A.mat_type == "s");
    v = normalize(v);
    vec diag;
    vec subdiag;
    diag.reserve(m);
    subdiag.reserve(m);

    std::vector<vec> V;
    V.push_back(v);
    vec w;

    for (size_t j = 0; j < m; j++)
    {
        w = vecprod(A, V[j]);
        if (j > 0)
        {
            w = vec_add(w, scale(V[j - 1], -subdiag[j - 1]));
        }
        diag[j] = dot(w, v);
        w = vec_add(w, scale(v, -diag[0]));
        subdiag[j] = norm(w);
        V.push_back(scale(w, 1 / subdiag[j]));
    }

    mat res;
    res.mat_type = "s";
    res.M = m;
    res.N = m;
    res.nnz = 2 * m - 1;
    res.I.assign(m + 1, 0);
    res.I[0] = 0;
    res.J.assign(res.nnz, 0);
    res.val.assign(res.nnz, 0);

    for (size_t i = 0; i < m; i++)
    {
        res.I[i + 1] = 2 * i + 1;
        if (i < m - 1)
        {
            res.J[2 * i] = i;
            res.J[2 * i + 1] = i;
            res.val[2 * i] = diag[i];
            res.val[2 * i + 1] = subdiag[i];
        }
    }
    res.J[2 * m - 2] = m - 1;
    res.val[2 * m - 2] = diag[m - 1];

    return std::make_tuple(res, V);
}

std::tuple<double, vec> power_iteration(mat A, vec q, double threshold, std::string name = "")
{

    auto start = std::chrono::high_resolution_clock::now();

    assert(norm(q) - 1.0 < 1e-10);
    int n = q.size();
    vec z;

    std::ofstream f;
    f.open(name);

    double lambda_old;
    double lambda = 0;

    double err = 1;
    for (size_t i = 0; i < 1e4; i++)
    {
        z = vecprod(A, q);
        q = normalize(z);

        if (i > 0)
            lambda_old = lambda;

        lambda = dot(q, vecprod(A, q));
        if (name != "")
        {
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
            f << i << " " << duration.count() << " " << std::setprecision(15) << lambda << "\n";
        }

        if (i > 0)
        {
            err = abs(lambda - lambda_old);
            if (err < threshold)
                break;
        }
    }
    f.close();
    return std::make_tuple(lambda, q);
}

std::tuple<double, vec> lanczos(mat A, vec v_0, int m, double err, std::string name)
{
    vec q_0 = normalize(v_0);
    auto [T_m, V_krylov] = getKrylov(A, q_0, m);

    q_0.assign(m, 1);
    q_0 = normalize(q_0);
    return power_iteration(T_m, q_0, err, name);
}