#include "matrix.hpp"
#include <vector>

struct krylov_old
{
    int m;
    mat H_bar;
    mat H;
    std::vector<std::vector<double>> V;
};


struct krylov
{

    std::vector<std::vector<double>> H_bar;
    std::vector<std::vector<double>> V;
};

std::vector<double> getKrylov(mat A, std::vector<std::vector<double>> &V)
{
    //std::cout << V[0][0] << std::endl;
    int sz = A.N;
    assert(sz == A.M);
    for (std::vector<double>v : V) assert(sz == v.size()); 

    //std::cout << "hallo" << std::endl;
    int j = V.size();

    std::vector<double> w(sz);
    w = vecprod(A,V[j-1]);
    
    //std::cout << "hallo 2" << std::endl;
    std::vector<double> h(j+1);


    for (size_t i = 0; i < j; i++)
    {
        h[i] = dot(V[i], w);

        //std::cout << i << " " << j << std::endl;
        for (size_t k = 0; k < sz; k++)
        {
            //if (k <= 10) printf("%d %d %.3f %.3f %.3f \n", (int) i, (int) k , w[k], h[i], V[i][k]);
            w[k] -= h[i]*V[i][k];
        }
    }
    h[j] = norm(w);
    //printf("% .3e \n", h[j]);
    std::vector<double> v(sz);

    for (size_t k = 0; k < sz; k++) v[k]=w[k]/h[j];

    V.push_back(v);
    
    return h;
}

krylov gramSchmidt(mat A, int m, std::vector<double> r0){

    krylov K;
    double r = norm(r0);
    int sz = r0.size();

    std::vector<double> v(sz);
    std::vector<std::vector<double>> H_bar;
    std::vector<std::vector<double>> V;

    for (size_t i = 0; i < sz; i++) v[i] = r0[i]/r;
    V.push_back(v);


    //std::vector<double> h(m+1);
    for (size_t i = 0; i < m; i++)
    {
        std::vector<double> h_temp(i+2);
        h_temp = getKrylov(A, V);
        H_bar.push_back(h_temp);
        //for (size_t j = 0; j <= j+1; j++) h[j]=h_temp;
    }

    /*
    for (size_t i = 0; i < m; i++)
    {
        printf("i= ");
        for (size_t j = 0; j < m; j++)
        {
            printf("% .3e ",dot(V[i],V[j]));
        }
        printf("\n");
    }
    */

    K.H_bar  = H_bar;
    K.V = V;
    return K;
}

std::vector<double> gmres(mat A, std::vector<double> x0, std::vector<double> b, int m)
{
    int sz = A.N;

    //impelemnt checks for dimenstions

    std::vector<double> v(sz);
    //std::vector<double> r0(sz,0.0);
    std::vector<double> g(m+1,0.0);

    //printf("gmres [1] r0: ");
    //for (auto &&e : r0) printf("% 3.3f, ",e);
    //printf("\n");

    //disp_full(A);

    std::vector<double> r0 = vecprod(A,x0);

    /*
    printf("gmres [2] r0: ");
    for (auto &&e : r0) printf("% 3.3f, ",e);
    printf("\n");
    */

    for (size_t i = 0; i < sz; i++) r0[i] = b[i]-r0[i];

    /*
    printf("gmres [3] r0: ");
    for (auto &&e : r0) printf("% 3.3f, ",e);
    printf("\n");
    
    printf("gmres b: ");
    for (auto &&e : b) printf("% 3.3f, ",e);
    printf("\n");
    */

    double r = norm(r0);
    for (size_t i = 0; i < sz; i++) v[i] = r0[i]/r;

    g[0]=r;

    krylov K;
    std::vector<std::vector<double>> H_bar;
    std::vector<std::vector<double>> V;
    V.push_back(v);

    std::vector<double> c(m);
    std::vector<double> s(m);
    for (size_t j = 0; j < m; j++)
    {
        //printf("m  %d \n", (int) j);
        std::vector<double> h(j+2,0.0);
        h = getKrylov(A, V);
        H_bar.push_back(h);

        for (size_t k_i = 0; k_i <= j; k_i++)
        {
            //printf("H_bar col[%d]: ",(int) k_i);
            for (size_t k_j = 0; k_j <= j+1; k_j++)
            {
                //printf("% 3.3f, ", H_bar[k_i][k_j]);
            }
            //printf("\n");
        }
        


        for (size_t k = 1; k <= j; k++)
        {
            //printf("k=%d \n",(int) k);
            double temp_h  =  c[k-1]*H_bar[j][k-1] + s[k-1]*H_bar[j][k];
            H_bar[j][k] = -s[k-1]*H_bar[j][k-1] + c[k-1]*H_bar[j][k];
            H_bar[j][k-1] = temp_h;
        }
        //printf("Hallo 2\n");

        double temp_sc = sqrt(pow(H_bar[j][j], 2)+pow(H_bar[j][j+1], 2));
        c[j] = H_bar[j][j]  /temp_sc;
        s[j] = H_bar[j][j+1]/temp_sc;
        //printf("j=%d \n",(int) j);
        //printf("H_j_j=% 3.3f, H_j+1_j=% 3.3f \n", H_bar[j][j], H_bar[j][j+1]);
        //printf("c_j=% 3.3f, s_j=% 3.3f \n", c[j], s[j]);
        //printf("sc= % 3.3f \n",temp_sc);

        //printf("Hallo 3\n");
        H_bar[j][j] = c[j]*H_bar[j][j] + s[j]*H_bar[j][j+1];
        g[j+1]= -s[j]*g[j];
        g[j] = c[j]*g[j];
        
        printf("\nm       = % 2d \n",(int)j+1);
        printf(  "abs_err = % 3.3f \n",abs(g[j+1]));
        printf(  "rel_err = % .3e \n",abs(g[j+1])/r);

        /*
        if (abs(g[j+1])/r<err)
        {
            m=j+1;
            break;
        }
        */
    }

    //std::vector<double> res(m);

    //printf("Result: ");
    for (int i = m - 1; i >= 0; i--)
    {
        for (int j = m - 1; j > i; j--)
        {
            g[i] -= g[j]*H_bar[j][i];
        }
        g[i] /= H_bar[i][i];
    }

    std::vector<double> res(sz,0.0);
    for (size_t i = 0; i < m; i++)
    {
        for (size_t k=0; k<sz; k++) x0[k] += g[i]*V[i][k];
        //printf("% .3e, ", g[i]);
    }

    //g.pop_back();

    x0.push_back(abs(g[m]));
    
    return x0;
    //printf("\n");
    
    //printf("\nres:\n");
    //for (double d:x0) printf("% 3.3f, ",d);
    //printf("\n");
}

std::vector<double> restarted_gmres(mat A, std::vector<double> x0, std::vector<double> b, int m)
{
    int sz = A.N;

    std::vector<double> v(sz);
    std::vector<double> g(m+1,0.0);

    std::vector<double> r0 = vecprod(A,x0);

    for (size_t i = 0; i < sz; i++) r0[i] = b[i]-r0[i];
    double r = norm(r0);
    double r_init = r;

    int i = 0;
    while (abs(r) > 1e-8)
    {
        x0 = gmres(A,x0,b,m);
        r=x0[sz]/r_init;
        i++;
        printf("#########\nFinished iteration %d with an error of % .5e \n#########\n",i,r);
        x0.pop_back();

    }
    return x0;
}
