#include "matrix.hpp"
#include "eigenvalues.hpp"

void task1()
{
    std::string input_path = "mat/power_test_csr.txt";
    mat A = read_mat(input_path);

    vec q_0;
    q_0.assign(A.N, 1.0);
    q_0 = normalize(q_0);

    std::string output_path = "output/task1_convergence.csv";

    auto [lambda_lanczos, v_lanczos] = power_iteration(A, q_0, 1e-8, output_path);
}

void task2()
{
    std::string input_path = "mat/cg_test_csr.txt";
    mat A = read_mat(input_path);
    double lambda_act = 9.5986080894852857e03;

    vec q_0;
    q_0.assign(A.N, 1.0);
    q_0 = normalize(q_0);

    std::string output_path = "output/task2_powerIt_convergence.csv";
    auto [lambda_lanczos, v_lanczos] = power_iteration(A, q_0, 1e-8, output_path);

    std::vector<int> ms = {30, 50, 75, 100};
    std::vector<double> thresholds = {1e-2, 1e-4, 1e-6, 1e-10};

    std::string out_temp = "output/temp.csv";

    for (size_t i = 0; i < 4; i++)
    {
        std::string out_path = "output/task2_lanczos_convergence_m" + std::to_string(ms[i]) + ".csv";
        auto [lambda_lanczos, v_lanczos] = lanczos(A, q_0, ms[i], thresholds[i], out_path);
    }
}

int main(int argc, char *argv[])
{
    task1();
    task2();
    /*
    assert(argc > 1 && "Path to matrix has to be passed as first argument");

    std::string path = argv[1];

    mat A = read_mat(path);

    std::vector<double> q_0(A.N, 1);

    q_0 = normalize(q_0);
    printf("%e\n", norm(q_0));

    std::string name = "test.txt";

    std::vector<int> ms = {30, 50, 75, 100};
    std::vector<double> errs = {1e-2, 1e-4, 1e-6, 1e-10};

    int m = 50;
    auto [lambda_lanczos, v_lanczos] = lanczos(A, q_0, m, 1e-4, name);

    printf("lanczos lambda = %e\n", lambda_lanczos);
    */
    return 0;
}