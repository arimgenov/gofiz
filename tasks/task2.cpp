#include <cmath>
#include <fstream>
#include <iostream>

const int step = 10;
const int N = 1000, K_x = 100, K_y = 1; // K - размер пространства
const int P_x = 10, P_y = 0; // размеры скорости
const int a = 20, c_x = 50, c_y = 0; // c-положение пика
const int p0_x = 5, p0_y = 0; // начальные скорости

const double dv_x = 5.0/P_x, dv_y= 1;

int index(int k_x, int k_y, int p_x, int p_y) {             
    return (p_y + P_y) * (2 * P_x + 1) * K_x * K_y
    + (p_x + P_x) * K_x * K_y
    + k_y * K_x
    + k_x;
}

void set_initials(double* data) {
    for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x <= P_x; p_x++) {for(int k_y = 0; k_y < K_y; k_y++) {for(int k_x = 0; k_x < K_x; k_x++) {
    {
        if (p_x == p0_x and p_y == p0_y) data[index(k_x,k_y,p_x,p_y)] = std::exp(-1. * ((k_x-c_x)*(k_x-c_x) + (k_y-c_y)*(k_y-c_y))/(1.*a*a));
        else data[index(k_x,k_y,p_x,p_y)] = 1e-9;
    }
} } } }
}

double calculate_denom() {
    double denom = 0;
    for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = 0; p_x <= P_x; p_x++) {
    {
        denom += std::exp(- (p_x * dv_x) * (p_x * dv_x) / 2.) * p_x;
    }
} }
    return denom;
}

double denom = calculate_denom();

void make_iter_x(double* next, double* data) {
    for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x <= P_x; p_x++) {for(int k_y = 0; k_y < K_y; k_y++) {for(int k_x = 0; k_x < K_x; k_x++) {
    {
        double gam = 1.* p_x /P_x;
        if (p_x > 0) {
            if (k_x == 0) continue; // граничное условие
            else next[index(k_x,k_y,p_x,p_y)] = data[index(k_x,k_y,p_x,p_y)] - gam * (data[index(k_x,k_y,p_x,p_y)] - data[index(k_x - 1,k_y,p_x,p_y)]);
        } else {
            if (k_x == K_x-1) continue;
            else next[index(k_x,k_y,p_x,p_y)] = data[index(k_x,k_y,p_x,p_y)] - gam * (data[index(k_x+1,k_y,p_x,p_y)] - data[index(k_x,k_y,p_x,p_y)]);
        }   
    }
} } } }

    //зеркальное граничное уловие
//     for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x <= P_x; p_x++) {for(int k_y = 0; k_y < K_y; k_y++) {
//     {
//         if (p_x > 0) next[index(0,k_y,p_x,p_y)] = next[index(0,k_y,-p_x,p_y)];
//         else next[index(K_x-1,k_y,p_x,p_y)] = data[index(K_x-1,k_y,-p_x,p_y)];
//     }
// } } }

    // диффузное граничное условие
    double nom0 = 0;
    double nomK = 0;

    for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x < 0; p_x++) {
    {
        nom0 += next[index(0, 0, p_x, 0)] * p_x;
    }
} }
    for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = 0; p_x <= P_x; p_x++) {
    {
        nomK += next[index(K_x-1, 0, p_x, 0)] * p_x;
    }
} }

    for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x <= P_x; p_x++) {for(int k_y = 0; k_y < K_y; k_y++) {
    {
        if (p_x > 0) next[index(0,k_y,p_x,p_y)] = - nom0 / denom * std::exp(- (p_x * dv_x) * (p_x * dv_x) / 2.);
        else next[index(K_x-1,k_y,p_x,p_y)] = nomK / denom * std::exp(- (p_x * dv_x) * (p_x * dv_x) / 2.);
    }
} } }
}

void write_to_file(double* data, char* filename) {
    std::ofstream out(filename);
    
    for (int k_y = 0; k_y < K_y; k_y++) {
        for (int k_x = 0; k_x < K_x; k_x++) {
            double con = 0;
            
            for (int p_y = -P_y; p_y <= P_y; p_y++) {
                for (int p_x = -P_x; p_x <= P_x; p_x++) {
                    con += data[index(k_x,k_y,p_x,p_y)];
            } }

            out << k_x << " " << dv_x * dv_y * con << std::endl;
    } }

    out.close();
}

int main() {
    auto size = K_x*K_y*(2*P_x+1)*(2*P_y+1);
    auto data = new double[size];
    auto next = new double[size];
    
    set_initials(data);

    for(int i=0; i<N; i++) {
        make_iter_x(next, data);
        std::swap(next, data);

        if (i % step == 0) {
            char filename[30];
            sprintf(filename, "D:/VisCode/gofiz/data/out_%03d.dat", i); // запись в файл
            write_to_file(data, filename);
        }
    }
    
    return 0;
}