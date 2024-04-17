#include <cmath>
#include <fstream>
#include <iostream>

const int step = 10;
const int N = 1000, K_x = 100, K_y = 100; // K - размер пространства
const int P_x = 10, P_y = 10; // размеры скорости
const int a = 20, c_x = 50, c_y = 50; // c-положение пика
const int p0_x = 5, p0_y = 5; // начальные скорости

const double dv_x = 5.0/P_x, dv_y= 5.0/P_y;

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

void print(double* data) {
    for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x <= P_x; p_x++) {for(int k_y = 0; k_y < K_y; k_y++) {for(int k_x = 0; k_x < K_x; k_x++) {
    {
        if (data[index(k_x,k_y,p_x,p_y)] > 0.5) std::cout << k_x << " " << k_y << " " << p_x << " " << p_y << " " << data[index(k_x,k_y,p_x,p_y)] << "\n";
    }
} } } }
}

void make_iter_x(double* next, double* data) {
    for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x <= P_x; p_x++) {for(int k_y = 0; k_y < K_y; k_y++) {for(int k_x = 0; k_x < K_x; k_x++) {
    {
        double gam = 1.* p_x /P_x;
        if (p_x > 0) {
            if (k_x == 0) continue; // граничное условие
            else {
                next[index(k_x,k_y,p_x,p_y)] = data[index(k_x,k_y,p_x,p_y)] - gam * (data[index(k_x,k_y,p_x,p_y)] - data[index(k_x - 1,k_y,p_x,p_y)]);
            }
        } else {
            if (k_x == K_x-1) continue;
            else {
                next[index(k_x,k_y,p_x,p_y)] = data[index(k_x,k_y,p_x,p_y)] - gam * (data[index(k_x+1,k_y,p_x,p_y)] - data[index(k_x,k_y,p_x,p_y)]);
            }
        }   
    }
} } } }
    //граничное уловие
    for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x <= P_x; p_x++) {for(int k_y = 0; k_y < K_y; k_y++) {
    {
        if (p_x > 0) next[index(0,k_y,p_x,p_y)] = data[index(0,k_y,p_x,p_y)];
        else next[index(K_x-1,k_y,p_x,p_y)] = data[index(K_x-1,k_y,p_x,p_y)];
    }
} } }
}

void make_iter_y(double* next, double* data) {
    for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x <= P_x; p_x++) {for(int k_y = 0; k_y < K_y; k_y++) {for(int k_x = 0; k_x < K_x; k_x++) {
    {
        double gam = 1.* p_y /P_y;
        if (p_y > 0) {
            if (k_y == 0) continue; // граничное условие
            else {
                next[index(k_x,k_y,p_x,p_y)] = data[index(k_x,k_y,p_x,p_y)] - gam * (data[index(k_x,k_y,p_x,p_y)] - data[index(k_x,k_y-1,p_x,p_y)]);
            }
        } else {
            if (k_y == K_y-1) continue;
            else {
                next[index(k_x,k_y,p_x,p_y)] = data[index(k_x,k_y,p_x,p_y)] - gam * (data[index(k_x,k_y+1,p_x,p_y)] - data[index(k_x,k_y,p_x,p_y)]);
            }
        }   
    }
} } } }
    //граничное уловие
    for (int p_y = -P_y; p_y <= P_y; p_y++) {for (int p_x = -P_x; p_x <= P_x; p_x++) {for(int k_x = 0; k_x < K_x; k_x++) {
    {
        if (p_x > 0) next[index(k_x,0,p_x,p_y)] = data[index(k_x,0,p_x,p_y)];
        else next[index(k_x,K_y-1,p_x,p_y)] = data[index(k_x,K_y-1,p_x,p_y)];
    }
} } }
}

int main() {
    auto size = K_x*K_y*(2*P_x+1)*(2*P_y+1);
    auto data = new double[size];
    auto next = new double[size];
    
    set_initials(data);

    make_iter_x(data, next);
    std::swap(data, next);
    make_iter_y(data, next);
    std::swap(data, next);

    print(data);
    
    return 0;
}