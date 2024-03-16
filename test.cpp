#include <cmath>
#include <fstream>

const int step = 10;
const int K = 100, N = 1000;
const double g = -0.2;
const int a = 20, c = 50;

double initial(int k) {
    return std::exp(-((k-c)*1.0/a) * ((k-c)*1.0/a)); // значение начального распределения в точке k
}

void set_initials(double* data) {
    for (int k=0; k < K; k++) data[k] = initial(k); // создание начального массива значений
}

double iter(double* data, int k) {
    if (k == K-1) return data[k];
    else return data[k] - (data[k+1]-data[k]) * g; // итерация одного элемента
}

void make_iter(double* next, double* data) {
    for (int k=0; k<K; k++) next[k] = iter(data, k); // итерация массива
}

void write_to_file(double* data, char* filename) {
    std::ofstream out(filename);
    for (int k=0; k<K; k++) out << data[k] << std::endl;
    out.close();
}

int main() {

    auto data = new double[K];
    auto next = new double[K];

    set_initials(data);

    for(int i=0; i<N; i++) {
        make_iter(next, data);
        std::swap(next, data);

        if (i % step == 0) {
            char filename[30];
            sprintf(filename, "D:/VisCode/gofiz/data/out_%03d.dat", i); // запись в файл
            write_to_file(data, filename);
        }
    }

  return 0;
}