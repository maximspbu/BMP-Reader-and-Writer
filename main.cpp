#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
/* Так не надо: https://stackoverflow.com/questions/4103086/defining-constants-in-c*/
# define M_PI           3.14159265358979323846
/* Так делать не стоит из-за конфликтов имен. Пусть у тебя есть функция с именем как у функции, 
 * которая орпеделена в стд. Тогда при вызове узнать, какая из эти двух функций на самом деле 
 * вызывается будет сложно. А вот если бы функцию из стд ты вызывал std::func(), то такой проблемы
 * возникнуть не может */
using namespace std;
/* Следует разбивать программу на файлы. Во первых, просто ради
 * организации кода. Во вторых, чтобы была возможность проводить
 * раздельную компиляцию. Стоит помещать каждый класс\структуру 
 * в отдельный .h файл, а определения их методов в .cpp. Тогда
 * при изменении реализации какой-то функции придется перекомпили-
 * ровать только этот .cpp файл. Если определение было в заголовке,
 * то придется перекомпилировать все файлы, в которые включен заголо-
 * вок (в крупном проекте это может сократить время компиляции на
 * несколько часов) */
#pragma pack(push, 1) // выравнивание размера полей для структуры
struct BMPFileHeader{
    uint16_t file_type{ 0x4D42 };
    uint32_t file_size{0};
    uint16_t reserved1{0};
    uint16_t reserved2{0};
    uint32_t offset_data{0};
};
#pragma pack(pop)
#pragma pack(push, 1)
struct BMPInfoHeader{
    uint32_t size{0};
    int32_t width{0};
    int32_t height{0};
    uint16_t planes{1};
    uint16_t bit_count{0};
    uint32_t compression{0};
    uint32_t size_image{0};
    int32_t x_pixels_per_meter{0};
    int32_t y_pixels_per_meter{0};
    uint32_t colors_used{0};
    uint32_t colors_important{0};
};
#pragma pack(pop)
/* Стоит сделать классом */
struct BMP{
    BMPFileHeader file_header;
    BMPInfoHeader info_header;
    vector<uint8_t> matrix;
    vector<uint8_t> padding(int* row_stride){
        *row_stride = info_header.width*3; // для считывания с учетом отступов
        uint32_t new_stride = *row_stride;
        while (new_stride%4!=0){
            new_stride++;
        }
        vector<uint8_t> padding_row(new_stride-*row_stride);
        return padding_row;
    }

    void read_bmp(const char* path){
        ifstream fs{path, ios_base::binary};
        if (fs.is_open()){
            fs.read((char*)&file_header, sizeof(file_header));
            fs.read((char*)&info_header, sizeof(info_header));
             /* Между заголовком и пикселями остались данные, которые ты не прочитала. 
              * Поэтому на моем компьютере твоя программа вовсе не работает. Эти данные
              * следовало считать, чтобы потом без изменений записать. А то потеряли*/
            fs.seekg(file_header.offset_data, fs.beg); // перемещаем указатель на место в файле, где записаны пиксели
            file_header.file_size = file_header.offset_data; // для подсчета размера файла мы пока считаем только размер заголовка
            matrix.resize(info_header.width*info_header.height*3);
            if (info_header.width%4==0){
                fs.read((char*)matrix.data(), matrix.size());
                file_header.file_size += (uint32_t)(matrix.size()); // размер файла уже с учетом вектора пикселей
            } else {
                int row_stride = 0;
                vector<uint8_t> padding_row = padding(&row_stride);
                for (int y = 0;y<info_header.height; ++y){
                    fs.read((char*)(matrix.data() + row_stride * y), row_stride);
                    fs.read((char*)padding_row.data(), padding_row.size());
                }
                file_header.file_size += (uint32_t)(matrix.size()) + info_header.height*(uint32_t)(padding_row.size());
            }
            fs.close();
        } else {
        cout << "No such file found" << endl;
        exit(0);
        }
    }
    void rotate_bmp_l(){
        vector<uint8_t> copy_matrix(matrix.size());
        int i, j, new_j, coords;
        for (int x = 0; x<matrix.size(); x++){
            i = x/(3*info_header.width);
            j = (x%(3*info_header.width))/3;
            new_j = info_header.height - 1 - i;
            coords = x%3 + 3*(j*info_header.height+new_j);
            copy_matrix[coords] = matrix[x];
        }
        matrix=copy_matrix;
        swap(info_header.height, info_header.width);
    }

    void rotate_bmp_r(){
        vector<uint8_t> copy_matrix(matrix.size());
        int i, j, new_i, coords;
        for (int x = 0; x<matrix.size(); x++){
            i = x/(3*info_header.width);
            j = (x%(3*info_header.width))/3;
            new_i = info_header.width - 1 - j;
            coords = x%3 + 3*(new_i*info_header.height+i);
            copy_matrix[coords] = matrix[x];
        }
        matrix = copy_matrix;
        swap(info_header.height, info_header.width);
    }

    void write_bmp(const char *fname){
        ofstream of{fname, ios_base::binary};
        if (info_header.width%4==0){
            of.write((const char*)&file_header, sizeof(file_header));
            of.write((const char*)&info_header, sizeof(info_header));
            of.write((const char*)matrix.data(), matrix.size());
        } else {
                int row_stride = 0;
                vector<uint8_t> padding_row = padding(&row_stride);
                of.write((const char*)&file_header, sizeof(file_header));
                of.write((const char*)&info_header, sizeof(info_header));
                for (int y = 0; y < info_header.height; ++y){
                    of.write((const char*)(matrix.data() + row_stride * y), row_stride);
                    of.write((const char*)padding_row.data(), padding_row.size());
                }
        }
        of.close();
    }

    double gaussianModel(double x, double y, double sigma){
        return 1. / exp(-(x*x+y*y)/(2*sigma*sigma));
    }

    vector<double> generate_coeff2(int radius, double sigma){
        vector<double> coeff(sizeof(double)*radius*radius);
        double sum = 0;
        int x;
        for (int i = -radius/2; i<radius/2+1; i++){
            for (int j = -radius/2; j<radius/2+1; j++){
                coeff[(i+radius/2)*radius+j+radius/2] = exp(-(i*i+j*j)/(2*sigma*sigma));
                sum+=coeff[(i+radius/2)*radius+j+radius/2];
            }
        }
        for (int i = 0; i<radius*radius; i++){
            coeff[i] /= sum;
        }
        return coeff;
    }

    void gauss(){
        vector<uint8_t> copy_matrix(matrix.size());
        double sigma = 2.0;
        int radius = 5;
        double b, g, r;
        vector<double> coeff = generate_coeff2(radius, sigma);
         /* Получилась довольно сильная вложенность циклов. Такое следует разбивать, возможно с использованием
          * вспомогательной функции. И вероятно, у нее будет много параметров. Однако нужно стремиться к тому,
          * чтобы у функции было не больше 5 параметров. Так что какие-то придется объединить в одну структуру */
        for (int i = radius/2; i < (info_header.height-radius/2); i++){
            for (int j = (radius/2); j < (info_header.width-radius/2); j++){
                b = g = r = 0;
                for (int m = -radius/2; m < radius/2+1; m++){
                    for (int n = -radius/2; n < radius/2+1; n++){
                        b += coeff[(m+radius/2) * radius + n+radius/2] * matrix[0 + 3*((i+m)*info_header.width + (j+n))];
                        g += coeff[(m+radius/2) * radius + n+radius/2] * matrix[1 + 3*((i+m)*info_header.width + (j+n))];
                        r += coeff[(m+radius/2) * radius + n+radius/2] * matrix[2 + 3*((i+m)*info_header.width + (j+n))];
                    }
                }
                copy_matrix[3*(i*info_header.width + j) + 0] = (uint8_t)b;
                copy_matrix[3*(i*info_header.width + j) + 1] = (uint8_t)g;
                copy_matrix[3*(i*info_header.width + j) + 2] = (uint8_t)r;
            }
        }
        matrix = copy_matrix;
    }
};

int main(){
    BMP bmp;
    bmp.read_bmp("test3.bmp");
    bmp.gauss();
    bmp.write_bmp("sample2_out.bmp");
    return 0;
}