#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
/* Так делать не стоит из-за конфликтов имен. Пусть у тебя есть функция с именем как у функции, 
 * которая орпеделена в стд. Тогда при вызове узнать, какая из эти двух функций на самом деле 
 * вызывается будет сложно. А вот если бы функцию из стд ты вызывал std::func(), то такой проблемы
 * возникнуть не может */
using namespace std;

#pragma pack(push, 1) // выравнивание размера полей для структуры
struct BMPFileHeader{
    uint16_t file_type{ 0x4D42 };
    uint32_t file_size{0};
    uint16_t reserved1{0};
    uint16_t reserved2{0};
    uint32_t offset_data{0};
};

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
/* Структуру стоило бы перенести в отдельный файл. Определения ее методов в соответствующий cpp */
struct BMP{
    BMPFileHeader fileHeader;
    BMPInfoHeader infoHeader;
    vector<char> matrix;

    void readBMP(const char* path){
        ifstream fs{path, ios_base::binary};
        if (fs){
            fs.read((char*)&fileHeader, sizeof(fileHeader));
            fs.read((char*)&infoHeader, sizeof(infoHeader));
            fs.seekg(fileHeader.offset_data, fs.beg); // перемещаем указатель на место в файле, где записаны пиксели
            fileHeader.file_size = fileHeader.offset_data; // для подсчета размера файла мы пока считаем только размер заголовка
            matrix.resize(infoHeader.width*infoHeader.height*3);
            if (infoHeader.width%4==0){
                fs.read((char*)matrix.data(), matrix.size());
                fileHeader.file_size += (uint32_t)(matrix.size()); // размер файла уже с учетом вектора пикселей
            } else {
                int row_stride = infoHeader.width*3; // для считывания с учетом отступов
                uint32_t new_stride = row_stride;
                while (new_stride%4!=0){
                    new_stride++;
                }
                vector<char> padding_row(new_stride-row_stride);
                for (int y = 0;y<infoHeader.height; ++y){
                    fs.read((char*)(matrix.data() + row_stride * y), row_stride);
                    fs.read((char*)padding_row.data(), padding_row.size());
                }
                fileHeader.file_size += (uint32_t)(matrix.size()) + infoHeader.height*(uint32_t)(padding_row.size());
            }
            fs.close();
        }
    }
    void rotateBMP_l(){
        vector<char> copy_matrix(matrix.size());
        for (int x = 0; x<matrix.size(); x++){
            //int i = x/(3*infoHeader.width);
            //int j = (x%(3*infoHeader.width))/3;
            //int new_j = i;
            //int new_i = infoHeader.width - 1 - j;
            //int coords = x%3 + 3*(new_i*infoHeader.height+new_j);
            //copy_matrix[coords] = matrix[x];
            /* Такие громоздкие выражения стоит разделять на несколько, чтобы было читабельнее */
            copy_matrix[x%3 + 3*(((x%(3*infoHeader.width))/3)*infoHeader.height+infoHeader.height - 1 - x/(3*infoHeader.width))] = matrix[x];
        }
        swap(infoHeader.height, infoHeader.width);
        matrix=copy_matrix;
    }

    void rotateBMP_r(){
        vector<char> copy_matrix(matrix.size());
        for (int x = 0; x<matrix.size(); x++){
            copy_matrix[x%3 + 3*((infoHeader.width - 1 - (x%(3*infoHeader.width))/3)*infoHeader.height+x/(3*infoHeader.width))] = matrix[x];
        }
        swap(infoHeader.height, infoHeader.width);
        matrix = copy_matrix;
    }

    void writeBMP(const char *fname){
        ofstream of{fname, ios_base::binary};
        if (infoHeader.width%4==0){
            of.write((const char*)&fileHeader, sizeof(fileHeader));
            of.write((const char*)&infoHeader, sizeof(infoHeader));
            of.write((const char*)matrix.data(), matrix.size());
        } else {
                int row_stride = infoHeader.width*3; // считаем отступ заново, т.к. размеры bmp могли поменяться при повороте
                uint32_t new_stride = row_stride;
                while (new_stride%4!=0){
                    new_stride++;
                }
                vector<char> padding_row(new_stride - row_stride);
                /* Если у тебя повторяются в коде несколько раз одни и те же несколько строчек, то стоит
                 * оформить их в отдельную функцию. В принципе это может стать приватной вспомогательной функцией
                 * записи загловка */
                of.write((const char*)&fileHeader, sizeof(fileHeader));
                of.write((const char*)&infoHeader, sizeof(infoHeader));
                for (int y = 0; y < infoHeader.height; ++y){
                    of.write((const char*)(matrix.data() + row_stride * y), row_stride);
                    of.write((const char*)padding_row.data(), padding_row.size());
                }
        }
        of.close();
    }

    double gaussianModel(double x, double y, double sigma){
        return 1. / exp(-(x*x+y*y)/(2*sigma*sigma));
    }

    double *generate_coeff(int radius, double sigma){
        double *coeff = new double[sizeof(double)*radius*radius];
        double sum = 0;
        for (int i = 0; i<radius; i++){
            for (int j = 0; j<radius; j++){
                coeff[i*radius+j] = gaussianModel(i-radius, j-radius/2, sigma);
                sum+=coeff[i*radius+j];
            }
        }
        for (int i = 0; i<radius*radius; i++){
            coeff[i] /= sum;
        }
        return coeff;
    }

    void gauss(){
        double sigma = 1.0;
        int radius = 5;
        int b, g, r;
        int coord=0;
        double* coeff = generate_coeff(radius, sigma);
        /* Довольно сильная вложенность, такое стоит как-то разделять. Введением вспомогательной функции
         * или еще как-то */
        for (int i = 0; i < 2*infoHeader.width; i++){
            for (int j = 0; j < 3*infoHeader.height; j+=3){
                b = g = r = 0;
                for (int m = 0; m < radius; m++){
                    for (int n = 0; n < radius; n++){
                        coord = infoHeader.width*(i+m)+j;
                        if (coord+3*n<0 || coord+2+3*n>matrix.size()) continue;
                        b += coeff[m*radius+n] * matrix[coord + 3*n];
                        g += coeff[m*radius+n] * matrix[coord+1 + 3*n];
                        r += coeff[m*radius+n] * matrix[coord+2 + 3*n];
                        //b += coeff[m*radius+n] * matrix[(i+m)*infoHeader.width + (j+n)];
                        //g += coeff[m*radius+n] * matrix[(i+m)*infoHeader.width + (j+1+n)];
                        //r += coeff[m*radius+n] * matrix[(i+m)*infoHeader.width + (j+2+n)];
                        //b += coeff[m * radius + n] * matrix[0 + (i+m)*infoHeader.width + (j + 3*n)];
                        //g += coeff[m * radius + n] * matrix[1 + (i+m)*infoHeader.width + (j + 3*n)];
                        //r += coeff[m * radius + n] * matrix[2 + (i+m)*infoHeader.width + (j + 3*n)];
                    }
                }
                matrix[i*infoHeader.width + j + 0] = b;
                matrix[i*infoHeader.width + j + 1] = g;
                matrix[i*infoHeader.width + j + 2] = r;
            }
        }
    }
};

int main(){
    BMP bmp;
    bmp.readBMP("dwsample-bmp-640.bmp");
    //bmp.rotateBMP_r();
    bmp.gauss();
    bmp.writeBMP("sample2_out.bmp");
    return 0;
}
