// Diplom.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
#include <iostream>
#include <math.h>
#include <stdlib.h> // нужен для вызова функций rand(), srand()
#include <time.h> // нужен для вызова функции time()
#include <string> // for string and to_string()
#include <fstream>
#include <cmath>
#include <utility> // для std::pair 

#define M_PI 3.14159265358979323846
const int n = 20; // размер сетки скоростей
const int Nx = 75;// размер сетки  координат nx= 3*N, ny = N
const int Ny = 25;
const  double T1 = 1;
const  double T2 = 2;
const  double u = 0.0;
const  double ms = 0.5; // масса моллекулы
// Функция распределения Максвелла через индексы скоростей

double maxwellDistFunc(double iksx, double iksy, double iksz, double kon, double T) {
    // T = 1;
    double v = (iksx - n / 2) + 0.5;
    double u = (iksy - n / 2) + 0.5;
    double w = (iksz - n / 2) + 0.5;
    double u0 = 0.0;
    double mnoj = (ms / (2 * M_PI * T));
    return kon * (pow(mnoj, 1.5) * exp(-ms * ((v - u0) * (v - u0) + u * u + w * w) / (2 * T))); // ix, iy - индексы скоростей
    //нормировка, чтобы суммы на единицу
}

// Функция распределения Максвелла через скорости
double distribution_speed(double v, double u, double w, double kon, double T) {
    // T = 1;
    double mnoj = (ms / (2 * M_PI * T));
    return kon * (pow(mnoj, 1.5) * exp(-ms * (v * v + u * u + w * w) / (2 * T))); // iksx, iksy - индексы скоростей
    //нормировка, чтобы суммы на единицу
}

double distribution_speed1(double iksx, double iksy, double iksz, double ux, double uy, double kon, double T) {
    // T = 1;
    double vx = (iksx - n / 2) + 0.5;
    double vy = (iksy - n / 2) + 0.5;
    double vz = (iksz - n / 2) + 0.5;
    double mnoj = (ms / (2 * M_PI * T));
    return kon * (pow(mnoj, 1.5) * exp(-ms * ((vx - ux) * (vx - ux) + (vy - uy) * (vy - uy) + vz * vz) / (2 * T))); // iksx, iksy - индексы скоростей
    //нормировка, чтобы суммы на единицу
}

// Функция для вычисления плотности
double computeDensity(double f[Nx + 1][Ny + 1][n][n][n], int ix, int iy) {
    double sum = 0;
    for (int iksx = 0; iksx < n; iksx++) {
        for (int iksy = 0; iksy < n; iksy++) {
            for (int iksz = 0; iksz < n; iksz++) {
                sum += f[ix][iy][iksx][iksy][iksz];
            }
        }
    }
    return sum;
}
std::pair<double, double> computeVelocity(double f[Nx + 1][Ny + 1][n][n][n], int ix, int iy) {
    double vx_sum = 0;
    double vy_sum = 0;
    double density = computeDensity(f, ix, iy);
    for (int iksx = 0; iksx < n; iksx++) {
        for (int iksy = 0; iksy < n; iksy++) {
            for (int iksz = 0; iksz < n; iksz++) {
                double vx = (iksx - n / 2) + 0.5;
                double vy = (iksy - n / 2) + 0.5;
                double vz = (iksy - n / 2) + 0.5;
                vx_sum += f[ix][iy][iksx][iksy][iksz] * vx;
                vy_sum += f[ix][iy][iksx][iksy][iksz] * vy;
            }
        }
    }
    return std::pair<double, double>(vx_sum / density, vy_sum / density);
}
// Функция для вычисления температуры
double computeTemperature(double f[Nx + 1][Ny + 1][n][n][n], int ix, int iy) {
    double sum = 0;
    for (int iksx = 0; iksx < n; iksx++) {
        for (int iksy = 0; iksy < n; iksy++) {
            for (int iksz = 0; iksz < n; iksz++) {
                sum += f[ix][iy][iksx][iksy][iksz] * ((iksx - n / 2 + 0.5) * (iksx - n / 2 + 0.5) +
                    (iksy - n / 2 + 0.5) * (iksy - n / 2 + 0.5) +
                    (iksz - n / 2 + 0.5) * (iksz - n / 2 + 0.5));

            }
        }
    }
    return sum * (ms / 3);


}



// Функция для вычисления теплового потока
std::pair<double, double> computeHeatFlux(double f[Nx + 1][Ny + 1][n][n][n], int ix, int iy) {
    double qx_sum = 0;
    double qy_sum = 0;
    for (int iksx = 0; iksx < n; iksx++) {
        for (int iksy = 0; iksy < n; iksy++) {
            for (int iksz = 0; iksz < n; iksz++) {
                double vx = (iksx - n / 2) + 0.5;
                double vy = (iksy - n / 2) + 0.5;
                double vz = (iksz - n / 2) + 0.5;
                double e = 0.5 * ms * (vx * vx + vy * vy + vz * vz); // кинетическая энергия
                qx_sum += f[ix][iy][iksx][iksy][iksz] * vx * e;
                qy_sum += f[ix][iy][iksx][iksy][iksz] * vy * e;
            }
        }
    }
    return  std::pair<double, double> { qx_sum, qy_sum };
}
// Функция для вычисления теплового потока с учетом гидродинамической скорости
std::pair<double, double> computeHeatFluxGydro(double f[Nx + 1][Ny + 1][n][n][n], int ix, int iy, double ux, double uy) {
    double qx_sum = 0;
    double qy_sum = 0;
    for (int iksx = 0; iksx < n; iksx++) {
        for (int iksy = 0; iksy < n; iksy++) {
            for (int iksz = 0; iksz < n; iksz++) {
                double vx = (iksx - n / 2) + 0.5 - ux;
                double vy = (iksy - n / 2) + 0.5 - uy;
                double vz = (iksz - n / 2) + 0.5;
                double e = 0.5 * ms * (vx * vx + vy * vy + vz * vz); // кинетическая энергия
                qx_sum += f[ix][iy][iksx][iksy][iksz] * vx * e;
                qy_sum += f[ix][iy][iksx][iksy][iksz] * vy * e;
            }
        }
    }
    return  std::pair<double, double> { qx_sum, qy_sum };
}


using namespace std;
int main()
{
    int i, j, k, l, m;
    int ix = 0, iy = 0, nx, ny, iksx, iksy, iksz, nksx, nksy, nksz;
    long ns, nt, kt, ks, nsl;
    double s1, s2, intkon, ksx, ksy, ksz, ks1x, ks1y, ks1z, dks = 1.0;
    double koef1 = 16.9/3.0, koef2;
    long ial, jal, kal, ibt, jbt, kbt; // координаты двух частиц начальных
    double bp, ep; // ep = азимутальный угол
    double dks_bolc = 0.001, dt = 0.0005, dx = 0.025, dy = 0.025;
    double dt2 = dt * koef1;//dlya int stolknovenii
    double koef = 1;
    double ff[Nx + 1][Ny + 1][n][n][n] = { 0 }; // Инициализация нулями
    double ff_12[Nx + 1][Ny + 1][n][n][n] = { 0 }; // Инициализация нулями
    double ints[n][n][n]; // интеграл столкновений
    double ints12[n][n][n];
    long ibl, jbl, kbl, i1bl, j1bl, k1bl;
    double rr[20]; // массив для нахождения кратчайшего расстояния, удвоил 22.05
    double rrp;
    double enn[10]; // энергии в этих точках
    double enpr[10]; // еще один массив для энергии всего 8 точек
    double s; // буферная переменная
    long lmbi, lmbj, lmbk, mui, muj, muk, di, dj, dk; // ближайшие точки после столкновения лямбда
    long lmpsi, lmpsj, lmpsk, mumsi, mumsj, mumsk; // вторая дискретная точка мю
    double mnr, xcm, ycm, zcm;
    double en2g, en1g, rgm; // rgamma
    long si, nsL; // переменные для дискретизации
    double ombgm; // omega bolshoe gamma

    double teta;
    double gx, gy, gz, gxy; // gxy = (gx^2 + gy^2)^0.5
    double g;
    double gshx, gshy, gshz;
    double kshx, kshy, kshz; // ks_h_x
    double ksh1x, ksh1y, ksh1z;
    double rpr;
    double gg; // модуль g
    ns = 10000;// /(40);//1200000; //количество точек выбирается вручную для интеграла столкновений
    nt = 60000; // количетсво циклов по времени для интеграла столкновений

    int ibL, jbL, kbL;
    double en0, enp; // enp энергия в ближайшей дискретной точке
    int iii;
    long en0cl, enpcl;
    bool fLr1;
    double ennb;
    double tprod, tpoper;
    string filname = "hello";
    string filname1, filname2, filename3;
    double dt1 = 0.0005;// 0.0005;
    double ron = 0;
    double vxn = 0;
    double vyn = 0;
    // double v_mean[Nx+1][Ny+1][2] = {0};
    double vx, vy;
    double spd_par = 0.5;

    double ss, pz, px, py;

    double vxnp = 0;
    double vynp = 0;

    int npr;
    //=================================26.05==================================================
    double f1excp1[n][n][n];
    double f1excp2[n][n][n];
    double f1excp3[n][n][n];
    double f1excp4[n][n][n];
    double f1excp5[n][n][n];
    double f1excp6[n][n][n];
    double f1excp7[n][n][n];
    double f1excp8[n][n][n];

    double fexcp1[n][n][n];
    double fexcp2[n][n][n];
    double fexcp3[n][n][n];
    double fexcp4[n][n][n];
    double fexcp5[n][n][n];
    double fexcp6[n][n][n];
    double fexcp7[n][n][n];
    double fexcp8[n][n][n];
    //=================================26.05=================================================


    double koefotr = 1;
    double sumpotok;
    double v2;
    double sumkonc;

    filname = "temp";
    ofstream outemp(filname);


    // Установить генератор случайных чисел
    srand(time(0));

    // Инициализация параметров
    nksx = n;
    nksy = n;
    nksz = n;
    nx = Nx;
    ny = Ny;
    int ratio = 5; // соотношение длины чипа к высоте
    int bias = 10; // смещение относительно левого 
    int height = 3; // высота чипа 
    double kon_1 = 1; // концентрация на входе
    double kon_2 = 1;// kon_1 * ny / (ny - height);// раньше было kon_2= kon_1* ny / (ny - height);
    double distribution;
    // Заполнение начальных условий
    for (iksz = 0; iksz < nksz; iksz++) {
        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {

                for (ix = 0; ix < bias; ix++) {
                    for (iy = 0; iy <= ny; iy++) {

                        //  часть трубы левее чипа
                        distribution = maxwellDistFunc(iksx, iksy, iksz, kon_1, T1);
                        ff[ix][iy][iksx][iksy][iksz] = distribution;
                    }
                }

                for (ix = bias; ix <= (bias + height * ratio); ix++) {
                    for (iy = height + 1; iy <= ny; iy++) {

                        // часть над чипом
                        distribution = maxwellDistFunc(iksx, iksy, iksz, kon_2, T1);
                        ff[ix][iy][iksx][iksy][iksz] = distribution;
                    }
                }

                for (ix = (bias + height * ratio + 1); ix <= nx; ix++) {
                    for (iy = 0; iy <= ny; iy++) {

                        //  часть трубы правее чипа
                        distribution = maxwellDistFunc(iksx, iksy, iksz, kon_1, T1);
                        ff[ix][iy][iksx][iksy][iksz] = distribution;
                    }
                }

                for (iy = 0; iy < height; iy++) {
                    // левая стенка чипа
                    distribution = maxwellDistFunc(iksx, iksy, iksz, kon_1, T2);
                    ff[bias][iy][iksx][iksy][iksz] = distribution;
                }
                for (ix = bias; ix <= (bias + height * ratio); ix++) {

                    // верхняя стенка чипа
                    distribution = maxwellDistFunc(iksx, iksy, iksz, kon_2, T2);
                    ff[ix][height][iksx][iksy][iksz] = distribution;
                }
                for (iy = 0; iy < height; iy++) {
                    // правая стенка чипа
                    distribution = maxwellDistFunc(iksx, iksy, iksz, kon_1, T2);
                    ff[bias + height * ratio][iy][iksx][iksy][iksz] = distribution;
                }

            }
        }
    }



    for (iksz = 0; iksz < nksz; iksz++) {
        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                for (ix = 0; ix < bias; ix++) {
                    for (iy = 0; iy <= ny; iy++) {
                        ff_12[ix][iy][iksx][iksy][iksz] = 0;

                    }
                }
            }
        }
    }
    for (iksz = 0; iksz < nksz; iksz++) {
        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                distribution = maxwellDistFunc(iksx, iksy, iksz, kon_1, T2);

                fexcp1[iksx][iksy][iksz] = distribution;
                fexcp2[iksx][iksy][iksz] = distribution;
                fexcp3[iksx][iksy][iksz] = distribution;
                fexcp4[iksx][iksy][iksz] = distribution;
                fexcp5[iksx][iksy][iksz] = distribution;
                fexcp6[iksx][iksy][iksz] = distribution;
                fexcp7[iksx][iksy][iksz] = distribution;
                fexcp8[iksx][iksy][iksz] = distribution;








                f1excp1[iksx][iksy][iksz] = 0;
                f1excp2[iksx][iksy][iksz] = 0;
                f1excp3[iksx][iksy][iksz] = 0;
                f1excp4[iksx][iksy][iksz] = 0;
                f1excp5[iksx][iksy][iksz] = 0;
                f1excp6[iksx][iksy][iksz] = 0;
                f1excp7[iksx][iksy][iksz] = 0;
                f1excp8[iksx][iksy][iksz] = 0;


            }
        }
    }












    for (kt = 0; kt < nt; kt++)
    {
        // массив средних скоростей
      /*  for (i = 0; i <= nx; i++) {
            for (j = 0; j <= ny; j++) {
                pair<double, double> velocity = computeVelocity(ff, ix, iy);
                double vx = velocity.first;
                double vy = velocity.second;
                v_mean[i][j][1] = vx;
                v_mean[i][j][2] = vy;
            }
        }*/
        //вывод в файл
        if (kt % 2 == 0) {
            std::ofstream data_file("box_bias10_koef_5.63_ns10000_nsl_2331.csv", std::ios::app); // Открываем файл в режиме добавления
            if (kt == 0) { // Записываем заголовки только при первом шаге
                data_file << "kt,x,y,temperature,density,vx,vy,heat_flux_x,heat_flux_y,heat_flux_x_gydro,heat_flux_y_gydro,\n";
            }
            for (ix = 0; ix <= Nx; ix++) {
                for (iy = 0; iy <= Ny; iy++) {
                    double density = computeDensity(ff, ix, iy);
                    double temperature = computeTemperature(ff, ix, iy);
                    if (density != 0)
                        temperature = temperature / density;


                    std::pair<double, double> velocity = computeVelocity(ff, ix, iy);
                    std::pair<double, double> heat_flux = computeHeatFlux(ff, ix, iy);
                    double vx = velocity.first;
                    double vy = velocity.second;
                    double heat_flux_x = heat_flux.first;
                    double heat_flux_y = heat_flux.second;
                    std::pair<double, double> heat_flux_gydro = computeHeatFluxGydro(ff, ix, iy, vx, vy);
                    double heat_flux_x_gydro = heat_flux_gydro.first;
                    double heat_flux_y_gydro = heat_flux_gydro.second;
                    data_file << kt << "," << ix << "," << iy << "," << temperature << "," << density << "," << vx << "," << vy << "," << heat_flux_x << "," << heat_flux_y << "," << heat_flux_x_gydro << "," << heat_flux_y_gydro << "\n";
                }
            }
            data_file.close();
            std::cout << "file" << kt << endl;


            sumkonc = 0;

            for (i = 1; i <= bias - 1; i++) {
                for (j = 1; j < ny; j++) {
                    for (iksx = 0; iksx < nksx; iksx++) {
                        for (iksy = 0; iksy < nksy; iksy++) {
                            for (iksz = 0; iksz < nksz; iksz++) {
                                sumkonc = sumkonc + ff_12[i][j][iksx][iksy][iksz];
                            }
                        }
                    }

                }
            }


            for (i = bias; i <= bias + height * ratio; i++) {
                for (j = height + 1; j < ny; j++) {
                    for (iksx = 0; iksx < nksx; iksx++) {
                        for (iksy = 0; iksy < nksy; iksy++) {
                            for (iksz = 0; iksz < nksz; iksz++) {
                                sumkonc = sumkonc + ff_12[i][j][iksx][iksy][iksz];
                            }
                        }
                    }

                }
            }

            for (i = bias + height * ratio + 1; i < nx; i++) {
                for (j = 1; j < ny; j++) {
                    for (iksx = 0; iksx < nksx; iksx++) {
                        for (iksy = 0; iksy < nksy; iksy++) {
                            for (iksz = 0; iksz < nksz; iksz++) {
                                sumkonc = sumkonc + ff_12[i][j][iksx][iksy][iksz];
                            }
                        }
                    }

                }
            }




            cout << sumkonc << endl;




            sumpotok = 0;

            for (j = 1; j <= 5; j++) {


                for (iksx = 0; iksx < nksx; iksx++) {
                    for (iksy = 0; iksy < nksy; iksy++) {
                        for (iksz = 0; iksz < nksz; iksz++) {
                            v2 = 0;
                            v2 = v2 + spd_par * (2 * (iksx - nksx / 2) + 1) * spd_par * (2 * (iksx - nksx / 2) + 1);
                            v2 = v2 + spd_par * (2 * (iksy - nksy / 2) + 1) * spd_par * (2 * (iksy - nksy / 2) + 1);
                            v2 = v2 + spd_par * (2 * (iksz - nksz / 2) + 1) * spd_par * (2 * (iksz - nksz / 2) + 1);
                            sumpotok = sumpotok + ff[bias - 2][iy][iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1) * v2;

                        }
                    }
                }
            }

            for (j = 1; j <= 5; j++) {

                for (iksx = 0; iksx < nksx; iksx++) {
                    for (iksy = 0; iksy < nksy; iksy++) {
                        for (iksz = 0; iksz < nksz; iksz++) {
                            v2 = 0;
                            v2 = v2 + spd_par * (2 * (iksx - nksx / 2) + 1) * spd_par * (2 * (iksx - nksx / 2) + 1);
                            v2 = v2 + spd_par * (2 * (iksy - nksy / 2) + 1) * spd_par * (2 * (iksy - nksy / 2) + 1);
                            v2 = v2 + spd_par * (2 * (iksz - nksz / 2) + 1) * spd_par * (2 * (iksz - nksz / 2) + 1);
                            sumpotok = sumpotok + ff[bias + height * ratio + 2][iy][iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1) * v2;

                        }
                    }
                }


            }

            for (i = bias - 2; i <= bias + height * ratio + 2; i++) {

                for (iksx = 0; iksx < nksx; iksx++) {
                    for (iksy = 0; iksy < nksy; iksy++) {
                        for (iksz = 0; iksz < nksz; iksz++) {
                            v2 = 0;
                            v2 = v2 + spd_par * (2 * (iksx - nksx / 2) + 1) * spd_par * (2 * (iksx - nksx / 2) + 1);
                            v2 = v2 + spd_par * (2 * (iksy - nksy / 2) + 1) * spd_par * (2 * (iksy - nksy / 2) + 1);
                            v2 = v2 + spd_par * (2 * (iksz - nksz / 2) + 1) * spd_par * (2 * (iksz - nksz / 2) + 1);
                            sumpotok = sumpotok + ff[i][5][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1) * v2;

                        }
                    }
                }


            }



            if (outemp.is_open()) {
                outemp << sumpotok << endl;
            }


        }


        // начало  дивергентной разностной схемы

        for (ix = 1; ix < nx; ix++) {
            for (iy = 1; iy < ny; iy++) {
                if (
                    ((iy >= 0 && iy <= ny) && ix == 0) || ((ix >= 0 && ix <= nx) && iy == ny) || ((iy >= 0 && iy <= ny) && ix == nx) ||
                    ((ix >= (bias + height * ratio) && ix <= nx) && iy == 0) || ((ix >= 0 && ix < bias) && iy == 0) ||
                    ((iy >= 0 && iy <= height) && ix == bias) || (ix >= bias && ix <= (bias + height * ratio) && iy == height) ||
                    ((iy >= 0 && iy <= height) && ix == bias + height * ratio)
                    )
                    continue;
                //  исключение периметра 
                if ((ix >= bias && ix <= bias + height * ratio) && (iy <= height))
                    continue;
                //область чипа
                if ((ix == (bias - 1) && iy == height) || (ix == (bias) && iy == height + 1) ||
                    (ix == (bias + height * ratio) && iy == height + 1) || (ix == (bias + height * ratio + 1) && iy == height))
                    continue;
                //область 4-ех точек исключений
                // Пропускаем итерацию для области чипа  и точек рядом с краями

                for (iksx = 0; iksx < nksx; iksx++) {
                    for (iksy = 0; iksy < nksy; iksy++) {
                        for (iksz = 0; iksz < nksz; iksz++) {

                            if ((iksx - n / 2) >= 0 && (iksy - n / 2) >= 0)
                                ff_12[ix][iy][iksx][iksy][iksz] = (ff[ix][iy][iksx][iksy][iksz] -
                                    spd_par * dt * (2 * (iksx - n / 2) + 1) * (ff[ix][iy][iksx][iksy][iksz] - ff[ix - 1][iy][iksx][iksy][iksz]) / dx -
                                    spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[ix][iy][iksx][iksy][iksz] - ff[ix][iy - 1][iksx][iksy][iksz]) / dy);
                            else if ((iksx - n / 2) >= 0 && (iksy - n / 2) < 0)
                                ff_12[ix][iy][iksx][iksy][iksz] = (ff[ix][iy][iksx][iksy][iksz] -
                                    spd_par * dt * (2 * (iksx - n / 2) + 1) * (ff[ix][iy][iksx][iksy][iksz] - ff[ix - 1][iy][iksx][iksy][iksz]) / dx -
                                    spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[ix][iy + 1][iksx][iksy][iksz] - ff[ix][iy][iksx][iksy][iksz]) / dy);

                            else if ((iksx - n / 2) < 0 && (iksy - n / 2) >= 0)
                                ff_12[ix][iy][iksx][iksy][iksz] = (ff[ix][iy][iksx][iksy][iksz] -
                                    spd_par * dt * (2 * (iksx - n / 2) + 1) * (ff[ix + 1][iy][iksx][iksy][iksz] - ff[ix][iy][iksx][iksy][iksz]) / dx -
                                    spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[ix][iy][iksx][iksy][iksz] - ff[ix][iy - 1][iksx][iksy][iksz]) / dy);
                            else if ((iksx - n / 2) < 0 && (iksy - n / 2) < 0)
                                ff_12[ix][iy][iksx][iksy][iksz] = (ff[ix][iy][iksx][iksy][iksz] -
                                    spd_par * dt * (2 * (iksx - n / 2) + 1) * (ff[ix + 1][iy][iksx][iksy][iksz] - ff[ix][iy][iksx][iksy][iksz]) / dx -
                                    spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[ix][iy + 1][iksx][iksy][iksz] - ff[ix][iy][iksx][iksy][iksz]) / dy);

                        }
                    }
                }//циклы по кси закрываются
            }
        }// закрытие циклов по координатам


        //=======================26.05==================точки прилегающие к верхним углам чипа============================


       // bias - 1, height
        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    if ((iksx - n / 2) >= 0 && (iksy - n / 2) >= 0)
                        ff_12[bias - 1][height][iksx][iksy][iksz] = (ff[bias - 1][height][iksx][iksy][iksz] -
                            spd_par * dt * (2 * (iksx - n / 2) + 1) * (ff[bias - 1][height][iksx][iksy][iksz] - ff[bias - 1 - 1][height][iksx][iksy][iksz]) / dx -
                            spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[bias - 1][height][iksx][iksy][iksz] - ff[bias - 1][height - 1][iksx][iksy][iksz]) / dy);



                    else if ((iksx - n / 2) >= 0 && (iksy - n / 2) < 0)
                        ff_12[bias - 1][height][iksx][iksy][iksz] = (ff[bias - 1][height][iksx][iksy][iksz] -
                            spd_par * dt * (2 * (iksx - n / 2) + 1) * (ff[bias - 1][height][iksx][iksy][iksz] - ff[bias - 1 - 1][height][iksx][iksy][iksz]) / dx -
                            spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[bias - 1][height + 1][iksx][iksy][iksz] - ff[bias - 1][height][iksx][iksy][iksz]) / dy);

                    else if ((iksx - n / 2) < 0 && (iksy - n / 2) >= 0)
                        ff_12[bias - 1][height][iksx][iksy][iksz] = (ff[bias - 1][height][iksx][iksy][iksz] -
                            spd_par * dt * (2 * (iksx - n / 2) + 1) * (fexcp1[iksx][iksy][iksz] - ff[bias - 1][height][iksx][iksy][iksz]) / dx -
                            spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[bias - 1][height][iksx][iksy][iksz] - ff[bias - 1][height - 1][iksx][iksy][iksz]) / dy);

                    else if ((iksx - n / 2) < 0 && (iksy - n / 2) < 0)
                        ff_12[bias - 1][height][iksx][iksy][iksz] = (ff[bias - 1][height][iksx][iksy][iksz] -
                            spd_par * dt * (2 * (iksx - n / 2) + 1) * (fexcp1[iksx][iksy][iksz] - ff[bias - 1][height][iksx][iksy][iksz]) / dx -
                            spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[bias - 1][height + 1][iksx][iksy][iksz] - ff[bias - 1][height][iksx][iksy][iksz]) / dy);

                }
            }
        }
        //bias,   height+1
        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    if ((iksx - n / 2) >= 0 && (iksy - n / 2) >= 0)
                        ff_12[bias][height + 1][iksx][iksy][iksz] = (ff[bias][height + 1][iksx][iksy][iksz] -
                            spd_par * dt * (2 * (iksx - n / 2) + 1) * (ff[bias][height + 1][iksx][iksy][iksz] - ff[bias - 1][height + 1][iksx][iksy][iksz]) / dx -
                            spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[bias][height + 1][iksx][iksy][iksz] - fexcp2[iksx][iksy][iksz]) / dy);

                    else if ((iksx - n / 2) >= 0 && (iksy - n / 2) < 0)
                        ff_12[bias][height + 1][iksx][iksy][iksz] = (ff[bias][height + 1][iksx][iksy][iksz] -
                            spd_par * dt * (2 * (iksx - n / 2) + 1) * (ff[bias][height + 1][iksx][iksy][iksz] - ff[bias - 1][height + 1][iksx][iksy][iksz]) / dx -
                            spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[bias][height + 1 + 1][iksx][iksy][iksz] - ff[bias][height + 1][iksx][iksy][iksz]) / dy);

                    else if ((iksx - n / 2) < 0 && (iksy - n / 2) >= 0)
                        ff_12[bias][height + 1][iksx][iksy][iksz] = (ff[bias][height + 1][iksx][iksy][iksz] -
                            spd_par * dt * (2 * (iksx - n / 2) + 1) * (ff[bias + 1][height + 1][iksx][iksy][iksz] - ff[bias][height + 1][iksx][iksy][iksz]) / dx -
                            spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[bias][height + 1][iksx][iksy][iksz] - fexcp2[iksx][iksy][iksz]) / dy);
                    else if ((iksx - n / 2) < 0 && (iksy - n / 2) < 0)
                        ff_12[bias][height + 1][iksx][iksy][iksz] = (ff[bias][height + 1][iksx][iksy][iksz] -
                            spd_par * dt * (2 * (iksx - n / 2) + 1) * (ff[bias + 1][height + 1][iksx][iksy][iksz] - ff[bias][height + 1][iksx][iksy][iksz]) / dx -
                            spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[bias][height + 1 + 1][iksx][iksy][iksz] - ff[bias][height + 1][iksx][iksy][iksz]) / dy);
                }
            }
        }
        //bias+height*ratio,   height+1
        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    if ((iksx - n / 2) >= 0 && (iksy - n / 2) >= 0)
                        ff_12[bias + height * ratio][height + 1][iksx][iksy][iksz] = (ff[bias + height * ratio][height + 1][iksx][iksy][iksz] -
                            spd_par * dt * (2 * (iksx - n / 2) + 1) * (ff[bias + height * ratio][height + 1][iksx][iksy][iksz] - ff[bias + height * ratio - 1][height + 1][iksx][iksy][iksz]) / dx -
                            spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[bias + height * ratio][height + 1][iksx][iksy][iksz] - fexcp3[iksx][iksy][iksz]) / dy);
                    else if ((iksx - n / 2) >= 0 && (iksy - n / 2) < 0)
                        ff_12[bias + height * ratio][height + 1][iksx][iksy][iksz] = (ff[bias + height * ratio][height + 1][iksx][iksy][iksz] -
                            spd_par * dt * (2 * (iksx - n / 2) + 1) * (ff[bias + height * ratio][height + 1][iksx][iksy][iksz] - ff[bias + height * ratio - 1][height + 1][iksx][iksy][iksz]) / dx -
                            spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[bias + height * ratio][height + 1 + 1][iksx][iksy][iksz] - ff[bias + height * ratio][height + 1][iksx][iksy][iksz]) / dy);

                    else if ((iksx - n / 2) < 0 && (iksy - n / 2) >= 0)
                        ff_12[bias + height * ratio][height + 1][iksx][iksy][iksz] = (ff[bias + height * ratio][height + 1][iksx][iksy][iksz] -
                            spd_par * dt * (2 * (iksx - n / 2) + 1) * (ff[bias + height * ratio + 1][height + 1][iksx][iksy][iksz] - ff[bias + height * ratio][height + 1][iksx][iksy][iksz]) / dx -
                            spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[bias + height * ratio][height + 1][iksx][iksy][iksz] - fexcp3[iksx][iksy][iksz]) / dy);
                    else if ((iksx - n / 2) < 0 && (iksy - n / 2) < 0)
                        ff_12[bias + height * ratio][height + 1][iksx][iksy][iksz] = (ff[bias + height * ratio][height + 1][iksx][iksy][iksz] -
                            spd_par * dt * (2 * (iksx - n / 2) + 1) * (ff[bias + height * ratio + 1][height + 1][iksx][iksy][iksz] - ff[bias + height * ratio][height + 1][iksx][iksy][iksz]) / dx -
                            spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[bias + height * ratio][height + 1 + 1][iksx][iksy][iksz] - ff[bias + height * ratio][height + 1][iksx][iksy][iksz]) / dy);
                }
            }
        }

        //bias+height*ratio+1,   height
        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    if ((iksx - n / 2) >= 0 && (iksy - n / 2) >= 0)
                        ff_12[bias + height * ratio + 1][height][iksx][iksy][iksz] = (ff[bias + height * ratio + 1][height][iksx][iksy][iksz] -
                            spd_par * dt * (2 * (iksx - n / 2) + 1) * (ff[bias + height * ratio + 1][height][iksx][iksy][iksz] - fexcp4[iksx][iksy][iksz]) / dx -
                            spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[bias + height * ratio + 1][height][iksx][iksy][iksz] - ff[bias + height * ratio + 1][height - 1][iksx][iksy][iksz]) / dy);
                    else if ((iksx - n / 2) >= 0 && (iksy - n / 2) < 0)
                        ff_12[bias + height * ratio + 1][height][iksx][iksy][iksz] = (ff[bias + height * ratio + 1][height][iksx][iksy][iksz] -
                            spd_par * dt * (2 * (iksx - n / 2) + 1) * (ff[bias + height * ratio + 1][height][iksx][iksy][iksz] - fexcp4[iksx][iksy][iksz]) / dx -
                            spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[bias + height * ratio + 1][height + 1][iksx][iksy][iksz] - ff[bias + height * ratio + 1][height][iksx][iksy][iksz]) / dy);
                    else if ((iksx - n / 2) < 0 && (iksy - n / 2) >= 0)
                        ff_12[bias + height * ratio + 1][height][iksx][iksy][iksz] = (ff[bias + height * ratio + 1][height][iksx][iksy][iksz] -
                            spd_par * dt * (2 * (iksx - n / 2) + 1) * (ff[bias + height * ratio + 1 + 1][height][iksx][iksy][iksz] - ff[bias + height * ratio + 1][height][iksx][iksy][iksz]) / dx -
                            spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[bias + height * ratio + 1][height][iksx][iksy][iksz] - ff[bias + height * ratio + 1][height - 1][iksx][iksy][iksz]) / dy);
                    else if ((iksx - n / 2) < 0 && (iksy - n / 2) < 0)
                        ff_12[bias + height * ratio + 1][height][iksx][iksy][iksz] = (ff[bias + height * ratio + 1][height][iksx][iksy][iksz] -
                            spd_par * dt * (2 * (iksx - n / 2) + 1) * (ff[bias + height * ratio + 1 + 1][height][iksx][iksy][iksz] - ff[bias + height * ratio + 1][height][iksx][iksy][iksz]) / dx -
                            spd_par * dt * (2 * (iksy - n / 2) + 1) * (ff[bias + height * ratio + 1][height + 1][iksx][iksy][iksz] - ff[bias + height * ratio + 1][height][iksx][iksy][iksz]) / dy);


                }
            }
        }

        //=======================26.05==================точки прилегающие к верхним углам чипа============================



















       // запуск интеграла столкновений по координатам
        for (ix = 1; ix < nx; ix++) {
            for (iy = 1; iy < ny; iy++) {
                if (
                    ((iy >= 0 && iy <= ny) && ix == 0) || ((ix >= 0 && ix <= nx) && iy == ny) || ((iy >= 0 && iy <= ny) && ix == nx) ||
                    ((ix >= (bias + height * ratio) && ix <= nx) && iy == 0) || ((ix >= 0 && ix < bias) && iy == 0) ||
                    ((iy >= 0 && iy <= height) && ix == bias) || (ix >= bias && ix <= (bias + height * ratio) && iy == height) ||
                    ((iy >= 0 && iy <= height) && ix == bias + height * ratio)
                    )
                    continue;
                //  исключение периметра 
                if ((ix >= bias && ix <= bias + height * ratio) && (iy <= height))
                    continue;



                // начало шага по интегралу столкновений
                // 
                nsl = 0;
                intkon = 0;
                rgm = 0.0;
                ////////////////////////////////////////////////начало вставки////////////////////////////////////
                dks = 1.0;

                koef = 1;
                //koef1 = 15;



                rgm = 0.0;
                //ns = 120000; //количество точек выбирается в ручную
                nt = 60000; // количетсво циклов по времени
                //srand(time(0));
                srand(time(0));


                //         for (kt = 1; kt <= nt; kt++) {
                nsl = 0;
                intkon = 0;
                rgm = 0.0;

                // Установить генератор случайных чисел
             //  srand(time(0));
                nsl = 0;

                // Обнуление ints
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        for (k = 0; k < n; k++) {
                            ints[i][j][k] = 0;
                            //cout<<ff[i][j][k]<<" ";
                        }
                    }
                }



                for (ks = 1; ks <= ns; ks++) {


                    // srand(time(0));


                    dks = 1.0;
                    bp = 0.123 * (rand() % 11);
                    bp = bp / (10);
                    bp = bp / (0.123);



                    ep = 2 * 3.14159265358979323846 * (rand() % 41) / 40; //случаный азимутальный угол


                    i = -((n) / 2) + rand() % ((n + 1));

                    npr = ceil(sqrt((n / 2) * (n / 2) - i * i));

                    j = -((npr) / 2) + rand() % ((npr + 1));

                    npr = ceil(sqrt((n / 2) * (n / 2) - i * i - j * j));
                    k = -((npr) / 2) + rand() % ((npr + 1));



                    ial = (2 * i + 1) * 0.5; jal = (2 * j + 1) * 0.5; kal = (2 * k + 1)*0.5; //необходимо для нечетной решетки


                    ksx = (ial)*dks;
                    ksy = (jal)*dks;
                    ksz = (kal)*dks;

                    i = -(n / 2) + rand() % (n + 1);
                    j = -(n / 2) + rand() % (n + 1);
                    k = -(n / 2) + rand() % (n + 1);

                    //  i=-((n-6) / 2)+rand()%((n-6)+1);
                     // j=-((n-6) / 2)+rand()%((n-6)+1);
                     // k=-((n-6) / 2)+rand()%((n-6)+1);



                    ibt =( 2 * i + 1) * 0.5; (jbt = 2 * j + 1) * 0.5; kbt = (2 * k + 1) * 0.5;



                    ks1x = (ibt)*dks;
                    ks1y = (jbt)*dks;
                    ks1z = (kbt)*dks;


                    if (((ial * ial + jal * jal + kal * kal) < 100) && ((ibt * ibt + jbt * jbt + kbt * kbt) < 100)) {       // {obrezanie po skorostyam}
                        nsl++;

                        koef = 1;
                        // koef1 = 15;


                        ksx = (ial)*dks;
                        ksy = (jal)*dks;
                        ksz = (kal)*dks;

                        ks1x = (ibt)*dks;
                        ks1y = (jbt)*dks;
                        ks1z = (kbt)*dks;


                        //cout<<ial<<" "<<jal<<" "<<kal<<" "<<ibt<<" "<<jbt<<" "<<kbt<<endl;
                        //cout<<ksx<<" "<<ks1x<<" "<<ksy<<" "<<ks1y<<" "<<ksz<<" "<<ks1z<<endl;


                        //cin>>iii;




                        teta = acos(bp);

                        gx = ks1x - ksx;
                        gy = ks1y - ksy;
                        gz = ks1z - ksz;
                        gxy = sqrt(gx * gx + gy * gy);
                        g = sqrt(gx * gx + gy * gy + gz * gz);


                        if (gxy != 0) {
                            gshx = gx * cos(teta) - (gx * gz / gxy) * cos(ep) * sin(teta) + (g * gy / gxy) * sin(ep) * sin(teta);
                            gshy = gy * cos(teta) - (gy * gz / gxy) * cos(ep) * sin(teta) - (g * gx / gxy) * sin(ep) * sin(teta);
                            gshz = gz * cos(teta) + gxy * cos(ep) * sin(teta);
                        }
                        else {
                            gshx = g * sin(ep) * sin(teta);
                            gshy = g * cos(ep) * sin(teta);
                            gshz = g * cos(teta);
                        }

                        kshx = 0.5 * (ksx + ks1x) - gshx / 2;
                        kshy = 0.5 * (ksy + ks1y) - gshy / 2;
                        kshz = 0.5 * (ksz + ks1z) - gshz / 2;

                        ksh1x = 0.5 * (ksx + ks1x) + gshx / 2;
                        ksh1y = 0.5 * (ksy + ks1y) + gshy / 2;
                        ksh1z = 0.5 * (ksz + ks1z) + gshz / 2;



                        gg = sqrt(gx * gx + gy * gy + gz * gz);


                        s1 = ksx * ksx + ksy * ksy + ksz * ksz + ks1x * ks1x + ks1y * ks1y + ks1z * ks1z; //для точной скорости
                        s2 = kshx * kshx + kshy * kshy + kshz * kshz + ksh1x * ksh1x + ksh1y * ksh1y + ksh1z * ksh1z; //для точных  разлетных скоростей

                        en0 = s2;


                        si = trunc((kshx) / 2);



                        if (kshx > 0)
                            ibl = 1 + si * 2;
                        else
                            ibl = -1 + si * 2;




                        si = trunc((kshy) / 2);

                        if (kshy > 0)
                            jbl = 1 + si * 2;
                        else
                            jbl = -1 + si * 2;



                        si = trunc((kshz) / 2);

                        if (kshz > 0)
                            kbl = 1 + si * 2;
                        else
                            kbl = -1 + si * 2;


                        si = trunc((ksh1x) / 2);


                        if (ksh1x > 0)
                            i1bl = 1 + si * 2;
                        else
                            i1bl = -1 + si * 2;

                        si = trunc((ksh1y) / 2);

                        if (ksh1y > 0)
                            j1bl = 1 + si * 2;
                        else
                            j1bl = -1 + si * 2;

                        si = trunc((ksh1z) / 2);

                        if (ksh1z > 0)
                            k1bl = 1 + si * 2;
                        else
                            k1bl = -1 + si * 2;
                        //нашли ближайшие точки

                        //cout<<"blizhaish k razletnym "<<ibl<<" "<<jbl<<" "<<kbl<<" "<<i1bl<<" "<<j1bl<<" "<<k1bl<<endl;
                        //cout<<ksx<<" "<<ks1x<<" "<<ksy<<" "<<ks1y<<" "<<ksz<<" "<<ks1z<<endl;

                        //cout<<nsl<<endl;
                        //cin>>iii;


                        fLr1 = 0;
                        //en0=1.2;
                        //enp=1.2;
                        en0 = (kshx * kshx + kshy * kshy + kshz * kshz + ksh1x * ksh1x + ksh1y * ksh1y + ksh1z * ksh1z);
                        enp = (ibl * ibl + jbl * jbl + kbl * kbl + i1bl * i1bl + j1bl * j1bl + k1bl * k1bl) * dks * dks;

                        en0cl = floor(s1);
                        enpcl = floor(enp);

                        if (en0cl == enpcl) {

                            fLr1 = 1;



                            //cout<<" en0=enp"<<endl;
                            //cin>>iii;

                            lmbi = ibl;
                            lmbj = jbl;
                            lmbk = kbl;

                            mui = i1bl;
                            muj = j1bl;
                            muk = k1bl;

                            lmpsi = lmbi;
                            lmpsj = lmbj;
                            lmpsk = lmbk;

                            mumsi = mui;
                            mumsj = muj;
                            mumsk = muk;

                            rgm = 1; //r ню ( не гаммма)
                            //cout<<" en0 "<<en0<<endl;
                            //cout<<" enp "<<enp<<endl;
                            //cout<<rgm<<endl;
                        }
                        //три случая страница 15, сначала 3, потом 2, первый случай выше
                        else {

                            if (enp < en0) {
                                lmbi = ibl;
                                lmbj = jbl;
                                lmbk = kbl;

                                mui = i1bl;
                                muj = j1bl;
                                muk = k1bl;

                                xcm = 0.5 * (lmbi + mui);
                                ycm = 0.5 * (lmbj + muj);
                                zcm = 0.5 * (lmbk + muk);

                                if (lmbi < xcm)
                                    di = -2;
                                else
                                    di = 2;
                                if (lmbi == xcm)
                                    di = 0;

                                if (lmbj < ycm)
                                    dj = -2;
                                else
                                    dj = 2;
                                if (lmbj == ycm)
                                    dj = 0;

                                if (lmbk < zcm)
                                    dk = -2;
                                else
                                    dk = 2;
                                if (lmbk == zcm)
                                    dk = 0;

                                rpr = (lmbi + di - xcm) * (lmbi + di - xcm) + (lmbj + dj - ycm) * (lmbj + dj - ycm) + (lmbk + dk - zcm) * (lmbk + dk - zcm); // расстояние до точки центра масс от точки блашей вне круга (скоростная сетка)
                                enn[1] = 2 * (xcm * xcm + ycm * ycm + zcm * zcm) + 2 * rpr;
                                rr[1] = (lmbi + di - kshx) * (lmbi + di - kshx) + (lmbj + dj - kshy) * (lmbj + dj - kshy) + (lmbk + dk - kshz) * (lmbk + dk - kshz);

                                lmpsi = lmbi + di;
                                lmpsj = lmbj + dj;
                                lmpsk = lmbk + dk;

                                mnr = rr[1]; //начальное значение минимума расстояния (пункт е первый метод)
                                //дпльше перебор  оставшихся 6 точек
                                rpr = (lmbi - xcm) * (lmbi - xcm) + (lmbj + dj - ycm) * (lmbj + dj - ycm) + (lmbk + dk - zcm) * (lmbk + dk - zcm);
                                enn[2] = 2 * rpr + 2 * (xcm * xcm + ycm * ycm + zcm * zcm);
                                s = enn[2];
                                rr[2] = (lmbi - kshx) * (lmbi - kshx) + (lmbj + dj - kshy) * (lmbj + dj - kshy) + (lmbk + dk - kshz) * (lmbk + dk - kshz);

                                if (s > en0 && rr[2] < mnr) {
                                    mnr = rr[2];
                                    lmpsi = lmbi;
                                    lmpsj = lmbj + dj;
                                    lmpsk = lmbk + dk;
                                }
                                rpr = (lmbi + di - xcm) * (lmbi + di - xcm) + (lmbj - ycm) * (lmbj - ycm) + (lmbk + dk - zcm) * (lmbk + dk - zcm);
                                enn[3] = 2 * rpr + 2 * (xcm * xcm + ycm * ycm + zcm * zcm);
                                s = enn[3];
                                rr[3] = (lmbi + di - kshx) * (lmbi + di - kshx) + (lmbj - kshy) * (lmbj - kshy) + (lmbk + dk - kshz) * (lmbk + dk - kshz);

                                if ((s > en0) && (rr[3] < mnr)) {
                                    mnr /*minivum rasst*/ = rr[3];
                                    lmpsi /*lambda + s*/ = lmbi + di;
                                    lmpsj /*lambda + s*/ = lmbj;
                                    lmpsk /*lambda + s*/ = lmbk + dk;
                                }

                                rpr = (lmbi + di - xcm) * (lmbi + di - xcm) + (lmbj + dj - ycm) * (lmbj + dj - ycm) + (lmbk - zcm) * (lmbk - zcm);
                                enn[4] = 2 * rpr + 2 * (xcm * xcm + ycm * ycm + zcm * zcm);
                                rr[4] = (lmbi + di - kshx) * (lmbi + di - kshx) + (lmbj + dj - kshy) * (lmbj + dj - kshy) + (lmbk - kshz) * (lmbk - kshz);
                                if ((enn[4] > en0) && (rr[4] < mnr)) {
                                    mnr /*minivum rasst*/ = rr[4];
                                    lmpsi /*lambda + s*/ = lmbi + di;
                                    lmpsj /*lambda + s*/ = lmbj + dj;
                                    lmpsk /*lambda + s*/ = lmbk;
                                }

                                rpr = (lmbi - xcm) * (lmbi - xcm) + (lmbj - ycm) * (lmbj - ycm) + (lmbk + dk - zcm) * (lmbk + dk - zcm);
                                enn[5] = 2 * rpr + 2 * (xcm * xcm + ycm * ycm + zcm * zcm);
                                rr[5] = (lmbi - kshx) * (lmbi - kshx) + (lmbj - kshy) * (lmbj - kshy) + (lmbk + dk - kshz) * (lmbk + dk - kshz);
                                if ((enn[5] > en0) && (rr[5] < mnr)) {
                                    mnr /*minivum rasst*/ = rr[5];
                                    lmpsi /*lambda + s*/ = lmbi;
                                    lmpsj /*lambda + s*/ = lmbj;
                                    lmpsk /*lambda + s*/ = lmbk + dk;
                                }

                                rpr = (lmbi - xcm) * (lmbi - xcm) + (lmbj + dj - ycm) * (lmbj + dj - ycm) + (lmbk - zcm) * (lmbk - zcm);
                                enn[6] = 2 * rpr + 2 * (xcm * xcm + ycm * ycm + zcm * zcm);
                                rr[6] = (lmbi - kshx) * (lmbi - kshx) + (lmbj + dj - kshy) * (lmbj + dj - kshy) + (lmbk - kshz) * (lmbk - kshz);
                                if ((enn[6] > en0) && (rr[6] < mnr)) {
                                    mnr /*minivum rasst*/ = rr[6];
                                    lmpsi /*lambda + s*/ = lmbi;
                                    lmpsj /*lambda + s*/ = lmbj + dj;
                                    lmpsk /*lambda + s*/ = lmbk;
                                }

                                rpr = (lmbi + di - xcm) * (lmbi + di - xcm) + (lmbj - ycm) * (lmbj - ycm) + (lmbk - zcm) * (lmbk - zcm);
                                enn[7] = 2 * rpr + 2 * (xcm * xcm + ycm * ycm + zcm * zcm);
                                rr[7] = (lmbi + di - kshx) * (lmbi + di - kshx) + (lmbj - kshy) * (lmbj - kshy) + (lmbk - kshz) * (lmbk - kshz);
                                if ((enn[7] > en0) && (rr[7] < mnr)) {
                                    mnr /*minivum rasst*/ = rr[7];
                                    lmpsi /*lambda + s*/ = lmbi + di;
                                    lmpsj /*lambda + s*/ = lmbj;
                                    lmpsk /*lambda + s*/ = lmbk;
                                }

                                /*posle nahodim parnyi usel punkt zh) str. 16 po 2.18*/

                                /*mumsi-> mu minus s i-toe*/
                                mumsi = lmbi + mui - lmpsi;
                                mumsj = lmbj + muj - lmpsj;
                                mumsk = lmbk + muk - lmpsk;

                                /*str 14 formula 2.20*/
                                en2g = lmpsi * lmpsi + lmpsj * lmpsj + lmpsk * lmpsk + mumsi * mumsi + mumsj * mumsj + mumsk * mumsk;
                                //if (en2g = enp)

                                //    break;
                                if (floor(en2g * dks) != floor(enp * dks)) {
                                    rgm = (en0 - enp) / (en2g - enp);
                                }
                                else
                                    rgm = 1;



                                /*enp=en1 if enp<en0 then*/
                               //cout<<" en2g "<<en2g<<endl;
                               //cout<<" en0 "<<en0<<endl;
                               //cout<<" enp "<<enp<<endl;

                               //cout<<rgm<<endl;
                            }
                            //второй случай пункта г
                            if (enp > en0) {
                                lmpsi = ibl;
                                lmpsj = jbl;
                                lmpsk = kbl;

                                mumsi = i1bl;
                                mumsj = j1bl;
                                mumsk = k1bl;

                                xcm = 0.5 * (lmpsi + mumsi);
                                ycm = 0.5 * (lmpsj + mumsj);
                                zcm = 0.5 * (lmpsk + mumsk);

                                if (lmpsi < xcm)
                                    di = 2;
                                else
                                    di = -2;
                                if (lmpsi == xcm)
                                    di = 0;

                                if (lmpsj < ycm)
                                    dj = 2;
                                else
                                    dj = -2;
                                if (lmpsj == ycm)
                                    dj = 0;

                                if (lmpsk < zcm)
                                    dk = 2;
                                else
                                    dk = -2;
                                if (lmpsk == zcm)
                                    dk = 0;

                                rpr = (lmpsi + di - xcm) * (lmpsi + di - xcm) + (lmpsj + dj - ycm) * (lmpsj + dj - ycm) + (lmpsk + dk - zcm) * (lmpsk + dk - zcm);
                                enn[1] = 2 * (xcm * xcm + ycm * ycm + zcm * zcm) + 2 * rpr;
                                rr[1] = (lmpsi + di - kshx) * (lmpsi + di - kshx) + (lmpsj + dj - kshy) * (lmpsj + dj - kshy) + (lmpsk + dk - kshz) * (lmpsk + dk - kshz);
                                lmbi = lmpsi + di;
                                lmbj = lmpsj + dj;
                                lmbk = lmpsk + dk;

                                mnr = rr[1];

                                rpr = (lmpsi - xcm) * (lmpsi - xcm) + (lmpsj + dj - ycm) * (lmpsj + dj - ycm) + (lmpsk + dk - zcm) * (lmpsk + dk - zcm);
                                enn[2] = 2 * rpr + 2 * (xcm * xcm + ycm * ycm + zcm * zcm);
                                s = enn[2];
                                rr[2] = (lmpsi - kshx) * (lmpsi - kshx) + (lmpsj + dj - kshy) * (lmpsj + dj - kshy) + (lmpsk + dk - kshz) * (lmpsk + dk - kshz);
                                if (s < en0 && rr[2] < mnr) {
                                    mnr = rr[2];
                                    lmbi = lmpsi;
                                    lmbj = lmpsj + dj;
                                    lmbk = lmpsk + dk;
                                }

                                rpr = (lmpsi + di - xcm) * (lmpsi + di - xcm) + (lmpsj - ycm) * (lmpsj - ycm) + (lmpsk + dk - zcm) * (lmpsk + dk - zcm);
                                enn[3] = 2 * rpr + 2 * (xcm * xcm + ycm * ycm + zcm * zcm);
                                s = enn[3];
                                rr[3] = (lmpsi + di - kshx) * (lmpsi + di - kshx) + (lmpsj - kshy) * (lmpsj - kshy) + (lmpsk + dk - kshz) * (lmpsk + dk - kshz);
                                if ((s < en0) && (rr[3] < mnr)) {
                                    mnr /*minivum rasst*/ = rr[3];
                                    lmbi /*lambda + s*/ = lmpsi + di;
                                    lmbj /*lambda + s*/ = lmpsj;
                                    lmbk /*lambda + s*/ = lmpsk + dk;
                                }

                                rpr = (lmpsi + di - xcm) * (lmpsi + di - xcm) + (lmpsj + dj - ycm) * (lmpsj + dj - ycm) + (lmpsk - zcm) * (lmpsk - zcm);
                                enn[4] = 2 * rpr + 2 * (xcm * xcm + ycm * ycm + zcm * zcm);
                                rr[4] = (lmpsi + di - kshx) * (lmpsi + di - kshx) + (lmpsj + dj - kshy) * (lmpsj + dj - kshy) + (lmpsk - kshz) * (lmpsk - kshz);

                                if ((enn[4] < en0) && (rr[4] < mnr)) {
                                    mnr /*minivum rasst*/ = rr[4];
                                    lmbi /*lambda + s*/ = lmpsi + di;
                                    lmbj /*lambda + s*/ = lmpsj + dj;
                                    lmbk /*lambda + s*/ = lmpsk;
                                }

                                rpr = (lmpsi - xcm) * (lmpsi - xcm) + (lmpsj - ycm) * (lmpsj - ycm) + (lmpsk + dk - zcm) * (lmpsk + dk - zcm);
                                enn[5] = 2 * rpr + 2 * (xcm * xcm + ycm * ycm + zcm * zcm);
                                rr[5] = (lmpsi - kshx) * (lmpsi - kshx) + (lmpsj - kshy) * (lmpsj - kshy) + (lmpsk + dk - kshz) * (lmpsk + dk - kshz);
                                if ((enn[5] < en0) && (rr[5] < mnr)) {
                                    mnr /*minivum rasst*/ = rr[5];
                                    lmbi /*lambda + s*/ = lmpsi;
                                    lmbj /*lambda + s*/ = lmpsj;
                                    lmbk /*lambda + s*/ = lmpsk + dk;
                                }

                                rpr = (lmpsi - xcm) * (lmpsi - xcm) + (lmpsj + dj - ycm) * (lmpsj + dj - ycm) + (lmpsk - zcm) * (lmpsk - zcm);
                                enn[6] = 2 * rpr + 2 * (xcm * xcm + ycm * ycm + zcm * zcm);
                                rr[6] = (lmpsi - kshx) * (lmpsi - kshx) + (lmpsj + dj - kshy) * (lmpsj + dj - kshy) + (lmpsk - kshz) * (lmpsk - kshz);
                                if ((enn[6] < en0) && (rr[6] < mnr)) {
                                    mnr /*minivum rasst*/ = rr[6];
                                    lmbi /*lambda + s*/ = lmpsi;
                                    lmbj /*lambda + s*/ = lmpsj + dj;
                                    lmbk /*lambda + s*/ = lmpsk;
                                }

                                rpr = (lmpsi + di - xcm) * (lmpsi + di - xcm) + (lmpsj - ycm) * (lmpsj - ycm) + (lmpsk - zcm) * (lmpsk - zcm);
                                enn[7] = 2 * rpr + 2 * (xcm * xcm + ycm * ycm + zcm * zcm);
                                rr[7] = (lmpsi + di - kshx) * (lmpsi + di - kshx) + (lmpsj - kshy) * (lmpsj - kshy) + (lmpsk - kshz) * (lmpsk - kshz);
                                if ((enn[7] < en0) && (rr[7] < mnr)) {
                                    mnr /*minivum rasst*/ = rr[7];
                                    lmbi /*lambda + s*/ = lmpsi + di;
                                    lmbj /*lambda + s*/ = lmpsj;
                                    lmbk /*lambda + s*/ = lmpsk;
                                }


                                /*posle nahodim parnyi usel punkt zh) str. 16 po 2.18*/
                                //относительно цетра масс по закону сохранения импульса
                                /*mumsi-> mu minus s i-toe*/
                                mui = lmpsi + mumsi - lmbi;
                                muj = lmpsj + mumsj - lmbj;
                                muk = lmpsk + mumsk - lmbk;

                                /*str 14 formula 2.20*/
                                en1g = lmbi * lmbi + lmbj * lmbj + lmbk * lmbk + mui * mui + muj * muj + muk * muk;
                                //if (en2g = enp)

                                //    break;

                                //obranotka redkogo sluchaya sonapravlennyh skorostey b azimut ugla =0;  2pi
                                if (floor(en1g * dks) != floor(enp * dks)) {
                                    rgm = (en0 - en1g) / (enp - en1g);
                                }
                                else
                                    rgm = 1;



                                /*enp=en1 if enp<en0 then*/


                                //cout<<" en1g "<<en1g<<endl;
                                //cout<<" en0 "<<en0<<endl;
                                //cout<<" enp "<<enp<<endl;

                                //cout<<rgm<<endl;

                            }// конец второго пунка г

                        /* =====================4.06===================== */
                        // nahodim omega bolshoe gamma

                        // nahozhu stepen ()^(1-r)
                        // promezhutocynaya peremenn s, s1


                          /* obrezanie po skorost */






                        }  // konec 2 & 3

                        //cout<<lmbi<<" "<<lmbj<<" "<<lmbk<<" "<<mui<<"  "<<muj<<" "<<muk<<endl;
                        //cout<<lmpsi<<" "<<lmpsj<<" "<<lmpsk<<" "<<mumsi<<" "<<mumsj<<" "<<mumsk<<endl;

                        //cout<<ksx<<" "<<ks1x<<" "<<ksy<<" "<<ks1y<<" "<<ksz<<" "<<ks1z<<endl;


                        //cin>>iii;




                        if (fLr1 == 0) {

                            s = (ff_12[ix][iy][(lmbi + n) / 2][(lmbj + n) / 2][(lmbk + n) / 2]) * (ff_12[ix][iy][(mui + n) / 2][(muj + n) / 2][(muk + n) / 2]);//первый множитель в первом слагаемом в первой скобке глобальной
                            if (abs(s) > 0) {
                                s = pow(abs(s), (1 - rgm));
                            }
                            else { s = 0; };


                            s1 = (ff_12[ix][iy][(lmpsi + n) / 2][(lmpsj + n) / 2][(lmpsk + n) / 2]) * (ff_12[ix][iy][(mumsi + n) / 2][(mumsj + n) / 2][(mumsk + n) / 2]); //второй множитель в первом слагаемом в первой скобке глобальной
                            if ((abs(s1) > 0)) {
                                s1 = pow(abs(s1), rgm);
                            }
                            else
                            {
                                s1 = 0;
                            };



                            s2 = sqrt((ial - ibt) * (ial - ibt) + (jal - jbt) * (jal - jbt) + (kal - kbt) * (kal - kbt));

                            //cout<<ial<<jal<<kal<<endl;

                            ombgm = (s * s1 - ff_12[ix][iy][(ial + n) / 2][(jal + n) / 2][(kal + n) / 2] * ff_12[ix][iy][(ibt + n) / 2][(jbt + n) / 2][(kbt + n) / 2]) * s2 * bp;

                        }

                        else {


                            //cout<<ial<<jal<<kal<<endl;


                            ombgm = (ff_12[ix][iy][(lmpsi + n) / 2][(lmpsj + n) / 2][(lmpsk + n) / 2] * ff_12[ix][iy][(mumsi + n) / 2][(mumsj + n) / 2][(mumsk + n) / 2]
                                - ff_12[ix][iy][(ial + n) / 2][(jal + n) / 2][(kal + n) / 2] * ff_12[ix][iy][(ibt + n) / 2][(jbt + n) / 2][(kbt + n) / 2]) * s2 * bp;


                            //cout<<ial<<jal<<kal<<endl;


                        }

                        if (ombgm != ombgm)//nan
                        {
                            continue;

                            cout << " omngm " << ombgm << endl;

                            cout << ial << " " << jal << "  " << kal << endl;
                            cout << ibt << "  " << jbt << " " << kbt << endl;
                            cout << lmbi << " " << lmbj << " " << lmbk << endl;
                            cout << mui << " " << muj << " " << muk << endl;
                            cout << " en0 " << en0 << endl;

                            cout << " enp " << enp << endl;
                            cout << " trunc en0 " << floor(en0) << endl;
                            cout << " en0cl " << en0cl << endl;
                            cout << " enpcl " << enpcl << endl;
                            cout << " fl " << fLr1 << endl;

                            cout << " rgm " << rgm << endl;

                            // cin >> iii;

                        }

                        koef = 1;


                        //cout<<"om "<<ombgm<<endl;

                        ints[(ial + n) / 2][(jal + n) / 2][(kal + n) / 2] = ints[(ial + n) / 2][(jal + n) / 2][(kal + n) / 2] + koef * ombgm;

                        ints[(ibt + n) / 2][(jbt + n) / 2][(kbt + n) / 2] = ints[(ibt + n) / 2][(jbt + n) / 2][(kbt + n) / 2] + koef * ombgm;


                        ints[(lmbi + n) / 2][(lmbj + n) / 2][(lmbk + n) / 2] = ints[(lmbi + n) / 2][(lmbj + n) / 2][(lmbk + n) / 2] - koef * (1 - rgm) * ombgm;

                        ints[(mui + n) / 2][(muj + n) / 2][(muk + n) / 2] = ints[(mui + n) / 2][(muj + n) / 2][(muk + n) / 2] - koef * (1 - rgm) * ombgm;


                        ints[(mumsi + n) / 2][(mumsj + n) / 2][(mumsk + n) / 2] = ints[(mumsi + n) / 2][(mumsj + n) / 2][(mumsk + n) / 2] - koef * (rgm)*ombgm;


                        ints[(lmpsi + n) / 2][(lmpsj + n) / 2][(lmpsk + n) / 2] = ints[(lmpsi + n) / 2][(lmpsj + n) / 2][(lmpsk + n) / 2] - koef * (rgm)*ombgm;
                        //cout<<"koef "<<koef<<" omgm "<<ombgm<<endl;

                        //cout<<(ial+n)/2<<"  "<<(jal+n)/2<<"  "<<(kal+n)/2<<" ints "<<ints[(ial+n)/ 2][(jal+n)/ 2][(kal+n)/ 2]<<endl;


                        //cin>>iii;
                        // konec  cikl obrezanie po skorostyam

                    }

                    //cout<<"     "<<ints[20][20][20]<<"  "<<endl;
                    //cout<<"  nsl   "<<nsl<<"  "<<endl;
                    //cin>>iii;

                    //cout<<" nsl"<<nsl<<" koef "<<koef<<endl;
                    //cin>>iii;



                    //for (i=0;i<=n;i++){
                    //cout<<endl;
                    //cout<<endl;
                    //for (j=0;j<=n;j++){
                    //cout<<endl;
                    //for (k=0;k<=n;k++){
                    //cout<<ints[i][j][k]<<" ";
                    //}
                    //}
                    //}




                }/*cikl ks*/


                //cout<<" nsl"<<nsl<<" koef "<<koef<<endl;;
                //cin>>iii;
               // koef1 = 15;//10.0;
                //                koef1 = 0.123 * nsl / ns;
                     //           koef1 = koef1 / 0.123;





                intkon = 0;//счетчик для проверки концентрации

                ennb = 0.0;
                tprod = 0.0;
                tpoper = 0.0;
                px = 0.0;
                py = 0.0;
                pz = 0.0;
                dks = 1.0;
                //cout<<"nsl "<<nsl<<endl;
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        for (k = 0; k < n; k++) {
                            intkon = intkon + ints[i][j][k];

                            ff_12[ix][iy][i][j][k] = ff_12[ix][iy][i][j][k] + dt2 * koef1 * ints[i][j][k];

                            ennb = ennb + ff_12[ix][iy][i][j][k] * ((2 * (i - 10) + 1) * (2 * (i - 10) + 1) + (2 * (j - 10) + 1) * (2 * (j - 10) + 1) + (2 * (k - 10) + 1) * (2 * (k - 10) + 1)) * dks * dks * 0.25;
                            tprod = tprod + ff_12[ix][iy][i][j][k] * ((2 * (i - 10) + 1) * (2 * (i - 10) + 1)) * dks * dks * 0.25;
                            tpoper = tpoper + 0.5 * (ff_12[ix][iy][i][j][k] * ((2 * (j - 10) + 1) * (2 * (j - 10) + 1) + (2 * (k - 10) + 1) * (2 * (k - 10) + 1))) * dks * dks * 0.25;
                            px = px + ff_12[ix][iy][i][j][k] * (2 * (i - 10) + 1) * dks * 0.5;
                            py = py + ff_12[ix][iy][i][j][k] * (2 * (j - 10) + 1) * dks * 0.5;
                            pz = pz + ff_12[ix][iy][i][j][k] * (2 * (k - 10) + 1) * dks * 0.5;

                        }
                    }
                }

            }
        }

        // закрытие цикла по координатам для интеграла столкновений
          // обновление границ должно быть тут
        ////////////////////////////////////////////////конец вставки////////////////////////////////////


        // обновление границ
              // левая стенка трубы начало/////
        for (iy = 0; iy <= ny; iy++) {
            // левая стенка трубы
            ix = 0;

            for (iksx = 0; iksx < nksx / 2; iksx++) {         //< 
                for (iksy = 0; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {
                        if (iy > 0 && iy < ny) {
                            //ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                            //ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;

                            if (iksy >= nksy / 2) {//iksx>nksx/2; ksx>0

                                ff_12[ix][iy][iksx][iksy][iksz] = ff[ix][iy][iksx][iksy][iksz]
                                    - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[ix + 1][iy][iksx][iksy][iksz] - ff[ix][iy][iksx][iksy][iksz])
                                    - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[ix][iy][iksx][iksy][iksz] - ff[ix][iy - 1][iksx][iksy][iksz]);
                            }

                            if (iksy < nksy / 2) {//iksx<nksx/2; ksx<0 нуля нет

                                ff_12[ix][iy][iksx][iksy][iksz] = ff[ix][iy][iksx][iksy][iksz]
                                    - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[ix][iy + 1][iksx][iksy][iksz] - ff[ix][iy][iksx][iksy][iksz])
                                    - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[ix + 1][iy][iksx][iksy][iksz] - ff[ix][iy][iksx][iksy][iksz]);
                            }
                        }
                        else
                            ff_12[ix][iy][iksx][iksy][iksz] = ff[ix][iy][iksx][iksy][iksz];
                    }
                }
            }
            //нашли для прилетающих распреление новое в конкретной точке x,y
  //находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn
            ron = 0;
            vxn = 0;
            vyn = 0;
            for (iksx = 0; iksx < nksx / 2; iksx++) {
                for (iksy = 0; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {
                        ron = ron + ff_12[ix][iy][iksx][iksy][iksz];
                        vxn = vxn + ff_12[ix][iy][iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);//ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                        vyn = vyn + ff_12[ix][iy][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                    }
                }
            }


            //тут надо подумать че будет если ron=0

            if ((ron) == 0) { vx = 0; vy = 0; }
            if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }

            if (abs(vxn) > 0)      ron = 1 * 2 * abs(vxn) * sqrt(ms * M_PI / (2 * T1));


            for (iksx = nksx / 2; iksx < nksx; iksx++) {
                for (iksy = 0; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        ff_12[ix][iy][iksx][iksy][iksz] = 1 * ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 1);
                    }
                }
            }

            vxnp = 0;
            // считаем какой будет поток при koefotr=1
            for (iksx = nksx / 2; iksx < nksx; iksx++) {
                for (iksy = 0; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        vxnp = vxnp + ff_12[ix][iy][iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);

                    }
                }
            }





            if (abs(vxnp) > 0) {
                koefotr = abs(vxn) / abs(vxnp);
            }
            else koefotr = 1;
            // задаем новое отражение для консерватизма


            if (abs(vxn) > 0)      ron = koefotr * 2 * abs(vxn) * sqrt(ms * M_PI / (2 * T1));


            for (iksx = nksx / 2; iksx < nksx; iksx++) {
                for (iksy = 0; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        ff_12[ix][iy][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 1);
                    }
                }
            }




            vxnp = 0;

        }


        ////////////////////////////////////////////конец левой стенки трубы///////////////

        // верхняя стенка трубы  здесь ksy вниз потому  производная по y (ff[ix][ny][iksx][iksy][iksz] - ff[ix][ny-1][iksx][iksy][iksz]) cм 2.3 zvmmf

        for (ix = 0; ix <= nx; ix++) {   //ix=0?1

            for (iksx = 0; iksx < nksx; iksx++) {         //< ? или <=
                for (iksy = nksy / 2; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        //ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                        //ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                        if (ix > 0 && ix < nx) {
                            if (iksx >= nksx / 2) {//iksx>nksx/2; ksx>0

                                ff_12[ix][ny][iksx][iksy][iksz] = ff[ix][ny][iksx][iksy][iksz] - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[ix][ny][iksx][iksy][iksz] - ff[ix][ny - 1][iksx][iksy][iksz])
                                    - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[ix][ny][iksx][iksy][iksz] - ff[ix - 1][ny][iksx][iksy][iksz]);
                            }

                            if (iksx < nksx / 2) {//iksx<nksx/2; ksx<0 нуля нет

                                ff_12[ix][ny][iksx][iksy][iksz] = ff[ix][ny][iksx][iksy][iksz] - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[ix][ny][iksx][iksy][iksz] - ff[ix][ny - 1][iksx][iksy][iksz])
                                    - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[ix + 1][ny][iksx][iksy][iksz] - ff[ix][ny][iksx][iksy][iksz]);
                            }
                        }
                        else
                            ff_12[ix][ny][iksx][iksy][iksz] = ff[ix][ny][iksx][iksy][iksz];
                    }
                }
            }
            //нашли для прилетающих распреление новое в конкретной точке x,y
  //находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn
            ron = 0;
            vxn = 0;
            vyn = 0;
            for (iksx = 0; iksx < nksx; iksx++) {
                for (iksy = nksy / 2; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {
                        ron = ron + ff_12[ix][ny][iksx][iksy][iksz];
                        vxn = vxn + ff_12[ix][ny][iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);//ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                        vyn = vyn + ff_12[ix][ny][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                    }
                }
            }





            //тут надо подумать че будет если ron=0


            if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }
            // здесь  пиши distribution ff_12[ix][ny][]


            //задаем сначала koefotr=1

            if (vyn > 0)      ron = 1 * 2 * vyn * sqrt(ms * M_PI / (2 * T1));



            for (iksx = 0; iksx < nksx; iksx++) {
                for (iksy = 0; iksy < nksy / 2; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {
                        ff_12[ix][ny][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 1);

                    }
                }
            }


            // считаем какой будет поток при koefotr=1
            vynp = 0;

            for (iksx = 0; iksx < nksx; iksx++) {
                for (iksy = 0; iksy < nksy / 2; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        vynp = vynp + ff_12[ix][ny][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);

                    }
                }
            }


            if (abs(vynp) > 0) {
                koefotr = abs(vyn) / abs(vynp);
            }
            else koefotr = 1;
            // задаем новое отражение для консерватизма


            if (vyn > 0)      ron = koefotr * 2 * vyn * sqrt(ms * M_PI / (2 * T1));

            for (iksx = 0; iksx < nksx; iksx++) {
                for (iksy = 0; iksy < nksy / 2; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {
                        ff_12[ix][ny][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 1);

                    }
                }
            }







            vynp = 0;


            // } //end  ix



         //    for (ix = 1; ix < nx; ix++) {
            for (iksx = 0; iksx < nksx; iksx++) {
                for (iksy = 0; iksy < nksy / 2; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {
                        ff_12[ix][ny][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 1);

                    }
                }
            }






            // } //end  ix





        }


        // end верхняя стенка трубы end








     // правая стенка трубы начало/////
        for (iy = 0; iy <= ny; iy++) {
            // правая стенка трубы
            ix = nx;
            for (iksx = nksx / 2; iksx < nksx; iksx++) {         //< 
                for (iksy = 0; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {
                        if (iy > 0 && iy < ny) {
                            if (iksy >= nksy / 2) {//iksx>nksx/2; ksx>0

                                ff_12[ix][iy][iksx][iksy][iksz] = ff[ix][iy][iksx][iksy][iksz]
                                    - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[ix][iy][iksx][iksy][iksz] - ff[ix - 1][iy][iksx][iksy][iksz])
                                    - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[ix][iy][iksx][iksy][iksz] - ff[ix][iy - 1][iksx][iksy][iksz]);
                            }

                            if (iksy < nksy / 2) {//iksx<nksx/2; ksx<0 нуля нет

                                ff_12[ix][iy][iksx][iksy][iksz] = ff[ix][iy][iksx][iksy][iksz]
                                    - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[ix][iy + 1][iksx][iksy][iksz] - ff[ix][iy][iksx][iksy][iksz])
                                    - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[ix][iy][iksx][iksy][iksz] - ff[ix - 1][iy][iksx][iksy][iksz]);
                            }

                        }
                        else  ff_12[ix][iy][iksx][iksy][iksz] = ff[ix][iy][iksx][iksy][iksz];
                    }

                }
            }  //нашли для прилетающих распреление новое в конкретной точке x,y
    //находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn
            ron = 0;
            vxn = 0;
            vyn = 0;
            for (iksx = nksx / 2; iksx < nksx; iksx++) {
                for (iksy = 0; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {
                        ron = ron + ff_12[ix][iy][iksx][iksy][iksz];
                        vxn = vxn + ff_12[ix][iy][iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);//ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                        vyn = vyn + ff_12[ix][iy][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                    }
                }
            }


            //тут надо подумать че будет если ron=0

            if ((ron) == 0) { vx = 0; vy = 0; }
            if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }

            if (abs(vxn) > 0)      ron = 1 * 2 * abs(vxn) * sqrt(ms * M_PI / (2 * T1));


            for (iksx = 0; iksx < nksx / 2; iksx++) {
                for (iksy = 0; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        ff_12[ix][iy][iksx][iksy][iksz] = 1 * ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 1);
                    }
                }
            }

            vxnp = 0;
            // считаем какой будет поток при koefotr=1
            for (iksx = 0; iksx < nksx / 2; iksx++) {
                for (iksy = 0; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        vxnp = vxnp + ff_12[ix][iy][iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);

                    }
                }
            }





            if (abs(vxnp) > 0) {
                koefotr = abs(vxn) / abs(vxnp);
            }
            else koefotr = 1;
            // задаем новое отражение для консерватизма


            if (abs(vxn) > 0)      ron = koefotr * 2 * abs(vxn) * sqrt(ms * M_PI / (2 * T1));


            for (iksx = 0; iksx < nksx / 2; iksx++) {
                for (iksy = 0; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        ff_12[ix][iy][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 1);
                    }
                }
            }




            vxnp = 0;

        }


        ////////////////////////////////////////////конец левой стенки трубы///////////////

        // нижняя стенка левее чипа
        for (ix = 0; ix < bias - 1; ix++) {
            // нижняя стенка левее чипа

            for (iksx = 0; iksx < nksx; iksx++) {         //< 
                for (iksy = 0; iksy < nksy / 2; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {
                        if (ix > 1) {

                            //ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                            //ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;

                            if (iksx >= nksx / 2) {//iksx>nksx/2; ksx>0

                                ff_12[ix][0][iksx][iksy][iksz] = ff[ix][0][iksx][iksy][iksz]
                                    - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[ix][0 + 1][iksx][iksy][iksz] - ff[ix][0][iksx][iksy][iksz])
                                    - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[ix][0][iksx][iksy][iksz] - ff[ix - 1][0][iksx][iksy][iksz]);
                            }

                            if (iksx < nksx / 2) {//iksx<nksx/2; ksx<0 нуля нет

                                ff_12[ix][0][iksx][iksy][iksz] = ff[ix][0][iksx][iksy][iksz]
                                    - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[ix][0 + 1][iksx][iksy][iksz] - ff[ix][0][iksx][iksy][iksz])
                                    - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[ix + 1][0][iksx][iksy][iksz] - ff[ix][0][iksx][iksy][iksz]);
                            }
                        }
                        else ff_12[ix][0][iksx][iksy][iksz] = ff[ix][0][iksx][iksy][iksz];
                    }
                }
            }
            //нашли для прилетающих распреление новое в конкретной точке x,y
 //находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn
            ron = 0;
            vxn = 0;
            vyn = 0;
            for (iksx = 0; iksx < nksx; iksx++) {
                for (iksy = 0; iksy < nksy / 2; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {
                        ron = ron + ff_12[ix][0][iksx][iksy][iksz];
                        vxn = vxn + ff_12[ix][0][iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);//ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                        vyn = vyn + ff_12[ix][0][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                    }
                }
            }


            //тут надо подумать че будет если ron=0

            if ((ron) == 0) { vx = 0; vy = 0; }
            if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }
            //}


          //  for (ix = 1; ix < bias; ix++) {


            if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }
            // здесь  пиши distribution ff_12[ix][ny][]

            if (abs(vyn) > 0)      ron = 1 * 2 * abs(vyn) * sqrt(ms * M_PI / (2 * T1));


            for (iksx = 0; iksx < nksx; iksx++) {
                for (iksy = nksy / 2; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        ff_12[ix][0][iksx][iksy][iksz] = 1 * ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 1);
                    }
                }
            }

            vynp = 0;
            // считаем какой будет поток при koefotr=1
            for (iksx = 0; iksx < nksx; iksx++) {
                for (iksy = nksy / 2; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        vynp = vynp + ff_12[ix][0][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);

                    }
                }
            }





            if (abs(vynp) > 0) {
                koefotr = abs(vyn) / abs(vynp);
            }
            else koefotr = 1;
            // задаем новое отражение для консерватизма


            if (abs(vyn) > 0)      ron = koefotr * 2 * abs(vyn) * sqrt(ms * M_PI / (2 * T1));


            for (iksx = 0; iksx < nksx; iksx++) {
                for (iksy = nksy / 2; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        ff_12[ix][0][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 1);
                    }
                }
            }




            vynp = 0;

        }


        // еще одна точка прлилегающая к левому нижнему углу


        for (iksx = 0; iksx < nksx; iksx++) {         //< 
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    //ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                    //ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;

                    if (iksx >= nksx / 2) {//iksx>nksx/2; ksx>0

                        ff_12[bias - 1][0][iksx][iksy][iksz] = ff[bias - 1][0][iksx][iksy][iksz] - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias - 1][0 + 1][iksx][iksy][iksz] - ff[bias - 1][0][iksx][iksy][iksz])
                            - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias - 1][0][iksx][iksy][iksz] - ff[bias - 1 - 1][0][iksx][iksy][iksz]);
                    }

                    if (iksx < nksx / 2) {//iksx<nksx/2; ksx<0 нуля нет

                        ff_12[bias - 1][0][iksx][iksy][iksz] = ff[bias - 1][0][iksx][iksy][iksz] - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias - 1][0 + 1][iksx][iksy][iksz] - ff[bias - 1][0][iksx][iksy][iksz])
                            - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (fexcp5[iksx][iksy][iksz] - ff[bias - 1][0][iksx][iksy][iksz]);
                    }

                }
            }
        }  //нашли для прилетающих распреление новое в конкретной точке x,y
//находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn
        ron = 0;
        vxn = 0;
        vyn = 0;
        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    ron = ron + ff_12[bias - 1][0][iksx][iksy][iksz];
                    vxn = vxn + ff_12[bias - 1][0][iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);//ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                    vyn = vyn + ff_12[bias - 1][0][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                }
            }
        }


        //тут надо подумать че будет если ron=0

        if ((ron) == 0) { vx = 0; vy = 0; }
        if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }
        //}


      //  for (ix = 1; ix < bias; ix++) {

        if ((ron) == 0) { vxn = 0; vyn = 0; }
        if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }
        // здесь  пиши distribution ff_12[ix][ny][]

        if (abs(vyn) > 0)      ron = koefotr * 2 * abs(vyn) * sqrt(ms * M_PI / (2 * T1));


        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = nksy / 2; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    ff_12[bias - 1][0][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 1);

                }
            }
        }






        //end нижняя стенка левее чипа

    // нижняя стенка правее чипа

        for (ix = (bias + height * ratio + 2); ix < nx; ix++) {
            // нижняя стенка правее чипа
            for (iksx = 0; iksx < nksx; iksx++) {         //< ? или <=
                for (iksy = 0; iksy < nksy / 2; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        //ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                        //ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;

                        if (iksx >= nksx / 2) {//iksx>nksx/2; ksx>0

                            ff_12[ix][0][iksx][iksy][iksz] = ff[ix][0][iksx][iksy][iksz] - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[ix][0 + 1][iksx][iksy][iksz] - ff[ix][0][iksx][iksy][iksz]) - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[ix][0][iksx][iksy][iksz] - ff[ix - 1][0][iksx][iksy][iksz]);
                        }

                        if (iksx < nksx / 2) {//iksx<nksx/2; ksx<0 нуля нет

                            ff_12[ix][0][iksx][iksy][iksz] = ff[ix][0][iksx][iksy][iksz] - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[ix][0 + 1][iksx][iksy][iksz] - ff[ix][0][iksx][iksy][iksz]) - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[ix + 1][0][iksx][iksy][iksz] - ff[ix][0][iksx][iksy][iksz]);
                        }

                    }
                }
            }  //нашли для прилетающих распреление новое в конкретной точке x,y
    //находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn
            ron = 0;
            vxn = 0;
            vyn = 0;
            for (iksx = 0; iksx < nksx; iksx++) {
                for (iksy = 0; iksy < nksy / 2; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {
                        ron = ron + ff_12[ix][0][iksx][iksy][iksz];
                        vxn = vxn + ff_12[ix][0][iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);//ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                        vyn = vyn + ff_12[ix][0][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                    }
                }
            }


            //тут надо подумать че будет если ron=0

            if ((ron) == 0) { vx = 0; vy = 0; }
            if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }


            if (abs(vyn) > 0)      ron = 1 * 2 * abs(vyn) * sqrt(ms * M_PI / (2 * T1));


            for (iksx = 0; iksx < nksx; iksx++) {
                for (iksy = nksy / 2; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        ff_12[ix][0][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 1);

                    }
                }
            }
            // for (ix = (bias + height * ratio + 1); ix < nx; ix++) {
            vynp = 0;
            // считаем какой будет поток при koefotr=1
            for (iksx = 0; iksx < nksx; iksx++) {
                for (iksy = nksy / 2; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        vynp = vynp + ff_12[ix][0][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);

                    }
                }
            }





            if (abs(vynp) > 0) {
                koefotr = abs(vyn) / abs(vynp);
            }
            else koefotr = 1;
            // задаем новое отражение для консерватизма


            if (abs(vyn) > 0)      ron = koefotr * 2 * abs(vyn) * sqrt(ms * M_PI / (2 * T1));


            for (iksx = 0; iksx < nksx; iksx++) {
                for (iksy = nksy / 2; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        ff_12[ix][0][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 1);
                    }
                }
            }




            vynp = 0;
            // нижняя стенка правее чипа













        }


        // еще одна точка прлилегающая к правомуму нижнему углу

        for (iksx = 0; iksx < nksx; iksx++) {         //< ? или <=
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    //ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                    //ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;

                    if (iksx >= nksx / 2) {//iksx>nksx/2; ksx>0

                        ff_12[bias + height * ratio + 1][0][iksx][iksy][iksz] = ff[bias + height * ratio + 1][0][iksx][iksy][iksz]
                            - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias + height * ratio + 1][0 + 1][iksx][iksy][iksz] - ff[bias + height * ratio + 1][0][iksx][iksy][iksz])
                            - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias + height * ratio + 1][0][iksx][iksy][iksz] - fexcp8[iksx][iksy][iksz]);
                    }

                    if (iksx < nksx / 2) {//iksx<nksx/2; ksx<0 нуля нет

                        ff_12[bias + height * ratio + 1][0][iksx][iksy][iksz] = ff[bias + height * ratio + 1][0][iksx][iksy][iksz]
                            - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias + height * ratio + 1][0 + 1][iksx][iksy][iksz] - ff[bias + height * ratio + 1][0][iksx][iksy][iksz])
                            - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias + height * ratio + 2][0][iksx][iksy][iksz] - ff[bias + height * ratio + 1][0][iksx][iksy][iksz]);
                    }
                }
            }
        }  //нашли для прилетающих распреление новое в конкретной точке x,y
//находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn
        ron = 0;
        vxn = 0;
        vyn = 0;
        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    ron = ron + ff_12[bias + height * ratio + 1][0][iksx][iksy][iksz];
                    vxn = vxn + ff_12[bias + height * ratio + 1][0][iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);//ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                    vyn = vyn + ff_12[bias + height * ratio + 1][0][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                }
            }
        }

        //тут надо подумать че будет если ron=0

        if ((ron) == 0) { vx = 0; vy = 0; }
        if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }


        if (abs(vyn) > 0)      ron = koefotr * 2 * abs(vyn) * sqrt(ms * M_PI / (2 * T1));


        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = nksy / 2; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    ff_12[bias + height * ratio + 1][0][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 1);
                }
            }
        }
        // for (ix = (bias + height * ratio + 1); ix < nx; ix++) {
             // нижняя стенка правее чипа












        // end нижняя стенка правее чипа

        // верхняя стенка чипа

        for (ix = bias + 1 + 1; ix < (bias + height * ratio) - 1; ix++) {



            for (iksx = 0; iksx < nksx; iksx++) {         //< ? или <=
                for (iksy = 0; iksy < nksy / 2; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        //ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                        //ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;

                        if (iksx >= nksx / 2) {//iksx>nksx/2; ksx>0

                            ff_12[ix][height][iksx][iksy][iksz] = ff[ix][height][iksx][iksy][iksz] - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[ix][height + 1][iksx][iksy][iksz] - ff[ix][height][iksx][iksy][iksz])
                                - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[ix][height][iksx][iksy][iksz] - ff[ix - 1][height][iksx][iksy][iksz]);
                        }

                        if (iksx < nksx / 2) {//iksx<nksx/2; ksx<0 нуля нет

                            ff_12[ix][height][iksx][iksy][iksz] = ff[ix][height][iksx][iksy][iksz] - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[ix][height + 1][iksx][iksy][iksz] - ff[ix][height][iksx][iksy][iksz])
                                - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[ix + 1][height][iksx][iksy][iksz] - ff[ix][height][iksx][iksy][iksz]);
                        }

                    }
                }
            }  //нашли для прилетающих распреление новое в конкретной точке x,y
    //находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn
            ron = 0;
            vxn = 0;
            vyn = 0;
            for (iksx = 0; iksx < nksx; iksx++) {
                for (iksy = 0; iksy < nksy / 2; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {
                        ron = ron + ff_12[ix][height][iksx][iksy][iksz];
                        vxn = vxn + ff_12[ix][height][iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);//ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                        vyn = vyn + ff_12[ix][height][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                    }
                }
            }

            //   }

            if ((ron) == 0) { vx = 0; vy = 0; }
            if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }

            // верхняя стенка чипа
            if (abs(vyn) > 0)      ron = 1 * 2 * abs(vyn) * sqrt(ms * M_PI / (2 * T2));


            for (iksx = 0; iksx < nksx; iksx++) {
                for (iksy = nksy / 2; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        ff_12[ix][height][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 2);

                    }
                }
            }


            vynp = 0;
            // считаем какой будет поток при koefotr=1
            for (iksx = 0; iksx < nksx; iksx++) {
                for (iksy = nksy / 2; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        vynp = vynp + ff_12[ix][height][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);

                    }
                }
            }





            if (abs(vynp) > 0) {
                koefotr = abs(vyn) / abs(vynp);
            }
            else koefotr = 1;
            // задаем новое отражение для консерватизма


            if (abs(vyn) > 0)      ron = koefotr * 2 * abs(vyn) * sqrt(ms * M_PI / (2 * T2));


            for (iksx = 0; iksx < nksx; iksx++) {
                for (iksy = nksy / 2; iksy < nksy; iksy++) {
                    for (iksz = 0; iksz < nksz; iksz++) {

                        ff_12[ix][height][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 2);
                    }
                }
            }




            vynp = 0;




        }


        //  еще одна точка прилегающая к левому верхнему на верхней стенке чипа


        for (iksx = 0; iksx < nksx; iksx++) {         //< ? или <=
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    //ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                    //ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;

                    if (iksx >= nksx / 2) {//iksx>nksx/2; ksx>0

                        ff_12[bias + 1][height][iksx][iksy][iksz] = ff[bias + 1][height][iksx][iksy][iksz]
                            - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias + 1][height + 1][iksx][iksy][iksz] - ff[bias + 1][height][iksx][iksy][iksz])
                            - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias + 1][height][iksx][iksy][iksz] - fexcp2[iksx][iksy][iksz]);
                    }

                    if (iksx < nksx / 2) {//iksx<nksx/2; ksx<0 нуля нет

                        ff_12[bias + 1][height][iksx][iksy][iksz] = ff[bias + 1][height][iksx][iksy][iksz]
                            - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias + 1][height + 1][iksx][iksy][iksz] - ff[bias + 1][height][iksx][iksy][iksz])
                            - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias + 1 + 1][height][iksx][iksy][iksz] - ff[bias + 1][height][iksx][iksy][iksz]);
                    }
                }
            }
        }  //нашли для прилетающих распреление новое в конкретной точке x,y
//находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn

        ix = bias + 1;

        ron = 0;
        vxn = 0;
        vyn = 0;
        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    ron = ron + ff_12[ix][height][iksx][iksy][iksz];
                    vxn = vxn + ff_12[ix][height][iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);//ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                    vyn = vyn + ff_12[ix][height][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                }
            }
        }

        //   }

        if ((ron) == 0) { vx = 0; vy = 0; }
        if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }

        // верхняя стенка чипа
        if (abs(vyn) > 0)      ron = 1 * 2 * abs(vyn) * sqrt(ms * M_PI / (2 * T2));


        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = nksy / 2; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    ff_12[ix][height][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 2);

                }
            }
        }


        vynp = 0;
        // считаем какой будет поток при koefotr=1
        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = nksy / 2; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    vynp = vynp + ff_12[ix][height][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);

                }
            }
        }





        if (abs(vynp) > 0) {
            koefotr = abs(vyn) / abs(vynp);
        }
        else koefotr = 1;
        // задаем новое отражение для консерватизма


        if (abs(vyn) > 0)      ron = koefotr * 2 * abs(vyn) * sqrt(ms * M_PI / (2 * T2));


        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = nksy / 2; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    ff_12[ix][height][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 2);
                }
            }
        }




        vynp = 0;


        //end  еще одна точка прилегающая к левому верхнему на верхней стенке чипа


        //  еще одна точка прилегающая к правому верхнему на верхней стенке чипа


        for (iksx = 0; iksx < nksx; iksx++) {         //< ? или <=
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    //ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                    //ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;

                    if (iksx >= nksx / 2) {//iksx>nksx/2; ksx>0

                        ff_12[bias + height * ratio - 1][height][iksx][iksy][iksz] = ff[bias + height * ratio - 1][height][iksx][iksy][iksz]
                            - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias + height * ratio - 1][height + 1][iksx][iksy][iksz] - ff[bias + height * ratio - 1][height][iksx][iksy][iksz])
                            - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias + height * ratio - 1][height][iksx][iksy][iksz] - ff[bias + height * ratio - 1 - 1][height][iksx][iksy][iksz]);
                    }

                    if (iksx < nksx / 2) {//iksx<nksx/2; ksx<0 нуля нет

                        ff_12[bias + height * ratio - 1][height][iksx][iksy][iksz] = ff[bias + height * ratio - 1][height][iksx][iksy][iksz]
                            - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias + height * ratio - 1][height + 1][iksx][iksy][iksz] - ff[bias + height * ratio - 1][height][iksx][iksy][iksz])
                            - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (fexcp3[iksx][iksy][iksz] - ff[bias + height * ratio - 1][height][iksx][iksy][iksz]);
                    }

                }
            }
        }  //нашли для прилетающих распреление новое в конкретной точке x,y
//находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn
        ix = bias + height * ratio - 1;
        ron = 0;
        vxn = 0;
        vyn = 0;
        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    ron = ron + ff_12[ix][height][iksx][iksy][iksz];
                    vxn = vxn + ff_12[ix][height][iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);//ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                    vyn = vyn + ff_12[ix][height][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                }
            }
        }

        //   }

        if ((ron) == 0) { vx = 0; vy = 0; }
        if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }

        // верхняя стенка чипа
        if (abs(vyn) > 0)      ron = 1 * 2 * abs(vyn) * sqrt(ms * M_PI / (2 * T2));


        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = nksy / 2; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    ff_12[ix][height][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 2);

                }
            }
        }


        vynp = 0;
        // считаем какой будет поток при koefotr=1
        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = nksy / 2; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    vynp = vynp + ff_12[ix][height][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);

                }
            }
        }





        if (abs(vynp) > 0) {
            koefotr = abs(vyn) / abs(vynp);
        }
        else koefotr = 1;
        // задаем новое отражение для консерватизма


        if (abs(vyn) > 0)      ron = koefotr * 2 * abs(vyn) * sqrt(ms * M_PI / (2 * T2));


        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = nksy / 2; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    ff_12[ix][height][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 2);
                }
            }
        }




        vynp = 0;


        //end  еще одна точка прилегающая к левому верхнему на верхней стенке чипа


        // end верхняя стенка чипа




          // левая стенка чипа


        // левая стенка чипа    здесь потому что вертикальная стенка  из 2-х точки, обрабатываем как прилегающие к угловым без цикла                     







        for (iksx = nksx / 2; iksx < nksx; iksx++) {  // скорость вправо
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {


                    if (iksy >= nksx / 2) {//iksy>nksx/2; ksy>0

                        ff_12[bias][1][iksx][iksy][iksz] = ff[bias][1][iksx][iksy][iksz]
                            - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias][1][iksx][iksy][iksz] - fexcp6[iksx][iksy][iksz])
                            - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias][1][iksx][iksy][iksz] - ff[bias - 1][1][iksx][iksy][iksz]);
                    }

                    if (iksy < nksx / 2) {//iksy<nksx/2; ksy<0 нуля нет

                        ff_12[bias][1][iksx][iksy][iksz] = ff[bias][1][iksx][iksy][iksz] - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias][1 + 1][iksx][iksy][iksz] - ff[bias][1][iksx][iksy][iksz])
                            - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias][1][iksx][iksy][iksz] - ff[bias - 1][1][iksx][iksy][iksz]);
                    }



                }
            }
        }


        //нашли для прилетающих распреление новое в конкретной точке x,y
    //находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn

        ron = 0;
        vxn = 0;
        vyn = 0;
        for (iksx = nksx / 2; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    ron = ron + ff_12[bias][1][iksx][iksy][iksz];
                    vxn = vxn + ff_12[bias][1][iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);//ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                    vyn = vyn + ff_12[bias][1][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                }
            }
        }




        if ((ron) == 0) { vx = 0; vy = 0; }
        if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }
        //здесь касательная vyn нормальная vxn

   // }
        if (abs(vxn) > 0)      ron = koefotr * 2 * abs(vxn) * sqrt(ms * M_PI / (2 * T2));


        for (iksx = 0; iksx < nksx / 2; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    ff_12[bias][1][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 2);
                }
            }
        }








        for (iksx = nksy / 2; iksx < nksx; iksx++) {  // скорость вправо
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {


                    if (iksy >= nksx / 2) {//iksy>nksx/2; ksy>0

                        ff_12[bias][2][iksx][iksy][iksz] = ff[bias][2][iksx][iksy][iksz]
                            - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias][2][iksx][iksy][iksz] - ff[bias][2 - 1][iksx][iksy][iksz])
                            - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias][2][iksx][iksy][iksz] - ff[bias - 1][2][iksx][iksy][iksz]);
                    }

                    if (iksy < nksx / 2) {//iksy<nksx/2; ksy<0 нуля нет

                        ff_12[bias][2][iksx][iksy][iksz] = ff[bias][2][iksx][iksy][iksz]
                            - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (fexcp1[iksx][iksy][iksz] - ff[bias][2][iksx][iksy][iksz])
                            - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias][2][iksx][iksy][iksz] - ff[bias - 1][2][iksx][iksy][iksz]);
                    }

                }

            }
        }


        //нашли для прилетающих распреление новое в конкретной точке x,y
    //находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn

        ron = 0;
        vxn = 0;
        vyn = 0;
        for (iksx = nksx / 2; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    ron = ron + ff_12[bias][2][iksx][iksy][iksz];
                    vxn = vxn + ff_12[bias][2][iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);//ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                    vyn = vyn + ff_12[bias][2][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                }
            }
        }




        if ((ron) == 0) { vx = 0; vy = 0; }
        if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }
        //здесь касательная vyn нормальная vxn

   // }
        if (abs(vxn) > 0)      ron = koefotr * 2 * abs(vxn) * sqrt(ms * M_PI / (2 * T2));


        for (iksx = 0; iksx < nksx / 2; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    ff_12[bias][2][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 2);
                }
            }
        }










        // end левая стенка чипа


        // end левая стенка чипа

















        //================================================26,05==========================================


            // правая стенка чипа    здесь потому что вертикальная стенка  из 2-х точки, обрабатываем как прилегающие к угловым без цикла                     


        for (iksx = 0; iksx < nksx / 2; iksx++) {  // скорость влево
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    if (iksy >= nksx / 2) {//iksy>nksx/2; ksy>0

                        ff_12[bias + height * ratio][2][iksx][iksy][iksz] = ff[bias + height * ratio][2][iksx][iksy][iksz] - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias + height * ratio][2][iksx][iksy][iksz] - ff[bias + height * ratio][2 - 1][iksx][iksy][iksz])
                            - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias + height * ratio + 1][2][iksx][iksy][iksz] - ff[bias + height * ratio][2][iksx][iksy][iksz]);
                    }

                    if (iksy < nksx / 2) {//iksy<nksx/2; ksy<0 нуля нет

                        ff_12[bias + height * ratio][2][iksx][iksy][iksz] = ff[bias + height * ratio][2][iksx][iksy][iksz] - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (fexcp4[iksx][iksy][iksz] - ff[bias + height * ratio][2][iksx][iksy][iksz])
                            - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias + height * ratio + 1][2][iksx][iksy][iksz] - ff[bias + height * ratio][2][iksx][iksy][iksz]);
                    }
                }


            }
        }


        //нашли для прилетающих распреление новое в конкретной точке x,y
    //находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn

        ron = 0;
        vxn = 0;
        vyn = 0;
        for (iksx = 0; iksx < nksx / 2; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    ron = ron + ff_12[bias + height * ratio][2][iksx][iksy][iksz];
                    vxn = vxn + ff_12[bias + height * ratio][2][iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);//ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                    vyn = vyn + ff_12[bias + height * ratio][2][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                }
            }
        }


        if ((ron) == 0) { vxn = 0; vyn = 0; }
        if (abs(ron) <= 0.000000000000001) { vxn = 0; vyn = 0; }

        if (abs(vxn) > 0)      ron = koefotr * 2 * abs(vxn) * sqrt(ms * M_PI / (2 * T2));


        for (iksx = nksx / 2; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    ff_12[bias + height * ratio][2][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 2);
                }
            }
        }






        for (iksx = 0; iksx < nksx / 2; iksx++) {  // скорость влево
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    if (iksy >= nksx / 2) {//iksy>nksx/2; ksy>0

                        ff_12[bias + height * ratio][1][iksx][iksy][iksz] = ff[bias + height * ratio][1][iksx][iksy][iksz]
                            - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias + height * ratio][1][iksx][iksy][iksz] - fexcp7[iksx][iksy][iksz])
                            - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias + height * ratio + 1][1][iksx][iksy][iksz] - ff[bias + height * ratio][1][iksx][iksy][iksz]);
                    }

                    if (iksy < nksx / 2) {//iksy<nksx/2; ksy<0 нуля нет

                        ff_12[bias + height * ratio][1][iksx][iksy][iksz] = ff[bias + height * ratio][1][iksx][iksy][iksz]
                            - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias + height * ratio][1 + 1][iksx][iksy][iksz] - ff[bias + height * ratio][1][iksx][iksy][iksz])
                            - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias + height * ratio + 1][1 + 1][iksx][iksy][iksz] - ff[bias + height * ratio][1][iksx][iksy][iksz]);
                    }



                }
            }
        }


        //нашли для прилетающих распреление новое в конкретной точке x,y
    //находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn

        ron = 0;
        vxn = 0;
        vyn = 0;
        for (iksx = 0; iksx < nksx / 2; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    ron = ron + ff_12[bias + height * ratio][1][iksx][iksy][iksz];
                    vxn = vxn + ff_12[bias + height * ratio][1][iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);//ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                    vyn = vyn + ff_12[bias + height * ratio][1][iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                }
            }
        }


        if (abs(ron) <= 0.000000000000001) { vxn = 0; vyn = 0; }

        if (abs(vxn) > 0)      ron = koefotr * 2 * abs(vxn) * sqrt(ms * M_PI / (2 * T2));


        for (iksx = nksx / 2; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    ff_12[bias + height * ratio][1][iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 2);

                }
            }
        }














        //end правая стенка чипа






         //левый верхний угол  1

        for (iksx = nksx / 2; iksx < nksx; iksx++) {
            for (iksy = nksy / 2; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    f1excp1[iksx][iksy][iksz] = fexcp1[iksx][iksy][iksz]
                        - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (fexcp1[iksx][iksy][iksz] - ff[bias - 1][height][iksx][iksy][iksz])
                        - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (fexcp1[iksx][iksy][iksz] - ff[bias][height - 1][iksx][iksy][iksz]);


                }
            }
        }


        for (iksx = nksx / 2; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    f1excp1[iksx][iksy][iksz] = fexcp1[iksx][iksy][iksz]
                        - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (fexcp1[iksx][iksy][iksz] - ff[bias - 1][height][iksx][iksy][iksz])
                        - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias][height + 1][iksx][iksy][iksz] - fexcp1[iksx][iksy][iksz]);

                }

            }
        }




        //нашли для прилетающих распреление новое в конкретной точке x,y
    //находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn

        ron = 0;
        vxn = 0;
        vyn = 0;
        for (iksx = nksx / 2; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    ron = ron + f1excp1[iksx][iksy][iksz];
                    vxn = vxn + f1excp1[iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);//ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                    vyn = vyn + f1excp1[iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                }
            }
        }


        if ((ron) == 0) { vx = 0; vy = 0; }
        if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }


        if (abs(vxn) > 0)      ron = koefotr * 2 * abs(vxn) * sqrt(ms * M_PI / (2 * T2));




        // fff[ix][height][iksx][iksy][iksz] = distribution_new;
        for (iksx = 0; iksx < nksx / 2; iksx++) {         //< ? или <=
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    f1excp1[iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 2);
                }
            }
        }


        //левый верхний угол  2

        for (iksx = nksx / 2; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    f1excp2[iksx][iksy][iksz] = fexcp2[iksx][iksy][iksz]
                        - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (fexcp2[iksx][iksy][iksz] - ff[bias - 1][height][iksx][iksy][iksz])
                        - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias][height + 1][iksx][iksy][iksz] - fexcp2[iksx][iksy][iksz]);
                }
            }
        }

        for (iksx = 0; iksx < nksx / 2; iksx++) {
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    f1excp2[iksx][iksy][iksz] = fexcp2[iksx][iksy][iksz]
                        - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias + 1][height][iksx][iksy][iksz] - fexcp2[iksx][iksy][iksz])
                        - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias][height + 1][iksx][iksy][iksz] - fexcp2[iksx][iksy][iksz]);
                }
            }
        }


        //нашли для прилетающих распреление новое в конкретной точке x,y
    //находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn

        ron = 0;
        vxn = 0;
        vyn = 0;
        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    ron = ron + f1excp2[iksx][iksy][iksz];
                    vxn = vxn + f1excp2[iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);//ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                    vyn = vyn + f1excp2[iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                }
            }
        }


        if ((ron) == 0) { vx = 0; vy = 0; }
        if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }

        if (abs(vyn) > 0)      ron = koefotr * 2 * abs(vyn) * sqrt(ms * M_PI / (2 * T2));

        for (iksx = 0; iksx < nksx; iksx++) {         //вверх  отражается
            for (iksy = nksy / 2; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    f1excp2[iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 2);

                }
            }
        }



        //правый верхний угол  3

        for (iksx = nksx / 2; iksx < nksx; iksx++) {         //< ? или <=
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    f1excp3[iksx][iksy][iksz] = fexcp3[iksx][iksy][iksz]
                        - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (fexcp3[iksx][iksy][iksz] - ff[bias + height * ratio - 1][height][iksx][iksy][iksz])
                        - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias + height * ratio][height + 1][iksx][iksy][iksz] - fexcp3[iksx][iksy][iksz]);

                }

            }
        }

        for (iksx = 0; iksx < nksx / 2; iksx++) {         //< ? или <=
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    f1excp3[iksx][iksy][iksz] = fexcp3[iksx][iksy][iksz]
                        - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias + height * ratio + 1][height][iksx][iksy][iksz] - fexcp3[iksx][iksy][iksz])
                        - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias + height * ratio][height + 1][iksx][iksy][iksz] - fexcp3[iksx][iksy][iksz]);
                }


            }
        }



        //нашли для прилетающих распреление новое в конкретной точке x,y
    //находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn

        ron = 0;
        vxn = 0;
        vyn = 0;
        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    ron = ron + f1excp3[iksx][iksy][iksz];
                    vxn = vxn + f1excp3[iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);//ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                    vyn = vyn + f1excp3[iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                }
            }
        }


        if ((ron) == 0) { vx = 0; vy = 0; }
        if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }

        if (abs(vyn) > 0)      ron = koefotr * 2 * abs(vyn) * sqrt(ms * M_PI / (2 * T2));

        for (iksx = 0; iksx < nksx; iksx++) {//вверх  отражается
            for (iksy = nksy / 2; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    f1excp3[iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 2);
                }
            }
        }













        //правый верхний угол  4

        for (iksx = 0; iksx < nksx / 2; iksx++) {         //

            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    f1excp4[iksx][iksy][iksz] = fexcp4[iksx][iksy][iksz]
                        - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias + 1 + height * ratio][height][iksx][iksy][iksz] - fexcp4[iksx][iksy][iksz])
                        - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias + height * ratio][height + 1][iksx][iksy][iksz] - fexcp4[iksx][iksy][iksz]);


                }
            }
        }

        for (iksx = 0; iksx < nksx / 2; iksx++) {         //

            for (iksy = nksy / 2; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    f1excp4[iksx][iksy][iksz] = fexcp4[iksx][iksy][iksz]
                        - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias + 1 + height * ratio][height][iksx][iksy][iksz] - fexcp4[iksx][iksy][iksz])
                        - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (fexcp4[iksx][iksy][iksz] - ff[bias + height * ratio][height - 1][iksx][iksy][iksz]);


                }
            }
        }


        //нашли для прилетающих распреление новое в конкретной точке x,y
    //находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn

        ron = 0;
        vxn = 0;
        vyn = 0;
        for (iksx = 0; iksx < nksx / 2; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    ron = ron + f1excp4[iksx][iksy][iksz];
                    vxn = vxn + f1excp4[iksx][iksy][iksz] * spd_par * (2 * (iksx - nksx / 2) + 1);//ksx=iksx-nksx/2  например  скорость ksx=-5;  iksx=-5+40/2=15 ;
                    vyn = vyn + f1excp4[iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                }
            }
        }


        if ((ron) == 0) { vx = 0; vy = 0; }
        if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }

        if (abs(vxn) > 0)      ron = koefotr * 2 * abs(vxn) * sqrt(ms * M_PI / (2 * T2));



        for (iksx = nksx / 2; iksx < nksx; iksx++) {                     //вправо  отражается
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    f1excp4[iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 2);
                }
            }
        }






        //левый нижний угол  5,6

        for (iksx = nksx / 2; iksx < nksx; iksx++) {         //< ? или <=
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {



                    f1excp5[iksx][iksy][iksz] = fexcp5[iksx][iksy][iksz]
                        - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (fexcp5[iksx][iksy][iksz] - ff[bias - 1][0][iksx][iksy][iksz])
                        - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias][0 + 1][iksx][iksy][iksz] - fexcp5[iksx][iksy][iksz]);


                }

            }
        }

        for (iksx = nksx / 2; iksx < nksx; iksx++) {         //< ? или <=
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {



                    f1excp6[iksx][iksy][iksz] = fexcp6[iksx][iksy][iksz]
                        - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias][0 + 1][iksx][iksy][iksz] - fexcp6[iksx][iksy][iksz])
                        - spd_par * (dt1 / dx) * (2 * (iksx - nksy / 2) + 1) * (fexcp6[iksx][iksy][iksz] - ff[bias - 1][0][iksx][iksy][iksz]);
                }
            }
        }




        //нашли для прилетающих распреление новое в конкретной точке x,y
    //находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn

        ron = 0;
        vxn = 0;
        vyn = 0;

        for (iksx = nksx / 2; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    vxn = vxn + f1excp5[iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                }
            }
        }


        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    vyn = vyn + f1excp6[iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                }
            }
        }








        if ((ron) == 0) { vx = 0; vy = 0; }
        if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }

        if (abs(vxn) > 0)      ron = koefotr * 2 * 2 * abs(vyn) * sqrt(ms * M_PI / (2 * T2));

        for (iksx = 0; iksx < nksx / 2; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {          // влево  отражается
                for (iksz = 0; iksz < nksz; iksz++) {
                    f1excp5[iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 2);
                }
            }
        }

        if (abs(vyn) > 0)      ron = koefotr * 2 * 2 * abs(vxn) * sqrt(ms * M_PI / (2 * T1));

        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = nksy / 2; iksy < nksy; iksy++) {          // вверх  отражается
                for (iksz = 0; iksz < nksz; iksz++) {
                    f1excp6[iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 1);
                }
            }
        }

        //правый нижний угол   7 , 8


        for (iksx = 0; iksx < nksx / 2; iksx++) {
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {
                    f1excp7[iksx][iksy][iksz] = fexcp7[iksx][iksy][iksz]
                        - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias + height * ratio][0 + 1][iksx][iksy][iksz] - fexcp7[iksx][iksy][iksz])
                        - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias + height * ratio + 1][0][iksx][iksy][iksz] - fexcp7[iksx][iksy][iksz]);
                }
            }
        }







        for (iksx = 0; iksx < nksx / 2; iksx++) {
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    f1excp8[iksx][iksy][iksz] = fexcp8[iksx][iksy][iksz]
                        - spd_par * (dt1 / dx) * (2 * (iksx - nksx / 2) + 1) * (ff[bias + height * ratio + 1][0][iksx][iksy][iksz] - fexcp8[iksx][iksy][iksz])
                        - spd_par * (dt1 / dx) * (2 * (iksy - nksy / 2) + 1) * (ff[bias + height * ratio][0 + 1][iksx][iksy][iksz] - fexcp8[iksx][iksy][iksz]);
                }
            }
        }




        //нашли для прилетающих распреление новое в конкретной точке x,y
    //находим параметры для отлетающих n,vx,vy я тут обозначил плотность ron ; vx*n обозначил vxn

        ron = 0;
        vxn = 0;
        vyn = 0;

        for (iksx = 0; iksx < nksx / 2; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    vxn = vxn + f1excp8[iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                }
            }
        }


        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy / 2; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    vyn = vyn + f1excp7[iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);//ksy=iksy-nksy/2  например  скорость ksy=-7;  iksy=-7+40/2=13 ;
                }
            }
        }








        if ((ron) == 0) { vx = 0; vy = 0; }
        if (abs(ron) <= 0.000000000000001) { vx = 0; vy = 0; }

        if (abs(vxn) > 0)      ron = koefotr * 2 * 2 * abs(vxn) * sqrt(ms * M_PI / (2 * T2));

        for (iksx = nksx / 2; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {          // вправо  отражается
                for (iksz = 0; iksz < nksz; iksz++) {
                    f1excp8[iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 2);
                }
            }
        }


        vxnp = 0;

        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = nksy / 2; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    vxnp = vxnp + f1excp8[iksx][iksy][iksz] * spd_par * (2 * (iksy - nksy / 2) + 1);
                }
            }
        }

        vxnp = 0;








        if (abs(vyn) > 0)      ron = koefotr * 2 * 2 * abs(vyn) * sqrt(ms * M_PI / (2 * T1));

        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = nksy / 2; iksy < nksy; iksy++) {          // вверх  отражается
                for (iksz = 0; iksz < nksz; iksz++) {
                    f1excp7[iksx][iksy][iksz] = ron * distribution_speed1(iksx, iksy, iksz, 0, 0, 1, 1);
                }
            }
        }






        //===========================================26.05================обновление здесь======================================





        //обновление массива ff
        for (ix = 0; ix <= nx; ix++) {
            for (iy = 0; iy <= ny; iy++) {
                for (iksx = 0; iksx < nksx; iksx++) {
                    for (iksy = 0; iksy < nksy; iksy++) {
                        for (iksz = 0; iksz < nksz; iksz++) {
                            ff[ix][iy][iksx][iksy][iksz] = ff_12[ix][iy][iksx][iksy][iksz];
                        }
                    }
                }
            }
        }






        for (iksx = 0; iksx < nksx; iksx++) {
            for (iksy = 0; iksy < nksy; iksy++) {
                for (iksz = 0; iksz < nksz; iksz++) {

                    fexcp1[iksx][iksy][iksz] = f1excp1[iksx][iksy][iksz];
                    fexcp2[iksx][iksy][iksz] = f1excp2[iksx][iksy][iksz];
                    fexcp3[iksx][iksy][iksz] = f1excp3[iksx][iksy][iksz];
                    fexcp4[iksx][iksy][iksz] = f1excp4[iksx][iksy][iksz];
                    fexcp5[iksx][iksy][iksz] = f1excp5[iksx][iksy][iksz];
                    fexcp6[iksx][iksy][iksz] = f1excp6[iksx][iksy][iksz];
                    fexcp7[iksx][iksy][iksz] = f1excp7[iksx][iksy][iksz];
                    fexcp8[iksx][iksy][iksz] = f1excp8[iksx][iksy][iksz];
                }
            }
        }


        //===========================================26.05======================================================




        //вставка под конец цикла по времени


    }
    // закрытие цикла по времени

    cout << "jjjj" << endl;
    return 0;



}









