#include <iostream>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

int calc_angle(int E) {    //rutherford formula

    int E = 1000000;

    int i=0;
    float k=0;

    float radian=0;
    float sin_4th=0;
    float cnst = 3.4e-14; //should be -54

    float verovatnoca = 0;

    for(i=0; i<90; i++) {

        radian=(i+1)*M_PI/180;
        sin_4th = pow(sin(radian), 4);

        k = cnst / sin_4th;
        verovatnoca = k / (pow(E,2)); //should be times 10e-40

        cout<<verovatnoca<<"   ";
    }
    return 0;
}


