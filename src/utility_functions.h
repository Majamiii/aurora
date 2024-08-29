#ifndef UTILITYFUNCTIONS
#define UTILITYFUNCTIONS

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <vector> 
#include <math.h>
#include <cmath>

using namespace std;



inline float gen_random() {
   // """Return a vector of (N/2 + N%2) elements sampled from a normal distribution """
    float rnd1 = (float)rand() / RAND_MAX;
    float rnd2 = (float)rand() / RAND_MAX;
    return sqrt( -2.0 * log(rnd1)) * cos( 2 * M_PI * rnd2) ;
//       ret.push_back( sqrt( -2.0 * log(rnd1)) * sin( 2 * M_PI * rnd2)   );
//    }
//    return ret;
}

#endif
