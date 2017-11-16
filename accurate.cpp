#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include <fstream>

using namespace std;

typedef numeric_limits< double > dbl;

long double t = 0;
const long double dt = 0.01;

const unsigned int N = 9;

struct reb_particle {
    long double x;
    long double y;
    long double z;
    long double vx;
    long double vy;
    long double vz;
    long double ax;
    long double ay;
    long double az;
    long double m;
};


// Initial conditions (in FP)
struct reb_particle p[9] = {
    {
        0.0021709922250528, 0.0057845061154043, -0.0001290326677066,
        -0.0003084904334499, 0.0003164862379414, 0.0000072860648107,
        0, 0, 0, 1.0000000000000000
    },
    {
        -0.1529074277548495, -0.4329649810759809, -0.0217536815870956,
        1.2130892048755062, -0.4636664872138580, -0.1492230266727991,
        0, 0, 0, 0.0000001660114153
    },
    {
        -0.7051385792282048, 0.1305062392874893, 0.0423980407616931,
        -0.2107118903711193, -1.1628741220859935, -0.0038067721592922,
        0, 0, 0, 0.0000024478382878
    },
    {
        0.8303864923760965, 0.5551748865431479, -0.0001556226179998,
        -0.5694403294004744, 0.8300359440285254, -0.0000250486216637,
        0, 0, 0, 0.0000030404326480
    },
    {
        -1.6007632981663540, 0.4507843866326728, 0.0485350310380760,
        -0.1874661855400607, -0.7140231189065021, -0.0103688562255236,
        0, 0, 0, 0.0000003227156038
    },
    {
        -4.5444724195553627, -2.9811209359531872, 0.1140115745580475,
        0.2354668506120313, -0.3459544002171689, -0.0038305410200901,
        0, 0, 0, 0.0009547919152112
    },
    {
        -0.2998316596246585, -10.0512228718170959, 0.1866942196718307,
        0.3063599906570191, -0.0107135147677418, -0.0120072161180579,
        0, 0, 0, 0.0002858856727222
    },
    {
        17.8418531053445939, 8.8433796310403689, -0.1982994964737093,
        -0.1032131635550300, 0.1941992816066720, 0.0020584917278455,
        0, 0, 0, 0.0000436624373583
    },
    {
        28.6228992820092181, -8.7910334836014847, -0.4786090163574258,
        0.0523633993793736, 0.1755278382196959, -0.0048214129381180,
        0, 0, 0, 0.0000515138377263
    },
};

static void drift(){
    for(unsigned int i=0; i<N; i++){
        p[i].x += (dt/2. * (long double)p[i].vx);
        p[i].y += (dt/2. * (long double)p[i].vy);
        p[i].z += (dt/2. * (long double)p[i].vz);
    }
}

static void kick(){
    for(unsigned int i=0; i<N; i++){
        p[i].vx += (dt*p[i].ax);
        p[i].vy += (dt*p[i].ax);
        p[i].vz += (dt*p[i].az);
    }
}

static void gravity(){
    for(unsigned int i=0; i<N; i++){
        p[i].ax = 0.;
        p[i].ay = 0.;
        p[i].az = 0.;
        for(unsigned int j=0; j<N; j++){
            if (i!=j){
                const long double dx = p[i].x - p[j].x;
                const long double dy = p[i].y - p[j].y;
                const long double dz = p[i].z - p[j].z;
                const long double _r = sqrt(dx*dx + dy*dy + dz*dz);
                const long double prefact = -1/(_r*_r*_r)*p[j].m;
                
                p[i].ax    += prefact*dx;
                p[i].ay    += prefact*dy;
                p[i].az    += prefact*dz;
            }
        }
    }
}

void janus_step(){

    drift();

    gravity();
    kick();

    drift();

    t += dt;

}

int main(){
    int i = 0;
    ofstream myfile ("accurate.txt");

    myfile.precision(dbl::max_digits10);
    while (t<2.*M_PI*1e1){ // 1 year
        i += 1;
        janus_step();
        myfile << "------------" << i << "th----------" << endl;
        for(int i=0;i<N;i++){
            myfile << p[i].vx << "," << p[i].vy << endl;
        }
    }
    myfile.close();
    return 1;
}