#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include <fstream>

using namespace std;

// Scale for conversion between FP and INT
const double scale_vel  = 1e-16; 
const double scale_pos  = 1e-16;
typedef numeric_limits< double > dbl;
// INT datatype used
#define REB_PARTICLE_INT_TYPE int64_t

struct reb_particle_int {
    REB_PARTICLE_INT_TYPE x;
    REB_PARTICLE_INT_TYPE y;
    REB_PARTICLE_INT_TYPE z;
    REB_PARTICLE_INT_TYPE vx;
    REB_PARTICLE_INT_TYPE vy;
    REB_PARTICLE_INT_TYPE vz;
};

double t = 0;
const double dt = 0.01;

const unsigned int N = 9;
struct reb_particle_int* p_int = NULL;

struct reb_particle {
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double ax;
    double ay;
    double az;
    double m;
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



static void to_int(){
    for(unsigned int i=0; i<N; i++){ 
        p_int[i].x = p[i].x/scale_pos; 
        p_int[i].y = p[i].y/scale_pos; 
        p_int[i].z = p[i].z/scale_pos; 
        p_int[i].vx = p[i].vx/scale_vel; 
        p_int[i].vy = p[i].vy/scale_vel; 
        p_int[i].vz = p[i].vz/scale_vel; 
    }
}
static void to_double(){
    for(unsigned int i=0; i<N; i++){ 
        p[i].x = ((double)p_int[i].x)*scale_pos; 
        p[i].y = ((double)p_int[i].y)*scale_pos; 
        p[i].z = ((double)p_int[i].z)*scale_pos; 
        p[i].vx = ((double)p_int[i].vx)*scale_vel; 
        p[i].vy = ((double)p_int[i].vy)*scale_vel; 
        p[i].vz = ((double)p_int[i].vz)*scale_vel; 
    }
}

static void drift(){
    for(unsigned int i=0; i<N; i++){
        p_int[i].x += (REB_PARTICLE_INT_TYPE)(dt / 2. * (double)p_int[i].vx * scale_vel / scale_pos);
        p_int[i].y += (REB_PARTICLE_INT_TYPE)(dt / 2. * (double)p_int[i].vy * scale_vel / scale_pos);
        p_int[i].z += (REB_PARTICLE_INT_TYPE)(dt / 2. * (double)p_int[i].vz * scale_vel / scale_pos);
    }
}

static void kick(){
    for(unsigned int i=0; i<N; i++){
        p_int[i].vx += (REB_PARTICLE_INT_TYPE)(dt * p[i].ax / scale_vel);
        p_int[i].vy += (REB_PARTICLE_INT_TYPE)(dt * p[i].ay / scale_vel);
        p_int[i].vz += (REB_PARTICLE_INT_TYPE)(dt * p[i].az / scale_vel);
    }
}

static void gravity(){
    for(unsigned int i=0; i<N; i++){
        p[i].ax = 0.;
        p[i].ay = 0.;
        p[i].az = 0.;
        for(unsigned int j=0; j<N; j++){
            if (i!=j){
                const double dx = p[i].x - p[j].x;
                const double dy = p[i].y - p[j].y;
                const double dz = p[i].z - p[j].z;
                const double _r = sqrt(dx*dx + dy*dy + dz*dz);
                const double prefact = -1/(_r * _r * _r)*p[j].m;
                
                p[i].ax    += prefact*dx;
                p[i].ay    += prefact*dy;
                p[i].az    += prefact*dz;
            }
        }
    }
}

void janus_step(){
    // One leapfrog step
    drift();
    
    to_double(); 
    gravity();
    kick();

    drift();

    t += dt;
    
    // Only needed for floating point outputs 
    to_double(); 
}

int main(){
    int i = 0;
    ofstream myfile ("original.txt");
    myfile.precision(dbl::max_digits10);
    // Setup initial integer coordinates
    if (p_int==NULL){
        p_int = (struct reb_particle_int *) malloc(sizeof(struct reb_particle_int)*N);
        to_int(); 
    }
    
    while (t<2.*M_PI*1e1){ // 1 year
        janus_step();
        i = i + 1;
        myfile << "------------" << i << "th----------" << endl;
        for(int i=0;i<N;i++){
            myfile << p[i].vx << "," << p[i].vy << endl;
        }
    }

    return 1;
}