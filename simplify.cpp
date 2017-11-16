#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include <fstream>

using namespace std;

typedef numeric_limits< double > dbl;

// Scale for conversion between FP and INT
const double scale_vel  = 1e-16; 
const double scale_pos  = 1e-16;

// INT datatype used
#define REB_PARTICLE_INT_TYPE int64_t

struct reb_particle_int {
    REB_PARTICLE_INT_TYPE x;
    REB_PARTICLE_INT_TYPE y;
    REB_PARTICLE_INT_TYPE z;
    REB_PARTICLE_INT_TYPE vx;
    REB_PARTICLE_INT_TYPE vy;
    REB_PARTICLE_INT_TYPE vz;
    REB_PARTICLE_INT_TYPE ax;
    REB_PARTICLE_INT_TYPE ay;
    REB_PARTICLE_INT_TYPE az;
    REB_PARTICLE_INT_TYPE m;
};

struct reb_particle_int_drift_temp {
    REB_PARTICLE_INT_TYPE x;
    REB_PARTICLE_INT_TYPE y;
    REB_PARTICLE_INT_TYPE z;
};

double t = 0;
const double dt = 0.01;

const unsigned int N = 9;
struct reb_particle_int* p_int = NULL;
struct reb_particle_int_drift_temp* p_int_drift_temp = NULL;

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
        p_int[i].ax = p[i].ax/scale_vel;
        p_int[i].ay = p[i].ay/scale_vel;
        p_int[i].az = p[i].az/scale_vel;
        p_int[i].m = p[i].m/scale_vel;
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
        // p_int_drift_temp[i].x = p_int[i].x;
        // p_int_drift_temp[i].y = p_int[i].y;
        // p_int_drift_temp[i].z = p_int[i].z;
        p_int[i].x += (REB_PARTICLE_INT_TYPE) (dt/2. * (double)p_int[i].vx);
        p_int[i].y += (REB_PARTICLE_INT_TYPE) (dt/2. * (double)p_int[i].vy);
        p_int[i].z += (REB_PARTICLE_INT_TYPE) (dt/2. * (double)p_int[i].vz);
    }
}

static void kick(){
    ofstream myfile ("simplify.txt", ios_base::app);
    for(unsigned int i=0; i<N; i++){
        myfile << "+++++++++++" << i << "th+++++++++++" << endl;
        p_int[i].vx += (REB_PARTICLE_INT_TYPE)(dt * p_int[i].ax);
        p_int[i].vy += (REB_PARTICLE_INT_TYPE)(dt * p_int[i].ay);
        p_int[i].vz += (REB_PARTICLE_INT_TYPE)(dt * p_int[i].az);
        myfile << "vx:" << p_int[i].vx << " vy:" << p_int[i].vy << " vz:" << p_int[i].vz << endl;
    }
}

static void gravity(){
    ofstream myfile ("simplify.txt", ios_base::app);
    for(unsigned int i=0; i<N; i++){
        p_int[i].ax = 0.;
        p_int[i].ay = 0.;
        p_int[i].az = 0.;
        for(unsigned int j=0; j<N; j++){
            if (i!=j){
                const double dx = (p_int[i].x - p_int[j].x)*1e-16;
                const double dy = (p_int[i].y - p_int[j].y)*1e-16;
                const double dz = (p_int[i].z - p_int[j].z)*1e-16;
                const double _r = (double)(sqrt(dx*dx + dy*dy + dz*dz)*1e16);
                myfile << "r is:" << _r << endl;
                const double prefact = -1/(_r*_r*_r)*p_int[j].m;
                
                p_int[i].ax    += prefact*dx;
                p_int[i].ay    += prefact*dy;
                p_int[i].az    += prefact*dz;
            }
        }
    }
}

void janus_step(){
    // One leapfrog step
    drift();
    
    //to_double(); 
    gravity();
    kick();

    drift();

    t += dt;
    
    // Only needed for floating point outputs 
    //to_double(); 
}

int main(){
    int i = 0;
    ofstream myfile ("simplify.txt");
    // Setup initial integer coordinates
    if (p_int==NULL){
        p_int = (struct reb_particle_int *) malloc(sizeof(struct reb_particle_int)*N);
        to_int();
        // myfile << "to_int()" << endl;
        // for (int i = 0; i < N; i++){
        //     myfile << "x:" << p_int[i].x << " y:" << p_int[i].y << " z:" << p_int[i].z << " vx:" << p_int[i].vx << " vy:" << p_int[i].vy << " vz:" << p_int[i].vz << endl;
        //     myfile << "ax:" << p_int[i].ax << " ay:" << p_int[i].ay << " az:" << p_int[i].az << " m:" << p_int[i].m << endl;;
        // }
    }
    
    myfile.precision(dbl::max_digits10);
    while (t<2.*M_PI*1e1){ // 1 year
        i += 1;
        janus_step();
        myfile << "------------" << i << "th----------" << endl;
        for(int i=0;i<N;i++){
            myfile << p[i].x << "," << p[i].y << endl;
        }
    }

    //janus_step();

    to_double();

    myfile.close();
    ofstream yourfile ("simplify.txt", ios_base::app);
    yourfile.precision(dbl::max_digits10);
    for(int i=0;i<N;i++){
        yourfile << p[i].x << "," << p[i].y << endl;
    }

    yourfile.close();
    return 1;
}