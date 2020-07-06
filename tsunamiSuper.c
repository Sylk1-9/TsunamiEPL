
/*
 * Hauteurs : Vanneste Sylvain (61851100) et Walsdorff Antoine (73571000)
 * Date de creation : Avril - Mai 2012
 */

// ---------------------------------
//------------- INCLUDE -----------
// ---------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

#include <GL/glfw.h> // A DECOMENTER POUR SOUMISSION !!!!!
#define GLFW_INCLUDE_GLU
//#include </Users/Antoine/Documents/Unif2eme/Q4/MECA1120ElemFinis/glfw-2.7.7/include/GL/glfw.h>
//#include </Users/sylvainvanneste/Dropbox/Sylvain_Antoine_Lu/Antoine-Sylvain/Element_fini/glfw-2.7.7/include/GL/glfw.h> // A EFFACER pour soummision

//#include <GLUT/glut.h> // erreur lors de la soumission sur serveur
//#include <GLUT/glut.h>

#include <time.h>

#define Error(a)   femError(a,__LINE__,__FILE__)
#define Warning(a) femWarning(a,  __LINE__, __FILE__)

// ---------------------------------
//------------ STRUCTURES -----------
// ---------------------------------
typedef enum {FEM_TRIANGLE,FEM_QUAD,FEM_EDGE} femElementType;
typedef enum {FEM_FULL,FEM_BAND,FEM_ITER} femSolverType;

typedef struct {
    int *elem;
    double *X;
    double *Y;
    double *H; 
    int nElem;
    int nNode;
    int nLocalNode;
} femMesh;

typedef struct {
    int elem[2];
    int node[2];
} femEdge;

typedef struct {
    femMesh *mesh;
    femEdge *edges;
    int nEdge;
    int nBoundary;
} femEdges;

typedef struct {
    int n;
    const double *xsi;
    const double *eta;
    const double *weight;
} femIntegration;

typedef struct {
    femSolverType type;
    void *solver;
} femSolver;

typedef struct {
    double *B;
    double **A;
    int size;
    int band;
} femBandSystem;

typedef struct {
    int n;
    int order;
    void (*x2)(double *xsi, double *eta);
    void (*phi2)(double xsi, double eta, double *phi);
    void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
    void (*x1)(double *xsi);
    void (*phi1)(double xsi, double *phi);
    void (*dphi1dx)(double xsi, double *dphidxsi);
    
} femDiscrete;

typedef struct {
    double *x;
    double *y;
    double *h;
    double *J;
    double *f;
    double *cT;
    double **dphidx;
    double **dphidy;
} femVectorElem;

typedef struct {
    double *h;
    double *J;
    double *nx;
    double *ny;
    double *cT;
    double ***Map;
} femVectorEdge;

typedef struct {
    femMesh *mesh;
    femEdges *edges;
    femDiscrete *space;
    femDiscrete *spaceSmall;
    femIntegration *rule1d;
    femIntegration *rule2d;
    femSolver *solver;
    femVectorElem *vectorElem;
    femVectorEdge *vectorEdge;
    int *map;
    int orderType;
    int size;
    double *X;
    double *Y;
    double *E;
    double *U;
    double *V;
    double *FE;
    double *FU;
    double *FV;
    double gravity;
    double Omega;
    double gamma;
    double R;
    double timeStep;

} femTsunami;


// ---------------------------------
//---------- PRE-DECLARATION --------
// ---------------------------------


void femError(char *text, int line, char *file);
int femEdgesCompare(const void *edgeOne, const void *edgeTwo);
void femTsunamiComputeMap(femTsunami *myTsunami, int orderType);

femSolver*           femSolverBandCreate(int size,int band);
femBandSystem*       femBandSystemCreate(int size, int band);
void                 femBandSystemFree(femBandSystem* myBandSystem);
void                 femBandSystemInit(femBandSystem *myBand);
double*              femBandSystemEliminate(femBandSystem *myBand);
void                 femBandSystemConstrain(femBandSystem *myBand, int myNode, double myValue);
void                 femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc);
double               femBandSystemGet(femBandSystem* myBandSystem, int i, int j);

void femTsunamiAddIntegralsElements(femTsunami *myTsunami);
void femTsunamiAddIntegralsEdges(femTsunami *myTsunami);
void femTsunamiMultiplyInverseMatrix(femTsunami *myTsunami);

double tsunamiInitialConditionOkada(double x, double y);
void tsunamiWriteFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem, int nsub);
int tsunamiReadFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem);

femVectorElem *femVectorElemCreate(femTsunami *myTsunami);
void femVectorElemFree(femVectorElem *theVectorElem, femTsunami *myTsunami);
void femVecteurEdgeCompute(femTsunami *myTsunami);
femVectorEdge *femVectorEdgeCreate(femTsunami *myTsunami);
void femVectorEdgeFree(femVectorEdge *theVectorEdge, femTsunami *myTsunami);


// ---------------------------------
//------- DONNEES INTEGRATION -------
// ---------------------------------

static const double _gaussEdge2Xsi[2]     = { 0.577350269189626,-0.577350269189626 };
static const double _gaussEdge2Weight[2]  = { 1.000000000000000, 1.000000000000000 };

static const double _gaussEdge3Xsi[3]     = {-0.774596669241483, 0.000000000000000, 0.774596669241483};
static const double _gaussEdge3Weight[3]  = { 0.555555555555556, 0.888888888888889,  0.555555555555556};

static const double _gaussQuad4Xsi[4]     = {-0.577350269189626,-0.577350269189626, 0.577350269189626, 0.577350269189626 };
static const double _gaussQuad4Eta[4]     = { 0.577350269189626,-0.577350269189626,-0.577350269189626, 0.577350269189626 };
static const double _gaussQuad4Weight[4]  = { 1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000 };
static const double _gaussQuad9Xsi[9]     = { 0.774596669241483, 0.000000000000000,-0.774596669241483,
    0.774596669241483, 0.000000000000000,-0.774596669241483,
    0.774596669241483, 0.000000000000000,-0.774596669241483 };
static const double _gaussQuad9Eta[9]     = { 0.774596669241483, 0.774596669241483, 0.774596669241483,
    0.000000000000000, 0.000000000000000, 0.000000000000000,
    -0.774596669241483,-0.774596669241483,-0.774596669241483 };
static const double _gaussQuad9Weight[9]  = { 0.308641975308642, 0.493827160493827, 0.308641975308642,
    0.493827160493827, 0.790123456790123, 0.493827160493827,
    0.308641975308642, 0.493827160493827, 0.308641975308642 };

static const double _gaussTri3Xsi[3]      = { 0.166666666666667, 0.666666666666667, 0.166666666666667 };
static const double _gaussTri3Eta[3]      = { 0.166666666666667, 0.166666666666667, 0.666666666666667 };
static const double _gaussTri3Weight[3]   = { 0.166666666666667, 0.166666666666667, 0.166666666666667 };
static const double _gaussTri12Xsi[12]    = { 0.249286745170910, 0.249286745170910, 0.501426509658179,
    0.063089014491502, 0.063089014491502, 0.873821971016996,
    0.310352451033785, 0.636502499121399, 0.053145049844816,
    0.310352451033785, 0.636502499121399, 0.053145049844816 };
static const double _gaussTri12Eta[12]    = { 0.249286745170910, 0.501426509658179, 0.249286745170910,
    0.063089014491502, 0.873821971016996, 0.063089014491502,
    0.636502499121399, 0.053145049844816, 0.310352451033785,
    0.053145049844816, 0.310352451033785, 0.636502499121399 };
static const double _gaussTri12Weight[12] = { 0.058393137863189, 0.058393137863189, 0.058393137863189,
    0.025422453185104, 0.025422453185104, 0.025422453185104,
    0.041425537809187, 0.041425537809187, 0.041425537809187,
    0.041425537809187, 0.041425537809187, 0.041425537809187 };


void _1c0_phi(double xsi, double *phi)
{
    phi[0] = (1.0 - xsi)/2.0;
    phi[1] = (1.0 + xsi)/2.0;
}

void _2c0_phi(double xsi, double *phi)
{
    phi[0] = (xsi - 1.0)*xsi/2.0;
    phi[1] = (xsi + 1.0)*xsi/2.0;
    phi[2] = (1.0 - xsi)*(1.0 + xsi);
}

void _3c0_phi(double xsi, double *phi)
{
    phi[0] =   9./16 * (-1./3 - xsi) * ( 1./3 - xsi) * (1.   - xsi);
    phi[1] = -27./16 * (-1.   - xsi) * ( 1./3 - xsi) * (1.   - xsi);
    phi[2] =  27./16 * (-1.   - xsi) * (-1./3 - xsi) * (1.   - xsi);
    phi[3] =  -9./16 * (-1.   - xsi) * (-1./3 - xsi) * (1./3 - xsi);
}

void _3c0_dphidx(double xsi, double *dphidxsi)
{
    dphidxsi[0] =   9./16 * ( - ( 1./3 - xsi) * (1.   - xsi) - (-1./3 - xsi) * (1.   - xsi) - (-1./3 - xsi) * ( 1./3 - xsi) ) ;
    dphidxsi[1] = -27./16 * ( - ( 1./3 - xsi) * (1.   - xsi) - (-1.   - xsi) * (1.   - xsi) - (-1.   - xsi) * ( 1./3 - xsi) ) ;
    dphidxsi[2] =  27./16 * ( - (-1./3 - xsi) * (1.   - xsi) - (-1.   - xsi) * (1.   - xsi) - (-1.   - xsi) * (-1./3 - xsi) ) ;
    dphidxsi[3] =  -9./16 * ( - (-1./3 - xsi) * (1./3 - xsi) - (-1.   - xsi) * (1./3 - xsi) - (-1.   - xsi) * (-1./3 - xsi) ) ;
}

void _q1c0_x(double *xsi, double *eta)
{
    xsi[0] =  1.0;  eta[0] =  1.0;
    xsi[1] = -1.0;  eta[1] =  1.0;
    xsi[2] = -1.0;  eta[2] = -1.0;
    xsi[3] =  1.0;  eta[3] = -1.0;
}

void _q1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = (1.0 + xsi) * (1.0 + eta) / 4.0;
    phi[1] = (1.0 - xsi) * (1.0 + eta) / 4.0;
    phi[2] = (1.0 - xsi) * (1.0 - eta) / 4.0;
    phi[3] = (1.0 + xsi) * (1.0 - eta) / 4.0;
}

void _q1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] =   (1.0 + eta) / 4.0;
    dphidxsi[1] = - (1.0 + eta) / 4.0;
    dphidxsi[2] = - (1.0 - eta) / 4.0;
    dphidxsi[3] =   (1.0 - eta) / 4.0;
    dphideta[0] =   (1.0 + xsi) / 4.0;
    dphideta[1] =   (1.0 - xsi) / 4.0;
    dphideta[2] = - (1.0 - xsi) / 4.0;
    dphideta[3] = - (1.0 + xsi) / 4.0;
}

void _q2c0_x(double *xsi, double *eta)
{
    xsi[0] =  1.0;  eta[0] =  1.0;
    xsi[1] = -1.0;  eta[1] =  1.0;
    xsi[2] = -1.0;  eta[2] = -1.0;
    xsi[3] =  1.0;  eta[3] = -1.0;
    xsi[4] =  0.0;  eta[4] =  1.0;
    xsi[5] = -1.0;  eta[5] =  0.0;
    xsi[6] =  0.0;  eta[6] = -1.0;
    xsi[7] =  1.0;  eta[7] =  0.0;
    xsi[8] =  0.0;  eta[8] =  0.0;
}

void _q2c0_phi(double xsi, double eta, double *phi)
{
    phi[0] =  xsi*(1.0+xsi)*eta*(1.0+eta)/4.0;
    phi[1] = -xsi*(1.0-xsi)*eta*(1.0+eta)/4.0;
    phi[2] =  xsi*(1.0-xsi)*eta*(1.0-eta)/4.0;
    phi[3] = -xsi*(1.0+xsi)*eta*(1.0-eta)/4.0;
    phi[4] =  (1.0-xsi*xsi)*eta*(1.0+eta)/2.0;
    phi[5] = -xsi*(1.0-xsi)*(1.0-eta*eta)/2.0;
    phi[6] = -(1.0-xsi*xsi)*eta*(1.0-eta)/2.0;
    phi[7] =  xsi*(1.0+xsi)*(1.0-eta*eta)/2.0;
    phi[8] =  (1.0-xsi*xsi)*(1.0-eta*eta);
}

void _q2c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] =  (1.0+2.0*xsi)*eta*(1.0+eta)/4.0;
    dphidxsi[1] = (-1.0+2.0*xsi)*eta*(1.0+eta)/4.0;
    dphidxsi[2] =  (1.0-2.0*xsi)*eta*(1.0-eta)/4.0;
    dphidxsi[3] = -(1.0+2.0*xsi)*eta*(1.0-eta)/4.0;
    dphidxsi[4] =       -2.0*xsi*eta*(1.0+eta)/2.0;
    dphidxsi[5] = (-1.0+2.0*xsi)*(1.0-eta*eta)/2.0;
    dphidxsi[6] =        2.0*xsi*eta*(1.0-eta)/2.0;
    dphidxsi[7] =  (1.0+2.0*xsi)*(1.0-eta*eta)/2.0;
    dphidxsi[8] =       -2.0*xsi*(1.0-eta*eta);
    dphideta[0] =  xsi*(1.0+xsi)*(1.0+2.0*eta)/4.0;
    dphideta[1] = -xsi*(1.0-xsi)*(1.0+2.0*eta)/4.0;
    dphideta[2] =  xsi*(1.0-xsi)*(1.0-2.0*eta)/4.0;
    dphideta[3] = -xsi*(1.0+xsi)*(1.0-2.0*eta)/4.0;
    dphideta[4] =  (1.0-xsi*xsi)*(1.0+2.0*eta)/2.0;
    dphideta[5] =  xsi*(1.0-xsi)*2.0*eta/2.0;
    dphideta[6] = -(1.0-xsi*xsi)*(1.0-2.0*eta)/2.0;
    dphideta[7] =  xsi*(1.0+xsi)*(-2.0*eta)/2.0;
    dphideta[8] =  (1.0-xsi*xsi)*(-2.0*eta);
}

void _p1c0_x(double *xsi, double *eta)
{
    xsi[0] =  0.0;  eta[0] =  0.0;
    xsi[1] =  1.0;  eta[1] =  0.0;
    xsi[2] =  0.0;  eta[2] =  1.0;
}

void _p1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = 1 - xsi - eta;
    phi[1] = xsi;
    phi[2] = eta;
}

void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] = -1.0;
    dphidxsi[1] =  1.0;
    dphidxsi[2] =  0.0;
    dphideta[0] = -1.0;
    dphideta[1] =  0.0;
    dphideta[2] =  1.0;
}

void _p2c0_x(double *xsi, double *eta)
{
    xsi[0] =  0.0;  eta[0] =  0.0;
    xsi[1] =  1.0;  eta[1] =  0.0;
    xsi[2] =  0.0;  eta[2] =  1.0;
    xsi[3] =  0.5;  eta[3] =  0.0;
    xsi[4] =  0.5;  eta[4] =  0.5;
    xsi[5] =  0.0;  eta[5] =  0.5;
}

void _p2c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = 1.0 - 3.0*(xsi+eta) + 2.0*(xsi+eta)*(xsi+eta);
    phi[1] = xsi*(2.0*xsi-1.0);
    phi[2] = eta*(2.0*eta-1.0);
    phi[3] = 4.0*xsi*(1.0-xsi-eta);
    phi[4] = 4.0*xsi*eta;
    phi[5] = 4.0*eta*(1.0-xsi-eta);
}

void _p2c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] = - 3.0 + 4.0*xsi + 4.0*eta;
    dphidxsi[1] = - 1.0 + 4.0*xsi          ;
    dphidxsi[2] =   0.0                    ;
    dphidxsi[3] =   4.0 - 8.0*xsi - 4.0*eta;
    dphidxsi[4] =                   4.0*eta;
    dphidxsi[5] =                 - 4.0*eta;
    dphideta[0] = - 3.0 + 4.0*xsi + 4.0*eta;
    dphideta[1] =   0.0                    ;
    dphideta[2] =  -1.0           + 4.0*eta;
    dphideta[3] =       - 4.0*xsi          ;
    dphideta[4] =         4.0*xsi          ;
    dphideta[5] =   4.0 - 4.0*xsi - 8.0*eta;
}

// ---------------------------------
//------------- SOLVEUR -------------
// ---------------------------------

femSolver *femSolverBandCreate(int size, int band)
{
    femSolver *mySolver = malloc(sizeof(femSolver));
    mySolver->type = FEM_BAND;
    mySolver->solver = (femSolver *)femBandSystemCreate(size,band);
    return(mySolver);
}

femBandSystem *femBandSystemCreate(int size, int band)
{
    femBandSystem *myBandSystem = malloc(sizeof(femBandSystem));
    myBandSystem->B = malloc(sizeof(double)*size*(band+1));
    myBandSystem->A = malloc(sizeof(double*)*size);
    myBandSystem->size = size;
    myBandSystem->band = band;
    myBandSystem->A[0] = myBandSystem->B + size;
    int i;
    for (i=1 ; i < size ; i++)
        myBandSystem->A[i] = myBandSystem->A[i-1] + band - 1;
    femBandSystemInit(myBandSystem);
    return(myBandSystem);
}

void femBandSystemFree(femBandSystem *myBandSystem)
{
    free(myBandSystem->B);
    free(myBandSystem->A);
    free(myBandSystem);
}

void femBandSystemInit(femBandSystem *myBandSystem)
{
    int i;
    int size = myBandSystem->size;
    int band = myBandSystem->band;
    for (i=0 ; i < size*(band+1) ; i++)
        myBandSystem->B[i] = 0;
}

double femBandSystemGet(femBandSystem* myBandSystem, int myRow, int myCol)
{
    double value = 0;
    if (myCol >= myRow && myCol < myRow+myBandSystem->band)  value = myBandSystem->A[myRow][myCol];
    return(value);
}

void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    int i,j;
    for (i = 0; i < nLoc; i++) {
        int myRow = map[i];
        for(j = 0; j < nLoc; j++) {
            int myCol = map[j];
            if (myCol >= myRow)  myBandSystem->A[myRow][myCol] += Aloc[i*nLoc+j]; }
        myBandSystem->B[myRow] += Bloc[i]; }
}

void femBandSystemConstrain(femBandSystem *myBand, int myNode, double myValue)
{
    double  **A, *B;
    int     i, size, band, ifirst, iend;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;
    
    ifirst = fmax(0,myNode - band + 1);
    iend   = myNode;
    for (i=ifirst; i < iend; i++) {
        B[i] -= myValue * A[i][myNode];
        A[i][myNode] = 0; }
    
    ifirst = myNode+1;
    iend = fmin(myNode + band,size);
    for (i=ifirst; i < iend; i++) {
        B[i] -= myValue * A[myNode][i];
        A[myNode][i] = 0; }
    
    A[myNode][myNode] = 1;
    B[myNode] = myValue;
}

double  *femBandSystemEliminate(femBandSystem *myBand)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;
    
    /* Incomplete Cholesky factorization */
    
    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-4 ) {
            Error("Cannot eleminate with such a pivot"); }
        jend = fmin(k + band,size);
        for (i = k+1 ; i <  jend; i++) {
            factor = A[k][i] / A[k][k];
            for (j = i ; j < jend; j++)
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
    
    /* Back-substitution */
    
    for (i = (size-1); i >= 0 ; i--) {
        factor = 0;
        jend = fmin(i + band,size);
        for (j = i+1 ; j < jend; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }
    
    return(myBand->B);
}


// ---------------------------------
//--------------- FEM ---------------
// ---------------------------------

femIntegration *femIntegrationCreate(int n, femElementType type)
{
    femIntegration *theRule = malloc(sizeof(femIntegration));
    if (type == FEM_QUAD && n == 4) {
        theRule->n      = 4;
        theRule->xsi    = _gaussQuad4Xsi;
        theRule->eta    = _gaussQuad4Eta;
        theRule->weight = _gaussQuad4Weight; }
    else if (type == FEM_QUAD && n == 9) {
        theRule->n      = 9;
        theRule->xsi    = _gaussQuad9Xsi;
        theRule->eta    = _gaussQuad9Eta;
        theRule->weight = _gaussQuad9Weight; }
    else if (type == FEM_TRIANGLE && n == 3) {
        theRule->n      = 3;
        theRule->xsi    = _gaussTri3Xsi;
        theRule->eta    = _gaussTri3Eta;
        theRule->weight = _gaussTri3Weight; }
    else if (type == FEM_TRIANGLE && n == 12) {
        theRule->n      = 12;
        theRule->xsi    = _gaussTri12Xsi;
        theRule->eta    = _gaussTri12Eta;
        theRule->weight = _gaussTri12Weight; }
    else if (type == FEM_EDGE && n == 2) {
        theRule->n      = 2;
        theRule->xsi    = _gaussEdge2Xsi;
        theRule->eta    = NULL;
        theRule->weight = _gaussEdge2Weight; }
    else if (type == FEM_EDGE && n == 3) {
        theRule->n      = 3;
        theRule->xsi    = _gaussEdge3Xsi;
        theRule->eta    = NULL;
        theRule->weight = _gaussEdge3Weight; }
    
    else Error("Cannot create such an integration rule !");
    return theRule;
}

void femIntegrationFree(femIntegration *theRule)
{
    free(theRule);
}

femDiscrete *femDiscreteCreate(int n, femElementType type)
{
    femDiscrete *theSpace = malloc(sizeof(femDiscrete));
    if (type == FEM_QUAD && n == 4) {
        theSpace->n       = 4;
        theSpace->order   = 1;
        theSpace->x2      = _q1c0_x;
        theSpace->phi2    = _q1c0_phi;
        theSpace->dphi2dx = _q1c0_dphidx;
        theSpace->phi1    = _1c0_phi;}
    else if (type == FEM_QUAD && n == 9) {
        theSpace->n       = 9;
        theSpace->order   = 2;
        theSpace->x2      = _q2c0_x;
        theSpace->phi2    = _q2c0_phi;
        theSpace->dphi2dx = _q2c0_dphidx;
        theSpace->phi1    = _2c0_phi;}
    else if (type == FEM_TRIANGLE && n == 3) {
        theSpace->n       = 3;
        theSpace->order   = 1;
        theSpace->x2      = _p1c0_x;
        theSpace->phi2    = _p1c0_phi;
        theSpace->dphi2dx = _p1c0_dphidx;
        theSpace->phi1    = _1c0_phi; }
    else if (type == FEM_TRIANGLE && n == 6) {
        theSpace->n       = 6;
        theSpace->order   = 2;
        theSpace->x2      = _p2c0_x;
        theSpace->phi2    = _p2c0_phi;
        theSpace->dphi2dx = _p2c0_dphidx;
        theSpace->phi1    = _2c0_phi;}
    else Error("Cannot create such a discrete space !");
    return theSpace;
}

void femDiscreteFree(femDiscrete *theSpace)
{
    free(theSpace);
}

void femTsunamiComputeMap(femTsunami *myTsunami, int orderType)
{
    femMesh  *theMesh = myTsunami->mesh;
    femEdges *theEdges = myTsunami->edges;
    
    int i,j,k, iElem;
    int nElem = theMesh->nElem;
    int nNode = theMesh->nNode;
    int nEdge = theEdges->nEdge;
    int nLocal = myTsunami->space->n;
    int nLocalNode = myTsunami->mesh->nLocalNode;
    double xsi[16],eta[16]; //,phi[4];
    myTsunami->space->x2(xsi,eta);
    femDiscrete *theLinearSpace;
    if (nLocalNode == 3) theLinearSpace = femDiscreteCreate(3,FEM_TRIANGLE);
    else        theLinearSpace = femDiscreteCreate(4,FEM_QUAD);
    
    if (orderType >= 1) {
        for (i=0; i < theMesh->nNode; i++) {
            myTsunami->X[i] = theMesh->X[i];
            myTsunami->Y[i] = theMesh->Y[i]; }
        for (i=0; i < nElem; i++)
            for (j=0; j < nLocalNode; j++)
                myTsunami->map[i*nLocal+j] = myTsunami->mesh->elem[i*nLocalNode+j]; }
    
    if (orderType == 2) {
        for (i=0; i < theEdges->nEdge; i++) {
            int node0 = theEdges->edges[i].node[0];
            int node1 = theEdges->edges[i].node[1];
            myTsunami->X[i+nNode] = (theMesh->X[theEdges->edges[i].node[0]]
                                      + theMesh->X[theEdges->edges[i].node[1]])/2.0;
            myTsunami->Y[i+nNode] = (theMesh->Y[theEdges->edges[i].node[0]]
                                      + theMesh->Y[theEdges->edges[i].node[1]])/2.0;
            
            for (j=0; j < 2; j++) {
                iElem = theEdges->edges[i].elem[j];
                if (iElem >= 0) {
                    int *elem = &myTsunami->mesh->elem[iElem*nLocalNode];
                    for (k=0; k < nLocalNode; k++) {
                        int k1 = k+1;
                        if (k1 == nLocalNode) k1 = 0;
                        if (node0 == elem[k] && node1 == elem[k1]) myTsunami->map[iElem*nLocal+nLocalNode+k] = i+nNode;
                        if (node1 == elem[k] && node0 == elem[k1]) myTsunami->map[iElem*nLocal+nLocalNode+k] = i+nNode; }}}}}
    
    if (orderType == 2 && nLocalNode ==  4) {
        for (i=0; i < nElem; i++) {
            myTsunami->map[i*nLocal+8] = i+nNode+nEdge;
            for (j=0; j < nLocalNode; j++) {
                myTsunami->X[i+nNode+nEdge] = 0.0;
                myTsunami->Y[i+nNode+nEdge] = 0.0; }
            for (j=0; j < nLocalNode; j++) {
                myTsunami->X[i+nNode+nEdge] += theMesh->X[theMesh->elem[i*nLocalNode+j]]/4.0;
                myTsunami->Y[i+nNode+nEdge] += theMesh->Y[theMesh->elem[i*nLocalNode+j]]/4.0; }}}
    
    femDiscreteFree(theLinearSpace);
}

femMesh *femMeshRead(femTsunami *myTsunami, const char *meshFileName, int orderType)
{
    femMesh *theMesh = malloc(sizeof(femMesh));
    myTsunami->orderType = orderType;

    int i,trash, *elem;
    int nodesOrderTriangle[] = {1,3,6,10};
    int integOrderTriangle[] = {1,3,12,12};
    int nodesOrderQuad[]     = {1,4,9,16};
    int integOrderQuad[]     = {1,4,9,9};
    int nElem, nNode, nEdge;
    FILE* file = fopen(meshFileName,"r");
    if (file == NULL) Error("No mesh file !");
    
    fscanf(file, "Number of nodes %d \n", &theMesh->nNode);
    nNode = theMesh->nNode;
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    theMesh->H = malloc(sizeof(double)*theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        fscanf(file,"%d : %le %le %le \n",&trash,&theMesh->X[i],&theMesh->Y[i], &theMesh->H[i]);
    }

    char str[256]; fgets(str, sizeof(str), file);
    if (!strncmp(str,"Number of triangles",19))  {
        sscanf(str,"Number of triangles %d \n", &theMesh->nElem);
        theMesh->elem = malloc(sizeof(int)*3*theMesh->nElem);
        theMesh->nLocalNode = 3;
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*3]);
            fscanf(file,"%d : %d %d %d\n", &trash,&elem[0],&elem[1],&elem[2]); }}
    
    else if (!strncmp(str,"Number of quads",15))  {
        sscanf(str,"Number of quads %d \n", &theMesh->nElem);
        theMesh->elem = malloc(sizeof(int)*4*theMesh->nElem);
        theMesh->nLocalNode = 4;
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*4]);
            fscanf(file,"%d : %d %d %d %d\n", &trash,&elem[0],&elem[1],&elem[2],&elem[3]); }}
    
    nElem = theMesh->nElem;
    nEdge = theMesh->nElem * theMesh->nLocalNode;
    
    if (theMesh->nLocalNode == 4) {
        switch (orderType) {
            case 0 : myTsunami->size = nElem; break;
            case 1 : myTsunami->size = nNode; break;
            case 2 : myTsunami->size = nNode + nEdge + nElem; break;
            case 3 : myTsunami->size = nNode + 2*nEdge + 4*nElem; break;  }
        myTsunami->space      = femDiscreteCreate(nodesOrderQuad[orderType],FEM_QUAD);
        myTsunami->spaceSmall = femDiscreteCreate(nodesOrderQuad[1],FEM_QUAD);
        myTsunami->rule2d     = femIntegrationCreate(integOrderQuad[orderType],FEM_QUAD);
    } 
    else if (theMesh->nLocalNode == 3) {
        switch (orderType) {
            case 0 : myTsunami->size = nElem; break;
            case 1 : myTsunami->size = nNode; break;
            case 2 : myTsunami->size = nNode + nEdge; break;
            case 3 : myTsunami->size = nNode + 2*nEdge + nElem; break;  }
        myTsunami->space      = femDiscreteCreate(nodesOrderTriangle[orderType],FEM_TRIANGLE);
        myTsunami->spaceSmall = femDiscreteCreate(nodesOrderTriangle[1],FEM_TRIANGLE);
        myTsunami->rule2d     = femIntegrationCreate(integOrderTriangle[orderType],FEM_TRIANGLE);
    }
    
    myTsunami->map     = malloc((myTsunami->space->n)*sizeof(int)*nElem);
    myTsunami->X       = malloc(sizeof(double)*myTsunami->size);
    myTsunami->Y       = malloc(sizeof(double)*myTsunami->size);
   
    fclose(file);
    return theMesh;
}

void femMeshFree(femMesh *theMesh)
{
    free(theMesh->X);
    free(theMesh->Y);
    free(theMesh->H);
    free(theMesh->elem);
    free(theMesh);
}

femEdges *femEdgesCreate(femMesh *theMesh)
{
    femEdges *theEdges = malloc(sizeof(femEdges));
    int nLoc = theMesh->nLocalNode;
    int i,j,n = theMesh->nElem * nLoc;
    femEdge* edges = malloc(n * sizeof(femEdge));
    theEdges->mesh  = theMesh;
    theEdges->edges = edges;
    theEdges->nEdge = n;
    theEdges->nBoundary = n;
    
    for (i = 0; i < theMesh->nElem; i++) {
        int *elem = &(theMesh->elem[i*nLoc]);
        for (j = 0; j < nLoc; j++) {
            int id = i * nLoc + j;
            edges[id].elem[0] = i;
            edges[id].elem[1] = -1;
            edges[id].node[0] = elem[j];
            edges[id].node[1] = elem[(j + 1) % nLoc]; }}
    
    qsort(theEdges->edges, theEdges->nEdge, sizeof(femEdge), femEdgesCompare);
    
    int index = 0;
    int nBoundary = 0;
    
    for (i=0; i < theEdges->nEdge; i++) {
        if (i == theEdges->nEdge - 1 || femEdgesCompare(&edges[i],&edges[i+1]) != 0) {
            edges[index] = edges[i];
            nBoundary++; }
        else {  edges[index] = edges[i];
            edges[index].elem[1] = edges[i+1].elem[0];
            i = i+1;}
        index++; }
    
    theEdges->edges = realloc(edges, index * sizeof(femEdge));
    theEdges->nEdge = index;
    theEdges->nBoundary = nBoundary;
    return theEdges;
}

void femEdgesFree(femEdges *theEdges)
{
    free(theEdges->edges);
    free(theEdges);
}

int femEdgesCompare(const void *edgeOne, const void *edgeTwo)
{
    int *nodeOne = ((femEdge*) edgeOne)->node;
    int *nodeTwo = ((femEdge*) edgeTwo)->node;
    int  diffMin = fmin(nodeOne[0],nodeOne[1]) - fmin(nodeTwo[0],nodeTwo[1]);
    int  diffMax = fmax(nodeOne[0],nodeOne[1]) - fmax(nodeTwo[0],nodeTwo[1]);
    
    if (diffMin < 0)    return  1;
    if (diffMin > 0)    return -1;
    if (diffMax < 0)    return  1;
    if (diffMax > 0)    return -1;
    return  0;
}

femTsunami *femTsunamiCreate(const char *meshFileName, int orderType)
{
    femTsunami *myTsunami = malloc(sizeof(femTsunami));
    myTsunami->orderType = orderType;
    myTsunami->Omega     = 2*M_PI/86400.0;
    myTsunami->gravity   = 9.81;
    myTsunami->gamma     = 0.0000001;
    myTsunami->R         = 6371220;
    myTsunami->mesh      = femMeshRead(myTsunami, meshFileName, orderType);
    myTsunami->edges     = femEdgesCreate(myTsunami->mesh);
    myTsunami->rule1d    = femIntegrationCreate(orderType+1,FEM_EDGE);

    femTsunamiComputeMap(myTsunami, orderType);

    int sizeLoc = myTsunami->space->n;
    int sizeGlo = myTsunami->mesh->nElem * sizeLoc + 1;
    myTsunami->solver = femSolverBandCreate(3*sizeLoc, sizeLoc);
    
    myTsunami->size = sizeGlo;
    myTsunami->E  = malloc(sizeof(double)*sizeGlo);
    myTsunami->U  = malloc(sizeof(double)*sizeGlo);
    myTsunami->V  = malloc(sizeof(double)*sizeGlo);
    myTsunami->FE = malloc(sizeof(double)*sizeGlo);
    myTsunami->FU = malloc(sizeof(double)*sizeGlo);
    myTsunami->FV = malloc(sizeof(double)*sizeGlo);
    
    myTsunami->vectorElem = femVectorElemCreate(myTsunami);
    myTsunami->vectorEdge = femVectorEdgeCreate(myTsunami);

    int i;
    for (i=0; i < sizeGlo; i++) {
        myTsunami->FE[i] = 0.0;
        myTsunami->FU[i] = 0.0;
        myTsunami->FV[i] = 0.0;
        myTsunami->E[i] = 0.0;
        myTsunami->U[i] = 0.0;
        myTsunami->V[i] = 0.0;  }
    
    return myTsunami;
}

void femTsunamiCompute(femTsunami *myTsunami)
{
    // RK : 
    int  size = myTsunami->size, i, j;
    double dt = myTsunami->timeStep;
    int sizeLoc = myTsunami->space->n;
    int sizeGlo = myTsunami->mesh->nElem * sizeLoc + 1;
    
    double* Eold = malloc(sizeof(double)*sizeGlo);
    double* Uold = malloc(sizeof(double)*sizeGlo);
    double* Vold = malloc(sizeof(double)*sizeGlo);
    double* Enew = malloc(sizeof(double)*sizeGlo);
    double* Unew = malloc(sizeof(double)*sizeGlo);
    double* Vnew = malloc(sizeof(double)*sizeGlo);

    for (i=0; i < size; i++) {
        Eold[i] = myTsunami->E[i];
        Enew[i] = myTsunami->E[i];
        Uold[i] = myTsunami->U[i];
        Unew[i] = myTsunami->U[i];
        Vold[i] = myTsunami->V[i];
        Vnew[i] = myTsunami->V[i];}
    
    //const int nStage   = 4;                               //  RK4
    const int nStage   = 2;                                 //  RK2
    //const double beta[4]  = {0.0,     0.5,   0.5, 1.0  }; //  RK4
    const double beta[4]  = {0.0, 1.0};                     //  RK2
    //const double gamma[4] = {1.0/6, 2.0/6, 2.0/6, 1.0/6}; //  RK4
    const double gamma[4] = {1.0/2.0, 1.0/2.0};
    
    for(j = 0; j < nStage; j++) {
        for (i=0; i < size; i++){
        	myTsunami->E[i] = Eold[i] + dt * beta[j] * myTsunami->FE[i];
            myTsunami->U[i] = Uold[i] + dt * beta[j] * myTsunami->FU[i];
            myTsunami->V[i] = Vold[i] + dt * beta[j] * myTsunami->FV[i];}

        femTsunamiAddIntegralsElements(myTsunami);
        femTsunamiAddIntegralsEdges(myTsunami);
        femTsunamiMultiplyInverseMatrix(myTsunami);
        
    	for (i=0; i < size; i++){
        	Enew[i] += dt * gamma[j] * myTsunami->FE[i];
            Unew[i] += dt * gamma[j] * myTsunami->FU[i];
            Vnew[i] += dt * gamma[j] * myTsunami->FV[i];}
    }
    
    for (i=0; i < size; i++) {
        myTsunami->E[i] = Enew[i];
        myTsunami->U[i] = Unew[i];
        myTsunami->V[i] = Vnew[i];}
    
    free(Enew);
    free(Unew);
    free(Vnew);
    free(Eold);
    free(Uold);
    free(Vold);

  
    /*
    
     // EE : 
    int  size = myTsunami->size, i;
    double *FE = myTsunami->FE;
    double *FU = myTsunami->FU;
    double *FV = myTsunami->FV;
    double *E = myTsunami->E;
    double *U = myTsunami->U;
    double *V = myTsunami->V;
    double theTimeStep = myTsunami->timeStep;
    
    
    for (i=0; i < size; i++) {
        FE[i] = 0.0;
        FU[i] = 0.0;
        FV[i] = 0.0; }
    femTsunamiAddIntegralsElements(myTsunami);
    femTsunamiAddIntegralsEdges(myTsunami);
    femTsunamiMultiplyInverseMatrix(myTsunami);
    for (i=0; i < size; i++) {
        E[i] += theTimeStep * FE[i];
        U[i] += theTimeStep * FU[i];
        V[i] += theTimeStep * FV[i]; }
     */
    
}

void femError(char *text, int line, char *file)
{
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);
}

void femTsunamiFree(femTsunami *myTsunami)
{
    free(myTsunami->FE);
    free(myTsunami->FU);
    free(myTsunami->FV);
    free(myTsunami->E);
    free(myTsunami->U);
    free(myTsunami->V);
    femVectorElemFree(myTsunami->vectorElem, myTsunami);
    femVectorEdgeFree(myTsunami->vectorEdge, myTsunami);
    femBandSystemFree(myTsunami->solver->solver);
    femIntegrationFree(myTsunami->rule1d);
    femIntegrationFree(myTsunami->rule2d);
    femDiscreteFree(myTsunami->space);
    femDiscreteFree(myTsunami->spaceSmall);
    femEdgesFree(myTsunami->edges);
    femMeshFree(myTsunami->mesh);
    free(myTsunami->map);
    free(myTsunami->X);
    free(myTsunami->Y);
    free(myTsunami);
}

// ---------------------------------
//------ FONCTION INTEGRALES----------
// ---------------------------------

void femInitialCondition(femTsunami *myTsunami)
{
    int iElem,k;
    for (iElem=0; iElem < myTsunami->mesh->nElem; iElem++)
    {
        for (k=0; k < myTsunami->space->n; k++){
            int ind = myTsunami->map[myTsunami->space->n*iElem+k];
            myTsunami->E[iElem*myTsunami->space->n+k] = tsunamiInitialConditionOkada(myTsunami->X[ind], myTsunami->Y[ind]);}
    }
}

void femVectorElemFree(femVectorElem *theVectorElem, femTsunami *myTsunami)
{
    free(theVectorElem->x);
    free(theVectorElem->y);
    free(theVectorElem->h);
    free(theVectorElem->f);
    free(theVectorElem->J);
    free(theVectorElem->cT);
    int nElem = myTsunami->mesh->nElem;
	int nbrInt = myTsunami->rule2d->n;
    int size = nbrInt*nElem;
    int i;
    for(i=0; i<size;i++){
        free(theVectorElem->dphidx[i]);
        free(theVectorElem->dphidy[i]);}
    free(theVectorElem->dphidx);
    free(theVectorElem->dphidy);
    free(theVectorElem);
}

femVectorElem *femVectorElemCreate(femTsunami *myTsunami)
{
    
    int nElem = myTsunami->mesh->nElem;
	int nbrInt = myTsunami->rule2d->n;
    int size = nbrInt*nElem;
    femVectorElem *theVectorElem = malloc(sizeof(femVectorElem));
    theVectorElem->h = malloc(sizeof(double)*size);
    theVectorElem->x = malloc(sizeof(double)*size);
    theVectorElem->y = malloc(sizeof(double)*size);
    theVectorElem->f = malloc(sizeof(double)*size);
    theVectorElem->J = malloc(sizeof(double)*size);
    theVectorElem->dphidx = malloc(sizeof(double)*size);
    theVectorElem->dphidy = malloc(sizeof(double)*size);
    theVectorElem->cT = malloc(sizeof(double)*size);
    int i;
    for(i=0; i<size;i++){
        theVectorElem->dphidx[i] = malloc(sizeof(double)*myTsunami->space->n);
        theVectorElem->dphidy[i] = malloc(sizeof(double)*myTsunami->space->n);
    }
    return theVectorElem;
}

void femVecteurElementsCompute(femTsunami *myTsunami)
{

    double  Omega = myTsunami->Omega;
    double  R     = myTsunami->R;
    int nElem = myTsunami->mesh->nElem;
	int nbrInt = myTsunami->rule2d->n;
    double Jack, dxdxsi, dxdeta, dydxsi, dydeta;
	double dxsidx, dxsidy, detadx, detady;
    
	double x,y, h;
    double dphidxsi[myTsunami->space->n];
    double dphideta[myTsunami->space->n];
    double dphidxsiSmall[myTsunami->spaceSmall->n];
    double dphidetaSmall[myTsunami->spaceSmall->n];
	double phi[myTsunami->space->n];
    double phiSmall[myTsunami->spaceSmall->n];
    
    const double *xsi = myTsunami->rule2d->xsi;
    const double *eta = myTsunami->rule2d->eta;
    
	int iElem, iInteg,m,k, ind;
    for(iElem=0;iElem<nElem; iElem++) {
        
		for(iInteg=0; iInteg<nbrInt; iInteg++) {
            double xsiLoc = xsi[iInteg], etaLoc = eta[iInteg];
            myTsunami->space->phi2(xsiLoc,etaLoc,phi);
            myTsunami->spaceSmall->phi2(xsiLoc,etaLoc,phiSmall);
            myTsunami->space->dphi2dx(xsiLoc, etaLoc, dphidxsi, dphideta);
            myTsunami->spaceSmall->dphi2dx(xsiLoc, etaLoc, dphidxsiSmall, dphidetaSmall);
            
			dxdxsi=0.0; dxdeta=0.0; dydxsi=0.0; dydeta=0.0;
			for(m=0; m<myTsunami->spaceSmall->n; m++){
				ind = myTsunami->mesh->elem[myTsunami->spaceSmall->n*iElem+m];
				dxdxsi += myTsunami->mesh->X[ind]*dphidxsiSmall[m];
				dxdeta += myTsunami->mesh->X[ind]*dphidetaSmall[m];
				dydxsi += myTsunami->mesh->Y[ind]*dphidxsiSmall[m];
				dydeta += myTsunami->mesh->Y[ind]*dphidetaSmall[m];}
            
			Jack=fabs(dxdxsi*dydeta-dxdeta*dydxsi);
            
            dxsidx= (1.0/Jack)*dydeta;
            dxsidy=-(1.0/Jack)*dxdeta;
            detadx=-(1.0/Jack)*dydxsi;
            detady= (1.0/Jack)*dxdxsi;
            
			for(m=0;m<myTsunami->space->n;m++){
				myTsunami->vectorElem->dphidx[iElem*nbrInt + iInteg][m]=dphidxsi[m]*dxsidx+dphideta[m]*detadx;
				myTsunami->vectorElem->dphidy[iElem*nbrInt + iInteg][m]=dphidxsi[m]*dxsidy+dphideta[m]*detady;}
            
            x=0.0; y=0.0; h=0.0;
            for (k=0 ; k<myTsunami->spaceSmall->n ; k++) {
                ind = myTsunami->mesh->elem[myTsunami->spaceSmall->n*iElem+k];
                x       += phiSmall[k]*myTsunami->mesh->X[ind];
				y       += phiSmall[k]*myTsunami->mesh->Y[ind];
                h       += phiSmall[k]*myTsunami->mesh->H[ind];}
            
            double z3d = R*(4.0*R*R - x*x - y*y) / (4.0*R*R + x*x + y*y);
            double lat = asin(z3d/R);
			double f = 2*Omega*sin(lat) ;
            double coefTerre = (4*R*R+ x*x + y*y  )/(4*R*R) ;

            
            myTsunami->vectorElem->x[iElem*nbrInt + iInteg] = x;
            myTsunami->vectorElem->y[iElem*nbrInt + iInteg] = y;
            myTsunami->vectorElem->h[iElem*nbrInt + iInteg] = h;
            myTsunami->vectorElem->f[iElem*nbrInt + iInteg] = f;
            myTsunami->vectorElem->J[iElem*nbrInt + iInteg] = Jack;
            myTsunami->vectorElem->cT[iElem*nbrInt + iInteg]= coefTerre;
        }
    }
    
    
}

void femTsunamiAddIntegralsElements(femTsunami *myTsunami)
{
    
    double *BE = myTsunami->FE;
    double *BU = myTsunami->FU;
    double *BV = myTsunami->FV;
    double  g     = myTsunami->gravity;
    double  R     = myTsunami->R;
    double  gamma = myTsunami->gamma;
    int nElem = myTsunami->mesh->nElem;
	int nbrInt = myTsunami->rule2d->n;    
	double u,v,x,y,f,coefTerre, etah,h, Jack;
    double dphidx[myTsunami->space->n];
    double dphidy[myTsunami->space->n];
	double phi[myTsunami->space->n];
    const double *w   = myTsunami->rule2d->weight;
    const double *xsi = myTsunami->rule2d->xsi;
    const double *eta = myTsunami->rule2d->eta;
    
	int iElem, iInteg,m,k, ind;
	for(iElem=0;iElem<nElem; iElem++) {
            
		for(iInteg=0; iInteg<nbrInt; iInteg++) {
            myTsunami->space->phi2(xsi[iInteg],eta[iInteg],phi);
            x    = myTsunami->vectorElem->x[iElem*nbrInt+iInteg];
            y    = myTsunami->vectorElem->y[iElem*nbrInt+iInteg];
            h    = myTsunami->vectorElem->h[iElem*nbrInt+iInteg];
            Jack = myTsunami->vectorElem->J[iElem*nbrInt+iInteg];
            f    = myTsunami->vectorElem->f[iElem*nbrInt+iInteg];
            coefTerre = myTsunami->vectorElem->cT[iElem*nbrInt+iInteg];
            
            for(m=0;m<myTsunami->space->n;m++){
                dphidx[m]= myTsunami->vectorElem->dphidx[iElem*nbrInt+iInteg][m];
                dphidy[m]= myTsunami->vectorElem->dphidy[iElem*nbrInt+iInteg][m];}

            u=0.0; v=0.0; etah=0.0;
            for (k=0 ; k<myTsunami->space->n ; k++) {
                ind = myTsunami->space->n*iElem+k;
				u       += phi[k]*myTsunami->U[ind];
				v       += phi[k]*myTsunami->V[ind];
                etah    += phi[k]*myTsunami->E[ind];}
            
            double JackW = Jack * w[iInteg];
			for (m=0 ; m<myTsunami->space->n ; m++){
                ind = myTsunami->space->n*iElem+m;
                BE[ind] += JackW * ((dphidx[m]*h*u + dphidy[m]*h*v) * coefTerre + phi[m]*h*(x*u + y*v)/(R*R));
                BU[ind] += JackW * (phi[m]*( f*v - gamma*u) + (dphidx[m] * g * etah) * coefTerre + phi[m]*g*x*etah/(2*R*R));
                BV[ind] += JackW * (phi[m]*(-f*u - gamma*v) + (dphidy[m] * g * etah) * coefTerre + phi[m]*g*y*etah/(2*R*R)); }
		}
	}
     
}

femVectorEdge *femVectorEdgeCreate(femTsunami *myTsunami)
{
    int nEdge = myTsunami->edges->nEdge;
	int nbrInt = myTsunami->rule1d->n;
    int size = nEdge*nbrInt;
    int n = myTsunami->orderType + 1;
    femVectorEdge *theVectorEdge = malloc(sizeof(femVectorEdge));
    theVectorEdge->h = malloc(sizeof(double)*size);
    theVectorEdge->cT = malloc(sizeof(double)*size);
    theVectorEdge->J = malloc(sizeof(double)*nEdge);
    theVectorEdge->nx = malloc(sizeof(double)*nEdge);
    theVectorEdge->ny = malloc(sizeof(double)*nEdge);
    theVectorEdge->Map = malloc(sizeof(double)*nEdge);
    
    int i, j;
    for(i=0;i<nEdge;i++) {
        theVectorEdge->Map[i] = malloc(sizeof(double)*2);
        for(j=0;j<2; j++) theVectorEdge->Map[i][j] = malloc(sizeof(double)*n);}
    return theVectorEdge;
}

void femVectorEdgeFree(femVectorEdge *theVectorEdge, femTsunami *myTsunami)
{
    int nEdge = myTsunami->edges->nEdge;
    free(theVectorEdge->h);
    free(theVectorEdge->J);
    free(theVectorEdge->nx);
    free(theVectorEdge->ny);
    free(theVectorEdge->cT);
    int i, j;
    for(i=0;i<nEdge;i++) {
        for(j=0;j<2; j++) free(theVectorEdge->Map[i][j]);
        free(theVectorEdge->Map[i]);}
    free(theVectorEdge->Map);
    free(theVectorEdge);   
}

void femVecteurEdgeCompute(femTsunami *myTsunami)
{
    int 	n = myTsunami->orderType + 1;
    double  xEdge[n],yEdge[n],phi[n], y, x, h;
    double  Jack;
    int     k,iEdge, iInteg, map[2][n];
    double phiSmall[2];
    double R = myTsunami->R;
    femIntegration *theRule = myTsunami->rule1d;
    femDiscrete *theSpace = myTsunami->space;
    int sizeGlo = myTsunami->mesh->nElem * theSpace->n + 1;
    for (iEdge=0; iEdge < myTsunami->edges->nEdge; iEdge++) {
        int nLocalNode = myTsunami->mesh->nLocalNode;
        int elem1 = myTsunami->edges->edges[iEdge].elem[0];
        int elem2 = myTsunami->edges->edges[iEdge].elem[1];
        int node1 = myTsunami->edges->edges[iEdge].node[0];
        int node2 = myTsunami->edges->edges[iEdge].node[1];
        int* nod = &(myTsunami->map[elem1*theSpace->n]);
        for(k=0; k<nLocalNode; k++) if (node1==nod[k]) break;
        if(theSpace->order==2) map[0][2] = elem1*theSpace->n + nLocalNode + k;
        map[0][0] = elem1*theSpace->n + k;
        map[0][1] = elem1*theSpace->n + (k + 1)%nLocalNode;
        if(elem2>=0)
        {
            int* nodes2 = &(myTsunami->map[elem2*theSpace->n]);
            for(k=0; k<nLocalNode; k++)
            {
                if (node2==nodes2[k]) break;
            }
            if(theSpace->order==2) map[1][2] = elem2*theSpace->n + nLocalNode + k;
            map[1][0] = elem2*theSpace->n + (k + 1)%nLocalNode;
            map[1][1] = elem2*theSpace->n + k;
        }
        else
        {
            if(theSpace->order==2) map[1][2] = sizeGlo-1;
            map[1][0] = sizeGlo-1;
            map[1][1] = sizeGlo-1;
        }
        for (k=0; k < 2; ++k){
        	int nod = myTsunami->edges->edges[iEdge].node[k];
        	xEdge[k] = myTsunami->X[nod];
        	yEdge[k] = myTsunami->Y[nod];}
        
        if(theSpace->order == 2){
            xEdge[2] = myTsunami->X[myTsunami->mesh->nNode+iEdge];
            yEdge[2] = myTsunami->Y[myTsunami->mesh->nNode+iEdge];}
        double dxdxsi = (xEdge[1] - xEdge[0]);
        double dydxsi = (yEdge[1] - yEdge[0]);
        double norm = sqrt(dxdxsi*dxdxsi + dydxsi*dydxsi);
        double nx =  dydxsi/norm;
        double ny = -dxdxsi/norm;
        Jack = norm / 2.0;
        for (iInteg=0; iInteg < theRule->n; iInteg++) {
            
            theSpace->phi1(theRule->xsi[iInteg],phi);
            myTsunami->spaceSmall->phi1(theRule->xsi[iInteg],phiSmall);
            
            x = 0.0; y = 0.0; h = 0.0;
            for (k=0 ; k<2 ; k++) {
                h += phiSmall[k]*myTsunami->mesh->H[myTsunami->edges->edges[iEdge].node[k]];
                x += phiSmall[k]*myTsunami->X[myTsunami->edges->edges[iEdge].node[k]];
                y += phiSmall[k]*myTsunami->Y[myTsunami->edges->edges[iEdge].node[k]];}
            
            double coefTerre = (4*R*R+ x*x + y*y  )/(4*R*R);
            myTsunami->vectorEdge->h[iEdge*theRule->n+iInteg] = h;
            myTsunami->vectorEdge->cT[iEdge*theRule->n+iInteg] = coefTerre;

        }
        myTsunami->vectorEdge->J[iEdge] = Jack;
        myTsunami->vectorEdge->ny[iEdge] = ny;
        myTsunami->vectorEdge->nx[iEdge] = nx;
        int j;
            for(j=0;j<2; j++)
                for(k=0;k<n;k++) myTsunami->vectorEdge->Map[iEdge][j][k] = map[j][k];

    }
    

}

void femTsunamiAddIntegralsEdges(femTsunami *myTsunami)
{
    
    int 	n = myTsunami->orderType + 1;
    double  phi[n], h, nx, ny, coefTerre;
    double  weight,Jack;
    double  eL,eR,uL,uR,vL,vR,unL,unR;
    double  qe,qu,qv;
    int     i,j,k,iEdge, iInteg, map[2][n];
    
    double *BE = myTsunami->FE;
    double *BU = myTsunami->FU;
    double *BV = myTsunami->FV;
    double *E = myTsunami->E;
    double *U = myTsunami->U;
    double *V = myTsunami->V;
    double g = myTsunami->gravity;
    femIntegration *theRule = myTsunami->rule1d;
    femDiscrete *theSpace = myTsunami->space;
    int sizeGlo = myTsunami->mesh->nElem * theSpace->n + 1;
    for (iEdge=0; iEdge < myTsunami->edges->nEdge; iEdge++) {
       
            for(j=0;j<2; j++)
                for(k=0;k<n;k++){ map[j][k] = myTsunami->vectorEdge->Map[iEdge][j][k];}
        nx = myTsunami->vectorEdge->nx[iEdge];
        ny = myTsunami->vectorEdge->ny[iEdge];
        Jack = myTsunami->vectorEdge->J[iEdge];
        
        int boundary = (map[1][0] == sizeGlo-1);
        for (iInteg=0; iInteg < theRule->n; iInteg++) {
            
            weight = theRule->weight[iInteg];
            theSpace->phi1(theRule->xsi[iInteg],phi);
            h = myTsunami->vectorEdge->h[iEdge*theRule->n+iInteg];
            coefTerre = myTsunami->vectorEdge->cT[iEdge*theRule->n+iInteg];
            
            uL = 0.0;vL = 0.0;uR = 0.0;vR = 0.0;eL = 0.0; eR = 0.0;
            for (k=0 ; k<n ; k++) {
                uL += phi[k]*U[map[0][k]];
                vL += phi[k]*V[map[0][k]];
                uR += phi[k]*U[map[1][k]];
                vR += phi[k]*V[map[1][k]];
                eL += phi[k]*E[map[0][k]];
                eR += phi[k]*E[map[1][k]];}
            
            eR = boundary ? eL : eR;
            unL = uL*nx+ vL*ny;
            unR = boundary ? -unL : uR*nx + vR*ny;
            qe =  0.5*h*   ( (unL+unR) + sqrt(g/h)*( eL-eR ) );
            qu =  0.5*g*nx*( ( eL+eR ) + sqrt(h/g)*(unL-unR) );
            qv =  0.5*g*ny*( ( eL+eR ) + sqrt(h/g)*(unL-unR) );

            for (i=0; i < n; i++) {
                BE[map[0][i]] -= qe*phi[i] * Jack * weight * coefTerre;
                BU[map[0][i]] -= qu*phi[i] * Jack * weight * coefTerre;
                BV[map[0][i]] -= qv*phi[i] * Jack * weight * coefTerre;
                BE[map[1][i]] += qe*phi[i] * Jack * weight * coefTerre;
                BU[map[1][i]] += qu*phi[i] * Jack * weight * coefTerre;
                BV[map[1][i]] += qv*phi[i] * Jack * weight * coefTerre;}}}

   }

void femTsunamiMultiplyInverseMatrix(femTsunami *myTsunami)
{
   
    double *BE = myTsunami->FE;
    double *BU = myTsunami->FU;
    double *BV = myTsunami->FV;
    femMesh *theMesh = myTsunami->mesh;
    femDiscrete *theSpace = myTsunami->space;
    femSolver *theSolver = myTsunami->solver;
    
    int n = theSpace->n;
    double Aloc[n*n];
    int iElem,i, j ,iInteg, mapElem[n],mapE[n],mapU[n],mapV[n];
    
	double Jack;
    const double *w   = myTsunami->rule2d->weight;
    const double *xsi = myTsunami->rule2d->xsi;
    const double *eta = myTsunami->rule2d->eta;
	double phi[myTsunami->space->n];
    
    for (i = 0; i < n; i++)   {
        mapE[i] = i;
        mapU[i] = i + n;
        mapV[i] = i + 2*n; }
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        femBandSystemInit((femBandSystem *)theSolver->solver);
        for (i = 0; i < n*n; i++)  Aloc[i] = 0;
        for (j=0; j < n; ++j) mapElem[j] = iElem*n + j;
        for(iInteg = 0; iInteg < myTsunami->rule2d->n; iInteg++)
        {
            myTsunami->space->phi2(xsi[iInteg],eta[iInteg],phi);            
            Jack = myTsunami->vectorElem->J[iElem*myTsunami->rule2d->n+iInteg];
            for (i=0; i < n; i++) {
                for (j=0; j < n; j++) {
                    Aloc[i*n + j] += phi[i]*phi[j]*w[iInteg]*Jack; }}
        }
        femBandSystemAssemble((femBandSystem *)theSolver->solver,Aloc,&BE[mapElem[0]],mapE,theSpace->n);
        femBandSystemAssemble((femBandSystem *)theSolver->solver,Aloc,&BU[mapElem[0]],mapU,theSpace->n);
        femBandSystemAssemble((femBandSystem *)theSolver->solver,Aloc,&BV[mapElem[0]],mapV,theSpace->n);
        double *soluce = femBandSystemEliminate((femBandSystem *)theSolver->solver);
        for(i = 0; i < theSpace->n; i++)
	    {
	    	BE[mapElem[i]] = soluce[i];
	    	BU[mapElem[i]] = soluce[i+n];
	    	BV[mapElem[i]] = soluce[i+2*n];
	    }
    }
}

// ---------------------------------
//------------- TSUNAMI -------------
// ---------------------------------

void tsunamiCompute(double dt, int nmax, int sub, int order, const char *meshFileName, const char *baseResultName)
{
    femTsunami *myTsunami = femTsunamiCreate(meshFileName, order);
    femInitialCondition(myTsunami);
    femVecteurElementsCompute(myTsunami);
    femVecteurEdgeCompute(myTsunami);
    
    double theDiscreteTime = 0.0;
    int    theIteration = 0;

    double theTimeStep = dt;
    myTsunami->timeStep = theTimeStep;
    double theStop = nmax*theTimeStep;
    clock_t tic = 0;
    int lastIteration;
    while (theStop >= theDiscreteTime) {
        theDiscreteTime += theTimeStep;

         if ((theIteration) % sub == 0 ) {
             tsunamiWriteFile(baseResultName, theIteration, myTsunami->U, myTsunami->V, myTsunami->E, myTsunami->mesh->nElem, myTsunami->space->n);
             clock_t toc = clock();
             printf("Time between the %d th and the %d th iteration : %f seconds\n", lastIteration ,theIteration , (double)(toc - tic) / CLOCKS_PER_SEC);
             tic = clock();
             lastIteration = theIteration;
             fflush(stdout);
         }
        
        femTsunamiCompute(myTsunami);
        theIteration += 1;
    }
    femTsunamiFree(myTsunami);

}

void tsunamiAnimate(double dt, int nmax, int sub, int order, const char *meshFileName, const char *baseResultName)
{
    
    int nElem,nNode,i,j,index,trash,*elem,xpos = 350,ypos = -36;
    double *X,*Y,*H,*E,*U,*V;
    int width,height;
    
    int nout = sub;
    double t;
    double R = 6371220;
    double BathMax = 9368;
    int pause = 0; int recule = 1; int avance = 0;
    //double zoom;
    GLfloat zoom = 2.f;

    
    char *forme = (char *) malloc(sizeof(char)*256);
    
    int cote = 0;
    int numberNodes = 0;
    double vagueFactor = 0.1;
    
    
    FILE* file = fopen(meshFileName,"r");
    if (file == NULL) {
        printf("Error : - cannot open mesh file :-) \n - \n");
        exit(0);
    }
    
    
    fscanf(file, "Number of nodes %d \n", &nNode);
    X = malloc(sizeof(double)*nNode);
    Y = malloc(sizeof(double)*nNode);
    H = malloc(sizeof(double)*nNode);
    
    for (i = 0; i < nNode; i++)
        fscanf(file,"%d : %le %le %le  \n",&trash,&X[i],&Y[i],&H[i]);
    
    fgets(forme, 256 , file);
    
    /// On detecte si on a affaire a des triangles
    if(strncmp(forme, "Number of triangles",19) == 0) {
        
        sscanf(forme, "Number of triangles %d \n", &nElem);
        cote = 3;
        switch (order) {
            case 0:
                numberNodes = 1;
                break;
            case 1:
                numberNodes = 3;
                break;
            case 2:
                numberNodes = 6;
                break;
        }
        
        elem = malloc(sizeof(int)*cote*nElem);
        U = malloc(sizeof(double)*numberNodes*nElem);
        V = malloc(sizeof(double)*numberNodes*nElem);
        E = malloc(sizeof(double)*numberNodes*nElem);
        
        for (i = 0; i < nElem; i++)
            fscanf(file,"%d : %d %d %d \n", &trash,&elem[i*cote],&elem[i*cote+1],&elem[i*cote+2]);
        fclose(file);
    }
    
    /// ou des carres
    else if(strncmp(forme, "Number of quads",15) == 0) {
        
        sscanf(forme, "Number of quads %d \n", &nElem);
        cote = 4;
        
        switch (order) {
            case 0:
                numberNodes = 1;
                break;
            case 1:
                numberNodes = 4;
                break;
            case 2:
                numberNodes = 9;
                break;
        }
        elem = malloc(sizeof(int)*cote*nElem);
        U = malloc(sizeof(double)*numberNodes*nElem);
        V = malloc(sizeof(double)*numberNodes*nElem);
        E = malloc(sizeof(double)*numberNodes*nElem);
        
        for (i = 0; i < nElem; i++)
            fscanf(file,"%d : %d %d %d %d \n", &trash,&elem[i*cote],&elem[i*cote+1],&elem[i*cote+2],&elem[i*cote+3]);
        fclose(file);
    }
    
    
    
    glfwInit();
    glfwOpenWindow(800,800,0,0,0,0,1,0,GLFW_WINDOW );
    glfwSetWindowTitle( "MECA1120 Tsunami 2013 - Sylvain V. et Antoine W." );
    
    
    glfwEnable( GLFW_STICKY_KEYS ); //add
    glfwSwapInterval( 1);
    GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 0.0 };
    GLfloat mat_shininess[] = { 50.0 };
    GLfloat light_position[] = { 8.0, 8.0, 8.0, 0.0 };
    
    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    GLfloat light_radiance[] = {1., 1., 1., 1.};
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_radiance);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_radiance);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_NORMALIZE);
    
    int xx = 345;
    int yy = -37;
    
    int coefXR = 0;
    int coefXL = 0;
    int coefYU = 0;
    int coefYD = 0;
    
    double t0 = 0;
    double frameInc = 0;
    int frame0,frame = -1;
    
    GLfloat colors[9], coord[3*cote], colorsS[9], coordS[3*cote];
    
    puts("\n \t Bonjour, bienvenu dans le guide de touches AZERTY (Soyez rapide pour appuyez sur les touche) :   \n \t\t- Pour tourner la terre, utilisez les touches directionelles haut/bas gauche/droite  \n \t\t- Pour agrandir/reduire la vague : W/Q \n \t\t- Pour zoomer/dezoomer : Z/A \n  \t\t- Pour mettre a simulation en pause/play : SPACE \n  \t\t- Pour reculer la simulation : R  \n  \t\t- Pour reeinitialiser la simulation : ENTER ");
    
    do {
        
        t = glfwGetTime() ;
        frame0 = frameInc;
        frameInc = (int) ((t-t0) * 4);
        
        if (frame0 != frameInc) {
            
            if(pause == 0){  if(avance == 0) frame +=1; if(recule == 0) frame -=1;}
            
            if (glfwGetKey(GLFW_KEY_ENTER) == GLFW_PRESS)
            {
                frame = 0; recule = 1; avance = 0; pause = 0; zoom = 2; vagueFactor = 0.1;
            }
            
            if (glfwGetKey(GLFW_KEY_SPACE) == GLFW_PRESS)
            {
                if (pause == 0)
                {
                    pause = 1;
                    avance = 1;
                    recule = 1;
                }
                else
                {
                    pause = 0;
                    avance = 0;
                    recule = 1;
                }
            }
            
            if (glfwGetKey('R') == GLFW_PRESS)
            {
                recule = 0;
                pause = 0;
                avance = 1;
            }
            
            if (glfwGetKey('Q') == GLFW_PRESS)
            {
                if(zoom>0.3)
                zoom -= 0.7;
            }
            if (glfwGetKey('W') == GLFW_PRESS)
            {
                if(zoom<5)
                zoom +=0.7;
            }
            
            //  glfwGetMousePos( &mouse, NULL );  //mouse = 389; // add
            glfwGetMousePos(&xpos, &ypos);
            
            char filename[256];
            char *chemin = (char *)malloc(sizeof(char)*256);
            strcpy(chemin, baseResultName);
            const char *basename = "-%08d.txt";
            strcat(chemin,  basename);
            sprintf(filename, chemin, frame * nout);
            if( (frame <= nmax/sub   &&  frame >=0 && pause == 0))
            { if( pause == 0 || ((recule == 0) )){
                if (access(filename, F_OK)) {
                    Error("File do not exist!");
                    glfwTerminate();
                    break;
                    exit( EXIT_SUCCESS );
                }
                
                tsunamiReadFile(baseResultName,frame*nout,U,V,E,nElem);
            }}
            
            glfwGetWindowSize( &width, &height );
            height = height > 0 ? height : 1;
            glViewport( 0, 0, width, height );
            
            glClearColor(0.027, 0.027, 0.171, 0); // Nice Dark blue for the sky
            
            /////////////////////////////////////////////
            //printf("string = %s \n", test);
            
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            
            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            gluPerspective(65.0f,(GLfloat)width/(GLfloat)height,1.0f,100.0f);
            
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            gluLookAt(0.0f,1.0f,0.0f,0.0f, 20.0f, 0.0f,0.0f,0.0f,1.0f);
            glTranslatef(0.0f,-zoom+14,0.0f);
            double tt = 1; // add
            
            glRotatef(0.3f*(GLfloat)389 + (GLfloat)tt*15.0f,0.0f,0.0f,1.0f);
            glScalef( 1, 1, 1);
            
            if (glfwGetKey(GLFW_KEY_RIGHT) == GLFW_PRESS) {
                coefXR +=GLFW_KEY_RIGHT;
                glRotated(yy + coefYU*0.005 - coefYD*0.005, 0.0, 1.0, 0.0);
                glRotated(xx - coefXL*0.005 + coefXR*0.005, 1.0, 0.0, 0.0);
                
            }
            else if (glfwGetKey(GLFW_KEY_LEFT) == GLFW_PRESS) {
                coefXL +=GLFW_KEY_LEFT;
                glRotated(yy + coefYU*0.005 - coefYD*0.005, 0.0, 1.0, 0.0);
                glRotated(xx - coefXL*0.005 + coefXR*0.005, 1.0, 0.0, 0.0);
            }
            else if (glfwGetKey(GLFW_KEY_UP) == GLFW_PRESS) {
                coefYU +=GLFW_KEY_UP;
                glRotated(yy + coefYU*0.005 - coefYD*0.005, 0.0, 1.0, 0.0);
                glRotated(xx - coefXL*0.005 + coefXR*0.005, 1.0, 0.0, 0.0);
                
                
            }
            
            else if (glfwGetKey(GLFW_KEY_DOWN) == GLFW_PRESS) {
                coefYD +=GLFW_KEY_DOWN;
                glRotated(yy + coefYU*0.005 - coefYD*0.005, 0.0, 1.0, 0.0);
                glRotated(xx - coefXL*0.005 + coefXR*0.005, 1.0, 0.0, 0.0);
                
            }
            
            else {
                glRotated(yy + coefYU*0.005 - coefYD*0.005, 0.0, 1.0, 0.0);
                glRotated(xx - coefXL*0.005 + coefXR*0.005, 1.0, 0.0, 0.0);
            }
            
            
            
            
            GLUquadricObj *quadratic = gluNewQuadric();
            gluQuadricNormals(quadratic, GLU_SMOOTH);
            glColor3f(0.906,0.82,0.527); // add better color for continent ;)
            
            gluSphere(quadratic,5.95,400,200);
     
            
            if (glfwGetKey('Z') == GLFW_PRESS && vagueFactor <= 0.1) vagueFactor += 0.01;
            else if (glfwGetKey('A') == GLFW_PRESS && vagueFactor >= 0.01 ) vagueFactor -= 0.01;
            
            for (i=0; i < nElem; ++i) {
                for (j=0; j < cote; ++j) {
                    index = elem[cote*i+j];
                    double value = H[index]/BathMax;
                    value = E[numberNodes*i+j]*10;
                    if (value < 0) value = 0;
                    if (value > 1) value = 1;
                    colors[j*3+0] = 3.5*(value)*(value);
                    colors[j*3+1] = (1-value)*(value)*3.5;
                    colors[j*3+2] = (1-value)*(1-value);
                    double x = X[index];
                    double y = Y[index];
                    double Factor = (4*R*R + x*x + y*y)*(R/6);
                    coord[j*3+0] = 4*R*R * x / Factor + (vagueFactor)*value*4*R*R * x / Factor;
                    coord[j*3+1] = 4*R*R * y / Factor + (vagueFactor)*value*4*R*R * y / Factor;
                    coord[j*3+2] = (4*R*R - x*x - y*y)*R/ Factor + (vagueFactor)*value*(4*R*R - x*x - y*y)*R/ Factor;
                    colorsS[j*3+0] = 3.5*(value)*(value);
                    colorsS[j*3+1] = (1-value)*(value)*3.5;
                    colorsS[j*3+2] = (1-value)*(1-value);
                    coordS[j*3+0] = 4*R*R * x / Factor ;
                    coordS[j*3+1] = 4*R*R * y / Factor ;
                    coordS[j*3+2] = (4*R*R - x*x - y*y)*R/ Factor ;
                }
                
                
                glEnableClientState(GL_VERTEX_ARRAY);
                glEnableClientState(GL_COLOR_ARRAY);
                glEnableClientState(GL_NORMAL_ARRAY);
                glVertexPointer(3, GL_FLOAT, 0, coord);
                glNormalPointer(GL_FLOAT, 0, coord);
                glColorPointer(3, GL_FLOAT, 0, colors);
                
                glDrawArrays(GL_POLYGON, 0, cote);
                
                glDisableClientState(GL_NORMAL_ARRAY);
                glDisableClientState(GL_COLOR_ARRAY);
                glDisableClientState(GL_VERTEX_ARRAY);
                
                glColor3f(0.0, 0.0, 0.0);
                glEnableClientState(GL_VERTEX_ARRAY);
                
                glEnableClientState(GL_VERTEX_ARRAY);
                glEnableClientState(GL_COLOR_ARRAY);
                glEnableClientState(GL_NORMAL_ARRAY);
                glVertexPointer(3, GL_FLOAT, 0, coordS);
                glNormalPointer(GL_FLOAT, 0, coordS);
                glColorPointer(3, GL_FLOAT, 0, colorsS);
                
                glDrawArrays(GL_POLYGON, 0, cote);
                
                glDisableClientState(GL_NORMAL_ARRAY);
                glDisableClientState(GL_COLOR_ARRAY);
                glDisableClientState(GL_VERTEX_ARRAY);
                
                glColor3f(0.0, 0.0, 0.0);
                glEnableClientState(GL_VERTEX_ARRAY);
                
            }
            
            
            glfwSwapBuffers();
            free(chemin);
            chemin=NULL;
            
        }
        
    }
    while( glfwGetKey( GLFW_KEY_ESC ) != GLFW_PRESS && glfwGetWindowParam( GLFW_OPENED ));
    
    free(X); X = NULL;
    free(Y); Y = NULL;
    free(H); H = NULL;
    free(U); U = NULL;
    free(V); V = NULL;
    free(E); E = NULL;
    free(elem); elem = NULL;
    free(forme); forme = NULL;
    glfwTerminate();
    exit( EXIT_SUCCESS );
    
    
}
double tsunamiInitialConditionOkada(double x, double y)
{
    double R = 6371220;
    double x3d = 4*R*R*x / (4*R*R + x*x + y*y);
    double y3d = 4*R*R*y / (4*R*R + x*x + y*y);
    double z3d = R*(4*R*R - x*x - y*y) / (4*R*R + x*x + y*y);
    double lat = asin(z3d/R)*180/M_PI;
    double lon = atan2(y3d,x3d)*180/M_PI;
    double lonMin = 142;
    double lonMax = 143.75;
    double latMin = 35.9;
    double latMax = 39.5;
    double olon = (lonMin+lonMax)/2;
    double olat = (latMin+latMax)/2;
    double angle = -12.95*M_PI/180;
    double lon2 = olon + (lon-olon)*cos(angle) + (lat-olat)*sin(angle);
    double lat2 = olat - (lon-olon)*sin(angle) + (lat-olat)*cos(angle);
    if ( lon2 <= lonMax && lon2 >= lonMin &&
         lat2 >= latMin && lat2 <= latMax ) 
            return 1.0;
    else    return 0.0; 
     
}

void tsunamiWriteFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem, int nsub)
{
    int i,j;
    const char *basename = "%s-%08d.txt";
    char filename[256];
    sprintf(filename,basename,baseResultName,iter);
    FILE* file = fopen(filename,"w");
    fprintf(file, "Number of elem %d \n", nelem);
    fprintf(file, "Number of local values per element %d \n", nsub);
    for (i = 0; i < nelem; ++i) {
    	for (j = 0; j < nsub; ++j) {
        	int index = i*nsub+j;
        	fprintf(file,"%d;%d;%le;%le;%le;\n",i,j,U[index],V[index],E[index]); }}
    fclose(file);
}

int tsunamiReadFile(const char *baseResultName, int iter, double *U, double *V, double *E, int nelem)
{
    int i,j,trash,nelemFile,nsub;
    const char *basename = "%s-%08d.txt";
    char filename[256];
      sprintf(filename,basename,baseResultName,iter);
      FILE* file = fopen(filename,"r");
      fscanf(file, "Number of elem %d \n", &nelemFile);
    fscanf(file, "Number of local values per element %d \n", &nsub);
    if (nelem != nelemFile) {
        printf("Error : wrong data file %d %d:-) \n",nelem,nelemFile);
        exit(0); }
    for (i = 0; i < nelem; ++i) {
    	for (j = 0; j < nsub; ++j) {
        	int index = i*nsub+j;
        	fscanf(file,"%d;%d;%le;%le;%le;\n",&trash,&trash,&U[index],&V[index],&E[index]); }}
    
    fclose(file);
    return nsub;
}

