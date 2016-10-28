#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <stdio.h>
#include <math.h>

const double kMicro = 1.0e-6;
double getTime()
{
        struct timeval TV;

        const int RC = gettimeofday(&TV, NULL);
        if(RC == -1)
        {
                printf("ERROR: Bad call to gettimeofday\n");
                return(-1);
        }

        return( ((double)TV.tv_sec) + kMicro * ((double)TV.tv_usec) );
}  

enum Direction{
    east=0, 
    west,
    north, 
    south, 
    up, 
    down 
};

double residual(double*** U, const int nx, const int ny, const int nz){
  double err = 0;
  for (int k=1; k<nz-1; k++)
    for (int j=1; j<ny-1; j++)
      for (int i=1; i<nx-1; i++){
        double du =  (U[k-1][j][i] + U[k+1][j][i] + U[k][j-1][i] + U[k][j+1][i] + U[k][j][i-1] + U[k][j][i+1] - 6 * U[k][j][i]);
        err = err + du * du;
      }
  return err;
}

int main(int argc, char* argv[])
{
    int myRank =0;
    int size =0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int N =128;
    int nIters=100;
    int px=1, py=1, pz=1;
    int rankx, ranky, rankz;
    int direct[6];     // east, west, north, south, up, down

    int argCount=0;
    while(++argCount <argc) {
   	if(!strcmp(argv[argCount], "-px")) px = atoi(argv[++argCount]);
	if(!strcmp(argv[argCount], "-py")) py = atoi(argv[++argCount]);
   	if(!strcmp(argv[argCount], "-pz")) pz = atoi(argv[++argCount]);
   	if(!strcmp(argv[argCount], "-N")) N = atoi(argv[++argCount]);
   	if(!strcmp(argv[argCount], "-i")) nIters = atoi(argv[++argCount]);
    }

    rankz= myRank/(px*py);
    ranky= (myRank%(px*py))/px;
    rankx= (myRank%(px*py))%px;

    direct[east] = rankx<(px-1)?myRank + 1:-1;
    direct[west] = rankx > 0 ?  myRank - 1:-1;
    direct[north] = ranky<(py-1)? myRank + px :-1;
    direct[south] = ranky > 0 ? myRank - px :-1;
    direct[up] = rankz<(pz-1)? myRank + px*py:-1;
    direct[down] = rankz > 0 ? myRank - px*py:-1;
    
//initialize the data
    double ***U, ***Un;
    double *messageEast=NULL, *messageWest=NULL, *messageNorth=NULL, *messageSouth=NULL, *messageUp=NULL, *messageDown=NULL;
    double *messageEast_recv=NULL, *messageWest_recv=NULL, *messageNorth_recv=NULL, *messageSouth_recv=NULL, *messageUp_recv=NULL, *messageDown_recv=NULL;
    int i, j, k;
    int nx, ny, nz;
    nx = N/px;
    ny = N/py;
    nz = N/pz;
    nx+=2;
    ny+=2;
    nz+=2;

    U = (double ***) malloc(nz*sizeof(double**));
    Un = (double ***) malloc(nz*sizeof(double**));

    for(k=0; k<nz; k++){
       U[k] = (double **) malloc(ny*sizeof(double*));
       Un[k] = (double **) malloc(ny*sizeof(double*));
    }

    for(k=0; k<nz; k++)
       for(j=0; j<ny; j++){
          U[k][j] = (double *) malloc(nx*sizeof(double));
          Un[k][j] = (double *) malloc(nx*sizeof(double));
    }

    for (k=0; k<nz; k++)
       for (j=0; j<ny; j++)
          for(i=0; i<nx; i++){
             if(((rankz==0)&&(k==0)) ||((rankz==pz-1)&&(k==nz-1)) ||((ranky==0)&&(j==0)) ||((ranky==py-1)&&(j==ny-1))||((rankx==0)&&(i==0)) ||((rankx==px-1)&&(i==nx-1))){
               U[k][j][i] = 0;
               Un[k][j][i] = 0;
             }else{
               U[k][j][i] = 1;
               Un[k][j][i] = 1;
             }
          }

    messageEast = (double *) malloc(ny*nz*sizeof(double));
    messageWest = (double *) malloc(ny*nz*sizeof(double));
    messageNorth = (double *) malloc(nx*nz*sizeof(double));
    messageSouth = (double *) malloc(nx*nz*sizeof(double));
    messageUp = (double *) malloc(nx*ny*sizeof(double));
    messageDown = (double *) malloc(nx*ny*sizeof(double));

    messageEast_recv = (double *) malloc(ny*nz*sizeof(double));
    messageWest_recv = (double *) malloc(ny*nz*sizeof(double));
    messageNorth_recv = (double *) malloc(nx*nz*sizeof(double));
    messageSouth_recv = (double *) malloc(nx*nz*sizeof(double));
    messageUp_recv = (double *) malloc(nx*ny*sizeof(double));
    messageDown_recv = (double *) malloc(nx*ny*sizeof(double));


    bzero(messageEast, ny*nz*sizeof(double));
    bzero(messageWest, ny*nz*sizeof(double));
    bzero(messageNorth, nx*nz*sizeof(double));
    bzero(messageSouth, nx*nz*sizeof(double));
    bzero(messageUp, nx*ny*sizeof(double));
    bzero(messageDown, nx*ny*sizeof(double));

    bzero(messageEast_recv, ny*nz*sizeof(double));
    bzero(messageWest_recv, ny*nz*sizeof(double));
    bzero(messageNorth_recv, nx*nz*sizeof(double));
    bzero(messageSouth_recv, nx*nz*sizeof(double));
    bzero(messageUp_recv, nx*ny*sizeof(double));
    bzero(messageDown_recv, nx*ny*sizeof(double));

    MPI_Request request[6];
    int southTAG=111, northTAG=112, eastTAG=113, westTAG=114, downTAG=115, upTAG=116;
    MPI_Barrier(MPI_COMM_WORLD);
    double execTime = -getTime();
     int request_count =0;

    double c=1.0/6.0;

    for (int iter=0; iter<nIters; iter++){
     for (int k=1; k< nz-1; k++) {
       for (int j=1; j<ny-1; j++) {
          double* Un0 = &Un[k][j][0];
          double *up  = &U[k-1][j][0];
          double *down = &U[k+1][j][0];
          double *east = &U[k][j][-1];
          double *west = &U[k][j][1];
          double *north = &U[k][j+1][0];
          double *south = &U[k][j-1][0];
          for (int i = 1; i < nx-1; i++){
            *(Un0+i) =  c*(*(up+i) + *(down+i) + *(east+i) + *(west+i) + *(north+i) + *(south+i));
          }
       }
     }

     double ***temp = NULL;
     temp = U;
     U = Un;
     Un = temp;

    request_count =0;
    if(direct[down] !=-1) {
       double* idx;
       double* ptr= messageDown;
       for(int y=0; y < ny; y++) {
           idx = &U[1][y][0];
           memcpy(ptr, idx, sizeof(double)*nx);
           ptr += nx;
       }
    }

    if(direct[up] !=-1) {
       double* idx;
       double* ptr= messageUp;
       for(int y=0; y < ny; y++) {
           idx = &U[nz-2][y][0];
           memcpy(ptr, idx, sizeof(double)*nx);
           ptr += nx;
       }
    }

    if(direct[east] !=-1) {
       double* idx;
       double* ptr= messageEast;
       for(int z=0; z < nz; z++) {
           for(int y=0; y < ny; y++) {
              idx = &U[z][y][nx-2];
              memcpy(ptr, idx, sizeof(double)*1);
              ptr ++;
           }
       }
    }

    if(direct[west] !=-1) {
       double* idx;
       double* ptr= messageWest;
       for(int z=0; z < nz; z++) {
           for(int y=0; y < ny; y++) {
              idx = &U[z][y][1];
              memcpy(ptr, idx, sizeof(double)*1);
              ptr ++;
           }
       }
    }

    if(direct[north] !=-1) {
       double* idx;
       double* ptr= messageNorth;
       for(int z=0; z < nz; z++) {
           idx = &U[z][ny-2][0];
           memcpy(ptr, idx, sizeof(double)*nx);
           ptr += nx;
       }
    }

    if(direct[south] !=-1) {
       double* idx;
       double* ptr= messageSouth;
       for(int z=0; z < nz; z++) {
           idx = &U[z][1][0];
           memcpy(ptr, idx, sizeof(double)*nx);
           ptr += nx;
       }
    }

        //communicate
     if(direct[down] !=-1){
        MPI_Isend(messageDown, nx*ny, MPI_DOUBLE, direct[down], downTAG, MPI_COMM_WORLD, &request[request_count]);
        request_count++;
     }

     if(direct[up] !=-1){
        MPI_Isend(messageUp, nx*ny, MPI_DOUBLE, direct[up], upTAG, MPI_COMM_WORLD, &request[request_count]);
        request_count++;
     }
     if(direct[east] !=-1){
        MPI_Isend(messageEast, nz*ny, MPI_DOUBLE, direct[east], eastTAG, MPI_COMM_WORLD, &request[request_count]);
        request_count++;
     }

     if(direct[west] !=-1){
        MPI_Isend(messageWest, nz*ny, MPI_DOUBLE, direct[west], westTAG, MPI_COMM_WORLD, &request[request_count]);
        request_count++;
     }
     if(direct[north] !=-1){
        MPI_Isend(messageNorth, nz*nx, MPI_DOUBLE, direct[north], northTAG, MPI_COMM_WORLD, &request[request_count]);
        request_count++;
     }

     if(direct[south] !=-1){
        MPI_Isend(messageSouth, nz*nx, MPI_DOUBLE, direct[south], southTAG, MPI_COMM_WORLD, &request[request_count]);
        request_count++;
     }

     if(direct[down] !=-1){
        MPI_Recv(messageDown_recv, nx*ny, MPI_DOUBLE, direct[down], upTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     }

     if(direct[up] !=-1){
        MPI_Recv(messageUp_recv, nx*ny, MPI_DOUBLE, direct[up], downTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     }
     if(direct[east] !=-1){
        MPI_Recv(messageEast_recv, nz*ny, MPI_DOUBLE, direct[east], westTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     }

     if(direct[west] !=-1){
        MPI_Recv(messageWest_recv, nz*ny, MPI_DOUBLE, direct[west], eastTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     }
     if(direct[north] !=-1){
        MPI_Recv(messageNorth_recv, nz*nx, MPI_DOUBLE, direct[north], southTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     }

     if(direct[south] !=-1){
        MPI_Recv(messageSouth_recv, nz*nx, MPI_DOUBLE, direct[south], northTAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
     }
     MPI_Waitall(request_count, request, MPI_STATUS_IGNORE);

    if(direct[down] !=-1) {
       double* idx;
       double* ptr= messageDown_recv;
       for(int y=0; y < ny; y++) {
           idx = &U[0][y][0];
           memcpy(idx, ptr, sizeof(double)*nx);
           ptr += nx;
       }
    }

    if(direct[up] !=-1) {
       double* idx;
       double* ptr= messageUp_recv;
       for(int y=0; y < ny; y++) {
           idx = &U[nz-1][y][0];
           memcpy(idx, ptr, sizeof(double)*nx);
           ptr += nx;
       }
    }

    if(direct[east] !=-1) {
       double* idx;
       double* ptr= messageEast_recv;
       for(int z=0; z < nz; z++) {
           for(int y=0; y < ny; y++) {
              idx = &U[z][y][nx-1];
              memcpy(idx, ptr, sizeof(double)*1);
              ptr ++;
           }
       }
    }

    if(direct[west] !=-1) {
       double* idx;
       double* ptr= messageWest_recv;
       for(int z=0; z < nz; z++) {
           for(int y=0; y < ny; y++) {
              idx = &U[z][y][0];
              memcpy(idx, ptr, sizeof(double)*1);
              ptr ++;
           }
       }
    }

    if(direct[north] !=-1) {
       double* idx;
       double* ptr= messageNorth_recv;
       for(int z=0; z < nz; z++) {
           idx = &U[z][ny-1][0];
           memcpy(idx, ptr, sizeof(double)*nx);
           ptr += nx;
       }
    }

    if(direct[south] !=-1) {
       double* idx;
       double* ptr= messageSouth_recv;
       for(int z=0; z < nz; z++) {
           idx = &U[z][0][0];
           memcpy(idx, ptr, sizeof(double)*nx);
           ptr += nx;
       }
    }
 }


 MPI_Barrier(MPI_COMM_WORLD);
 execTime += getTime();

 double res = residual(U, nx, ny, nz);
 double  global_err=0;
 MPI_Reduce (&res, &global_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
 if(myRank ==0) {
   printf("Residual: %e  \n ", sqrt(global_err/((double)(N+1)*(double)(N+1)*(double)(N+1))));
 }

 if(myRank ==0) {
   double gflops = nIters*(double)N*(double)N*(double)N*8/(1.0e9);
   printf("Exec time collected by user: %e  seconds\n ",execTime);
   printf("GFLOP/S: %e  \n ",gflops/execTime);
 }
 MPI_Finalize();
 return 0;
}

