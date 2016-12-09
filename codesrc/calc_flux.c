#include <stdio.h>
#include <math.h>
#include "time.h"
#include <unistd.h>
#include <stdlib.h>
#include "consts.h"
#include "get_freqs.c"
#include <getopt.h>

/* A cleaned-up version of calcFlux.c.
 * Started 5/25/2016, Austin Sousa -- asousa@stanford.edu
 *  -Added code to get frequencies from directory listing
 *  -Deleted alpha file stuff (no resuming)
 */ 

#define   NAPe    2.71828182845905

// -----------------------------------------------
// GLOBAL VARIABLES:
// -----------------------------------------------

double      L_TARG;
int nancounter;


// -----------------------------------------------
// Constants to perform Gauss quadrature integration
// -----------------------------------------------
// double t5[] =   {   -0.9061798459, 
//             -0.5384693101, 
//             0, 
//             0.9061798459,  
//             0.5384693101    };

// double beta5[] = {  0.2369268851,
//             0.4786286745, 
//             0.56888888889, 
//             0.2369268851,
//             0.4786286745    };

double beta5[] ={0.56888889, 0.47862867,  0.47862867, 0.23692689,  0.23692689};
double t5[] =   {0.,         0.53846931, -0.53846931, 0.90617985, -0.90617985};


// -----------------------------------------------
// PROTOTYPES
// -----------------------------------------------

float *getArr(void);
void updateArr(float *arr1, float *arr2);
void compFlux(float *arr, float out_lat, float out_lon, int k, char *dir, char *flux_filename);
void readJ(float J[][100], char *filename);
double getJdiff(float J[][100], double E, double alpha_lc);



/* -----------------------------------------------
 * calc_flux Main
 *
 * - Calculates differential electron flux due to
 *   the deflection matrices (pN*, pS* files), based
 *   on an initial magnetosphere population (flux_filename).
 *   
 * Inputs:
 *  -dir:    Directory to folder of pN, pS files
 *  -L_TARG: L-shell to compute on
 *  -flux_filename: Path to the flux file (EQFLUXMA.dat)
 */
int main(int argc, char *argv[])
{
  FILE *inPtr, *alphaPtr;
  int numFreqs, i, k, m, nin, ei, ti;
  char sysCmd[512], filename[64], *NS, alphaFile[128];
  // char *flux_filename;
  float L, *arr, *arrarr[128]; // <- array of pointers to precip arrays
  float J[100][100];
  double Jext;
  int *freqs;
  int num_freqs;

  nancounter = 0;
  int opt =0;
  int opt_index = 0;

  // Default inputs:
  char *dir = "outputs";
  char *out_dir = "outputs";
  char *flux_filename = "EQFLUXMA.dat";
  float out_lat = 0;
  float out_lon = 0;

  // Parse input arguments:
  static struct option long_options[] =
  {
      {"p_dir",       required_argument,    0, 'd'},
      {"out_dir",     required_argument,    0, 'e'},
      {"out_lat",     required_argument,    0, 'f'},
      {"out_lon",     required_argument,    0, 'g'},
      {"flux_file",   required_argument,    0, 'h'},
      {0, 0, 0, 0}
  };

  while (opt != -1) {
      opt = getopt_long (argc, argv, "d:e:f:g:h:", long_options, &opt_index);
      // cout << "opt is " << opt << "\n";
      switch(opt) {
          case 0:
          if (long_options[opt_index].flag != 0)      break;
          case 'd':   // p_dir
            dir = optarg;                             break;
          case 'e':
            out_dir = optarg;                         break;
          case 'f':
            out_lat  = atof(optarg);                  break;
          case 'g':
            out_lon = atof(optarg);                   break;
          case 'h':
            flux_filename = optarg;                   break;
          case '?':
               printf("\nUnknown option: %s\n",opt);  break;
      }
  }




  // Get L-shell at output:
  L_TARG = round(100.*((R_E + H_IONO)/(R_E*pow(cos(D2R*out_lat),2))))/100.;


  printf(" ------ 2D WIPP Code ------ \n");
  printf("         calc_flux          \n");
  printf(" -------------------------- \n\n");

  printf("input directory:\t%s\n",dir);
  printf("output directory:\t%s\n",out_dir);
  printf("out lat:\t %g\tlon:\t%g\n",out_lat, out_lon);
  printf("L-shell:\t%g\n",L_TARG);
  printf("flux file:\t%s\n",flux_filename);




  printf("\n\aWill do %d time steps\n", NUM_STEPS);


  freqs = get_freqs_at(dir, out_lon, &num_freqs);

  printf("\n\aWill do %d frequencies\n",num_freqs);  


  sprintf(alphaFile, "%s/alpha_%g_%s", out_dir, L_TARG, "N" );
  printf("alphaFile: %s\n",alphaFile);
  arr = getArr();


  for(i=0; i<num_freqs; i++) {
    // printf("Frequency: %d \n",freqs[i]);
    for(k=0; k<2; k++) {

        NS = (k==0) ? "N" : "S";

          // Initialize the array (if we're the first)
          if(i==0)  arrarr[k]=getArr();

              sprintf(filename,"%s/p%s_%g_%g_%d.dat",dir, NS, out_lat, out_lon, freqs[i]);
              // sprintf(filename,"%s/p%s%d_%g.dat",dir,NS,freqs[i],L_TARG);
              // printf("i: %d, k: %d, filename: %s\n", i, k, filename);
              inPtr = fopen(filename, "r");
              if (inPtr != NULL) {
                printf("opened %s, pointer is %d\n",filename,inPtr);
                nin = fread(arr, sizeof(float), NUM_E*NUM_STEPS, inPtr);
                fclose(inPtr);
                // printf("read %d values\n",nin);
                // Add to rolling sum from previous files
                updateArr( arrarr[k] , arr );
              } else {
                printf("could not open file at %s\n",filename);
              }
            } // for(k ... )  N/S - hemisphere

      } // freqs
  

  // Compute flux for north and south
  for(k=0; k<2; k++) {
    NS = (k == 0) ? "N" : "S";
    printf("Calling compFlux... hemisphere = %s\n",NS);
    compFlux(arrarr[k], out_lat, out_lon, k, out_dir, flux_filename);
  } // N/S - hemisphere


  return 0;
}



/*
 * FUNCTION: compFlux
 * ------------------
 * This function will open a file with the appropriate hemisphere and
 * L-shell inserted into the name, calculate the precipitated flux  
 * due to each cell in the dAlpha_RMS matrix and write it into the 
 * file.
 *
 */

void compFlux(float *arr, float out_lat, float out_lon, int k, char *dir, char *flux_filename)
{
  FILE *phiPtr, *QPtr, *NPtr, *alphaPtr;
  int ei, ti, i, nout;
  double da_pk, P, Q, epsm, alpha_eq, v, crunch;
  double I=0.0, x, g, field, Phi_p, b=1.0, vcm, fcgs;
  double v_tot_arr[NUM_E], E_tot_arr[NUM_E], dE_arr[NUM_E], gamma;
  double Jdiff[NUM_E];
  float Phi_float, alpha, J[100][100];
  char *NS, PhiFile[128], QFile[128], NFile[128], alphaFile[128];

  double t1, t2, t3, t4, t6;
  double max_alpha = 0;

  // Open up Phi file for writing
  if(k==0) {NS = "N";} else {NS="S";}  

  sprintf( PhiFile, "%sphi_%g_%g_%s.dat", dir, out_lat, out_lon, NS );
  sprintf( QFile,   "%sQ_%g_%g_%s.dat", dir, out_lat, out_lon, NS );
  sprintf( NFile,   "%sN_%g_%g_%s.dat", dir, out_lat, out_lon, NS );
  sprintf( alphaFile,   "%salpha_%g_%g_%s.dat", dir, out_lat, out_lon, NS );  

  printf("writing %s\n", PhiFile);

  if( (phiPtr=fopen(PhiFile, "w"))==NULL ) {
    printf("\nProblem opening %s\n", PhiFile);
    exit(0);
   }

  if( (alphaPtr=fopen(alphaFile, "w"))==NULL ) {
    printf("\nProblem opening %s\n", alphaFile);
   exit(0);
  }

  
  epsm = (1/L_TARG)*(R_E+H_IONO)/R_E;
  // (1/sin(alpha_lc)^2) -- geometric focusing term (Bortnik 5.2)
  crunch  = sqrt(1+3*(1-epsm))/pow(epsm,3); 
  // Loss-cone angle at equator:
  alpha_eq  = asin(sqrt( 1.0/crunch ));
  
  printf("alpha_eq: %g\n",R2D*alpha_eq);
  printf("crunch: %g\n",crunch);
  // Load flux file:
  readJ(J, flux_filename);

  // Precalculate energy and velocity values
  for(i=0; i<NUM_E; i++) {
    E_tot_arr[i] = pow(10, (E_EXP_BOT+(DE_EXP/2)+DE_EXP*i) ); //E in eV
    Jdiff[i] = getJdiff( J, E_tot_arr[i], alpha_eq );
    v_tot_arr[i] = C*sqrt(1 - pow( (E_EL/(E_EL+E_tot_arr[i])) ,2) );

    // Energy differential dE in keV
    dE_arr[i] = 1e-3 * pow(10, (E_EXP_BOT+ (DE_EXP/2))) * 
      exp(DE_EXP*i / log10(NAPe)) * DE_EXP / log10(NAPe);
  }

  // Loop over energies:
  for(ei=0; ei<NUM_E; ei++) {

    if(ALPHA_DISTRIBUTION) {
      // Suprathermal velocity distribution
      // (I forget where this one is from -- should be consistent
      //  with the damping code. --aps 12.2016)

      v = v_tot_arr[ei];
      vcm = v*100;      // v in cm for distrib fn calculation
      gamma = 1.0/sqrt( 1 - v*v/(C*C) );
      fcgs =  4.9e5/pow( (vcm*gamma) ,4) - 
              8.3e14/pow( (vcm*gamma) ,5) + 
              5.4e23/pow( (vcm*gamma) ,6);
     

      b = (v*v/M_EL)*pow( sqrt(1 - (v*v)/(C*C)), 3) * 1.6e-8 * fcgs;

      // b = 1e8 / pow(E_tot_arr[i],2);
    } else {
      b = Jdiff[ei]*1000;
    }




    // Loop over timesteps:
    for(ti=0; ti<NUM_STEPS; ti++) {
      
      alpha = sqrt( arr[ei*NUM_STEPS+ti] );

      if (alpha > max_alpha) { max_alpha = alpha; }

      if (alpha > 0) {

        // printf("ei: %d ti: %d alpha: %g\n",ei, ti, alpha);
      }
      nout=fwrite(&alpha, sizeof(float), 1, alphaPtr);      

      // Peak change in pitch-angle at this time and energy:
      da_pk = sqrt(2.0)*alpha;  // sqrt(2)*alpha_RMS = peak

      // Integrate the distribution from (alpha = 0 to alpha_lc) 
      // using Gaussian Quadrature (5th order)

      // First we need to convert our limits from (alpha_lc - da_pk,  alpha_lc) to (0, 1):
      P = da_pk/2.;             //[ alpha_lc - (alpha_lc-da_pk) ] /2
      Q = alpha_eq - da_pk/2.;  //[ alpha_lc + (alpha_lc-da_pk) ] /2

      I = 0.0;


      if(da_pk != 0) {
        // Integrate wrt alpha:
        for(i=0; i<5; i++) {
          x = P*t5[i] + Q ;

          if(ALPHA_DISTRIBUTION) {

            // Square distribution
            // g = (P/PI) * sin(2.0*x) * (asin((x-alpha_eq)/da_pk)+ (PI/2.0) );

            // same thing - but avoids NaNs when da_pk is very very small:
            g = (P/PI) * sin(2.0*x) * (asin(t5[i]/2. - 0.5) + PI/2.0);          

          } else {
            // S
            // g = (P/PI)*sin(2.0*x)*((x - alpha_eq)*( asin((x-alpha_eq)/da_pk)+ (PI/2.0) ) +
            //   sqrt(da_pk*da_pk-pow((x-alpha_eq),2)));

            // same thing, but simplified to ditch numerical errors when t1 is small:
            t1 = (x - alpha_eq);
            t2 = fabs(da_pk)*sqrt(1 -0.25*pow(t5[i] - 1, 2));
            g = (P/PI)*sin(2.0*x)*(t1*( asin(t5[i]/2. - 0.5)+ (PI/2.0) ) + t2);
          }; 

      // if isnan(g) {
      //   // printf("G ISNAN at t: %d e: %d alpha: %g g: %g i=%d\n", ti, ei, alpha, g, i);
      //   // printf("x: %g alpha_eq: %g da_pk: %g P %g\n",x,alpha_eq, da_pk, P);
      //   printf("t2: %g t3: %g diff: %g \n",t2,t3, fabs(t2 - t3));
      //      };
        I += ( beta5[i]*g );

        } // for(i ... ) -> Gauss quad integration
      } // if da_pk != 0

      Phi_p = PI*crunch*b*I;

      // printf("I: %g\n",I);
      Phi_float = ( float ) Phi_p;
      
      // if isnan(I) { printf("I ISNAN at t: %d e: %d alpha: %g g: %g i=%d\n", ti, ei, alpha, g, i); };



      if isnan(Phi_p) { 
        nancounter=nancounter + 1;
        //printf("Total NaNs: %i\n",nancounter);
      };

      nout=fwrite(&Phi_float, sizeof(float), 1, phiPtr);
      
    } // for(ti ... )
  } // for(ei ... )
  
  fclose(phiPtr);
  fclose(alphaPtr);

  printf("Total NaNs: %i\n",nancounter);
  printf("max pitch-angle deflection: %g deg\n",R2D*max_alpha);

//   /// HEY AUSTIN, YOU HAVEN'T TESTED THIS PART YET (6.2.2016)
  // Now calculate Q and N
  // Need to integrate over E
  // float Qarr[NUM_E][NUM_TIMES];
  // float Qarr[NUM_E][NUM_TIMES];
  // for(ti=0; ti<NUM_TIMES; ti++) {
  //  for(ei=0; ei<NUM_E; ei++) {
  //     Qarr[ti] += (float)( Phi_p * E_tot_arr[ei] * dE_arr[ei] * 1.602e-9); // joules to milliergs
  //     Narr[ti] += (float)( Phi_p * dE_arr[ei]); //eV->keV
  //    } // ei
  // }  // ti

//   nout=fwrite(Qarr, sizeof(float), (NUM_TIMES), QPtr);
//   if(nout!=NUM_LONGS*NUM_TIMES) printf("\n\aProblem writing Q\n");
//   nout=fwrite(Narr, sizeof(float), (NUM_TIMES), NPtr);
//   if(nout!=NUM_LONGS*NUM_TIMES) printf("\n\aProblem writing N\n");

//   fclose(QPtr);
//   fclose(NPtr);

}

/*
 * FUNCTION: getJdiff
 * ------------------
 * Using the AE8 data stored in J, calculate the differential flux 
 * by taking the (energy) derivative of the integral flux, dividing
 * by the total solid angle and extrapolating down in energy if 
 * need be.
 *
 */
double  getJdiff(float J[][100], double E, double alpha_lc)
{
  int row, i, topCol, botE;
  double J1, J2, I, x1, x2, y1, y2, m, c, x_ext, y_ext, J_ext;

  row = (int)floor((L_TARG+0.11 - J[1][0])/0.1); // to make sure! 
  
  // if(  fabs((double)J[row][0]-L_TARG) > 1e-3   ) 
  //   printf("\nL-shell not matching data\n\a");

  I = PI * cos(alpha_lc) * (PI - 2*alpha_lc);

  // Find column corresponding to highest energy value
  for(i=0; i<100; i++) {
    if(J[0][i+1] < 0.01) { 
      topCol = i; 
      break; 
    }
  }



  // Case 1. E < 100 keV
  // -------------------

  if( E <= 1e5 ) {
 
    // diff flux @ 100 keV and 200 keV
    J1 = 1e-6*fabs(J[row][2] - J[row][1]) / (J[0][2] - J[0][1]); 
    J2 = ((1e-6*fabs(J[row][3] - J[row][2]) / (J[0][3] - J[0][2])) 
      + J1 )/2; // central difference

    // do extrapolation in log-log space for best fit 
    x1 = log10( J[0][1]*1e6 );
    x2 = log10( J[0][2]*1e6 );
    y1 = log10( J1 );
    y2 = log10( J2 );

    m = (y2-y1)/(x2-x1);            // gradient of line
    c = (y1*x2 - y2*x1)/(x2-x1) ;   // offset of line, i.e.
    
    // y = m*x + c
    x_ext = log10( E );
    y_ext = m*x_ext + c;
    J_ext = pow(10, y_ext);

    return (J_ext/I);

  }


  

  // Case 2. E > 7 MeV
  // -----------------

  if( E >= 7e6 ) {
  
    // If flux at 7 Mev = 0, flux above it is zero too
    if( J[row][topCol]==0 )  return 0;

    // Otherwise need to extrapolate as in case 1.
    // diff flux @ 6.5 MeV and 7 MeV
    J2 = 1e-6*fabs( J[row][topCol] - J[row][topCol-1] ) 
      / (J[0][topCol] - J[0][topCol-1]); 

    J1 = ((1e-6*fabs( J[row][topCol-1] - J[row][topCol-2]) / 
       (J[0][topCol-1] - J[0][topCol-2]) ) + J2 )/2; // cdiff

    // do extrapolation in log-log space for best fit 
    x1 = log10( J[0][topCol-1]*1e6 );
    x2 = log10( J[0][topCol]*1e6 );
    y1 = log10( J1 );
    y2 = log10( J2 );

    m = (y2-y1)/(x2-x1);        // gradient of line
    c = (y1*x2 - y2*x1)/(x2-x1) ;   // offset of line, i.e.
                    // y = m*x + c
    x_ext = log10( E );
    y_ext = m*x_ext + c;
    J_ext = pow(10, y_ext);

    if(J_ext < 1e-10 ) J_ext = 0.0;

    return (J_ext/I);
  }


  // Case 3. 100 keV < E < 7 MeV
  if( E<7e6 && E>1e5 ) {


    // Find column corresponding lower energy value
    for(i=1; i<100; i++) {
      if( (J[0][i+1]*1e6) > E ) { 
    botE = i; 
    break; 
      }
    }


    // central diff flux @ lower and higher energies
    J1 = ( (1e-6 * fabs( J[row][botE] - J[row][botE-1] )
        / ( J[0][botE] - J[0][botE-1] ) ) + 
       (1e-6 * fabs( J[row][botE+1] - J[row][botE] )
        / ( J[0][botE+1] - J[0][botE] ) )  ) / 2;

    J2 = ( (1e-6 * fabs( J[row][botE+1] - J[row][botE] )
        / ( J[0][botE+1] - J[0][botE] ) ) + 
       (1e-6 * fabs( J[row][botE+2] - J[row][botE+1] )
        / ( J[0][botE+2] - J[0][botE+1] ) )  ) / 2;

    if(botE == 1)
      J1 =  (1e-6 * fabs( J[row][botE+1] - J[row][botE] )
          / ( J[0][botE+1] - J[0][botE] ) );
    
    if(botE == (topCol-1))
      J2 = (1e-6 * fabs( J[row][botE+1] - J[row][botE] )
        / ( J[0][botE+1] - J[0][botE] ) );
    



    // If J1 = J2 = 0, interpolated value also 0
    if( J1==0 && J2==0 ) return 0;



    // If only J2 = 0, do linear interpolation
    if( J2 == 0 ) {
      J_ext = J1*( ( J[0][botE+1]-(E*1e-6) )/
           ( J[0][botE+1] - J[0][botE] ) );
      return (J_ext/I);
    }



    // Otherwise interpolate as in case 1 (log-log space)

    x1 = log10( J[0][botE]*1e6 );
    x2 = log10( J[0][botE+1]*1e6 );
    y1 = log10( J1 );
    y2 = log10( J2 );

    m = (y2-y1)/(x2-x1);        // gradient of line
    c = (y1*x2 - y2*x1)/(x2-x1) ;   // offset of line, i.e.
                    // y = m*x + c
    x_ext = log10( E );
    y_ext = m*x_ext + c;
    J_ext = pow(10, y_ext);

    return (J_ext/I);
  }

}


// -----------------------------------------------
// UTILITIES
// -----------------------------------------------

/*
 * FUNCTION: readJ
 * ---------------
 * This function simply looks for a file 'filename' and reads 
 * it in.  The columns are energies and the rows are L-shells.
 * The first column is just a list of L-shells and the first row 
 * is just a list of energies (in MeV).
 *
 */
void readJ(float J[][100], char *filename)
{
  // char *filename;
  FILE *filePtr;
  int i,j;

  // filename = "EQFLUXMA.dat";

  if( (filePtr = fopen( filename ,"r")) == NULL ) {
    printf("Error opening the flux file! path: %s\n",filename);
    exit(1);
  }
  

  // INITIALIZE
  for(i=0; i<100; i++) {
    for(j=0; j<100; j++) {
      J[i][j] = 0.0;
    }
  }

  // READ IN VALUES
  for(i=0; i<47; i++) {
    for(j=0; j<31; j++) {
      fscanf(filePtr, "%e", &(J[i][j]));
    }
  }

  fclose(filePtr);
}




/* 
 * FUNCTION: getArr
 * ----------------
 * This function simply allocates dynamically a block of memory the 
 * size of NUM_E * NUM_STEPS of type float, and initializes it.
 * It returns the pointer to the memory.
 *
 */
float *getArr(void)
{
  float *arr;
  int ei, ti;
    printf("calling getArr...\n");

  arr = (float *) malloc( NUM_E * NUM_STEPS * sizeof(float) );
  if(arr == NULL) {
    printf("\nProb assigning mem in calcFlux\n");
    exit(0);
  }

  for(ei=0; ei<NUM_E; ei++) {
    for(ti=0; ti<NUM_STEPS; ti++) {
      arr[ei*NUM_STEPS+ti] = 0.0;
    }
  }
  //printf("Finishing getArr...\n");
  return arr;

}


/* 
 * FUNCTION: updateArr
 * -------------------
 * This function updates the values of arr1 with those of arr2.
 *
 */
void updateArr(float *arr1, float *arr2)
{
  int ei, ti;
  //printf("starting updateArr\n");
  for(ei=0; ei<NUM_E; ei++) {
    for(ti=0; ti<NUM_STEPS; ti++) {
      if(arr2[ei*NUM_STEPS+ti]>0.0)
        arr1[ei*NUM_STEPS+ti] += arr2[ei*NUM_STEPS+ti];
    }
  }
}
