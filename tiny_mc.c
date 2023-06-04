/* Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)"
 * 1 W Point Source Heating in Infinite Isotropic Scattering Medium
 * http://omlc.ogi.edu/software/mc/tiny_mc.c
 *
 * Adaptado para CP2014, Nicolas Wolovick
 */

#define _XOPEN_SOURCE 500  // M_PI

#define VectSize 8 //Tamano del vector para SIMD

#include "params.h"
#include "wtime.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

char t1[] = "Tiny Monte Carlo by Scott Prahl (http://omlc.ogi.edu)";
char t2[] = "1 W Point Source Heating in Infinite Isotropic Scattering Medium";
char t3[] = "CPU version, adapted for PEAGPGPU by Gustavo Castellano"
            " and Nicolas Wolovick";


// global state, heat and heat square in each shell
static float heat[SHELLS];
static float heat2[SHELLS];

const float albedo = MU_S / (MU_S + MU_A);
const float shells_per_mfp = 1e4 / MICRONS_PER_SHELL / (MU_A + MU_S);
const float tol = 0.001f; //tolerancia para inicio de roulette
//const float n_roulette = logf(0.001f / 1.0f) / logf (20.0f/22.0f);
unsigned int n_roulette = 73;
static unsigned int g_seed;
unsigned int i;
unsigned int j;
float pi = (float) M_PI;
float dos_pi = 2* M_PI;



//################  FUNCIONES PARA FAST_RAND() ###########

// Used to seed the generator. 
          
void fast_srand(int seed) {
    g_seed = seed;
    }

// Compute a pseudorandom integer.
// Output value in range [0, 32767]

int fast_rand(void) {
    g_seed = (214013*g_seed+2531011);
    return (g_seed>>16)&0x7FFF;
    }    
//##########################################################

//################  FUNCION FAST_SIN() y FAST_COS() ###########

float fast_sin(float x) {
	const float B = 4.0f / pi;
	const float C = -4.0f / (pi * pi);
	const float P = 0.225f;
	
	float y = B * x + C * x * fabs(x);
		
	y = P * (y * fabs(y) - y) + y;
	
	return y; 
	}
	
float fast_cos(float x) {
	const float B = 4.0f / pi;
	const float C = -4.0f / (pi * pi);
	const float pi_2  = pi/2.0f;
	const float pi_32 = 3.0f*pi/2.0f;
	const float P = 0.225f;

	float y = (B * (x + pi_2) + C * (x + pi_2) * fabs(x + pi_2))*(x<=pi_2) + (B * (x - pi_32) + C * (x - pi_32) * fabs(x - pi_32))*(x>pi_2);
		
	y = P * (y * fabs(y) - y) + y;
	
	return y; 
	}	
	
	
//##########################################################

/***
 * Photon
 ***/

static void photon(void)
{


    /* launch */
    float x[VectSize] = {0.0f};
    float y[VectSize] = {0.0f};
    float z[VectSize] = {0.0f};
    float r[VectSize] = {0.0f};
    float u[VectSize] = {0.0f};
    float v[VectSize] = {0.0f};
    float w[VectSize] = {1.0f,1.0f,1.0f,1.0f ,1.0f,1.0f,1.0f,1.0f};//  ,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,  
    			// 1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f};
    
    float t[VectSize];
    float xi1[VectSize], xi2[VectSize];
    //float Rand[VectSize];
    
    //unsigned int n_rands = 219; //438; // 3 * n_roulette * 2;
    //float Rand[n_rands];
	float Rand1[VectSize];
	float Rand2[VectSize];
	float Rand3[VectSize];
	//float Rand4[VectSize];
    float rad[VectSize];
	float phi[VectSize];
    
    /* CUIDADO: Cambiar tambien en variable global "n_roulette" implicita en caso de tocar "weight" */
    float weight = 1.0f;
    
    unsigned int shell_2D[73][VectSize]; // Matriz para guarda posicion de interaccion de los fotones mientras no haya roulette
	unsigned int shell_1D[VectSize]; // Vector para guardar posicion de interaccion paso a paso con roulette



	/* #################  LOOP SIN ROULETTE  ####################*/
	for(j=0;j<n_roulette;j++){ // Iterar sobre los 73 pasos SIN ROULETTE
	
		
		/* ##### Creacion de 3 array de tamaño [VectSize] ##### */
		for (i=0;i<VectSize;i++) {
			Rand1[i] = fast_rand() / (float)32767.0 ;
			Rand2[i] = fast_rand() / (float)32767.0 ;
			Rand3[i] = fast_rand() / (float)32767.0 ;
			} 
		
		/* #### Loop principal de Photon() Vectrorizado sin Roulete #### */	
		for (i=0;i<VectSize;i++) {
			t[i] = -logf( Rand1[i] ) / (MU_A + MU_S);	// Primer Random ###		
		    x[i] += t[i] * u[i];
		    y[i] += t[i] * v[i];
		    z[i] += t[i] * w[i];

		    shell_2D[j][i] = sqrtf(x[i] * x[i] + y[i] * y[i] + z[i] * z[i] + r[i] * r[i]) * shells_per_mfp;
		    if (shell_2D[j][i] > SHELLS - 1) {
		        shell_2D[j][i] = SHELLS - 1;
		    	}
		
        	/* New direction - Distribution - Vectorized*/
	  	    rad[i] = sqrtf(Rand2[i]);				// Segundo Random ###
	  	    phi[i] = pi * (2.0f * Rand3[i] - 1.0f);	// Tercer Random ###
	  	    // Coseno Funcion Libreria
	  	    xi1[i] = rad[i] * cosf(phi[i]);	  	    
	  	    // Seno Implementacion in=situ de Aprox Trigonometricas
	  	    const float B = 4.0 / pi;
			const float C = -4.0 / (pi * pi);
			const float P = 0.225;
			float y;
			y = B * phi[i] + C * phi[i] * fabs(phi[i]);
			xi2[i] = rad[i] * ( P * (y * fabs(y) - y) + y);
			t[i] = xi1[i] * xi1[i] + xi2[i] * xi2[i];
 			u[i] = 2.0f * t[i] - 1.0f;
			v[i] = xi1[i] * sqrtf((1.0f - u[i] * u[i]) / t[i]);
			w[i] = xi2[i] * sqrtf((1.0f - u[i] * u[i]) / t[i]);
			}

        }

	//####################################################################
	
	/* ##### Depositar Energia de los 73 pasos SIN ROULETTE ##### */
	
	for(j=0;j<n_roulette;j++){ //La energia que se deposita depende del paso j
	
		for (i=0;i<VectSize;i++) { //Cada foton del vector contribuye con la misma energia
		
		    heat[shell_2D[j][i]]  += (1.0f - albedo) * weight;
		    heat2[shell_2D[j][i]] += (1.0f - albedo) * (1.0f - albedo) * weight * weight;
			
			}
		
		weight *= albedo; // Actualizo energia para el proximo j
		
		}	
		    
	//####################################################################
	        

	/* #################  LOOP final de ROULETTE  ####################*/

	for(;;){ // Loop abierto de ROULETTE
	
		for (i=0;i<VectSize;i++) {
			Rand1[i] = fast_rand() / (float)32767.0 ;
			Rand2[i] = fast_rand() / (float)32767.0 ;
			Rand3[i] = fast_rand() / (float)32767.0 ;

			} 		
	
		for (i=0;i<VectSize;i++) {
			t[i] = -logf( Rand1[i] ) / (MU_A + MU_S);			
		    x[i] += t[i] * u[i];
		    y[i] += t[i] * v[i];
		    z[i] += t[i] * w[i];
		    
		    /* Calculo shell de interaccion para los VectSize fotones */
		   	shell_1D[i] = sqrtf(x[i] * x[i] + y[i] * y[i] + z[i] * z[i] + r[i] * r[i]) * shells_per_mfp;
		    if (shell_1D[i] > SHELLS - 1) {
		        shell_1D[i] = SHELLS - 1;
		    	}	
		    }
		/* Como loop general es abierto, deposito energia en cada paso  */    		    
		for (i=0;i<VectSize;i++) { 
		    heat[shell_1D[i]]  += (1.0f - albedo) * weight;
		    heat2[shell_1D[i]] += (1.0f - albedo) * (1.0f - albedo) * weight * weight;
			}		    
		    
        weight *= albedo;  

		/* Roulette para el paquete de [VectSize] fotones */
        if (weight < tol) { /* roulette */
            if (fast_rand() / (float)32767.0 > 0.1f)
                break;
           weight /= 0.1f;}

		
        /* New direction - Distribution - Vectorized*/
	  	for (i=0;i<VectSize;i++) {	
	  	    rad[i] = sqrtf(Rand2[i]);
	  	    phi[i] = pi * (2.0f * Rand3[i] - 1.0f);
	  	    
	  	    // Coseno Funcion Libreria
	  	    xi1[i] = rad[i] * cosf(phi[i]);	  	    
	  	    // Seno Implementacion in=situ de Aprox Trigonometricas
	  	    const float B = 4.0 / pi;
			const float C = -4.0 / (pi * pi);
			const float P = 0.225;
			float y;
			y = B * phi[i] + C * phi[i] * fabs(phi[i]);
			xi2[i] = rad[i] * ( P * (y * fabs(y) - y) + y);
			
	  	    t[i] = xi1[i] * xi1[i] + xi2[i] * xi2[i];
 			u[i] = 2.0f * t[i] - 1.0f;
			v[i] = xi1[i] * sqrtf((1.0f - u[i] * u[i]) / t[i]);
			w[i] = xi2[i] * sqrtf((1.0f - u[i] * u[i]) / t[i]);
			}

        }
        
        
        
}



/***
 * Main matter
 ***/

int main(void)
{
    //heading
    
    printf("# %s\n# %s\n# %s\n", t1, t2, t3);
    printf("# Scattering = %8.3f/cm\n", MU_S);
    printf("# Absorption = %8.3f/cm\n", MU_A);
    printf("# Photons    = %8d\n", PHOTONS);
	
	
    // configure RNG
    fast_srand(SEED); // Seeder para fast_rand()
    // start timer
    double start = wtime();
    // simulation
    unsigned int PHOTONS_M = (unsigned int) PHOTONS / VectSize;
    for (unsigned int i = 0; i < PHOTONS_M; ++i) {
        photon();
    }
    // stop timer
    double end = wtime();
    assert(start <= end);
    double elapsed = end - start;

    printf("# %lf seconds\n", elapsed);
    printf("# %lf K photons per second\n#\n", 1e-3 * PHOTONS / elapsed);

	
    printf("# Radius\tHeat\n");
    printf("# [microns]\t[W/cm^3]\tError\n");
    float t = 4.0f * M_PI * powf(MICRONS_PER_SHELL, 3.0f) * PHOTONS / 1e12;
//    for (unsigned int i = 0; i < SHELLS - 1; ++i) {
//        printf("%6.0f\t%12.5f\t%12.5f\n", i * (float)MICRONS_PER_SHELL,
//               heat[i] / t / (i * i + i + 1.0 / 3.0),
//               sqrt(heat2[i] - heat[i] * heat[i] / PHOTONS) / t / (i * i + i + 1.0f / 3.0f));
//    }
//    printf("# extra\t%12.5f\n", heat[SHELLS - 1] / PHOTONS);
    

    return 0;
}
