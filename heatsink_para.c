#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <mpi.h>
#include <string.h>

/* Auteur : Charles Bouillaguet <charles.bouillaguet@univ-lille.fr>
   USAGE  : compiler avec -lm (et -O3 tant qu'à faire)
            rediriger la sortie standard vers un fichier texte
            gcc heatsink.c -O3 -lm -o heatsink
            ./heatsink > steady_state.txt
            [lancer le script python indiqué pour le rendu graphique]
  
   DISCLAIMER : ce code ne prétend pas à un réalisme physique absolu.
                ce code est améliorable de façon évidente (mais il a 
                été écrit de façon à rendre le plus évident possible le
                principe physique de la simulation).
*/

/* on peut changer la matière du dissipateur, sa taille, la puissance du CPU, etc. */
#define ALUMINIUM
#define NORMAL			/* MEDIUM est plus rapide, FAST est encore plus rapide (debuging) */
#define DUMP_STEADY_STATE

const double L = 0.15;		/* largeur (x) du dissipateur thermique (m) */
const double l = 0.12;		/* hauteur (y) du dissipateur thermique (m) */
const double E = 0.008;		/* épaisseur (z) du dissipateur thermique (m) */
const double watercooling_T = 20;	/* température du fluide de water-cooling, (°C) */
const double CPU_TDP = 280;	/* puissance dissipée par le CPU (W) */

#ifdef FAST
double dl = 0.004;		/* pas de simulation spatial (m) */
double dt = 0.004;		/* pas de simulation temporel (s) */
#endif

#ifdef MEDIUM
double dl = 0.002;
double dt = 0.002;
#endif

#ifdef NORMAL
double dl = 0.001;
double dt = 0.001;
#endif

#ifdef ALUMINIUM
double sink_heat_capacity = 897;	/* Capacité thermique massique du dissipateur (J / kg / K) */
double sink_density = 2710;		/* densité du dissipateur (kg / m^3) */
double sink_conductivity = 237;		/* conductivité thermique du dissipateur (W / m / K) */
double euros_per_kg = 1.594;		/* prix de la matière au kilo */
#endif

#ifdef COPPER
double sink_heat_capacity = 385;
double sink_density = 8960;
double sink_conductivity = 390;
double euros_per_kg = 5.469;
#endif

#ifdef GOLD
double sink_heat_capacity = 128;
double sink_density = 19300;
double sink_conductivity = 317;
double euros_per_kg = 47000;
#endif

#ifdef IRON
double sink_heat_capacity = 444;
double sink_density = 7860;
double sink_conductivity = 80;
double euros_per_kg = 0.083;
#endif

const double Stefan_Boltzmann = 5.6703e-8;	/* (W / m^2 / K^4), rayonnement du corps noir */
const double heat_transfer_coefficient = 10;	/* Coefficient de convection thermique (W / m^2 / K) */
double CPU_surface;

/* renvoie True si le CPU est en contact avec le dissipateur au point (x,y).
   Ceci décrit un AMD EPYC "Rome". */
static inline bool CPU_shape(double x, double y)
{
	x -= (L - 0.0754) / 2;
	y -= (l - 0.0585) / 2;
	bool small_y_ok = (y > 0.015 && y < 0.025) || (y > 0.0337 && y < 0.0437);
	bool small_x_ok = (x > 0.0113 && x < 0.0186) || (x > 0.0193 && x < 0.0266)
	    || (x > 0.0485 && x < 0.0558) || (x > 0.0566 && x < 0.0639);
	bool big_ok = (x > 0.03 && x < 0.045 && y > 0.0155 && y < 0.0435);
	return big_ok || (small_x_ok && small_y_ok);
}

/* renvoie la surface totale de contact entre le CPU et le radiateur (en m^2) */
double CPU_contact_surface()
{
	double S = 0;
	for (double x = dl / 2; x < L; x += dl)
		for (double y = dl / 2; y < l; y += dl)
			if (CPU_shape(x, y))
				S += dl * dl;
	return S;
}

/* Renvoie la nouvelle température de la cellule (i, j, k). Pour ce faire, accède aux températures
   des cellules voisines (gauche, droite, haut, bas, avant, arrière), sauf si on est sur le bord. */
static inline double update_temperature(const double *T, int u, int n, int m, int o, int i, int j, int k, int rank)
{
	/* quantité d'energie thermique qu'il faut apporter à une cellule pour la faire chauffer de 1°C */
	const double cell_heat_capacity = sink_heat_capacity * sink_density * dl * dl * dl;	/* J.K */
	const double dl2 = dl * dl;
	double thermal_flux = 0;

	if (i > 0)
		thermal_flux += (T[u - 1] - T[u]) * sink_conductivity * dl;	/* voisin x-1 */
	else {
		thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
		thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
	}

	if (i < n - 1)
		thermal_flux += (T[u + 1] - T[u]) * sink_conductivity * dl;	/* voisin x+1 */
	else {
		thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
		thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
	}
	if (j > 0)
		thermal_flux += (T[u - n] - T[u]) * sink_conductivity * dl;	/* voisin y-1 */
	else {
		/* Cellule du bas: reçoit-elle de la chaleur du CPU ? */
		if (CPU_shape(i * dl, k * dl))
			thermal_flux += CPU_TDP / CPU_surface * dl2;
		else {
			thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
			thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
		}
	}

	if (j < m - 1)
		thermal_flux += (T[u + n] - T[u]) * sink_conductivity * dl;	/* voisin y+1 */
	else {
		thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
		thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
	}

	if (k > 0){
		thermal_flux += (T[u - n * m] - T[u]) * sink_conductivity * dl;	/* voisin z-1 */
	} else {
		thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
		thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
	}

	if (k < o - 1)
		thermal_flux += (T[u + n * m] - T[u]) * sink_conductivity * dl;	/* voisin z+1 */
	else {
		thermal_flux -= Stefan_Boltzmann * dl2 * pow(T[u], 4);
		thermal_flux -= heat_transfer_coefficient * dl2 * (T[u] - watercooling_T);
	}
	/* ajuste la température avec le flux thermique */
	return T[u] + thermal_flux * dt / cell_heat_capacity;
}

/* Effectue la simulation sur le k-ème plan xy.
   v est l'indice du début du k-ème plan xy dans T et R */
static inline void do_xy_plane(const double *T, double *R, int v, int n, int m, int o, int k, int rank,int i1)
{
	if (k == 0)
		// on ne touche pas au plan z = 0 : il est maintenu à température constante par le watercooling
		return;
	for (int j = 0; j < m; j++) {	// y
		for (int i = 0; i < n; i++) {	// x
			int u = v + j * n + i;
			int u1 = (i1*n*m) + j * n + i;
			R[u1] =  update_temperature(T, u, n, m, o, i, j, k,rank);
		}
	}
}

int main(int argc, char *argv[]) {

	CPU_surface = CPU_contact_surface();
	double V_global = L * l * E;
	int n = ceil(L / dl);
	int m = ceil(E / dl);
	int o_global = ceil(l / dl);

	int rank;
	int size;
	MPI_Status status;
 
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int o_local = o_global/size;
	float epsilon = 0.1;
	//float epsilon = 200;

	if (rank == 0){
		fprintf(stderr, "DISSIPATEUR\n");
		fprintf(stderr, "\tDimension (cm) [x,y,z] = %.1f x %.1f x %.1f\n", 100 * L, 100 * E, 100 * l);
		fprintf(stderr, "\tVolume = %.1f cm^3\n", V_global * 1e6);
		fprintf(stderr, "\tMasse = %.2f kg\n", V_global * sink_density);
		fprintf(stderr, "\tPrix = %.2f €\n", V_global * sink_density * euros_per_kg);
		fprintf(stderr, "SIMULATION\n");
		fprintf(stderr, "\tGrille (x,y,z) = %d x %d x %d (%.1fMo)\n", n, m, o_global, 7.6293e-06 * n * m * o_global);
		fprintf(stderr, "\tdt = %gs\n", dt);
		fprintf(stderr, "CPU\n");
		fprintf(stderr, "\tPuissance = %.0fW\n", CPU_TDP);
		fprintf(stderr, "\tSurface = %.1f cm^2\n", CPU_surface * 10000);

		/* température de chaque cellule, en degré Kelvin. */
		double *T_global = malloc(n * m * o_global * sizeof(*T_global));
		double *R_global = malloc(n * m * o_global * sizeof(*R_global));
		if (T_global == NULL || R_global == NULL) {
			perror("Impossible d'allouer T_global et R_global");
			exit(1);
		}

		double *T_local = malloc(n * m * (o_local +1) * sizeof(*T_local));
		double *R_local = malloc(n * m * (o_local+1) * sizeof(*R_local));
		if (T_local == NULL || R_local == NULL) {
			perror("Impossible d'allouer T_global et R_global");
			exit(1);
		}

		/* initialement le radiateur est à la température du fluide de watercooling */
		for (int u = 0; u < n * m * (o_local+1); u++){
			T_local[u] = watercooling_T + 273.15;
			R_local[u] = watercooling_T + 273.15;
		}

		/* c'est parti, on allume le CPU et on simule jusqu'à avoir atteint le régime stationnaire. */
		double t = 0;
		int n_steps = 0;
		int convergence = 0;

		/* simulation des pas de temps */
		while (convergence == 0) {

			// demande les halos
			double haloLeftSd[n * m * sizeof(double)];
			int k=0;
			for (int u = n * m * (o_local-1); u < n * m * (o_local); u++){
				haloLeftSd[k] = T_local[u];
				k++;
			}
			double haloRightRcv[n * m * sizeof(double)];
			MPI_Send(&haloLeftSd, m*n*sizeof(double), MPI_DOUBLE, 1, 12 , MPI_COMM_WORLD);
			MPI_Recv(&haloRightRcv, m*n*sizeof(double), MPI_DOUBLE, 1, 99 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			k =0;
			for (int u = n*m*(o_local); u < n*m*(o_local+1); u++){
				T_local[u] = haloRightRcv[k];
				k++;
			}

			/* Met à jour toutes les cellules. On traite les plans xy par z croissant. */
			for (int k = 0; k < o_local; k++) {	// z
				int v = (k) * n * m;
				do_xy_plane(T_local, R_local, v, n, m, o_global, k, 0, k);
			}

			/* toutes les secondes, on teste la convergence et on affiche un petit compte-rendu */
			if (n_steps % ((int)(1 / dt)) == 0) {
				double delta_T = 0;
				double max = -INFINITY;
				
				for (int u = 0; u < n * m * o_local; u++) {
					delta_T += (R_local[u] - T_local[u]) * (R_local[u] - T_local[u]);
					if (R_local[u] > max)
						max = R_local[u];
				}
				double Tmax=0;
				MPI_Reduce(&max, &Tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
				max = Tmax;
				double eps;
				MPI_Allreduce(&delta_T, &eps, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				delta_T = eps;
				delta_T = sqrt(delta_T) / dt;

				fprintf(stderr, "t = %.1fs ; T_max = %.1f°C ; convergence = %g\n", t, max - 273.15, delta_T);
				if (delta_T < epsilon){
					convergence = 1;
				}
			}

			/* les nouvelles températures sont dans R */
			double * tmp = R_local;
			R_local = T_local;
			T_local = tmp;

			t += dt;
			n_steps += 1;
		}
		/* garde les valeurs sans halo */
		double *T_local1 = malloc(n * m * o_local * sizeof(*T_local1));
		for (int u = 0; u < n * m * o_local; u++) {
			memcpy(&T_local1[u], &T_local[u], sizeof(double));
		}
		MPI_Gather(T_local1, n*m*o_local, MPI_DOUBLE, T_global, n*m*o_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		#ifdef DUMP_STEADY_STATE
			printf("###### STEADY STATE; t = %.1f\n", t);
			for (int k = 0; k < o_global; k++) {	// z
				printf("# z = %g\n", k * dl);
				for (int j = 0; j < m; j++) {	// y
					for (int i = 0; i < n; i++) {	// x
						printf("%.1f ", T_global[k * n * m + j * n + i] - 273.15);
					}
					printf("\n");
				}
			}
			printf("\n");
			fprintf(stderr, "Rendu graphique : python3 rendu_picture_steady.py [filename.txt] %d %d %d\n", n, m, o_global);
		#endif

	} else{
		int nb;
		if (rank != size-1){
			nb=2;
		} else{
			nb=1;
		}
		double *T_local = malloc(n * m * (o_local+2) * sizeof(*T_local));
		double *R_local = malloc(n * m * (o_local+2) * sizeof(*R_local));
		if (T_local == NULL || R_local == NULL) {
			perror("Impossible d'allouer T_global et R_global");
			exit(1);
		}
		
		for (int u = 0; u < n * m * (o_local+2); u++){
			T_local[u] = watercooling_T + 273.15;
			R_local[u] = watercooling_T + 273.15;
		}

		/* c'est parti, on allume le CPU et on simule jusqu'à avoir atteint le régime stationnaire. */
		double t = 0;
		int n_steps = 0;
		int convergence = 0;

		/* simulation des pas de temps */
		while (convergence == 0) {

			// demande les halos
			double haloRightSd[n * m * sizeof(double)];
			int k=0;
			for (int u = n*m; u < n * m*2; u++){
				haloRightSd[k] = T_local[u];
				k++;
			}
			double haloLeftRcv[n * m * sizeof(double)];
			if (rank%2==0){
					MPI_Send(&haloRightSd, m*n*sizeof(double), MPI_DOUBLE, rank-1, 99 , MPI_COMM_WORLD);
					MPI_Recv(&haloLeftRcv, m*n*sizeof(double), MPI_DOUBLE, rank-1, 12 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}else{
					MPI_Recv(&haloLeftRcv, m*n*sizeof(double), MPI_DOUBLE, rank-1, 12 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
					MPI_Send(&haloRightSd, m*n*sizeof(double), MPI_DOUBLE, rank-1, 99 , MPI_COMM_WORLD);
			}
			for (int u = 0; u < n * m ; u++){
				T_local[u] = haloLeftRcv[u];
			}
			if (rank != size-1){ // si je suis pas le dernier processeur
				// demande les halos
				double haloLeftSd[n * m * sizeof(double)];
				int k=0;
				for (int u = n * m * (o_local); u < n * m * (o_local+1); u++){
					haloLeftSd[k] = T_local[u];
					k++;
				}
				double haloRightRcv[n * m * sizeof(double)];
				if (rank%2==0){
					MPI_Send(&haloLeftSd, m*n*sizeof(double), MPI_DOUBLE, rank+1, 12 , MPI_COMM_WORLD);
					MPI_Recv(&haloRightRcv, m*n*sizeof(double), MPI_DOUBLE, rank+1, 99 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				}else{
					MPI_Recv(&haloRightRcv, m*n*sizeof(double), MPI_DOUBLE, rank+1, 99 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					MPI_Send(&haloLeftSd, m*n*sizeof(double), MPI_DOUBLE, rank+1, 12 , MPI_COMM_WORLD);
				}
				k =0;
				for (int u = n*m*(o_local+1); u<n*m*(o_local+2); u++){
					T_local[u] = haloRightRcv[k];
					k++;
				}
			}

			/* Met à jour toutes les cellules. On traite les plans xy par z croissant. */
			int i=0;
			for (int k = rank*o_local; k < (rank*o_local)+o_local; k++) {	// z
				int v = (i+1) * n * m;
				do_xy_plane(T_local, R_local, v, n, m, o_global, k,rank,i+1);
				i++;
			}

			/* toutes les secondes, on teste la convergence et on affiche un petit compte-rendu */
			if (n_steps % ((int)(1 / dt)) == 0) {
				double delta_T = 0;
				double max = -INFINITY;
				int u1 = n*m;
				for (int u = 0; u < n * m * o_local; u++) {
					delta_T += (R_local[u1] - T_local[u1]) * (R_local[u1] - T_local[u1]);
					u1+=1;
					if (R_local[u1] > max){
						max = R_local[u1];
					}
				}
				MPI_Reduce(&max, NULL, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
				double eps;
				MPI_Allreduce(&delta_T, &eps, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
				delta_T = eps;
				delta_T = sqrt(delta_T) / dt;
				if (delta_T < epsilon)
					convergence = 1;
			}
			
			/* les nouvelles températures sont dans R */
			double * tmp = R_local;
			R_local = T_local;
			T_local = tmp;

			t += dt;
			n_steps += 1;
		}

		/* garde les valeurs sans halo */
		double *T_local1 = malloc(n * m * o_local * sizeof(*T_local1));
		int u1 = 0;
		for (int u = n*m; u < n * m * (o_local+1); u++) {
			memcpy(&T_local1[u1], &T_local[u], sizeof(double));
			u1++;
		}
		MPI_Gather(T_local1, n*m*o_local, MPI_DOUBLE, NULL, n*m*o_local, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	exit(EXIT_SUCCESS);
}
