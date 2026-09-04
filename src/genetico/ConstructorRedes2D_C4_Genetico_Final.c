// MinGW-w64 only exposes its rand_r() compatibility macro under this guard
// (must be defined before <stdlib.h> is included). On Windows it expands to
// a call to the regular rand(), whose per-thread state is already
// thread-local in the UCRT, so the decorrelated seeds computed in
// mutarCromosoma() are ignored on this platform but the parallel exchanges
// still draw independent sequences; on glibc (Linux/macOS) rand_r() is the
// real reentrant implementation and does use the seed as intended.
#define _POSIX_THREAD_SAFE_FUNCTIONS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

// Some MinGW header versions only define M_PI under _USE_MATH_DEFINES; the
// fallback below keeps randomNormal() portable without relying on that.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============================================================================
 * Constraint-preserving Genetic Algorithm for the construction of 2D porous
 * networks under the Dual Site-Bond Model (DSBM), coordination number four.
 *
 * This implementation follows:
 *   Moreno-Montiel, Cordero-Sanchez, Villegas-Cortez, "A Constraint-Preserving
 *   Genetic Algorithm for the Construction of Porous Networks Based on the
 *   Dual Site-Bond Model".
 *
 * Section map (paper -> code):
 *   Sec. 3.1  Genetic encoding / chromosome structure  -> struct NODO_BSM
 *   Sec. 3.2  Population initialization (Algorithm 1)  -> inicializarPoblacion
 *   Sec. 3.3  Fitness function (Eq. 2-3)               -> evaluarFitness
 *   Sec. 3.4  Constraint-preserving crossover           -> cruzarCromosomas /
 *                                                          cruzarPorCategoriaEnlace
 *   Sec. 3.5  Constraint-preserving mutation             -> mutarCromosoma
 *   Sec. 2.2  Legacy relaxation fallback                 -> relajaSitios
 *   Fig. 5    General evolutionary cycle                 -> main()
 * ==========================================================================*/

// A gene = one lattice node. It stores the site radius together with the
// LEFT and UPPER bond radii only; the right/lower bonds are inferred from
// neighboring nodes through periodic (modular) indexing (Sec. 3.1), so every
// bond is represented exactly once in the whole chromosome.
typedef struct Nodo_Red {
    double r_Sitio;     // Site radius, r_Site
    double r_EIzq;      // Left bond radius, r_BLeft
    double r_EArr;      // Upper bond radius, r_BUp
    int errT1;          // Number of Type 1 Construction Principle violations at this node
    int errG;           // Legacy geometric-error counter from 3D DSBM models; the present
                         // 2D formulation only optimizes Type 1 violations (Sec. 2.1), so
                         // this field is kept for struct compatibility but is not used.
    int tipo;            // Site size category: 0 Small, 1 Medium, 2 Large (type_S)
    int tipoEIzq;         // Left bond size category: 0 Small, 1 Medium, 2 Large (type_BLeft)
    int tipoEArr;         // Upper bond size category: 0 Small, 1 Medium, 2 Large (type_BUp)
} NODO_BSM;

// Global site-size (x1,x2) and bond-size (e1,e2) category thresholds, computed
// once from the baseline network M_base (Sec. 2.1, Fig. 1: small/medium/large
// classification). Every chromosome in the population shares these thresholds
// because they are all spatial rearrangements of the very same statistical
// inventory (Sec. 3.2).
double x1, x2; // Site category thresholds (small|medium, medium|large)
double e1, e2; // Bond category thresholds (small|medium, medium|large)

// Custom "large sentinel" used to seed the search for the minimum fitness
// value in main(); NOT the standard <limits.h> INT_MAX.
#define INT_MAX 50000000

void inicializarPoblacion(NODO_BSM ***, int , int , double , double , double );
double randomUniform(double , double );
double randomNormal(double , double );
int evaluarFitness(NODO_BSM **, int );
void cruzarCromosomas(int , int , NODO_BSM ***);
void mutarCromosoma(int, int , NODO_BSM *** );
void relajaSitios(int , int , NODO_BSM **);
void copiarRed(NODO_BSM **, NODO_BSM **, int );
void exportarRedConColores(NODO_BSM **, int , const char *);

int main(int argc, char **argv) {
    int L;

    printf("Proporciona el tamaño de la red: ");
    scanf("%d", &L);

    // Reserved for a future distribution-selection feature (e.g. uniform vs.
    // normal sampling); currently unused by the algorithm itself.
    char *distribucion = argv[2];
    double mediaS, mediaE, desviacion;

    printf("Proporciona la media de los Sitios (50 - 150): ");
    scanf("%lf", &mediaS);
    printf("Proporciona la media de los Enlaces < media de los Sitios: ");
    scanf("%lf", &mediaE);
    printf("Proporciona el valor de la desviación: ");
    scanf("%lf", &desviacion);

    srand(time(NULL));

    // Configuración del algoritmo genético
    int numCromosomas = 8; // Population size P (Sec. 4.2 found P = 64 to be the
                            // empirically optimal trade-off between minimum error
                            // and convergence speed on a 50x50 lattice; adjust here).
    int gen = 0; // Número de generaciones

    // Iniciar medición de tiempo
    clock_t start = clock();

    // Reserve one L x L chromosome (NODO_BSM matrix) per individual.
    NODO_BSM ***Poblacion = malloc(numCromosomas * sizeof(NODO_BSM **));
    for (int i = 0; i < numCromosomas; i++) {
        Poblacion[i] = malloc(L * sizeof(NODO_BSM *));
        for (int j = 0; j < L; j++) {
            Poblacion[i][j] = malloc(L * sizeof(NODO_BSM));
        }
    }

    // Reservar memoria para el mejor cromosoma
    NODO_BSM **RED_2D = malloc(L * sizeof(NODO_BSM *));
    for (int i = 0; i < L; i++) {
        RED_2D[i] = malloc(L * sizeof(NODO_BSM));
    }

    // Sec. 3.2 / Algorithm 1: build M_base and derive the rest of the
    // population as uniform Fisher-Yates permutations of it.
    inicializarPoblacion(Poblacion, numCromosomas, L, mediaS, mediaE, desviacion);

    // ---- Evolutionary cycle (Fig. 5): Evaluate Fitness (+Elitism) -> Converged? ->
    // Selection/Variation (Crossover/Mutation) -> repeat until f(Mp) = 0. ----
    int mejorFitness = INT_MAX;
    while (mejorFitness != 0) {

        // Evaluar fitness (Sec. 3.3, Eq. 3) and keep the best chromosome (elitism).
        int mejorCromosoma = 0;
        for (int i = 0; i < numCromosomas; i++) {
            int fitness = evaluarFitness(Poblacion[i], L);
            if (fitness < mejorFitness) {
                mejorFitness = fitness;
                mejorCromosoma = i;
                printf("Optimo[%d] = %d\n", gen, mejorFitness);
                // Copiar el mejor cromosoma a RED_2D
                for (int x = 0; x < L; x++) {
                    for (int y = 0; y < L; y++) {
                        RED_2D[x][y] = Poblacion[mejorCromosoma][x][y];
                    }
                }

            }
        }

        if (mejorFitness == 0) {
            printf("Solución óptima encontrada en la generación %d.\n", gen);
            break;
        }else if(((double)evaluarFitness(RED_2D, L) * 100) / (L * L) < 0.001){ //TEN CUIDADO CON DIVISIONES ENTERAS
            //printf("Entro aqui = %d, %.6f\n", evaluarFitness(RED_2D, L), ((double)evaluarFitness(RED_2D, L) * 100) / (L * L));
            // Fallback from the pre-GA relaxation strategy (Sec. 2.2): once the
            // population is extremely close to a feasible solution, slightly grow
            // the still-violating sites instead of relying on crossover/mutation.
            relajaSitios(numCromosomas, L, Poblacion[mejorCromosoma]);
        }else{

            // Sec. 3.4: constraint-preserving crossover (recombines bonds by category).
            cruzarCromosomas(numCromosomas, L, Poblacion);
            // Sec. 3.5: constraint-preserving mutation (relocates sites by category).
            mutarCromosoma(numCromosomas, L, Poblacion);
        }

        gen++;
    }

    exportarRedConColores(RED_2D, L, "red_colores.csv");

    // Liberar memoria
    for (int i = 0; i < numCromosomas; i++) {
        for (int j = 0; j < L; j++) {
            free(Poblacion[i][j]);
        }
        free(Poblacion[i]);
    }
    free(Poblacion);

    // Fin de medición de tiempo
    clock_t end = clock();
    double tiempoTotal = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Tiempo total de ejecución: %.2f segundos\n", tiempoTotal);

    return 0;
}

// ============================================================================
// Sec. 3.2 - Population Initialization with Statistical Conservation
// (Algorithm 1: "Population Initialization via Uniform Mapping")
//
// Individual k = 0 is the baseline network M_base: its site and bond radii are
// drawn directly from the prescribed Gaussian/uniform distributions and its
// size-category thresholds (x1, x2 for sites; e1, e2 for bonds) are computed
// once and shared by the whole population.
//
// Every other individual k = 1..P-1 is an INDEPENDENT, uniformly random
// spatial rearrangement of M_base (never of the previously built individual):
// the chromosome is flattened into a 1D array of genes, shuffled with a full
// Fisher-Yates permutation (P(pi) = 1/N!, Eq. 1), and reshaped back into the
// L x L periodic lattice. Because whole gene structures are swapped, site and
// bond radii (and their category labels) always travel together, so the
// prescribed statistical inventories F_S(R_S) and F_B(R_B) are preserved
// exactly - only the spatial arrangement changes.
// ============================================================================
void inicializarPoblacion(NODO_BSM ***Poblacion, int numCromosomas, int L, double mediaS, double mediaE, double desviacion) {
    int N = L * L;

    for (int k = 0; k < numCromosomas; k++) {
        if (k == 0) {
            // ---- Baseline Synthesis: build M_base from the prescribed distributions ----
            for (int i = 0; i < L; i++) {
                for (int j = 0; j < L; j++) {
                    Poblacion[k][i][j].r_Sitio = randomUniform(mediaS - desviacion, mediaS + desviacion);
                    //Poblacion[k][i][j].r_Sitio = randomNormal(mediaS, desviacion);
                    Poblacion[k][i][j].r_EIzq = randomUniform(mediaE - desviacion, mediaE + desviacion);
                    //Poblacion[k][i][j].r_EIzq = randomNormal(mediaE, desviacion);
                    Poblacion[k][i][j].r_EArr = randomUniform(mediaE - desviacion, mediaE + desviacion);
                    //Poblacion[k][i][j].r_EArr = randomNormal(mediaE, desviacion);
                    Poblacion[k][i][j].errT1 = 0;
                    Poblacion[k][i][j].errG = 0;
                    Poblacion[k][i][j].tipo = -1;
                    Poblacion[k][i][j].tipoEIzq = -1;
                    Poblacion[k][i][j].tipoEArr = -1;
                }
            }

            // ---- Site-size categorization: Small / Medium / Large (Sec. 2.1, Fig. 1) ----
            int hist[3] = {0, 0, 0};
            // Encontrar los valores mínimo y máximo de r_Sitio
            double minSitio = Poblacion[k][0][0].r_Sitio;
            double maxSitio = Poblacion[k][0][0].r_Sitio;

            for (int i = 0; i < L; i++) {
                for (int j = 0; j < L; j++) {
                    if (Poblacion[k][i][j].r_Sitio < minSitio) {
                        minSitio = Poblacion[k][i][j].r_Sitio;
                    }
                    if (Poblacion[k][i][j].r_Sitio > maxSitio) {
                        maxSitio = Poblacion[k][i][j].r_Sitio;
                    }
                }
            }
            // Calcular los límites de los intervalos
            double intervalo = (maxSitio - minSitio) / 3.0;
            x1 = minSitio + intervalo;
            x2 = maxSitio - intervalo;
            printf("x1 =  %f, x2 = %f\n", x1, x2);
            for (int i = 0; i < L; i++) {
                for (int j = 0; j < L; j++) {
                    if (Poblacion[k][i][j].r_Sitio <= x1) {
                        Poblacion[k][i][j].tipo = 0;
                        hist[0]++;
                    } else if (Poblacion[k][i][j].r_Sitio <= x2) {
                        hist[1]++;
                        Poblacion[k][i][j].tipo = 1;
                    } else {
                        hist[2]++;
                        Poblacion[k][i][j].tipo = 2;
                    }
                }
            }

            printf("Histograma de colores (sitios): Amarillo(Pequeños) = %d, Azul(Medianos) = %d, Rojo(Grandes) = %d\n", hist[0], hist[1], hist[2]);

            // ---- Bond-size categorization: Small / Medium / Large ----
            // Required by the constraint-preserving crossover operator (Sec. 3.4),
            // which recombines only bonds sharing the same category via the
            // type_BLeft / type_BUp attributes. Left and upper bonds are both
            // drawn from the same distribution F_B(R_B), so they share thresholds.
            int histB[3] = {0, 0, 0};
            double minEnlace = Poblacion[k][0][0].r_EIzq;
            double maxEnlace = Poblacion[k][0][0].r_EIzq;

            for (int i = 0; i < L; i++) {
                for (int j = 0; j < L; j++) {
                    double izq = Poblacion[k][i][j].r_EIzq;
                    double arr = Poblacion[k][i][j].r_EArr;
                    if (izq < minEnlace) minEnlace = izq;
                    if (izq > maxEnlace) maxEnlace = izq;
                    if (arr < minEnlace) minEnlace = arr;
                    if (arr > maxEnlace) maxEnlace = arr;
                }
            }

            double intervaloB = (maxEnlace - minEnlace) / 3.0;
            e1 = minEnlace + intervaloB;
            e2 = maxEnlace - intervaloB;
            printf("e1 =  %f, e2 = %f\n", e1, e2);

            for (int i = 0; i < L; i++) {
                for (int j = 0; j < L; j++) {
                    double izq = Poblacion[k][i][j].r_EIzq;
                    double arr = Poblacion[k][i][j].r_EArr;

                    if (izq <= e1) Poblacion[k][i][j].tipoEIzq = 0;
                    else if (izq <= e2) Poblacion[k][i][j].tipoEIzq = 1;
                    else Poblacion[k][i][j].tipoEIzq = 2;

                    if (arr <= e1) Poblacion[k][i][j].tipoEArr = 0;
                    else if (arr <= e2) Poblacion[k][i][j].tipoEArr = 1;
                    else Poblacion[k][i][j].tipoEArr = 2;

                    histB[Poblacion[k][i][j].tipoEIzq]++;
                    histB[Poblacion[k][i][j].tipoEArr]++;
                }
            }

            printf("Histograma de colores (enlaces): Amarillo(Pequeños) = %d, Azul(Medianos) = %d, Rojo(Grandes) = %d\n", histB[0], histB[1], histB[2]);

        } else {
            // ---- Individuals 1..P-1: uniform Fisher-Yates permutation of M_base ----

            // Structural Replication + Linear Projection (flatten): copy M_base
            // (chromosome 0, NOT the previously generated individual) into a 1D
            // working buffer of length N = L^2.
            NODO_BSM *V = malloc(N * sizeof(NODO_BSM));
            for (int i = 0; i < L; i++) {
                for (int j = 0; j < L; j++) {
                    V[i * L + j] = Poblacion[0][i][j];
                }
            }

            // Fisher-Yates Permutation (Alg. 1): every permutation of the N genes
            // is equally likely, P(pi) = 1/N!. Complete gene structures are
            // exchanged, so r_Sitio, r_EIzq, r_EArr and their category labels
            // always move together.
            for (int i = N - 1; i >= 1; i--) {
                int j = rand() % (i + 1);
                NODO_BSM temp = V[i];
                V[i] = V[j];
                V[j] = temp;
            }

            // Spatial Mapping (reshape): restore the L x L lattice. Periodic
            // boundaries need no extra bookkeeping here because right/lower
            // bonds are always derived from neighboring nodes' left/upper
            // fields at evaluation time (Sec. 3.1), not stored redundantly.
            for (int i = 0; i < L; i++) {
                for (int j = 0; j < L; j++) {
                    Poblacion[k][i][j] = V[i * L + j];
                }
            }

            free(V);
        }
    }
}

// Función para generar números aleatorios uniformes en un rango dado
double randomUniform(double min, double max) {
    return min + (max - min) * ((double)rand() / RAND_MAX);
}

// Función para generar números aleatorios con distribución normal
double randomNormal(double media, double desviacion) {
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2); // Distribución normal estándar
    return z * desviacion + media; // Transformación a la distribución deseada
}

// ============================================================================
// Sec. 3.3 - Fitness function (Eq. 2-3): counts Type 1 Construction Principle
// violations, i.e. bonds whose radius exceeds the radius of a site they
// connect. For each node (i,j), its own left/upper bonds are checked against
// its own site, and the bonds "borrowed" from the right/lower neighbors
// (their left/upper fields, via periodic modular indexing) are also checked
// against this same site. Summed over the whole lattice, every physical bond
// therefore gets evaluated once from EACH of its two endpoint sites, so a
// single bond that violates both neighbors can contribute up to two units to
// the total (a stricter accounting than Eq. 2's per-bond OR indicator, but
// f(Mp) still reaches exactly 0 iff the lattice is fully consistent).
// ============================================================================
int evaluarFitness(NODO_BSM **RED_2D, int L) {
    int errores = 0;

    #pragma omp parallel for reduction(+:errores)
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            RED_2D[i][j].errT1 = 0;

            if (RED_2D[i][j].r_Sitio < RED_2D[i][j].r_EIzq)
                RED_2D[i][j].errT1++;

            if (RED_2D[i][j].r_Sitio < RED_2D[i][(j + 1) % L].r_EIzq)
                RED_2D[i][j].errT1++;

            if (RED_2D[i][j].r_Sitio < RED_2D[i][j].r_EArr)
                RED_2D[i][j].errT1++;

            if (RED_2D[i][j].r_Sitio < RED_2D[(i + 1) % L][j].r_EArr)
                RED_2D[i][j].errT1++;

            // Sumar errores al total
            errores += RED_2D[i][j].errT1;
        }
    }

    return errores;
}

// ----------------------------------------------------------------------------
// Small accessors used by the bond-category crossover below. campo = 0 selects
// the left bond (r_EIzq / tipoEIzq, i.e. type_BLeft); campo = 1 selects the
// upper bond (r_EArr / tipoEArr, i.e. type_BUp).
// ----------------------------------------------------------------------------
static int obtenerTipoEnlace(NODO_BSM *nodo, int campo) {
    return (campo == 0) ? nodo->tipoEIzq : nodo->tipoEArr;
}

// Exchanges only the bond radius (and its category label) stored in the
// requested slot between two nodes; site radii and site categories are left
// untouched, exactly as illustrated in Fig. 10 (only the highlighted bonds move).
static void intercambiarRadioEnlace(NODO_BSM *a, NODO_BSM *b, int campo) {
    if (campo == 0) {
        double tempR = a->r_EIzq;
        a->r_EIzq = b->r_EIzq;
        b->r_EIzq = tempR;

        int tempT = a->tipoEIzq;
        a->tipoEIzq = b->tipoEIzq;
        b->tipoEIzq = tempT;
    } else {
        double tempR = a->r_EArr;
        a->r_EArr = b->r_EArr;
        b->r_EArr = tempR;

        int tempT = a->tipoEArr;
        a->tipoEArr = b->tipoEArr;
        b->tipoEArr = tempT;
    }
}

// ============================================================================
// Sec. 3.4 - Constraint-preserving crossover, single bond field (left or
// upper). Coordinates are bucketed into Small/Medium/Large by bond category
// (type_BLeft or type_BUp) for each offspring, then bond radii are exchanged
// only between matching categories. Because donor and receiver always belong
// to the same category, the prescribed bond-size distribution F_B(R_B) stays
// exactly invariant; only the spatial arrangement of bonds changes. Exchanges
// are focused on positions that currently violate the Construction Principle
// (errT1 > 0), so recombination concentrates on defective regions of the lattice.
// ============================================================================
static void cruzarPorCategoriaEnlace(NODO_BSM **Hijo1, NODO_BSM **Hijo2, int L, int campo) {
    int N = L * L;
    int *idxCat1[3], *idxCat2[3];
    int cnt1[3] = {0, 0, 0}, cnt2[3] = {0, 0, 0};

    for (int c = 0; c < 3; c++) {
        idxCat1[c] = malloc(N * sizeof(int));
        idxCat2[c] = malloc(N * sizeof(int));
    }

    // Classify every lattice position by bond category, for both offspring.
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            int idx = i * L + j;
            int t1 = obtenerTipoEnlace(&Hijo1[i][j], campo);
            int t2 = obtenerTipoEnlace(&Hijo2[i][j], campo);
            idxCat1[t1][cnt1[t1]++] = idx;
            idxCat2[t2][cnt2[t2]++] = idx;
        }
    }

    // Exchange bond radii within each category (small/medium/large).
    for (int c = 0; c < 3; c++) {
        int numIntercambios = (cnt1[c] < cnt2[c]) ? cnt1[c] : cnt2[c];
        for (int n = 0; n < numIntercambios; n++) {
            int idx1 = idxCat1[c][rand() % cnt1[c]];
            int idx2 = idxCat2[c][rand() % cnt2[c]];
            int i1 = idx1 / L, j1 = idx1 % L;
            int i2 = idx2 / L, j2 = idx2 % L;

            if (Hijo1[i1][j1].errT1 > 0) {
                intercambiarRadioEnlace(&Hijo1[i1][j1], &Hijo2[i2][j2], campo);
            }
        }
    }

    for (int c = 0; c < 3; c++) {
        free(idxCat1[c]);
        free(idxCat2[c]);
    }
}

// ============================================================================
// Sec. 3.4 - Constraint-preserving crossover operator (top level).
// For every parent pair (ii, ii+1): offspring start as exact replicas of their
// parents; category-preserving recombination is then applied independently to
// left bonds and to upper bonds (Fig. 10). Parents and offspring are jointly
// evaluated and the two fittest chromosomes survive, implementing the elitist
// replacement strategy described in Sec. 3.4.
// ============================================================================
void cruzarCromosomas(int tamPoblacion, int L, NODO_BSM ***Poblacion) {
    // Matriz auxiliar para los nuevos hijos
    NODO_BSM ***Hijos = malloc(2 * sizeof(NODO_BSM **));
    for (int i = 0; i < 2; i++) {
        Hijos[i] = malloc(L * sizeof(NODO_BSM *));
        for (int j = 0; j < L; j++) {
            Hijos[i][j] = malloc(L * sizeof(NODO_BSM));
        }
    }

    // Scratch buffers used to snapshot the two winning chromosomes before
    // writing them back into Poblacion[ii]/Poblacion[ii+1]. This avoids a
    // data-aliasing hazard: whenever a parent itself is one of the two
    // winners, writing the first winner into Poblacion[ii] would otherwise
    // silently corrupt the value later read for the second winner if that
    // winner happens to alias Poblacion[ii].
    NODO_BSM **Temp1 = malloc(L * sizeof(NODO_BSM *));
    NODO_BSM **Temp2 = malloc(L * sizeof(NODO_BSM *));
    for (int j = 0; j < L; j++) {
        Temp1[j] = malloc(L * sizeof(NODO_BSM));
        Temp2[j] = malloc(L * sizeof(NODO_BSM));
    }

    for (int ii = 0; ii < tamPoblacion; ii += 2) {

        // Offspring are initialized as exact replicas of their parents.
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                Hijos[0][i][j] = Poblacion[ii][i][j];
                Hijos[1][i][j] = Poblacion[ii + 1][i][j];
            }
        }

        // Category-preserving recombination on left bonds, then upper bonds.
        cruzarPorCategoriaEnlace(Hijos[0], Hijos[1], L, 0); // left bonds  (type_BLeft)
        cruzarPorCategoriaEnlace(Hijos[0], Hijos[1], L, 1); // upper bonds (type_BUp)

        // Arreglo para anotar los fitness de los 4 seleccionados
        int EV_C[4];
        EV_C[0] = evaluarFitness(Poblacion[ii], L);       // Fitness del padre 1
        EV_C[1] = evaluarFitness(Poblacion[ii + 1], L);   // Fitness del padre 2
        EV_C[2] = evaluarFitness(Hijos[0], L);            // Fitness del hijo 1
        EV_C[3] = evaluarFitness(Hijos[1], L);            // Fitness del hijo 2

        // Indices de los mejores cromosomas
        int indices[4] = {0, 1, 2, 3};

        // Ordenar los índices por el fitness correspondiente
        for (int i = 0; i < 4; i++) {
            for (int j = i + 1; j < 4; j++) {
                if (EV_C[indices[i]] >= EV_C[indices[j]]) {
                    int temp = indices[i];
                    indices[i] = indices[j];
                    indices[j] = temp;
                }
            }
        }

        // Seleccionar los dos mejores cromosomas (elitist replacement, Sec. 3.4).
        NODO_BSM **candidatos[4] = { Poblacion[ii], Poblacion[ii + 1], Hijos[0], Hijos[1] };
        NODO_BSM **mejor1 = candidatos[indices[0]];
        NODO_BSM **mejor2 = candidatos[indices[1]];

        // Snapshot both winners first (aliasing-safe), then commit to the population.
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                Temp1[i][j] = mejor1[i][j];
                Temp2[i][j] = mejor2[i][j];
            }
        }
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                Poblacion[ii][i][j] = Temp1[i][j];
                Poblacion[ii + 1][i][j] = Temp2[i][j];
            }
        }
    }

    // Liberar memoria de los buffers temporales
    for (int j = 0; j < L; j++) {
        free(Temp1[j]);
        free(Temp2[j]);
    }
    free(Temp1);
    free(Temp2);

    // Liberar memoria de los hijos
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < L; j++) {
            free(Hijos[i][j]);
        }
        free(Hijos[i]);
    }
    free(Hijos);
}

// ============================================================================
// Sec. 3.5 - Constraint-preserving mutation operator.
// Mutation is a spatial relocation operator: it only exchanges the POSITIONS
// of two sites that already belong to the same size category (type_S), so the
// prescribed Gaussian site-size distribution F_S(R_S) is left exactly
// invariant; bonds are never touched here. Several random candidate pairs are
// tried per chromosome (multi-point mutation) and each candidate exchange is
// accepted independently with probability 15% (Sec. 3.5).
// ============================================================================
void mutarCromosoma(int numCromosomas, int L, NODO_BSM ***POBLACION) {
    const double probMutacion = 0.15; // Acceptance probability per candidate exchange

    #pragma omp parallel for collapse(3)
    for (int k = 0; k < numCromosomas; k++) {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                // Each (k,i,j) iteration needs its own decorrelated seed: reusing a
                // seed tied only to the thread id (and never updating it) would make
                // every iteration handled by the same thread draw the identical
                // "random" partner and the identical probability roll every time.
                // omp_get_thread_num() is only called under _OPENMP (defined by the
                // compiler exclusively when -fopenmp is passed) so that the sequential
                // build - compiled from this very same source without -fopenmp - does
                // not need to link the OpenMP runtime just for this one symbol; with a
                // single thread the k/i/j terms alone already decorrelate iterations.
#ifdef _OPENMP
                unsigned int tid = (unsigned int)omp_get_thread_num();
#else
                unsigned int tid = 0u;
#endif
                unsigned int seed = (unsigned int)(tid * 104729u
                                                    + (unsigned int)k * 131u
                                                    + (unsigned int)i * 1000003u
                                                    + (unsigned int)j * 7919u + 1u);

                int i2 = rand_r(&seed) % L;
                int j2 = rand_r(&seed) % L;
                double roll = (double)rand_r(&seed) / (double)RAND_MAX;

                if (roll < probMutacion && POBLACION[k][i2][j2].tipo == POBLACION[k][i][j].tipo) {
                    // Relocate: swap only the site radius between two same-category
                    // sites; the site category itself is unaffected by construction,
                    // since both positions already belonged to the same category.
                    double temp = POBLACION[k][i][j].r_Sitio;
                    POBLACION[k][i][j].r_Sitio = POBLACION[k][i2][j2].r_Sitio;
                    POBLACION[k][i2][j2].r_Sitio = temp;
                }
            }
        }
    }
}


// ============================================================================
// Sec. 2.2 - Legacy relaxation fallback (not part of the GA proper).
// Once the population is extremely close to a feasible solution (see the
// convergence check in main()), this slightly grows the radius of any site
// still involved in a violation (+0.1%) instead of relying on crossover and
// mutation, mirroring the relaxation-based correction strategy used by
// earlier sequential DSBM reconstruction methods before this GA was proposed.
// ============================================================================
void relajaSitios(int numCromosomas, int L, NODO_BSM **POBLACION) {
    #pragma omp parallel for collapse(3)
    for (int k = 0; k < numCromosomas; k++) {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                if (POBLACION[i][j].errT1 > 0) {
                    POBLACION[i][j].r_Sitio += (POBLACION[i][j].r_Sitio * 0.001);
                }
            }
        }
    }

}

// Función para copiar una red
void copiarRed(NODO_BSM **origen, NODO_BSM **destino, int L) {
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            destino[i][j] = origen[i][j];
        }
    }
}

// ============================================================================
// Sec. 4.5 - Exports the final lattice as a CSV with the same small/medium/
// large (yellow/blue/red) site-size color coding used throughout the paper's
// figures (Fig. 1, 4, 15-17), based on the site thresholds x1/x2.
// ============================================================================
void exportarRedConColores(NODO_BSM **RED_2D, int L, const char *filename) {
    // Encontrar los valores mínimo y máximo de r_Sitio
    int hist[3] = {0, 0, 0};

    printf("x1 =  %f, x2 = %f\n", x1, x2);

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            if (RED_2D[i][j].r_Sitio <= x1) {
                hist[0]++;
            } else if (RED_2D[i][j].r_Sitio <= x2) {
                hist[1]++;
            } else {
                hist[2]++;
            }
        }
    }

    printf("Histograma de colores: Amarillo(Pequeños) = %d, Azul(Medianos) = %d, Rojo(Grandes) = %d\n", hist[0], hist[1], hist[2]);

    // Abrir archivo para escritura
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error al abrir el archivo para escritura.\n");
        return;
    }

    // Escribir encabezado
    fprintf(file, "x,y,r_Sitio,color\n");

    // Asignar colores basados en los intervalos dinámicos
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            const char *color;
            if (RED_2D[i][j].r_Sitio <= x1) {
                color = "yellow"; // Sitios pequeños
            } else if (RED_2D[i][j].r_Sitio <= x2) {
                color = "blue"; // Sitios medianos
            } else {
                color = "red"; // Sitios grandes
            }
            fprintf(file, "%d,%d,%f,%s\n", i, j, RED_2D[i][j].r_Sitio, color);
        }
    }

    fclose(file);
    printf("Archivo %s generado con éxito.\n", filename);
}
