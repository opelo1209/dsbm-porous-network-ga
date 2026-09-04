#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>

/* ============================================================================
 * Pure Monte Carlo construction of a 2D DSBM porous network (coordination
 * number four), used as the classical baseline the paper contrasts the
 * genetic algorithm against (Sec. 1.2/2.2): "Early reconstruction algorithms
 * relied on Monte Carlo procedures, where lattice elements were repeatedly
 * exchanged through random site-site and bond-bond permutations until the
 * Construction Principle was satisfied [8]. [...] this stochastic strategy
 * behaves as a blind search process, requiring a very large number of random
 * exchanges before feasible lattice configurations can be obtained."
 *
 * Unlike ConstructorRedes2D_C4_Genetico_Final.c, this file has:
 *   - No population, no crossover, no mutation, no size-category restriction
 *     on which elements may be exchanged (a move can pair a small site with a
 *     large site, etc.) - moves are chosen uniformly at random over the whole
 *     lattice, which is what makes this "pure" Monte Carlo rather than the
 *     constraint-preserving GA.
 *   - No relaxation fallback: convergence relies exclusively on random
 *     site-site and bond-bond exchanges, exactly as historically described.
 *     As a direct consequence, and unlike the GA version, this program is NOT
 *     guaranteed to reach zero error within a practical time budget for large
 *     lattices; a maximum-attempts safety cap is provided purely so the
 *     program terminates (see "Proporciona el numero maximo de intentos"
 *     below), not as an algorithmic correction.
 *
 * Acceptance rule: a candidate exchange is accepted only if it does not
 * increase the total number of Type 1 Construction Principle violations
 * (Eq. 2-3 of the paper); otherwise it is reverted. The exchanged elements
 * themselves are still chosen completely at random - this is the standard
 * way Monte Carlo reconstruction is implemented in the porous-media
 * literature this paper cites (greedy/random-walk hybrid, not simulated
 * annealing: there is no temperature schedule and no acceptance of
 * error-increasing moves).
 *
 * Efficiency note: evaluating the FULL lattice fitness after every single
 * exchange (as the naive description might suggest) would cost O(L^2) per
 * attempt, which is intractable for the L values used in the paper (up to
 * 500, i.e. 250,000 nodes) over the millions of attempts blind search
 * typically needs. Instead, only the (at most four) lattice positions whose
 * Type 1 error count can possibly change as a result of a given exchange are
 * recomputed - see "Incremental (local) fitness evaluation" below - which
 * matches the same per-node accounting used by evaluarFitness() exactly, so
 * results are directly comparable to the genetic algorithm's reported error.
 * ==========================================================================*/

// Same node layout as ConstructorRedes2D_C4_Genetico_Final.c: a gene stores a
// site radius plus its LEFT and UPPER bond radii only; right/lower bonds are
// inferred from neighboring nodes via periodic (modular) indexing.
typedef struct Nodo_Red {
    double r_Sitio;     // Site radius, r_Site
    double r_EIzq;      // Left bond radius, r_BLeft
    double r_EArr;      // Upper bond radius, r_BUp
    int errT1;          // Number of Type 1 Construction Principle violations at this node
    int errG;           // Legacy geometric-error counter from 3D DSBM models; unused here (2D, Type 1 only)
    int tipo;            // Site size category: 0 Small, 1 Medium, 2 Large - informational only,
                          // kept so the exported CSV uses the same small/medium/large color coding
                          // as the genetic algorithm's output (see exportarRedConColores). Pure
                          // Monte Carlo moves are NOT restricted by category.
    int tipoEIzq;
    int tipoEArr;
} NODO_BSM;

// Site-size (x1,x2) and bond-size (e1,e2) category thresholds, computed once
// from the randomly generated network, purely for the CSV color export -
// they play no role in how moves are selected or accepted.
double x1, x2;
double e1, e2;

double randomUniform(double, double);
void generarRedBase(NODO_BSM **, int, double, double, double);
int calcularErrT1Local(NODO_BSM **, int, int, int);
int evaluarFitness(NODO_BSM **, int);
static int intentarIntercambioSitios(NODO_BSM **, int, int *);
static int intentarIntercambioEnlaces(NODO_BSM **, int, int *);
void exportarRedConColores(NODO_BSM **, int, const char *);

int main(void) {
    int L;
    double mediaS, mediaE, desviacion;
    long maxIntentos;

    printf("=== Construccion por Monte Carlo puro (baseline, Sec. 2.2 del articulo) ===\n");
    printf("Proporciona el tamano de la red: ");
    scanf("%d", &L);
    printf("Proporciona la media de los Sitios: ");
    scanf("%lf", &mediaS);
    printf("Proporciona la media de los Enlaces < media de los Sitios: ");
    scanf("%lf", &mediaE);
    printf("Proporciona el valor de la desviacion: ");
    scanf("%lf", &desviacion);
    printf("Proporciona el numero maximo de intentos (0 = sin limite, usar con cuidado): ");
    scanf("%ld", &maxIntentos);
    if (maxIntentos <= 0) {
        maxIntentos = LONG_MAX;
    }

    srand((unsigned int)time(NULL));

    clock_t start = clock();

    NODO_BSM **RED = malloc(L * sizeof(NODO_BSM *));
    for (int i = 0; i < L; i++) {
        RED[i] = malloc(L * sizeof(NODO_BSM));
    }

    // Baseline Synthesis: same random generation + small/medium/large
    // categorization used as the GA's Mbase (Sec. 3.2, step 1). The pure
    // Monte Carlo search then works directly on this single lattice.
    generarRedBase(RED, L, mediaS, mediaE, desviacion);

    int totalError = evaluarFitness(RED, L);
    printf("Error inicial (Optimo[0]) = %d\n", totalError);

    int ultimoError = totalError;
    long iter = 0;
    const long progressInterval = 1000000L;

    while (totalError > 0 && iter < maxIntentos) {
        iter++;

        // "site-site and bond-bond permutations" (Sec. 1.2/2.2): each attempt
        // picks one of the two move families with equal probability, and the
        // two elements being exchanged are chosen uniformly at random over
        // the WHOLE lattice - no category restriction, unlike the GA.
        if (rand() % 2 == 0) {
            intentarIntercambioSitios(RED, L, &totalError);
        } else {
            intentarIntercambioEnlaces(RED, L, &totalError);
        }

        if (totalError < ultimoError) {
            ultimoError = totalError;
            printf("Optimo[%ld] = %d\n", iter, totalError);
        } else if (iter % progressInterval == 0) {
            printf("... intento %ld, error actual sin mejora reciente = %d\n", iter, totalError);
        }
    }

    if (totalError == 0) {
        printf("Solucion optima encontrada en el intento %ld.\n", iter);
    } else {
        printf("Maximo de intentos alcanzado (%ld) sin converger. Error final = %d.\n", iter, totalError);
        printf("Esto es esperado del Monte Carlo puro para redes grandes: es busqueda ciega,\n");
        printf("sin la garantia de convergencia que si tiene la version genetica con relajacion.\n");
    }

    exportarRedConColores(RED, L, "red_colores.csv");

    for (int i = 0; i < L; i++) {
        free(RED[i]);
    }
    free(RED);

    clock_t end = clock();
    double tiempoTotal = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Intentos totales: %ld\n", iter);
    printf("Tiempo total de ejecución: %.2f segundos\n", tiempoTotal);

    return 0;
}

// Función para generar números aleatorios uniformes en un rango dado
double randomUniform(double min, double max) {
    return min + (max - min) * ((double)rand() / RAND_MAX);
}

// ============================================================================
// Baseline Synthesis (Sec. 3.2, step 1 - shared with the genetic algorithm):
// samples site and bond radii uniformly from the prescribed distributions and
// labels every site/bond as Small/Medium/Large for CSV export purposes only.
// ============================================================================
void generarRedBase(NODO_BSM **RED, int L, double mediaS, double mediaE, double desviacion) {
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            RED[i][j].r_Sitio = randomUniform(mediaS - desviacion, mediaS + desviacion);
            RED[i][j].r_EIzq = randomUniform(mediaE - desviacion, mediaE + desviacion);
            RED[i][j].r_EArr = randomUniform(mediaE - desviacion, mediaE + desviacion);
            RED[i][j].errT1 = 0;
            RED[i][j].errG = 0;
            RED[i][j].tipo = -1;
            RED[i][j].tipoEIzq = -1;
            RED[i][j].tipoEArr = -1;
        }
    }

    int hist[3] = {0, 0, 0};
    double minSitio = RED[0][0].r_Sitio, maxSitio = RED[0][0].r_Sitio;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            if (RED[i][j].r_Sitio < minSitio) minSitio = RED[i][j].r_Sitio;
            if (RED[i][j].r_Sitio > maxSitio) maxSitio = RED[i][j].r_Sitio;
        }
    }
    double intervalo = (maxSitio - minSitio) / 3.0;
    x1 = minSitio + intervalo;
    x2 = maxSitio - intervalo;
    printf("x1 =  %f, x2 = %f\n", x1, x2);

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            if (RED[i][j].r_Sitio <= x1) { RED[i][j].tipo = 0; hist[0]++; }
            else if (RED[i][j].r_Sitio <= x2) { RED[i][j].tipo = 1; hist[1]++; }
            else { RED[i][j].tipo = 2; hist[2]++; }
        }
    }
    printf("Histograma de colores (sitios): Amarillo(Pequeños) = %d, Azul(Medianos) = %d, Rojo(Grandes) = %d\n",
           hist[0], hist[1], hist[2]);

    double minEnlace = RED[0][0].r_EIzq, maxEnlace = RED[0][0].r_EIzq;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            double izq = RED[i][j].r_EIzq, arr = RED[i][j].r_EArr;
            if (izq < minEnlace) minEnlace = izq;
            if (izq > maxEnlace) maxEnlace = izq;
            if (arr < minEnlace) minEnlace = arr;
            if (arr > maxEnlace) maxEnlace = arr;
        }
    }
    double intervaloB = (maxEnlace - minEnlace) / 3.0;
    e1 = minEnlace + intervaloB;
    e2 = maxEnlace - intervaloB;

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            double izq = RED[i][j].r_EIzq, arr = RED[i][j].r_EArr;
            RED[i][j].tipoEIzq = (izq <= e1) ? 0 : (izq <= e2) ? 1 : 2;
            RED[i][j].tipoEArr = (arr <= e1) ? 0 : (arr <= e2) ? 1 : 2;
        }
    }
}

// ============================================================================
// Per-node Type 1 violation count (Eq. 2), identical to the four checks in
// ConstructorRedes2D_C4_Genetico_Final.c's evaluarFitness(): a site's own
// left/upper bonds are checked against it, and the bonds "borrowed" from the
// right/lower neighbors (their left/upper fields) are checked against it too.
// Used both for the one-off full-lattice pass (evaluarFitness) and for the
// incremental local updates below - keeping a single shared implementation
// guarantees both are always exactly consistent with each other.
// ============================================================================
int calcularErrT1Local(NODO_BSM **RED, int L, int i, int j) {
    int e = 0;
    if (RED[i][j].r_Sitio < RED[i][j].r_EIzq) e++;
    if (RED[i][j].r_Sitio < RED[i][(j + 1) % L].r_EIzq) e++;
    if (RED[i][j].r_Sitio < RED[i][j].r_EArr) e++;
    if (RED[i][j].r_Sitio < RED[(i + 1) % L][j].r_EArr) e++;
    RED[i][j].errT1 = e;
    return e;
}

// One-off full-lattice fitness evaluation (Eq. 3), used only once at start-up
// to obtain the initial error count; O(L^2), never called again afterwards.
int evaluarFitness(NODO_BSM **RED, int L) {
    int total = 0;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            total += calcularErrT1Local(RED, L, i, j);
        }
    }
    return total;
}

// ============================================================================
// Incremental (local) fitness evaluation.
//
// r_Sitio(i,j) is read ONLY by errT1(i,j) itself (both of its own-bond checks
// compare against bonds stored at the same node), so a site-only exchange
// only ever changes errT1 at the two swapped positions - no neighbors need to
// be touched.
//
// r_EIzq(i,j) is read by errT1(i,j) (its own check) AND by errT1 of its LEFT
// neighbor, (i, j-1 mod L) (whose "right bond" check reads this same stored
// value). Symmetrically, r_EArr(i,j) is read by errT1(i,j) and by errT1 of
// its UP neighbor, (i-1 mod L, j). A bond exchange can therefore change errT1
// at up to four positions: the two swapped bond slots and their respective
// neighbors (left neighbor for a left-bond slot, up neighbor for an
// upper-bond slot) - fewer if any of those positions coincide.
// ============================================================================

static int intentarIntercambioSitios(NODO_BSM **RED, int L, int *totalError) {
    int i1 = rand() % L, j1 = rand() % L;
    int i2 = rand() % L, j2 = rand() % L;

    int before = calcularErrT1Local(RED, L, i1, j1) + calcularErrT1Local(RED, L, i2, j2);

    double temp = RED[i1][j1].r_Sitio;
    RED[i1][j1].r_Sitio = RED[i2][j2].r_Sitio;
    RED[i2][j2].r_Sitio = temp;

    int after = calcularErrT1Local(RED, L, i1, j1) + calcularErrT1Local(RED, L, i2, j2);
    int delta = after - before;

    if (delta <= 0) {
        *totalError += delta;
        return 1;
    }

    // Revert: not an improvement, so the "pure" blind search does not keep it.
    temp = RED[i1][j1].r_Sitio;
    RED[i1][j1].r_Sitio = RED[i2][j2].r_Sitio;
    RED[i2][j2].r_Sitio = temp;
    calcularErrT1Local(RED, L, i1, j1);
    calcularErrT1Local(RED, L, i2, j2);
    return 0;
}

static int intentarIntercambioEnlaces(NODO_BSM **RED, int L, int *totalError) {
    int i1 = rand() % L, j1 = rand() % L, campo1 = rand() % 2; // 0 = left bond, 1 = upper bond
    int i2 = rand() % L, j2 = rand() % L, campo2 = rand() % 2;

    int ni1 = (campo1 == 0) ? i1 : (i1 - 1 + L) % L;
    int nj1 = (campo1 == 0) ? (j1 - 1 + L) % L : j1;
    int ni2 = (campo2 == 0) ? i2 : (i2 - 1 + L) % L;
    int nj2 = (campo2 == 0) ? (j2 - 1 + L) % L : j2;

    int pi[4] = {i1, ni1, i2, ni2};
    int pj[4] = {j1, nj1, j2, nj2};
    int ui[4], uj[4], uniqueCount = 0;
    for (int k = 0; k < 4; k++) {
        int dup = 0;
        for (int m = 0; m < uniqueCount; m++) {
            if (ui[m] == pi[k] && uj[m] == pj[k]) { dup = 1; break; }
        }
        if (!dup) { ui[uniqueCount] = pi[k]; uj[uniqueCount] = pj[k]; uniqueCount++; }
    }

    int before = 0;
    for (int k = 0; k < uniqueCount; k++) before += calcularErrT1Local(RED, L, ui[k], uj[k]);

    double *ref1 = (campo1 == 0) ? &RED[i1][j1].r_EIzq : &RED[i1][j1].r_EArr;
    double *ref2 = (campo2 == 0) ? &RED[i2][j2].r_EIzq : &RED[i2][j2].r_EArr;
    double temp = *ref1;
    *ref1 = *ref2;
    *ref2 = temp;

    int after = 0;
    for (int k = 0; k < uniqueCount; k++) after += calcularErrT1Local(RED, L, ui[k], uj[k]);
    int delta = after - before;

    if (delta <= 0) {
        *totalError += delta;
        return 1;
    }

    // Revert
    temp = *ref1;
    *ref1 = *ref2;
    *ref2 = temp;
    for (int k = 0; k < uniqueCount; k++) calcularErrT1Local(RED, L, ui[k], uj[k]);
    return 0;
}

// Same CSV format as ConstructorRedes2D_C4_Genetico_Final.c (x,y,r_Sitio,color)
// so the existing analysis/plot_site_distribution.py and plot_correlation.py
// scripts work on Monte Carlo output without any changes.
void exportarRedConColores(NODO_BSM **RED, int L, const char *filename) {
    int hist[3] = {0, 0, 0};
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            if (RED[i][j].r_Sitio <= x1) hist[0]++;
            else if (RED[i][j].r_Sitio <= x2) hist[1]++;
            else hist[2]++;
        }
    }
    printf("Histograma de colores: Amarillo(Pequeños) = %d, Azul(Medianos) = %d, Rojo(Grandes) = %d\n",
           hist[0], hist[1], hist[2]);

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        printf("Error al abrir el archivo para escritura.\n");
        return;
    }

    fprintf(file, "x,y,r_Sitio,color\n");
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            const char *color;
            if (RED[i][j].r_Sitio <= x1) color = "yellow";
            else if (RED[i][j].r_Sitio <= x2) color = "blue";
            else color = "red";
            fprintf(file, "%d,%d,%f,%s\n", i, j, RED[i][j].r_Sitio, color);
        }
    }
    fclose(file);
    printf("Archivo %s generado con éxito.\n", filename);
}
