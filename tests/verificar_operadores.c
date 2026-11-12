/* ============================================================================
 * Arnes de verificacion de los operadores del GA (NO modifica el algoritmo).
 *
 * Incluye el archivo fuente original tal cual (renombrando su main) y llama
 * directamente a inicializarPoblacion / evaluarFitness / cruzarCromosomas /
 * mutarCromosoma para comprobar cuatro afirmaciones:
 *
 *   T0  RAND_MAX de la plataforma (afecta Fisher-Yates y el cruce).
 *   T1  (A1) El mapa de categorias de sitios nunca cambia: todo cromosoma
 *        final tiene exactamente el mapa de categorias de algun cromosoma
 *        inicial, y la fraccion de vecinos de misma categoria no se mueve.
 *   T2  (A4) Inventario por cromosoma: radios de sitios y de enlaces
 *        (multiconjunto ordenado) comparados contra M_base.
 *   T3  (B1) La mutacion es determinista entre llamadas: aplicada dos veces
 *        al mismo estado inicial, ¿produce el mismo resultado?
 *
 * Compilar (desde la raiz del repo):
 *   gcc -O2 -std=c99 -fopenmp -DGA_SRC='"../src/genetico/ConstructorRedes2D_C4_Genetico_Final.c"' \
 *       tests/verificar_operadores.c -o verificar_ga -lm
 * Ejecutar:
 *   ./verificar_ga [L] [P] [G] [mediaS] [mediaE] [semiancho] [semilla]
 *   (por defecto: 50 16 300 410 400 50 12345)
 * ==========================================================================*/
#ifndef GA_SRC
#define GA_SRC "../src/genetico/ConstructorRedes2D_C4_Genetico_Final.c"
#endif

#define main ga_main_original
#include GA_SRC
#undef main

#include <string.h>

static NODO_BSM **nuevaRed(int L) {
    NODO_BSM **r = malloc(L * sizeof(NODO_BSM *));
    for (int i = 0; i < L; i++) r[i] = malloc(L * sizeof(NODO_BSM));
    return r;
}
static void liberarRed(NODO_BSM **r, int L) {
    for (int i = 0; i < L; i++) free(r[i]);
    free(r);
}
static int cmpDouble(const void *a, const void *b) {
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}
/* Categoria recalculada desde r_Sitio (no se confia en el campo .tipo). */
static unsigned char catSitio(double r) {
    return (r <= x1) ? 0 : (r <= x2) ? 1 : 2;
}
static void mapaCategorias(NODO_BSM **red, int L, unsigned char *mapa) {
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++)
            mapa[i * L + j] = catSitio(red[i][j].r_Sitio);
}
static double fraccionMismaCategoria(const unsigned char *m, int L) {
    long same = 0;
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++) {
            same += m[i * L + j] == m[i * L + (j + 1) % L];
            same += m[i * L + j] == m[((i + 1) % L) * L + j];
        }
    return (double)same / (2.0 * L * L);
}
static void sitiosOrdenados(NODO_BSM **red, int L, double *v) {
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++) v[i * L + j] = red[i][j].r_Sitio;
    qsort(v, (size_t)L * L, sizeof(double), cmpDouble);
}
static void enlacesOrdenados(NODO_BSM **red, int L, double *v) {
    int N = L * L;
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++) {
            v[i * L + j] = red[i][j].r_EIzq;
            v[N + i * L + j] = red[i][j].r_EArr;
        }
    qsort(v, (size_t)2 * N, sizeof(double), cmpDouble);
}
static double media(const double *v, int n) {
    double s = 0; for (int i = 0; i < n; i++) s += v[i]; return s / n;
}

int main(int argc, char **argv) {
    int L        = argc > 1 ? atoi(argv[1]) : 50;
    int P        = argc > 2 ? atoi(argv[2]) : 16;
    int G        = argc > 3 ? atoi(argv[3]) : 300;
    double mS    = argc > 4 ? atof(argv[4]) : 410;
    double mE    = argc > 5 ? atof(argv[5]) : 400;
    double a     = argc > 6 ? atof(argv[6]) : 50;
    unsigned sem = argc > 7 ? (unsigned)strtoul(argv[7], NULL, 10) : 12345u;
    if (P % 2) P++;                       /* cruzarCromosomas necesita P par */
    int N = L * L;
    setvbuf(stdout, NULL, _IONBF, 0);
    srand(sem);

    printf("=== Verificacion de operadores ===\nFuente: %s\n", GA_SRC);
    printf("L=%d  P=%d  G=%d  mediaS=%.1f  mediaE=%.1f  semiancho=%.1f  semilla=%u\n\n",
           L, P, G, mS, mE, a, sem);

    /* ---------------- T0 ---------------- */
    printf("[T0] RAND_MAX = %d\n", RAND_MAX);
    if ((long)RAND_MAX < (long)N - 1)
        printf("     ATENCION: N-1 = %d > RAND_MAX -> Fisher-Yates con rand()%%(i+1) NO es uniforme a este L,\n"
               "     y el cruce solo puede elegir las primeras %d posiciones de cada categoria.\n",
               N - 1, RAND_MAX + 1);
    else
        printf("     OK para este L (problema aparece si L^2-1 > RAND_MAX, es decir L > %d).\n",
               (int)sqrt((double)RAND_MAX + 1.0));

    /* ---------------- Poblacion inicial ---------------- */
    NODO_BSM ***pob = malloc(P * sizeof(NODO_BSM **));
    for (int k = 0; k < P; k++) pob[k] = nuevaRed(L);
    printf("\n--- salida de inicializarPoblacion() ---\n");
    inicializarPoblacion(pob, P, L, mS, mE, a);
    printf("----------------------------------------\n");

    NODO_BSM **base = nuevaRed(L);
    for (int i = 0; i < L; i++) memcpy(base[i], pob[0][i], L * sizeof(NODO_BSM));

    double *sBase = malloc(N * sizeof(double)), *eBase = malloc(2 * N * sizeof(double));
    double *sTmp  = malloc(N * sizeof(double)), *eTmp  = malloc(2 * N * sizeof(double));
    sitiosOrdenados(base, L, sBase);
    enlacesOrdenados(base, L, eBase);
    double mediaEBase = media(eBase, 2 * N);

    unsigned char *mapasIni = malloc((size_t)P * N), *mapaTmp = malloc(N);
    double fracIni = 0;
    for (int k = 0; k < P; k++) {
        mapaCategorias(pob[k], L, mapasIni + (size_t)k * N);
        fracIni += fraccionMismaCategoria(mapasIni + (size_t)k * N, L);
    }
    fracIni /= P;

    /* ---------------- Evolucion: SOLO cruce + mutacion ----------------
     * Es el nucleo del ciclo de main() sin relajacion ni inyeccion de
     * diversidad, para aislar el efecto de los dos operadores. */
    int *fit = malloc(P * sizeof(int));
    printf("\nEvolucion (solo cruzarCromosomas + mutarCromosoma):\n");
    printf("  gen   mejor_fitness   media_enlaces_poblacion (M_base = %.3f)\n", mediaEBase);
    for (int g = 0; g <= G; g++) {
        int mejor = INT_MAX;
        for (int k = 0; k < P; k++) { fit[k] = evaluarFitness(pob[k], L); if (fit[k] < mejor) mejor = fit[k]; }
        if (g % (G / 5 > 0 ? G / 5 : 1) == 0 || g == G) {
            double me = 0;
            for (int k = 0; k < P; k++) { enlacesOrdenados(pob[k], L, eTmp); me += media(eTmp, 2 * N); }
            printf("  %4d   %13d   %.3f\n", g, mejor, me / P);
        }
        if (g == G) break;
        cruzarCromosomas(P, L, pob, fit);
        mutarCromosoma(P, L, pob);
    }

    /* ---------------- T1 ---------------- */
    int coinciden = 0; double fracFin = 0;
    for (int k = 0; k < P; k++) {
        mapaCategorias(pob[k], L, mapaTmp);
        fracFin += fraccionMismaCategoria(mapaTmp, L);
        for (int q = 0; q < P; q++)
            if (memcmp(mapaTmp, mapasIni + (size_t)q * N, N) == 0) { coinciden++; break; }
    }
    fracFin /= P;
    printf("\n[T1] Mapa de categorias de sitios (A1)\n");
    printf("     Cromosomas finales cuyo mapa es IDENTICO a algun mapa inicial: %d de %d\n", coinciden, P);
    printf("     Vecinos de misma categoria: inicial %.1f%%  ->  final %.1f%%  (aleatorio ~33.3%%)\n",
           100 * fracIni, 100 * fracFin);
    printf("     %s\n", coinciden == P
           ? "CONFIRMADO: los operadores no crean ningun arreglo nuevo de categorias de sitio."
           : "Los operadores SI generan arreglos nuevos de categorias.");

    /* ---------------- T2 ---------------- */
    int sitiosIguales = 0, enlacesIguales = 0;
    double meMin = 1e300, meMax = -1e300, meSum = 0;
    for (int k = 0; k < P; k++) {
        sitiosOrdenados(pob[k], L, sTmp);
        enlacesOrdenados(pob[k], L, eTmp);
        sitiosIguales  += memcmp(sTmp, sBase, N * sizeof(double)) == 0;
        enlacesIguales += memcmp(eTmp, eBase, 2 * N * sizeof(double)) == 0;
        double me = media(eTmp, 2 * N);
        meSum += me; if (me < meMin) meMin = me; if (me > meMax) meMax = me;
    }
    printf("\n[T2] Inventario por cromosoma vs M_base (A4)\n");
    printf("     Sitios  con multiconjunto identico a M_base: %d de %d\n", sitiosIguales, P);
    printf("     Enlaces con multiconjunto identico a M_base: %d de %d\n", enlacesIguales, P);
    printf("     Media de enlaces: M_base %.3f | poblacion final min %.3f / prom %.3f / max %.3f\n",
           mediaEBase, meMin, meSum / P, meMax);
    printf("     %s\n", enlacesIguales == P
           ? "OK: el inventario de enlaces se conserva."
           : "CONFIRMADO: el cruce NO conserva el inventario de enlaces de cada cromosoma.");

    /* ---------------- T3 ---------------- */
    NODO_BSM **A = nuevaRed(L), **B = nuevaRed(L);
    for (int i = 0; i < L; i++) { memcpy(A[i], base[i], L * sizeof(NODO_BSM)); memcpy(B[i], base[i], L * sizeof(NODO_BSM)); }
    NODO_BSM ***pa = &A, ***pb = &B;
    mutarCromosoma(1, L, pa);          /* "generacion 1" */
    mutarCromosoma(1, L, pb);          /* "generacion 2" sobre el MISMO estado de partida */
    int cambiosA = 0, iguales = 0;
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++) {
            cambiosA += A[i][j].r_Sitio != base[i][j].r_Sitio;
            iguales  += A[i][j].r_Sitio == B[i][j].r_Sitio;
        }
    mutarCromosoma(1, L, pa);          /* segunda aplicacion consecutiva sobre A */
    int regresan = 0, movidos = 0;
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++)
            if (base[i][j].r_Sitio != B[i][j].r_Sitio) { movidos++; regresan += A[i][j].r_Sitio == base[i][j].r_Sitio; }
    printf("\n[T3] Determinismo de la mutacion (B1)\n");
    printf("     Posiciones alteradas por una llamada: %d de %d\n", cambiosA, N);
    printf("     Dos llamadas desde el MISMO estado coinciden en %d de %d posiciones\n", iguales, N);
    if (iguales == N && cambiosA > 0)
        printf("     CONFIRMADO: la mutacion produce exactamente lo mismo en cada generacion\n"
               "     (la semilla de rand_r no depende de la generacion). Tras 2 llamadas seguidas,\n"
               "     %d de %d sitios movidos regresaron a su valor original.\n", regresan, movidos);
    else
        printf("     Las dos llamadas difieren: en esta plataforma rand_r ignora la semilla\n"
               "     (MinGW lo convierte en rand()); el comportamiento difiere del de Linux.\n");

    printf("\nFin.\n");
    return 0;
}
