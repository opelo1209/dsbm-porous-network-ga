/* ============================================================================
 * Redes porosas 2D del DSBM (conectividad 4) - VERSION 3: POBLACIONES POR TESELA
 *
 * Idea (B. Moreno-Montiel, 2026-11-12): se construye UNA sola configuracion
 * L x L. En lugar de una poblacion de redes completas (v2), la red se parte
 * en teselas S x S y CADA TESELA tiene su propia poblacion de K versiones.
 * La funcion objetivo es local (cada sitio depende de sus 4 enlaces), asi que
 * una tesela se puede optimizar casi independientemente de las demas.
 *
 * Algoritmo por epocas:
 *   1. CORTAR: la configuracion actual se parte en T x T teselas (T = L/S),
 *      con un corrimiento (oi, oj) de la rejilla. Cada tesela recibe K
 *      versiones, todas copia de su pedazo. Todas las versiones de una tesela
 *      comparten EXACTAMENTE el mismo inventario de radios (se guardan como
 *      permutaciones de ids locales), asi que el cruce entre versiones
 *      conserva el inventario sin desplazar genes a otras teselas.
 *   2. EVOLUCIONAR (Gep generaciones): en cada generacion, para cada tesela:
 *        - cada version se evalua en el contexto de las versiones instaladas
 *          (las mejores) de las teselas vecinas (coevolucion cooperativa);
 *        - cruce PMX de sub-bloques entre parejas de versiones + seleccion
 *          elitista 2 de 4 (opcional);
 *        - mutacion guiada (intercambios sitio-sitio y enlace-enlace dentro
 *          de la tesela, acepta si no empeora; recocido opcional);
 *        - se instala en la red la mejor version.
 *      Las teselas se procesan en 4 fases (coloreo 2x2 por paridad de su
 *      fila y columna de tesela). Dos teselas del mismo color no comparten
 *      sitios, enlaces ni incidencias, asi que se procesan en paralelo y el
 *      resultado no depende del numero de hilos.
 *   3. MIGRAR: intercambios guiados de radios entre teselas vecinas sobre la
 *      configuracion ensamblada (mueve inventario entre teselas).
 *   4. DESPLAZAR la rejilla S/2 en diagonal y volver a 1, para que las
 *      fronteras no queden fijas (evita costuras).
 *
 * Con K = 1 el metodo es busqueda local restringida a cada tesela + migracion
 * + desplazamiento de rejilla (aisla el efecto de la poblacion por tesela).
 *
 * Se conservan de v2: inventario inmutable (permutaciones de M_base),
 * evaluacion local O(1), generador xoshiro256**, costo medido en llamadas a
 * errSitio(), y el mismo formato de salida (RESUMEN, traza, CSV).
 * Con la misma semilla, M_base y la configuracion inicial son las mismas que
 * en v2, asi que las corridas son comparables semilla por semilla.
 *
 * Uso:
 *   ga_v3 L mediaS mediaE semiancho [K=4] [Gmax=0] [semilla] [prefijo=red]
 *         [--tesela=S] [--epoca=G] [--migracion=f] [--sin-cruce]
 *         [--recocido=T0,Tf] [--presupuesto=EVALS] [--traza=f.csv]
 *         [--variante=nombre]
 *   --tesela=S     lado de la tesela (por defecto L/4). L/S debe ser par >= 2.
 *   --epoca=G      generaciones por epoca antes de migrar y desplazar (20).
 *   --migracion=f  intentos de migracion por epoca, en barridos de la red
 *                  (f*N de sitios y f*2N de enlaces; por defecto 0.15).
 *   --sin-cruce    sin cruce entre versiones (solo mutacion y seleccion).
 *   --recocido     como en v2: aceptacion de Metropolis con temperatura que
 *                  baja geometricamente de T0 a Tf con el presupuesto.
 * ==========================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define PROB_MUTACION 0.15
#define KMAX 64

/* ------------------------------------------------------------------------- */
/* Generador xoshiro256** (igual que v2)                                     */
/* ------------------------------------------------------------------------- */
typedef struct { uint64_t s[4]; } Rng;
static uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static void rngSembrar(Rng *r, uint64_t semilla) { for (int i = 0; i < 4; i++) r->s[i] = splitmix64(&semilla); }
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t rngNext(Rng *r) {
    uint64_t *s = r->s;
    uint64_t res = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t; s[3] = rotl(s[3], 45);
    return res;
}
static inline double rngUnif(Rng *r) { return (rngNext(r) >> 11) * (1.0 / 9007199254740992.0); }
static inline int rngMenor(Rng *r, int n) { return (int)(((unsigned __int128)rngNext(r) * (uint64_t)n) >> 64); }

/* ------------------------------------------------------------------------- */
/* Red global                                                                */
/* ------------------------------------------------------------------------- */
static int L, N;
static double *sVal, *bVal, x1, x2;
typedef struct { int *sid; int *bid; int fit; } Crom;   /* igual que v2 */

static inline int errSitio(const Crom *c, int p) {
    int i = p / L, j = p % L;
    double s = sVal[c->sid[p]];
    int e = 0;
    e += s < bVal[c->bid[p]];
    e += s < bVal[c->bid[i * L + (j + 1) % L]];
    e += s < bVal[c->bid[N + p]];
    e += s < bVal[c->bid[N + ((i + 1) % L) * L + j]];
    return e;
}
static int evaluar(const Crom *c) { int t = 0; for (int p = 0; p < N; p++) t += errSitio(c, p); return t; }
static inline void extremosEnlace(int u, int *a, int *b) {
    if (u < N) { int i = u / L, j = u % L; *a = u; *b = i * L + (j - 1 + L) % L; }
    else { int p = u - N, i = p / L, j = p % L; *a = p; *b = ((i - 1 + L) % L) * L + j; }
}

static double g_T = 0;
static inline int rechazar(int delta, Rng *r) {
    if (delta <= 0) return 0;
    if (g_T > 0 && rngUnif(r) < exp(-delta / g_T)) return 0;
    return 1;
}
static inline int intentoSitio(Crom *c, int p, int q, long *evals, Rng *r) {
    if (p == q) return 0;
    int antes = errSitio(c, p) + errSitio(c, q);
    int t = c->sid[p]; c->sid[p] = c->sid[q]; c->sid[q] = t;
    int delta = errSitio(c, p) + errSitio(c, q) - antes;
    *evals += 4;
    if (rechazar(delta, r)) { t = c->sid[p]; c->sid[p] = c->sid[q]; c->sid[q] = t; return 0; }
    c->fit += delta;
    return 1;
}
static inline int intentoEnlace(Crom *c, int u, int v, long *evals, Rng *r) {
    if (u == v) return 0;
    int s[4], ns = 0, a, b;
    extremosEnlace(u, &a, &b); s[ns++] = a; if (b != a) s[ns++] = b;
    extremosEnlace(v, &a, &b);
    int dupA = 0, dupB = 0;
    for (int k = 0; k < ns; k++) { dupA |= s[k] == a; dupB |= s[k] == b; }
    if (!dupA) s[ns++] = a;
    if (!dupB && b != a) s[ns++] = b;
    int antes = 0, despues = 0;
    for (int k = 0; k < ns; k++) antes += errSitio(c, s[k]);
    int t = c->bid[u]; c->bid[u] = c->bid[v]; c->bid[v] = t;
    for (int k = 0; k < ns; k++) despues += errSitio(c, s[k]);
    int delta = despues - antes;
    *evals += 2 * ns;
    if (rechazar(delta, r)) { t = c->bid[u]; c->bid[u] = c->bid[v]; c->bid[v] = t; return 0; }
    c->fit += delta;
    return 1;
}

/* ------------------------------------------------------------------------- */
/* Teselas                                                                   */
/* Una tesela k = ti*T + tj tiene S2 = S*S sitios. Es duena del sitio p y de  */
/* sus ranuras de enlace izquierdo (p) y superior (N+p). Ranuras locales:    */
/*   m < S2  -> enlace izquierdo del sitio local m                           */
/*   m >= S2 -> enlace superior del sitio local m - S2                       */
/* Una version es un bloque de 3*S2 ids locales: [0,S2) sitios, [S2,3S2)      */
/* enlaces. inv*[] traduce id local -> indice global de valor (en sVal/bVal). */
/* ------------------------------------------------------------------------- */
static int S, T, S2, K;
static int oi = 0, oj = 0;           /* corrimiento de la rejilla */
static int *tpos;                    /* T*T*S2: sitio global del sitio local l de la tesela k */
static int *dueno;                   /* N: tesela de cada sitio */
static int *invS, *invB;             /* T*T*S2 y T*T*2*S2 */
static int *locS, *locB;             /* N y 2N: valor global -> id local en su tesela */
static int **vers;                   /* T*T*K bloques de 3*S2 */
static int *actual;                  /* version instalada de cada tesela */

static inline int ranura(const int *pos, int m) { return m < S2 ? pos[m] : N + pos[m - S2]; }

static void cortar(const Crom *G) {
    for (int ti = 0; ti < T; ti++)
        for (int tj = 0; tj < T; tj++) {
            int k = ti * T + tj;
            int *pos = tpos + (long)k * S2, *iS = invS + (long)k * S2, *iB = invB + 2L * k * S2;
            for (int di = 0; di < S; di++)
                for (int dj = 0; dj < S; dj++) {
                    int l = di * S + dj;
                    int p = ((oi + ti * S + di) % L) * L + (oj + tj * S + dj) % L;
                    pos[l] = p; dueno[p] = k;
                }
            for (int l = 0; l < S2; l++) { iS[l] = G->sid[pos[l]]; locS[iS[l]] = l; }
            for (int m = 0; m < 2 * S2; m++) { iB[m] = G->bid[ranura(pos, m)]; locB[iB[m]] = m; }
            for (int v = 0; v < K; v++) {
                int *blk = vers[(long)k * K + v];
                for (int l = 0; l < S2; l++) blk[l] = l;
                for (int m = 0; m < 2 * S2; m++) blk[S2 + m] = m;
            }
            actual[k] = 0;
        }
}
static void instalar(Crom *G, int k, const int *blk) {
    const int *pos = tpos + (long)k * S2, *iS = invS + (long)k * S2, *iB = invB + 2L * k * S2;
    for (int l = 0; l < S2; l++) G->sid[pos[l]] = iS[blk[l]];
    for (int m = 0; m < 2 * S2; m++) G->bid[ranura(pos, m)] = iB[blk[S2 + m]];
}
static void guardar(const Crom *G, int k, int *blk) {
    const int *pos = tpos + (long)k * S2;
    for (int l = 0; l < S2; l++) blk[l] = locS[G->sid[pos[l]]];
    for (int m = 0; m < 2 * S2; m++) blk[S2 + m] = locB[G->bid[ranura(pos, m)]];
}
/* Costo de la tesela k en el contexto actual: incidencias de sus sitios mas
   las de la columna oeste y la fila norte vecinas (las unicas otras que
   dependen de sus enlaces). Incluye terminos constantes iguales para todas
   las versiones, por lo que sirve para compararlas y para el delta global. */
static int costo(const Crom *G, int k, long *evals) {
    const int *pos = tpos + (long)k * S2;
    int e = 0;
    for (int l = 0; l < S2; l++) e += errSitio(G, pos[l]);
    for (int d = 0; d < S; d++) {
        int p = pos[d * S], i = p / L, j = p % L;
        e += errSitio(G, i * L + (j - 1 + L) % L);
        p = pos[d]; i = p / L; j = p % L;
        e += errSitio(G, ((i - 1 + L) % L) * L + j);
    }
    *evals += S2 + 2 * S;
    return e;
}
static long mutarTesela(Crom *g, int k, Rng *r) {
    const int *pos = tpos + (long)k * S2;
    long ev = 0;
    for (int l = 0; l < S2; l++)
        if (rngUnif(r) < PROB_MUTACION) intentoSitio(g, pos[l], pos[rngMenor(r, S2)], &ev, r);
    for (int m = 0; m < 2 * S2; m++)
        if (rngUnif(r) < PROB_MUTACION) intentoEnlace(g, ranura(pos, m), ranura(pos, rngMenor(r, 2 * S2)), &ev, r);
    return ev;
}

/* PMX en ids locales (todas las versiones de una tesela son permutaciones del
   mismo conjunto, asi que el hijo sigue siendo permutacion). */
static void pmxLocal(int *H, const int *D, int *posm, int n, const int *idx, int nb) {
    for (int p = 0; p < n; p++) posm[H[p]] = p;
    for (int k = 0; k < nb; k++) {
        int p = idx[k], w = D[p];
        if (H[p] == w) continue;
        int q = posm[w], v = H[p];
        H[p] = w; H[q] = v; posm[w] = p; posm[v] = q;
    }
}

typedef struct { int *h1, *h2, *posm, *bl, *blB; } TrabT;

static long procesarTesela(Crom *G, int k, Rng *r, TrabT *w, int conCruce,
                           int *delta, long *hEval, long *hAcept, long *hMejora) {
    Crom g = *G;                         /* vista: mismos arreglos, fit propio */
    int **V = vers + (long)k * K;
    int E[KMAX];
    long ev = 0;
    int Eantes = costo(&g, k, &ev);
    for (int v = 0; v < K; v++) {
        if (v == actual[k]) E[v] = Eantes;
        else { instalar(&g, k, V[v]); E[v] = costo(&g, k, &ev); }
    }
    if (conCruce && K >= 2) {
        int ord[KMAX];
        for (int v = 0; v < K; v++) ord[v] = v;
        for (int i = K - 1; i >= 1; i--) { int j = rngMenor(r, i + 1); int t = ord[i]; ord[i] = ord[j]; ord[j] = t; }
        int minLado = S / 4 > 1 ? S / 4 : 1, maxLado = S / 2 > minLado ? S / 2 : minLado;
        for (int q = 0; q < K / 2; q++) {
            int a = ord[2 * q], b = ord[2 * q + 1];
            int h = minLado + rngMenor(r, maxLado - minLado + 1);
            int ww = minLado + rngMenor(r, maxLado - minLado + 1);
            int i0 = rngMenor(r, S - h + 1), j0 = rngMenor(r, S - ww + 1), nb = 0;
            for (int di = 0; di < h; di++)
                for (int dj = 0; dj < ww; dj++) {
                    int l = (i0 + di) * S + j0 + dj;
                    w->bl[nb] = l; w->blB[2 * nb] = l; w->blB[2 * nb + 1] = S2 + l; nb++;
                }
            memcpy(w->h1, V[a], 3L * S2 * sizeof(int));
            memcpy(w->h2, V[b], 3L * S2 * sizeof(int));
            pmxLocal(w->h1, V[b], w->posm, S2, w->bl, nb);
            pmxLocal(w->h1 + S2, V[b] + S2, w->posm, 2 * S2, w->blB, 2 * nb);
            pmxLocal(w->h2, V[a], w->posm, S2, w->bl, nb);
            pmxLocal(w->h2 + S2, V[a] + S2, w->posm, 2 * S2, w->blB, 2 * nb);
            instalar(&g, k, w->h1); int e1 = costo(&g, k, &ev);
            instalar(&g, k, w->h2); int e2 = costo(&g, k, &ev);
            /* Seleccion 2 de 4; en empate se prefiere a los hijos. */
            int *c[4] = { w->h1, w->h2, V[a], V[b] };
            int e[4] = { e1, e2, E[a], E[b] }, id[4] = {0, 1, 2, 3};
            for (int x = 1; x < 4; x++)
                for (int y = x; y > 0 && e[id[y]] < e[id[y - 1]]; y--) { int t = id[y]; id[y] = id[y - 1]; id[y - 1] = t; }
            *hAcept += (id[0] < 2) + (id[1] < 2);
            /* Estricto: hijo que entra y es mejor que el peor padre (excluye clones
               de padres identicos, frecuentes al inicio de cada epoca). */
            int peor = E[a] > E[b] ? E[a] : E[b];
            *hMejora += (id[0] < 2 && e[id[0]] < peor) + (id[1] < 2 && e[id[1]] < peor);
            *hEval += 2;
            V[a] = c[id[0]]; E[a] = e[id[0]];
            V[b] = c[id[1]]; E[b] = e[id[1]];
            w->h1 = c[id[2]]; w->h2 = c[id[3]];
        }
    }
    for (int v = 0; v < K; v++) {
        instalar(&g, k, V[v]);
        g.fit = E[v];
        ev += mutarTesela(&g, k, r);
        E[v] = g.fit;
        guardar(&g, k, V[v]);
    }
    int b = 0;
    for (int v = 1; v < K; v++) if (E[v] < E[b]) b = v;
    instalar(&g, k, V[b]);
    actual[k] = b;
    *delta = E[b] - Eantes;
    return ev;
}

/* Migracion: intercambios guiados entre un sitio (o enlace) y otro de una
   tesela vecina, sobre la configuracion ensamblada. */
static inline int vecina(int k, int d) {
    int ti = k / T, tj = k % T;
    if (d == 0) ti = (ti - 1 + T) % T; else if (d == 1) ti = (ti + 1) % T;
    else if (d == 2) tj = (tj - 1 + T) % T; else tj = (tj + 1) % T;
    return ti * T + tj;
}
static long migrar(Crom *G, Rng *r, double f) {
    long ev = 0, nS = (long)(f * N), nB = (long)(f * 2 * N);
    for (long it = 0; it < nS; it++) {
        int p = rngMenor(r, N), kk = vecina(dueno[p], rngMenor(r, 4));
        intentoSitio(G, p, tpos[(long)kk * S2 + rngMenor(r, S2)], &ev, r);
    }
    for (long it = 0; it < nB; it++) {
        int u = rngMenor(r, 2 * N), p = u < N ? u : u - N, kk = vecina(dueno[p], rngMenor(r, 4));
        intentoEnlace(G, u, ranura(tpos + (long)kk * S2, rngMenor(r, 2 * S2)), &ev, r);
    }
    return ev;
}

/* ------------------------------------------------------------------------- */
/* Metricas y exportacion (iguales a v2)                                     */
/* ------------------------------------------------------------------------- */
static double correlacionVecinos(const Crom *c) {
    double m = 0, v = 0, s = 0;
    for (int p = 0; p < N; p++) m += sVal[c->sid[p]];
    m /= N;
    for (int p = 0; p < N; p++) { double z = sVal[c->sid[p]] - m; v += z * z; }
    v /= N;
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++) {
            double z = sVal[c->sid[i * L + j]] - m;
            s += z * (sVal[c->sid[i * L + (j + 1) % L]] - m);
            s += z * (sVal[c->sid[((i + 1) % L) * L + j]] - m);
        }
    return s / (2.0 * N * v);
}
static int esPermutacion(const int *v, int n) {
    char *visto = calloc(n, 1); int ok = 1;
    for (int i = 0; i < n && ok; i++) { if (v[i] < 0 || v[i] >= n || visto[v[i]]) ok = 0; else visto[v[i]] = 1; }
    free(visto); return ok;
}
static int categoria(double r) { return r <= x1 ? 0 : (r <= x2 ? 1 : 2); }
static double mismaCategoria(const Crom *c) {
    long same = 0;
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++) {
            int k = categoria(sVal[c->sid[i * L + j]]);
            same += k == categoria(sVal[c->sid[i * L + (j + 1) % L]]);
            same += k == categoria(sVal[c->sid[((i + 1) % L) * L + j]]);
        }
    return (double)same / (2.0 * N);
}
static void exportar(const Crom *c, const char *prefijo) {
    char ruta[512];
    const char *nombre[3] = {"yellow", "blue", "red"};
    snprintf(ruta, sizeof ruta, "%s_colores.csv", prefijo);
    FILE *f = fopen(ruta, "w");
    if (!f) { printf("No se pudo escribir %s\n", ruta); return; }
    fprintf(f, "x,y,r_Sitio,color\n");
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++) {
            double s = sVal[c->sid[i * L + j]];
            fprintf(f, "%d,%d,%f,%s\n", i, j, s, nombre[categoria(s)]);
        }
    fclose(f);
    snprintf(ruta, sizeof ruta, "%s_completa.csv", prefijo);
    f = fopen(ruta, "w");
    if (!f) return;
    fprintf(f, "x,y,r_Sitio,r_BIzq,r_BArr,errT1\n");
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++) {
            int p = i * L + j;
            fprintf(f, "%d,%d,%f,%f,%f,%d\n", i, j, sVal[c->sid[p]], bVal[c->bid[p]], bVal[c->bid[N + p]], errSitio(c, p));
        }
    fclose(f);
    printf("Archivos %s_colores.csv y %s_completa.csv generados.\n", prefijo, prefijo);
}
static double ahora(void) {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    return (double)clock() / CLOCKS_PER_SEC;
#endif
}

int main(int argc, char **argv) {
    const char *pos[8] = {0};
    int npos = 0, conCruce = 1, Gep = 20;
    double presupuesto = 0, fMig = 0.15, T0 = 0, Tf = 0;
    const char *rutaTraza = NULL, *variante = NULL;
    S = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--sin-cruce") == 0) conCruce = 0;
        else if (strncmp(argv[i], "--presupuesto=", 14) == 0) presupuesto = atof(argv[i] + 14);
        else if (strncmp(argv[i], "--traza=", 8) == 0) rutaTraza = argv[i] + 8;
        else if (strncmp(argv[i], "--variante=", 11) == 0) variante = argv[i] + 11;
        else if (strncmp(argv[i], "--tesela=", 9) == 0) S = atoi(argv[i] + 9);
        else if (strncmp(argv[i], "--epoca=", 8) == 0) Gep = atoi(argv[i] + 8);
        else if (strncmp(argv[i], "--migracion=", 12) == 0) fMig = atof(argv[i] + 12);
        else if (strncmp(argv[i], "--recocido=", 11) == 0) {
            if (sscanf(argv[i] + 11, "%lf%*1[,/]%lf", &T0, &Tf) != 2 || T0 <= 0 || Tf <= 0 || Tf > T0) {
                fprintf(stderr, "Use --recocido=T0,Tf (o T0/Tf) con 0 < Tf <= T0\n"); return 1;
            }
        }
        else if (strncmp(argv[i], "--", 2) == 0) { fprintf(stderr, "Opcion desconocida: %s\n", argv[i]); return 1; }
        else if (npos < 8) pos[npos++] = argv[i];
    }
    if (npos < 4) {
        fprintf(stderr, "Uso: %s L mediaS mediaE semiancho [K=4] [Gmax=0] [semilla] [prefijo=red]\n"
                        "        [--tesela=S] [--epoca=G] [--migracion=f] [--sin-cruce] [--recocido=T0,Tf]\n"
                        "        [--presupuesto=EVALS] [--traza=archivo.csv] [--variante=nombre]\n", argv[0]);
        return 1;
    }
    L = atoi(pos[0]);
    double mS = atof(pos[1]), mE = atof(pos[2]), a = atof(pos[3]);
    K = npos > 4 ? atoi(pos[4]) : 4;
    long Gmax = npos > 5 ? atol(pos[5]) : 0;
    uint64_t semilla = npos > 6 ? strtoull(pos[6], NULL, 10) : (uint64_t)time(NULL);
    const char *prefijo = npos > 7 ? pos[7] : "red";
    if (S <= 0) S = L / 4;
    if (L < 4 || a <= 0 || K < 1 || K > KMAX || S < 2 || L % S || (L / S) % 2 || Gep < 1 || fMig < 0) {
        fprintf(stderr, "Parametros invalidos (se requiere L/S par >= 2, 1 <= K <= %d).\n", KMAX); return 1;
    }
    if (K % 2 && conCruce && K > 1) K++;
    if (K > KMAX) K = KMAX;
    if (K == 1) conCruce = 0;
    if (Gmax <= 0) Gmax = 2147483647L;
    if (T0 > 0 && presupuesto <= 0 && Gmax == 2147483647L) {
        fprintf(stderr, "El recocido necesita --presupuesto o Gmax finito.\n"); return 1;
    }
    if (!variante) variante = "teselas_v3";
    N = L * L; T = L / S; S2 = S * S;
    setvbuf(stdout, NULL, _IONBF, 0);

    double omega = 1.0 - fabs(mS - mE) / (2.0 * a);
    if (omega < 0) omega = 0;
    if (omega > 1) omega = 1;
    printf("v3 teselas [%s] | L=%d S=%d (%dx%d teselas) K=%d epoca=%d migracion=%.3g cruce=%s | presupuesto=%.3g semilla=%llu | Omega=%.3f\n",
           variante, L, S, T, T, K, Gep, fMig, conCruce ? "PMX" : "no", presupuesto, (unsigned long long)semilla, omega);
    if (T0 > 0) printf("RECOCIDO,T0=%g,Tf=%g (geometrico en la fraccion de presupuesto)\n", T0, Tf);
    int nHilos = 1;
#ifdef _OPENMP
    nHilos = omp_get_max_threads();
    printf("Hilos OpenMP: %d\n", nHilos);
#endif
    double t0 = ahora();

    /* ---- Inventario M_base: misma secuencia que v2 ---- */
    Rng maestro; rngSembrar(&maestro, semilla);
    sVal = malloc(N * sizeof(double));
    bVal = malloc(2 * N * sizeof(double));
    for (int p = 0; p < N; p++) {
        sVal[p] = mS - a + 2 * a * rngUnif(&maestro);
        bVal[p] = mE - a + 2 * a * rngUnif(&maestro);
        bVal[N + p] = mE - a + 2 * a * rngUnif(&maestro);
    }
    double mn = sVal[0], mx = sVal[0];
    for (int p = 1; p < N; p++) { if (sVal[p] < mn) mn = sVal[p]; if (sVal[p] > mx) mx = sVal[p]; }
    x1 = mn + (mx - mn) / 3.0; x2 = mx - (mx - mn) / 3.0;

    /* ---- Configuracion unica (identidad = individuo 0 de v2) ---- */
    Crom G; G.sid = malloc(N * sizeof(int)); G.bid = malloc(2 * N * sizeof(int));
    for (int p = 0; p < N; p++) { G.sid[p] = p; G.bid[p] = p; G.bid[N + p] = N + p; }
    G.fit = evaluar(&G);
    double evals = N;
    Crom mejor; mejor.sid = malloc(N * sizeof(int)); mejor.bid = malloc(2 * N * sizeof(int));
    memcpy(mejor.sid, G.sid, N * sizeof(int)); memcpy(mejor.bid, G.bid, 2 * N * sizeof(int)); mejor.fit = G.fit;
    printf("Optimo[0] = %d\n", G.fit);

    /* ---- Estructuras de teselas ---- */
    int nT = T * T;
    tpos = malloc((long)nT * S2 * sizeof(int));
    dueno = malloc(N * sizeof(int));
    invS = malloc((long)nT * S2 * sizeof(int));
    invB = malloc(2L * nT * S2 * sizeof(int));
    locS = malloc(N * sizeof(int));
    locB = malloc(2 * N * sizeof(int));
    actual = malloc(nT * sizeof(int));
    vers = malloc((long)nT * K * sizeof(int *));
    for (long v = 0; v < (long)nT * K; v++) vers[v] = malloc(3L * S2 * sizeof(int));
    TrabT *trab = malloc(nHilos * sizeof(TrabT));
    for (int h = 0; h < nHilos; h++) {
        trab[h].h1 = malloc(3L * S2 * sizeof(int));
        trab[h].h2 = malloc(3L * S2 * sizeof(int));
        trab[h].posm = malloc(2L * S2 * sizeof(int));
        trab[h].bl = malloc(S2 * sizeof(int));
        trab[h].blB = malloc(2L * S2 * sizeof(int));
    }
    uint64_t *semillas = malloc(nT * sizeof(uint64_t));
    int *lista = malloc(nT * sizeof(int));
    cortar(&G);

    FILE *traza = rutaTraza ? fopen(rutaTraza, "w") : NULL;
    if (rutaTraza && !traza) fprintf(stderr, "No se pudo abrir %s\n", rutaTraza);
    if (traza) {
        fprintf(traza, "generacion,evaluaciones,tiempo_s,mejor_error\n");
        fprintf(traza, "0,%.0f,%.4f,%d\n", evals, ahora() - t0, mejor.fit);
    }

    long g = 0, epocas = 0, hijosEvaluados = 0, hijosAceptados = 0, hijosMejoran = 0;
    int ge = 0, fitG = G.fit, falloDelta = 0;
    while (g < Gmax && mejor.fit > 0 && (presupuesto <= 0 || evals < presupuesto)) {
        if (T0 > 0) {
            double frac = presupuesto > 0 ? evals / presupuesto : (double)g / Gmax;
            if (frac > 1) frac = 1;
            g_T = T0 * pow(Tf / T0, frac);
        }
        for (int k = 0; k < nT; k++) semillas[k] = rngNext(&maestro);
        long ev = 0;
        for (int color = 0; color < 4; color++) {
            int nl = 0;
            for (int k = 0; k < nT; k++)
                if (((k / T) % 2) * 2 + (k % T) % 2 == color) lista[nl++] = k;
            long evC = 0, he = 0, ha = 0, hm = 0; int dsum = 0;
            #pragma omp parallel for schedule(dynamic) reduction(+:evC,he,ha,hm,dsum)
            for (int x = 0; x < nl; x++) {
                int tid = 0;
#ifdef _OPENMP
                tid = omp_get_thread_num();
#endif
                int k = lista[x], d = 0;
                Rng r; rngSembrar(&r, semillas[k]);
                evC += procesarTesela(&G, k, &r, &trab[tid], conCruce, &d, &he, &ha, &hm);
                dsum += d;
            }
            ev += evC; hijosEvaluados += he; hijosAceptados += ha; hijosMejoran += hm; fitG += dsum;
        }
        evals += (double)ev;
        g++; ge++;
        if (fitG < mejor.fit) {
            memcpy(mejor.sid, G.sid, N * sizeof(int)); memcpy(mejor.bid, G.bid, 2 * N * sizeof(int)); mejor.fit = fitG;
            printf("Optimo[%ld] = %d\n", g, mejor.fit);
            if (traza) fprintf(traza, "%ld,%.0f,%.4f,%d\n", g, evals, ahora() - t0, mejor.fit);
        }
        if (ge >= Gep && fitG > 0) {
            /* Fin de epoca: verificar, migrar, desplazar la rejilla y volver a cortar. */
            int f = evaluar(&G); evals += N;
            if (f != fitG) { printf("ERROR: fitness acumulado %d != recalculado %d (epoca %ld)\n", fitG, f, epocas); falloDelta = 1; fitG = f; }
            G.fit = fitG;
            evals += (double)migrar(&G, &maestro, fMig);
            fitG = G.fit;
            oi = (oi + S / 2) % L; oj = (oj + S / 2) % L;
            cortar(&G);
            ge = 0; epocas++;
            if (fitG < mejor.fit) {
                memcpy(mejor.sid, G.sid, N * sizeof(int)); memcpy(mejor.bid, G.bid, 2 * N * sizeof(int)); mejor.fit = fitG;
                printf("Optimo[%ld] = %d (migracion)\n", g, mejor.fit);
                if (traza) fprintf(traza, "%ld,%.0f,%.4f,%d\n", g, evals, ahora() - t0, mejor.fit);
            }
        }
    }

    double t = ahora() - t0;
    if (traza) { fprintf(traza, "%ld,%.0f,%.4f,%d\n", g, evals, t, mejor.fit); fclose(traza); }
    if (mejor.fit == 0) printf("Solucion sin violaciones en la generacion %ld.\n", g);
    else printf("Fin. Error residual = %d.\n", mejor.fit);

    int okS = esPermutacion(mejor.sid, N), okB = esPermutacion(mejor.bid, 2 * N);
    int recalculado = evaluar(&mejor);
    double misma = mismaCategoria(&mejor), c1 = correlacionVecinos(&mejor);
    printf("--- Verificacion ---\n");
    printf("Inventario de sitios conservado:  %s\n", okS ? "SI (permutacion exacta de M_base)" : "NO");
    printf("Inventario de enlaces conservado: %s\n", okB ? "SI (permutacion exacta de M_base)" : "NO");
    printf("Fitness recalculado = %d (registrado = %d)%s\n", recalculado, mejor.fit, falloDelta ? " [hubo un error de delta]" : "");
    printf("Error residual por incidencia = %.4f%%\n", 100.0 * mejor.fit / (4.0 * N));
    printf("Vecinos de misma categoria (crudo) = %.1f%%\n", 100.0 * misma);
    printf("C(1) vecinos inmediatos = %.4f\n", c1);
    printf("Generaciones = %ld, epocas = %ld\n", g, epocas);
    printf("Evaluaciones de sitio = %.0f\n", evals);
    printf("TESELAS,S=%d,T=%d,K=%d,epoca=%d,migracion=%g,cruce=%d,hijos_evaluados=%ld,hijos_aceptados=%ld,tasa_aceptacion=%.4f,hijos_que_mejoran=%ld,tasa_estricta=%.4f\n",
           S, T, K, Gep, fMig, conCruce, hijosEvaluados, hijosAceptados,
           hijosEvaluados ? (double)hijosAceptados / hijosEvaluados : 0.0, hijosMejoran,
           hijosEvaluados ? (double)hijosMejoran / hijosEvaluados : 0.0);
    printf("Tiempo total de ejecucion: %.2f segundos\n", t);
    printf("RESUMEN,variante,L,omega,semilla,P,cruce,generaciones,evaluaciones,tiempo_s,error,error_pct,misma_cat,C1,inventario_ok,fitness_ok\n");
    printf("RESUMEN,%s,%d,%.3f,%llu,%d,%d,%ld,%.0f,%.3f,%d,%.5f,%.4f,%.4f,%d,%d\n",
           variante, L, omega, (unsigned long long)semilla, K, conCruce, g, evals, t, mejor.fit,
           100.0 * mejor.fit / (4.0 * N), misma, c1, okS && okB, recalculado == mejor.fit && !falloDelta);
    exportar(&mejor, prefijo);

    for (long v = 0; v < (long)nT * K; v++) free(vers[v]);
    for (int h = 0; h < nHilos; h++) { free(trab[h].h1); free(trab[h].h2); free(trab[h].posm); free(trab[h].bl); free(trab[h].blB); }
    free(vers); free(trab); free(tpos); free(dueno); free(invS); free(invB); free(locS); free(locB);
    free(actual); free(semillas); free(lista); free(G.sid); free(G.bid); free(mejor.sid); free(mejor.bid);
    free(sVal); free(bVal);
    return 0;
}
