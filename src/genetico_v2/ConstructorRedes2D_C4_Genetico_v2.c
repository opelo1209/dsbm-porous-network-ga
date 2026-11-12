/* ============================================================================
 * GA para redes porosas 2D del Modelo Dual Sitio-Enlace (DSBM), conectividad 4
 * VERSION 2 - corrige los hallazgos A1 y A4 de la revision del repositorio.
 *
 * Diferencias respecto a src/genetico/ConstructorRedes2D_C4_Genetico_Final.c:
 *
 *  A4  Conservacion EXACTA del inventario en cada cromosoma.
 *      Cada cromosoma se representa como dos permutaciones de indices:
 *        sid[p]  -> indice del radio de sitio (en sVal[]) ubicado en la posicion p
 *        bid[u]  -> indice del radio de enlace (en bVal[]) en la ranura u
 *                   (u < N: enlace izquierdo de p=u; u >= N: enlace superior de p=u-N)
 *      Los valores sVal/bVal son los de M_base y NUNCA cambian. Todos los
 *      operadores son intercambios dentro del mismo cromosoma, de modo que
 *      cada cromosoma conserva exactamente el multiconjunto de M_base.
 *      El cruce es un PMX por bloques: el hijo 1 recibe un bloque rectangular
 *      del padre 2 (sitios y enlaces) haciendo intercambios internos, y
 *      viceversa. Ya no se intercambian valores ENTRE cromosomas.
 *
 *  A1  La mutacion ya no exige "misma categoria": puede mover cualquier
 *      sitio (y cualquier enlace) a cualquier posicion del mismo cromosoma,
 *      por lo que el mapa de categorias puede reorganizarse. Por defecto
 *      (MUTACION_GUIADA=1) una propuesta se acepta solo si no aumenta el
 *      numero de violaciones (evaluacion local O(1)); con -DMUTACION_GUIADA=0
 *      se acepta a ciegas con probabilidad pm, como en el original.
 *
 *  Ademas (necesario para que lo anterior funcione igual en Linux y Windows):
 *   - Generador xoshiro256** propio (sin rand()/rand_r(): B1, B2), con
 *     semilla en la linea de comandos y resultados independientes del
 *     numero de hilos.
 *   - El mejor cromosoma se rastrea correctamente y no se muta (B3, B4).
 *   - Sin relajacion de radios: los radios jamas se modifican (B5).
 *   - Sin scanf ni argv[2] (B6).
 *
 * Uso:
 *   ./ga_v2 L mediaS mediaE semiancho [P=64] [Gmax=13000] [semilla=time] [prefijo=red]
 * Salida: <prefijo>_colores.csv (x,y,r_Sitio,color: mismo formato que el
 * original) y <prefijo>_completa.csv (incluye enlaces).
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

#ifndef MUTACION_GUIADA
#define MUTACION_GUIADA 1
#endif
#define PROB_MUTACION 0.15
#define ESTANCAMIENTO 50

/* ------------------------------------------------------------------------- */
/* Generador xoshiro256** (Blackman & Vigna), sembrado con splitmix64.       */
/* ------------------------------------------------------------------------- */
typedef struct { uint64_t s[4]; } Rng;

static uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static void rngSembrar(Rng *r, uint64_t semilla) {
    for (int i = 0; i < 4; i++) r->s[i] = splitmix64(&semilla);
}
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t rngNext(Rng *r) {
    uint64_t *s = r->s;
    uint64_t res = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t; s[3] = rotl(s[3], 45);
    return res;
}
static inline double rngUnif(Rng *r) {            /* [0,1) con 53 bits */
    return (rngNext(r) >> 11) * (1.0 / 9007199254740992.0);
}
static inline int rngMenor(Rng *r, int n) {       /* entero uniforme en [0,n) */
    return (int)(((unsigned __int128)rngNext(r) * (uint64_t)n) >> 64);
}

/* ------------------------------------------------------------------------- */
/* Datos globales: inventario de M_base (inmutable)                          */
/* ------------------------------------------------------------------------- */
static int L, N;
static double *sVal;          /* N radios de sitio  */
static double *bVal;          /* 2N radios de enlace */
static double x1, x2;         /* umbrales de categoria de sitio */

typedef struct {
    int *sid;                 /* N  */
    int *bid;                 /* 2N */
    int fit;
    /* Conjunto de sitios con al menos una violacion (solo con --focal > 0).
       vl[0..nv) lista los sitios; vpos[p] es su indice en vl o -1.
       valid = 0 obliga a reconstruirlo (tras cruce, copia o reinicio). */
    int *vl, *vpos, nv, valid;
} Crom;

static Crom cromNuevo(void) {
    Crom c; c.sid = malloc(N * sizeof(int)); c.bid = malloc(2 * N * sizeof(int)); c.fit = 0;
    c.vl = malloc(N * sizeof(int)); c.vpos = malloc(N * sizeof(int)); c.nv = 0; c.valid = 0;
    return c;
}
static void cromCopiar(Crom *d, const Crom *o) {
    memcpy(d->sid, o->sid, N * sizeof(int));
    memcpy(d->bid, o->bid, 2 * N * sizeof(int));
    d->fit = o->fit;
    d->valid = 0;
}
static void cromLiberar(Crom *c) { free(c->sid); free(c->bid); free(c->vl); free(c->vpos); }

/* ------------------------------------------------------------------------- */
/* Fitness (Ec. 2-3): violaciones por incidencia sitio-enlace                */
/* ------------------------------------------------------------------------- */
static inline int errSitio(const Crom *c, int p) {
    int i = p / L, j = p % L;
    double s = sVal[c->sid[p]];
    int e = 0;
    e += s < bVal[c->bid[p]];                              /* izquierdo propio   */
    e += s < bVal[c->bid[i * L + (j + 1) % L]];            /* derecho (vecino)   */
    e += s < bVal[c->bid[N + p]];                          /* superior propio    */
    e += s < bVal[c->bid[N + ((i + 1) % L) * L + j]];      /* inferior (vecino)  */
    return e;
}
static int evaluar(const Crom *c) {
    int t = 0;
    for (int p = 0; p < N; p++) t += errSitio(c, p);
    return t;
}
/* Sitios cuyo error depende de la ranura de enlace u (sus dos extremos). */
static inline void extremosEnlace(int u, int *a, int *b) {
    if (u < N) { int i = u / L, j = u % L; *a = u; *b = i * L + (j - 1 + L) % L; }
    else { int p = u - N, i = p / L, j = p % L; *a = p; *b = ((i - 1 + L) % L) * L + j; }
}

/* ------------------------------------------------------------------------- */
/* Inicializacion: M_base + permutaciones Fisher-Yates de genes completos    */
/* ------------------------------------------------------------------------- */
static void permutacionGenes(Crom *c, Rng *r) {
    c->valid = 0;
    int *pi = malloc(N * sizeof(int));
    for (int p = 0; p < N; p++) pi[p] = p;
    for (int i = N - 1; i >= 1; i--) { int j = rngMenor(r, i + 1); int t = pi[i]; pi[i] = pi[j]; pi[j] = t; }
    for (int p = 0; p < N; p++) { c->sid[p] = pi[p]; c->bid[p] = pi[p]; c->bid[N + p] = N + pi[p]; }
    free(pi);
}
static void identidad(Crom *c) {
    c->valid = 0;
    for (int p = 0; p < N; p++) { c->sid[p] = p; c->bid[p] = p; c->bid[N + p] = N + p; }
}

/* ------------------------------------------------------------------------- */
/* Intentos de intercambio dentro de un cromosoma con evaluacion local.      */
/* Mantienen c->fit exacto y suman a *evals las llamadas a errSitio().       */
/* guiado = 1: se revierte si aumenta el error; 0: se acepta siempre.        */
/* Recocido simulado (--recocido=T0,Tf): con g_T > 0, un movimiento guiado   */
/* que empeora en delta se acepta con probabilidad exp(-delta/T) (Metropolis).*/
/* Con g_T = 0 no se consume ningun numero aleatorio extra, asi que el       */
/* comportamiento por defecto es identico bit a bit al de la evidencia 01.   */
/* ------------------------------------------------------------------------- */
static double g_T = 0;          /* temperatura actual (0 = sin recocido) */
static double g_focal = 0;      /* --focal: prob. de elegir el primer elemento entre sitios violados */
static inline void marcar(Crom *c, int p, int e) {
    if (e > 0) { if (c->vpos[p] < 0) { c->vpos[p] = c->nv; c->vl[c->nv++] = p; } }
    else if (c->vpos[p] >= 0) {
        int k = c->vpos[p], ult = c->vl[--c->nv];
        c->vl[k] = ult; c->vpos[ult] = k; c->vpos[p] = -1;
    }
}
static inline int rechazar(int guiado, int delta, Rng *r) {
    if (!guiado || delta <= 0) return 0;
    if (g_T > 0 && rngUnif(r) < exp(-delta / g_T)) return 0;
    return 1;
}
static inline int intentoSitio(Crom *c, int p, int q, int guiado, long *evals, Rng *r) {
    if (p == q) return 0;
    int antes = errSitio(c, p) + errSitio(c, q);
    int t = c->sid[p]; c->sid[p] = c->sid[q]; c->sid[q] = t;
    int ep = errSitio(c, p), eq = errSitio(c, q);
    int delta = ep + eq - antes;
    *evals += 4;
    if (rechazar(guiado, delta, r)) { t = c->sid[p]; c->sid[p] = c->sid[q]; c->sid[q] = t; return 0; }
    c->fit += delta;
    if (c->valid) { marcar(c, p, ep); marcar(c, q, eq); }
    return 1;
}
static inline int intentoEnlace(Crom *c, int u, int v, int guiado, long *evals, Rng *r) {
    if (u == v) return 0;
    int s[4], ns = 0, a, b;
    extremosEnlace(u, &a, &b); s[ns++] = a; if (b != a) s[ns++] = b;
    extremosEnlace(v, &a, &b);
    int dupA = 0, dupB = 0;
    for (int k = 0; k < ns; k++) { dupA |= s[k] == a; dupB |= s[k] == b; }
    if (!dupA) s[ns++] = a;
    if (!dupB && b != a) s[ns++] = b;
    int antes = 0, despues = 0, es[4];
    for (int k = 0; k < ns; k++) antes += errSitio(c, s[k]);
    int t = c->bid[u]; c->bid[u] = c->bid[v]; c->bid[v] = t;
    for (int k = 0; k < ns; k++) { es[k] = errSitio(c, s[k]); despues += es[k]; }
    int delta = despues - antes;
    *evals += 2 * ns;
    if (rechazar(guiado, delta, r)) { t = c->bid[u]; c->bid[u] = c->bid[v]; c->bid[v] = t; return 0; }
    c->fit += delta;
    if (c->valid) for (int k = 0; k < ns; k++) marcar(c, s[k], es[k]);
    return 1;
}

/* ------------------------------------------------------------------------- */
/* Cruce. Dos modos (opcion --cruce=):                                        */
/*  pmx   : PMX por bloques (evidencia 01). El hijo recibe los valores del    */
/*          bloque del donante; cada valor desplazado va a la posicion donde  */
/*          estaba el valor recibido (fuera del bloque, en cualquier lugar).  */
/*  rango : transferencia por rangos. El hijo conserva SUS valores dentro del */
/*          bloque, pero los reacomoda imitando el patron del donante: su     */
/*          k-esimo valor mas pequeno va donde el donante tiene su k-esimo    */
/*          mas pequeno. No desplaza nada fuera del bloque.                   */
/* En ambos el multiconjunto de cada hijo es exactamente el de M_base (A4).   */
/* Opcion --reparar=R: R barridos de busqueda local guiada restringida a la   */
/* franja de frontera del bloque (+ posiciones desplazadas en pmx).           */
/* ------------------------------------------------------------------------- */
static int g_modoCruce = 0;     /* 0 = pmx, 1 = rango */
static int g_reparar = 0;       /* barridos de reparacion de frontera */

typedef struct { double v; int i; } Par;
static int cmpPar(const void *a, const void *b) {
    double x = ((const Par *)a)->v, y = ((const Par *)b)->v;
    return (x > y) - (x < y);
}

/* Memoria de trabajo por hilo para el cruce. */
typedef struct {
    int *pos, *bloqS, *bloqB, *zona, *zonaB, *desp;
    unsigned char *enBloque, *marca;
    Par *pa, *pb;
} Trabajo;
static Trabajo trabajoNuevo(void) {
    Trabajo w;
    w.pos = malloc(2 * N * sizeof(int));
    w.bloqS = malloc(N * sizeof(int));
    w.bloqB = malloc(2 * N * sizeof(int));
    w.zona = malloc(N * sizeof(int));
    w.zonaB = malloc(2 * N * sizeof(int));
    w.desp = malloc(4 * N * sizeof(int));
    w.enBloque = calloc(N, 1);
    w.marca = calloc(N, 1);
    w.pa = malloc(2 * N * sizeof(Par));
    w.pb = malloc(2 * N * sizeof(Par));
    return w;
}
static void trabajoLiberar(Trabajo *w) {
    free(w->pos); free(w->bloqS); free(w->bloqB); free(w->zona); free(w->zonaB);
    free(w->desp); free(w->enBloque); free(w->marca); free(w->pa); free(w->pb);
}

/* PMX: registra en desp[] las posiciones de sitio afectadas fuera del bloque. */
static void pmxBloque(int *H, const int *D, int *pos, int n, const int *bloque, int nb, int *desp, int *nd) {
    for (int p = 0; p < n; p++) pos[H[p]] = p;
    for (int k = 0; k < nb; k++) {
        int p = bloque[k], w = D[p];
        if (H[p] == w) continue;
        int q = pos[w], v = H[p];
        H[p] = w; H[q] = v; pos[w] = p; pos[v] = q;
        desp[(*nd)++] = q < N ? q : q - N;
    }
}
/* Rango: reacomoda los valores propios de H en las ranuras del bloque segun el orden de D. */
static void rangoBloque(int *H, const int *D, const double *val, const int *ranuras, int nb, Par *pa, Par *pb) {
    for (int k = 0; k < nb; k++) {
        int u = ranuras[k];
        pa[k].v = val[H[u]]; pa[k].i = H[u];
        pb[k].v = val[D[u]]; pb[k].i = u;
    }
    qsort(pa, nb, sizeof(Par), cmpPar);
    qsort(pb, nb, sizeof(Par), cmpPar);
    for (int k = 0; k < nb; k++) H[pb[k].i] = pa[k].i;
}

/* Reparacion local de la franja de frontera (anillo de ancho 1 dentro y fuera
   del bloque) y de las posiciones desplazadas. Siempre guiada. */
static long repararZona(Crom *H, Trabajo *w, int i0, int j0, int h, int ww, const int *desp, int nd, Rng *r) {
    if (g_reparar <= 0) return 0;
    int nz = 0;
    for (int di = -1; di <= h; di++)
        for (int dj = -1; dj <= ww; dj++) {
            int borde = di <= 0 || dj <= 0 || di >= h - 1 || dj >= ww - 1;
            if (!borde) continue;
            int p = ((i0 + di + L) % L) * L + (j0 + dj + L) % L;
            if (!w->marca[p]) { w->marca[p] = 1; w->zona[nz++] = p; }
        }
    for (int k = 0; k < nd; k++)
        if (!w->marca[desp[k]]) { w->marca[desp[k]] = 1; w->zona[nz++] = desp[k]; }
    for (int k = 0; k < nz; k++) w->marca[w->zona[k]] = 0;
    int nzb = 0;
    for (int k = 0; k < nz; k++) { w->zonaB[nzb++] = w->zona[k]; w->zonaB[nzb++] = N + w->zona[k]; }
    long evals = 0;
    for (int it = 0; it < g_reparar; it++) {
        for (int k = 0; k < nz; k++) intentoSitio(H, w->zona[k], w->zona[rngMenor(r, nz)], 1, &evals, r);
        for (int k = 0; k < nzb; k++) intentoEnlace(H, w->zonaB[k], w->zonaB[rngMenor(r, nzb)], 1, &evals, r);
    }
    return evals;
}

/* Devuelve las evaluaciones de errSitio() usadas. */
static long cruzarPar(Crom *P1, Crom *P2, Crom *H1, Crom *H2, Rng *r, Trabajo *w) {
    cromCopiar(H1, P1);
    cromCopiar(H2, P2);
    /* Bloque rectangular periodico de lado entre L/8 y L/2. */
    int minLado = L / 8 > 1 ? L / 8 : 1, maxLado = L / 2 > minLado ? L / 2 : minLado;
    int h = minLado + rngMenor(r, maxLado - minLado + 1);
    int ww = minLado + rngMenor(r, maxLado - minLado + 1);
    int i0 = rngMenor(r, L), j0 = rngMenor(r, L), nb = 0;
    for (int di = 0; di < h; di++)
        for (int dj = 0; dj < ww; dj++) {
            int p = ((i0 + di) % L) * L + (j0 + dj) % L;
            w->bloqS[nb] = p; w->bloqB[2 * nb] = p; w->bloqB[2 * nb + 1] = N + p; nb++;
        }
    int nd1 = 0, nd2 = 0;
    int *desp1 = w->desp, *desp2 = w->desp + 2 * N;
    if (g_modoCruce == 0) {
        pmxBloque(H1->sid, P2->sid, w->pos, N, w->bloqS, nb, desp1, &nd1);
        pmxBloque(H1->bid, P2->bid, w->pos, 2 * N, w->bloqB, 2 * nb, desp1, &nd1);
        pmxBloque(H2->sid, P1->sid, w->pos, N, w->bloqS, nb, desp2, &nd2);
        pmxBloque(H2->bid, P1->bid, w->pos, 2 * N, w->bloqB, 2 * nb, desp2, &nd2);
    } else {
        rangoBloque(H1->sid, P2->sid, sVal, w->bloqS, nb, w->pa, w->pb);
        rangoBloque(H1->bid, P2->bid, bVal, w->bloqB, 2 * nb, w->pa, w->pb);
        rangoBloque(H2->sid, P1->sid, sVal, w->bloqS, nb, w->pa, w->pb);
        rangoBloque(H2->bid, P1->bid, bVal, w->bloqB, 2 * nb, w->pa, w->pb);
    }
    H1->fit = evaluar(H1);
    H2->fit = evaluar(H2);
    long evals = 2L * N;
    evals += repararZona(H1, w, i0, j0, h, ww, desp1, nd1, r);
    evals += repararZona(H2, w, i0, j0, h, ww, desp2, nd2, r);
    return evals;
}

/* ------------------------------------------------------------------------- */
/* Mutacion (A1): intercambio de posiciones dentro del cromosoma SIN          */
/* restriccion de categoria, para sitios y para enlaces. Mantiene c->fit     */
/* exacto mediante evaluacion local.                                         */
/* Devuelve el numero de evaluaciones de errSitio() realizadas (medida de costo). */
/* ------------------------------------------------------------------------- */
/* Mutacion focalizada (--focal=f): mismo numero de intentos que la normal en
   promedio (pm*N de sitios y pm*2N de enlaces), pero con probabilidad f el
   primer elemento sale de los sitios con violaciones: el sitio mismo, o uno
   de sus 4 enlaces. La pareja sigue siendo uniforme en toda la red. */
static long mutarFocal(Crom *c, Rng *r) {
    long evals = 0;
    if (!c->valid) {
        c->nv = 0;
        for (int p = 0; p < N; p++) { c->vpos[p] = -1; if (errSitio(c, p) > 0) { c->vpos[p] = c->nv; c->vl[c->nv++] = p; } }
        c->valid = 1;
        evals += N;
    }
    long nS = (long)(PROB_MUTACION * N + 0.5), nB = (long)(PROB_MUTACION * 2 * N + 0.5);
    for (long it = 0; it < nS && c->fit > 0; it++) {
        int p = (c->nv > 0 && rngUnif(r) < g_focal) ? c->vl[rngMenor(r, c->nv)] : rngMenor(r, N);
        intentoSitio(c, p, rngMenor(r, N), MUTACION_GUIADA, &evals, r);
    }
    for (long it = 0; it < nB && c->fit > 0; it++) {
        int u;
        if (c->nv > 0 && rngUnif(r) < g_focal) {
            int p = c->vl[rngMenor(r, c->nv)], i = p / L, j = p % L;
            switch (rngMenor(r, 4)) {
                case 0: u = p; break;                                   /* izquierdo */
                case 1: u = i * L + (j + 1) % L; break;                 /* derecho   */
                case 2: u = N + p; break;                               /* superior  */
                default: u = N + ((i + 1) % L) * L + j; break;          /* inferior  */
            }
        } else u = rngMenor(r, 2 * N);
        intentoEnlace(c, u, rngMenor(r, 2 * N), MUTACION_GUIADA, &evals, r);
    }
    return evals;
}
static long mutar(Crom *c, Rng *r) {
    if (g_focal > 0) return mutarFocal(c, r);
    long evals = 0;
    for (int p = 0; p < N; p++)
        if (rngUnif(r) < PROB_MUTACION) intentoSitio(c, p, rngMenor(r, N), MUTACION_GUIADA, &evals, r);
    for (int u = 0; u < 2 * N; u++)
        if (rngUnif(r) < PROB_MUTACION) intentoEnlace(c, u, rngMenor(r, 2 * N), MUTACION_GUIADA, &evals, r);
    return evals;
}

/* Correlacion de radios de sitio entre vecinos inmediatos (4 vecinos). */
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

/* ------------------------------------------------------------------------- */
/* Verificaciones y exportacion                                              */
/* ------------------------------------------------------------------------- */
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

/* ------------------------------------------------------------------------- */
/* Uso:                                                                      */
/*   ga_v2 L mediaS mediaE semiancho [P=64] [Gmax=13000] [semilla] [prefijo] */
/*         [--sin-cruce] [--presupuesto=EVALS] [--traza=archivo.csv]         */
/*         [--variante=nombre] [--cruce=pmx|rango] [--reparar=R]             */
/* --sin-cruce     : desactiva el cruce (solo mutacion; P puede ser 1).      */
/* --cruce=MODO    : pmx (por defecto, evidencia 01) o rango.                */
/* --reparar=R     : R barridos de busqueda local en la frontera del bloque  */
/*                   tras el cruce (0 = sin reparacion, por defecto).        */
/* --focal=f       : mutacion focalizada en sitios con violaciones (0..1).   */
/* --recocido=T0,Tf: recocido simulado en la mutacion guiada. La temperatura */
/*                   baja geometricamente de T0 a Tf segun la fraccion del   */
/*                   presupuesto usada (o de Gmax si no hay presupuesto).    */
/* --presupuesto=E : se detiene al alcanzar E evaluaciones de errSitio()     */
/*                   (costo independiente del hardware y del n. de hilos).   */
/* --traza=f.csv   : escribe generacion,evaluaciones,tiempo,mejor_error en   */
/*                   cada mejora (curva de convergencia).                    */
/* Gmax = 0 significa sin limite de generaciones.                            */
/* Al final imprime una linea "RESUMEN,..." en formato CSV.                  */
/* ------------------------------------------------------------------------- */
int main(int argc, char **argv) {
    const char *pos[8] = {0};
    int npos = 0, conCruce = 1;
    double presupuesto = 0;
    const char *rutaTraza = NULL, *variante = NULL;
    double T0 = 0, Tf = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--sin-cruce") == 0) conCruce = 0;
        else if (strncmp(argv[i], "--presupuesto=", 14) == 0) presupuesto = atof(argv[i] + 14);
        else if (strncmp(argv[i], "--traza=", 8) == 0) rutaTraza = argv[i] + 8;
        else if (strncmp(argv[i], "--variante=", 11) == 0) variante = argv[i] + 11;
        else if (strcmp(argv[i], "--cruce=pmx") == 0) g_modoCruce = 0;
        else if (strcmp(argv[i], "--cruce=rango") == 0) g_modoCruce = 1;
        else if (strncmp(argv[i], "--reparar=", 10) == 0) g_reparar = atoi(argv[i] + 10);
        else if (strncmp(argv[i], "--focal=", 8) == 0) {
            g_focal = atof(argv[i] + 8);
            if (g_focal < 0 || g_focal > 1) { fprintf(stderr, "Use --focal=f con 0 <= f <= 1\n"); return 1; }
        }
        else if (strncmp(argv[i], "--recocido=", 11) == 0) {
            if (sscanf(argv[i] + 11, "%lf%*1[,/]%lf", &T0, &Tf) != 2 || T0 <= 0 || Tf <= 0 || Tf > T0) {
                fprintf(stderr, "Use --recocido=T0,Tf (o T0/Tf) con 0 < Tf <= T0\n"); return 1;
            }
        }
        else if (strncmp(argv[i], "--", 2) == 0) { fprintf(stderr, "Opcion desconocida: %s\n", argv[i]); return 1; }
        else if (npos < 8) pos[npos++] = argv[i];
    }
    if (npos < 4) {
        fprintf(stderr, "Uso: %s L mediaS mediaE semiancho [P=64] [Gmax=13000] [semilla] [prefijo=red]\n"
                        "        [--sin-cruce] [--presupuesto=EVALS] [--traza=archivo.csv] [--variante=nombre]\n"
                        "        [--cruce=pmx|rango] [--reparar=R] [--recocido=T0,Tf] [--focal=f]\n", argv[0]);
        return 1;
    }
    L = atoi(pos[0]);
    double mS = atof(pos[1]), mE = atof(pos[2]), a = atof(pos[3]);
    int P = npos > 4 ? atoi(pos[4]) : 64;
    long Gmax = npos > 5 ? atol(pos[5]) : 13000;
    uint64_t semilla = npos > 6 ? strtoull(pos[6], NULL, 10) : (uint64_t)time(NULL);
    const char *prefijo = npos > 7 ? pos[7] : "red";
    if (!variante) variante = conCruce ? "GA" : "sin_cruce";
    if (L < 2 || P < 1 || a <= 0 || (conCruce && P < 2)) { fprintf(stderr, "Parametros invalidos.\n"); return 1; }
    if (conCruce && P % 2) P++;
    if (Gmax <= 0) Gmax = 2147483647L;
    if (T0 > 0 && presupuesto <= 0 && Gmax == 2147483647L) {
        fprintf(stderr, "El recocido necesita --presupuesto o Gmax finito para el programa de temperatura.\n"); return 1;
    }
    N = L * L;
    setvbuf(stdout, NULL, _IONBF, 0);

    double omega = 1.0 - fabs(mS - mE) / (2.0 * a);
    if (omega < 0) omega = 0;
    if (omega > 1) omega = 1;
    printf("GA v2 [%s] | L=%d P=%d Gmax=%ld presupuesto=%.3g semilla=%llu | mediaS=%.2f mediaE=%.2f semiancho=%.2f -> Omega=%.3f\n",
           variante, L, P, Gmax, presupuesto, (unsigned long long)semilla, mS, mE, a, omega);
    printf("Mutacion %s, pm=%.2f, sin restriccion de categoria; cruce %s; reparacion de frontera: %d barridos.\n",
           MUTACION_GUIADA ? "GUIADA (acepta si no empeora)" : "CIEGA", PROB_MUTACION,
           !conCruce ? "DESACTIVADO" : (g_modoCruce ? "por RANGOS en bloque" : "PMX por bloques"),
           conCruce ? g_reparar : 0);
    if (T0 > 0) printf("RECOCIDO,T0=%g,Tf=%g (geometrico en la fraccion de presupuesto)\n", T0, Tf);
    if (g_focal > 0) printf("FOCAL,f=%g (primer elemento de cada intento tomado de los sitios violados)\n", g_focal);
#ifdef _OPENMP
    printf("Hilos OpenMP: %d\n", omp_get_max_threads());
#endif
    double t0 = ahora();

    /* ---- Inventario M_base ---- */
    Rng maestro; rngSembrar(&maestro, semilla);
    sVal = malloc(N * sizeof(double));
    bVal = malloc(2 * N * sizeof(double));
    for (int p = 0; p < N; p++) {
        sVal[p] = mS - a + 2 * a * rngUnif(&maestro);
        bVal[p] = mE - a + 2 * a * rngUnif(&maestro);       /* izquierdo */
        bVal[N + p] = mE - a + 2 * a * rngUnif(&maestro);   /* superior  */
    }
    double mn = sVal[0], mx = sVal[0];
    for (int p = 1; p < N; p++) { if (sVal[p] < mn) mn = sVal[p]; if (sVal[p] > mx) mx = sVal[p]; }
    x1 = mn + (mx - mn) / 3.0; x2 = mx - (mx - mn) / 3.0;
    double medS = 0, medB = 0;
    for (int p = 0; p < N; p++) medS += sVal[p];
    for (int u = 0; u < 2 * N; u++) medB += bVal[u];
    printf("M_base: media sitios=%.3f  media enlaces=%.3f  x1=%.3f x2=%.3f\n", medS / N, medB / (2 * N), x1, x2);

    /* ---- Poblacion ---- */
    int nScr = P < 2 ? 2 : P;
    Crom *pob = malloc(P * sizeof(Crom)), *scr = malloc(nScr * sizeof(Crom));
    for (int k = 0; k < P; k++) pob[k] = cromNuevo();
    for (int k = 0; k < nScr; k++) scr[k] = cromNuevo();
    identidad(&pob[0]);
    for (int k = 1; k < P; k++) permutacionGenes(&pob[k], &maestro);
    for (int k = 0; k < P; k++) pob[k].fit = evaluar(&pob[k]);
    double evals = (double)P * N;

    Crom mejor = cromNuevo();
    int iMejor = 0;
    for (int k = 1; k < P; k++) if (pob[k].fit < pob[iMejor].fit) iMejor = k;
    cromCopiar(&mejor, &pob[iMejor]);
    printf("Optimo[0] = %d\n", mejor.fit);

    FILE *traza = rutaTraza ? fopen(rutaTraza, "w") : NULL;
    if (rutaTraza && !traza) fprintf(stderr, "No se pudo abrir %s\n", rutaTraza);
    if (traza) {
        fprintf(traza, "generacion,evaluaciones,tiempo_s,mejor_error\n");
        fprintf(traza, "0,%.0f,%.4f,%d\n", evals, ahora() - t0, mejor.fit);
    }

    int *orden = malloc(P * sizeof(int));
    uint64_t *semillas = malloc(P * sizeof(uint64_t));
    long g = 0; int sinMejora = 0;
    long hijosEvaluados = 0, hijosAceptados = 0;
    int nHilos = 1;
#ifdef _OPENMP
    nHilos = omp_get_max_threads();
#endif
    Trabajo *trab = NULL;
    if (conCruce) {
        trab = malloc(nHilos * sizeof(Trabajo));
        for (int k = 0; k < nHilos; k++) trab[k] = trabajoNuevo();
    }

    while (g < Gmax && mejor.fit > 0 && (presupuesto <= 0 || evals < presupuesto)) {
        int fitAntes = mejor.fit;
        if (T0 > 0) {
            double frac = presupuesto > 0 ? evals / presupuesto : (double)g / Gmax;
            if (frac > 1) frac = 1;
            g_T = T0 * pow(Tf / T0, frac);
        }

        /* Inyeccion de diversidad: reemplaza el peor cuarto (nunca al mejor actual). */
        if (sinMejora >= ESTANCAMIENTO && P >= 4) {
            for (int k = 0; k < P; k++) orden[k] = k;
            for (int a2 = 0; a2 < P; a2++)
                for (int b2 = a2 + 1; b2 < P; b2++)
                    if (pob[orden[b2]].fit > pob[orden[a2]].fit) { int t = orden[a2]; orden[a2] = orden[b2]; orden[b2] = t; }
            int reemplazos = 0;
            for (int r = 0; r < P && reemplazos < P / 4; r++) {
                int k = orden[r];
                if (k == iMejor) continue;
                permutacionGenes(&pob[k], &maestro);
                pob[k].fit = evaluar(&pob[k]);
                evals += N;
                reemplazos++;
            }
            printf("Estancamiento: %d individuos reemplazados en la generacion %ld.\n", reemplazos, g);
            sinMejora = 0;
        }

        if (conCruce) {
            /* Emparejamiento aleatorio y semillas por pareja (independiente del n. de hilos). */
            for (int k = 0; k < P; k++) orden[k] = k;
            for (int i = P - 1; i >= 1; i--) { int j = rngMenor(&maestro, i + 1); int t = orden[i]; orden[i] = orden[j]; orden[j] = t; }
            for (int k = 0; k < P; k++) semillas[k] = rngNext(&maestro);

            /* Cruce + seleccion elitista por pareja, en paralelo (parejas disjuntas). */
            long evCruce = 0, acept = 0;
            #pragma omp parallel reduction(+:evCruce,acept)
            {
                int tid = 0;
#ifdef _OPENMP
                tid = omp_get_thread_num();
#endif
                Trabajo *w = &trab[tid];
                #pragma omp for schedule(dynamic)
                for (int q = 0; q < P / 2; q++) {
                    int p1 = orden[2 * q], p2 = orden[2 * q + 1];
                    Rng r; rngSembrar(&r, semillas[q]);
                    evCruce += cruzarPar(&pob[p1], &pob[p2], &scr[2 * q], &scr[2 * q + 1], &r, w);
                    /* Los 2 mejores de {P1,P2,H1,H2} quedan en la poblacion; en empate
                       se prefiere a los hijos. Intercambio de punteros, sin copias. */
                    Crom c[4] = { scr[2 * q], scr[2 * q + 1], pob[p1], pob[p2] };
                    int id[4] = {0, 1, 2, 3};
                    for (int x = 1; x < 4; x++)
                        for (int y = x; y > 0 && c[id[y]].fit < c[id[y - 1]].fit; y--) { int t = id[y]; id[y] = id[y - 1]; id[y - 1] = t; }
                    acept += (id[0] < 2) + (id[1] < 2);   /* hijos que entran a la poblacion */
                    pob[p1] = c[id[0]]; pob[p2] = c[id[1]];
                    scr[2 * q] = c[id[2]]; scr[2 * q + 1] = c[id[3]];
                }
            }
            evals += (double)evCruce;
            hijosEvaluados += 2L * (P / 2);
            hijosAceptados += acept;
        }

        /* Mutacion en paralelo. Con mutacion guiada el error nunca aumenta, asi que
           se muta a todos; con mutacion ciega se protege al mejor (elitismo).
           Con recocido el error puede subir, pero el mejor historico vive en
           'mejor', que es una copia aparte. */
        iMejor = 0;
        for (int k = 1; k < P; k++) if (pob[k].fit < pob[iMejor].fit) iMejor = k;
        for (int k = 0; k < P; k++) semillas[k] = rngNext(&maestro);
        long evMut = 0;
        #pragma omp parallel for schedule(dynamic) reduction(+:evMut)
        for (int k = 0; k < P; k++) {
            if (!MUTACION_GUIADA && k == iMejor) continue;
            Rng r; rngSembrar(&r, semillas[k]);
            evMut += mutar(&pob[k], &r);
        }
        evals += (double)evMut;

        iMejor = 0;
        for (int k = 1; k < P; k++) if (pob[k].fit < pob[iMejor].fit) iMejor = k;
        g++;
        if (pob[iMejor].fit < mejor.fit) {
            cromCopiar(&mejor, &pob[iMejor]);
            printf("Optimo[%ld] = %d\n", g, mejor.fit);
            if (traza) fprintf(traza, "%ld,%.0f,%.4f,%d\n", g, evals, ahora() - t0, mejor.fit);
        }
        sinMejora = (mejor.fit < fitAntes) ? 0 : sinMejora + 1;
    }

    double t = ahora() - t0;
    if (traza) { fprintf(traza, "%ld,%.0f,%.4f,%d\n", g, evals, t, mejor.fit); fclose(traza); }
    if (mejor.fit == 0) printf("Solucion sin violaciones en la generacion %ld.\n", g);
    else if (presupuesto > 0 && evals >= presupuesto) printf("Presupuesto de %.3g evaluaciones agotado. Error residual = %d.\n", presupuesto, mejor.fit);
    else printf("Presupuesto de %ld generaciones agotado. Error residual = %d.\n", Gmax, mejor.fit);

    /* ---- Verificaciones finales ---- */
    int okS = esPermutacion(mejor.sid, N), okB = esPermutacion(mejor.bid, 2 * N);
    int recalculado = evaluar(&mejor);
    double misma = mismaCategoria(&mejor), c1 = correlacionVecinos(&mejor);
    printf("--- Verificacion ---\n");
    printf("Inventario de sitios conservado:  %s\n", okS ? "SI (permutacion exacta de M_base)" : "NO");
    printf("Inventario de enlaces conservado: %s\n", okB ? "SI (permutacion exacta de M_base)" : "NO");
    printf("Fitness recalculado = %d (registrado = %d)\n", recalculado, mejor.fit);
    printf("Error residual por incidencia = %.4f%%\n", 100.0 * mejor.fit / (4.0 * N));
    printf("Vecinos de misma categoria (crudo) = %.1f%%\n", 100.0 * misma);
    printf("C(1) vecinos inmediatos = %.4f\n", c1);
    printf("Generaciones = %ld\n", g);
    printf("Evaluaciones de sitio = %.0f\n", evals);
    if (conCruce)
        printf("CRUCE,modo=%s,reparar=%d,hijos_evaluados=%ld,hijos_aceptados=%ld,tasa_aceptacion=%.4f\n",
               g_modoCruce ? "rango" : "pmx", g_reparar, hijosEvaluados, hijosAceptados,
               hijosEvaluados ? (double)hijosAceptados / hijosEvaluados : 0.0);
    printf("Tiempo total de ejecucion: %.2f segundos\n", t);
    printf("RESUMEN,variante,L,omega,semilla,P,cruce,generaciones,evaluaciones,tiempo_s,error,error_pct,misma_cat,C1,inventario_ok,fitness_ok\n");
    printf("RESUMEN,%s,%d,%.3f,%llu,%d,%d,%ld,%.0f,%.3f,%d,%.5f,%.4f,%.4f,%d,%d\n",
           variante, L, omega, (unsigned long long)semilla, P, conCruce, g, evals, t, mejor.fit,
           100.0 * mejor.fit / (4.0 * N), misma, c1, okS && okB, recalculado == mejor.fit);
    exportar(&mejor, prefijo);

    if (trab) { for (int k = 0; k < nHilos; k++) trabajoLiberar(&trab[k]); free(trab); }
    for (int k = 0; k < P; k++) cromLiberar(&pob[k]);
    for (int k = 0; k < nScr; k++) cromLiberar(&scr[k]);
    cromLiberar(&mejor);
    free(pob); free(scr); free(orden); free(semillas); free(sVal); free(bVal);
    return 0;
}
