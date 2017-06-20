/**
 * Implementación de cuerpos finitos de Galois de tamanho 2^m
 * y polinomios sobre ellos.
 */
#ifndef GALOIS_H
#define GALOIS_H

/* Maximo tamanio de un cuerpo finito */
#define MAX_Q 256


/* CUERPOS GF(q) -------------------------------------------- */

typedef unsigned short gf_t;

/* ciclo de vida */
/*
 * Inicializacion del cuerpo por defecto
 * Esto require los siguientes pasos:
 *
 * 1-Encontrar un polinomio de F2[x] irreducible de orden m
 * 2-Tomar el elemento [x] (siempre primitivo) y multiplicarlo por si
 *   mismo hasta encontrar todos los elementos del cuerpo y asi 
 *   construir la tabla de logaritmos.
 *
 * El primer paso para m > 1 es relativamente complicado y se hace en una 
 * funcion aparte. Para m = 1 es trivial y este caso se usa para operar 
 * con los polinomios de F2 usando las funciones de polinomios.
 */
void gf_ini(unsigned int m);

int gf_get_m();

int gf_get_q();

int gf_get_n();

int gf_get_mascara();

/* operaciones */
gf_t gf_antilog(unsigned int L);

gf_t gf_sum( gf_t a, gf_t b);

/*
gf_t gf_neg( gf_t a);

gf_t gf_sub( gf_t a, gf_t b);
*/

gf_t gf_mul( gf_t a, gf_t b);

gf_t gf_div( gf_t a, gf_t b);

gf_t gf_inv( gf_t a);

int gf_log( gf_t a);

gf_t gf_pot( gf_t a, int e);

char* gf_as_bits(gf_t a);

char* gf_to_string(const int maxlen,char* buffer);

void gf_print(const int level);

/* POLINOMIOS SOBRE GF(q) ---------------------------------- */

/* Maximo orden de polinomios */
#define MAX_N 512

typedef struct {
  /* orden */
  int n;
  gf_t coefs[MAX_N];
} gf_pol_t;

/* ciclo de vida */

/* 
 * crea un polinomio de orden n sobre el cuerpo especificado y 
 * con los coeficientes especificados o todo 0 si coefs = NULL
 */
gf_pol_t* gf_pol_new(int n, const gf_t* pcoefs);

gf_pol_t* gf_pol_ini(gf_pol_t* pp, int n, const gf_t* pcoefs);

/*
 * borra un polinomio dado
 */
void gf_pol_del(gf_pol_t* pf);

gf_t gf_pol_get_coef(gf_pol_t* pp, int l);

void gf_pol_set_coef(gf_pol_t* pp, int l,gf_t a);

gf_t gf_pol_eval(const gf_pol_t* pa, const gf_t x);

/* operaciones: c <- a op b (binarias) b <- op a */ 
gf_pol_t* gf_pol_sum(const gf_pol_t* pa, const gf_pol_t* pb, gf_pol_t* pc);

/*
gf_pol_t* gf_pol_neg(const gf_pol_t* pa, gf_pol_t* pb);

gf_pol_t* gf_pol_sub(const gf_pol_t* pa, const gf_pol_t* pb, gf_pol_t* pc);
*/

gf_pol_t* gf_pol_mul(const gf_pol_t* pa, const gf_pol_t* pb, gf_pol_t* pc);

gf_pol_t* gf_pol_scale(const gf_t k, const gf_pol_t* pa, gf_pol_t* pb);

/* multiplica al polinomio por x^e */
void gf_pol_shift(gf_pol_t* pp, unsigned int e);

gf_pol_t* gf_pol_div(const gf_pol_t* pa, const gf_pol_t* pb, 
		     gf_pol_t* pq, gf_pol_t* pr);

gf_pol_t* gf_pol_copy(gf_pol_t* pa, const gf_pol_t* pb);

/* devuelve: 0 si son iguales, 1 si A(x) mayor que B(X) o -1 en caso contrario */
int gf_pol_comp(const gf_pol_t* pa, const gf_pol_t* pb);

int gf_pol_cero(const gf_pol_t* pa);

gf_pol_t* gf_pol_der(const gf_pol_t* pa, gf_pol_t* pd);

char* gf_pol_to_string(const gf_pol_t* pp, const int maxlen,char* buffer);

void gf_pol_print(const int level, const gf_pol_t* pp);

int gf_pol_raices(const gf_pol_t* pp);

/* algoritmos */

/*
 * Dado un polinomio irreducible de orden n genera la tabla
 * de logaritmos para el cuerpo finito de tamanio 2^n
 */
void gf_ini_log();

/* algoritmo de euclides */

#endif
