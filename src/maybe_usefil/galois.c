#include "galois.h"

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <memory.h>

#include "rscd.h"
#include "util.h"

/* Cuerpos y polinomios predefinidos */

/* polinomios irreducibles sobre GF(2) */
static unsigned long irreducibles[13]={ 
  // en octal!
  07, //2
  013, //3
  023, //4 010011
  045, //5 100101
  0103, //6 001000011
  0211,
  0435, 
  01021,
  02011,
  04005,
  010123,
  020033
};

/* Atributos de GF(q) */
/* potencia de 2 */
static  int m;
/* tamanio del cuerpo */
static  int q; 
/* orden multiplicativo */
static  int n;
/* tabla precalculada de antilogaritmos discretos */
static  gf_t mascara;
/* */
static  gf_t logaritmos[MAX_Q];
static  gf_t antilogaritmos[MAX_Q];

int gf_get_m() { return m; }
int gf_get_q() { return q; }
int gf_get_n() { return n; }
int gf_get_mascara() { return mascara; }

void gf_ini(unsigned int _m) {

  rscd_trace("gf_ini()");
  m=_m;
  q = 1<<_m;
  n = (1<<_m)-1; /* orden multiplicativo */
  mascara=n;
  if (_m == 1) {
    antilogaritmos[0]=1;
    logaritmos[0]=0;
    logaritmos[1]=0;
  } else {
    gf_ini_log();
  }
}

gf_t gf_antilog(unsigned int L){ return antilogaritmos[L]; }

/* 
 * se trata de tomar la raiz [x] del polinomio irreducible 
 * y multiplicarla por si misma hasta desbordar en x^8. Ahi se sustituye
 * x^8 por los otros coefs. del polinomio (en este caso 0,0,0,1,1,1,0,1)
 * y se sigue multiplicando por [x] hasta que un nuevo elemento desborde.
 * Nuevamente en el caso de gf8 se tiene (0,0,
 */
void gf_ini_log() {
  int i;
  int elem;
  unsigned long mascara_msb,mascara_lsb;
  unsigned long coefs_irr,coefs_lsb;
  //gf_t* plog;
  char tmp[128];

  rscd_trace("gf_ini_log()");
  //plog = &antilogaritmos[0];
  mascara_msb=q; // un 1 en el bit m (2^m)
  mascara_lsb=q-1; // todos 1s hasta bit m-1 
  coefs_irr=irreducibles[m-2];
  coefs_lsb=coefs_irr & mascara_lsb;
  antilogaritmos[0]=1; 
  logaritmos[0]=0;
  elem=02; // el elemento [x] = ..000010
  snprintf(tmp,127,"M_msb=%04lX, M_lsb=%04lX, C=%04lX, C_lsb=%04lX",mascara_msb,mascara_lsb,coefs_irr,coefs_lsb);
  rscd_debug(tmp);
  for (i=1; i < n; i++) {
    antilogaritmos[i]=elem;
    logaritmos[elem]=i;
    // elem[i+1]=elem[i]*[x] ...
    elem<<=1; 
    // ... hasta que desborde, 
    if (elem & mascara_msb) { 
      // ahi se sustituye x^m por P(x)-x^m=coefs_lsb
      // y se suma a elem. Esto es un xor  con coefs_lsb
      elem  = (elem ^ coefs_lsb) & mascara_lsb;      
    }
  }
}

/* operaciones */
gf_t gf_sum( gf_t a, gf_t b) {
  return (a ^ b); 
}

int gf_log( gf_t a) {
  if ((a==0)|| (a>n)) {
    return -1;
  } else {
    return logaritmos[a];
  }
}

gf_t gf_div( gf_t a, gf_t b) {
  if (a==0) return 0;
  else if (b==0) {
    rscd_error("Division por cero");
    exit(1);
    return 0; // no se como representar
  } else
    return antilogaritmos[(logaritmos[a]+n-logaritmos[b])%n]; 
}

gf_t gf_inv( gf_t a) {
  if (a ==0) return 0;
  if (a==1) return 1;
  return antilogaritmos[n-gf_log(a)]; 
}

gf_t gf_mul( gf_t a, gf_t b) {
  if ((a == 0) || (b==0)) 
    return 0;
  return antilogaritmos[(logaritmos[a]+logaritmos[b])%n]; 
}

gf_t gf_pot( gf_t a, int e) {
  if (a==0) return 0;
  return antilogaritmos[(gf_log(a)*e)%(n)];
}

char* gf_as_bits(gf_t a) {
  static char buf[35]; 
  int i;
  int mask = 1<<(m-1);
  buf[0]='[';
  for (i=1; i <= m; i++, mask>>=1) 
    buf[i] =(a & mask) ? '1':'0';
  buf[m+1]=']';
  buf[m+2]=0;
  return buf;
}

void gf_print(const int level) {
  char out[2048];
  gf_to_string(2047,out);
  rscd_log(level,out);
}

char* gf_to_string(const int maxlen,char* buffer) {
  int i;
  char tmp[128];
  buffer[0]=0;
  snprintf(buffer,maxlen,"GF(%d), orden=%d, mascara=%08X\n",q,n,mascara);
  for (i=0; i < n; i++) {
    snprintf(tmp,127,"%03d:%04X\t",i,antilogaritmos[i]);
    strncat(buffer,tmp,maxlen);
  }
  for (i=1; i <= n; i++) {
    snprintf(tmp,127,"%03d:%04X\t",i,logaritmos[i]);
    strncat(buffer,tmp,maxlen);
  }
  return buffer;
}

/* POLINOMIOS SOBRE GF(q) ---------------------------------- */

/* ciclo de vida */
gf_pol_t* gf_pol_new(int n, const gf_t* pcoefs) {
  gf_pol_t* pp;
  rscd_trace("gf_pol_new()");
  pp=malloc(sizeof(gf_pol_t));
  return gf_pol_ini(pp,n,pcoefs);
}

gf_pol_t* gf_pol_ini(gf_pol_t* pp, int n, const gf_t* pcoefs) {
  memset(pp->coefs,0,sizeof(gf_t)*MAX_N);
  if ((n!=0) && (pcoefs != NULL)) {
    /* copiamos todos los coeficientes */
    memcpy((void*)pp->coefs,(const void*)pcoefs,sizeof(gf_t)*(n+1));
    /* pero el orden definitivo lo da el ultimo elemento no nulo  */
    while (n && !pcoefs[n])
      n--;
    pp->n=n;
  } else {
    pp->n=0;
  }
  return pp;
}

/*
 * borra un polinomio dado
 */
void gf_pol_del(gf_pol_t* pp) {
  rscd_trace("gf_pol_del()");
  free(pp);
}

gf_t gf_pol_get_coef(gf_pol_t* pp, int l) {
  return pp->coefs[l];
}

void gf_pol_set_coef(gf_pol_t* pp, int l,gf_t a) {
  pp->coefs[l]=a;
}

gf_t gf_pol_eval(const gf_pol_t* pa, const gf_t x) {
  register int i;
  gf_t y;
  const gf_t* pca;
  pca=&pa->coefs[pa->n];
  y = 0;
  for (i=pa->n;i>=0;i--) {
    // printf("x * y +a =%d * %d + %d\n",x,y,*pca);
    y = gf_sum(gf_mul(y,x),*(pca--));
  }
  return y;
}


/* operaciones: c <- a op b (binarias) b <- op a */ 
gf_pol_t* gf_pol_sum(const gf_pol_t* pa, const gf_pol_t* pb, gf_pol_t* pc) {
  register int i;
  register int na,nb,nc;
  na=pa->n; nb=pb->n;
  nc=(na > nb)?na:nb;  
  if (nc>nb) {
    for(i=nc; i > nb; i--) {
      pc->coefs[i]=pa->coefs[i];
    }
    for(; i >= 0; i--) {
      pc->coefs[i]=gf_sum(pa->coefs[i],pb->coefs[i]);
    }
  } else {
    for(i=nc; i > na; i--) {
      pc->coefs[i]=pb->coefs[i];
    }
    for(; i >= 0; i--) {
      pc->coefs[i]=gf_sum(pa->coefs[i],pb->coefs[i]);
    }
  }
  while(nc && !pc->coefs[nc])
    nc--;
  pc->n=nc;
  return pc; 
}


gf_pol_t* gf_pol_mul(const gf_pol_t* pa, const gf_pol_t* pb, gf_pol_t* pc) {
  rscd_trace("gf_pol_mul()");
  /* convolucion c[i]=sum_{j=0}^{i}{a[j]b[i-j]}*/
  int na,nb,nc;
  register int i,j;
  register gf_t tmp;
  na=pa->n;
  nb=pb->n;
  nc=na+nb;
  for (i=0; i <= n; i++) {
    tmp=0;
    //printf("c[%d]=",i);
    for (j=0; j <=i; j++) {
      if (j > na) continue;
      if ((i-j) > nb) continue;
      tmp = gf_sum(tmp,gf_mul(pa->coefs[j],pb->coefs[i-j]));
      //printf("\ta[%d]=%d, b[%d]=%d, a[%d]b[%d]=%d\n",j,pa->coefs[j],i-j,pb->coefs[i-j],j,i-j,tmp);
      //printf("%04x ",tmp);
    }
    pc->coefs[i]=tmp;
  }  
  while(nc && !pc->coefs[nc])
    nc--;
  pc->n=nc;  
  return pc;
}

gf_pol_t* gf_pol_scale(const gf_t k, const gf_pol_t* pa, gf_pol_t* pb) {
  const register gf_t *pca;
  register gf_t *pcb;
  register int i,n;
  pca=&pa->coefs[0];
  pcb=&pb->coefs[0];
  n=pa->n;
  pb->n=n;
  for (i=0; i<=n; i++)  {
    pcb[i]=gf_mul(k,pca[i]);
  } 
  return pb;
}

void gf_pol_shift(gf_pol_t* pp, unsigned int e) {
  register int i;
  gf_t* pc;
  if (!e) 
    return;
  pc=&pp->coefs[0];
  for (i = pp->n; i>=0; i++)
    pc[i+e]=pc[i];
  for (i = e-1; i>=0; i++)
    pc[i]=0;  
}

gf_pol_t* gf_pol_copy(gf_pol_t* pa, const gf_pol_t* pb) {
  memcpy(pa,pb,sizeof(gf_pol_t));
  return pa;
}

gf_pol_t* gf_pol_div(const gf_pol_t* pa, const gf_pol_t* pb, 
		     gf_pol_t* pq, gf_pol_t* pr) 
{
  /* 
   * Algoritmo a/b: a=qb+r
   * sea a = a_nx^n+ ... + a_0 y b=b_nx^m+ ... b_0
   * 1. r:=a q:=0
   */
  /* 2. si m > n */  
  int m,n;
  int cero_izq;
  register int i;
  gf_t c;
  gf_t *pcr;
  const gf_t *pcb;
  rscd_trace("gf_pol_div()");

  gf_pol_copy(pr,pa);
  gf_pol_ini(pq,0,NULL);
  m=pb->n;
  n=pr->n;
  pcr=&pr->coefs[n];
  pcb=&pb->coefs[m];
  while (m <= n) {
    /* 2.b.1 k=m-n, c = (a_n/b_m) */
    c = gf_div(*pcr,*pcb);     
    /*printf("coef r=%x coef b=%x c=%x\n",*pcr,*pcb,c); */
    /* 2.b.2 q := q + c*x^k */
    pq->coefs[n-m]=c;
    if (!pq->n)
      pq->n=n-m;
    /* 2.b.3 r := r - c*x^k*b */
    cero_izq=1;
    /* printf("pcr:"); */
    for (i = m; i >= 0; i--,pcr--,pcb--) {
      *pcr=gf_sum(*pcr,gf_mul(c,*pcb));
      /* printf(" %x",*pcr); */
      if (cero_izq) {
	if (!*pcr) 
	  pr->n--;
	else
	  cero_izq=0;
      }
      if (pr->n <0) {
	pr->n=0;
	return pq;
      }
    }
    /* printf(" pr->n:%d\n",pr->n); */
    m=pb->n;
    n=pr->n;
    pcr=&pr->coefs[n];
    pcb=&pb->coefs[m];
  }
  return pq;
}

void gf_pol_print(const int log_level, const gf_pol_t* pp) {
  char out[512];
  rscd_trace("gf_pol_print()");
  gf_pol_to_string(pp,511,out);
  puts(out); /* rscd_log(log_level,out); */
}

char* gf_mon_to_string(const gf_t a, const int exp) {
  static char monstr[15];
  char tmp[10];
  monstr[0]=0;
  strcat(monstr,gf_as_bits(a));
  if (exp >= 1) {
    strcat(monstr,"x");
    if (exp > 1) {
      snprintf(tmp,9,"^%d",exp);
      strcat(monstr,tmp);
    }
  }
  return monstr;
}

char* gf_pol_to_string(const gf_pol_t* pp, const int maxlen,char* buffer){
  int i;
  int primero;
  primero=1;
  buffer[0]=0;
  if (pp->n == 0) {
    strncat(buffer,gf_as_bits(pp->coefs[0]),maxlen);
    return buffer;
  }
  for (i=pp->n;i >=0; i--) {
    //    if (pp->coefs[i] != 0) {
      if (!primero) {
	strncat(buffer, " + ",maxlen);
      } else {
	primero=0;
      }
      strncat(buffer,gf_mon_to_string(pp->coefs[i],i),maxlen);
      //}
  }
  return buffer;
}

int gf_pol_comp(const gf_pol_t* pa, const gf_pol_t* pb) {
  rscd_trace("gf_pol_comp()");
  if (pa->n > pb->n)
    return 1;
  else if (pa->n < pb->n)
    return -1;
  else
    return memcmp(pa->coefs,pb->coefs,(pa->n+1)*sizeof(gf_t));
}

int gf_pol_cero(const gf_pol_t* pa) {
  rscd_trace("gf_pol_cero()");
  register int i;
  i=pa->n;
  if (i == 0)
    return 1;

  while (i>=0)
    if (pa->coefs[i]) return 0;
  return 1;
}

gf_pol_t* gf_pol_der(const gf_pol_t* pa, gf_pol_t* pd) {
  /* d(x)=a'(x)=(a_nx^n+...+a_0)'=n*a_nx^(n-1)+...+a_(n-1)*/
  /* o sea que d.coefs[i-1]=i*a.coefs[i], i=1,...,n-1 */
  int i;
  int na,nd;
  const gf_t* pca;
  gf_t* pcd;
  na=pa->n;
  if (na==0) {
    pd->n=0;
    pd->coefs[0]=0;
    return pd;
  } else {
    nd=na-1;
    pca=pa->coefs + 1;
    pcd=pd->coefs;
    for (i=1; i<=na; i++,pcd++,pca++) {
      // bestialidad: tendria que verificarse que i < tamanio del cuerpo
      // no lo hago porque siempre voy a trabajar con polinomios de orden bastante
      // menor que el tamanio del cuerpo
      //  *pcd = gf_mul(*pca,i); 
      if (i%2==0) // si i es par entonces (a+a+...+a)=0 porque la caracteristica del cuerpo es 2
	*pcd=0;
      else // si es impar es = *pca
	*pcd = *pca;
    }
    while(nd && !pd->coefs[nd])
      nd--;
    pd->n=nd;
    return pd;
  }
}


int gf_pol_raices(const gf_pol_t* pp) {
  gf_t a;
  gf_t nr=0;
  char tmp[512];
  printf("Raices de %s: ",gf_pol_to_string(pp,511,tmp));
  for (a=0; a < gf_get_q(); a++) {
    if (!gf_pol_eval(pp,a)) {
      printf("%0X ",a);
      nr++;
    }    
  } 
  printf("\nTotal %d raices distintas. ",nr);
  if (nr == pp->n) { 
    printf(" Todas simples.\n"); 
  } else {
    printf(" Hay raices con multiplicidad >1.\n");
  }
  return nr;
}
