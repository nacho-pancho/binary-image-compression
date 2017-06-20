#include <cstdio>
#include <cctype>
#include "pnm.h"

static int skip_comments(FILE *fp)
{
    int ch;
    char line[100];

    while ((ch = fgetc(fp)) != EOF && isspace(ch)) 
        ;
    if (ch == '#') {
        fgets(line, sizeof(line), fp);
        skip_comments(fp);
    } else
        fseek(fp, -1, SEEK_CUR);
    return 0;
}

int read_pnm_header(FILE* fhandle, int& tipo, int& ancho, int& alto, int& maxval) {
  int res;
  /* ya leemos toda la metadata */    
  if (fgetc(fhandle) != 'P') { 
    printf("El archivo no es un PNM valido.\n");
    fclose(fhandle); return -1;
  }
  tipo = fgetc(fhandle)-'0';
  if ( (tipo != 2) && (tipo != 5) && (tipo != 6))  {
    printf("Error: el archivo no es un PNM soportado (2,5 o 6, pero tipo=%d)\n",tipo);
    fclose(fhandle); return -1;
  }
  res = skip_comments(fhandle); if (res) return -1;
  res = fscanf(fhandle, "%d", &ancho); if (res <=0) return -1;
  res = skip_comments(fhandle); if (res) return -1;
  res = fscanf(fhandle, "%d", &alto); if (res <= 0) return -1;
  res = skip_comments(fhandle); if (res) return -1;
  res = fscanf(fhandle, "%d", &maxval); if (res <= 0) return -1;
  fgetc(fhandle);
  /* res = skip_comments(fhandle); if (res) return -1; */
  /* fgetc(fhandle); fscanf ya lee uno mas de la cuenta*/ 
  return 0;
}


int read_pgm_p2_data(FILE* fhandle, int ancho, int alto, int maxval, pixel_t* buf) {
  int i,n,res = 1;
  n = ancho*alto;
  for (i = 0; (i < n) && (res > 0); ++i) {
    res = fscanf(fhandle,"%u",&buf[i]);
  }
  return 0;
}

int read_pgm_p5_data(FILE* fhandle, int ancho, int alto, int maxval, pixel_t* buf) {
  int nread,nexpected;
  nexpected = ancho*alto;
  if (maxval < 256) {
    nread = 0;
    for (int i = 0; i < nexpected; ++i) {
      int pix = fgetc(fhandle);
      if (pix < 0) break;
      buf[i] = pix;
      nread ++;
    }
  } else {
    unsigned char hilo[2];
    nread = 0;
    for (int i = 0; i < nexpected; ++i) {
      if ( fread(&hilo[0],1,2,fhandle)==2 )
	nread ++;
      buf[i] = (hilo[0]<<8)+hilo[1];
    }
  }
  if (nread < nexpected) {
    printf("Solo se leyeron %d bytes pero se esperaban %d: feof()=%d\tferror()=%d\n",nread,nexpected,feof(fhandle),ferror(fhandle));
  }
  return nexpected-nread;
}

int read_pgm_data(FILE* fhandle, int tipo, int ancho, int alto, int maxval, pixel_t* buf) {
  switch (tipo) {
  case 2:
    return read_pgm_p2_data(fhandle, ancho,  alto, maxval, buf);
  case 5:
    return read_pgm_p5_data(fhandle, ancho,  alto, maxval, buf);
  default:
    return 0;
  }
}

int write_ppm_header(int tipo,
		     int ancho, 
		     int alto, 
		     int maxval, 
		     FILE* fw) {
  fprintf(fw, "P%c\n",tipo+'0'); 
  fprintf(fw, "%d %d\n", ancho, alto);
  fprintf(fw, "%d\n", maxval);
  return ferror(fw);
}

int write_p2_data(const pixel_t* pixels, const int npixels, const int maxval, FILE* fw) {
  register int i;
  for (i = 0; i < npixels; ++i) {
    fprintf(fw,"%d\t",pixels[i]);
    if (!((i+1) % 20)) fprintf(fw,"\n");
  }
  return ferror(fw);
}

int write_p5_data(const pixel_t* pixels, const int npixels, const int maxval, FILE* fw) {
  register int i;
  if (maxval < 256) {
    for (i = 0; (i < npixels) && !ferror(fw); ++i) {
      fputc((unsigned char) pixels[i],fw);
    }
  } else {
    for (i = 0; (i < npixels) && !ferror(fw) ; ++i) {
      const pixel_t unsigned pix = (pixel_t)pixels[i];
      fputc((pix >> 8) & 0xff,fw);
      fputc(pix & 0xff,fw);
    }
  }
  return ferror(fw);
}

int write_pgm_p2(const pixel_t* pixels, 
		 int ancho, 
		 int alto, 
		 int maxval, 
		 const char* ruta_archivo) {
  int i;
  FILE* fw;

  fw = fopen(ruta_archivo, "wb"); 
  if (fw == NULL) {
    i = ferror(fw);
    fclose(fw);
    printf("Error al intentar escribir archivo %s\n",ruta_archivo);
    return i;
  }
  write_ppm_header(2,ancho,alto,maxval,fw);
  write_p2_data(pixels,ancho*alto,maxval,fw);
  i = fclose(fw);
  return i;
}

int write_pgm_p5(const pixel_t* pixels, 
		 int ancho, 
		 int alto, 
		 int maxval, 
		 const char* ruta_archivo) {
  int i;
  FILE* fw;
  fw = fopen(ruta_archivo, "wb"); 
  if (fw == NULL) {
    i = ferror(fw);
    fclose(fw);
    printf("Error al intentar escribir archivo %s\n",ruta_archivo);
    return i;
  }
  write_ppm_header(5,ancho,alto,maxval,fw);
  write_p5_data(pixels,ancho*alto,maxval,fw);
  i = fclose(fw);
  return i;
}


int write_pgm(const pixel_t* pixels, 
	      int tipo,
	      int ancho, 
	      int alto, 
	      int maxval, 
	      const char* ruta_archivo) {
  switch (tipo) {
  case 2:
    return write_pgm_p2(pixels, 		
			ancho, 
			alto, 
			maxval, 
			ruta_archivo);
  case 5:
    return write_pgm_p5(pixels, 
			ancho, 
			alto, 
			maxval, 
			ruta_archivo);
  default:
    return -1;
  }
}


int read_ppm_data(FILE* fimg, int tipo, int ancho, int alto, int maxval, pixel_t* buf) {
  register int i, aux;
  unsigned char r,g,b;
  const int n = ancho*alto;
  aux = 3;
  for (i = 0; (i < n) && (aux == 3); ++i) {
    aux  = fread(&r,sizeof(char),1,fimg);  /* maximo asumido de 8 bits por canal! */
    aux += fread(&g,sizeof(char),1,fimg);
    aux += fread(&b,sizeof(char),1,fimg);
    buf[i] = (r << 16) + (g << 8) + b;
  }  
  if (i < n)
    return -1;
  else
    return 0;
}


int write_ppm(const pixel_t* pixels,
	      int tipo,
	      int ancho, 
	      int alto, 
	      int maxval, 
	      const char* ruta_archivo) {
  const int n = ancho*alto;   
  register int i;
  FILE* fimg;
  fimg = fopen(ruta_archivo, "wb"); 
  if (fimg == NULL) {
    i = ferror(fimg);
    fclose(fimg);
    printf("Error writing to file %s\n",ruta_archivo);
    return i;
  }
  fprintf(fimg, "P6\n");
  fprintf(fimg, "%d %d\n", ancho, alto);
  fprintf(fimg, "%d\n", maxval);
  for (i = 0; i < n; ++i,++pixels) {
    int p = *pixels;
    fputc((p >> 16) & 0xff, fimg); /* r */
    fputc((p  >> 8) & 0xff, fimg); /* g */
    fputc((p      ) & 0xff, fimg); /* b */
  }
  fclose(fimg);
  return 0;
}
