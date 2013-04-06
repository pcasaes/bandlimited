/**
 
 
Apache License 2.0

bandlimited~
    Copyright [2010] Paulo Casaes

      This product includes software developed at
      Gitorious (https://gitorious.org/bandlimited).
 
 -- 
 http://http://gitorious.org/~saturno
 mailto:irmaosaturno@gmail.com
 
 */
 
#ifndef BANDLIMITED_DEFS_H_
#define BANDLIMITED_DEFS_H_


 
 
#define BANDLIMITED_PI     3.14159265358979323846
#define BANDLIMITED_PISQ   9.8696044010893586188

#ifdef BANDLIMITED_MAXHARMONICS
#else
#define BANDLIMITED_MAXHARMONICS 1104 // down to about 20hz at 44.1khz, 40hz at 88.2khz...
#endif

#define BANDLIMITED_TABSIZE 2048						//2048
#define BANDLIMITED_FINVNPOINTS	0.00048828125				//0.00048828125   // 1.0 / BANDLIMITED_TABSIZE

#define GETSTRING(s) (s)->s_name
#define ISFLOAT(a) (a.a_type==A_FLOAT)
#define ISSYMBOL(a) (a.a_type==A_SYMBOL)

#define UNITBIT32 1572864.  /* 3*2^19; bit 32 has place value 1 */

/* machine-dependent definitions.  These ifdefs really
 should have been by CPU type and not by operating system! */
#ifdef IRIX
/* big-endian.  Most significant byte is at low address in memory */
#define HIOFFSET 0    /* word offset to find MSB */
#define LOWOFFSET 1    /* word offset to find LSB */
#define int32 long  /* a data type that has 32 bits */
#endif /* IRIX */

#ifdef MSW
/* little-endian; most significant byte is at highest address */
#define HIOFFSET 1
#define LOWOFFSET 0
#define int32 long
#endif



#if defined(__FreeBSD__) || defined(__APPLE__)
#include <machine/endian.h>
#endif

#ifdef __linux__
#include <endian.h>
#endif

#if defined(__unix__) || defined(__APPLE__)
#if !defined(BYTE_ORDER) || !defined(LITTLE_ENDIAN)                         
#if defined(__USE_BSD)
#error No byte order defined                                                    
#else
# define LITTLE_ENDIAN  __LITTLE_ENDIAN
# define BIG_ENDIAN     __BIG_ENDIAN
# define PDP_ENDIAN     __PDP_ENDIAN
# define BYTE_ORDER     __BYTE_ORDER
#endif
#endif                                                                          

#if BYTE_ORDER == LITTLE_ENDIAN                                             
#define HIOFFSET 1                                                              
#define LOWOFFSET 0                                                             
#else                                                                           
#define HIOFFSET 0    /* word offset to find MSB */                             
#define LOWOFFSET 1    /* word offset to find LSB */                            
#endif /* __BYTE_ORDER */                                                       
#include <sys/types.h>
#define int32 int32_t
#endif /* __unix__ or __APPLE__*/


#define BANDLIMITED_INCREMENT 8						//16		4
#define BANDLIMITED_HAMSTART 1104					//1104
#define BANDLIMITED_HAMSIZE 138						//69			276  BANDLIMITED_HAMSTART / BANDLIMITED_INCREMENT

#define DEBUG 0
#ifdef DEBUG
#define debug(x) x
#else
#define debug(x)
#endif


union tabfudge
{
    double tf_d;
    int32 tf_i[2];
};

static long bandlimited_count=0l;
static float *bandlimited_sin_table=0;
static float **bandlimited_triangle_table=0;
static float **bandlimited_sawwave_table=0;
static float **bandlimited_sawtriangle_table=0;
static float **bandlimited_square_table=0;


#endif /*BANDLIMITED_DEFS_H_*/
