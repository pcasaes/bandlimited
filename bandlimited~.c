/**


Copyright (C) 2010 Paulo Casaes

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


-- 
http://http://gitorious.org/~saturno
mailto:irmaosaturno@gmail.com


 
 Code was liberally borrowed from d_osc.c

*/

#include "m_pd.h"
#include <math.h>
#include <string.h>

#define BANDLIMITED_PI     3.14159265358979323846
#define BANDLIMITED_PISQ   9.8696044010893586188

#ifdef BANDLIMITED_MAXHARMONICS
#else
#define BANDLIMITED_MAXHARMONICS 734 // about 24hz at 44.1khz
#endif

#define GETSTRING(s) (s)->s_name

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
#error No byte order defined                                                    
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

union tabfudge
{
    double tf_d;
    int32 tf_i[2];
};

static t_class *bandlimited_class;

typedef struct _bandlimited
	{
		t_object x_obj;
		
		//phasor
		double x_phase;
		float x_conv;
		float x_f;      /* scalar frequency */
		
		//bandlimited
		float s_nq;
		float cutoff;
		int max_harmonics_set;
		int max_harmonics_mod;
		int max_harmonics;
		

		
		//type
		t_float (*generator)(int, t_float, t_float);
		t_float (*dutycycle)(t_float);

		
	} t_bandlimited;


static t_float bandlimite_nodutycycle(t_float in) {
	return 1.0f;
}

static t_float bandlimite_pulsedutycycle(t_float in) {
	t_float dc;
	dc = in;
	if(dc < 0.0)
		dc *= -1.0f;
	dc = dc - (int)(dc/1.0f);
	dc = 1.0f - dc;
	if(dc <= 0.0f)
		dc=1.0f;
	return dc;
}

static inline t_float bandlimited_sin(t_float in) {
    double dphase;
    int normhipart;
    union tabfudge tf;
    float *tab = cos_table, *addr, f1, f2, frac;
	
    tf.tf_d = UNITBIT32;
    normhipart = tf.tf_i[HIOFFSET];
	
	in = 0.25f - in;
	
	dphase = (double)(in * (float)(COSTABSIZE)) + UNITBIT32;
	tf.tf_d = dphase;
	addr = tab + (tf.tf_i[HIOFFSET] & (COSTABSIZE-1));
	tf.tf_i[HIOFFSET] = normhipart;
	frac = tf.tf_d - UNITBIT32;
	f1 = addr[0];
	f2 = addr[1];
	return (f1 + frac * (f2 - f1));
	
	
	
}

static t_float bandlimited_square(int max_harmonics, t_float p, t_float dutycycle) {
	int i;
	t_float sum=0.0f;
	
	
	for(i = 1; i <= max_harmonics; i += 2) {
		
		sum += bandlimited_sin(p * i )/i;
	}
	
	return  4.0f *  sum / BANDLIMITED_PI;
}

static t_float bandlimited_pulse(int max_harmonics, t_float p, t_float dutycycle) {
	int i;
	t_float sum=0.0f;
	
	
	for(i = 1; i <= max_harmonics; i += 2) {
		
		sum += bandlimited_sin(p * i * dutycycle)/i;
	}
	
	return  4.0f *  sum / BANDLIMITED_PI;
}


static t_float bandlimited_triangle(int max_harmonics, t_float p, t_float dutycycle) {
	int i;
	t_float sum=0.0f;
	
	
	for(i = 1; i <= max_harmonics; i += 2) {
		
		//sum += (powf(-1.0f, (i-1)/2.0f) * bandlimited_sin(p * i))/powf(i, 2.0f) ;
		sum += (bandlimited_sin(p * i)/powf(i, 2.0f)) * (i%4==3 ? -1 : 1 );
		
	}
	
	return  8.0f *  sum / BANDLIMITED_PISQ;
}

static inline t_float bandlimited_sawwave(int max_harmonics, t_float p, t_float dutycycle) {
	int i;
	t_float sum=0.0f;
	
	for(i = 1; i <= max_harmonics; i++) {
		
		sum += bandlimited_sin(p * i)/i;
	}
	
	return   sum / BANDLIMITED_PI;
}

static t_float bandlimited_saw(int max_harmonics, t_float p, t_float dutycycle) {
	
	return -2.0f * bandlimited_sawwave(max_harmonics, p, dutycycle);

}

static t_float bandlimited_rsaw(int max_harmonics, t_float p, t_float dutycycle) {
	return 2.0f * bandlimited_sawwave(max_harmonics, p, dutycycle);
}


static inline int bandlimited_typeset(t_bandlimited *x, t_symbol *type) {
	if(strcmp(GETSTRING(type), "saw") == 0) {
		x->generator=  &bandlimited_saw;
		x->dutycycle= &bandlimite_nodutycycle;
		x->max_harmonics_mod=1;
	} else if(strcmp(GETSTRING(type), "rsaw") == 0) {
		x->generator=  &bandlimited_rsaw;
		x->dutycycle= &bandlimite_nodutycycle;
		x->max_harmonics_mod=1;
	} else if(strcmp(GETSTRING(type), "square") == 0) {
		x->generator=  &bandlimited_square;
		x->dutycycle= &bandlimite_nodutycycle;
		x->max_harmonics_mod=2;
	} else if(strcmp(GETSTRING(type), "triangle") == 0) {
		x->generator=  &bandlimited_triangle;
		x->dutycycle= &bandlimite_nodutycycle;
		x->max_harmonics_mod=2;
	} else if(strcmp(GETSTRING(type), "pulse") == 0) {
		x->generator=  &bandlimited_pulse;
		x->dutycycle= &bandlimite_pulsedutycycle;
		x->max_harmonics_mod=2;
	} else {
		goto type_unknown;

	}
	x->max_harmonics = x->max_harmonics * x->max_harmonics_mod;
	return 0;
type_unknown:
	return 1;
}

static void *bandlimited_new(t_symbol *type, t_floatarg f) {
    t_bandlimited *x = (t_bandlimited *)pd_new(bandlimited_class);
	x->cutoff=0;
    x->x_f = f;
	x->s_nq=0;
	x->max_harmonics=BANDLIMITED_MAXHARMONICS;
	if(bandlimited_typeset(x, type) == 1) {
		error("bandlimited~: Uknown type %s, using saw", GETSTRING(type));
		x->generator=  &bandlimited_saw;
	}
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_float, gensym("ft1"));
	inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);	
    x->x_phase = 0;
    x->x_conv = 0;	
    outlet_new(&x->x_obj, gensym("signal"));

    return (x);
}

static void bandlimited_ft1(t_bandlimited *x, t_float f)
{
    x->x_phase = f;
}

static void bandlimited_cutoff(t_bandlimited *x, t_float f)
{
	if(x->s_nq != 0 && f > x->s_nq) 
		error("bandlimited~: %f is greater than the nyquist limit %f, ignoring", f, x->s_nq);
	else if(f < 1 )
		x->cutoff = x->s_nq-1;
	else 
	x->cutoff = f;
}

static void bandlimited_max(t_bandlimited *x, t_float f)
{
	int val = (int)f;
	if(val < 1) {
		val = BANDLIMITED_MAXHARMONICS;
	}
	else if(val > BANDLIMITED_MAXHARMONICS) 
		post("bandlimited~: maximum number of harmonics %d might be too high. you are warned", val);
	x->max_harmonics_set = val;
	x->max_harmonics = x->max_harmonics_set * x->max_harmonics_mod;
	
}





static void bandlimited_type(t_bandlimited *x, t_symbol *type)
{
	if(bandlimited_typeset(x, type) == 1) {
		error("bandlimited~: Uknown type %s, leaving as is", GETSTRING(type));
	}
}




static inline t_float bandlimited_phasor(t_bandlimited *x, t_float in)
{
    double dphase = x->x_phase + UNITBIT32;
    union tabfudge tf;
    int normhipart;
    float conv = x->x_conv;
	t_float out;
	
    tf.tf_d = UNITBIT32;
    normhipart = tf.tf_i[HIOFFSET];
    tf.tf_d = dphase;
	

        tf.tf_i[HIOFFSET] = normhipart;
        dphase += in * conv;
        out = tf.tf_d - UNITBIT32;
        tf.tf_d = dphase;

    tf.tf_i[HIOFFSET] = normhipart;
    x->x_phase = tf.tf_d - UNITBIT32;
    return out;
}





static t_int *bandlimited_perform(t_int *w)
{
    t_bandlimited *x = (t_bandlimited *)(w[1]);
    t_float *in = (t_float *)(w[2]);
    t_float *dutycycle = (t_float *)(w[3]);
    t_float *out = (t_float *)(w[4]);
    int n = (int)(w[5]);
	t_float p;


	int max_harmonics = (int)( x->cutoff / *in);
	if(max_harmonics > x->max_harmonics)
		max_harmonics = x->max_harmonics;
	
    while (n--)
    {

		p = bandlimited_phasor(x, *in++);
		*out++ = x->generator(max_harmonics, p, x->dutycycle(*dutycycle++));
		

		
    }
    return (w+6);
}

static void bandlimited_dsp(t_bandlimited *x, t_signal **sp)
{
	x->x_conv = 1./sp[0]->s_sr;
	x->s_nq = sp[0]->s_sr / 2.0f;
	if(x->cutoff == 0)
		x->cutoff = ( x->s_nq)  - 1.0;
	else if(x->cutoff > x->s_nq) {
		error("bandlimited~: %f is greater than the nyquist limit %f, you are warned", x->cutoff, x->s_nq);
	}
    
	dsp_add(bandlimited_perform, 5, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_n);
	//dsp_add(bandlimited_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

static void bandlimite_dmaketable(void)
{
    int i;
    float *fp, phase, phsinc = (2. * 3.14159) / COSTABSIZE;
    union tabfudge tf;
    
    if (cos_table) return;
    cos_table = (float *)getbytes(sizeof(float) * (COSTABSIZE+1));
    for (i = COSTABSIZE + 1, fp = cos_table, phase = 0; i--;
		 fp++, phase += phsinc)
		*fp = cos(phase);
	
	/* here we check at startup whether the byte alignment
	 is as we declared it.  If not, the code has to be
	 recompiled the other way. */
    tf.tf_d = UNITBIT32 + 0.5;
    if ((unsigned)tf.tf_i[LOWOFFSET] != 0x80000000)
        bug("cos~: unexpected machine alignment");
}


extern void bandlimited_tilde_setup(void)
{
    bandlimited_class = class_new(gensym("bandlimited~"), (t_newmethod)bandlimited_new, 0,
						  sizeof(t_bandlimited), 0, A_DEFSYMBOL, A_DEFFLOAT, 0);
    CLASS_MAINSIGNALIN(bandlimited_class, t_bandlimited, x_f);
    class_addmethod(bandlimited_class, (t_method)bandlimited_dsp, gensym("dsp"), 0);
    class_addmethod(bandlimited_class, (t_method)bandlimited_ft1,
					gensym("ft1"), A_FLOAT, 0);	
    class_addmethod(bandlimited_class, (t_method)bandlimited_type,
					gensym("type"), A_SYMBOL, 0);	
    class_addmethod(bandlimited_class, (t_method)bandlimited_cutoff,
					gensym("cutoff"), A_FLOAT, 0);		
    class_addmethod(bandlimited_class, (t_method)bandlimited_max,
					gensym("max"), A_FLOAT, 0);		
		
    bandlimite_dmaketable();
	post("bandlimited~: band limited signal generator. Using %d as the default maximum harmonics (to redefine compile with -DBANDLIMITED_MAXHARMONICS=x flag).", BANDLIMITED_MAXHARMONICS);


}