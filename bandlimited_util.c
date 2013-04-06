/**

This software is copyrighted by Miller Puckette and others.  The following
terms (the "Standard Improved BSD License") apply to all files associated with
the software unless explicitly disclaimed in individual files:

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above  
   copyright notice, this list of conditions and the following 
   disclaimer in the documentation and/or other materials provided
   with the distribution.
3. The name of the author may not be used to endorse or promote
   products derived from this software without specific prior 
   written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,   
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "m_pd.h"
#include "bandlimited_defs.h"

/*
 * This function performs a 4 point interpolation lookup on the table.
 * code borrowed from tabread4~
 *
 * param float *  pointer to wavetable
 * param t_float phase to lookup
 *
 * return t_float result of the table lookup
 */
t_float bandlimited_read4(float *table, t_float p) {
	
    double dphase;
    int normhipart;
    union tabfudge tf;
    float *tab = table, *addr, a, b, c,d,cminusb, frac;
	
	
    tf.tf_d = UNITBIT32;
    normhipart = tf.tf_i[HIOFFSET];
	
	
	
	dphase = (double)(p * (float)(BANDLIMITED_TABSIZE)) + UNITBIT32;
	tf.tf_d = dphase;
	addr = tab + (tf.tf_i[HIOFFSET] & (BANDLIMITED_TABSIZE-1)) + 1;
	tf.tf_i[HIOFFSET] = normhipart;
	frac = tf.tf_d - UNITBIT32;
	a = addr[-1];
	b = addr[0];
	c = addr[1];
	d = addr[2];
	
	
	cminusb = c-b;
	return b + frac * (
					   cminusb - 0.1666667f * (1.-frac) * (
														   (d - a - 3.0f * cminusb) * frac + (d + 2.0f*a - 3.0f*b)
														   )
					   );
}

#ifdef DEBUG
/*
 * This function performs sin(2pi * x) using linear interpolation on top of a wavetable.
 * It's here for testing purposes.
 *
 * param t_float phase
 *
 * return double evaluation of sin function
*/
double bandlimited_sin_lin(t_float p) {
    double dphase;
    int normhipart;
    union tabfudge tf;
    float *tab = bandlimited_sin_table, *addr, f1, f2, frac;
	
    tf.tf_d = UNITBIT32;
    normhipart = tf.tf_i[HIOFFSET];
	
	
	
	dphase = (double)(p * (float)(BANDLIMITED_TABSIZE)) + UNITBIT32;
	tf.tf_d = dphase;
	addr = tab + (tf.tf_i[HIOFFSET] & (BANDLIMITED_TABSIZE-1))+1;
	tf.tf_i[HIOFFSET] = normhipart;
	frac = tf.tf_d - UNITBIT32;
	f1 = addr[0];
	f2 = addr[1];
	return (f1 + frac * (f2 - f1));
}
#endif

