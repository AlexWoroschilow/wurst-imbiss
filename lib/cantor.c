/*
 * cantor.c
 * coding cantor style
 *
 * @author Steve Hoffmann
 * @date Tue 23 Nov 2006
 * @remark handle with care Uint has a limited range
 */
 
 #include "mathematics.h"
 #include "cantor.h"
 #include "memory.h"
 #include "basic-types.h"

 Uint code_2tupel_cantor(Uint i, Uint j) {
	Uint code;
	
	code =((i+j)*(i+j+1))/2+j;
	return code;
 }


 Uint codeCantor(vector_t *v) {
	Uint i;
	Uint code=0;
	
	if (LENGTHVEC(v) == 1) return VECTOR(v,0);
	code = VECTOR(v,LENGTHVEC(v)-1);
	
	for(i=LENGTHVEC(v)-1; i > 0; i--) {		
			code = code_2tupel_cantor(VECTOR(v,i-1), code);
	}	
	return code;
 }


Uint f_cantor(Uint w) {

	return (w*(w+1)*0.5);
}


Uint q_cantor(Uint z) {
	Uint v=0;

	while (f_cantor(v)<=z) v+=1;
	return (v-1);
}


vector_t *decode_2tupel_cantor(void *space, Uint i) {
	Uint j,y,x;
	vector_t *v=NULL;
	
	v=ALLOCMEMORY(space, v, vector_t, 1);
	INITVECTOR(v);

	j = q_cantor(i);
	y = i-f_cantor(j);
	x = j-y;
	
	APPENDVEC(space, v, x);
	APPENDVEC(space, v, y);

	return v;
}


vector_t *decodeCantor(void *space, Uint code, Uint n) {
	Uint i;
	vector_t *v = NULL;
	vector_t *r = NULL;
	
	v=ALLOCMEMORY(space, v, vector_t, 1);
	INITVECTOR(v);

	for (i=0; i < (n-1); i++) {
		r = decode_2tupel_cantor(space, code);		
		APPENDVEC(space, v, VECTOR(r,0));
		code = VECTOR(r,1);

		FREEMEMORY(space, r);
	}

	APPENDVEC(space, v, VECTOR(r,1));
	return (v);
}

