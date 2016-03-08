#ifndef LIST_H
#define LIST_H

/*
 *
 *	list.h
 *  declarations for a linearly linked list
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 12/27/06 14:00:46 CET  
 *
 */

 #include "basic-types.h"

 #define LISTNILVALUE -1

 typedef struct {
 	void *data;
	Uint next,
		 prev;
 } Listelem;

#define GETLISTELEM(L, A)  (L)->nodes[A].data

typedef struct {
	Listelem *nodes;
	Uint allocated,
		 nextfree,
		 firstNode,
		 lastNode,
		 length;
} List;

 List* initList (void *, Uint elems);
 void insertBefore (void *, List *, Uint, void *);
 void insertAfter (void *, List *, Uint, void *);
 void unlinkListElem (List *, Uint);
 void wrapList (void *, List *l, void (*rmv)(void *, void *));
 void **dataList (void *, List *);
 void sweepList (void *, List *l);
 
#endif
 
