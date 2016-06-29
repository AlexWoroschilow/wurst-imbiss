
/*
 *  list.c
 *  implementation for a linearly linked list
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/27/06 13:41:58 CET
 *  
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include "basic-types.h"
 #include "memory.h"
 #include "list.h"

/*--------------------------------- initList ---------------------------------
 *    
 * initializes a list allocating the first list element
 * 
 */
 
List*
initList (void *space, Uint elems)
{
    List *l; 
    
	if (elems == 0)
	  return NULL;
	
	l = ALLOCMEMORY(space, NULL, List, 1);	
	l->nodes = ALLOCMEMORY(space, NULL, Listelem, elems);
	l->allocated = elems;
	l->length = 0;
	l->firstNode = 0;
    l->nextfree = 0;
    l->lastNode = 0;

	return l;
}

void
appendlist(void *space, List *l, void *data) {

	if (l->nextfree >= l->allocated-2) {
		l->nodes = ALLOCMEMORY(space, l->nodes, Listelem, l->allocated*2);
		l->allocated *= 2;
	}
	  
    if (l->length==0) {
        l->firstNode = l->nextfree;
        l->lastNode = LISTNILVALUE;
    } else {
        l->nodes[l->lastNode].next = l->nextfree;   
    }
    
    l->nodes[l->nextfree].data = data;
	l->nodes[l->nextfree].next = LISTNILVALUE;
    l->nodes[l->nextfree].prev = l->lastNode; 
    l->lastNode = l->nextfree; 
    l->nextfree++;
    l->length++;
       
    return;
}

void
insertBefore(void *space, List *l, Uint cur, void *data) {	
	
   if (cur > l->nextfree && cur != LISTNILVALUE)
	  return ;
	
   if (l->length==0) {
/*        appendlist(space, l, data);
          return;*/
    } 
  
	
	if (l->nextfree >= l->allocated-2) {
		l->nodes = ALLOCMEMORY(space, l->nodes, Listelem, l->allocated*2);
		l->allocated *= 2;
	}
	  
	l->nodes[l->nextfree].data = data;
	l->nodes[l->nextfree].next = cur;
	
	if (cur != LISTNILVALUE) {
	  l->nodes[l->nextfree].prev = l->nodes[cur].prev;
	  l->nodes[cur].prev = l->nextfree;
	} else {
	  l->nodes[l->nextfree].prev = LISTNILVALUE;
	}

	if (l->nodes[l->nextfree].prev == LISTNILVALUE)
	  l->firstNode = l->nextfree;
	else
	  l->nodes[l->nodes[l->nextfree].prev].next = l->nextfree;   

	l->nextfree++;
	l->length++;
}

void
insertAfter(void *space, List *l, Uint cur, void *data) { 
	
    if (cur > l->length && cur != LISTNILVALUE)
	  return ;
	
   if (l->nextfree >= l->allocated-2) {
		l->nodes = ALLOCMEMORY(space, l->nodes,Listelem,l->allocated*2);
		l->allocated *= 2;
	}
	  
	l->nodes[l->nextfree].data = data;
	l->nodes[l->nextfree].prev = cur;
	
	if (cur != LISTNILVALUE) {
	  l->nodes[l->nextfree].next = l->nodes[cur].next;
	  l->nodes[cur].next = l->nextfree;	
	} else {
	  l->nodes[l->nextfree].next = LISTNILVALUE;
	}

	if (l->nodes[l->nextfree].next == LISTNILVALUE) 
	  l->lastNode = l->nextfree;
	else
	  l->nodes[l->nodes[l->nextfree].next].prev = l->nextfree;  
	
	l->nextfree++;
	l->length++;
}


/*-------------------------------- removeElem --------------------------------
 *    
 * removes an element from a list
 * but doesn't free
 * 
 */
 
void
unlinkListElem (List *l, Uint cur)
{
    if (l->nodes[cur].prev != LISTNILVALUE)
    	l->nodes[l->nodes[cur].prev].next = l->nodes[cur].next;
	else
	    l->firstNode = l->nodes[cur].next;
	
	if (l->nodes[cur].next != LISTNILVALUE)
		l->nodes[l->nodes[cur].next].prev = l->nodes[cur].prev;
	else
	  	l->lastNode = l->nodes[cur].prev;
	
	l->length--;
}


/*--------------------------------- wrapList ---------------------------------
 *    
 * destructs a list
 * 
 */
 
void
wrapList (void *space, List *l, 
		  void (*rmv)(void *, void *))
{
    Uint cur;

	if (rmv != NULL) {
		for (cur = l->firstNode; cur != LISTNILVALUE; cur=l->nodes[cur].next) 
		{
	   		rmv(space, l->nodes[cur].data);	 
		}
	}

	FREEMEMORY(space, l->nodes);
	FREEMEMORY(space, l);
}


/*--------------------------------- dataList ---------------------------------
 *    
 * returns a pointer to a sorted array of list data
 * 
 */
 
void**
dataList (void *space, List *l)
{
    void **data; 
	Uint cur, k=0;

    data = ALLOCMEMORY(space, NULL, void*, l->length);
	
	for (cur = l->firstNode; cur != LISTNILVALUE; cur = l->nodes[cur].next) {
        data[k] = l->nodes[cur].data; 
		k++;
	}

	return data;
}

/*-------------------------------- sweepList ---------------------------------
 *    
 * cleans a list of all unlinked elements and sorts the nodes 
 * array implicitly
 * 
 */
 
void sweepList(void *space, List *l) {
	Uint cur, k=0;
	Listelem *buf;
	 
	buf = ALLOCMEMORY(space, NULL, Listelem, l->length+1);
	
	for (cur = l->firstNode; cur != LISTNILVALUE; cur = l->nodes[cur].next) {
	    buf[k].data = l->nodes[cur].data;
		
		if (k > 0)
			buf[k].prev = k-1;
		else
		    buf[k].prev = LISTNILVALUE;
		
		if (l->nodes[cur].next != LISTNILVALUE)
			buf[k].next = k+1;
		else
		    buf[k].next = LISTNILVALUE;	
		
		k++;
	}

	FREEMEMORY(space, l->nodes);

	l->nodes = buf;
	l->firstNode = 0;
	l->lastNode  = l->length-1;
	l->allocated = l->length+1;
	l->nextfree  = l->length;
}

