/*
	  This file is part of gdub.
	  (C) 2006 Steve Hoffmann 
 
	  gdub is free software; you can redistribute it and/or modify
	  it under the terms of the GNU General Public License as published
	  by the Free Software Foundation; either version 2, or (at your
	  option) any later version.
 
	  gdub is distributed in the hope that it will be useful, but
	  WITHOUT ANY WARRANTY; without even the implied warranty of
	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	  General Public License for more details.
 
	  You should have received a copy of the GNU General Public License
	  along with gdub; see the file COPYING.  If not, write to the
	  Free Software Foundation, Inc., 59 Temple Place - Suite 330,
	  Boston, MA 02111-1307, USA.	
 
 */

/**
 * @file stack.c
 * @author Steve Hoffmann
 * @brief implementation of a simple stack
 */

/*
 * $Log$
 *
 */


 #include <stdio.h>
 #include <stdlib.h>
 #include "stack.h"


  void initStack(void *spacetab, Stack *stack, int size) {
 
	stack->stackelements = ALLOCMEMORY(spacetab, NULL, int, size);
	stack->size=size;
	stack->top=-1;

 }

  int stackisempty(Stack *stack) {
		return (stack->top < 0);
  }
		
  void stackpush(void* spacetab, Stack *stack, Stackelement elem) {
		
		if(stack->top >= stack->size-1) {
			
			stack->stackelements = ALLOCMEMORY(spacetab, stack->stackelements, Stackelement, ((stack->size+1)+STACKINCREMENT));
			stack->size=(stack->size+1)+STACKINCREMENT; 
		}

		stack->top++;
		stack->stackelements[stack->top] = elem;
  }

 Stackelement stackpop(Stack *stack){
 

	if(!stackisempty(stack)) {
		return stack->stackelements[stack->top--];	
	}
	
	return STACK_NULL_TYPE;
 }

 Stackelement stacktop(Stack *stack){
	
	if(!stackisempty(stack)) {
		return stack->stackelements[stack->top];
	}

	return STACK_NULL_TYPE;
 }
 
 Stackelement stacktopn(Stack *stack, Uint n){
	
	if(!stackisempty(stack) && stack->top >= n) {
		return stack->stackelements[stack->top-n];
	}

	return STACK_NULL_TYPE;
 }

 void destructStack(void *spacetab, Stack *stack) {
	FREEMEMORY(spacetab, stack->stackelements);
 }
