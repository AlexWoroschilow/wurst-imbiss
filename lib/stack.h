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
 * @file stack.h
 * @author Steve Hoffmann
 * @brief header file for a simple stack
 */

/*
 * $Log$
 *
 */

#ifndef STACK_H
#define STACK_H

#include <stdio.h>
#include <stdlib.h>
#include "basic-types.h"
#include "memory.h"

#define STACK_NULL_TYPE 0
#define STACKINCREMENT 10000

typedef int Stackelement;


typedef struct{
	Stackelement* stackelements;
	int size;
	int top;
} Stack;

void  initStack(void *spacetab, Stack* stack, int size);
int  stackisempty(Stack *stack);
void stackpush(void* spacetab, Stack *stack, Stackelement elem);
Stackelement stackpop(Stack *stack);
Stackelement stacktop(Stack *stack);
void destructStack(void* spacetab, Stack *stack);
Stackelement stacktopn(Stack *stack, Uint n);

#endif
