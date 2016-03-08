
/*
 *  vtprogressbar.c
 *  implementation for a very simple
 *  progress bar
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/05/06 02:11:30 CET
 *  
 */
 #include <stdio.h>
 #include <stdlib.h>
 #include "basic-types.h"
 #include "vtprogressbar.h"


void
cursorInvisible() {
  fprintf(stderr, "%c%c%c%d%c", 27, '[', '?', 25, 'l');
}

void
cursorVisible() {
  fprintf(stderr, "%c%c%c%d%c", 27, '[', '?', 25, 'h');
}



/*---------------------------- initProgressBarVT -----------------------------
 *    
 * initializes a progress bar for VT
 * 
 */
 
void
initProgressBarVT ()
{	
  	fprintf(stderr, "%c%c%c", 27, '[', 's');
    fprintf(stderr, "%c%c%c", 27, '[', 'K');
    return ;
}

/*------------------------------ progressBarVT -------------------------------
 *    
 * a simple progress bar for VT terminals
 * 
 */

void
progressBarVT (char *message, Uint complete, Uint processed, Uint size)
{
    Uint i, percent, bar;
	char cur;
    
    if (complete == 0) complete = 1;
	bar=(size*processed)/(complete);
	percent =(processed*100)/(complete);
	fprintf(stderr, "[");
	for(i=0; i < size; i++) {
		if(i<=bar) fprintf(stderr, "=");
		else fprintf(stderr," ");
	}
	i = processed % 30;
	if (i<=10) cur = '/'; 
	else if (i<=20)cur='\\';
	else cur='-';
	fprintf(stderr,"]   %d%c   %s   %c\n", percent, '%', message, cur);
	fprintf(stderr,"%c%c%c", 27, '[', 'A');
	return;
}


