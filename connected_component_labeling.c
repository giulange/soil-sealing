/*
	Object:		Raster-scan and label-equivalence-based algorithm.
	Authors:	Massimo Nicolazzo & Giuliano Langella
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

// I: "urban" 	--> [0,0] shifted
// O: "polygon" --> [1,1] shifted
/*
	The "forward scan mask" for eight connected connectivity is the following:
		nw		nn		ne
		ww		cc		xx
		xx		xx		xx
	assuming that:
		> cc is is the background(=0)/foreground(=1) pixel at (r,c),
		> nw, nn, ne, ww are the north-west, north, north-east and west pixels in the eight connected connectivity,
		> xx are skipped pixels.
	Therefore the mask has 4 active pixels with(out) object pixels (that is foreground pixels).
*/
#define durban(c,r)		urban[		(c)		+	(r)		*ncols		] // I: scan value at current [r,c] 
#define nw_pol(c,r)		polygon[	(c)		+	(r)		*(ncols+1)	] // O: scan value at North-West
#define nn_pol(c,r)		polygon[	(c+1)	+	(r)		*(ncols+1)	] // O: scan value at North
#define ne_pol(c,r)		polygon[	(c+2)	+	(r)		*(ncols+1)	] // O: scan value at North-East
#define ww_pol(c,r)		polygon[	(c)		+	(r+1)	*(ncols+1)	] // O: scan value at West
#define cc_pol(c,r)		polygon[	(c+1)	+	(r+1)	*(ncols+1)	] // O: scan value at current [r,c] which is shifted by [1,1] in O

// MIN & MAX
#define max(val1,val2)		(val1)>(val2)?(val1):(val2)
#define min(val1,val2)		(val1)<(val2)?(val1):(val2)

unsigned char Vb	= 0;
unsigned char Vo	= 1;

//----------------------------
// to-be-deleted GLOBALS:
unsigned char 

Urban[]={
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,
0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,
0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,
1,1,0,0,1,0,0,0,0,0,1,0,1,0,0,
1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,
0,1,1,1,0,0,0,0,0,0,1,1,1,1,1,
1,0,0,0,0,0,0,1,1,0,1,1,1,0,0,
0,0,0,0,0,0,1,1,1,0,1,1,0,1,1,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
		};
/*
Urban[]={
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,1,0,1,1,1,1,1,1,1,1,1,1,1,0,
0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,
0,1,0,1,0,1,1,1,1,1,1,1,0,1,0,
0,1,0,1,0,1,0,0,0,0,0,1,0,1,0,
0,1,0,1,0,1,1,1,1,1,0,1,0,1,0,
1,1,0,1,0,0,0,0,0,0,0,1,0,1,0,
0,1,0,1,1,1,1,1,1,1,1,1,0,1,0,
0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,
0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
		};
*/
unsigned int dc		= 15;
unsigned int dr		= 11;
unsigned int nID	= 300;
//----------------------------

unsigned int first_scan(	unsigned char	*urban,
							unsigned int	nrows,
							unsigned int	ncols,
							unsigned char	*polygon,
							unsigned int	*count,
							unsigned int	*PARENT		)
{
	unsigned int r			= 0;
	unsigned int c			= 0;
	unsigned int maxcount	= 0;	
	for(r=0; r<nrows; r++)
	{
		for(c=0; c<ncols; c++)
		{
			/* IMPORTANT !!!!!!!!!!
				Multiply every assignment xx_pol(c,r) by durban(c,r) to avoid the following if.
				This is important when passing to CUDA!
				But take care that the input raster must be binary: 0:=BackGround; 1:=Object).
			*/
			if(durban(c,r)==Vo) // if (r,c) is object pixel 
			{
				// NORTH:
				if(nn_pol(c,r)!=Vb) cc_pol(c,r) = nn_pol(c,r);
				
				// WEST:
				else if(ww_pol(c,r)!=Vb)
				{
					cc_pol(c,r) = ww_pol(c,r);
					if (ne_pol(c,r)!=Vb) record_equivalence( ne_pol(c,r), ww_pol(c,r), PARENT );
				}
				
				// NORTH-WEST:
				else if(nw_pol(c,r)!=Vb)
				{
					cc_pol(c,r) = nw_pol(c,r);
					if (ne_pol(c,r)!=Vb) record_equivalence( ne_pol(c,r), nw_pol(c,r), PARENT );
				}
				
				// NORTH-EAST:
				else if(ne_pol(c,r)!=Vb) cc_pol(c,r) = ne_pol(c,r);
				
				//none object pixels in mask:
				else cc_pol(c,r) = ++maxcount;
				
				// I am not sure that we should count right NOW.
				// Maybe It could be better to performe this is a separate phase when all 
				// equivalent labels are already known and applied.
				count[cc_pol(c,r)]++;
			}
		}
	}	
	return maxcount;
}

record_equivalence(unsigned int val1, unsigned int val2, unsigned int *PARENT)
{
	PARENT[max(val1,val2)] = min(val1,val2);
}

union_equivalence(unsigned int maxcount, unsigned int *PARENT)
{
	unsigned int j,k;
	for(j=1;j<=maxcount;j++)
	{
		if(PARENT[j]==0) continue;
		k = j;
		while(PARENT[k] != 0) k = PARENT[k];
		PARENT[j] = k;
	}
}

relabel_equivalence(unsigned int maxcount, unsigned int *PARENT)
{
	unsigned int j;
	for(j=1;j<=maxcount;j++) if(PARENT[j]==0) PARENT[j] = j;
}

second_scan(	
		unsigned char	*polygon,
		unsigned int	nrows,
		unsigned int	ncols,
		unsigned int	*PARENT		)
{
	unsigned int i,j;
    for(i=0;i<nrows;i++) for(j=0;j<ncols;j++) cc_pol(j,i) = PARENT[cc_pol(j,i)];
}

void print_mat(unsigned char *u,unsigned int nrows,unsigned int ncols, char *Label)
{
	unsigned int i,j;
	printf("\n%s:\n", Label);
	for(i=0;i<nrows;i++)
	{
		for(j=0;j<ncols;j++) printf("%d ",u[(j)+(i)*ncols]);
		printf("\n");
	}
	printf("\n");
}

void print_vec( unsigned int *vec, unsigned int numel, unsigned char *Label )
{
	unsigned int i;
	printf("\n%s:\n", Label);
	for(i=0;i<numel;i++) printf("%3.0d ",vec[i]);
	printf("\n");
}


int main()
{
	// DECLARATION:
	unsigned char * pol;
	unsigned int * cont,*rcont;
	unsigned int i,mc;
	unsigned int *PARENT;
	
	// INITIALIZATION:
	pol		= (unsigned char*)calloc((dc+1)*(dr+1),sizeof(unsigned char));// la pol deve avere 1^ riga e colonna di zeri in più per consentire loop (i,j):
	cont 	= (unsigned int*)malloc(nID*sizeof(unsigned int));
	PARENT	= (unsigned int*)calloc(nID,sizeof(unsigned int));
	
	// KERNELs INVOCATION:
	mc		= first_scan(Urban,dr,dc,pol,cont,PARENT);	//	(1) 1st SCAN
	union_equivalence(		mc, PARENT);				//	(2) UNION
	relabel_equivalence(	mc, PARENT);				//	(3) RELABEL
	// ??												//	(?)	COMPACT ==> no gaps in 
	second_scan(pol,dr,dc, PARENT);						//	(4) 2nd SCAN
	
	// PRINTS:
	print_mat(	Urban,	dr,		dc,		"Binary Raster"							);
	print_mat(	pol,	dr+1,	dc+1, 	"CCL (Connected-Components Labeling)"	);
	printf("Indexes:\n"); for(i=0;i<=mc;i++) {printf("%3.0d ",i);};  printf("\n");
	print_vec(	PARENT,	mc+1,			"Labels"								);
	print_vec(	cont,	mc+1,			"Counts"								);

	// FREE MEMORY:
	free(cont);
	free(pol);
	free(PARENT);
	
	// RETURN:
	return 0;
}