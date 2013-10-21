/*
	Object:		Raster-scan and label-equivalence-based algorithm.
	Authors:	Massimo Nicolazzo & Giuliano Langella
	email:		gyuliano@libero.it
*/

/*
|------------------------------------|
|                  |                 |
|        I        2|1     II         |      20001=10002         I  II III IV
|      10000      0|0   20000        |      20001=10004		1   1   1   0  1
|                 4|1                |      30002=10004     2   2   1   0  2
|                 4|                 |                      3   1   1
|------------------------------------|                      4   4   1
|                 2|                 |					 2500	0
|                  |                 |
|       III        |      IV         |
|                  |                 |
|                  |                 |
|------------------------------------|

how many rows: (ntiles * 2 + ntilesX + ntilesY) *  ((tiledimX+tiledimY)/2)

*/


/*
 *
--------------+
0 1 2 3 4 5  6 |7 8 9 10   11| 12 13 14 15
            |----------------|


0 0 0 0 1 0  0 0 0 0 1
1 1 0 0 1 1  1 1 0 0 1
1 0 1 0 0 1  1 0 1 0 0
0 0 0 0 1 0  0 0 0 0 1
0 0 0 0 0 0  0 0 0 0 0
0 0 0 0 1 0


s_1
0 0 0 0 1 0  0 0 0 0 1		if(1) --> equiv(t1,t2)
2 2 0 0 1 1  2 2 0 0 1
2 0 2 0 0 1  2 0 2 0 0
0 0 0 0 1 0  0 0 0 0 3
0 0 0 0 0 0  0 0 0 0 0


mc1=2, mc2=2 => 6 x 2

  t1 t2   t1*mc2+t2
0  0  0      0
1  0  1      1
2  1  0      2
4  1  1      3



s_1-nt
0 0 0 0 1 0  0 0 0 0 1		if(1) --> equiv(t1,t2)
2 2 0 0 1 1  2 2 0 0 1
2 0 2 0 0 0  0 2 2 0 0
0 0 0 0 3 3  2 0 0 0 3
0 0 0 0 0 0  0 0 0 0 0

s_2
 0  0  0  0 x1    0  0  0 0 y1
x2 x2  0  0 x1   x1 x1  0 0 y1
x2  0 x2  0  0   x1  0 x1 0  0
 0  0  0  0 x1    0  0  0 0 y3
 0  0  0  0  0    0  0  0 0  0

s_3
s_1
0 0 0 0 1  0 0 0 0 3
2 2 0 0 1  1 1 0 0 3
2 0 2 0 0  1 0 1 0 0
0 0 0 0 1  0 0 0 0 4
0 0 0 0 0  0 0 0 0 0


 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

// I: "urban" 	--> [0,0] shifted
// O: "lab_mat" --> [1,1] shifted
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
#define nw_pol(c,r)		lab_mat[	(c)		+	(r)		*(ncols+1)	] // O: scan value at North-West
#define nn_pol(c,r)		lab_mat[	(c+1)	+	(r)		*(ncols+1)	] // O: scan value at North
#define ne_pol(c,r)		lab_mat[	(c+2)	+	(r)		*(ncols+1)	] // O: scan value at North-East
#define ww_pol(c,r)		lab_mat[	(c)		+	(r+1)	*(ncols+1)	] // O: scan value at West
#define cc_pol(c,r)		lab_mat[	(c+1)	+	(r+1)	*(ncols+1)	] // O: scan value at current [r,c] which is shifted by [1,1] in O

// MIN & MAX
#define max(val1,val2)		(val1)>(val2)?(val1):(val2)
#define min(val1,val2)		(val1)<(val2)?(val1):(val2)

unsigned char 	Vb			= 0;	// background value
unsigned char 	Vo			= 1;	// object value
unsigned int	cross_cols	= 4;	// number of columns of the cross_parent matrix

//----------------------------

/*[]={
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
		*/
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

0.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0.0,1,0,1,1,1,1,1,1,1,1,1,1,1,0,
0.0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,
0.0,1,0,1,0,1,1,1,1,1,1,1,0,1,0,
0.0,1,0,1,0,1,0,0,0,0,0,1,0,1,0,
0.0,1,0,1,0,1,1,1,1,1,0,1,0,1,0,
0.1,1,0,1,0,0,0,0,0,0,0,1,0,1,0,
0.0,1,0,1,1,1,1,1,1,1,1,1,0,1,0,
0.0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,
0.0,1,1,1,1,1,1,1,1,1,1,1,1,1,0,
0.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0

*/

//----------------------------

unsigned int first_scan(	unsigned char	*urban,
							unsigned int	nrows,
							unsigned int	ncols,
							unsigned int	*lab_mat,
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

			if(durban(c,r)==Vo) // if (r,c) is object pixel 
			{
				/*
				 * 	We use the so called "forward scan mask"
				 * 	in which only four adjacent pixels are considered:
				 * 		{ ww, nw, nn, ne }
				 */
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
	if(val1!=val2) PARENT[max(val1,val2)] = min(val1,val2);
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
		unsigned int	*lab_mat,
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

void read_mat(unsigned char *urban, unsigned int nrows, unsigned int ncols, char *filename)
{
	unsigned int rr,cc;
	FILE *fid ;
	int a;
	fid= fopen(filename,"rt");
	for(rr=0;rr<nrows;rr++) for(cc=0;cc<ncols;cc++) { fscanf(fid, "%d",&a);	durban(cc,rr)=(unsigned char)a;}
	fclose(fid);
}

// to be developed [MAYBE]:
void record_cross_equivalence()
{
	/* Record in cross_parent in the first rows the root equivalence for the first
	 * tile (i.e. tile[0] ) of the grid of tiles. That is I have to write the following two times:
	 * 		{ parent tile number, parent ID of parent tile number }
	 * having:
	 * 		{ parent tile number, parent ID of parent tile number, parent tile number, parent ID of parent tile number }
	 * which is able to give the ROOT to the union_cross_equivalence function.
	 */

	// here I must call the objects_stintching_XX functions!
	// first:
	//objects_stitching_nn
	// second:
	//objects_stitching_ww
	// third [maybe not useful]:
	//objects_stitching_nw
}

void objects_stitching_ww(
		unsigned int *lm_ww,		//label matrix of western tile in the mask
		unsigned int *lm_cc,		//label matrix of target (=centre) tile in the mask
		unsigned int nr,			// number of rows
		unsigned int nc,			// number of columns
		unsigned int ntile_ww,		// number of ww tile in the grid of tiles
		unsigned int ntile_cc,		// number of cc tile in the grid of tiles
		unsigned int *cross_parent,	// pointer to the first row in cross_parent
		unsigned int *first_empty	// the first empty row in cross_parent SCALAR
						)
{
	unsigned int r;
	for(r=0;r<nr;r++)
	{
		if (lm_ww(nc*(r+1)-1)!=0)
		{
			cross_parent[&first_empty*cross_cols+0]	= ntile_cc;				// parent tile number
			cross_parent[&first_empty*cross_cols+1]	= lm_cc(nc*r);			// parent ID of parent tile number
			cross_parent[&first_empty*cross_cols+2]	= ntile_ww;				// root tile number
			cross_parent[&first_empty*cross_cols+3]	= lm_ww(nc*(r+1)-1);	// root ID of root tile number
			first_empty +=1;
		}
	}
}

void objects_stitching_nn(
		unsigned int *lm_nn,		//label matrix of northern tile in the mask
		unsigned int *lm_cc,		//label matrix of target (=centre) tile in the mask
		unsigned int nr,			// number of rows
		unsigned int nc,			// number of columns
		unsigned int ntile_nn,		// number of nn tile in the grid of tiles
		unsigned int ntile_cc,		// number of cc tile in the grid of tiles
		unsigned int *cross_parent,	// pointer to the first row in cross_parent
		unsigned int *first_empty	// the first empty row in cross_parent SCALAR
						)
{
	unsigned int c;
	for(c=0;c<nc;c++)
	{
		if (lm_nn(nc*(nr-1)+c)!=0)
		{
			cross_parent[&first_empty*cross_cols+0]	= ntile_cc;				// parent tile number
			cross_parent[&first_empty*cross_cols+1]	= lm_cc(c);				// parent ID of parent tile number
			cross_parent[&first_empty*cross_cols+2]	= ntile_nn;				// root tile number
			cross_parent[&first_empty*cross_cols+3]	= lm_nn(nc*(nr-1)+c);	// root ID of root tile number
			first_empty +=1;
		}
	}
}

void objects_stitching_nw(
		unsigned int *lm_nw,		//label matrix of north-western tile in the mask
		unsigned int *lm_cc,		//label matrix of target (=centre) tile in the mask
		unsigned int nr,			// number of rows
		unsigned int nc,			// number of columns
		unsigned int ntile_nw,		// number of nw tile in the grid of tiles
		unsigned int ntile_cc,		// number of cc tile in the grid of tiles
		unsigned int *cross_parent,	// pointer to the first row in cross_parent
		unsigned int *first_empty	// the first empty row in cross_parent SCALAR
						)
{
	if (lm_nw(nc*nr-1)!=0)
	{
		cross_parent[&first_empty*cross_cols+0]	= ntile_cc;				// parent tile number
		cross_parent[&first_empty*cross_cols+1]	= lm_cc(0);				// parent ID of parent tile number
		cross_parent[&first_empty*cross_cols+2]	= ntile_nw;				// root tile number
		cross_parent[&first_empty*cross_cols+3]	= lm_nw(nc*nr-1);		// root ID of root tile number
		first_empty +=1;
	}
}

void union_cross_equivalence(unsigned int first_empty, unsigned int *cross_parent)
{
	unsigned int i,j;
	// I have to check whether one for loop suffices to solve the whole ROOT tree
	for(i=0;i<first_empty;first_empty++)
	{
		j = 0;
		// DO (j=j+1) Until (I will find the (first occurrence of) ROOT cross equivalence row recorded in cross_parent):
		while( (cross_parent[j*cross_cols+0]!=cross_parent[i*cross_cols+2]) && (cross_parent[j*cross_cols+1]!=cross_parent[i*cross_cols+3]) && j<i) j+=1;
		if(j==i) // the while loop should return j=i if any previous ROOT was found!
		{
			cross_parent[i*cross_cols+2] = cross_parent[j*cross_cols+0];
			cross_parent[i*cross_cols+3] = cross_parent[j*cross_cols+1];
		}
	}
}

int main()
{
	// DECLARATION:
	unsigned int nTiles = 2,iTile;
	unsigned int * (pol[nTiles]);
	unsigned int * (cont[nTiles]);
	unsigned int i,	mc[nTiles];
	unsigned int *(PARENT[nTiles]);
	unsigned int dc		= 6;
	unsigned int dr		= 5;
	unsigned int nID	= dc*dr/4 +1;
	unsigned char *(urban[nTiles]);
	char s[255];

	for(iTile = 0;iTile<nTiles;iTile++)
	{
		// INITIALIZATION:
		pol[iTile]		= (unsigned int*)calloc((dc+1)*(dr+1),sizeof(unsigned int));// la pol deve avere 1^ riga e colonna di zeri in più per consentire loop (i,j):
		cont[iTile]		= (unsigned int*)malloc(nID*sizeof(unsigned int));
		PARENT[iTile]	= (unsigned int*)calloc(nID,sizeof(unsigned int));
		urban[iTile] 	= (unsigned char*)calloc((dc)*(dr),sizeof(unsigned int));

		// READ
		sprintf(s,"/home/giuliano/git/soil-sealing/data/%c.txt",65+iTile);
		printf("here!\n");
		puts(s);
		read_mat(urban[iTile], dr, dc, s);

		// KERNELs INVOCATION:
		mc[iTile]		= first_scan(urban[iTile],dr,dc,pol[iTile],cont[iTile],PARENT[iTile]);	//	(1) 1st SCAN
		union_equivalence(		mc[iTile], PARENT[iTile]);				//	(2) UNION
		relabel_equivalence(	mc[iTile], PARENT[iTile]);				//	(3) RELABEL
		// ??												//	(?)	COMPACT ==> no gaps in
		second_scan(pol[iTile],dr,dc, PARENT[iTile]);						//	(4) 2nd SCAN

	}
	/*
	// PRINTS:
	print_mat(	urban,	dr,		dc,		"Binary Raster"							);
	print_mat(	pol,	dr+1,	dc+1, 	"CCL (Connected-Components Labeling)"	);
	printf("Indexes:\n"); for(i=0;i<=mc;i++) {printf("%3.0d ",i);};  printf("\n");
	print_vec(	PARENT,	mc+1,			"Labels"								);
	print_vec(	cont,	mc+1,			"Counts"								);
*/




	// FREE MEMORY:
	for(iTile = 0;iTile<nTiles;iTile++)
	{
		free(cont[iTile]);
		free(pol[iTile]);
		free(PARENT[iTile]);
		free(urban[iTile]);
	}
	
	// RETURN:
	return 0;
}
