/*
	Object:		Raster-scan and label-equivalence-based algorithm.
	Authors:	Massimo Nicolazzo & Giuliano Langella
	email:		gyuliano@libero.it
*/

/*
how many rows: (ntiles * 2 + ntilesX + ntilesY) *  ((tiledimX+tiledimY)/2)
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <errno.h>        /* errno */
#include <string.h>       /* strerror */

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
#define durban(c,r)		urban[		(c)	+	(r)	*ncols		] // I: scan value at current [r,c] 
#define nw_pol(c,r)		lab_mat[	(c-1)	+	(r-1)	*(ncols)	] // O: scan value at North-West
#define nn_pol(c,r)		lab_mat[	(c+0)	+	(r+0)	*(ncols)	] // O: scan value at North
#define ne_pol(c,r)		lab_mat[	(c+1)	+	(r-1)	*(ncols)	] // O: scan value at North-East
#define ww_pol(c,r)		lab_mat[	(c-1)	+	(r+0)	*(ncols)	] // O: scan value at West
#define cc_pol(c,r)		lab_mat[	(c+0)	+	(r+0)	*(ncols)	] // O: scan value at current [r,c] which is shifted by [1,1] in O
// MIN & MAX
#define max(val1,val2)		(val1)>(val2)?(val1):(val2)
#define min(val1,val2)		(val1)<(val2)?(val1):(val2)

extern int	errno;
unsigned char 	Vb			= 0;	// background value
unsigned char 	Vo			= 1;	// object value
unsigned int	cross_cols	= 4;	// number of columns of the cross_parent matrix

//---------------------------- FUNCTIONS PROTOTYPES
// 	FIRST STAGE
unsigned int first_scan( unsigned char *urban, unsigned int nrows, unsigned int ncols,unsigned int *lab_mat,unsigned int *count,unsigned int *PARENT);
void record_equivalence(unsigned int val1, unsigned int val2, unsigned int *PARENT);
void union_equivalence(unsigned int maxcount, unsigned int *PARENT);
void relabel_equivalence(unsigned int maxcount, unsigned int *PARENT);
void second_scan( unsigned int *lab_mat,unsigned int nrows,unsigned int ncols,unsigned int *PARENT );
void print_mat(unsigned char *u,unsigned int nrows,unsigned int ncols, char *Label);
void print_vec( unsigned int *vec, unsigned int numel, unsigned char *Label );
void read_mat(unsigned char *urban, unsigned int nrows, unsigned int ncols, char *filename);

// 	SECOND STAGE
unsigned int *objects_stitching_nn(unsigned int *lm_nn,unsigned int *lm_cc,unsigned int nr,unsigned int nc,unsigned int ntile_nn,unsigned int ntile_cc,unsigned int *cross_parent);
unsigned int *objects_stitching_ww(unsigned int *lm_ww,unsigned int *lm_cc,unsigned int nr,unsigned int nc,unsigned int ntile_ww,unsigned int ntile_cc,unsigned int *cross_parent);
unsigned int *objects_stitching_cc(unsigned int *lm_cc,unsigned int nr,unsigned int nc,unsigned int ntile_cc,unsigned int *cross_parent,unsigned int *cross_parent_ii,unsigned int first_empty,unsigned int *PARENT,unsigned int maxcount);
unsigned int *record_cross_equivalence(unsigned int **lm,unsigned int *cross_parent,unsigned int nr,unsigned int nc,unsigned int ntile_cc,int ntile_nn,int ntile_ww,unsigned int *PARENT,unsigned int mc);
void union_cross_equivalence(unsigned int first_empty, unsigned int *cross_parent);
void relabel_cross_equivalence(unsigned int first_empty, unsigned int *cross_parent);
//---------------------------- FUNCTIONS PROTOTYPES

unsigned int first_scan(	unsigned char	*urban,
				unsigned int	nrows,
				unsigned int	ncols,
				unsigned int	*lab_mat,
				unsigned int	*count,
				unsigned int	*PARENT		)
{
	unsigned int r		= 0;
	unsigned int c		= 0;
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
				if(r>0 && nn_pol(c,r)!=Vb) cc_pol(c,r) = nn_pol(c,r);
				
				// WEST:
				else if(c>0 && ww_pol(c,r)!=Vb)
				{
					cc_pol(c,r) = ww_pol(c,r);
					if (ne_pol(c,r)!=Vb) record_equivalence( ne_pol(c,r), ww_pol(c,r), PARENT );
				}
				
				// NORTH-WEST:
				else if(r>0 && c>0 && nw_pol(c,r)!=Vb)
				{
					cc_pol(c,r) = nw_pol(c,r);
					if (ne_pol(c,r)!=Vb) record_equivalence( ne_pol(c,r), nw_pol(c,r), PARENT );
				}
				
				// NORTH-EAST:
				else if(r>0 && c<ncols && ne_pol(c,r)!=Vb) cc_pol(c,r) = ne_pol(c,r);
				
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

void record_equivalence(unsigned int val1, unsigned int val2, unsigned int *PARENT)
{
	if(val1!=val2) PARENT[max(val1,val2)] = min(val1,val2);
}

void union_equivalence(unsigned int maxcount, unsigned int *PARENT)
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

void relabel_equivalence(unsigned int maxcount, unsigned int *PARENT)
{
	/*
	 * 		Assign to PARENT[j] j itself only if it is ROOT id.
	 * 		This is useful during the second_scan because we avoid
	 * 		an if statement.
	 */
	unsigned int j;
	for(j=1;j<=maxcount;j++) if(PARENT[j]==0) PARENT[j] = j;
}

void second_scan(
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

unsigned int *objects_stitching_nn(
		unsigned int *lm_nn,		// label matrix of northern tile in the mask
		unsigned int *lm_cc,		// label matrix of target (=centre) tile in the mask
		unsigned int nr,		// number of rows
		unsigned int nc,		// number of columns
		unsigned int ntile_nn,		// number of nn tile in the grid of tiles
		unsigned int ntile_cc,		// number of cc tile in the grid of tiles
		unsigned int *cross_parent	// pointer to the first row in cross_parent
						)
{
	unsigned int c;
	for(c=0;c<nc;c++)
	{
		if (lm_nn[nc*(nr-1)+c]!=0) // usare la lm_cc
		{
			*cross_parent++	= ntile_cc;		// parent tile number
			*cross_parent++	= lm_cc[c];		// parent ID of parent tile number
			*cross_parent++	= ntile_nn;		// root tile number
			*cross_parent++	= lm_nn[nc*(nr-1)+c];	// root ID of root tile number
		}
	}
	return cross_parent;
}

unsigned int *objects_stitching_ww(
		unsigned int *lm_ww,		// label matrix of western tile in the mask
		unsigned int *lm_cc,		// label matrix of target (=centre) tile in the mask
		unsigned int nr,		// number of rows
		unsigned int nc,		// number of columns
		unsigned int ntile_ww,		// number of ww tile in the grid of tiles
		unsigned int ntile_cc,		// number of cc tile in the grid of tiles
		unsigned int *cross_parent	// pointer to the first row in cross_parent
						)
{
	unsigned int r;
	for(r=0;r<nr;r++)
	{
		if (lm_ww[nc*(r+1)-1]!=0)
		{
			*cross_parent++	= ntile_cc;		// parent tile number
			*cross_parent++	= lm_cc[nc*r];		// parent ID of parent tile number
			*cross_parent++	= ntile_ww;		// root tile number
			*cross_parent++	= lm_ww[nc*(r+1)-1];	// root ID of root tile number
		}
	}
	return cross_parent;
}

unsigned int *objects_stitching_cc(
		unsigned int *lm_cc,		//label matrix of target (=centre) tile in the mask
		unsigned int nr,		// number of rows
		unsigned int nc,		// number of columns
		unsigned int ntile_cc,		// number of cc tile in the grid of tiles
		unsigned int *cross_parent,	// pointer to the first row in cross_parent
		unsigned int *cross_parent_ii,	// pointer to the first element before stitching by {nn,ww}
		unsigned int first_empty,	// the first empty row in cross_parent SCALAR
		unsigned int *PARENT,		// the PARENT vector of target tile within tiles-mask
		unsigned int maxcount		// number of labels in target tile stored in PARENT
						)
{
	unsigned int i,j;
	unsigned char found=0;
	for(i=1;i<=maxcount;i++) // start from 1, because 0 is background
	{
		found = 0;
		for(j=0;j<first_empty;j++) if( (cross_parent_ii[j*cross_cols+0]==ntile_cc) && (cross_parent_ii[j*cross_cols+1]==PARENT[i]) ) {found = 1; break;}

		if (found==0) // then record current label as ROOT label! (i.e. elements {1,2} are equal to elements {3,4} of cross_parent
		{
			*cross_parent++	= ntile_cc;	// parent tile number
			*cross_parent++	= PARENT[i];	// parent ID of parent tile number
			*cross_parent++	= ntile_cc;	// root tile number
			*cross_parent++	= PARENT[i];	// root ID of root tile number
		}
	}
	return cross_parent;
}
unsigned int *record_cross_equivalence(
		unsigned int **lm,
		unsigned int *cross_parent,
		unsigned int nr,
		unsigned int nc,
		unsigned int ntile_cc,
		int ntile_nn,
		int ntile_ww,
		unsigned int *PARENT,
		unsigned int mc)
{
	unsigned int *lm_cc;
	unsigned int *lm_nn;
	unsigned int *lm_ww;
	lm_cc=lm[ntile_cc];
	/*
	 * For every tile, record equivalences in cross_parent in following order:
	 *
	 *	(1) nn [only if existent] 	---> objects_stitching_nn
	 *	(2) ww [only if existent]	---> objects_stitching_ww
	 *	(3) cc [always  existent]	---> objects_stitching_cc
	 *
	 * Point (3) records non-crossing objects and then it can build the complete and continuous
	 * (i.e. without jumping) set of labels for the whole image.
	 *
	 */
	unsigned int *cross_parent_ii,first_empty;
	cross_parent_ii = cross_parent;
	// here I must call the objects_stintching_XX functions!
	// (1)
	//		objects_stitching_nn()
	if( ntile_nn>=0 ) {
		lm_nn=lm[ntile_nn];
		cross_parent=objects_stitching_nn(lm_nn,lm_cc, nr, nc,ntile_nn, ntile_cc, cross_parent);
	}
	// (2)
	if( ntile_ww>=0 ) {
		lm_ww=lm[ntile_ww];
		cross_parent=objects_stitching_ww(lm_ww,lm_cc, nr, nc,ntile_ww, ntile_cc, cross_parent);
	}
	first_empty = cross_parent-cross_parent_ii;
	// (3)
	cross_parent=objects_stitching_cc(lm_cc,nr,nc,ntile_cc,cross_parent, cross_parent_ii, first_empty, PARENT,mc);
	return cross_parent;
}

void union_cross_equivalence(unsigned int first_empty, unsigned int *cross_parent)
{
	unsigned int i,j;
	// I have to check whether one for loop suffices to solve the whole ROOT tree
	for(i=0;i<first_empty;first_empty++)
	{
		j = 0;
		// DO (j=j+1) Until (I will find the (first occurrence of) ROOT cross equivalence row recorded in cross_parent):
		while( (cross_parent[j*cross_cols+0]!=cross_parent[i*cross_cols+2]) && (cross_parent[j*cross_cols+1]!=cross_parent[i*cross_cols+3]) && j<=i) j+=1;
		/*
		 *	The while loop returns j=i if any previous ROOT was found!
		 *	Hence if j<i I have to update equivalence changing current label with ROOT label.
		 */
		if(j<i)
		{
			cross_parent[i*cross_cols+2] = cross_parent[j*cross_cols+2];
			cross_parent[i*cross_cols+3] = cross_parent[j*cross_cols+3];
		}
	}
}

void relabel_cross_equivalence(unsigned int first_empty, unsigned int *cross_parent)
{
	unsigned int *final_parent;
	final_parent = (unsigned int*)calloc(first_empty,sizeof(unsigned int));
	unsigned int i,j,final_count=0;
	// I have to check whether one for loop suffices to solve the whole ROOT tree
	for(i=0;i<first_empty;first_empty++)
	{
		if( (cross_parent[i*cross_cols+0]==cross_parent[i*cross_cols+2]) && (cross_parent[i*cross_cols+1]==cross_parent[i*cross_cols+3]) ) final_parent[i] = final_count++;
		else
		{
			j = 0;
			// DO (j=j+1) Until (I will find the (first occurrence of) ROOT cross equivalence row recorded in cross_parent):
			while( (cross_parent[j*cross_cols+0]!=cross_parent[i*cross_cols+2]) && (cross_parent[j*cross_cols+1]!=cross_parent[i*cross_cols+3]) && j<=i) j+=1;
			if(j<i) final_parent[i] = final_parent[j];
			else puts("Error in cross_parent!\nCheck one of {union_cross_equivalence, record_cross_equivalence}"); puts(strerror(errno));
		}
	}
}

int main()
{
	// DECLARATION:
	unsigned int nTiles,iTile;
	unsigned int ntilesX=2,ntilesY=2;
	nTiles = ntilesX*ntilesY;
	unsigned int * (lab_mat[nTiles]);
	unsigned int * (cont[nTiles]);
	unsigned int i,	mc[nTiles];
	unsigned int *(PARENT[nTiles]);
	unsigned int nc		= 6;
	unsigned int nr		= 5;
	unsigned int nID	= nc*nr/4 +1;
	unsigned char *(urban[nTiles]);
	char s[255];

	for(iTile = 0;iTile<nTiles;iTile++)
	{
		// INITIALIZATION:
		lab_mat[iTile]	= (unsigned int*)calloc((nc)*(nr),sizeof(unsigned int));// la lab_mat deve avere 1^ riga e colonna di zeri in più per consentire loop (i,j):
		cont[iTile]		= (unsigned int*)malloc(nID*sizeof(unsigned int));
		PARENT[iTile]	= (unsigned int*)calloc(nID,sizeof(unsigned int));
		urban[iTile] 	= (unsigned char*)calloc((nc)*(nr),sizeof(unsigned int));

		// READ
		sprintf(s,"/home/giuliano/git/soil-sealing/data/%c.txt",65+iTile);
		printf("here!\n");
		puts(s);
		read_mat(urban[iTile], nr, nc, s);

		// KERNELs INVOCATION:
		mc[iTile]		= first_scan(urban[iTile],nr,nc,lab_mat[iTile],cont[iTile],PARENT[iTile]);	//	(1) 1st SCAN
		union_equivalence(		mc[iTile], PARENT[iTile]);				//	(2) UNION
		relabel_equivalence(	mc[iTile], PARENT[iTile]);				//	(3) RELABEL
		// ??												//	(?)	COMPACT ==> no gaps in
		second_scan(lab_mat[iTile],nr,nc, PARENT[iTile]);						//	(4) 2nd SCAN

	}
	/*
	// PRINTS:
	print_mat(	urban,	nr,		nc,		"Binary Raster"							);
	print_mat(	lab_mat,	nr+1,	nc+1, 	"CCL (Connected-Components Labeling)"	);
	printf("Indexes:\n"); for(i=0;i<=mc;i++) {printf("%3.0d ",i);};  printf("\n");
	print_vec(	PARENT,	mc+1,			"Labels"								);
	print_vec(	cont,	mc+1,			"Counts"								);
*/

	unsigned int *cross_parent,*first_pos_cp,dim=0,rr,cc;
	int nn,ww;
	for(iTile = 0;iTile<nTiles;iTile++) 	{	dim += mc[iTile]; 	}
	cross_parent=calloc(dim,sizeof(int));
	first_pos_cp=cross_parent;

	//kernel 2:
	for(iTile = 0;iTile<nTiles;iTile++)
	{
		rr = iTile / ntilesX; 		// CURRENT ROW (quoziente intero)
		nn=iTile - ntilesX;			// tile index of nn
		ww=((rr*ntilesX)==iTile)?-1:iTile-1;	// tile index of ww
		cross_parent=record_cross_equivalence(lab_mat,cross_parent,nr,nc,iTile, nn, ww, PARENT[iTile],mc[iTile]);
	}




	// FREE MEMORY:
	for(iTile = 0;iTile<nTiles;iTile++)
	{
		free(cont[iTile]);
		free(lab_mat[iTile]);
		free(PARENT[iTile]);
		free(urban[iTile]);
	}
	
	// RETURN:
	return 0;
}
