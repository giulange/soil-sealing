/*
	Object:		Raster-scan and label-equivalence-based algorithm.
	Authors:	Massimo Nicolazzo & Giuliano Langella
	email:		gyuliano@libero.it


-----------
DESCRIPTION:
-----------

 I: "urban"		--> [0,0] shifted
 O: "lab_mat"	--> [1,1] shifted

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

//	INCLUDES
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <errno.h>        /* errno */
#include <string.h>       /* strerror */
#include <math.h>			// ceil

// DEFINES
//	-indexes
#define durban(c,r)		urban[		(c)	+	(r)	*ncols		] // I: scan value at current [r,c] 
#define nw_pol(c,r)		lab_mat[	(c-1)	+	(r-1)	*(ncols)	] // O: scan value at North-West
#define nn_pol(c,r)		lab_mat[	(c+0)	+	(r-1)	*(ncols)	] // O: scan value at North
#define ne_pol(c,r)		lab_mat[	(c+1)	+	(r-1)	*(ncols)	] // O: scan value at North-East
#define ww_pol(c,r)		lab_mat[	(c-1)	+	(r+0)	*(ncols)	] // O: scan value at West
#define cc_pol(c,r)		lab_mat[	(c+0)	+	(r+0)	*(ncols)	] // O: scan value at current [r,c] which is shifted by [1,1] in O
//	-min/max
#define _max(val1,val2)		(val1)>(val2)?(val1):(val2)
#define _min(val1,val2)		(val1)<(val2)?(val1):(val2)
//#define tiledimX 12
//#define tiledimY 12

// GLOBAL VARIABLES
unsigned char 	Vb			= 0;	// background value
unsigned char 	Vo			= 1;	// object value
unsigned int	cross_cols	= 4;	// number of columns of the cross_parent matrix
unsigned char	buffer[255];

//---------------------------- FUNCTIONS PROTOTYPES
// 	FIRST STAGE
unsigned int first_scan( unsigned char *urban, unsigned int nrows, unsigned int ncols,unsigned int *lab_mat,unsigned int *count,unsigned int nID,unsigned int *PARENT);
unsigned int *record_equivalence(unsigned int val1, unsigned int val2, unsigned int nID, unsigned int *PARENT);
void union_equivalence(unsigned int nID, unsigned int *PARENT);
unsigned int *relabel_equivalence(unsigned int nID, unsigned int maxcount, unsigned int *PARENT);
void second_scan( unsigned int *lab_mat,unsigned int nrows,unsigned int ncols,unsigned int *PARENT );
void print_mat(unsigned char *u,unsigned int nrows,unsigned int ncols, char *Label);
void print_vec( unsigned int *vec, unsigned int numel, unsigned char *Label );
void read_mat(unsigned char *urban, unsigned int nrows, unsigned int ncols, char *filename);
void write_mat(unsigned int **lab_mat, unsigned int nr, unsigned int nc, unsigned int ntilesX, unsigned int ntilesY, char *filename);

// 	SECOND STAGE
unsigned int *objects_stitching_nn(unsigned int *lm_nn,unsigned int *lm_cc,unsigned int nr,unsigned int nc,unsigned int ntile_nn,unsigned int ntile_cc,unsigned int *cross_parent);
unsigned int *objects_stitching_ww(unsigned int *lm_ww,unsigned int *lm_cc,unsigned int nr,unsigned int nc,unsigned int ntile_ww,unsigned int ntile_cc,unsigned int *cross_parent);
unsigned int *objects_stitching_cc(unsigned int *lm_cc,unsigned int nr,unsigned int nc,unsigned int ntile_cc,unsigned int *cross_parent,unsigned int *cross_parent_ii,unsigned int first_empty,unsigned int *PARENT,unsigned int maxcount);
unsigned int *record_cross_equivalence(unsigned int **lm,unsigned int *cross_parent,unsigned int nr,unsigned int nc,unsigned int ntile_cc,int ntile_nn,int ntile_ww,unsigned int *PARENT,unsigned int mc);
unsigned int  union_cross_equivalence(unsigned int first_empty, unsigned int *cross_parent);
unsigned int *relabel_cross_equivalence(unsigned int *final_parent, unsigned int *cross_parent,unsigned int dim_cross,unsigned int nTiles,unsigned int *mc);
// 	THIRD STAGE
unsigned int *third_scan(unsigned int	nrows, unsigned int	ncols, unsigned int	*lab_mat, unsigned int	*cur_final_parent);
//---------------------------- FUNCTIONS PROTOTYPES


/*
	how many rows: (ntiles * 2 + ntilesX + ntilesY) *  ((tiledimX+tiledimY)/2)
*/

unsigned int first_scan(
				unsigned char	*urban,
				unsigned int	nrows,
				unsigned int	ncols,
				unsigned int	*lab_mat,
				unsigned int	*count,
				unsigned int 	nID,
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
					if (c<ncols-1 && r>0 && ne_pol(c,r)!=Vb) PARENT = record_equivalence( ne_pol(c,r), ww_pol(c,r), nID, PARENT );
				}
				
				// NORTH-WEST:
				else if(r>0 && c>0 && nw_pol(c,r)!=Vb)
				{
					cc_pol(c,r) = nw_pol(c,r);
					if ( c<ncols-1 && ne_pol(c,r)!=Vb ) PARENT = record_equivalence( ne_pol(c,r), nw_pol(c,r), nID, PARENT );
				}
				
				// NORTH-EAST:
				else if(r>0 && c<ncols-1 && ne_pol(c,r)!=Vb) cc_pol(c,r) = ne_pol(c,r);
				
				//none object pixels in mask:
				else
					cc_pol(c,r) = ++maxcount;
				
				// I am not sure that we should count right NOW.
				// Maybe It could be better to performe this is a separate phase when all 
				// equivalent labels are already known and applied.
				count[cc_pol(c,r)]++;
			}
		}
	}	
	return maxcount;
}

unsigned int * record_equivalence(unsigned int val1, unsigned int val2, unsigned int nID, unsigned int *PARENT)
{
	int skip = 0;
	int i,k;
	if(val1!=val2) // write equivalence only if values are different
	{	// check if equivalence was already written
		for (i=nID-1;i>=0;i--)
		{
			//printf("PARENT[%d,1:2]=(%d,%d)\t\tequivalence=(%d,%d)\n",i,PARENT[i*2+0],PARENT[i*2+1],_max(val1,val2),_min(val1,val2));
			if( (PARENT[i*2+0] == (int)fmax((float)val1,(float)val2)) && (PARENT[i*2+1] == (int)fmin((float)val1,(float)val2)) )
			{
				skip = 1;
			}
			if(PARENT[i*2+0]==0) k=i;
		}
		if( skip==0 ) // write new equivalence only if it is not duplicated
		{
			/*
			*PARENT++ = (int)fmax((float)val1,(float)val2);
			*PARENT++ = (int)fmin((float)val1,(float)val2);
			*/
			PARENT[k*2+0] = (int)fmax((float)val1,(float)val2);
			PARENT[k*2+1] = (int)fmin((float)val1,(float)val2);
			//printf("PARENT[%d,1:2]=(%d,%d)\n",k,PARENT[k*2+0],PARENT[k*2+1]);
		}
	}
	return PARENT;
}

void union_equivalence(unsigned int nID, unsigned int *PARENT)
{
	unsigned int j,k,*list_convergent_equivalence;
	int found_combined_equivalence=1;

	//printf("\n\nUNION_EQUIVALENCE:\n");

	// (A) compare left-hand side on "j" with left-hand side on "k"  ==> SPLIT COMBINED EQUIVALENCES
	while(found_combined_equivalence==1)
	{
		found_combined_equivalence = 0;
		for(j=0;j<nID;j++)
		{
			if(PARENT[j*2+0]==0) break; // It stops when the first zero on the left-hand side on "j" is found.

			list_convergent_equivalence = (unsigned int*)calloc(nID*2,sizeof(unsigned int));
			int cRow = 0;
			int root_eq = nID;
			for(k=j;k<nID;k++) // create a list of ID's with same left-hand side ID.
			{
				if(PARENT[k*2+0]==0) break;
				if( (PARENT[j*2+0]==PARENT[k*2+0]) ) // If a same ID is found on the left column on "j" and "k"...
				{
					if( (j!=k) && (PARENT[j*2+1]!=PARENT[k*2+1]) ) found_combined_equivalence = 1; // It is true only if I'm on different rows of PARENT and right-hand side is different (otherwise endless while loop)
					list_convergent_equivalence[2*cRow +0] = k;				// ...record "k" position and
					list_convergent_equivalence[2*cRow +1] = PARENT[k*2+1]; // ...record right-hand side ID on "k"
					root_eq = fmin((float)root_eq,(float)PARENT[k*2+1]);
					cRow+=1;
				}
			}
			for(k=0;k<cRow;k++) // split equivalences
			{
				if(list_convergent_equivalence[2*k +1]!=root_eq) // decode equivalence only if it does not contain the "root_eq" of the current list
				{
					//printf("PARENT[%d,1:2]=(%d,%d)\t--->\t",k,PARENT[(list_convergent_equivalence[2*k +0])*2+0],PARENT[(list_convergent_equivalence[2*k +0])*2+1]);
					PARENT[(list_convergent_equivalence[2*k +0])*2+0] = list_convergent_equivalence[2*k +1];
					PARENT[(list_convergent_equivalence[2*k +0])*2+1] = root_eq;
					//printf("PARENT[%d,1:2]=(%d,%d)\n",k,PARENT[(list_convergent_equivalence[2*k +0])*2+0],PARENT[(list_convergent_equivalence[2*k +0])*2+1]);
				}
			}
			if(found_combined_equivalence==1) break;
		}
	}
	// ----(A)----
}

unsigned int *relabel_equivalence(unsigned int nID, unsigned int maxcount, unsigned int *PARENT)
{
	/*
	 *
	 */
	int cc;
	int j,k,dependent_id;
	unsigned int *new_PARENT;
	new_PARENT = (unsigned int*)calloc(nID*2,sizeof(unsigned int));

	//printf("\n\nRE-LABEL EQUIVALENCE:\n");

	// (B) compare left-hand side on "j" with right-hand side on "k" ==> ASSIGN ROOT EQUIVALENCE
	cc=1;
	while (cc)
	{	cc=0;
		for(j=0;j<nID;j++)
		{
			if(PARENT[j*2+0]==0) break; // It stops when the first zero on the left-hand side on "j" is found.
			for(k=0;k<nID;k++)
			{
				if(PARENT[j*2+1]==PARENT[k*2+0])
				{
					PARENT[j*2+1]=PARENT[k*2+1];
					cc=1;
				}
			}
	// ----(B)----
		}
	}
	// (C) delete duplicates, sort IDs, define also ROOT equivalences!
	for(j=0;j<=maxcount;j++)
	{
		dependent_id = 0;
		for(k=0;k<nID;k++)
		{
			if( PARENT[k*2+0]==j )
			{
				dependent_id = 1;
				break;
			}
		}
		if( dependent_id==1 )
		{
			new_PARENT[j*2+0] = j;
			new_PARENT[j*2+1] = PARENT[k*2+1];
		}
		else
		{
			new_PARENT[j*2+0] = j;
			new_PARENT[j*2+1] = j;
		}
		//printf("PARENT[%d,1:2]=(%d,%d)\n",j,new_PARENT[j*2+0],new_PARENT[j*2+1]);
	}
	// ----(C)----
	return new_PARENT;
}

void second_scan(
		unsigned int	*lab_mat,
		unsigned int	nrows,
		unsigned int	ncols,
		unsigned int	*PARENT		)
{
	unsigned int i,j;
    for(i=0;i<nrows;i++) for(j=0;j<ncols;j++) cc_pol(j,i) = PARENT[2*cc_pol(j,i)+1];
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
void print_mat_int(unsigned int *u,unsigned int nrows,unsigned int ncols, char *Label)
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
	if (fid == NULL) { printf("Error opening file!\n"); exit(1); }
	for(rr=1;rr<nrows-1;rr++) for(cc=1;cc<ncols-1;cc++) { fscanf(fid, "%d",&a);	durban(cc,rr)=(unsigned char)a; }
	fclose(fid);
}

void write_mat(unsigned int **lab_mat, unsigned int nr, unsigned int nc, unsigned int ntilesX, unsigned int ntilesY, char *filename)
{
	unsigned int rr,cc,ntX,ntY;
	FILE *fid ;
	fid = fopen(filename,"w");
	if (fid == NULL) { printf("Error opening file %s!\n",filename); exit(1); }

	for(ntY=0;ntY<ntilesY;ntY++)
	{
		for(rr=1;rr<nr;rr++)
		{
			for(ntX=0;ntX<ntilesX;ntX++)
			{
				for(cc=1;cc<nc;cc++)
				{
					if( !(((cc==nc-1) && ((ntilesX*ntY+ntX+1)%ntilesX)==0))	&&	// do not print last column
						!(((rr==nr-1) && (ntY==ntilesY-1))) 					// do not print last row
					)
					{
						fprintf(fid, "%d ",lab_mat[ntilesX*ntY+ntX][nc*rr+cc]);
						//printf("%d ",lab_mat[ntilesX*ntY+ntX][nc*rr+cc]);
					}
				}
			}
			fprintf(fid,"\n");
			//printf("\n");
		}
	}
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
		if (lm_nn[nc*(nr-1)+c]!=0)
		{
			*cross_parent++	= ntile_cc;				// parent tile number
			*cross_parent++	= lm_cc[c];				// parent ID of parent tile number
			*cross_parent++	= ntile_nn;				// root tile number
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
		unsigned int *lm_cc,			// label matrix of target (=centre) tile in the mask
		unsigned int nr,				// number of rows
		unsigned int nc,				// number of columns
		unsigned int ntile_cc,			// number of cc tile in the grid of tiles
		unsigned int *cross_parent,		// pointer to the first row in cross_parent
		unsigned int *cross_parent_ii,	// pointer to the first element before stitching by {nn,ww}
		unsigned int first_empty,		// the first empty row in cross_parent [SCALAR]
		unsigned int *PARENT,			// the PARENT vector of target tile within tiles-mask
		unsigned int maxcount			// number of labels in target tile stored in PARENT
						)
{
	unsigned int i,j;
	unsigned char found=0;
	for(i=1;i<=maxcount;i++) // start from 1, because 0 is background
	{
		// The current ntile_cc is on the left-hand side? If yes then the equivalence with (tile) + (intra-tile ID) is not a ROOT equivalence, because cross_parent only stores non-ROOT equivalences!
		found = 0;
		for(j=0;j<first_empty;j++) if( (cross_parent_ii[j*cross_cols+0]==ntile_cc) && (cross_parent_ii[j*cross_cols+1]==PARENT[i*2+1]) ) {found = 1; break;}

		// If the current tile was not found in cross_parent then it gives a new ROOT equivalence.
		if (found==0) // then record current label as ROOT label! (i.e. elements {1,2} are equal to elements {3,4} of cross_parent
		{
			*cross_parent++	= ntile_cc;			// parent tile number
			*cross_parent++	= PARENT[i*2+0];	// parent ID of parent tile number
			*cross_parent++	= ntile_cc;			// root tile number
			*cross_parent++	= PARENT[i*2+1];	// root ID of root tile number
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

unsigned int union_cross_equivalence(unsigned int effective_nr, unsigned int *cross_parent)
{
	unsigned int i,j,curr_r_newrule=0;
	// I have to check whether one for loop suffices to solve the whole ROOT tree
	for(i=0;i<effective_nr;i++)
	{
		j = 0;
		// DO (j=j+1) Until (I will find the (first occurrence of) ROOT cross equivalence row recorded in cross_parent):
		while(!((j>i) || ((cross_parent[j*cross_cols+0]==cross_parent[i*cross_cols+2]) && (cross_parent[j*cross_cols+1]==cross_parent[i*cross_cols+3])))) j+=1;
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

	int aa,x=0;
	while(1)
	{	x++;
		aa=0;
		curr_r_newrule =0;
		for(i=0;i<effective_nr;i++)
		{
			for(j = i+1;j<effective_nr;j++)
			{
				if  (
					((cross_parent[j*cross_cols+0]==cross_parent[i*cross_cols+0]) && (cross_parent[j*cross_cols+1]==cross_parent[i*cross_cols+1]))
					&&
					((cross_parent[j*cross_cols+2]!=cross_parent[i*cross_cols+2]) || (cross_parent[j*cross_cols+3]!=cross_parent[i*cross_cols+3]))
					)
				{
/*
					printf("Found [j] %d %d %d %d \t<--> [i] %d %d %d %d\n",
							cross_parent[j*cross_cols+0],cross_parent[j*cross_cols+1],cross_parent[j*cross_cols+2],cross_parent[j*cross_cols+3],
							cross_parent[i*cross_cols+0],cross_parent[i*cross_cols+1],cross_parent[i*cross_cols+2],cross_parent[i*cross_cols+3]
					);
*/
					aa=1;
					if( (cross_parent[j*cross_cols+2]>cross_parent[i*cross_cols+2]) ||
						((cross_parent[j*cross_cols+2]==cross_parent[i*cross_cols+2]) && ( cross_parent[j*cross_cols+3]>cross_parent[i*cross_cols+3]))	)
					{
						cross_parent[(effective_nr+curr_r_newrule)*cross_cols+0]=cross_parent[j*cross_cols+2];
						cross_parent[(effective_nr+curr_r_newrule)*cross_cols+1]=cross_parent[j*cross_cols+3];
						cross_parent[(effective_nr+curr_r_newrule)*cross_cols+2]=cross_parent[i*cross_cols+2];
						cross_parent[(effective_nr+curr_r_newrule)*cross_cols+3]=cross_parent[i*cross_cols+3];
						cross_parent[j*cross_cols+2] = cross_parent[i*cross_cols+2];
						cross_parent[j*cross_cols+3] = cross_parent[i*cross_cols+3];
					}
					else
					{
						cross_parent[(effective_nr+curr_r_newrule)*cross_cols+0]=cross_parent[i*cross_cols+2];
						cross_parent[(effective_nr+curr_r_newrule)*cross_cols+1]=cross_parent[i*cross_cols+3];
						cross_parent[(effective_nr+curr_r_newrule)*cross_cols+2]=cross_parent[j*cross_cols+2];
						cross_parent[(effective_nr+curr_r_newrule)*cross_cols+3]=cross_parent[j*cross_cols+3];
						cross_parent[i*cross_cols+2] = cross_parent[j*cross_cols+2];
						cross_parent[i*cross_cols+3] = cross_parent[j*cross_cols+3];
					}
					curr_r_newrule += 1;
				}
			}
		}
		effective_nr += curr_r_newrule;
		if (aa==0) break;
	}
	return effective_nr;
}

unsigned int *relabel_cross_equivalence( unsigned int *final_parent, unsigned int *cross_parent,unsigned int dim_cross,unsigned int nTiles,unsigned int *mc)
{
	unsigned int i,j,k,iTile,final_count=1;
	// It reuses the 4th column writing the overall image IDs
	for(i=0;i<dim_cross;i++)
	{
		if( (cross_parent[i*cross_cols+0]==cross_parent[i*cross_cols+2]) && (cross_parent[i*cross_cols+1]==cross_parent[i*cross_cols+3]) )
		{
			cross_parent[i*cross_cols+3] = final_count++;
		}
		else
		{
			j = 0;
			while(!((j>i) || ((cross_parent[j*cross_cols+0]==cross_parent[i*cross_cols+2]) && (cross_parent[j*cross_cols+1]==cross_parent[i*cross_cols+3])))) j+=1;
			if (j<=i) cross_parent[i*cross_cols+3]=cross_parent[j*cross_cols+3];
			else puts("Error in cross_parent!\nCheck one of {union_cross_equivalence, record_cross_equivalence}");
		}
	}
	//costruzione di final_parente: semplificato la cross eliminando i doppioni considerando solo le prime due colonne.
	//printf("cross_parent -- nicolazzo abrupt code\n");
	k=0;
	for(iTile = 0;iTile<nTiles;iTile++)
	{
		for (j=1;j<=mc[iTile];j++)
		{
			for (i=0;i<dim_cross;i++)
			{
				if ((cross_parent[i*cross_cols+0]==iTile) && (cross_parent[i*cross_cols+1]==j)) break;
			}
			final_parent[k++]=cross_parent[i*cross_cols+3];
			//printf("%d %d %d\n",iTile,j,cross_parent[i*cross_cols+3]);
		}
	}
	return final_parent;
}


unsigned int *third_scan(
				unsigned int	nrows,
				unsigned int	ncols,
				unsigned int	*lab_mat,
				unsigned int	*cur_final_parent	)
{
	unsigned int r;
	unsigned int c;
	for(r=0; r<nrows; r++)
		for(c=0; c<ncols; c++)
			if(cc_pol(c,r)!=0) // if (r,c) is object pixel
			{	//printf("%d<-%d\n",cc_pol(c,r), cur_final_parent[cc_pol(c,r)-1]);
				cc_pol(c,r) = cur_final_parent[cc_pol(c,r)-1];
			}
	return lab_mat;
}



int main(int argc, char **argv)
{
	unsigned int tiledimX 	= atoi( argv[1] );
	unsigned int tiledimY 	= atoi( argv[2] );
	unsigned int NC1 		= atoi( argv[3] );//98; // passed by JAI
	unsigned int NR1 		= atoi( argv[4] );//98; // passed by JAI

	// DECLARATION:
	unsigned int nID		= tiledimX*tiledimY/4 +1;
	unsigned int ntilesX,ntilesY,nTiles,iTile;
	// X dir
	ntilesX 				= ceil( (NC1+2-1) / (tiledimX-1)  );
	unsigned int NC 		= ntilesX*(tiledimX-1) +1;
	// Y dir
	ntilesY 				= ceil( (NR1+2-1) / (tiledimY-1)  );
	unsigned int NR 		= ntilesY*(tiledimY-1) +1;
	//
	nTiles 					= ntilesX*ntilesY;

	unsigned int * (lab_mat[nTiles]);
	unsigned int * (cont[nTiles]);
	unsigned int i,j,k;
	unsigned int mc[nTiles];
	unsigned int *(PARENT[nTiles]);
	unsigned char *(urban[nTiles]);
	unsigned int rr,cc,nn,ww;
	unsigned int dim=0;
	unsigned int *dim_cum;
	unsigned int *cross_parent;
	unsigned int *final_parent;
	unsigned int *first_pos_cp;
	unsigned char *urban_gl;

	dim_cum = (unsigned int*)malloc(nTiles*sizeof(int));


	urban_gl	= (unsigned char*)calloc((NC)*(NR),sizeof(unsigned char));
	sprintf(buffer,"/home/giuliano/git/soil-sealing/data/ALL.txt");
	read_mat(urban_gl, NR, NC, buffer);

	// 1st KERNEL INVOCATION: intra-tile labelingt-
	for(iTile = 0;iTile<nTiles;iTile++)
	{
		// INITIALIZATION:
		lab_mat[iTile]	= (unsigned int*)calloc((tiledimX)*(tiledimY),sizeof(unsigned int));
		cont[iTile]		= (unsigned int*)malloc(nID*sizeof(unsigned int));
		PARENT[iTile]	= (unsigned int*)calloc(nID*2,sizeof(unsigned int)); //critico! controllare quanto e' la dim max possibile!!
		urban[iTile] 	= (unsigned char*)calloc((tiledimX)*(tiledimY),sizeof(unsigned char));

		for(i=0;i<tiledimY;i++)
			for(j=0;j<tiledimX;j++)
				urban[iTile][j+tiledimX*i]=urban_gl[+ i*(NC)+j 						// dentro la 1° tile
				                              + (tiledimX-1)*(iTile%ntilesX)		// itile in orizzontale
				                              + (iTile/ntilesX)*(tiledimY-1)*NC];	// itile in verticale
/*
		sprintf(buffer,"iTile=%d",iTile);
		print_mat(	urban[iTile],	tiledimY,		tiledimX, 	buffer	);
*/

		// KERNELs INVOCATION:
		mc[iTile] = first_scan(urban[iTile],tiledimY,tiledimX,lab_mat[iTile],cont[iTile],nID,PARENT[iTile]);	//	(1) 1st SCAN
		//sprintf(buffer,"iTile=%d -- 1scan",iTile);
		//print_mat_int(	lab_mat[iTile],	tiledimY,		tiledimX, 	buffer	);

		union_equivalence(		nID, PARENT[iTile]);										//	(2) UNION
		//sprintf(buffer,"ueq iTile=%d",iTile);
		//print_mat_int(	lab_mat[iTile],	tiledimY,		tiledimX, 	buffer	);

		PARENT[iTile] = relabel_equivalence(	nID, mc[iTile], PARENT[iTile]);										//	(3) RELABEL
		//sprintf(buffer,"PARENT[%d]",iTile);
		//print_mat_int(	PARENT[iTile],	nID,		2, 	buffer	);

		//																						//	(?)	COMPACT ==> no gaps in
		second_scan(lab_mat[iTile],tiledimY,tiledimX, PARENT[iTile]);										//	(4) 2nd SCAN
		//sprintf(buffer,"iTile=%d -- 2scan",iTile);
		//print_mat_int(	lab_mat[iTile],	tiledimY,		tiledimX, 	buffer	);
	}

	dim_cum[0] = 0;
	dim = mc[0];
	for(iTile = 1;iTile<nTiles;iTile++)
	{
		dim += mc[iTile];
		dim_cum[iTile] = dim_cum[iTile-1] + mc[iTile-1];
	}

	// 2nd KERNEL INVOCATION: inter-tile labeling
	cross_parent = (unsigned int*)calloc(((tiledimY+tiledimX)*nTiles+dim)*cross_cols,sizeof(int));
	final_parent = (unsigned int*)calloc(dim,sizeof(unsigned int));
	first_pos_cp = cross_parent;
	for(iTile = 0;iTile<nTiles;iTile++)
	{
		rr = iTile / ntilesX;					// CURRENT ROW (quoziente intero)
		nn = iTile - ntilesX;					// tile index of nn
		ww = ((rr*ntilesX)==iTile)?-1:iTile-1;	// tile index of ww
		cross_parent = record_cross_equivalence(lab_mat,cross_parent,tiledimY,tiledimX,iTile, nn, ww, PARENT[iTile],mc[iTile]);
	}

	unsigned int effective_nr;
	effective_nr = (cross_parent - first_pos_cp)/cross_cols;
	cross_parent = first_pos_cp;
	//sprintf(buffer,"cross_parent -- after record_cross_equivalence");
	//print_mat_int(	cross_parent,	effective_nr,		cross_cols, 	buffer	);

	// UNION CROSS...
	effective_nr = union_cross_equivalence( effective_nr, cross_parent);
	//sprintf(buffer,"cross_parent -- after union_cross_equivalence");
	//print_mat_int(	cross_parent,	effective_nr,		cross_cols, 	buffer	);

	// RELABEL CROSS...
	cross_parent = relabel_cross_equivalence(  final_parent,cross_parent,effective_nr,nTiles,mc);
	//sprintf(buffer,"cross_parent -- after relabel_cross_equivalence");
	//print_mat_int(	cross_parent,	effective_nr,		cross_cols, 	buffer	);

	// FINAL SCAN
	for(iTile = 0;iTile<nTiles;iTile++)
	{
		lab_mat[iTile] = third_scan(tiledimY, tiledimX, lab_mat[iTile], final_parent+dim_cum[iTile]);
		//sprintf(buffer,"iTile[%4.0d ] -- 3scan",iTile);
		//print_mat_int(	lab_mat[iTile],	tiledimY,		tiledimX, 	buffer	);
	}

	// SAVE lab_mat to file and compare with MatLab
	sprintf(buffer,"/home/giuliano/git/soil-sealing/data/Ccode.txt");
	write_mat(lab_mat, tiledimX, tiledimY, ntilesX, ntilesY, buffer);

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
