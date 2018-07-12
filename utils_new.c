/******************************************************************************/
/***  3D Unsteady Navier-Stokes solver                                      ***/
/***   Version 1.0                                                          ***/
/***   Version 2.0 - March 1998; P. Minev                                   ***/
/***                 fill_bb_info, find_min_max, sort2, print_*vector,      ***/
/***                 print_matrix added                                     ***/
/***   Version 2.1 - June 1998; P. Minev                                    ***/
/***                 compute_substeps added                                 ***/
/***   Version 2.2 - March 10; J. Myers                                     ***/
/***                 lag_uvw_in added                                       ***/
/***                 Courant number calcualtion added                       ***/
/***                 to the compute_substeps function                       ***/
/***   Version 2.3 - March 22; J.Myers                                      ***/
/***                 scaled alpha in get_uvw_in and lag_uvw_in              ***/
/***                 by radius to insure womersley profile correct          ***/
/***                 for inlets with radius other than 1                    ***/
/***                 Changed the way zt,zf10, and zu calculated within      ***/
/***                 get_uvw_in and lag_uvw_in (fixed error in              ***/
/***                 Womersley flow calculation )                           ***/
/***                 Added pass and check of structure bc_table to          ***/
/***                 get_uvw_in                                             ***/
/***                 (newtetr_v3.2)                                         ***/
/***   Version 2.4 - April 16; J.Myers                                      ***/
/***                 get_uvw_in and lag_uvw_in (fixed error in              ***/
/***                 Analytical flow calculation)                           ***/
/***                                                                        ***/
/***                                                                        ***/
/***                                                                        ***/
/***                                                                        ***/
/***   Copyright 1994 - C. Ross Ethier, D.A. Steinman, X. Zhang             ***/
/***                                                                        ***/
/***  FILE NAME:   utils.c                                                  ***/
/***                                                                        ***/
/***  DESCRIPTION:     This file contains assorted utility functions        ***/
/***                     add_vectors                                        ***/
/***                     add_vectors_nonD                                   ***/
/***                     bin_search                                         ***/
/***              Cabs, Cadd, Cdiv, Cexp, Cmul, Complex, Csqrt              ***/
/***                     compute_substeps                                   ***/
/***                     copy_vector                                        ***/
/***                     copy_vector_D                                      ***/
/***                     dot_product                                        ***/
/***                     dump_inlet_vel                                     ***/ 
/***                     fill_bb_info                                       ***/
/***                     find_min_max                                       ***/
/***                     generate_tau                                       ***/
/***                     get_uvw_in                                         ***/
/***                     lag_uvw_in                                         ***/  
/***                     L2_norm                                            ***/
/***                     Message                                            ***/
/***                     pack_vector                                        ***/
/***                     print_ivector                                      ***/
/***                     print_matrix                                       ***/
/***                     print_vector                                       ***/
/***                     RCmul                                              ***/
/***                     shell                                              ***/
/***                     sort2                                              ***/
/***                     unpack_vector                                      ***/
/***                     usage                                              ***/
/***                     zbes                                               ***/
/***                     zero_vector                                        ***/
/***                     zero_vector_nonD                                   ***/
/***                                                                        ***/
/******************************************************************************/

#include "../include/newtetr.h"

/************************************************************************
FUNCTION: add_vectors
 PURPOSE: adds factor1*vector1 + factor2*vector 2, overwrites sum vector with answer
PARAMETERS:
   Input:
     - sum: vector to have sum placed in it.
     - vector1, vector2: input vectors
     - factor1, factor2: factors multiplying input vectors
     - size: vector length
   Output:
     - sum: overwritten with answer
 RETURNS: void
************************************************************************/

void add_vectors(VECTOR sum, VECTOR vector1, VECTOR vector2, FLOAT factor1,
   FLOAT factor2, NODE size)
{
   FLOAT
     *last_row_ptr;

   last_row_ptr = sum + size;
   for (; sum <= last_row_ptr; *sum++ = factor1* *vector1++ + factor2 * *vector2++);
   return;
}

/************************************************************************
FUNCTION: add_vectors_nonD
 PURPOSE: adds factor1*vector1 + factor2*vector 2, overwrites sum vector with answer
          IMPORTANT: Loops over only non-Dirichlet entries.
PARAMETERS:
   Input:
     - sum: vector to have sum placed in it.
     - vector1, vector2: input vectors
     - factor1, factor2: factors multiplying input vectors
     - size: # of entries in vector to be updated
   Output:
     - sum: overwritten with answer
 RETURNS: void
************************************************************************/

void add_vectors_nonD(VECTOR sum, VECTOR vector1, VECTOR vector2, FLOAT factor1,
   FLOAT factor2, CONDENSE_INFO *Condense)
{
   FLOAT
      *last_row_ptr;

   char
      *non_D_ptr = Condense->non_D;

   last_row_ptr = sum + Condense->rowmap[Condense->tot_rows];
   while(sum <= last_row_ptr)
   {
      if (*non_D_ptr++)
         *sum = factor1 * *vector1 + factor2 * *vector2;
      sum++;
      vector1++;
      vector2++;
   }
   return;
}

/************************************************************************
FUNCTIONS: bin_search
PURPOSE: performs a binary search in an ordered array of integers 
PARAMETERS:
   Input: array_ptr: the pointer to the 0-th element of the array
          element: the integer to be searched for
 RETURNS: the position of element if it is found in the array;
          the position of the right closest entry of the array if element
          is not found in the array(with - sign); if it's to the right of all entries
          of the array it returns the size of the array + 1.
************************************************************************/
long int bin_search( NODE *array_ptr, NODE element, NODE l_count, NODE r_count)
{
NODE left, right, m_count1, m_count2, middle1, middle2, shift;

     left = array_ptr[l_count]; right = array_ptr[r_count];

     if( element > right ) return -((long int)r_count+1);
       else if( element < left ) return -(long int)l_count;
     
            shift = (NODE)(((float)(r_count-l_count))/2e0+1e-1);
            m_count1=l_count+shift;
                    middle1=array_ptr[m_count1];
            m_count2=r_count-shift;
                    middle2=array_ptr[m_count2];

                if( element == left ) return (long int)l_count;
                  else if( element == right ) return (long int)r_count;
                       else if( element < middle2 && element > middle1 ) return -(long int)m_count2;

                       else
                       {
                       while( shift!=0 )
                       {
                       if( element < middle1 )
                         return bin_search( array_ptr, element, l_count, m_count1);
                     
                        else if( middle2 < element ) 
                                return bin_search( array_ptr, element, m_count2, r_count);
       
                                else if( middle1 == element ) return (long int)m_count1;
                                        else if( middle2 == element ) return (long int)m_count2;
                        } /* END while */
                       
                        return -(long int)r_count;
                        } /* END else */
}

/************************************************************************
FUNCTIONS: Cabs, Cadd, Cdiv, Cexp, Cmul, Complex, Csqrt (from Numerical Recipes)
 PURPOSE: handle operations with complex numbers, namely
         modulus, addition, division, exponentiation, multiplication,
         formation of a complex # from 2 reals, square root, subtraction (respectively)
PARAMETERS:
   Input: various; see functions for self-explanatory description
 RETURNS: ditto
************************************************************************/

FLOAT Cabs(COMPLEX z)
{
   FLOAT x,y,ans,temp;

   x=fabs(z.re);
   y=fabs(z.im);
   if (x == 0.0)
      ans=y;
   else if (y == 0.0)
      ans=x;
   else if (x > y) {
      temp=y/x;
      ans=x*sqrt(1.0+temp*temp);
   } else {
      temp=x/y;
      ans=y*sqrt(1.0+temp*temp);
   }
   return (FLOAT) ans;
}

/******************************************************************************/
COMPLEX Cadd(COMPLEX a, COMPLEX b)
{
   COMPLEX
     c;
   c.re=a.re+b.re;
   c.im=a.im+b.im;
   return c;
}

/******************************************************************************/
COMPLEX Cdiv(COMPLEX a, COMPLEX b)
{
   COMPLEX
     c;
   FLOAT
     r,
     den;

   if (fabs(b.re) >= fabs(b.im)) {
      r=b.im/b.re;
      den=b.re+r*b.im;
      c.re=(a.re+r*a.im)/den;
      c.im=(a.im-r*a.re)/den;
   } else {
      r=b.re/b.im;
      den=b.im+r*b.re;
      c.re=(a.re*r+a.im)/den;
      c.im=(a.im*r-a.re)/den;
   }
   return c;
}

/******************************************************************************/
COMPLEX Cexp(COMPLEX a)
{
   FLOAT
     im = exp(a.re),
     re = exp(a.re);

   re *= cos(a.im);
   im *= sin(a.im);

   return Complex(re, im);
}

/******************************************************************************/
COMPLEX Cmul(COMPLEX a, COMPLEX b)
{
   COMPLEX
     c;
   c.re=a.re*b.re-a.im*b.im;
   c.im=a.im*b.re+a.re*b.im;
   return c;
}

/******************************************************************************/
COMPLEX Complex(FLOAT re, FLOAT im)
{
   COMPLEX c;
   c.re=re;
   c.im=im;
   return c;
}

/******************************************************************************/
COMPLEX Csqrt(COMPLEX z)
{
   COMPLEX
     c;
   FLOAT
     x,
     y,
     w,
     r;

   if ((z.re == 0.0) && (z.im == 0.0)) {
      c.re=0.0;
      c.im=0.0;
      return c;
   } else {
      x=fabs(z.re);
      y=fabs(z.im);
      if (x >= y) {
         r=y/x;
         w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
      } else {
         r=x/y;
         w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
      }
      if (z.re >= 0.0) {
         c.re=w;
         c.im=z.im/(2.0*w);
      } else {
         c.im=(z.im >= 0) ? w : -w;
         c.re=z.im/(2.0*c.im);
      }
      return c;
   }
}

/************************************************************************
FUNCTION: compute_substeps
 PURPOSE: computes the number of substeps for the current step according to the 
          critical Courant number
PARAMETERS:
   Input:
     - tot_nodes: total number of nodes
     - Stokes:  structure with the data for the timestepping
     - courant: contains the minimum distance around each node
     - uvw_old1: the velocity components
     - Re: Reynolds number
     - alpha: alpha parameter
     - crcrit: critical Courant number
   Output:
     - Stokes->tot_substeps appropriately filled
 RETURNS: void
************************************************************************/
void compute_substeps( NODE tot_nodes, short *substep_ptr, FLOAT dt, 
                       FLOAT *courant, VECTOR u_old1, VECTOR v_old1, 
                       VECTOR w_old1, FLOAT Re, FLOAT alpha, FLOAT crcrit)
{
NODE i,   /*  loop variables */
  cu_n=0; /* small Courant number counter */
FLOAT cr_max=0.,   /* maximum Courant number */
      cr_min =1000., /* minimum Courant number */
      cr_temp;

/*    compute the maximum Courant number    */
      for( i=1; i <= tot_nodes; i++ )
      {
      cr_temp = Re*sqrt(u_old1[i]*u_old1[i]+v_old1[i]*v_old1[i]+w_old1[i]*w_old1[i])
                *dt/(alpha*alpha*courant[i]);
      
      if( cr_temp > cr_max ) 
      cr_max=cr_temp;
      if(cr_temp < cr_min)
      cr_min=cr_temp;
      if (cr_temp < DEFAULT_CRMIN)
      cu_n++;
      }

/*       compute the number of substeps   */
         cr_temp = (short) (cr_max/crcrit) + 1;
         if( cr_temp > *substep_ptr )
         *substep_ptr = cr_temp;
         Message(1,"  Number of substeps for this step: %d \n", *substep_ptr);
         Message(1,"  Maximum Courant number is: %lf \n", cr_max); 
         Message(1,"  Minimum Courant number is: %lf \n", cr_min);
         Message(1,"  %d locations found where local Courant numbers < %lf (compute_substeps)\n", cu_n, DEFAULT_CRMIN); 
 
}

/************************************************************************
FUNCTION: copy_vector
 PURPOSE: copies vector source to target
PARAMETERS:
   Input:
     - source vector
     - target vector
     - vector length
   Output:
     - target vector, appropriately filled
 RETURNS: void
************************************************************************/
void copy_vector(VECTOR target, VECTOR source, NODE length)
{

   NODE
     iloop;

   for (iloop=1; iloop <= length; target[iloop] = source[iloop], iloop++);
}

/************************************************************************
FUNCTION: copy_vector_D
 PURPOSE: copies vector source to target, looping only over Dirichlet nodes
PARAMETERS:
   Input:
     - target vector
     - source vector
     - Condese: info about which nodes are Dirichlet
   Output:
     - target vector, appropriately filled
 RETURNS: void
************************************************************************/
void copy_vector_D(VECTOR target, VECTOR source, CONDENSE_INFO *Condense)
{

   NODE
     iloop,
     length = Condense->tot_rows;

   char
      *non_D_ptr = Condense->non_D + 1;  

   for (iloop = 1; iloop <= length; iloop++)
   {
      if (!(*non_D_ptr++))
         *target = *source;
      target++;
      source++;
   }
}


/************************************************************************
FUNCTION: dot_product
 PURPOSE: computes dot product of Vector1 and Vector 2, looping over nonD entries only
PARAMETERS:
   Input:
     - Vector1, Vector2: 2 vectors to use
     - Condense: info about non-D entries in vectors
 RETURNS: value of dot product
************************************************************************/
FLOAT dot_product(VECTOR Vector1, VECTOR Vector2, CONDENSE_INFO *Condense)
{
   FLOAT
      product;               /* value of dot product */

   NODE
      *last_row_ptr,         /* pointer to last non-Dirichlet row number in vectors */
      row,                   /* current row in vectors */
      *row_ptr;              /* points to row */

   last_row_ptr = Condense->rowmap + Condense->tot_rows;
   product = 0.0;

   /* Loop over all non-Dirichlet rows */
   for (row_ptr = Condense->rowmap + 1; row_ptr <= last_row_ptr; row_ptr++)
   {
      row = *row_ptr;
      product += *(Vector1 + row) * *(Vector2 + row);
   }
   return product;
}


/**************************************************
FUNCTION: dump_inlet_vel

****************************************************/
void dump_inlet_vel(INLET_INFO *Inlets, FLOAT **nodal_xyz, NODE tot_nodes,
                    FLOAT *u,FLOAT *v,FLOAT *w,FLOAT current_time, 
		    short timestep, short timestep_start )
{
   NODE  i,node ;
   static NODE center_node;
   FLOAT dist=100.0, tmp,vel_mag;
   FLOAT xc,yc,zc, x,y,z;
   /*** assuming that there is only one inlet **/
   if( timestep == timestep_start + 1 )  {
      xc = Inlets->origin[1][1];
      yc = Inlets->origin[1][2];
      zc = Inlets->origin[1][3];
      
      for(i=1;i<=Inlets->nodes_per_inlet[1];i++)  {
         node = Inlets->nodes[1][i] ;
         x = nodal_xyz[node][1];
	 y = nodal_xyz[node][2];
	 z = nodal_xyz[node][3];
	 tmp = sqrt( pow(x-xc,2.0) + pow(y-yc,2.0 ) + pow(z-zc,2.0) ) ;
         
	 if( tmp < dist )  {
	    dist = tmp ;
	    center_node = node ;
	    }
	       
         }
      }
   
   vel_mag = sqrt( pow(u[center_node],2.0) + pow(v[center_node],2.0) 
             + pow(w[center_node],2.0) );

   printf("Time: %f Velocity at an inlet node: u= %f v= %f w= %f , Mag= %f \n",
           current_time, u[center_node], v[center_node],w[center_node],vel_mag); 
   return ;
}

/************************************************************************
FUNCTION: fill_bb_info
PURPOSE: fills the information about a grid of bounding boxes
PARAMETERS:
   Input:
     - nodal_xyz: the array of the coordinates of the nodal points
     - tot_nodes: the number of nodes in the mesh
     - bb: the structure for the bounding box grid which is partially filled
   Output:
     - the remaining entries of the structure bb are appropriately
       filled
 RETURNS: void
************************************************************************/
void fill_bb_info( FLOAT **nodal_xyz, NODE tot_nodes, BBOX_INFO *bb)
{
NODE count=0, i, iloop, j, k, l; /* loop variables and counters */
NODE 
   *pos, /* store the coordinates of the point in the grid of bounding boxes */
   *sb, *sa,  /* temporary rearrangement of the node numbers in order to speed up filling the struct  */
   store1, store2; /* temporary storage variables */

NODE slice, row, column; /* store the slice, row and column number of the current point */
NODE offset1, offset2;  /* store the offsets in column_info */

long int pos_row, pos_column, pos_point;

short unsigned calloc_err = FALSE;

          /* decide which is slice direction, which is row direction & which is column direction */
          bb->n_slices = bb->size_x; bb->d_slices=1; bb->d_rows=2; bb->d_columns=3;
               if ( bb->size_y > bb->n_slices ) 
                  {bb->n_slices = bb->size_y; bb->d_slices=2;
                   bb->d_rows=1; bb->d_columns=3;}
                  if ( bb->size_z > bb->n_slices ) 
                     {bb->n_slices = bb->size_z; bb->d_slices=3;
                      bb->d_rows=1; bb->d_columns=2;}
          bb->n_bb = 0;

/***  Sort the points of the mesh with increasing slice-coordinate; the convention
      is that the three principle directions are called slice, row and column
      directions;  the slice direction is chosen to be the one containing the 
      largest number of voxels (bounding boxes); it is previously given in
      bb->d_slices.  sa has the "voxel co-ordinates" of each node in the mesh.
      sb has the corresponding node #.  After the shell sort is called, sa is
      ordered in increasing order.  sb is permuted in the same way as sa.
      Why do we do this?  Having the nodes sorted  in increase slice co-ordinate order means
      that we can step through all nodes in an orderly way.  This avoids the use of a linked
      list data structure, and means that the arrays in the bb structure can be "grown"
      always from the end, without having to insert elements into the middle of an existing
      array.  ***/

      if((sa = (NODE *) CALLOC(tot_nodes+1, sizeof(NODE))) == NULL)
         Message(FATAL," in calloc for sa.");
      if((sb = (NODE *) CALLOC(tot_nodes+1, sizeof(NODE))) == NULL)
         Message(FATAL," in calloc for sb.");

      for( i=1; i<=tot_nodes; i++) {
         if( bb->d_slices==1 )
            sa[i]=(NODE)((nodal_xyz[i][1]-bb->min_x)/bb->step+1);
         else if( bb->d_slices==2 )
            sa[i]=(NODE) ((nodal_xyz[i][2]-bb->min_y)/bb->step+1);
         else if( bb->d_slices==3 )
            sa[i]=(NODE) ((nodal_xyz[i][3]-bb->min_z)/bb->step+1);
         sb[i] = i;
      }
      shell( tot_nodes, sa, sb); /* Call sorting routine */

/*** Check the sorting  ***/
      for( iloop=1; iloop<tot_nodes; iloop++ ) {
         if( sa[iloop] > sa[iloop+1] )
            Message(FATAL," In fill_bb _info; wrong sorting; points %d %d.", sb[iloop], sb[iloop+1]);
      }

      FREE( sa );

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/

   if((pos = (NODE *) CALLOC(4, sizeof(NODE))) == NULL)
      Message(FATAL," in calloc for pos.");

/***  Alocate memory for the arrays of the structure bb ***/
/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
   if((bb->slice_info1 = (NODE **) CALLOC(bb->n_slices + 1, sizeof(NODE *))) == NULL)
      Message(FATAL," in calloc for bb->slice_info1.");
   for (i=1; i <= bb->n_slices; i++)
       if ((bb->slice_info1[i] = (NODE*) CALLOC(2, sizeof(NODE ))) == NULL)
          calloc_err = TRUE;
   if (calloc_err)
      Message(FATAL," In fill_bb_info; calloc for rows of slice_info1\n");

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
      if((bb->slice_info2 = (NODE **) CALLOC(bb->n_slices + 1, sizeof(NODE *))) == NULL)
         Message(FATAL," In fill_bb_info; in calloc for bb->slice_info2.");
      for (i=1; i <= bb->n_slices; i++)
      {
      if ((bb->slice_info2[i] = (NODE*) CALLOC(2, sizeof(NODE ))) == NULL)
      calloc_err = TRUE;

      bb->slice_info2[i][0] = 0;
      }     
      if (calloc_err)
         Message(FATAL,"In fill_bb_info; in calloc for rows of slice_info2\n");

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
   if((bb->points_in_bb = (NODE **) CALLOC(2, sizeof(NODE *))) == NULL)            
      Message(FATAL," In fill_bb_info; calloc for bb->points_in_bb.");
       if ((bb->points_in_bb[1] = (NODE *) CALLOC(1, sizeof(NODE ))) == NULL)
          calloc_err = TRUE;

   if (calloc_err)
       Message(FATAL," In fill_bb_info; in calloc for rows of slice_info1\n");

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
      if((bb->column_info = (NODE *) CALLOC(1, sizeof(NODE))) == NULL)            
         Message(FATAL," In fill_bb_info; calloc for bb->column_info.");

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
/*   Initialize the offsets for the slices  */
     for( i=1; i<= bb->n_slices; i++)
        bb->slice_info1[i][0]=1;

/*** Now do the hard job - filling the arrays of the structure ***/

  for ( iloop=1; iloop<=tot_nodes; iloop++ ) {
     i = sb[iloop];

     pos[1]=(NODE) ((nodal_xyz[i][1]-bb->min_x)/bb->step+1);
     pos[2]=(NODE) ((nodal_xyz[i][2]-bb->min_y)/bb->step+1);
     pos[3]=(NODE) ((nodal_xyz[i][3]-bb->min_z)/bb->step+1);

     slice = pos[bb->d_slices];
     row =  pos[bb->d_rows];
     column = pos[bb->d_columns];
     if( slice > bb->n_slices )
        Message(FATAL," In fill_bb_info;  slice # larger than the total number of slices.\n");


/**  Find the position of row in slice_info2 **/
    if( bb->slice_info2[slice][0] < 1 ) /* No entries yet in this row */
       pos_row = -1; 
    else
       pos_row = bin_search( bb->slice_info2[slice], row, 1, bb->slice_info2[slice][0] );

    if( pos_row < 0 ) {
/*  Point not found; we must realloc all the arrays of the structure;
    it is supposed that the points are sorted
    into ascending slice coordinate thus only the current row offsets must be updated.
    See documentation for bin_search to understand possible returned values  */

     pos_row = fabs(pos_row);
     bb->slice_info2[slice][0] += 1;
     bb->n_bb += 1;

     if (( bb->slice_info2[slice] = (NODE *) REALLOC (bb->slice_info2[slice], 
         ( bb->slice_info2[slice][0]+1)*sizeof(NODE))) == NULL)
        Message(FATAL," In fill_bb_info; realloc for rows of bb->slice_info2.\n");

     if (( bb->slice_info1[slice] = (NODE *) REALLOC (bb->slice_info1[slice], 
         ( bb->slice_info2[slice][0]+1)*sizeof(NODE))) == NULL)
        Message(FATAL," In fill_bb_info; realloc for rows of bb->slice_info1.\n");

     if (( bb->column_info = (NODE *) REALLOC (bb->column_info, 
         (bb->n_bb+1)*sizeof(NODE))) == NULL)
        Message(FATAL," In fill_bb_info; realloc for bb->column_info.\n");

     if((bb->points_in_bb = (NODE **) REALLOC(bb->points_in_bb,
        (bb->n_bb+1)*sizeof(NODE *))) == NULL)
        Message(FATAL," In fill_bb_info; realloc for bb->points_in_bb.");
     if((bb->points_in_bb[bb->n_bb] = (NODE *) CALLOC(1, sizeof(NODE))) == NULL)
        Message(FATAL," In fill_bb_info;calloc for the last row of bb->points_in_bb.");
            
/*  Shift the entries in points_in_bb and column_info so that we can fit the new entry into the
    middle somewhere */

    if( pos_row == bb->slice_info2[slice][0] ) {
        /*  The new row is the rightmost row in the current slice; we do not shift it; update the offset for it */
        offset1 = bb->n_bb;
        bb->slice_info1[slice][pos_row] = bb->n_bb - bb->slice_info1[slice][0];
    }
    else
       offset1 = bb->slice_info1[slice][0]+bb->slice_info1[slice][pos_row];

  for( j=bb->n_bb; j>=offset1+1; j-- ) 
     bb->column_info[j] = bb->column_info[j-1];

  for( j=bb->n_bb; j>=offset1+1; j-- ) {
     if((bb->points_in_bb[j] = (NODE *) REALLOC(bb->points_in_bb[j],
        (bb->points_in_bb[j-1][0]+1)*sizeof(NODE))) == NULL)
        Message(FATAL," In fill_bb_info; realloc for columns of bb->points_in_bb.");

           for( k=0; k<=bb->points_in_bb[j-1][0]; k++ )
           bb->points_in_bb[j][k]=bb->points_in_bb[j-1][k];
  }
  /*  Put the number of the current point into corresponding entry of points_in_bb  */
      if((bb->points_in_bb[offset1] = (NODE *) REALLOC(bb->points_in_bb[offset1],
         2*sizeof(NODE))) == NULL)
         Message(FATAL," In fill_bb_info; realloc for row of bb->points_in_bb.");
       bb->points_in_bb[offset1][0] = 1;
       bb->points_in_bb[offset1][1] = i;
       bb->column_info[offset1] = column;

   /* Shift entries in the arrays info1 and info2      */
    
           for( j=bb->slice_info2[slice][0]; j>=pos_row+1; j-- )
           {
           bb->slice_info2[slice][j] = bb->slice_info2[slice][j-1];
           bb->slice_info1[slice][j] = bb->slice_info1[slice][j-1]+1;
           }
           bb->slice_info2[slice][pos_row] = row;

   /*  Update the offsets of the slices below the current one  */

           for( j=slice+1; j<=bb->n_slices; j++ )
           bb->slice_info1[j][0] = bb->n_bb + 1;
     } /* END if; case point is not found in the search through slices and rows */            
 
    else {
    /**  The slice and row are not empty; check if the column is empty  **/
          offset1 = bb->slice_info1[slice][0]+bb->slice_info1[slice][pos_row];
          if( pos_row < bb->slice_info2[slice][0] ) 
            offset2 = bb->slice_info1[slice][0]+bb->slice_info1[slice][pos_row+1]-1;
            else
                offset2 = bb->n_bb;

             if( offset1>offset2 )
               Message(FATAL," In fill_bb_info; offset1 > offset2.");

             pos_column = bin_search( bb->column_info, column, offset1, offset2 );

     if( pos_column<0 )
/*  slice and row  already exist but the column not */
     {
     pos_column = fabs(pos_column);
/**  Realloc column_info and points_in_bb  **/
             
     bb->n_bb += 1;

     if (( bb->column_info = (NODE *) REALLOC (bb->column_info,
         (bb->n_bb+1)*sizeof(NODE))) == NULL)
        Message(FATAL," In fill_bb_info; realloc for bb->column_info.\n");

     if((bb->points_in_bb = (NODE **) REALLOC(bb->points_in_bb,
        (bb->n_bb+1)*sizeof(NODE *))) == NULL)
        Message(FATAL," In fill_bb_info; realloc for bb->points_in_bb.");
     if((bb->points_in_bb[bb->n_bb] = (NODE *) CALLOC(1, sizeof(NODE))) == NULL)
        Message(FATAL," In fill_bb_info;calloc for the last row of bb->points_in_bb.");

/*  Shift the entries in points_in_bb and column_info */

    offset1 = pos_column;
    for( j=bb->n_bb; j>=offset1+1; j-- ) {
       bb->column_info[j] = bb->column_info[j-1];
    }

       for( j=bb->n_bb; j>=offset1+1; j-- ) {
          if((bb->points_in_bb[j] = (NODE *) REALLOC(bb->points_in_bb[j],
             (bb->points_in_bb[j-1][0]+1)*sizeof(NODE))) == NULL)
          Message(FATAL," In fill_bb_info; realloc for columns of bb->points_in_bb.");

          for( k=0; k<=bb->points_in_bb[j-1][0]; k++ )
             bb->points_in_bb[j][k]=bb->points_in_bb[j-1][k];
       }
/*  Put the number of the current point into corresponding entry of points_in_bb  */
      bb->points_in_bb[offset1][0] = 1;
      bb->column_info[offset1] = column;
      if((bb->points_in_bb[offset1] = (NODE *) REALLOC(bb->points_in_bb[offset1],
         2*sizeof(NODE))) == NULL)
         Message(FATAL," In fill_bb_info; realloc for row of bb->points_in_bb.");
         bb->points_in_bb[offset1][1] = i;
         bb->column_info[offset1] = column;

/*  Update the offsets of the rows to the right of the current one */
    
           for( j=bb->slice_info2[slice][0]; j>=pos_row+1; j-- )
              bb->slice_info1[slice][j] += 1;
           
/*  Update the offsets of the slices below the current one  */

           for( j=slice+1; j<=bb->n_slices; j++ )
              bb->slice_info1[j][0] = bb->n_bb+1;

     } /* END if; the bounding box was empty */
           
      else {
/**  The bounding box contains an entry already  **/
/*  Put the number of the current point into corresponding entry of points_in_bb  */
          pos_point = bb->points_in_bb[pos_column][0] += 1;
          if((bb->points_in_bb[pos_column] = (NODE *) REALLOC(bb->points_in_bb[pos_column],
             (bb->points_in_bb[pos_column][0]+1)*sizeof(NODE))) == NULL)
             Message(FATAL," In fill_bb_info; realloc for row of bb->points_in_bb.");
          bb->points_in_bb[pos_column][pos_point] = i;
      } /* END else; the bounding box was not empty */
             } /* END else; the slice  and row were not empty */
/*
Message(1,"\n Point : %d", i);
Message(1,"\n Number of bounding boxes: %d", bb->n_bb);
for( j=1; j<=bb->n_bb; j++ )
print_ivector( bb->points_in_bb[j], bb->points_in_bb[j][0], "\n bounding box %d \n",j);

for( j=1; j<=bb->n_slices; j++ )
{
Message(1,"\n Offset for slice: %d - %d \n", j, bb->slice_info1[j][0]);
print_ivector( bb->slice_info1[j], bb->slice_info2[j][0], "\n offsets for slice %d \n",j);
}
*/

   } /* END loop over the nodes of the mesh; index iloop */

   Message(1,"\n Number of bounding boxes: %d \n", bb->n_bb);
   /*for( i=1; i<=bb->n_bb; i++ )
   print_ivector( bb->points_in_bb[i], bb->points_in_bb[i][0], "\n bounding box %d \n",i);

   for( i=1; i<=bb->n_slices; i++ )
   print_ivector( bb->slice_info1[i], bb->slice_info2[i][0], "\n offsets for slice %d \n",i);  
   */

   /***   Now check if the info is correct.  We do this by computing the voxel co-ordinates of
    each mesh point, and by then searching for the point in the bb structure.  If the voxel
    co-ordinates returned by the bb structure search match the ones known beforehand, then
    everything is OK  ***/

   for ( i=1; i<=bb->n_slices; i++ ) {
       if( bb->slice_info2[i][0]==0 ) 
         count++; /* count is # empty slices, should be zero */
       for( j=1; j<=bb->slice_info2[i][0]; j++ ) {
           offset1 = bb->slice_info1[i][0] + bb->slice_info1[i][j];
           if( j==bb->slice_info2[i][0] ) {/* if true, last row in this slice */
             if( i==bb->n_slices )
               offset2 = bb->n_bb;
               else
                   offset2 = bb->slice_info1[i+1][0]-1;
        }
         else
             offset2 = bb->slice_info1[i][0] + bb->slice_info1[i][j+1]-1;
          
         for( k=offset1; k<=offset2; k++ )
             for( iloop=1; iloop<=bb->points_in_bb[k][0]; iloop++ ) {
                l = bb->points_in_bb[k][iloop]; /* Nodal point # */
                pos[1]=(nodal_xyz[l][1]-bb->min_x)/bb->step+1;
                pos[2]=(nodal_xyz[l][2]-bb->min_y)/bb->step+1;
                pos[3]=(nodal_xyz[l][3]-bb->min_z)/bb->step+1;
                slice = pos[bb->d_slices];
                row =  pos[bb->d_rows];
                column = pos[bb->d_columns];
                if( slice!=i || row!=bb->slice_info2[i][j] || column!=bb->column_info[k] ) {

                   Message(1,"\n Directions: slices %d, rows %d, columns %d; step: %lf\n",
                           bb->d_slices, bb->d_rows, bb->d_columns, 
                           bb->step);
 
                   Message(FATAL,"\n Point %d (%d, %d, %d) is in a wrong box; bb-coord: %d, %d, %d", 
                           bb->points_in_bb[k][iloop], slice, row, column,
                           i, bb->slice_info2[i][j], bb->column_info[k]);
             }
          } /* END loop; index iloop */
       } /* END loop; index j */
    }/* END loop; index i */
    Message(1,"\n %d empty slices in the bb grid.\n", count);

   FREE( sb );
   FREE( pos );
   return;
}

/************************************************************************
FUNCTION: find_min_max
PURPOSE: computes the bounding box of the domain, the minimum distance between
         the points and the maximum distance in each direction.  This is used
	 to partially fill the bb1 & bb2 structures, as well as the courant vector.
PARAMETERS:
   Input:
     - nodal_xyz: the array of the coordinates of the nodal points
     - con_table: connectivity table
     - tot_nodes: the number of nodes in the mesh
     - tot_elements: the  number of elements in the mesh
     - node_el: the invert of the connectivity table
     - node_el_fill: the fill of node_el
   Output:
     - the corresponding entries of the structures bb1 and bb2 are appropriately
       filled; array courant filled with the minimum distance around each point
 RETURNS: void
************************************************************************/
void find_min_max( FLOAT **nodal_xyz, NODE **con_table, NODE **node_el, 
                   short *node_el_fill, NODE tot_nodes, NODE tot_elements, 
                   BBOX_INFO *bb1, BBOX_INFO *bb2, FLOAT *courant)
{
NODE i, j, k, node1, el,  /* loop variables */
     node1_store, node2_store; /*  store the points where the minimum distance is achieved */
FLOAT dist, dist_min;     /* distance between two points */

/* Initialization of certain entries of the structures  */

      bb1->min_x = bb1->min_y = bb1->min_z = 10e6;
      bb1->max_x = bb1->max_y = bb1->max_z = -10e6;
      bb1->step = 0.;
      bb2->step = 10e6;

      for ( i=1; i<= tot_nodes; i++ )
      {
          if ( nodal_xyz[i][1] < bb1->min_x ) bb1->min_x = nodal_xyz[i][1];
          if ( nodal_xyz[i][2] < bb1->min_y ) bb1->min_y = nodal_xyz[i][2];
          if ( nodal_xyz[i][3] < bb1->min_z ) bb1->min_z = nodal_xyz[i][3];

       if ( nodal_xyz[i][1] > bb1->max_x ) bb1->max_x = nodal_xyz[i][1];
       if ( nodal_xyz[i][2] > bb1->max_y ) bb1->max_y = nodal_xyz[i][2];
       if ( nodal_xyz[i][3] > bb1->max_z ) bb1->max_z = nodal_xyz[i][3];
       }

/*   Now find the minimum distance within the points of the mesh.
 Note from CRE:  shold rewrite this at some point to use node-node connectivity table
 stored in MATRIX struct   */

           for ( i=1; i<=tot_nodes; i++ )
           {
           dist_min = 1e6;

               for ( j=1; j<=node_el_fill[i]; j++ )
               {
               el=node_el[i][j];
                    for ( k=1; k<= NODES_PER_EL; k++ )
                    {
                    if ( con_table[el][k] != i )
                       {
                       node1 = con_table[el][k]; 
                       dist = sqrt((nodal_xyz[node1][1]-nodal_xyz[i][1])*
                                   (nodal_xyz[node1][1]-nodal_xyz[i][1])+
                                   (nodal_xyz[node1][2]-nodal_xyz[i][2])*
                                   (nodal_xyz[node1][2]-nodal_xyz[i][2])+
                                   (nodal_xyz[node1][3]-nodal_xyz[i][3])*
                                   (nodal_xyz[node1][3]-nodal_xyz[i][3]));
                       if( dist < dist_min ) 
                       dist_min = dist; 

                    } /*END if */
                    } /*END for loop; index k */
                } /*END for loop; index j */

/*  for each point load the minimum distance to any other point into courant[i]  */

           courant[i] = dist_min;

/*  find the minimum and maximum of all local minimums; they are used as voxel sizes
    of the "large" and "small" structured grids;  in the current version only
    the "large" structured grid is used because it seems to be enough for a fast
    search but if necessary the bb2 structure can also be used    */

           if ( dist_min < bb2->step )  bb2->step = dist_min;
           if ( dist_min > bb1->step )  bb1->step = dist_min;
           } /* END for loop; index i */ 
 
        bb2->step = bb2->step; /*From CRE: Nees`ds to be updated if we ever decide to use bb2.  To date, never
					used in actual code, so don't worry too much about it. */

        bb2->min_x = bb1->min_x; bb2->min_y = bb1->min_y; bb2->min_z = bb1->min_z;
        bb2->max_x = bb1->max_x; bb2->max_y = bb1->max_y; bb2->max_z = bb1->max_z;

        bb2->size_x = (bb2->max_x-bb2->min_x)/(bb2->step)+1;
        bb2->size_y = (bb2->max_y-bb2->min_y)/(bb2->step)+1;
        bb2->size_z = (bb2->max_z-bb2->min_z)/(bb2->step)+1;

             bb1->step = 2.*bb1->step+1e-5; /* bb1->step is modified so that we can guarantee
                that all voxels lying inside the original mesh have at least one point in them.  In other
		words, if a voxel is empty, then we can be sure that it is outside the mesh.
		1e-5 is a safety factor for numerical roundoff reasons. We may want to later multiply this by some 
                number > 1 if we decide that we want larger voxels.*/
             bb1->size_x = (bb1->max_x-bb1->min_x)/(bb1->step)+1;
             bb1->size_y = (bb1->max_y-bb1->min_y)/(bb1->step)+1;
             bb1->size_z = (bb1->max_z-bb1->min_z)/(bb1->step)+1;
return;
}
           

/************************************************************************
FUNCTION: generate_tau
 PURPOSE: calculates coefficients for closed Gear scheme, for variable timestep
PARAMETERS:
   Input:
     - nsolve: the order (type) of the solver
     - uvw_old2_exists: == TRUE if uvw_tilde_2 was successfully read from .ini file
     - tconst: scaling factor for times
     - timestep: the current timestep number
     - times: vector of timestep times
   Output:
     - tau0, tau1, tau2: coefficients multiplying {u} vectors in timestepping.
       Specifically, du/dt = tau_0*{u}^n+1 + tau_1*{u}^n + tau_2*{u}^n-1
 RETURNS: void
************************************************************************/
void generate_tau(short nsolve, short uvw_old2_exists, FLOAT tconst, short timestep,
   FLOAT *times, FLOAT *tau0_ptr, FLOAT *tau1_ptr, FLOAT *tau2_ptr)
{
   FLOAT
     dt1,
     dt2;

   /*** Calculate timesteps ***/
   dt1 = tconst*(times[timestep]-times[timestep-1]);
   dt2 = tconst*(times[timestep-1]-times[timestep-2]);
   /* Note that the times vector effectively goes from -1 to Stokes->tot_timesteps.
      This allows the timestep from the last step of a complete cycle to go into
      the -1 spot in the vector.  See newtetr.c (after read_input_pass_0) to see
      how this is done */

   /*** Check for invalid timesteps ***/
   if (dt1 <= 0.0) 
       Message(FATAL,"Computed invalid dt1 = %f\n", dt1);
   if (dt2 <= 0.0 && abs(nsolve) != 1 && uvw_old2_exists) 
       Message(FATAL,"Computed invalid dt2 = %f\n", dt2);

   /*** Calculate coefficients ***/

   if (abs(nsolve) == 1 || !uvw_old2_exists)
   {
      *tau0_ptr = 1.0/dt1;
      *tau1_ptr = -1.0/dt1;
      *tau2_ptr = 0.0;
   }
   else
   {
      *tau0_ptr = (1.0/dt1)*(1.0+2.0*dt1/dt2)/(1.0+dt1/dt2);
      *tau1_ptr = -(1.0/dt1)*(1.0+dt1/dt2);
      *tau2_ptr = (1.0/dt1)*(dt1/dt2)*(dt1/dt2)/(1.0+dt1/dt2);
   }

   return;
}

/************************************************************************
FUNCTION: get_uvw_in
 PURPOSE: obtainss fully-developed inlet velocity conditions.  Note that this
     is not implemented in a very efficient way, but this approach saves a
     great deal of storage in the convective marching routines.
PARAMETERS:
   Input:
     - t: time at which velocity is to be calculated 
     - alpha: Womersley parameter
     - nodal_xyz: positions of nodes
     - velocity: vector to be filled with inlet velocities (in appropriate locations)
     - dvel_dt: vector to be filled with time deriv of inlet velocities (in appropriate locations)
     - Inlets: information about inlets in mesh (number, location, node #s, etc).
     - compute_derivative: TRUE if dvel_dt is to be filled
     - gen_uvw_in_flag: if TRUE, generate velocities analytically
     - coord_direction: indicates if x, y or z component of velocity is to be calculated
     - bc_table    :info about boundary conditions
     - tconst: time scaling factor
   Output:
     - appropriately filled velocity vectors u, v, w
 RETURNS: void
************************************************************************/
void get_uvw_in(FLOAT t, FLOAT alpha, FLOAT **nodal_xyz, VECTOR velocity,
     VECTOR dvel_dt, INLET_INFO *Inlets, short compute_derivative,
     short gen_uvw_in_flag, GLOBAL_TRIPLET coord_direction,
     BC_TYPE *bc_table, FLOAT tconst)
{
   COMPLEX
     za,
     zar,
     zconst,
     zf10,
     zi,      /* = sqrt(-1.0) */
     zt,
     zu;

   NODE
     iloop,
     node;    /* global node number for inlet node */

   FLOAT
     cosa,    /* cosine of angle inlet w.r.t. horizontal axis ?*/
     r,       /* distance from inlet origin to node */
     sina,    /* sine of angle inlet w.r.t. horizontal axis ?*/
     U,       /* inlet velocity magnitude */
     dU_dt,   /* magnitude of d(inlet velocity)/dt */
     x,       /* x displacement from inlet origin to node */
     y,       /* y displacement from inlet origin to node */
     z,       /* z displacement from inlet origin to node */
     norm; /* store the unit normal */

   short
     inlet,   /* inlet number */
     k;       /* Fourier mode of inlet flow waveform */

   #if defined(ANALYTIC)
   FILE
      *file_cfg;

   FLOAT
      a,
      d,
      ax, ay, az,
      dx, dy, dz,
      et;

   if ((file_cfg = fopen("ethier.cfg", "r")) == NULL)
      Message(FATAL," in opening ethier.cfg.");
   fscanf(file_cfg, "%f %f", &a, &d);
   fclose(file_cfg);
   #endif

   t = t - WAVEFORM_PHASE * tconst ;

   if (gen_uvw_in_flag)
   {
      zi = Complex(0.0,1.0);

      for (inlet = 1; inlet <= Inlets->tot_inlets; inlet++)
      {
/*  First normalize the normal   */

       norm = sqrt(Inlets->normal[inlet][1]*Inlets->normal[inlet][1] + 
                   Inlets->normal[inlet][2]*Inlets->normal[inlet][2] +
                   Inlets->normal[inlet][3]* Inlets->normal[inlet][3]); 
       Inlets->normal[inlet][1] /= norm;
       Inlets->normal[inlet][2] /= norm;
       Inlets->normal[inlet][3] /= norm;

         for (iloop = 1; iloop <= Inlets->nodes_per_inlet[inlet]; iloop++)
         {
             node = Inlets->nodes[inlet][iloop];
             x = nodal_xyz[node][1] - Inlets->origin[inlet][1];
             y = nodal_xyz[node][2] - Inlets->origin[inlet][2];
             z = nodal_xyz[node][3] - Inlets->origin[inlet][3];

             #if defined(ANALYTIC)
             x = nodal_xyz[node][1];
             y = nodal_xyz[node][2];
             z = nodal_xyz[node][3];
 
             et = exp(-d*d*t);
             ax = a*x;
             ay = a*y;
             az = a*z;
             dx = d*x;
             dy = d*y;
             dz = d*z;

             if (coord_direction == 1)   /* x direction */
                velocity[node] = -a*(exp(ax)*sin(ay+dz) + exp(az)*cos(ax+dy))*et;
             else if (coord_direction == 2) /* y direction */
                velocity[node] = -a*(exp(ay)*sin(az+dx) + exp(ax)*cos(ay+dz))*et;
             else
                velocity[node] = -a*(exp(az)*sin(ax+dy) + exp(ay)*cos(az+dx))*et;
             /*p = -a*a*et*et*(exp(2.0*ax) + exp(2.0*ay) + exp(2.0*az)
                + 2.0*sin(ax+dy)*cos(az+dx)*exp(ay+az)
                + 2.0*sin(ay+dz)*cos(ax+dy)*exp(az+ax)
                + 2.0*sin(az+dx)*cos(ay+dz)*exp(ax+ay))/2.0 */

             #else

             r = sqrt(x*x+y*y+z*z)/Inlets->radius[inlet];
             U = 2.*(Inlets->zFour[inlet][0].re)*(1.-r*r);
             dU_dt = 0.0;

             for (k = 1; k <= Inlets->Four_per_inlet[inlet]; k++)
             {
             
               /* Note that because the inlet radius may not be equal to the */
               /* radius used in the normalization of the model, we must     */
               /* scale the alpha value appropriately such that the inlet    */
               /* womersley profile is correctly sized for the inlet         */
               /* We do this by multiplying the alpha value by the           */
               /* normalized radius of the inlet                             */
               
               za   = RCmul(alpha*Inlets->radius[inlet]*sqrt(1.0*k), Cmul(zi, Csqrt(zi)));
               zar = RCmul(r, za);
               
               /* calcualting zt in this manner assumes that the t value is  */
               /* to be scaled by tconst (default tconst is 2PI but can be   */
               /* different)                                                 */
               
               zt = Complex(0, 2.*PI*k*t/tconst);
               zf10 = RCmul(2., Cdiv(zbes(1,za), Cmul(zbes(0,za),za)));
               zconst = Cdiv(Cmul(Inlets->zFour[inlet][k], Cexp(zt)),
                  Complex(1. - zf10.re, 0. - zf10.im));
               zu = Cmul(zconst, Complex(1.0 - Cdiv(zbes(0,zar), zbes(0,za)).re,
                                         0.0 - Cdiv(zbes(0,zar), zbes(0,za)).im));
               U += zu.re;
               dU_dt += zu.im * k;
             } 

             /* Normalize wrt inlet radius */
             U /= pow(Inlets->radius[inlet], 2);
             dU_dt /= pow(Inlets->radius[inlet], 2);

            /* To insure that the correct velocity is assigned to the node coordinate */
            /* the boundary conditon must also be checked.  This is done using the    */
            /* bc_table structure.  Note that if bc_table == 0, then the point is     */
            /* a wall point and should be assigned a zero velocity                    */
            
            
             if (coord_direction == 1)        /* x direction */
                {   
                  if (bc_table[node].x == 1)  /* inlet node */
                     velocity [node] = U*Inlets->normal[inlet][1];
                  else                        /*  Wall node or some other BC condition */    
                     velocity[node] = 0.0;   
                } /* end if coord_direction == 1 */
                
             else if (coord_direction == 2)   /* y direction */
               {
                 if (bc_table[node].y == 1)   /* inlet node */
                    velocity[node] = U*Inlets->normal[inlet][2];
                 else                        /*  Wall node or some other BC condition */    
                    velocity[node] = 0.0;   
                     
                } /* end else if coord_direction == 2 */
               
             else                             /* z direction */
               {
                 if (bc_table[node].z == 1)   /* inlet node */
                    velocity[node] = U*Inlets->normal[inlet][3];
                 else                        /*  Wall node or some other BC condition */    
                    velocity[node] = 0.0;   
                    
               } /* end else coord_direction == 3  */
               
             #endif

             if (compute_derivative)
             {
                #if defined(ANALYTIC)
                   dvel_dt[node] = -d*d*velocity[node];
                #else
                   if (coord_direction == 1)
                     dvel_dt[node] = dU_dt*Inlets->normal[inlet][1];
                   else if (coord_direction == 2)
                     dvel_dt[node] = dU_dt*Inlets->normal[inlet][2];
                   else
                     dvel_dt[node] = dU_dt*Inlets->normal[inlet][3];
                #endif
             }

         } /* end loop over inlet nodes */
      } /* end loop over inlets */
   }
   else
   {
      /* iterpolate current and previous inlet velocities if read from uin file.
          Note that non-Dirichlet velocities must be set equal to zero. */
      /*vadd(3*nn,1.0*itss/nsub,uin,1.0-1.0*itss/nsub,uinold, uins) */
      Message(FATAL,".  Inlet velocity file reading not yet implemented.\n");
   }

   return;
}

/************************************************************************
FUNCTION: lag_uvw_in 
 PURPOSE: obtains fully-developed inlet velocity conditions at none nodal postions
    on an inlet plane. (primarily for lag points)
    Note that this is not implemented in a very efficient way, but this 
    approach saves a great deal of storage in the convective marching routines.
PARAMETERS:
   Input:
     - t: time at which velocity is to be calculated 
     - alpha: Womersley parameter
     - lag_x: positions of nodes (x position)
     - lax_y:
     - lag_z:
     - vel_comp: velocity componentto be filled with  inlet velocities (in appropriate locations)
     - Inlets: information about inlets in mesh (number, location, node #s, etc).
     - gen_uvw_in_flag: if TRUE, generate velocities analytically
     - tconst:  time scaling factor used to normalize t
   Output:
     - appropriately filled velocity vectors u, v, w
 RETURNS: void
************************************************************************/
void lag_uvw_in(FLOAT t, FLOAT alpha, FLOAT lag_x, FLOAT lag_y, FLOAT lag_z, NODE iloop,
      VECTOR u, VECTOR v, VECTOR w, NODE inlet, INLET_INFO *Inlets,
      short gen_uvw_in_flag, FLOAT tconst)
{
   COMPLEX
     za,
     zar,
     zconst,
     zf10,
     zi,      /* = sqrt(-1.0) */
     zt,
     zu;

   FLOAT
     cosa,    /* cosine of angle inlet w.r.t. horizontal axis ?*/
     r,       /* distance from inlet origin to node */
     sina,    /* sine of angle inlet w.r.t. horizontal axis ?*/
     U,       /* inlet velocity magnitude */
     x,       /* x displacement from inlet origin to node */
     y,       /* y displacement from inlet origin to node */
     z,       /* z displacement from inlet origin to node */
     norm;    /* store the unit normal */

   short
     k;       /* Fourier mode of inlet flow waveform */

   if (gen_uvw_in_flag)
   {
      zi = Complex(0.0,1.0);

      /*  First normalize the normal   */

       norm = sqrt(Inlets->normal[inlet][1]*Inlets->normal[inlet][1] + 
                   Inlets->normal[inlet][2]*Inlets->normal[inlet][2] +
                   Inlets->normal[inlet][3]*Inlets->normal[inlet][3]); 
       Inlets->normal[inlet][1] /= norm;
       Inlets->normal[inlet][2] /= norm;
       Inlets->normal[inlet][3] /= norm;

             x = lag_x - Inlets->origin[inlet][1];
             y = lag_y - Inlets->origin[inlet][2];
             z = lag_z - Inlets->origin[inlet][3];

             r = sqrt(x*x+y*y+z*z)/Inlets->radius[inlet];
             U = 2.*(Inlets->zFour[inlet][0].re)*(1.-r*r);

             for (k = 1; k <= Inlets->Four_per_inlet[inlet]; k++)
             {
               /* Note that because the inlet radius may not be equal to the */
               /* radius used in the normalization of the model, we must     */
               /* scale the alpha value appropriately such that the inlet    */
               /* womersley profile is correctly sized for the inlet         */
               /* We do this by multiplying the alpha value by the           */
               /* normalized radius of the inlet                             */
               
               za   = RCmul(alpha*Inlets->radius[inlet]*sqrt(1.0*k), Cmul(zi, Csqrt(zi)));
               zar = RCmul(r, za);
               
               /* calcualting zt in this manner assumes that the t value is  */
               /* normalized by tconst (default tconst is 2PI but can be     */
               /* different)                                                 */
               
               zt   = Complex(0, 2.*PI*k*t/tconst);
               zf10 = RCmul(2., Cdiv(zbes(1,za), Cmul(zbes(0,za),za)));
               zconst = Cdiv(Cmul(Inlets->zFour[inlet][k], Cexp(zt)),
                  Complex(1. - zf10.re, 0. - zf10.im));
               zu = Cmul(zconst, Complex(1.0 - Cdiv(zbes(0,zar), zbes(0,za)).re,
                                         0.0 - Cdiv(zbes(0,zar), zbes(0,za)).im));
               U += zu.re;
             } 

             /* Normalize wrt inlet radius */
             U /= pow(Inlets->radius[inlet], 2);

             /* Should not need to check the wall boundary conditon here since */
             /* by design this function is not called unless the lag point is  */
             /* the inlet radius.  This may need checked in the future         */

             u[iloop] = U*Inlets->normal[inlet][1];   /* x direction */
             v[iloop] = U*Inlets->normal[inlet][2];   /* y direction */
             w[iloop] = U*Inlets->normal[inlet][3];   /* z direction */

   }
   else
   {
      /* iterpolate current and previous inlet velocities if read from uin file.
          Note that non-Dirichlet velocities must be set equal to zero. */
      /*vadd(3*nn,1.0*itss/nsub,uin,1.0-1.0*itss/nsub,uinold, uins) */
      Message(FATAL,".  Inlet velocity file reading not yet implemented.\n");
   }

   return;
}



/************************************************************************
FUNCTION: L2_norm
 PURPOSE: computes L2 norm of a vector, looping over non-Dirichlet entries only
PARAMETERS:
   Input:
      - Vector: vector for which nowm is to be computed
      - Condense: info about non-Dirichlet nodes in Vector
  RETURNS:
      - value of L2 norm
************************************************************************/
FLOAT L2_norm(VECTOR Vector, CONDENSE_INFO *Condense)
{

   FLOAT
      norm,
      value;                 /* value of current entry of Vector */

   NODE
      *last_row_ptr,         /* pointer to last non-Dirichlet row number in vectors */
      *row_ptr;              /* points to row */

   last_row_ptr = Condense->rowmap + Condense->tot_rows;
   norm = 0.0;

   /* Loop over all non-Dirichlet rows */
   for (row_ptr = Condense->rowmap + 1; row_ptr <= last_row_ptr; row_ptr++)
   {
      value = *(Vector + *row_ptr);
      norm += value*value;
   }
   return sqrt(norm);
}

/************************************************************************
FUNCTION: Message
 PURPOSE: displays an error/debug message, possibly aborting
PARAMETERS:
   input: message_type (see newtetr.h for macros)
  RETURNS: void
************************************************************************/

void Message(short message_type, char *format, ...)
{
   va_list ap;
   int i;

   /* don't print anything if the debug level is too low */
   if (message_type>=0 && message_type>Debug_flag) return;

   /* initialize variable argument list */
   va_start(ap, format);

   /* print message preamble based on requested message type */
   switch(message_type) {
     case FATAL:
      printf("FATAL ERROR! ");
      break;
     case ERROR:
      printf("Error! ");
      break;
     case WARNING:
      printf("Warning! ");
      break;
     default:
        for (i=1; i <= message_type; i++) printf("   ");
      break;
   }

   /* print out main message */
   (void) vprintf(format, ap);

   /* end argument list */
   va_end(ap);

   /* post process errors */
   switch(message_type) {
     case FATAL:
        printf("...aborting immediately...\n");
      exit(-1);
      break;
     case ERROR:
      Error_count++;
      break;
     default:
      break;
   }
}

/************************************************************************
FUNCTION: pack_vector
 PURPOSE: packs non-Dirichlet entries of a source vctor into target vector
PARAMETERS:
   Input:
     - target: vector to be filled
     - source: vector which supplies entries to target
     - Condense: info about Dirichlet structure of source vector
   Output:
     - appropriately filled target vector
 RETURNS: void
************************************************************************/

void pack_vector(FLOAT *target, FLOAT *source, CONDENSE_INFO *Condense)
{
   NODE
      i;                     /* location in target corresponding to *row_ptr */

   /* Loop over all non-Dirichlet rows */
   for (i=1; i<=Condense->tot_rows; i++)
      target[i] = source[Condense->rowmap[i]];
}

/************************************************************************
FUNCTION: print_ivector
 PURPOSE: prints a NODE type vector
PARAMETERS:
   Input:
     - vec_ptr: pointer to the vector
     - size: size of the vector
 RETURNS:
     - void
************************************************************************/
void print_ivector( NODE *vec_ptr, NODE size, char *title, ... )
{
NODE i; /*  loop variable */
va_list ap;

   /* initialize variable argument list */
   va_start(ap, title);

   (void) vprintf(title, ap);

   /* end argument list */
   va_end(ap);

for( i=1; i<=size; i++ )
printf(" %d \n", vec_ptr[i]);
} 

/************************************************************************
FUNCTION: print_matrix
 PURPOSE: prints a FLOAT type matrix
PARAMETERS:
   Input:
     - mat: pointer to the matrix
     - m_info: info structure for the matrix
 RETURNS:
     - void
************************************************************************/
void print_matrix( MATRIX mat, MATRIX_INFO *m_info, char *title, ... )
{
NODE row, column, /*  loop variables */
     first_column, end_column;
va_list ap;

   /* initialize variable argument list */
   va_start(ap, title);

   (void) vprintf(title, ap);

   /* end argument list */
   va_end(ap);

   for( row=1; row<=m_info->tot_rows; row++ )
   {
   first_column = m_info->entry_offset[row]+1;
   end_column = m_info->entry_offset[row+1];
        for( column=first_column; column<=end_column; column++)
        printf(" [%d][%d] = %f \n", row, m_info->columns[column], mat[column]);
    }     
} 

/************************************************************************
FUNCTION: print_vector
 PURPOSE: prints a NODE type vector
PARAMETERS:
   Input:
     - vec_ptr: pointer to the vector
     - size: size of the vector
 RETURNS:
     - void
************************************************************************/
void print_vector( FLOAT *vec_ptr, NODE size, char *title, ... )
{
NODE i; /*  loop variable */
va_list ap;
 
   /* initialize variable argument list */
   va_start(ap, title);
 
   (void) vprintf(title, ap);
 
   /* end argument list */
   va_end(ap);
 
for( i=1; i<=size; i++ )
printf(" %f \n", vec_ptr[i]);
}

/************************************************************************
FUNCTION: RCmul  (from Numerical Recipes)
 PURPOSE: multiplies a real and a complex #
PARAMETERS:
   Input:
     - x, real #
     - a, complex #
 RETURNS:
     - the product (COMPLEX)
************************************************************************/
COMPLEX RCmul(FLOAT x, COMPLEX a)
{
   COMPLEX
     c;
   c.re=x*a.re;
   c.im=x*a.im;
   return c;
}

/************************************************************************
FUNCTION: shell
 PURPOSE: sorts an array into ascending numerical order while making the corresponding rearrangment
          of another, integer array
PARAMETERS:
   Input:
     - n: the size of the arrays
     - arr: NODE array to be sorted
     - brr: NODE array to be rearranged
   Output:
     - appropriately filled arr and brr
 RETURNS: void
************************************************************************/
void shell(NODE n, NODE *arr, NODE *brr)
{
        int nn,m,j,i,lognb2;
        NODE t, u;
        double ALN2I, TINY;
 
ALN2I = 1.442695022;
TINY = 1.0e-5;
        lognb2=(log((double) n)*ALN2I+TINY);
        m=n;
        for (nn=1;nn<=lognb2;nn++) {
                m >>= 1;
                for (j=m+1;j<=n;j++) {
                        i=j-m;
                        t=arr[j];
                        u=brr[j];
                        while (i >= 1 && arr[i] > t) {
                                arr[i+m]=arr[i];
                                brr[i+m]=brr[i];
                                i -= m;
                        }
                        arr[i+m]=t;
                        brr[i+m]=u;
                }
        }
}



/************************************************************************
FUNCTION: sort2
 PURPOSE: sorts an array into ascending numerical order while making the corresponding rearrangment
          of another, integer array
PARAMETERS:
   Input:
     - n: the size of the arrays
     - arr: FLOAT array to be sorted
     - brr: NODE array to be rearranged
   Output:
     - appropriately filled arr and brr
 RETURNS: void
************************************************************************

void sort2(NODE n, FLOAT arr[], NODE brr[])
{
	NODE i, ir=n, j, k, l=1;
        int *istack, jstack=0;
	FLOAT a, b, temp;

istack = ivector(1, NSTACK);

for ( ; ; ){
    if (ir-1 < M_NUMR ) {
        for (j=l+1; j<=ir; j++) {
             a=arr[j];
             b=brr[j];
             for (i=j-1; i>=l; i--) {
                  if (arr[i] <= a) break;
                  arr[i+1]=arr[i];
                  brr[i+1]=brr[i];
             }
             arr[i+1]=a;
             brr[i+1]=b;
        }
        if (!jstack) {
            free_ivector( istack, 1, NSTACK);
            return;
        }
        ir=istack[jstack];
        l=istack[jstack-1];
        jstack-=2;
     } else {
        k=(l+ir) >> 1;
        SWAP(arr[k], arr[l+1])        
        SWAP(brr[k], brr[l+1])        
        if (arr[l] > arr[ir]) {
            SWAP( arr[l], arr[ir] )
            SWAP( brr[l], brr[ir] )
        }
        if ( arr[l+1] > arr[ir] ) {
             SWAP( arr[l+1], arr[ir] )
             SWAP( brr[l+1], brr[ir] )
        }
        if ( arr[l] > arr[l+1] ) {
             SWAP( arr[l], arr[l+1] )
             SWAP( brr[l], brr[l+1] )
        }
        i=l+1;
        j=ir;
        a=arr[l+1];
        b=brr[l+1];
        for ( ; ; ) {
            do i++; while ( i<=n && arr[i] < a );
            do j--; while ( j>=1 && arr[j] > a );
            if ( j<i ) break;
            SWAP( arr[i], arr[j] ) 
            SWAP( brr[i], brr[j] ) 
        }   
        arr[l+1] = arr[j];
        arr[j] = a;
        brr[l+1] = brr[j];
        brr[j] = b;
        jstack += 2;
        if ( jstack > NSTACK ) nrerror("NSTACK too small in sort2.");
        if ( ir-i+1 >= j-1 ) {
             istack[jstack] = ir;
             istack[jstack-1] = i;
             ir = j-1;
        } else {
            istack[jstack] = j-1;
            istack[jstack-1] = l;
            l= i;
        }
             }
         }
    }
*/
/************************************************************************
FUNCTION: unpack_vector
 PURPOSE: unpacks non-Dirichlet entries of a source vctor into target vector
PARAMETERS:
   Input:
     - target: vector to be filled
     - source: vector which supplies entries to target
     - Condense: info about Dirichlet structure of target vector
   Output:
     - appropriately filled source vector
 RETURNS: void
************************************************************************/

void unpack_vector(FLOAT *target, FLOAT *source, CONDENSE_INFO *Condense)
{
   NODE
      i,                     /* location in target corresponding to *row_ptr */
      *last_row_ptr,         /* pointer to last non-Dirichlet row number in vectors */
      *row_ptr;              /* points to row */

   last_row_ptr = Condense->rowmap + Condense->tot_rows;

   for (row_ptr = Condense->rowmap + 1, i=1; row_ptr <= last_row_ptr; row_ptr++, i++)
      target[*row_ptr] = source[i];

}

/************************************************************************
FUNCTION: usage
 PURPOSE: displays summary of program and options
PARAMETERS:
  RETURNS: void
************************************************************************/

void usage(void)
{
  printf("\nUsage: newtetr [-dn] [-fn] [-l] [-pC|J] [-tn] input_file\n\n"
         "-dn: dump debugging info, where n controls level of debugging:\n"
         "\tn = 1 - milestones in problem setup, time marching\n"
         "\tn = 2 - level 1 plus extra time marching information\n"
         "\tn = 3 - level 2 plus iterative solver information\n"
         "\tWarning: amount of info at -d2 and -d3 is potentially HUGE\n"
         "-fn: ILUT fill level [implicit scheme] (default = %d)\n"
         "-lfilename: file for logging screen i/o (useful for Windows;\n"
         "\tfor Unix use \"newtetr > filename\").  If no filename\n"
         "\tspecified, default dumpfile is tetr.log\n"
         "-pC: preconditioner type for mass matrix\n"
         "\tJ = Jacobi (default), C = Incomplete Cholesky order 0\n"
         "-tn: ILUT threshold level [implicit scheme] (default = %.8g)\n",
        DEFAULT_ILU_FILL_LEVEL, DEFAULT_ILU_TOL);
  exit(0);
}

/************************************************************************
FUNCTION: zbes
 PURPOSE: calculates the complex Bessel function Jn(y)
PARAMETERS:
   Input:
     - n, the order of the Bessel fcn
     - the argument (COMPLEX)
 RETURNS:
     - the value of the bessel fcn (COMPLEX)
************************************************************************/
COMPLEX zbes(short n, COMPLEX argument)
{
   COMPLEX
     z = {1.0, 0.0},
     zanswer,
     zarg,
     zproduct = {1.0, 0.0};

   unsigned
     i;

   zarg = RCmul(-0.25, Cmul(argument, argument));
   zanswer.re = 1.0;
   zanswer.im = 0.0;

   for (i=1; i <= 10000; i++)
   {
     z = RCmul(1./((float)i * (float)(i+n)), Cmul(z, zarg));
     if (Cabs(z) <= 1.e-20)
       break;
     zanswer = Cadd(zanswer,z);
   }

   for (i=1; i <= n; zproduct = Cmul(zproduct, RCmul(0.5, argument)), i++);

   zanswer = Cmul(zanswer, zproduct);

   return zanswer;
}

/************************************************************************
FUNCTION: zero_vector
 PURPOSE: fills a vector with zeroes
PARAMETERS:
   Input:
     - the vector
     - the length of the vector
   Output:
     - appropriately zeroed vector
 RETURNS: void
************************************************************************/

void zero_vector(FLOAT *vector, NODE length)
{
   NODE iloop;
   for (iloop=1; iloop <= length; vector[iloop++] = 0.0);
}

/************************************************************************
FUNCTION: zero_vector_nonD
 PURPOSE: fills non-Dirichlet entries of a vector with zeroes
PARAMETERS:
   Input:
     - vector: vector to be zeroed
     - Condense: info about Dirichlet structure of vector
   Output:
     - appropriately zeroed vector
 RETURNS: void
************************************************************************/

void zero_vector_nonD(FLOAT *vector, CONDENSE_INFO *Condense)
{
   NODE
      *last_row_ptr,         /* pointer to last non-Dirichlet row number in vectors */
      *row_ptr;              /* points to row */

   last_row_ptr = Condense->rowmap + Condense->tot_rows;

   /* Loop over all non-Dirichlet rows */
   for (row_ptr = Condense->rowmap + 1; row_ptr <= last_row_ptr; row_ptr++)
      *(vector + *row_ptr) = 0.0;

}


