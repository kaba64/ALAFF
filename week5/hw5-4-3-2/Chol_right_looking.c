#include <math.h>
#include <FLAME.h>
#include "G_VScale.h"
#include "M_Sym_Rank1_Update.h"

/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/

int Chol_right_looking_unb( FLA_Obj A )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );
    
    /*------------------------------------------------------------*/
    *( double *) FLA_Obj_buffer_at_view(alpha11) = sqrt(*( double *) FLA_Obj_buffer_at_view(alpha11));
    G_VScale_unb(a21,alpha11,1,'Y',a21);
    M_Sym_Rank1_Update_unb(A22,a21,a21,-1);
    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );

  }

  return FLA_SUCCESS;
}
