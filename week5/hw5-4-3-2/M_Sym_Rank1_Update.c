#include <FLAME.h>
#include "G_AXPY.h"
/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/
int M_Sym_Rank1_Update_unb( FLA_Obj A, FLA_Obj x, FLA_Obj y, int sign)
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj xT,              x0,
          xB,              chi1,
                           x2;

  FLA_Obj yT,              y0,
          yB,              psi1,
                           y2;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_BR );

  FLA_Part_2x1( x,    &xT, 
                      &xB,            0, FLA_BOTTOM );

  FLA_Part_2x1( y,    &yT, 
                      &yB,            0, FLA_BOTTOM );

  while ( FLA_Obj_length( ABR ) < FLA_Obj_length( A ) ){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  &a01,     /**/ &A02,
                                                &a10t, &alpha11, /**/ &a12t,
                        /* ************* */   /* ************************** */
                           ABL, /**/ ABR,       &A20,  &a21,     /**/ &A22,
                           1, 1, FLA_TL );

    FLA_Repart_2x1_to_3x1( xT,                &x0, 
                                              &chi1, 
                        /* ** */            /* **** */
                           xB,                &x2,        1, FLA_TOP );

    FLA_Repart_2x1_to_3x1( yT,                &y0, 
                                              &psi1, 
                        /* ** */            /* **** */
                           yB,                &y2,        1, FLA_TOP );

    /*------------------------------------------------------------*/
    *(double*)FLA_Obj_buffer_at_view(alpha11)+=sign*((*(double*)FLA_Obj_buffer_at_view(psi1))*(*(double*)FLA_Obj_buffer_at_view(chi1)));
    G_AXPY_unb(a21,x2,psi1,sign);
    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  /**/ a01,     A02,
                            /* ************** */  /* ************************ */
                                                     a10t, /**/ alpha11, a12t,
                              &ABL, /**/ &ABR,       A20,  /**/ a21,     A22,
                              FLA_BR );

    FLA_Cont_with_3x1_to_2x1( &xT,                x0, 
                            /* ** */           /* **** */
                                                  chi1, 
                              &xB,                x2,     FLA_BOTTOM );

    FLA_Cont_with_3x1_to_2x1( &yT,                y0, 
                            /* ** */           /* **** */
                                                  psi1, 
                              &yB,                y2,     FLA_BOTTOM );

  }

  return FLA_SUCCESS;
}
