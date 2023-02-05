#include <FLAME.h>
#include "Dot_Column_v_unb.h"
/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/

int Mvmult_A_Tr_unb( FLA_Obj A, FLA_Obj x, FLA_Obj y )
{
  FLA_Obj AL,    AR,       A0,  a1,  A2;

  FLA_Obj yT,              y0,
          yB,              psi1,
                           y2;

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  FLA_Part_2x1( y,    &yT, 
                      &yB,            0, FLA_TOP );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &a1, &A2,
                           1, FLA_RIGHT );

    FLA_Repart_2x1_to_3x1( yT,                &y0, 
                        /* ** */            /* **** */
                                              &psi1, 
                           yB,                &y2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/
    Dot_Column_v_unb(a1,x,psi1);
    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, a1, /**/ A2,
                              FLA_LEFT );

    FLA_Cont_with_3x1_to_2x1( &yT,                y0, 
                                                  psi1, 
                            /* ** */           /* **** */
                              &yB,                y2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}
