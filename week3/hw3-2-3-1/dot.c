#include <FLAME.h>
#include "dot.h"
/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/
int Dot_unb( FLA_Obj alpha, FLA_Obj xt, FLA_Obj y )
{
  FLA_Obj xLt,    xRt,       x0t,  chi1,  x2t;

  FLA_Obj yT,              y0,
          yB,              psi1,
                           y2;
  double* alphain;
  alphain = ( double * ) FLA_Obj_buffer_at_view(alpha);
  
  FLA_Part_1x2( xt,    &xLt,  &xRt,      0, FLA_LEFT );

  FLA_Part_2x1( y,    &yT, 
                      &yB,            0, FLA_TOP );

  while ( FLA_Obj_width( xLt ) < FLA_Obj_width( xt ) ){

    FLA_Repart_1x2_to_1x3( xLt,  /**/ xRt,        &x0t, /**/ &chi1, &x2t,
                           1, FLA_RIGHT );

    FLA_Repart_2x1_to_3x1( yT,                &y0, 
                        /* ** */            /* **** */
                                              &psi1, 
                           yB,                &y2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/
    *alphain += (*(double*)FLA_Obj_buffer_at_view(chi1))*(*(double*)FLA_Obj_buffer_at_view(psi1));
    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &xLt,  /**/ &xRt,        x0t, chi1, /**/ x2t,
                              FLA_LEFT );

    FLA_Cont_with_3x1_to_2x1( &yT,                y0, 
                                                  psi1, 
                            /* ** */           /* **** */
                              &yB,                y2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}
