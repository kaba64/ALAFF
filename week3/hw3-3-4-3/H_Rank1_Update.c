#include "FLAME.h"
#include "Rank1-update.h"

/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/

int H_Rank1_Update_unb( FLA_Obj A, FLA_Obj b, FLA_Obj wt )
{
  FLA_Obj AL,    AR,       A0,  a1,  A2;

  FLA_Obj wLt,    wRt,       w0t,  omega1,  w2t;

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  FLA_Part_1x2( wt,    &wLt,  &wRt,      0, FLA_LEFT );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &a1, &A2,
                           1, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( wLt,  /**/ wRt,        &w0t, /**/ &omega1, &w2t,
                           1, FLA_RIGHT );

    /*------------------------------------------------------------*/
    Rank1_update_unb(omega1,b,a1);
    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, a1, /**/ A2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &wLt,  /**/ &wRt,        w0t, omega1, /**/ w2t,
                              FLA_LEFT );

  }

  return FLA_SUCCESS;
}
