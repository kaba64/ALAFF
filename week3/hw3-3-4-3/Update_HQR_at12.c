#include "FLAME.h"

/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/

int Update_HQR_at12_unb( FLA_Obj at, FLA_Obj wt )
{
  FLA_Obj aLt,    aRt,       a0t,  alpha1,  a2t;

  FLA_Obj wLt,    wRt,       w0t,  omega1,  w2t;

  FLA_Part_1x2( at,    &aLt,  &aRt,      0, FLA_LEFT );

  FLA_Part_1x2( wt,    &wLt,  &wRt,      0, FLA_LEFT );

  while ( FLA_Obj_width( aLt ) < FLA_Obj_width( at ) ){

    FLA_Repart_1x2_to_1x3( aLt,  /**/ aRt,        &a0t, /**/ &alpha1, &a2t,
                           1, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( wLt,  /**/ wRt,        &w0t, /**/ &omega1, &w2t,
                           1, FLA_RIGHT );

    /*------------------------------------------------------------*/
    *( double * ) FLA_Obj_buffer_at_view(alpha1)-= *( double * ) FLA_Obj_buffer_at_view(omega1); 
    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &aLt,  /**/ &aRt,        a0t, alpha1, /**/ a2t,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &wLt,  /**/ &wRt,        w0t, omega1, /**/ w2t,
                              FLA_LEFT );

  }

  return FLA_SUCCESS;
}
