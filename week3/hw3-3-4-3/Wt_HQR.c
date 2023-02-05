#include "FLAME.h"
#include "Dot_Column_v_unb.h"
/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/

int Wt_HQR_unb( FLA_Obj at, FLA_Obj b, FLA_Obj C, FLA_Obj alpha, FLA_Obj wt )
{
  FLA_Obj aLt,    aRt,       a0t,  alpha1,  a2t;

  FLA_Obj CL,    CR,       C0,  c1,  C2;

  FLA_Obj wLt,    wRt,       w0t,  omega1,  w2t;

  FLA_Obj temp_alpha;
  double* ptr_temp;
  double alphain;
  FLA_Obj_create(FLA_DOUBLE,1,1,0,0,&temp_alpha);
  alphain = 1.0/(*( double * ) FLA_Obj_buffer_at_view(alpha));
  ptr_temp = ( double * ) FLA_Obj_buffer_at_view(temp_alpha);
  
  FLA_Part_1x2( at,    &aLt,  &aRt,      0, FLA_LEFT );

  FLA_Part_1x2( C,    &CL,  &CR,      0, FLA_LEFT );

  FLA_Part_1x2( wt,    &wLt,  &wRt,      0, FLA_LEFT );

  while ( FLA_Obj_width( aLt ) < FLA_Obj_width( at ) ){

    FLA_Repart_1x2_to_1x3( aLt,  /**/ aRt,        &a0t, /**/ &alpha1, &a2t,
                           1, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, /**/ &c1, &C2,
                           1, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( wLt,  /**/ wRt,        &w0t, /**/ &omega1, &w2t,
                           1, FLA_RIGHT );

    /*------------------------------------------------------------*/
    Dot_Column_v_unb(b,c1,temp_alpha);
    *( double * ) FLA_Obj_buffer_at_view(omega1) = (*( double * ) FLA_Obj_buffer_at_view(alpha1)+*ptr_temp)*alphain;
    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &aLt,  /**/ &aRt,        a0t, alpha1, /**/ a2t,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, c1, /**/ C2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &wLt,  /**/ &wRt,        w0t, omega1, /**/ w2t,
                              FLA_LEFT );

  }
  FLA_Obj_free(&temp_alpha);
  return FLA_SUCCESS;
}
