#include "FLAME.h"

/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/

int AXPY_Inverse_unb( FLA_Obj a, FLA_Obj alpha, FLA_Obj q )
{
  FLA_Obj aT,              a0,
          aB,              alpha1,
                           a2;

  FLA_Obj qT,              q0,
          qB,              chi1,
                           q2;
  double* alphain;
  alphain = ( double * ) FLA_Obj_buffer_at_view(alpha);
  
  FLA_Part_2x1( a,    &aT, 
                      &aB,            0, FLA_TOP );

  FLA_Part_2x1( q,    &qT, 
                      &qB,            0, FLA_TOP );

  while ( FLA_Obj_length( aT ) < FLA_Obj_length( a ) ){

    FLA_Repart_2x1_to_3x1( aT,                &a0, 
                        /* ** */            /* ****** */
                                              &alpha1, 
                           aB,                &a2,        1, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( qT,                &q0, 
                        /* ** */            /* **** */
                                              &chi1, 
                           qB,                &q2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/
    (*(double*)FLA_Obj_buffer_at_view(chi1))= (*(double*)FLA_Obj_buffer_at_view(alpha1))/(*alphain);
    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &aT,                a0, 
                                                  alpha1, 
                            /* ** */           /* ****** */
                              &aB,                a2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &qT,                q0, 
                                                  chi1, 
                            /* ** */           /* **** */
                              &qB,                q2,     FLA_TOP );

  }

  return FLA_SUCCESS;
}
