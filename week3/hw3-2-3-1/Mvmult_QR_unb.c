#include "FLAME.h"
#include "dot.h"

/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/
int Mvmult_QR_unb( FLA_Obj Q, FLA_Obj r, FLA_Obj a, FLA_Obj b )
{
  FLA_Obj QT,              Q0,
          QB,              q1t,
                           Q2;

  FLA_Obj aT,              a0,
          aB,              alpha1,
                           a2;

  FLA_Obj bT,              b0,
          bB,              beta1,
                           b2;
  /* A temporarily alpha */
  FLA_Obj alpha;
  FLA_Obj_create(FLA_DOUBLE,1,1,0,0,&alpha);
  *(double*) FLA_Obj_buffer_at_view(alpha) = 0.0;
  /*end*/
  
  FLA_Part_2x1( Q,    &QT, 
                      &QB,            0, FLA_TOP );

  FLA_Part_2x1( a,    &aT, 
                      &aB,            0, FLA_TOP );

  FLA_Part_2x1( b,    &bT, 
                      &bB,            0, FLA_TOP );

  while ( FLA_Obj_length( QT ) < FLA_Obj_length( Q ) ){

    FLA_Repart_2x1_to_3x1( QT,                &Q0, 
                        /* ** */            /* *** */
                                              &q1t, 
                           QB,                &Q2,        1, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( aT,                &a0, 
                        /* ** */            /* ****** */
                                              &alpha1, 
                           aB,                &a2,        1, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( bT,                &b0, 
                        /* ** */            /* ***** */
                                              &beta1, 
                           bB,                &b2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/
    Dot_unb(alpha,q1t,r);
    *(double*) FLA_Obj_buffer_at_view(beta1)=*(double*) FLA_Obj_buffer_at_view(alpha1)-*(double*) FLA_Obj_buffer_at_view(alpha);
    *(double*) FLA_Obj_buffer_at_view(alpha) = 0.0;
    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &QT,                Q0, 
                                                  q1t, 
                            /* ** */           /* *** */
                              &QB,                Q2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &aT,                a0, 
                                                  alpha1, 
                            /* ** */           /* ****** */
                              &aB,                a2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &bT,                b0, 
                                                  beta1, 
                            /* ** */           /* ***** */
                              &bB,                b2,     FLA_TOP );

  }
  FLA_Obj_free(&alpha);
  return FLA_SUCCESS;
}
