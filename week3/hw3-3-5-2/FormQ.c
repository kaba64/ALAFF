#include <FLAME.h>
#include "XT_A_S_N.h"
#include "M_Rank1_Update.h"
#include "G_VScale.h"

/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/

int FormQ_unb( FLA_Obj A, FLA_Obj t )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj tT,              t0,
          tB,              tau1,
                           t2;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_BR );

  FLA_Part_2x1( t,    &tT, 
                      &tB,            0, FLA_BOTTOM );

  while ( FLA_Obj_width(ATL)>0){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  &a01,     /**/ &A02,
                                                &a10t, &alpha11, /**/ &a12t,
                        /* ************* */   /* ************************** */
                           ABL, /**/ ABR,       &A20,  &a21,     /**/ &A22,
                           1, 1, FLA_TL );

    FLA_Repart_2x1_to_3x1( tT,                &t0, 
                                              &tau1, 
                        /* ** */            /* **** */
                           tB,                &t2,        1, FLA_TOP );

    /*------------------------------------------------------------*/
    *(double*)FLA_Obj_buffer_at_view(alpha11) = 1.0 - (1.0/(*(double*)FLA_Obj_buffer_at_view(tau1)));
    XT_A_S_N_unb(a12t,a21,A22,tau1);
    M_Rank1_Update_unb(A22,a21,a12t,1);
    G_VScale_unb(a21,tau1,-1,'Y',a21);
    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  /**/ a01,     A02,
                            /* ************** */  /* ************************ */
                                                     a10t, /**/ alpha11, a12t,
                              &ABL, /**/ &ABR,       A20,  /**/ a21,     A22,
                              FLA_BR );

    FLA_Cont_with_3x1_to_2x1( &tT,                t0, 
                            /* ** */           /* **** */
                                                  tau1, 
                              &tB,                t2,     FLA_BOTTOM );

  }

  return FLA_SUCCESS;
}
