#include <FLAME.h>
#include "Housev.h"
#include "Update_HQR_at12.h"
#include "H_Rank1_Update.h"
#include "Wt_HQR.h"

/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/

int HQR_unb_var1( FLA_Obj A, FLA_Obj t )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj tT,              t0,
          tB,              tau1,
                           t2;
  FLA_Obj wt;
  int n_w;
  
  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x1( t,    &tT, 
                      &tB,            0, FLA_TOP );

  while ( FLA_Obj_width( ATL ) < FLA_Obj_width( A ) ){

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    FLA_Repart_2x1_to_3x1( tT,                &t0, 
                        /* ** */            /* **** */
                                              &tau1, 
                           tB,                &t2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/
    Housev_unb(alpha11,a21,tau1);
    n_w = FLA_Obj_width (a12t);
    FLA_Obj_create(FLA_DOUBLE,1,n_w,0,0,&wt);
    Wt_HQR_unb(a12t,a21,A22,tau1,wt);
    Update_HQR_at12_unb(a12t,wt);
    H_Rank1_Update_unb(A22,a21,wt);
    FLA_Obj_free(&wt);
    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &tT,                t0, 
                                                  tau1, 
                            /* ** */           /* **** */
                              &tB,                t2,     FLA_TOP );
    
  }

  return FLA_SUCCESS;
}
