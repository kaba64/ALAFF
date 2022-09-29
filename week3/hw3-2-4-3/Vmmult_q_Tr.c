#include <FLAME.h>
#include "Dot_Column_v_unb.h"

/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/

int Vmmult_q_Tr_unb( FLA_Obj a, FLA_Obj A, FLA_Obj rt )
{
  FLA_Obj AL,    AR,       A0,  a1,  A2;
  
  FLA_Obj rLt,    rRt,       r0t,  rho1,  r2t;
  
  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );
  
  FLA_Part_1x2( rt,    &rLt,  &rRt,      0, FLA_LEFT );
  
  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &a1, &A2,
                           1, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( rLt,  /**/ rRt,        &r0t, /**/ &rho1, &r2t,
                           1, FLA_RIGHT );

    /*------------------------------------------------------------*/
    Dot_Column_v_unb(a,a1,rho1);
    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, a1, /**/ A2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &rLt,  /**/ &rRt,        r0t, rho1, /**/ r2t,
                              FLA_LEFT );

  }

  return FLA_SUCCESS;
}
