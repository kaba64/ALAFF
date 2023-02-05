#include <FLAME.h>
#include "Norm_2_unb.h"
#include <math.h>
#include "axpy_inverse.h"
/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/

int Housev_unb( FLA_Obj alpha, FLA_Obj x, FLA_Obj beta, FLA_Obj u, FLA_Obj tau )
{
  FLA_Obj norm2_rest, visc1;
  double a1, a2;
  FLA_Obj_create(FLA_DOUBLE,1,1,0,0,&norm2_rest);
  FLA_Obj_create(FLA_DOUBLE,1,1,0,0,&visc1);
  /*------------------------------------------------------------*/
  Norm_2_unb(x,norm2_rest);
  a1 = (*(double*)FLA_Obj_buffer_at_view(norm2_rest))*(*(double*)FLA_Obj_buffer_at_view(norm2_rest));
  a2 = (*(double*)FLA_Obj_buffer_at_view(alpha))*(*(double*)FLA_Obj_buffer_at_view(alpha));
  if(*(double*)FLA_Obj_buffer_at_view(alpha)>=0.0){
    *(double*)FLA_Obj_buffer_at_view(beta) = -sqrt(a1+a2);
  }else{
    *(double*)FLA_Obj_buffer_at_view(beta) = sqrt(a1+a2);
  }
  *(double*)FLA_Obj_buffer_at_view(visc1) = *(double*)FLA_Obj_buffer_at_view(alpha) - *(double*)FLA_Obj_buffer_at_view(beta);
  AXPY_Inverse_unb(x,visc1,u);
  *(double*)FLA_Obj_buffer_at_view(norm2_rest)/=abs(*(double*)FLA_Obj_buffer_at_view(visc1));
  a1 = *(double*)FLA_Obj_buffer_at_view(norm2_rest)*(*(double*)FLA_Obj_buffer_at_view(norm2_rest));
  *(double*)FLA_Obj_buffer_at_view(tau) = (1+a1)*0.5;
  /*------------------------------------------------------------*/
  FLA_Obj_free(&norm2_rest);
  FLA_Obj_free(&visc1);
  //}

  return FLA_SUCCESS;
}
