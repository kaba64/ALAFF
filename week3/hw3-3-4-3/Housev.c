#include <FLAME.h>
#include "Norm_2_unb.h"
#include <math.h>
#include "axpy_inverse.h"
#include "zero_vector.h"
/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/

int Housev_unb(FLA_Obj beta, FLA_Obj u, FLA_Obj tau )
{
  FLA_Obj norm2_rest, visc1;
  FLA_Obj  alpha;
  double a1, a2;
  FLA_Obj_show( " beta = [ ", beta, "%le", " ];" );
  FLA_Obj_show( " u = [ ", u, "%le", " ];" );
  FLA_Obj_create(FLA_DOUBLE,1,1,0,0,&alpha);
  FLA_Obj_create(FLA_DOUBLE,1,1,0,0,&norm2_rest);
  FLA_Obj_create(FLA_DOUBLE,1,1,0,0,&visc1);
  *(double*)FLA_Obj_buffer_at_view(alpha) = *(double*)FLA_Obj_buffer_at_view(beta);
  /*------------------------------------------------------------*/
  
  Norm_2_unb(u,norm2_rest);
  a1 = (*(double*)FLA_Obj_buffer_at_view(norm2_rest))*(*(double*)FLA_Obj_buffer_at_view(norm2_rest));
  a2 = (*(double*)FLA_Obj_buffer_at_view(alpha))*(*(double*)FLA_Obj_buffer_at_view(alpha));
  if(*(double*)FLA_Obj_buffer_at_view(alpha)>=0.0){
    *(double*)FLA_Obj_buffer_at_view(beta) = -sqrt(a1+a2);
  }else{
    *(double*)FLA_Obj_buffer_at_view(beta) = sqrt(a1+a2);
  }
  *(double*)FLA_Obj_buffer_at_view(visc1) = *(double*)FLA_Obj_buffer_at_view(alpha) - *(double*)FLA_Obj_buffer_at_view(beta);
  if(*(double*)FLA_Obj_buffer_at_view(visc1)==0.0){
    Zero_Vector_unb(u);
    a1 = 0.0;
  }else{
    *(double*)FLA_Obj_buffer_at_view(norm2_rest)/=fabs(*(double*)FLA_Obj_buffer_at_view(visc1));
    AXPY_Inverse_unb(u,visc1,u);
    a1 = *(double*)FLA_Obj_buffer_at_view(norm2_rest)*(*(double*)FLA_Obj_buffer_at_view(norm2_rest));
  }
  *(double*)FLA_Obj_buffer_at_view(tau) = (1.0+a1)*0.5;
  /*------------------------------------------------------------*/
  FLA_Obj_free(&norm2_rest);
  FLA_Obj_free(&visc1);
  FLA_Obj_free(&alpha);
  
  return FLA_SUCCESS;
}
