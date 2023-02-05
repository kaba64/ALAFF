#include <stdio.h>
#include <sys/time.h>
#include <FLAME.h>
#include "Mvmult_A_Tr_unb.h"
#include "Mvmult_QR_unb.h"
#include "Norm_2_unb.h"
#include "axpy_inverse.h"

#define PRINT_DATA
// Declaration of local prototypes.
void matrix_generate( FLA_Obj A, char type,char name);
double spendtimeinsecond();
int CGS_QR( FLA_Obj A, FLA_Obj Q, FLA_Obj R );

int main( int argc, char *argv[] ) {
  int m_A,n_A;
  FLA_Obj A, Q,R;
  double t1,t2;
  printf("Please enter the size of the matrix A : m_A and n_A \n");
  scanf("%d%d",&m_A,&n_A);
  FLA_Init();
  FLA_Obj_create(FLA_DOUBLE,m_A,n_A,0,0,&A);
  FLA_Obj_create(FLA_DOUBLE,m_A,n_A,0,0,&Q);
  FLA_Obj_create(FLA_DOUBLE,n_A,n_A,0,0,&R);
  matrix_generate(A,'F','A');
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
  FLA_Obj_show( " Q = [ ", Q, "%le", " ];" );
  FLA_Obj_show( " R = [ ", R, "%le", " ];" );
#endif
  t1=spendtimeinsecond();
  CGS_QR(A,Q,R);
  t2=spendtimeinsecond();
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
  FLA_Obj_show( " Q = [ ", Q, "%le", " ];" );
  FLA_Obj_show( " R = [ ", R, "%le", " ];" );
#endif
  printf("The spent time for computing the QR factorizatio of A in second : %lf\n",t2-t1);
  FLA_Obj_free(&A);
  FLA_Obj_free(&Q);
  FLA_Obj_free(&R);
  FLA_Finalize();
  return 0;
}

/* Copyright 2022 The University of Texas at Austin  
   
   For licensing information see
   http://www.cs.utexas.edu/users/flame/license.html 
   
   Programmed by: Kazem Bazesefidpar
   kazemba@kth.se
*/

int CGS_QR( FLA_Obj A, FLA_Obj Q, FLA_Obj R )
{
  FLA_Obj AL,    AR,       A0,  a1,  A2;

  FLA_Obj QL,    QR,       Q0,  q1,  Q2;

  FLA_Obj RTL,   RTR,      R00,  r01,   R02, 
          RBL,   RBR,      r10t, rho11, r12t,
                           R20,  r21,   R22;
  /* A temporarily column vector */
  FLA_Obj a_orth;
  int  m_orth;
  m_orth = FLA_Obj_length( A );
  FLA_Obj_create(FLA_DOUBLE,m_orth,1,0,0,&a_orth);
  /*end*/
  
  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  FLA_Part_1x2( Q,    &QL,  &QR,      0, FLA_LEFT );

  FLA_Part_2x2( R,    &RTL, &RTR,
                      &RBL, &RBR,     0, 0, FLA_TL );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &a1, &A2,
                           1, FLA_RIGHT );

    FLA_Repart_1x2_to_1x3( QL,  /**/ QR,        &Q0, /**/ &q1, &Q2,
                           1, FLA_RIGHT );

    FLA_Repart_2x2_to_3x3( RTL, /**/ RTR,       &R00,  /**/ &r01,   &R02,
                        /* ************* */   /* ************************ */
                                                &r10t, /**/ &rho11, &r12t,
                           RBL, /**/ RBR,       &R20,  /**/ &r21,   &R22,
                           1, 1, FLA_BR );

    /*------------------------------------------------------------*/
    Mvmult_A_Tr_unb(Q0,a1,r01);
    Mvmult_QR_unb(Q0,r01,a1,a_orth);
    Norm_2_unb(a_orth,rho11);
    AXPY_Inverse_unb(a_orth,rho11,q1);
    /*------------------------------------------------------------*/

    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, a1, /**/ A2,
                              FLA_LEFT );

    FLA_Cont_with_1x3_to_1x2( &QL,  /**/ &QR,        Q0, q1, /**/ Q2,
                              FLA_LEFT );

    FLA_Cont_with_3x3_to_2x2( &RTL, /**/ &RTR,       R00,  r01,   /**/ R02,
                                                     r10t, rho11, /**/ r12t,
                            /* ************** */  /* ********************** */
                              &RBL, /**/ &RBR,       R20,  r21,   /**/ R22,
                              FLA_TL );

  }
  FLA_Obj_free(&a_orth);
  return FLA_SUCCESS;
}
//
void matrix_generate( FLA_Obj A, char type,char name) {
  double  * buff_A;
  int  m_A, n_A, ldim_A;
  int  num;
  double xin;
  
  buff_A = ( double * ) FLA_Obj_buffer_at_view( A );
  m_A    = FLA_Obj_length( A );
  n_A    = FLA_Obj_width ( A );
  ldim_A = FLA_Obj_col_stride( A );
  num = 1;
  if(type=='F'){
    for (int j = 0; j < n_A; j++ ) {
      for (int i = 0; i < m_A; i++ ) {
	printf("Please enter %c(%d,%d) \t",name,i,j);
	scanf("%lf",&buff_A[ i + j * ldim_A ]);
        //buff_A[ i + j * ldim_A ] = ( double ) num;
        //num++;
      }
    }
  }else if(type=='S'){
    for (int j = 0; j < n_A; j++ ) {
      for (int i = 0; i < m_A; i++ ) {
	if(i>=j){
	  buff_A[ i + j * ldim_A ] = ( double ) num;
	  num++;
	}else{
	  buff_A[ i + j * ldim_A ]=buff_A[j+i*ldim_A];
	}
      }
    }
  }else if(type=='U'){
    for (int j = 0; j < n_A; j++ ) {
      for (int i = 0; i < j+1; i++ ) {
        buff_A[ i + j * ldim_A ] = ( double ) num;
        num++;
      }
    }
  }else if(type=='L'){
    for (int j = 0; j < n_A; j++ ) {
      for (int i = j; i < m_A; i++ ) {
        buff_A[ i + j * ldim_A ] = ( double ) num;
        num++;
      }
    }
  }else if(type=='I'){
    for (int j = 0; j < n_A; j++ ) {
      for (int i = 0; i < m_A; i++ ) {
	if(i==j){
	  buff_A[ i + j * ldim_A ] = 1.0;
        }
      }
    }
  }
}
//
double spendtimeinsecond()
{
  struct timeval tp;
  struct timezone tzp;
  int i;

  i = gettimeofday(&tp,&tzp);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}
