#include <stdio.h>
#include <sys/time.h>
#include <FLAME.h>
#include "Housev.h"
#include "Update_HQR_at12.h"
#include "H_Rank1_Update.h"
#include "Wt_HQR.h"

#define PRINT_DATA
// Declaration of local prototypes.
void matrix_generate( FLA_Obj A, char type,char name);
double spendtimeinsecond();
int HQR_unb_var1( FLA_Obj A, FLA_Obj t );

int main( int argc, char *argv[] ) {
  int m_A,n_A;
  FLA_Obj A,t;
  double t1,t2;
  printf("Please enter the size of the matrix A : m_A and n_A \n");
  scanf("%d%d",&m_A,&n_A);
  FLA_Init();
  FLA_Obj_create(FLA_DOUBLE,m_A,n_A,0,0,&A);
  FLA_Obj_create(FLA_DOUBLE,n_A,1,0,0,&t);
  matrix_generate(A,'F','A');
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
  FLA_Obj_show( " t = [ ", t, "%le", " ];" );
#endif
  t1=spendtimeinsecond();
  HQR_unb_var1(A,t);
  t2=spendtimeinsecond();
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
  FLA_Obj_show( " t = [ ", t, "%le", " ];" );
#endif
  printf("The spent time for computing the QR factorizatio of A in second : %lf\n",t2-t1);
  FLA_Obj_free(&A);
  FLA_Obj_free(&t);
  FLA_Finalize();
  return 0;
}
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
