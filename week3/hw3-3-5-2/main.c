#include <stdio.h>
#include <sys/time.h>
#include <FLAME.h>
#include "Housev.h"
#include "Update_HQR_at12.h"
#include "H_Rank1_Update.h"
#include "Wt_HQR.h"
#include "HQR.h"
#include "FormQ.h"


#define PRINT_DATA
// Declaration of local prototypes.
void matrix_generate( FLA_Obj A, char type,char name);
double spendtimeinsecond();

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
  FormQ_unb(A,t);
#ifdef PRINT_DATA
  FLA_Obj_show( " Q = [ ", A, "%le", " ];" );
#endif
  printf("The spent time for computing the QR factorizatio of A in second : %lf\n",t2-t1);
  FLA_Obj_free(&A);
  FLA_Obj_free(&t);
  FLA_Finalize();
  return 0;
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
