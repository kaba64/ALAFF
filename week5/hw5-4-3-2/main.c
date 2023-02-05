#include <sys/time.h>
#include <stdio.h>
#include <FLAME.h>
#include "Chol_right_looking.h"

#define PRINT_DATA
// Declaration of local prototypes.
void matrix_generate( FLA_Obj A, char type,char name);
double spendtimeinsecond();

int main( int argc, char *argv[] ) {
  int m_A,n_A;
  FLA_Obj A;
  double t1,t2;
  printf("Please enter the size of the matrix A : m_A and n_A \n");
  scanf("%d%d",&m_A,&n_A);
  FLA_Init();
  FLA_Obj_create(FLA_DOUBLE,m_A,n_A,0,0,&A);
  matrix_generate(A,'F','A');
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
#endif
  t1=spendtimeinsecond();
  Chol_right_looking_unb(A);
  t2=spendtimeinsecond();
#ifdef PRINT_DATA
  FLA_Obj_show( " A = [ ", A, "%le", " ];" );
#endif
  printf("The spent time for computing the LU factorizatio of A in second : %lf\n",t2-t1);
  FLA_Obj_free(&A);
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
