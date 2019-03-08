#include <Rcpp.h>
using namespace Rcpp;


//[[Rcpp::export]]
NumericVector distances_icells(NumericMatrix expr_driving_norm, NumericMatrix log_scores, IntegerVector ind_A, NumericVector lib_size ) {
  
  int gene_num = expr_driving_norm.nrow(), sample_num = expr_driving_norm.ncol();
  NumericVector distances(sample_num-1);
  int log_num = log_scores.nrow();
  int i,j,count;
  
  // fos the calculation of a single distance
  int A,B;
  double log_score;
  int k,pos_a,pos_b;
  double mean_lib_size;
  //printf("gene_num=%d,sample_num=%d",gene_num,sample_num);
  
  count=0;
  
  i=0;
    for (j=i+1; j<sample_num; j++)
    {
      
      pos_a=i*gene_num;
      pos_b=j*gene_num;
      //printf("i=%d,j=%d,lib_size[col_a]=%f, lib_size[col_b]=%f\n",i,j,lib_size[i],lib_size[j]);
      mean_lib_size=(lib_size[i]+lib_size[j])/2*10;
      
      log_score=0;
      for (k=0; k<gene_num; k++) 
      {
        A=expr_driving_norm[pos_a+k]*mean_lib_size;
        B=expr_driving_norm[pos_b+k]*mean_lib_size;
        log_score = log_score + log_scores[(ind_A[B])*log_num+ind_A[A]];
        //printf("k=%d) col_a=%d,col_b=%d, A=%d B=%d, ind_A(A)=%d, ind_A(B)=%d, log_score=%f\n",k, i , j, A, B, ind_A[A], ind_A[B],log_scores[(ind_A[B])*log_num+ind_A[A]]);
      }
      //mexPrintf("col_a=%d,col_b=%d pos_A=%d pos_b=%d conta=%d log_score=%f\n",col_a,col_b,pos_a,pos_b,conta,log_score);
      distances[count]=log_score;
      count=count+1;
      
    }
    
    return(distances);
}




