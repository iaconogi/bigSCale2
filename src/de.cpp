#include <Rcpp.h>
using namespace Rcpp;


//[[Rcpp::export]]
NumericVector DE_Rcpp(NumericMatrix expr_norm_1,NumericMatrix expr_norm_2, NumericMatrix log_scores, IntegerVector ind_A, NumericVector lib_size_1, NumericVector lib_size_2 ) {
  
  int gene_num = expr_norm_1.nrow(), sample_num_1 = expr_norm_1.ncol(),sample_num_2 =expr_norm_2.ncol();
  NumericVector results_DE(gene_num*2);
  int log_num = log_scores.nrow();
  int i,j;
  
  // fos the calculation of a single distance
  int A,B;
  int k,pos_a,pos_b;
  double mean_lib_size;
  //printf("Rcpp: gene_num=%d,sample_num_1=%d,sample_num_2=%d\n",gene_num,sample_num_1,sample_num_2);

  
  for (k=0; k<(gene_num*2); k++)
    {
    results_DE[k]=0;
    }


  for (i=0; i<sample_num_1; i++)//(i=0; i<1; i++) //(i=0; i<sample_num_1; i++)
    for (j=0; j<sample_num_2; j++) //(j=0; j<1; j++) //(j=0; j<sample_num_2; j++)
    {

      pos_a=i*gene_num;
      pos_b=j*gene_num;
      //printf("i=%d,j=%d,lib_size[i]=%f, lib_size[j]=%f\n",i,j,lib_size_1[i],lib_size_2[j]);
      mean_lib_size=(lib_size_1[i]+lib_size_2[j])/2*10;

      for (k=0; k<gene_num; k++)//gene_num
      {
        A=expr_norm_1[pos_a+k]*mean_lib_size;
        B=expr_norm_2[pos_b+k]*mean_lib_size;
        results_DE[k] = results_DE[k] + log_scores[(ind_A[B])*log_num+ind_A[A]];
        results_DE[k+gene_num] = results_DE[k+gene_num] + ((A>0) | (B>0));
        //merda=(int)((A+B)>0);
        //if ( (k==35859) & (results_DE[k+gene_num]>previous_value) )
        //  {
        //  previous_value=results_DE[k+gene_num];
        //  printf("k=%d) col_a=%d,col_b=%d, A=%d B=%d, ind_A(A)=%d, ind_A(B)=%d, log_score=%f,results_DE[k+gene_num]=%f\n",k, i , j, A, B, ind_A[A], ind_A[B],log_scores[(ind_A[B])*log_num+ind_A[A]],results_DE[k+gene_num]);
        //  }
      }
      //mexPrintf("col_a=%d,col_b=%d pos_A=%d pos_b=%d conta=%d log_score=%f\n",col_a,col_b,pos_a,pos_b,conta,log_score);

    }

return(results_DE);
  
}





// count=0;
// i=0;
// for (j=i+1; j<sample_num; j++) 
// {
//   score(expr_driving_norm,lib_size,ind_A,log_scores, i,j,gene_num,log_num,distances,count);
//   count=count+1;
// }
