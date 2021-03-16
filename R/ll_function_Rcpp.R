### Rcpp Functions ###

#require(Rcpp)

### function to use in cont_normal3 (creating pvals matrix)
sourceCpp( code = '#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec create_pval(arma::vec x, arma::vec means, arma::vec means2, arma::vec sds, arma::uvec log, arma::vec sel, arma::uvec lowerT) {
  int n = x.size();
  arma::vec res(n);
  for(int i = 0; i < n; i++) {
    if(sel[i] = 1){
      res[i] = R::dnorm(x[i], means[i], sds[i], log[i]);
    } else{ 
            if(log[i]){
                res[i] = log(R::qnorm(x[i], means[i], sds[i], lowerT[i], false) -  R::qnorm(x[i], means2[i], sds[i], lowerT[i], false)) ;
            } else{
                res[i] = R::qnorm(x[i], means[i], sds[i], lowerT[i], false) -  R::qnorm(x[i], means2[i], sds[i], lowerT[i], false) ;
            }
    }
  }
  return res;
}')


#create_pval( x = d_df$Xb,
#             means = d_df$internal_lb,
#             means2 = d_df$internal_ub, 
#             sds = d_df$sigma_est,
#             log = rep(log.p,length(d_df$Xb)),
#             sel = d_df$internal_type,
#             lowerT = rep(TRUE,length(d_df$Xb)))






### USE THIS IN FOR LOOP to get Xb
#cppFunction( '
#NumericVector getXb(ExpressionVector exp, DataFrame df, NumericVector x){
#  int n = x.size();
#  NumericVector out(n);
#  out = exp.eval(df);
#  return(out);
#}
#
#')
#getXb(exp = parse(text = formulas[[1]]), df = d_df, x = d_df$d0)



###
### GIBT DATA FRAME (wie lapply) aus
#cppFunction( '
#DataFrame getXb(List formulas, DataFrame df){
#  int size = df.nrow();
#  int lenL = formulas.size();
#  List out(lenL);

#  for(int i = 0; i < lenL; i++) {
#    ExpressionVector exp = formulas[i];
#    out[i] = exp.eval(df);
#  }
#    out.attr("class") = "data.frame";
#  out.attr("row.names") = seq(1, size);
#  out.attr("names") = CharacterVector::create("main","SIGMA");
#return(out);
#}
#')
#head(getXb(formulas = lapply(formulas,function(x)parse(text = x)), df = d_df))