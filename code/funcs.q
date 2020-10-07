// These are functions that are not currently being used within the ARIMA toolkit but might be useful for future

// Durbin Levinson function to calculate the ocefficients in a pure AR model with no trend for a univariate dataset
// Implementation can be found here https://www.stat.purdue.edu/~zhanghao/STAT520/handout/DurbinLevHandout.pdf
/. r - returns coefficients for lag values
durbin_lev:{[data;p]
 mat:(1+p;1+p)#0f;
 v:(1+p)#0f;
 mat[1;1]:.tm.i.acf[data;1];
 v[1]:var[data]*(1-xexp[mat[1;1];2]);
 reverse 1_last first(p-1){[data;d]
 mat:d[0];v:d[1];n:d[2];
 k:n+1;
 mat[k;k]:(.tm.i.lagcov[data;k]-sum mat[n;1+til n]mmu .tm.i.lagcov[data]each k-1+til n)%v[n];
 upd:{[data;n;mat;j]
  mat[n;j]-(mat[n+1;n+1]*mat[n;1+n-j])
 }[data;n;mat]each 1+til n;
 mat[k;1+til n]:upd;
 v[k]:v[n]*(1-xexp[mat[k;k];2]);
 (mat;v;n+1)}[data]/(mat;v;1)
 }

// Finding the params for pure MA models
/. r - returns innovations matrix to calculate past errors
i.innovations:{[data;q]
 v:(q+1)#0f;
 alpha:(q+1;q+1)#0f;
 v[0]:i.gamma[data;0]; 
 d:`alpha`v`n!(alpha;v;1);
 d:q{[data;d]
    d[`k]:0;
    d:d[`n]{[data;n;d]k:d[`k];
      s:sum d[`alpha][k;k-til k]*d[`alpha][n;n-til k]*d[`v][til k];
      d[`alpha;n;n-k]:(i.gamma[data;n-k]-s)%d[`v][k];
      l:sum(d[`alpha][n;n-til n] xexp 2)*d[`v][til n];
      d[`v;n]:v[0]-l;d[`k]+:1;d}[data;d`n]/d;
    d[`n]+:1;
    d}[data]/d;
 d`alpha
 }
// Find residual errors for pure ma model
/*mat - innovations matrix
/. r - returns residual errors from model
i.residuals:{[data;q;mat]
 qdata:neg[q]#data;
 est:(q)#0f;
 est:first q-1{[mat;qdata;en]
  est:en[0];n:en[1];
  j:1+til n;
  s:sum(mat[n;j]*(qdata[j-1]-est[j-1]));
  est[n]:s;
  (est;n+1)}[mat;qdata]/(est;1);
 qdata - est}

