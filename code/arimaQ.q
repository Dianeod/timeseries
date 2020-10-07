// This code is copied from the thelasians

// Moving Average
gamma:{[data;h](neg[h]_data)cov h _data}


innovations:{[data;q]
// data -- historical data (list of floating numbers) used to estimatethe MA(q) process
// q -- order of the MA(q) process 
 p:10*q; // until convergence 
 v:(p+1)#0f;
 alpha:(p+1;p+1)#0f; 
 v[0]:gamma[0;data]; 
 n:1;
 k:0;
 while[n<=p;
  k:0; 
  while[k<n;
   j:0; s:0f; 
   while[j<k;
    s:s+(alpha[k;k-j]*alpha[n;n-j]*v[j]);
    j:j+1; ];
   alpha[n;n-k]:(gamma[n-k;data] - s)%v[k]; k:k+1;
   s:0f; j:0;
   while[j<n;
   s:s+(alpha[n;n-j] xexp 2)*v[j];
  j:j+1; ];
  v[n]:gamma[0;data]-s; ];
  n:n+1; ];
   :alpha[p-1;1+ til q]; 
 };


// Thesasians
.quantQ.ts.HannanRissanen:{[data;p;q]
// data -- historical data (list of floating numbers) used to estimate the ARMA(p,q) process
// p -- order of the AR(p) process // q -- order of the MA(q) process 
 // Step 1 - Estimate AR(max(p+q)+1) 
 lags:1+(p|q); 
 arEst:.quantQ.ts.DurbinLevinson[data;lags]; 
  endog:data til[count[data]-lags]+\:til (lags);
  ests:endog mmu arEst;
  res:((lags)_data)-ests;
/ res:.quantQ.ts.residualsAR[arEst;data];
  res:(lags#0f),res;
// Step 2 - regress data[t] on (data[t-1],...,data[t-p],res[t-1],...,res [t-q])
// for t = lag+q, ..., n by OLS
 t:(lags+q) _ til count data;
 x:(reverse each .quantQ.ts.matrix[t;p;data]),'(reverse each .quantQ.ts.matrix[t;q;res]);
 :(inv[(flip x) mmu x] mmu flip x) mmu data[t];
 };

.quantQ.ts.DurbinLevinson:{[data;p]
// data -- historical data (list of floating numbers) used to estimate the AR(p) process.
// p -- order of the AR(p) model.
 p:p+1;
 phi:(p;p)#0f;
 v:p#0f;
  phi[1;1]:.quantQ.ts.gamma[1;data]%.quantQ.ts.gamma[0;data];
  v[0]:.quantQ.ts.gamma[0;data];
 n:1;
 s:0f;
 while[n<p;
 j:1; while[j<n;
 s:s+phi[n-1;j]*.quantQ.ts.gamma[n-j;data];
 j:j+1]; phi[n;n]:(.quantQ.ts.gamma[n;data]-s)%v[n-1]; s:0f;
 v[n]:v[n-1]*(1f-(phi[n;n] xexp 2));
 j:1;
 while[j<n;
 phi[n;j]:phi[n-1;j]-phi[n;n]*phi[n-1;n-j];
 j:j+1]; n:n+1];
 :1_phi[p-1;]; };

.quantQ.ts.gamma:{[data;h](neg[h]_data)cov h _data}

.quantQ.ts.matrix:{[t;q;d]
// t -- index
// q -- order of the MA(q) process
// d -- series
 :$[q>1; d[t-q],'.quantQ.ts.matrix[t;q-1;d];d[t-1]];
 };

.quantQ.ts.toeplitz:{[data;p]
// data -- historical data
// p -- order of the AR(p)
 t:abs (til p)-/:(til p); :(p;p)#{[x;data].quantQ.ts.r[data;x]}[;data] each raze t;
 };


.quantQ.ts.r:{[data;h]
// h -- order of autocorrelation
// data -- data :
 :.quantQ.ts.gamma[h;data]%.quantQ.ts.gamma[0;data];};

.quantQ.ts.forecastMA:{[data;q]
    // data --  historical data (list of floating) used to forecast the AR(lags) model
    // q -- degree of MA
    est:(q+2)#0f;
    v:(q+1)#0f;
    alpha:(q+1;q+1)#0f;
    v[0]:.quantQ.ts.gamma[0;data];
    n:1;
    k:0;
    while[n<=q;
        k:0;
        while[k<n;
            j:0; s:0f;
            while[j<k;
                s:s+(alpha[k;k-j]*alpha[n;n-j]*v[j]);
                j:j+1;
            ];
        alpha[n;n-k]:(.quantQ.ts.gamma[n-k;data] - s)%v[k];
        k:k+1;
        s:0f; j:0;
            while[j<n;
                s:s+(alpha[n;n-j] xexp 2)*v[j];
                j:j+1;
            ];
            v[n]:.quantQ.ts.gamma[0;data]-s;
        ];
        n:n+1;
    ];
    SFV;
    // one-step predictor of MA(q)
    n:1;
    while[n <= q; 
        j:1;s:0f;
        while[j <= n;
            s:s+(alpha[n;j]*(data[j]-est[j]));
            j:j+1];
        est[n+1]:s;
        n:n+1;
        s:0f
    ];
    :1_est;
 };

.quantQ.ts.innovations:{[data;q]
// data -- historical data (list of floating numbers) used to estimate
// q -- order of the MA(q) process 
 p:10*q; // until convergence 
 v:(p+1)#0f;
 alpha:(p+1;p+1)#0f; v[0]:.quantQ.ts.gamma[0;data]; n:1;
 k:0;
 while[n<=p;
 k:0; while[k<n;
 j:0; s:0f; while[j<k;
 s:s+(alpha[k;k-j]*alpha[n;n-j]*v[j]);
 j:j+1; ];
 alpha[n;n-k]:(.quantQ.ts.gamma[n-k;data] - s)%v[k]; k:k+1;
 s:0f; j:0;
 while[j<n;
 s:s+(alpha[n;n-j] xexp 2)*v[j];
 j:j+1; ];
 v[n]:.quantQ.ts.gamma[0;data]-s; ];
 n:n+1; ];
 :alpha[p-1;1+ til q]; };


.quantQ.ts.yuleWalker:{[data;p]
// data -- historical data (list of floating numbers)
// p -- order of the AR(p) model
 :inv[.quantQ.ts.toeplitz[data;p]] mmu (.quantQ.ts.r[data;] each 1+til p); };
