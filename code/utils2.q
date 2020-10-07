// utility funcs for time series mdls

\d .tm

/*endog - variable to be predicted by model
/*exog - additional data table to be included when predicting data 
/*p - number of lags
/*d - number of differences to apply
/*q - number of residual errors
/*tr - include trend or not

// ARMA/AR models utils
// Fit ARMA mdl using Hannan-Rissanen 
/. r - model parameters used for predictions
i.ARMAmdl:{[endog;exog;p;q;tr]
 // convert exon table to matrix 
 if[98h~type exog;exog:"f"$i.mat[exog]];
 errs:i.esterrs[endog;exog;n:1+(p|q)];
 // Using the resid errorrs calculate coefficients for ARMA model
 coef:i.estparam[endog;exog;errs`errors;`p`q!p,q;tr];
 // return dictionary with required values for forecasting
 `params`tr_param`exog_param`p_param`q_param`lags`resid`estresid!
   (coef;coef[tr-1];coef[tr+til count exog[0]];p#neg[q+p]#coef;
        neg[q]#coef;(neg n)#endog;(neg q)#errs`errors;errs`params)
 }

// Estimate parameters of ARMA mdl
/*errors - errors from AR mdl
/. r - estimated params for ARMA mdl
i.estparam:{[endog;exog;errors;d;tr]
 // create lag matrix for endow values
 endogm:i.lagmat[endog;d`p];
 // create lag matrix for resid errors
 resid:i.lagmat[errors;d`q];
 // decide how many values to use from training
 m:neg min raze(count[endog]-d[`p`P]),count[errors]-d[`q`Q];
 // join exog, endow and resid values
 x:(m #exog),'(m #endogm),'m #resid;
 // add seasonality components
 if[not 0N~d[`P];x:x,'(m #flip[d[`P]xprev\:endog])];
 if[not 0N~d[`Q];x:x,'(m #flip[d[`Q]xprev\:errors])];
 // add trend line if specified
 if[tr;x:1f,'x];
 // values to predict
 y:m #endog;
 // use least squared error method to get coefficient values
 inv[fx mmu x]mmu(fx:flip x)mmu y
 }

// Predict single ARMA/AR/SARIMA value
/* d - list of p lag values, q resid errors and the previous predicted values
/*estresid - The model params used for estimating resin errors
/ . r - list of p lag values, q resid errors and list of predicted values
i.sngpred:{[params;exog;dict;d;estresid;typ]
 // check if trend is present
 pred:$[count[params]~count[m:get[".tm.i.",typ,"val"][exog;d;dict]];
 // no trend, so multiply coefficients by values
 params mmu m;
 // add trend value
 params[0]+((1_params) mmu m)];
 // if MA>0, estimate error values
 if[count[d[1]];
   estvals:exog[count d[2]],$["sar"~typ;#[neg[dict`est];];]d[0];
  d[1]:(1_d[1]),pred-mmu[estresid;estvals]];
 // append new lag values, and resid errors for next step calculations
 ((1_d[0]),pred;d[1];d[2],pred)}

i.sarval:{[exog;d;dict];
  exog[count d[2]],raze#[neg[dict`p];d[0]],#[neg[dict`q];d[1]],
   d[0][dict[`P]],d[1][dict[`Q]]}

i.armaval:{[exog;d;dict]
 exog[count d[2]],dict[`p] _raze d[0 1]}

// Estimate errors for Hannan Riessanan method
/*n - AR param to use
/.r - the residual errors and paramaters used to calculate them
i.esterrs:{[endog;exog;n]
 // build AR model to estimate resid errors
 estresid:ARfit[endog;exog;n;0b];
 endogm:i.lagmat[endog;n];
 pred:((neg[count endogm]#exog),'endogm) mmu estresid`params;
 errors:(n _endog)-pred;
 `params`errors!(estresid`params;errors)}

// SARIMAX model utils

// Revert season differenced data
/*origd - original data before being differenced
/*dfdata - differenced data
/.r - the data reverted back to its original format before differencing 
i.revseasdf:{[origd;dfdata]
 seasd:origd,dfdata;
 n:count origd;
 [n]_first{x[1]<y}[;count[seasd]]{[n;sdi]
 sd:sdi[0];i:sdi[1];
 sd[i]:sd[i-n]+sd[i];
 (sd;i+1)}[n]/(seasd;n)}

i.SARIMAmdl:{[endog;exog;p;q;tr;s]
 if[98h~type exog;exog:"f"$.tm.i.mat[exog]];
 errs:i.esterrs[endog;exog;n:1+(p|q)];
 dict:`p`q`P`Q!p,q,(1+til each s[`P`Q])*s[`m];
 coef:i.estparam[endog;exog;errs`errors;dict;tr];
 ns:sum s`P`Q;
 coefn:neg[ns]_coef;coefs:neg[ns]#coef;
 i.SARIMAcols!(coef;coefn[tr-1];coefn[tr+til count exog[0]];p#neg[q+p]#coefn;
   neg[q]#coefn;count[dict`P] #coefs;neg count[dict`Q] #coefs;
   (neg n|max[dict`P])#endog;(neg q|max dict`Q)#errs`errors;errs`params;2_dict)}

// Seasonal differencing
/*m - order of the seasonal component
/*d - data to apply differencing on
/.r - seasonal diffenced data
i.sdiff:{[m;d][m]_ d-(m xprev d)}

i.SARIMAp:{[mdl;exog;len]  // convert exog table to matrix
 if[98h~type exog;exog:"f"$.tm.i.mat[exog]];
 // if MA=0 then use ARpred
 $[count mdl`resid;
 [cnt:count each mdl;
 dict:`q`p`est!cnt[`q_param`p_param],cnt[`estresid]-cnt[`exog_param];
 dict:dict,enlist each mdl[`S]-min each mdl[`S];
 preds:{x>count[y[2]]}[len;]i.sngpred[mdl`params;exog;
     dict;;mdl`estresid;"sar"]/(mdl`lags;
  mdl`resid;());last preds];
 ARpred[mdl;exog;len]]}

// AIC utils

// Fit the model with the params and score it using AIC
/*train - training data
/*test - testing data
/*param - list of parameters used (q,p,trend)
/. r - the AIC score to the corresponding params
i.fitscore:{[train;test;len;param]
 // create model using specified params
 mdl:ARIMAfit[train`endog;train`exog;param`p;param`d;param`q;param`tr];
 // predict values using model coeffs
 preds:ARIMApred[mdl;test`exog;len];
 // calculate aic score using predictions and true values
 i.aic[len# test`endog;preds;param]}

// Apply AIC scoring
/*true - true values
/*pred - predicted values
/. r - returns AIC score
i.aic:{[true;pred;params]
 // get residual sum of squares
 rss:sqr wsum sqr:true-pred;
 // get aic score
 sc:(2*k:sum params)+n*log(rss%n:count[pred]);
 //If k<40, use altered aic score
 $[k<40;sc+(2*k*(k-1))%n-k-1;sc]}

// Auto correlation function
/*h - order of autocorrelation
i.acf:{[data;h]
 i.gamma[h;data]%i.gamma[0;data]}


// General Utils

// Create a lagmatrix
/*data - data to create matrix 
/*lag - number of lags
i.lagmat:{[data;lag]
 n:count data;
 data til[n-lag]+\:til lag
 }

// Create matrix from table
i.mat:{flip value flip x}

/ Indicates if the data is stationary
/*r - returns boolean
i.stat:{[data]
 adfull:.ml.fresh.i.adfuller[data]`;
 $[adfull[1]<0.05;1b;0b]}

// Differencing to produce stationary time series
/d* - order of differencing
/ . r - list of differences in data
i.diff:{[data;d]
 d _d{deltas x}/data}

// AutoCorrelation function
/h* - order of autocorrelation
/. r -  autocorrelation of order h 
i.gamma:{[h;data]
 (neg[h]_data) cov h _data}

// Columns of ARMA mdl
i.ARMAcols:`params`tr_param`exog_param`p_param`q_param`lags`resid`estresid
// Columns of SARIMA mdl
i.SARIMAcols:`params`tr_param`exog_param`p_param`q_param`P_param`Q_param`lags`resid`estresid`S 

// Error calls

i.err.steps:{'`$"Exog length not long enough"} 
i.err.stat:{'`$"Time series not stationary, try another value of d"}
i.err.len:{'`$"Endog and Exog not the same length"}
