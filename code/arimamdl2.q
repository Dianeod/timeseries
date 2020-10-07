//The AR models were compared against statsmodel AutoReg

//Ordinary Least squares method is used to predict parameters

\d .tm

/*endog - variable to be predicted by model
/*exog - additional data to be included when predicting data
/*p - number of lags
/*d - number of differences to apply
/*q - number of residual errors
/*tr - include trend or not

// ARIMA mdls

// fit the ARIMA mdl
/. r the model parameters and data needed for future predictions
ARIMAfit:{[endog;exog;p;d;q;tr]
 // chk length
 if[not count[endog]~count[exog];i.err.len];
 // diff model
 I:i.diff[endog;d];
 // do stationary check
 if[not i.stat[I];i.err.stat];
 // run ARMA model
 ARMAfit[I;exog;p;q;tr],(enlist[`origd]!enlist d{deltas x}/neg[d] #endog)}

// predict future values
/*mdl - the fitted mdl
/*len - length of future values to predict
/ . r predicted values
ARIMApred:{[mdl;exog;len]
 // predict values
 preds:ARMApred[mdl;exog;len];
 // if differenced, revert back to original
 $[dval:count mdl`origd;dval _dval{count[x] msum x}/mdl[`origd],preds;preds]}

// ARMA mdls

// fit the autoregressive moving average model
/. r the model parameters and data needed for future predictions
ARMAfit:{[endog;exog;p;q;tr]
 // if ma=0, use AR model
 $[0~q;ARfit[endog;exog;p;tr],`q_coeff`resid!(();());
  i.ARMAmdl[endog;exog;p;q;tr]]}

// predict future values
/ . r predicted values
ARMApred:{[mdl;exog;len]
  // convert exog table to matrix
 if[98h~type exog;exog:"f"$i.mat[exog]];
 // if MA=0 then use ARpred
 $[count mdl`resid;
 [dict:enlist[`p]!enlist count[mdl`lags]-count[mdl`p_param];
  preds:{x>count[y[2]]}[len;]i.sngpred[mdl`params;exog;
     dict;;mdl`estresid;"arma"]/(mdl`lags;mdl`resid;());
 last preds];
 ARpred[mdl;exog;len]] 
 }

// PURE AR models

// fit the autoregressive mdl
/. r the model parameters and data needed for future predictions
ARfit:{[endog;exog;p;tr]
 // convert exog table to matrix
 if[98h~type exog;exog:"f"$i.mat[exog]];
 // estimate coefficients
 params:i.estparam[endog;exog;endog;`p`q!p,0;tr];
 // get lag values
 lagd:neg[p]#endog;
 // return dictionary with required info for predictions
 `params`tr_coeff`exog_coeff`p_param`lags!(params;tr#params;
     params[tr _til count[exog]];neg[p]#params;lagd)}

// predict future values
/ . r predicted values 
ARpred:{[mdl;exog;len]
 // convert exog to matrix
 if[98h~type exog;exog:"f"$i.mat[exog]];
 // predict future values
 preds:{x>count[y[2]]}[len;]i.sngpred[mdl`params;exog;enlist[`p]!enlist 0;;();"arma"]/(mdl`lags;();());
 last preds}


// PURE MA models

// Fit the moving average model 
/. r - returns the coefficients and residual errors	
MAfit:{[data;q]	
 // use innovation method to calculate params	
 innmat:i.innovations[data;q];	
 // get resid errors	
 resids:i.residuals[data;q;innmat];	
 `params`errors!(1_last innmat;resids)}	

// One step MA pred	
MApred:{[mdl]	
 // calculate 1 step ahead	
 mdl[`params] mmu mdl`errors}	


/ SARIMAX models

// Fit the SARIMAX model
/*s - a dictionary of seasonal components
/.r the model parameters and data needed for future predictions
SARIMAfit:{[endog;exog;p;d;q;tr;s]
 // chk length
 if[not count[endog]~count[exog];i.err.len];
 // diff model
 I:i.diff[endog;d];
 // Seasonality Diff
 if[s`D;I:s[`D]i.sdiff[s`m]/I];
 // do stationary check
 if[not i.stat[I];i.err.stat];
 // run ARMA model
 i.SARIMAmdl[I;exog;p;q;tr;s],`origd`origs!(d{deltas x}/neg[d] #endog;neg[s[`D]*s`m]#endog)}

// Predict future values
/.r predicted values
SARIMApred:{[mdl;exog;len]
 // predict values
 preds:i.SARIMAp[mdl;exog;len];
 // if seasonal differenced, revert to original
 if[count mdl[`origs];preds:i.revseasdf[mdl`origd;preds]];
 // if differenced, revert back to original
 $[dval:count[mdl`origd];dval _dval{count[x] msum x}/mdl[`origd],preds;preds]}

// AIC  Best model Prediction

/Return best parameters based on lowes AIC score
/*train - training data dictionary with keys `endog`exog
/*test - testing data dictionry  with keys `engod`exog
/*len - how many steps ahead to predict
/* dict - dictionary of parameters to try
/. r - params that scored the best
aicparam:{[train;test;len;dict]
 // get aic scores for each set of params
 scores:i.fitscore[train;test;len;]each flip dict;
 // return best value
 dict[;scores?bsc],enlist[`score]!enlist bsc:min scores}
