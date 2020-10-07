\l ml/ml.q
.ml.loadfile`:init.q

plt:.p.import[`matplotlib.pyplot]

printfacts:{[qsar;pysar;typ]
 -1"\nTotal ",typ," models for each q and py implementations: ",string count qsar;
 -1"\nTotal time taken for q ",typ," models is ",string sum qsar`overall_tm;
 -1"\nTotal error for q ",typ," models is ",string sum qsar`score;
 -1"\nTotal time taken for py ",typ," models is ",string sum pysar`overall_tm;
 -1"\nTotal error for py ",typ," models is ",string sum pysar`score;}

perc:{[qtab;pytab;n]count[where qtab[`score]<(pytab[`score]*1+n)]%count qtab}

perc_print:{[qdata;pydata;typ]
 -1"\nPercentage of times q models recieved lower error scores than py models ",typ,": ",string count[where qdata[`score]<pydata`score]%count qdata;}

compare:{[qscore;pyscore;x;y]
 plt[`:plot][qscore;`label pykw "q_ARIMA"];
 plt[`:plot][pyscore;`label pykw "py_ARIMA"];
 plt[`:legend][];
 plt[`:xlabel][x];
 plt[`:ylabel][y];
 plt[`:title][y," vs ",x];
 plt[`:show][];}

tmplot:{[tab;comp;descr;typ] 
 plt[`:plot][.ml.minmaxscaler[tab`overall_tm];`label pykw "Time"];
 plt[`:plot][.ml.minmaxscaler[comp];`label pykw descr];
 plt[`:legend][];
 plt[`:xlabel]["Run"];
 plt[`:ylabel]["Scale"];
 plt[`:title][typ," Scale vs run"];
 plt[`:show][];}

scorevar:{[qdata;pydata;vars]
 p:perc[qdata;pydata]each vars;
 plt[`:bar]["f"$vars;p;`width pykw 0.03];
 plt[`:xlabel][`variance];
 plt[`:ylabel][`percentage];
 plt[`:title]["Percentage vs Variance"];
 plt[`:show][];}

time_ratio:{[qdata;pydata;typ]
 -1"On average, the Python vs q time ratio for the ",typ," model was: ",string avg pydata[`overall_tm]%qdata[`overall_tm];
 plt[`:plot][pydata[`overall_tm]%qdata[`overall_tm]];
 plt[`:xlabel]["py/q"];
 plt[`:ylabel]["run"];
 plt[`:title][typ," py/q vs run"];
 plt[`:show][];}

scoring:{[qdata;pydata;depen;descr;typ]
 plt[`:plot][.ml.minmaxscaler qdata[`score]-pydata`score;`label pykw "Python-q error score"];
 plt[`:plot][.ml.minmaxscaler depen;`label pykw descr]; 
 plt[`:legend][];
 plt[`:xlabel]["Run"];
 plt[`:ylabel]["Scale"];
 plt[`:title][typ," Error Differencing Scale vs run"];
 plt[`:show][];
 }
