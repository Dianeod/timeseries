# Example Time Series Notebooks

Time series is defined as a sequence of information collected sequentially in time. Unlike most supervised machine learning methods which takes every row in a data frame as an independent entry, time series forecasting takes into account trends observed in historical data in order to forecast future observations. This includes examining attributes such as seasonality components and lagged values among others. 

Throughout these notebooks various methods of time series forecasting are implemented on a range of datasets in order to compare implementation and results of models.


The following time series libraries are demonstrated in the notebooks:

1. [ARIMA](https://machinelearningmastery.com/arima-for-time-series-forecasting-with-python/): ARIMA is a auto regression model that uses past lagged values to forecast future values. Both a [statsmodel](https://www.statsmodels.org/stable/generated/statsmodels.tsa.arima_model.ARIMA.html) and kdb implementation are available.

1. [SARIMA](https://machinelearningmastery.com/sarima-for-time-series-forecasting-in-python/#:~:text=An%20extension%20to%20ARIMA%20that,data%20containing%20trends%20and%20seasonality.): A SARIMA model is an extension of the ARIMA model with the addition of seasonal components.


In each notebook, two datasets are tested:

- Daily Temp: Daily Minimum temperatures reached in Melbourne over a 10 year period
- Bike Rental : The number of bikes rented every hour recorded by TFL

These two datasets were chosen to give a well rounded view of how to build both a simple and more complex model. 
`Daily Temp` represents a simple time series dataframe which requires very little data preparation and has no additional data columns to be included in the training of the model. `Bike Rental` however, is a more complex time series which requires additional data preparation and exogenous variables to be added to the models. 

**Notes**
These notebooks implement vanilla models of each time series forcasters with no feature extraction performed on the dataset. For a more in detailed description of these models and how to achieve best results, please follow the links provided in the notebooks.  

## Requirements

- kdb+>=? v3.5 64-bit
- Python 3.x
- [embedPy](https://github.com/KxSystems/embedPy)
- [JupyterQ](https://github.com/KxSystems/jupyterq)
- [ML-Toolkit](https://github.com/KxSystems/ml) (v0.3.x)

## Dependencies

Install the Python dependencies with

pip
```bash
pip install -r requirements.txt
```

or with conda
```bash
conda install --file requirements.txt
```
**N.B.** Additionally [graphviz](http://www.graphviz.org/download/) must be installed on the system running the notebooks.
