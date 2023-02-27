# Cullen and Frey Space

Python implementation of Cullen and Frey space based on fitdistrplus package.

## Example

    import cullenFrey
    import numpy as np
    import scipy.stats as stats
    import matplotlib.pyplot as plt 
    
    fig = plt.figure(figsize=(10,10))
    data = stats.beta(2.31, 0.127).rvs(100000)
    cullenFrey.cullenFrey(data,plotData=True,resample=True)
    plt.show()
  

![alt text](https://github.com/rsautter/CullenFrey/blob/master/Figures/example.png)

## Requirements
  
  - Scipy
  - Matplotlib

## Install

    !pip install git+https://github.com/rsautter/CullenFrey/
    
 ## References:
 
 [1] Cullen, Alison C., H. Christopher Frey, and Christopher H. Frey. Probabilistic techniques in exposure assessment: a handbook for dealing with variability and uncertainty in models and inputs. Springer Science & Business, 1999.
 
 [2] Delignette-Muller, Marie Laure, and Christophe Dutang. "fitdistrplus: An R package for fitting distributions." Journal of statistical software 64 (2015): 1-34.
 
 
