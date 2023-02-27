import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

def cullenFrey(data=np.random.rand(100000),plotData=False,resample=False,samplePercentage=0.85,npoints=10000,label='Sample'):
  '''
  Cullen-Frey plot
  --------------------------------------------------
  data - the dataset (it will be normalized)
  plotData - plot the data in the cullen-frey space
  resample - if true, then a set of kutosis and skewness will be generated from the dataset sample 
  samplePercentage - how many points should be sampled? (only used when resample=True)
  npoints - how many resamples should be performed? (only used when resample=True)
  --------------------------------------------------
  Rubens A. Sautter (2023)
  Based on fitdistrplus (https://cran.r-project.org/web/packages/fitdistrplus/index.html)
  '''
  dataN = data-np.average(data)
  dataN = dataN / np.std(dataN) 

  # measuring the data cullen-frey
  if resample:
    kurt, skew = [], []
    n = int(np.ceil(len(dataN)*samplePercentage))
    for i in range(npoints):
      sample = np.random.choice(dataN,n)
      kurt.append(stats.kurtosis(sample,fisher=False))
      skew.append(stats.skew(sample)**2)
  else:
    kurt, skew = stats.kurtosis(dataN,fisher=False), stats.skew(dataN)**2
  kurtmax = np.max(np.append(kurt,19))
  skewmax = np.max(np.append(skew,9))
  

  # estimating the beta distribution region according to fitdistrplus
  p = np.exp(-70)
  q = np.exp(np.linspace(-100,100,1000))
  
  skewBeta1 = (4*((q-p)**2)*(p+q+1))/(((p+q+2)**2)*p*q)
  kurtBeta1 = 3*(p+q+1)*(p*q*(p+q-6)+2*((p+q)**2))/(p*q*(p+q+2)*(p+q+3) )
  filter = ((skewBeta1<1e10) & (skewBeta1>1e-4))
  skewBeta1, kurtBeta1 = skewBeta1[filter],kurtBeta1[filter]
  fitBetaInf = np.polyfit(skewBeta1,kurtBeta1,deg=1)

  p = np.exp(70)
  skewBeta2 = (4*((q-p)**2)*(p+q+1))/(((p+q+2)**2)*p*q)
  kurtBeta2 = 3*(p+q+1)*(p*q*(p+q-6)+2*((p+q)**2))/(p*q*(p+q+2)*(p+q+3)) 
  filter = ((skewBeta2<1e10) & (skewBeta2>1e-4))
  skewBeta2, kurtBeta2 = skewBeta2[filter],kurtBeta2[filter]
  fitBetaSup = np.polyfit(skewBeta2,kurtBeta2,deg=1)

  #gamma
  gammaX = 4/q
  gammaY = 3 + 6/q

  #lognorm
  q = np.exp(np.linspace(-2,2,1000))
  q2 = np.exp(q**2)
  q2 = q2[q2<1e10]
  logNormX = ((q2+2)**2)*(q2-1)
  logNormY = ((q2**4)+2*(q2**3)+3*(q2**2)-3)

  # plot
  x = np.sort(skewBeta1) 
  plt.fill_between(x = x,y1=np.polyval(fitBetaInf, x),y2=np.polyval(fitBetaSup, x),facecolor = "lightgray",label='Beta')
  plt.plot(gammaX,gammaY,ls='--',color='k',label='Gamma')
  plt.plot(logNormX,logNormY,ls='dotted',color='k',label='LogNorm')
  plt.scatter([0],[3],s=200,label='Normal',marker='*',color='k')
  plt.scatter([0],[9/5],s=100,label='Uniform',marker='^',color='k')
  plt.scatter([4],[9],s=100,label='Exponential',marker='s',color='k')
  plt.scatter([0],[4.2],s=100,label='Logistic',marker='+',color='k')
  plt.xlabel(r"Skewness$^2$",fontsize=15)
  plt.ylabel(r"Kurtosis (Pearson)",fontsize=15)
  plt.title("Cullen and Frey",fontsize=17)
  
  if plotData:
  	plt.plot(skew, kurt,marker = '.', color='blue')
  	plt.plot(np.average(skew),np.average(kurt),marker = '.', color='red',label=label)

  plt.xlim(-0.25,skewmax+1)
  plt.ylim(0,kurtmax+1)
  plt.gca().invert_yaxis()
  plt.legend()


if __name__ == '__main__':
	fig = plt.figure(figsize=(10,10))
	cullenFrey(stats.beta(2.31, 0.127).rvs(100000),plotData=True)
	fig = plt.figure(figsize=(10,10))
	cullenFrey(stats.beta(2.31, 0.127).rvs(100000),plotData=True,resample=True)
	plt.show()
