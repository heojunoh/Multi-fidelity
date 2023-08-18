import GPy
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.mlab as ml
import matplotlib.patches as mpatches
import scipy.stats as stats

import time

l2error= []
meanscore= []
meancrps= []
comptime= []

for kk in range(1, 101): 
 np.random.seed(kk)

 def high(x):
  x1 = x[:,0]
  x2 = x[:,1]
  return (1 - np.exp(-1/(2*x2))) * (2300*x1**3 + 1900*x1**2 + 2092*x1 + 60) / (100*x1**3 + 500*x1**2 + 4*x1 +20)

 def low(x):
  x1 = x[:,0]
  x2 = x[:,1]
  xx11 = x1+0.05
  xx12 = x1-0.05
  xx21 = x2+0.05
  xx22 = np.maximum([0], [x2-0.05])
  return 1/4*(high(np.vstack((xx11,xx21)).T)+high(np.vstack((xx11,xx22)).T)+high(np.vstack((xx12,xx21)).T)+high(np.vstack((xx12,xx22)).T))

 def scale_range(x,ub,lb):
  Np = x.shape[0]
  dim = x.shape[1]
  for i in range(0,Np):
   for j in range(0,dim):
    tmp = ub[j] -lb[j]
    x[i][j] = tmp*x[i][j] + lb[j]
  return x

 ''' Define training and test points '''

 ''' Create training set '''
 N1 = 20
 N2 = 10

 plot = 1
 save = 0

 dim = 2
 lb = np.array([0, 0])
 ub = np.array([1, 1])

 tmp = np.random.rand(1000,dim)
 Xtrain = scale_range(tmp,ub,lb)
 idx = np.random.permutation(1000)
 X1 = Xtrain[idx[0:N1], :]
 X2 = Xtrain[idx[0:N2], :]

 Y1 = low(X1)[:,None]
 Y2 = high(X2)[:,None]

 nn = 40
 x1 = np.linspace(lb[0], ub[0], 50)
 x2 = np.linspace(lb[1], ub[1], 50)
 X, Y = np.meshgrid(x1, x2)

 tmp = np.random.rand(1000,dim)
 Xtest = scale_range(tmp,ub,lb)

 Exact = high(Xtest)
 Low = low(Xtest)
 
 active_dimensions = np.arange(0,dim)

 ''' Train level 1 '''
 start = time.time()
 k1 = GPy.kern.RBF(dim, ARD = True)
 m1 = GPy.models.GPRegression(X=X1, Y=Y1, kernel=k1)

 m1[".*Gaussian_noise"] = m1.Y.var()*0.01
 m1[".*Gaussian_noise"].fix()

 m1.optimize(max_iters = 500)

 m1[".*Gaussian_noise"].unfix()
 m1[".*Gaussian_noise"].constrain_positive()

 m1.optimize_restarts(30, optimizer = "bfgs",  max_iters = 1000)

 mu1, v1 = m1.predict(X2)


 ''' Train level 2 '''
 XX = np.hstack((X2, mu1))

 k2 = GPy.kern.RBF(1, active_dims = [dim])*GPy.kern.RBF(dim, active_dims = active_dimensions, ARD = True) \
 + GPy.kern.RBF(dim, active_dims = active_dimensions, ARD = True)

 m2 = GPy.models.GPRegression(X=XX, Y=Y2, kernel=k2)

 m2[".*Gaussian_noise"] = m2.Y.var()*0.01
 m2[".*Gaussian_noise"].fix()

 m2.optimize(max_iters = 500)

 m2[".*Gaussian_noise"].unfix()
 m2[".*Gaussian_noise"].constrain_positive()

 m2.optimize_restarts(30, optimizer = "bfgs",  max_iters = 1000)


 ''' Predict at test points '''
 # sample f_1 at xtest
 nsamples = 100
 ntest = Xtest.shape[0]
 mu1, C1 = m1.predict(Xtest, full_cov=True)
 Z = np.random.multivariate_normal(mu1.flatten(),C1,nsamples)

 # push samples through f_2
 tmp_m = np.zeros((nsamples,ntest))
 tmp_v = np.zeros((nsamples,ntest))
 for i in range(0,nsamples):
  mu, v = m2.predict(np.hstack((Xtest, Z[i,:][:,None])))
  tmp_m[i,:] = mu.flatten()
  tmp_v[i,:] = v.flatten()

 # get posterior mean and variance
 mean = np.mean(tmp_m, axis = 0)[:,None]
 var = np.mean(tmp_v, axis = 0)[:,None]+ np.var(tmp_m, axis = 0)[:,None]
 var = np.abs(var)
 end = time.time()
 
 Exact = Exact[:,None]
 
 error = np.sqrt(np.mean((mean-Exact)**2))
 score = np.mean(-(Exact-mean)**2/var-np.log(var))
 crps = np.mean(-np.sqrt(var)*(1/np.sqrt(np.pi)-2*stats.norm.pdf((Exact-mean)/np.sqrt(var))-(Exact-mean)/np.sqrt(var)*(2*stats.norm.cdf((Exact-mean)/np.sqrt(var))-1)))
 ctime = (end - start)
 # print( "N1 = %d, N2 = %d, sample = %d, error = %e" % (N1, N2[ii], jj+1, error))

 l2error.append(error)
 meanscore.append(score)
 meancrps.append(crps)
 comptime.append(ctime)


l2error
np.mean(l2error) # 0.7549168
np.sort(l2error) 

meanscore 
np.mean(meanscore) 
np.sort(meanscore) 

meancrps
np.mean(meancrps) # 0.33174714
np.sort(meancrps) 

comptime


### RMSE ###
c(0.29442103, 0.3186567 , 0.32653252, 0.35884138, 0.3612266 ,
       0.3844722 , 0.39507451, 0.39712256, 0.41062758, 0.41489027,
       0.41807915, 0.4235495 , 0.42748424, 0.42926974, 0.43253144,
       0.43507348, 0.43554247, 0.43916288, 0.45281416, 0.45451204,
       0.45868008, 0.45940901, 0.46398707, 0.47300148, 0.47714258,
       0.47728713, 0.48612575, 0.49195389, 0.49343751, 0.50160201,
       0.50723251, 0.50786866, 0.5086829 , 0.51449708, 0.52858122,
       0.53043482, 0.53224675, 0.53875863, 0.54796368, 0.54964504,
       0.56915951, 0.58196077, 0.59245869, 0.60428343, 0.609985  ,
       0.63588918, 0.64744829, 0.65831741, 0.66844321, 0.6746748 ,
       0.69451237, 0.69489135, 0.70941937, 0.71535954, 0.71961617,
       0.72219468, 0.74525468, 0.74751615, 0.75195187, 0.76328598,
       0.78247198, 0.78269446, 0.78646322, 0.80000871, 0.81419144,
       0.82354549, 0.82600036, 0.83651875, 0.8392843 , 0.84059   ,
       0.87707428, 0.89894201, 0.90958787, 0.93962372, 0.9621871 ,
       0.96368175, 1.00470257, 1.0155331 , 1.04454381, 1.05469583,
       1.06000599, 1.06246045, 1.07694135, 1.08235461, 1.08749761,
       1.10841457, 1.12305717, 1.14630055, 1.15815202, 1.15942962,
       1.19638412, 1.27647976, 1.40148636, 1.44381612, 1.47485264,
       1.48467365, 1.48758505, 1.58189863, 1.70460261, 2.00790359)

### mean CRPS result.park.meancrps ### The smaller, the better
c(0.14984793, 0.15505572, 0.15752919, 0.1698199 , 0.17900319,
       0.18479192, 0.18752176, 0.19301432, 0.19358862, 0.19858049,
       0.19965511, 0.20126463, 0.206678  , 0.2098936 , 0.21043123,
       0.21064509, 0.21116791, 0.21476699, 0.21666309, 0.21724425,
       0.21825302, 0.21901363, 0.22262649, 0.2237749 , 0.22923687,
       0.22969995, 0.23134981, 0.23409576, 0.23670036, 0.23862148,
       0.2401337 , 0.24342123, 0.24376436, 0.24652922, 0.24822088,
       0.25209849, 0.25473364, 0.25586932, 0.25760047, 0.25864618,
       0.26220888, 0.26239111, 0.26497999, 0.26526329, 0.26629332,
       0.27691318, 0.27783033, 0.28691924, 0.28803394, 0.28989036,
       0.29701956, 0.29913689, 0.29994336, 0.29995128, 0.30109846,
       0.30355179, 0.30439197, 0.31807817, 0.31906798, 0.33022416,
       0.33680343, 0.33803476, 0.34622141, 0.34853128, 0.35212875,
       0.36136308, 0.37145629, 0.3787886 , 0.38130153, 0.38332794,
       0.38349665, 0.38485145, 0.38525533, 0.39298484, 0.39986938,
       0.40135883, 0.40558233, 0.40605587, 0.40960297, 0.41125834,
       0.42218459, 0.43216168, 0.43272242, 0.43966233, 0.47053592,
       0.47574082, 0.48443594, 0.48444645, 0.49302779, 0.51296516,
       0.53790712, 0.55106439, 0.56819775, 0.59170665, 0.61832829,
       0.64467375, 0.67897264, 0.72814085, 0.76438428, 0.80044238)

### computation time result.park.comptime ### The smaller, the better
c(38.46457505226135, 45.29613900184631, 41.67369318008423, 13.852711915969849, 66.1725401878357, 
48.59166121482849, 42.97861981391907, 53.27743482589722, 42.89934515953064, 42.774128913879395, 
46.519551038742065, 36.85790181159973, 41.81056785583496, 40.66720390319824, 41.6387300491333, 
35.29694604873657, 45.70790982246399, 45.42031407356262, 56.07017803192139, 60.64559721946716, 
26.092360973358154, 53.423269271850586, 38.765719175338745, 29.527088165283203, 38.00006413459778, 
53.50367999076843, 45.0527720451355, 73.1283929347992, 29.322190046310425, 47.27030396461487, 
62.09379696846008, 24.999520301818848, 49.110761880874634, 18.03965473175049, 58.86049675941467, 
66.25800490379333, 37.40067386627197, 33.23986196517944, 44.91863489151001, 67.02728796005249, 
63.02077388763428, 33.2865469455719, 41.235963106155396, 50.644490003585815, 44.38124108314514, 
52.33662271499634, 25.506235122680664, 41.71459484100342, 71.1365180015564, 51.028106927871704, 
51.71041417121887, 36.98497986793518, 34.90028500556946, 37.8846960067749, 52.581767082214355, 
48.512667179107666, 36.47699785232544, 42.46441912651062, 49.92925405502319, 49.961191177368164, 
35.59922981262207, 44.32139492034912, 65.25653123855591, 51.57457995414734, 31.872881174087524, 
36.46962308883667, 44.602237939834595, 30.923030853271484, 30.89026689529419, 40.68335199356079, 
55.80264902114868, 29.661072254180908, 47.33366298675537, 34.87639594078064, 46.38042712211609, 
42.20359015464783, 44.61259698867798, 26.800437927246094, 66.36775994300842, 32.11775302886963, 
42.078288316726685, 27.806385040283203, 56.80334281921387, 44.93955183029175, 36.94848990440369, 
38.07759118080139, 31.344873905181885, 60.28279399871826, 77.21438503265381, 16.452447175979614, 
41.06525111198425, 42.43532395362854, 53.989501953125, 49.655259132385254, 49.13196682929993, 
34.81497597694397, 55.88296604156494, 46.69751811027527, 35.98944282531738, 32.41061568260193)



