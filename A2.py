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
  return -np.exp(medium(x))*np.sin(medium(x))

 def medium(x):
  return np.log(low(x)**2+np.sqrt(x+3))

 def low(x):
  return np.exp(-1.4*x) * np.cos(3.5*np.pi*x)

 def scale_range(x,ub,lb):
  Np = x.shape[0]
  dim = x.shape[1]
  for i in range(0,Np):
   for j in range(0,dim):
    tmp = ub[j] -lb[j]
    x[i][j] = tmp*x[i][j] + lb[j]
  return x

 # def rmse(pred, truth):
 #  pred = pred.flatten()
 #  truth = truth.flatten()
 #  return np.sqrt(np.mean((pred-truth)**2))


 ''' Create training set '''
 N1 = 15
 N2 = 10
 N3 = 5

 plot = 1
 save = 0

 dim = 1
 lb = np.array([0])
 ub = np.array([1])

 tmp = np.random.rand(1000,dim)
 Xtrain = scale_range(tmp,ub,lb)
 idx = np.random.permutation(1000)
 X1 = Xtrain[idx[0:N1], :]
 X2 = Xtrain[idx[0:N2], :]
 X3 = Xtrain[idx[0:N3], :]

 Y1 = low(X1)
 Y2 = medium(X2)
 Y3 = high(X3)

 tmp = np.random.rand(100,dim)
 Xtest = scale_range(tmp,ub,lb)

 Exact = high(Xtest)
 Medium = medium(Xtest)
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


 # Prepare for level 3: sample f_1 at X3
 nsamples = 100
 ntest = X3.shape[0]
 mu0, C0 = m1.predict(X3, full_cov=True)
 Z = np.random.multivariate_normal(mu0.flatten(),C0,nsamples)
 tmp_m = np.zeros((nsamples,ntest))
 tmp_v = np.zeros((nsamples,ntest))

 # push samples through f_2
 for i in range(0,nsamples):
  mu, v = m2.predict(np.hstack((X3, Z[1,:][:,None])))
  tmp_m[i,:] = mu.flatten()
  tmp_v[i,:] = v.flatten()

 # get mean and variance at X3
 mu2 = np.mean(tmp_m, axis = 0)
 v2 = np.mean(tmp_v, axis = 0) + np.var(tmp_m, axis = 0)
 mu2 = mu2[:,None]
 v3 = np.abs(v2[:,None])


 ''' Train level 3 '''
 XX = np.hstack((X3, mu2))

 k3 = GPy.kern.RBF(1, active_dims = [dim])*GPy.kern.RBF(dim, active_dims = active_dimensions, ARD = True) \
 + GPy.kern.RBF(dim, active_dims = active_dimensions, ARD = True)

 m3 = GPy.models.GPRegression(X=XX, Y=Y3, kernel=k3)

 m3[".*Gaussian_noise"] = m3.Y.var()*0.01
 m3[".*Gaussian_noise"].fix()

 m3.optimize(max_iters = 500)

 m3[".*Gaussian_noise"].unfix()
 m3[".*Gaussian_noise"].constrain_positive()

 m3.optimize_restarts(30, optimizer = "bfgs",  max_iters = 1000)

 # Compute posterior mean and variance for level 3 evaluated at the test points

 # sample f_1 at Xtest
 nsamples = 100
 ntest = Xtest.shape[0]
 mu0, C0 = m1.predict(Xtest, full_cov=True)
 Z = np.random.multivariate_normal(mu0.flatten(),C0,nsamples)

 # push samples through f_2 and f_3
 tmp_m = np.zeros((nsamples**2,ntest))
 tmp_v = np.zeros((nsamples**2,ntest))
 cnt = 0
 for i in range(0,nsamples):
  mu, C = m2.predict(np.hstack((Xtest, Z[i,:][:,None])), full_cov=True)
  try: 
    Q = np.random.multivariate_normal(mu.flatten(),C,nsamples)
    for j in range(0,nsamples):
     mu, v = m3.predict(np.hstack((Xtest, Q[j,:][:,None])))
     tmp_m[cnt,:] = mu.flatten()
     tmp_v[cnt,:] = v.flatten()
     cnt = cnt + 1
  except:
    pass

 # get f_2 posterior mean and variance at Xtest
 mu3 = np.mean(tmp_m, axis = 0)
 v3 = np.mean(tmp_v, axis = 0) + np.var(tmp_m, axis = 0)
 mu3 = mu3[:,None]
 # v3 = np.abs(v3[:,None])
 end = time.time()
 

 # Exact = Exact[:,None]
 error = np.sqrt(np.mean((mu3-Exact)**2))
 score = np.mean(-(Exact-mu3)**2/v3-np.log(v3))
 crps = np.mean(-np.sqrt(v3)*(1/np.sqrt(np.pi)-2*stats.norm.pdf((Exact-mu3)/np.sqrt(v3))-(Exact-mu3)/np.sqrt(v3)*(2*stats.norm.cdf((Exact-mu3)/np.sqrt(v3))-1)))
 ctime = (end - start)
 
 # print ("N1 = %d, N2 = %d, N3 = %d, error = %e, GP error = %e" % (N1, N2, N3, error, gpe))

 l2error.append(error)
 meanscore.append(score)
 meancrps.append(crps)
 comptime.append(ctime)


l2error
np.mean(l2error) 
np.sort(l2error) 

meanscore 
np.mean(meanscore) 
np.sort(meanscore) 

meancrps
np.mean(meancrps) 
np.sort(meancrps) 

comptime



### RMSE ###
c(0.00088608, 0.00111335, 0.00121508, 0.00125022, 0.00168015,
       0.00178216, 0.00212914, 0.00212981, 0.00225237, 0.00237729,
       0.00238091, 0.00279292, 0.00350591, 0.00364954, 0.00377493,
       0.00435657, 0.00480207, 0.00522593, 0.00524829, 0.00545775,
       0.00580225, 0.00612863, 0.00633811, 0.00655344, 0.00690516,
       0.00716252, 0.00735164, 0.00740154, 0.00760074, 0.00782773,
       0.00791001, 0.0081752 , 0.00821942, 0.0084137 , 0.00847452,
       0.00883437, 0.00883448, 0.00893649, 0.00949677, 0.00963008,
       0.00963531, 0.01016337, 0.010431  , 0.01174373, 0.01185404,
       0.01241278, 0.01264896, 0.01325998, 0.01338268, 0.01355914,
       0.01381407, 0.01409632, 0.01427213, 0.01445461, 0.01468199,
       0.01623944, 0.01642372, 0.01755051, 0.01798454, 0.01853255,
       0.01922675, 0.01932157, 0.01948749, 0.02104365, 0.02124347,
       0.02209883, 0.0234252 , 0.02442722, 0.02642747, 0.0273333 ,
       0.02825075, 0.0293954 , 0.02959933, 0.03031772, 0.03080215,
       0.03151073, 0.03299759, 0.03299915, 0.03432835, 0.03449503,
       0.04029537, 0.04086719, 0.04380565, 0.04470441, 0.0474623 ,
       0.04773728, 0.04802601, 0.05221179, 0.05398925, 0.05495056,
       0.05743543, 0.05884322, 0.06235784, 0.06300685, 0.08541811,
       0.08842227, 0.09524274, 0.11083501, 0.1231543 , 0.13667536)

### mean CRPS result.nonlinear.meancrps ### The smaller, the better
c(0.00055287, 0.00056893, 0.00058467, 0.00061104, 0.00076101,
       0.00082206, 0.00087378, 0.00098689, 0.00107653, 0.00127548,
       0.00145001, 0.00152505, 0.00159384, 0.00177264, 0.00203542,
       0.0024308 , 0.00245266, 0.00251307, 0.00255224, 0.00266618,
       0.00304524, 0.00317158, 0.00319169, 0.00323388, 0.00335426,
       0.0035017 , 0.00353358, 0.00362131, 0.0036398 , 0.00371769,
       0.0037292 , 0.00398891, 0.0040474 , 0.00421169, 0.00441127,
       0.00446011, 0.00452807, 0.00454843, 0.00454899, 0.00457419,
       0.00510512, 0.00517273, 0.00539537, 0.00569706, 0.00587295,
       0.00613868, 0.00615771, 0.00619656, 0.0062838 , 0.00633786,
       0.00636156, 0.00636392, 0.00657937, 0.00671933, 0.00680745,
       0.0069364 , 0.00700807, 0.00716673, 0.00755597, 0.00763527,
       0.00776563, 0.00803802, 0.00814938, 0.00822005, 0.00837698,
       0.00860951, 0.00868789, 0.00906137, 0.00907762, 0.00911888,
       0.01007754, 0.01097358, 0.01118291, 0.01139766, 0.01287789,
       0.01294118, 0.01313971, 0.01376261, 0.01410549, 0.01435616,
       0.01484915, 0.01508538, 0.0155802 , 0.01637103, 0.01734518,
       0.01800809, 0.01811862, 0.01828745, 0.01876778, 0.01954165,
       0.02236439, 0.02317984, 0.02412885, 0.02461406, 0.02762206,
       0.03129925, 0.033886  , 0.04575798, 0.04626372, 0.05836373)

### computation time result.park.comptime ### The smaller, the better
c(39.624094009399414, 46.95026159286499, 34.916425943374634, 47.21524000167847, 37.400604248046875, 
45.498799085617065, 33.87744998931885, 39.53625512123108, 45.21037697792053, 49.92695999145508, 
56.98320007324219, 44.459715127944946, 43.311001777648926, 47.756458044052124, 41.915776014328, 
44.746052980422974, 36.52470779418945, 42.316728830337524, 60.811194896698, 50.32084798812866, 
37.72214698791504, 37.757896900177, 56.45163297653198, 57.74765682220459, 52.71590280532837, 
43.58770418167114, 45.15324425697327, 43.24698090553284, 41.527403593063354, 47.31889009475708, 
48.130861043930054, 45.93967390060425, 38.87239599227905, 49.24992799758911, 48.55687594413757, 
43.193994998931885, 51.445817708969116, 41.52069973945618, 35.637669801712036, 42.80051398277283, 
35.81785297393799, 45.82896900177002, 46.663573026657104, 42.933902740478516, 45.37310218811035, 
51.94463586807251, 51.03771901130676, 46.58602023124695, 45.73778510093689, 48.49455118179321, 
50.580835819244385, 52.98732018470764, 50.787434816360474, 53.64008188247681, 56.803860902786255, 
48.620413303375244, 49.04798722267151, 44.372803926467896, 55.25683808326721, 59.887545108795166, 
39.63615584373474, 53.97800302505493, 66.13159012794495, 44.871362924575806, 42.267254114151, 
45.71021819114685, 42.43766903877258, 62.80208396911621, 55.91352295875549, 43.39822316169739, 
43.84833216667175, 43.77510976791382, 39.46296668052673, 46.66157412528992, 44.343899965286255, 
41.69769883155823, 44.23315191268921, 50.49132323265076, 46.413718938827515, 42.80970001220703, 
46.939780712127686, 40.08373284339905, 49.70857214927673, 64.62264800071716, 55.84169101715088, 
54.77261686325073, 37.368059158325195, 44.3578040599823, 49.684186935424805, 50.14914584159851, 
55.8731791973114, 43.82171392440796, 47.45457100868225, 36.53011894226074, 50.965299129486084, 
54.529321908950806, 52.29884910583496, 51.083418130874634, 43.7864670753479, 36.60940074920654)


