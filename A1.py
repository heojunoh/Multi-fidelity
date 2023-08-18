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
  return np.exp(medium(x))

 def medium(x):
  return (x+low(x)**2) + (np.cos(low(x))) 

 def low(x):
  return np.sin(2*np.pi*(x-0.1))


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

 end = time.time()
 print ("Training done in %f seconds" % (end - start))

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
c(0.01079722, 0.01080525, 0.01225556, 0.01253609, 0.01391071,
       0.01404068, 0.01496848, 0.01619847, 0.01628984, 0.01652756,
       0.01683593, 0.01728447, 0.01902683, 0.02035204, 0.02118697,
       0.0224731 , 0.0267823 , 0.02834109, 0.0288071 , 0.02903569,
       0.02980256, 0.02987355, 0.0305842 , 0.03145871, 0.03411151,
       0.03520274, 0.03530964, 0.0360884 , 0.0362846 , 0.03749156,
       0.04024009, 0.04042458, 0.04286345, 0.04348038, 0.04564838,
       0.04700627, 0.04736182, 0.04763847, 0.04792603, 0.04916653,
       0.04940488, 0.05045759, 0.05271451, 0.05562384, 0.05942964,
       0.06473726, 0.06935024, 0.06981348, 0.07063607, 0.0768579 ,
       0.07700743, 0.07761212, 0.0787409 , 0.08136152, 0.08773271,
       0.09748464, 0.09880646, 0.09907652, 0.09963727, 0.09989578,
       0.10034275, 0.10226903, 0.10497224, 0.10718671, 0.10747865,
       0.11320633, 0.1149581 , 0.11800641, 0.11905073, 0.12550582,
       0.13504052, 0.13515264, 0.13534362, 0.13902403, 0.1442205 ,
       0.15176257, 0.158932  , 0.16211878, 0.1638927 , 0.16396067,
       0.17097982, 0.17680585, 0.18305469, 0.20266778, 0.20440453,
       0.22506543, 0.2328039 , 0.24165253, 0.27264706, 0.27887932,
       0.28872882, 0.29440095, 0.30338803, 0.31811753, 0.38688249,
       0.47306996, 0.59381205, 0.78583299, 1.68426123, 2.29701339)

### mean CRPS result.nonlinear.meancrps ### The smaller, the better
c(0.0061601 , 0.00675008, 0.00699821, 0.00743983, 0.00850194,
       0.00864816, 0.00872859, 0.00899436, 0.00943959, 0.00954674,
       0.00966113, 0.01019109, 0.01194476, 0.01304387, 0.01445312,
       0.01460863, 0.01478206, 0.01574597, 0.01587721, 0.01600961,
       0.01695539, 0.01704405, 0.01737952, 0.01745323, 0.01758647,
       0.01779312, 0.01896038, 0.01911534, 0.01931503, 0.02030299,
       0.0212123 , 0.02273054, 0.02303463, 0.02327831, 0.02351862,
       0.02484613, 0.02620303, 0.02709762, 0.02821576, 0.02922276,
       0.0299305 , 0.03058252, 0.03084539, 0.03351106, 0.03363218,
       0.03386175, 0.03439767, 0.03719221, 0.03931102, 0.04078921,
       0.04096701, 0.04175526, 0.04338947, 0.04575712, 0.04844803,
       0.04870093, 0.04906037, 0.04911261, 0.05167955, 0.05268473,
       0.0537914 , 0.05401252, 0.05484361, 0.0549566 , 0.05551573,
       0.05765073, 0.06008867, 0.06027573, 0.06051757, 0.06333809,
       0.06724103, 0.06761747, 0.07354435, 0.07365155, 0.07459737,
       0.08270232, 0.08555094, 0.0867788 , 0.08971259, 0.09247778,
       0.09480289, 0.09819905, 0.10156655, 0.10686493, 0.11909083,
       0.12436003, 0.12591607, 0.1337602 , 0.13973111, 0.1401706 ,
       0.14055677, 0.18127244, 0.19265215, 0.19624693, 0.24537706,
       0.26311173, 0.28835533, 0.3967282 , 1.1177908 , 1.47679853)

### computation time result.park.comptime ### The smaller, the better
c(47.58649206161499, 43.9633629322052, 58.49593496322632, 54.456708908081055, 56.52161478996277, 
49.368226051330566, 55.05099678039551, 53.668368101119995, 53.84085416793823, 51.4822211265564, 
68.60556888580322, 47.54049468040466, 64.74662709236145, 52.18308091163635, 54.054673194885254, 
57.47967576980591, 52.729902029037476, 59.617583990097046, 79.26344513893127, 54.88638114929199, 
51.17211699485779, 43.48248600959778, 53.7282280921936, 57.232295989990234, 50.06611084938049, 
50.28274202346802, 54.38559603691101, 49.00454592704773, 64.1526210308075, 59.42630410194397, 
52.358598709106445, 53.64052414894104, 48.4995436668396, 52.42168188095093, 52.4921989440918, 
67.74826121330261, 53.58530807495117, 48.47695302963257, 70.34393882751465, 53.54081678390503, 
50.88684916496277, 51.01155090332031, 60.98721218109131, 59.63371229171753, 50.3912410736084, 
47.958996057510376, 57.67331409454346, 59.09324073791504, 58.864184856414795, 61.172019958496094, 
52.51533603668213, 60.630939245224, 58.41107892990112, 57.09679293632507, 73.43942928314209, 
73.84559488296509, 64.96402168273926, 61.80919289588928, 63.782763957977295, 50.60345387458801, 
79.61962389945984, 48.932740926742554, 69.91919684410095, 58.39850425720215, 48.368841886520386, 
62.49673581123352, 69.85646390914917, 56.79120206832886, 57.94804406166077, 69.74532198905945, 
54.605026960372925, 64.59307718276978, 65.60859799385071, 65.09089708328247, 80.25910091400146, 
51.68548512458801, 64.37776708602905, 62.56523394584656, 51.06946301460266, 52.336650133132935, 
63.450737714767456, 61.173306703567505, 62.934667110443115, 52.35030913352966, 62.97181010246277, 
60.45982480049133, 77.5310890674591, 56.71697783470154, 68.59799885749817, 60.86932611465454, 
54.839762926101685, 66.61005401611328, 60.643391132354736, 54.99477410316467, 49.10493493080139, 
51.39285087585449, 56.866870164871216, 47.940619230270386, 43.08665084838867, 50.482587814331055)


