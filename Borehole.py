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
  x3 = x[:,2]
  x4 = x[:,3]
  x5 = x[:,4]
  x6 = x[:,5]
  x7 = x[:,6]
  x8 = x[:,7]
  return (2*np.pi*x3*(x4-x6))/(np.log(x2/x1)*(1+(2*x7*x3)/(np.log(x2/x1)*x1**2*x8)+(x3/x5)))

 def low(x):
  x1 = x[:,0]
  x2 = x[:,1]
  x3 = x[:,2]
  x4 = x[:,3]
  x5 = x[:,4]
  x6 = x[:,5]
  x7 = x[:,6]
  x8 = x[:,7]
  return (5*x3*(x4-x6))/(np.log(x2/x1)*(1.5+(2*x7*x3)/(np.log(x2/x1)*x1**2*x8)+(x3/x5)))

 def scale_range(x,ub,lb):
  Np = x.shape[0]
  dim = x.shape[1]
  for i in range(0,Np):
   for j in range(0,dim):
    tmp = ub[j] -lb[j]
    x[i][j] = tmp*x[i][j] + lb[j]
  return x

 ''' Create training set '''
 N1 = 80
 N2 = 40

 plot = 1
 save = 0

 dim = 8
 lb = np.array([0.05, 100, 63070, 990, 63.1, 700, 1120, 9855])
 ub = np.array([0.15, 50000, 115600, 1110, 116, 820, 1680, 12045])

 tmp = np.random.rand(1000,dim)
 Xtrain = scale_range(tmp,ub,lb)
 idx = np.random.permutation(1000)
 X1 = Xtrain[idx[0:N1], :]
 X2 = Xtrain[idx[0:N2], :]

 Y1 = low(X1)[:,None]
 Y2 = high(X2)[:,None]

 # nn = 40
 # x1 = np.linspace(lb[0], ub[0], 10)
 # x2 = np.linspace(lb[1], ub[1], 10)
 # x3 = np.linspace(lb[2], ub[2], 10)
 # x4 = np.linspace(lb[3], ub[3], 10)
 # x5 = np.linspace(lb[4], ub[4], 10)
 # x6 = np.linspace(lb[5], ub[5], 10)
 # x7 = np.linspace(lb[6], ub[6], 10)
 # x8 = np.linspace(lb[7], ub[7], 10)
 # X, Y = np.meshgrid(x1, x2)#, x3, x4, x5, x6, x7, x8)

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
 mu0, C0 = m1.predict(Xtest, full_cov=True)
 Z = np.random.multivariate_normal(mu0.flatten(),C0,nsamples)

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
c(85.64001612, 85.80437339, 86.73002939, 86.95017722, 86.98712327,
       87.12769297, 87.45519789, 87.48634341, 87.56742474, 87.78279148,
       87.85621479, 87.86694683, 87.9233014 , 87.96753205, 88.17682491,
       88.21380928, 88.26297619, 88.26526763, 88.26864269, 88.29456226,
       88.38521236, 88.39347625, 88.40050458, 88.43921543, 88.50386939,
       88.53470627, 88.62088881, 88.68516991, 88.69037   , 88.72536774,
       88.74113825, 88.74784762, 88.79661657, 88.88463065, 88.88683606,
       88.93382358, 88.96772845, 89.13823852, 89.24070832, 89.3073373 ,
       89.31514742, 89.32450153, 89.34169547, 89.3947988 , 89.45621231,
       89.47074324, 89.47469893, 89.53404981, 89.64991893, 89.7310376 ,
       89.77329303, 89.77857304, 89.78571346, 89.83626268, 89.92640951,
       90.01750269, 90.03437723, 90.07276807, 90.21247281, 90.23524965,
       90.27048604, 90.27910259, 90.3586534 , 90.3712481 , 90.44051929,
       90.46377021, 90.4640873 , 90.49262182, 90.5057352 , 90.51683443,
       90.58513476, 90.69702102, 90.71305564, 90.7579166 , 90.79728355,
       90.80506742, 90.97474317, 91.01390924, 91.0202576 , 91.08809224,
       91.24277077, 91.24576713, 91.30474663, 91.31146244, 91.39138725,
       91.47915909, 91.56630331, 91.67621847, 91.69788027, 91.73715048,
       91.75316665, 91.77531024, 91.83136602, 91.94769331, 92.14544891,
       92.51136494, 92.76423062, 93.40757245, 94.59379251, 95.43260832)

### mean CRPS result.park.meancrps ### The smaller, the better
c(49.29462471, 49.2980267 , 50.09820259, 50.10611749, 50.13905544,
       50.2642037 , 50.31809112, 50.45530659, 50.47234833, 50.56785915,
       50.63670308, 50.66977955, 50.71428898, 50.73148708, 50.82184035,
       50.88616105, 50.88804366, 50.89153971, 50.92732821, 50.93055938,
       50.95578799, 50.96262694, 50.9823663 , 51.0253229 , 51.05430568,
       51.07780624, 51.0780425 , 51.10396648, 51.1840877 , 51.20991604,
       51.28421588, 51.31527123, 51.33451128, 51.34107754, 51.37791952,
       51.37958484, 51.40765674, 51.44229555, 51.46384142, 51.59682893,
       51.61899525, 51.6746597 , 51.69808218, 51.7027772 , 51.70894838,
       51.73378427, 51.76525051, 51.78308338, 51.78315322, 51.80762076,
       51.82194404, 51.83817118, 51.84682794, 51.86668421, 51.88695548,
       51.89809504, 51.91261588, 51.98959469, 52.00891176, 52.02293739,
       52.05633049, 52.0896567 , 52.10920556, 52.12855853, 52.13015747,
       52.13356444, 52.2188267 , 52.28982923, 52.2953657 , 52.37094534,
       52.38177441, 52.39335955, 52.39580686, 52.47252351, 52.48761522,
       52.51523912, 52.52257224, 52.69457141, 52.69949172, 52.77416841,
       52.7806585 , 52.7959367 , 52.80766792, 52.89438096, 52.90452621,
       52.94153982, 52.94372474, 52.98543036, 52.99177638, 53.05195537,
       53.0609831 , 53.08978817, 53.15678613, 53.16424325, 53.38292602,
       53.88583814, 53.99809494, 54.41812143, 55.08558337, 56.51043403)

### computation time result.park.comptime ### The smaller, the better
c(13.21924090385437, 13.53073501586914, 13.371935844421387, 11.861296892166138, 11.814447164535522, 
12.085968017578125, 10.727544069290161, 11.397588014602661, 12.062791347503662, 12.70897102355957, 
11.754456043243408, 12.361602067947388, 12.184675216674805, 11.925798892974854, 12.289186000823975, 
12.273238182067871, 12.223232984542847, 10.867138862609863, 10.797057151794434, 10.527328968048096, 
8.187646865844727, 7.72758674621582, 7.996772050857544, 9.626137018203735, 12.089679956436157, 
12.500027179718018, 11.977385997772217, 11.9743070602417, 12.053840160369873, 12.146662950515747, 
12.169353008270264, 12.173588037490845, 12.008018016815186, 12.097518920898438, 9.47256088256836, 
8.371622800827026, 7.960578203201294, 8.297514915466309, 10.958580017089844, 10.838160276412964, 
11.078902959823608, 12.161823987960815, 12.194358587265015, 11.948587894439697, 12.152981042861938, 
12.213150978088379, 12.224562883377075, 11.742861270904541, 11.80446982383728, 10.729051113128662, 
10.69556188583374, 10.297934770584106, 10.584675073623657, 10.599033832550049, 10.759460926055908, 
10.844261169433594, 11.092674016952515, 12.070417881011963, 11.965091943740845, 12.085356950759888, 
11.874701023101807, 11.900409936904907, 11.35857629776001, 10.798681020736694, 10.756439208984375, 
10.980239868164062, 11.817746877670288, 10.73458480834961, 10.9587881565094, 10.77275276184082, 
10.927067041397095, 10.698134899139404, 11.161715984344482, 12.062072038650513, 12.030118942260742, 
12.011702060699463, 11.213928937911987, 10.675285816192627, 11.264014959335327, 12.211453914642334, 
12.011328935623169, 11.928683042526245, 12.068357229232788, 10.677475214004517, 10.824889898300171, 
10.739693880081177, 10.844388961791992, 10.546398878097534, 10.858871936798096, 9.812835216522217, 
10.622009992599487, 11.292539358139038, 12.025216102600098, 12.089556217193604, 11.908856868743896, 
12.000484943389893, 12.011229038238525, 11.76927399635315, 12.344007015228271, 12.410010814666748)



