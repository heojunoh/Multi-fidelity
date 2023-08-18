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

dim = 4

X1 = np.array([[.000550, 315.00, 310.00, 365.00],
            [.000552, 293.53, 318.63, 388.29],
            [.000566, 285.77, 266.71, 367.27],
            [.000578, 302.17, 358.13, 343.72],
            [.000580, 272.26, 211.71, 333.65],
            # [.000589, 278.16, 225.78, 351.83],
            [.000612, 280.83, 291.53, 394.72],
            [.000626, 284.89, 350.46, 352.29],
            [.000639, 270.45, 241.21, 341.81],
            [.000643, 276.17, 216.99, 371.60],
            [.000652, 298.04, 303.96, 361.58],
            # [.000700, 288.15, 300.00, 400.00],
            [.000751, 287.99, 326.02, 354.08],
            [.000780, 292.73, 267.84, 369.00],
            [.000800, 303.15, 250.00, 350.00],
            [.000814, 286.39, 339.92, 332.40],
            [.000842, 294.39, 203.45, 346.05],
            # [.000850, 301.31, 317.85, 341.00],
            [.000857, 282.12, 262.30, 350.10],
            [.000903, 284.25, 290.90, 364.99],
            [.000910, 248.87, 206.74, 398.00],
            [.000940, 271.32, 362.73, 400.00],
            [.000950, 280.00, 270.00, 330.00],
            # [.001000, 293.15, 202.40, 373.15],
            [.000669, 296.33, 343.16, 385.81],
            [.000670, 303.07, 321.41, 370.48],
            [.000683, 287.05, 227.31, 358.24],
            [.000689, 272.70, 260.91, 355.37],
            [.000694, 278.35, 212.79, 376.24],
            # [.000698, 277.52, 299.39, 338.40],
            [.000711, 292.26, 273.31, 392.54],
            [.000714, 283.08, 306.69, 344.34],
            [.000722, 276.53, 353.75, 374.41],
            [.000730, 285.51, 217.74, 383.92],
            [.000738, 295.01, 295.02, 347.22]#,
            # [.000741, 270.95, 275.19, 356.87]
            ])
            
X2 = np.array([[.000550, 315.00, 310.00, 365.00],
            [.000552, 293.53, 318.63, 388.29],
            [.000566, 285.77, 266.71, 367.27],
            [.000578, 302.17, 358.13, 343.72],
            [.000580, 272.26, 211.71, 333.65],
            # [.000589, 278.16, 225.78, 351.83],
            [.000612, 280.83, 291.53, 394.72],
            [.000626, 284.89, 350.46, 352.29],
            [.000639, 270.45, 241.21, 341.81],
            [.000643, 276.17, 216.99, 371.60],
            [.000652, 298.04, 303.96, 361.58],
            # [.000700, 288.15, 300.00, 400.00],
            [.000751, 287.99, 326.02, 354.08],
            [.000780, 292.73, 267.84, 369.00],
            [.000800, 303.15, 250.00, 350.00],
            [.000814, 286.39, 339.92, 332.40],
            [.000842, 294.39, 203.45, 346.05],
            # [.000850, 301.31, 317.85, 341.00],
            [.000857, 282.12, 262.30, 350.10],
            [.000903, 284.25, 290.90, 364.99],
            [.000910, 248.87, 206.74, 398.00],
            [.000940, 271.32, 362.73, 400.00],
            [.000950, 280.00, 270.00, 330.00]#,
            # [.001000, 293.15, 202.40, 373.15]
            ])
            
            
Y1 = np.array([7.02, 25.61, 21.23, 11.44, 15.03, #18.55,
        30.22, 18.13, 17.92, 24.20, 17.47, #30.90,
        18.17, 20.92, 13.08, 12.68, 13.75, #11.30, 
        18.25, 22.22, 36.56, 35.53, 13.54, #21.60,
        25.07, 18.93, 18.61, 21.31, 25.11, #16.02, 
        27.47, 16.43, 26.50, 25.88, 14.37#, #22.36
        ])[:,None]

Y2 = np.array([7.48, 23.54, 20.15, 10.17, 15.29, #18.39, 
        30.12, 18.17, 19.05, 24.96, 16.95, # 34.45,
        19.57, 21.97, 14.83, 14.36, 15.12, #11.92, 
        21.31, 25.37, 47.05, 42.93, 17.41, #22.89
        ])[:,None]

Xtest = np.array([[.000500, 293.15, 362.73, 393.15],
               [.000560, 277.01, 354.98, 374.00],
               [.000594, 279.54, 258.51, 360.13],
               [.000620, 275.00, 225.00, 340.00],
               [.000627, 287.60, 243.96, 382.54],
               [.000657, 294.24, 330.63, 375.53],
               [.000680, 313.28, 259.12, 350.00],
               [.000763, 292.82, 254.84, 373.38],
               [.000850, 270.00, 325.00, 385.00],
               [.000851, 273.71, 315.27, 381.14],
               [.000874, 282.50, 253.25, 396.36],
               [.000882, 299.22, 288.45, 385.07]])

Exact = np.array([25.82, 19.77, 20.52, 18.78, 24.68, 22.30, 4.55, 23.33, 32.85, 34.80, 36.11, 27.36])[:,None]
 
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
 mean
 
error = np.linalg.norm(Exact - mean)/np.linalg.norm(Exact)
# score = np.mean(-(Exact-mean)**2/var-np.log(var))
crps = np.mean(-np.sqrt(var)*(1/np.sqrt(np.pi)-2*stats.norm.pdf((Exact-mean)/np.sqrt(var))-(Exact-mean)/np.sqrt(var)*(2*stats.norm.cdf((Exact-mean)/np.sqrt(var))-1)))
ctime = (end - start)
# print( "N1 = %d, N2 = %d, sample = %d, error = %e" % (N1, N2[ii], jj+1, error))

l2error.append(error)
meanscore.append(score)
meancrps.append(crps)
comptime.append(ctime)


l2error

meanscore 

meancrps

comptime



### RMSE ###
2.842415

### mean score result.park.meanscore ### The larger, the better


### mean CRPS result.park.meancrps ### The smaller, the better
8.562537

### computation time result.park.comptime ### The smaller, the better




