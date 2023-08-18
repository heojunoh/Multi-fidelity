import GPy
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.mlab as ml
import matplotlib.patches as mpatches
import pandas as pd
import scipy.stats as stats

import time

l2error= []
meanscore= []
meancrps= []
comptime= []

dim=2
active_dimensions = np.arange(0,dim)
 
Xtest = np.array(pd.read_table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/Xtest.txt", sep=","))
Exact = np.array(pd.read_table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/ytest.txt", sep=","))

#########################################

X1 = np.array(pd.read_table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_X.txt", sep=","))
Y1 = np.array(pd.read_table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep=",", header=None)[3])[:,None]
Y1


X2 = np.array(pd.read_table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_X.txt", sep=","))
Y2 = np.array(pd.read_table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep=",", header=None)[3])[:,None]
Y2


X3 = np.array(pd.read_table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_X.txt", sep=","))
Y3 = np.array(pd.read_table("/Users/junoh/Downloads/StackingDesign-Reproducibility/Rmatlab_files/generate_text/temp_to_r.txt", sep=",", header=None)[3])[:,None]
Y3


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
 mu, v = m2.predict(np.hstack((X3, Z[i,:][:,None])))
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
 Q = np.random.multivariate_normal(mu.flatten(),C,nsamples)
 for j in range(0,nsamples):
  mu, v = m3.predict(np.hstack((Xtest, Q[j,:][:,None])))
  tmp_m[cnt,:] = mu.flatten()
  tmp_v[cnt,:] = v.flatten()
  cnt = cnt + 1


 # get f_2 posterior mean and variance at Xtest
mu3 = np.mean(tmp_m, axis = 0)
v3 = np.mean(tmp_v, axis = 0) + np.var(tmp_m, axis = 0)
mu3 = mu3[:,None]
v3 = np.abs(v3[:,None])
end = time.time()

 # Exact = Exact[:,None]
error = np.linalg.norm(Exact - mu3)/np.linalg.norm(Exact)
score = np.mean(-(Exact-mu3)**2/v3-np.log(v3))
crps = np.mean(-np.sqrt(v3)*(1/np.sqrt(np.pi)-2*stats.norm.pdf((Exact-mu3)/np.sqrt(v3))-(Exact-mu3)/np.sqrt(v3)*(2*stats.norm.cdf((Exact-mu3)/np.sqrt(v3))-1)))
ctime = (end - start)

l2error.append(error)
meanscore.append(score)
meancrps.append(crps)
comptime.append(ctime)





l2error

meanscore 

meancrps

comptime



np.mean(l2error) 
np.sort(l2error) 

np.mean(meanscore) 
np.sort(meanscore) 

np.mean(meancrps) 
np.sort(meancrps) 


### L2 error ###
c(0.08296562284056933, 0.26944568322745277, 0.08139932349505961, 0.06591543174151877, 0.07539813412693266, 
0.054232110710972156, 0.08305997653664889, 0.06354733784199235, 0.05857402773073935, 0.1358282065717591, 
0.05184968080311192, 0.060033594336085724, 0.08171029623631187, 0.17766831135547617, 0.05291498799731878, 
0.06584816603092461, 0.03009986302911407, 0.04810813255207206, 0.06714902793781717, 0.10171316513508435, 
0.06671993366295795, 0.050836904522805114, 0.12107033428616551, 0.06044577038368387, 0.0677842654704768, 
0.0712270546508961, 0.08555776795597886, 0.06678370359631516, 0.07219886395976191, 0.055628482275447574, 
0.08072438995677511, 0.055625115918103755, 0.07098293549593825, 0.06434466047336887, 0.10888285474073596, 
0.07422725981139436, 0.0633772920358081, 0.05796373984467013, 0.05492804850918673, 0.08209933537618808,
0.05932549535668908, 0.05262276438485215, 0.05718896275464744, 0.06434278236699577, 0.0738492087860961, 
0.06421319849590207, 0.11957146209403301, 0.12055272961881312, 0.07923707571466478, 0.0432476866983063)

### Mean score ###
c(-1.226594881411466, -5.918071411935286, -1.9558134642277685, -1.1085329528957575, -1.488964865587706,
-1.0080137923033425, -1.4315922877597598, -1.7087502834053412, -1.2037667689509974, -1.972561304609182, 
-1.4881905697693818, -1.6943623807970218, -14.979284910707147, -2.1177263505148507, -2.569032939900389,
-1.4842397062145218, -0.8799405048802597, -1.090908447988584, -1.4026922180858028, -7.702982712178904,
-1.6542554315350002, -1.2317690422106842, -1.774163331111413, -1.3164878700815088, -2.016471303590664,
-1.7111590264163186, -1.7097943540423484, -1.815623748849739, -1.8387515849737113, -1.3239455410075098, 
-2.9060598636177826,  -1.409988196868963, -1.7499975034338784, -1.2701768252674852,-4.132453999219285, 
-1.6906501203035058, -1.3586512654534149, -1.4651167782278274, -1.1106552422294056, -2.0169536276943574, 
-1.4044022126524454, -1.1321504460605327, -1.46435129580732, -1.7986160831840154, -1.0989744158653132, 
-1.501738937874928, -1.5566837480191227, -2.2544837541746277, -1.6985180936658764, -1.6379162095425286)

### Mean CRPS ###
c(0.8665093009231732, 3.996773266725916, 0.9802456550518821, 0.6696942386202817, 0.8310129570466523, 
0.6079964990070106, 0.8242488435136438, 0.7511318502757607, 0.6535410076500228, 1.103256421420058, 
0.6275321686864954, 0.7299826779332969, 1.0389419776713766, 2.146308417291548, 0.7090765223645336, 
0.7745347977930385, 0.4509185498132722, 0.5678350412189771, 0.7818611893871226, 1.0411421605781788,
0.7248813654899693, 0.6415190011288245, 1.316842718739142, 0.6781357965493838, 0.8546743101841922, 
0.8404619653981069, 0.8762054608835957, 0.855108086986757, 0.8582863606674513, 0.6732108401997121, 
0.967623455683443, 0.6701625579856461, 0.8315109786094997, 0.6687164971636793, 1.1600744842502486, 
0.8909211338068274, 0.7130320995639049, 0.6955417205478769, 0.6198301378860349, 0.9090131335823368,
0.6687303161557062, 0.6335639385687759, 0.6685754180500859, 0.7814467739024099, 0.7310736180116385, 
0.7622633520528582, 0.9875835837418133, 1.1314801482714438, 0.8567294136300516, 0.6545641497112342)

### computation time ###
c(59.946934938430786, 52.38716220855713, 45.430317878723145, 53.316006898880005, 48.21130299568176, 
56.10589289665222, 44.731395959854126, 47.18866419792175, 60.56380605697632, 57.43849778175354, 
48.32838797569275, 53.80115723609924, 45.838510274887085, 46.60700488090515, 45.262287855148315, 
40.9620361328125, 34.07575297355652, 57.633909940719604, 43.846577167510986, 44.71290373802185, 
57.686196088790894, 50.56412315368652, 62.51956605911255, 51.00797200202942, 53.17219305038452, 
51.76465892791748, 58.04509210586548, 109.99326395988464, 47.08927607536316, 47.45210385322571, 
41.366633892059326, 49.39511013031006, 57.91219902038574, 45.636048316955566, 39.643718242645264, 
50.07098174095154, 46.20361399650574, 63.02182722091675, 60.79853701591492, 44.08073091506958,
48.222208976745605, 53.435014963150024, 41.8203980922699, 50.117005825042725, 52.72865009307861, 
52.23375916481018, 33.32086992263794, 55.11390280723572, 46.73026132583618, 51.14283609390259)








