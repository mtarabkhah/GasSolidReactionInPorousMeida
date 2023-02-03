from joblib import Parallel, delayed
import numpy as np
import time

def testfcn(i):
    time.sleep(0.0000005)
    return i*i,i

H = 200
W = 200
a = np.random.randint(1,100,[H,W])
#print(a)
count = 0
Sim = np.zeros_like(a)
tic = time.time()
for i in range(H):
    par_list = Parallel(n_jobs=4, prefer="threads")(delayed(testfcn)(a[i,j]) for j in range(W))
    count = count+sum([item[1] for item in par_list]) #Just for test
    Sim[i][:] = [item[0] for item in par_list]

print(time.time()-tic)
#print(Sim)
print(count)

count = 0
Sim = np.zeros_like(a)
tic = time.time()
for i in range(H):
    par_list = Parallel(n_jobs=1, prefer="threads")(delayed(testfcn)(a[i,j]) for j in range(W))
    count = count+sum([item[1] for item in par_list]) #Just for test
    Sim[i][:] = [item[0] for item in par_list]

print(time.time()-tic)
#print(Sim)
print(count)

count = 0
b = a.flatten()
Sim = np.zeros_like(b)
tic = time.time()
par_list = Parallel(n_jobs=4, prefer="threads")(delayed(testfcn)(b[j]) for j in range(H*W))
count = count+sum([item[1] for item in par_list]) #Just for test
#print(par_list)
Sim = [item[0] for item in par_list]
#for i in range(H):
#    Sim[i] = b[i*W:(i+1)*W]

print(time.time()-tic)
#print(Sim)
print(count)