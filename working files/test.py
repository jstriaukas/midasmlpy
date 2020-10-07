import numpy as np
import midasml2.midasml2 as ml
nrow=np.array([2000],dtype='int32')
eta=[]
y=[]
ldot=[]
for i in range(0,2000):
    eta.append(i**3)#some alternate numbers
    y.append((i**2)+i**2)# some alternate numbers
    ldot.append(0) #empty
eta=np.array(eta,dtype='float64')
y=np.array(y,dtype='float64')
ldot=np.array(ldot,dtype='float64')
ml.lin_grad_calc(nrow,eta,y,ldot)
print(ldot)