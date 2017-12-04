
import numpy as np
import matplotlib
matplotlib.use('Agg')
from scipy import *
from matplotlib import pyplot as plt
import math
import cmath
import sys
#import random
from pylab import *
from scipy import integrate
from scipy.integrate import ode
import matplotlib.pyplot as pyplot
#from adj import adj_mat1
from scipy import linalg as LA
from numpy.linalg import *
#import stulanN3

import matplotlib.cm as cm
import collections
from random import randint
from scipy.linalg import eig
from scipy.sparse.linalg.eigen.arpack import eigsh
from matplotlib import ticker


try:

     beta = float(sys.argv[1])
     N = int(sys.argv[2])

except:
     print "Usage:",sys.argv[1], "Not enough entries"; sys.exit(1)

#print n, b, inp, beta
sys.stdout.flush()


def things():
     return beta

def vals():
     beta=things()
     return beta
beta=vals()
#adj_mat = np.ones([2*int(N)-1,2*int(N)-1])
adj_mat = np.loadtxt('adj_mat.txt')
'''for i in range(int(N)-1):
     for j in range(int(N)):
         adj_mat[i][j]=adj_mat[i+N][j+N-1]=adj_mat1[i][j]
for i in range(int(N)-1):
    adj_mat[i][N-1]=adj_mat[i+N][N-1]=adj_mat1[i][N-1]
adj_mat[N-1]=adj_mat[N-1]/2'''
#np.savetxt('adj_mat.txt',adj_mat)
#N, adj_mat = stulanN3.value()

#mu=eiges()
#mu.sort()
#mu=[1,cmath.exp((1j*2*math.pi)/5),cmath.exp((1j*4*math.pi)/5),cmath.exp((1j*6*math.pi)/5),cmath.exp((1j*8*math.pi)/5)] #unidire one neighbor coupling
#mu=[2.0/2.0,(cmath.exp((1j*2*math.pi)/5)+cmath.exp((1j*8*math.pi)/5))/2,(cmath.exp((1j*4*math.pi)/5)+cmath.exp((1j*16*math.pi)/5))/2,(cmath.exp((1j*6*math.pi)/5)+cmath.exp((1j*24*math.pi)/5))/2,(cmath.exp((1j*8*math.pi)/5)+cmath.exp((1j*32*math.pi)/5))/2] #this is for R=1 bidriectional
#mu = [-2.0/3.0, -1.0/3.0, 0.0, 0.0, 1.0]
#print len(mu)
#zz=[cmath.exp((1j*6*math.pi)/5)]


mu,eigen_vec=LA.eig(adj_mat)
print"eigenvalues: ", mu
print "eigenvectors: " ,eigen_vec

np.savetxt('eigenvalue.txt',mu)
np.savetxt('eigenvector.txt',eigen_vec)
#this is for complete synchronization solution

for eta in mu:
    for lam in lambd:
        for sig in sigma:
            count=0
            for i in range(len(mu)):
                #x= sig*(eta-mu[i])*(2*lam*math.cos(beta) + sig*(eta*math.cos(2*beta) + 2*eta - mu[i])) > 0
                #print x
                if mu[i]== eta:#
                     
                    continue
                elif  (-2*(lam + sig*math.cos(beta)*(2*eta - mu[i])) < 0 and sig*(eta-mu[i])*(2*lam*math.cos(beta) + sig*(eta*math.cos(2*beta) + 2*eta - mu[i])) > 0 and lam > -eta*sig*math.cos(beta)) :
                    count+=1
               # print count    
                
            if count == (len(mu)-1):
                lambdasol.append(lam)
                sigmasol.append(sig)


# now for steady state trivial solution


for lam in lambd:
    for sig in sigma:
        count=0
        for i in range(len(mu)):
                #x= sig*(eta-mu[i])*(2*lam*math.cos(beta) + sig*(eta*math.cos(2*beta) + 2*eta - mu[i])) > 0
                #print x
            
            if  (lam<-mu[i]*sig*math.cos(beta)) :
                count+=1
               # print count    
                
            if count == (len(mu)):
                lambdasol1.append(lam)
                sigmasol1.append(sig)


#now for solutions containing amplitude death


I=np.eye(int(N))
Jb=np.array([[3,0],[0,1]])
R=np.array([[math.cos(beta),-math.sin(beta)],[math.sin(beta),math.cos(beta)]])


#eigzz=np.array([[1.0/math.sqrt(5),cmath.exp((1j*6*math.pi)/5)/math.sqrt(5),cmath.exp((1j*12*math.pi)/5)/math.sqrt(5),cmath.exp((1j*18*math.pi)/5)/math.sqrt(5),cmath.exp((1j*24*math.pi)/5)/math.sqrt(5)], [1,1,1,1,1]])
n=[]
hi=[]
V1=[]
V=[]
for i in range(int(N)):
    V1.append(eigen_vec[i]*eigen_vec[i])
#V1=np.array([1,1,0,1,1])
#print V1
for i in range(int(N)):
    V.append(np.diag(V1[i]))

C=np.zeros((int(N),int(N),2,2))
m=np.zeros((int(N),int(N),2,2))
n=np.zeros((int(N),int(N),2,2))
o=np.zeros((int(N),int(N),2,2))

lambdasol2=[]
sigmasol2=[]
count=0
#zz=[eigen_vec[3]]
for eta, ve in zip(mu,V):
#for eta,ve in zip(zz,eigzz): #ad,V[1:3]): 
#considering only eigenv to see it's stability
    for lam in lambd:
        for sig in sigma:
           Ja=np.array([[lam,eta*sig*math.sin(beta)],[-eta*sig*math.sin(beta),lam]])
           lam1=lam+eta*sig*math.cos(beta)
           c=[]
           for i in range(int(N)):
               for j in range(int(N)):
                   for k in range(len(Ja)):
                       for l in range(len(Ja)):
                            m[i][j][k][l]=I[i][j]*Ja[k][l]
                            n[i][j][k][l]=lam1*ve[i][j]*Jb[k][l]
                            o[i][j][k][l]=sig*adj_mat[i][j]*R[k][l]
           C=m-n+o
           #C=C.reshape((5,5,2,2))
           for i in range(int(N)):
               for k in range(len(Ja)):
                   for j in range(int(N)):
                       for l in range(len(Ja)):
                           c.append(C[i][j][k][l])
    
           c=np.array(c).reshape((int(N)*2,int(N)*2))
          # C=C.reshape((int(N)*2,int(N)*2))
           #print c
           #evals_large,evecs = eigsh(c, 1, which='LA', return_eigenvectors='False')
           evals, evecs=eig(c)
           #print evals_large
           evals.sort()
           #print evals
           hi.append(evals[(len(evals)-1)])
           #count+=1
           if (hi[count]<0):
               lambdasol2.append(lam)
               sigmasol2.append(sig)
           count+=1


#now for solutions for ring topology (no amplitude death) we use the matrix for synchronization (general) analysi find the eigee matrix

lambdasol3=np.empty((int(N),g*g))
sigmasol3=np.empty((int(N),g*g))

  

I=np.eye(int(N))
r =np.linspace(0.6,2.5,5)
#J2=np.array([[-2*ro*ro,0],[0,0]])
hello=[]
cos1=np.zeros((int(N),int(N)))
sin1=np.zeros((int(N),int(N)))
One=np.array([[1,0],[0,1]])
One1=np.array([[0,1],[1,0]])
final=np.zeros((int(N),int(N),2,2))
e=np.zeros((int(N),int(N),2,2))
f=np.zeros((int(N),int(N),2,2))
g=np.zeros((int(N),int(N),2,2))
h=np.zeros((int(N),int(N),2,2))
#M=np.zeros((int(N),int(N),2,2))
C=np.zeros((int(N),int(N)))
B=np.zeros((int(N),int(N)))
colors=cm.rainbow(np.linspace(0, 1, N))
#sigmm=[-2]

for kay in range(int(N)):
   # sigmasol3f=[]
   # lambdasol3f=[]    
    for i in range(int(N)):
        for j in range(int(N)):
            cos1[i][j]=math.cos(beta-2*math.pi*kay/int(N)*(i-j))
            sin1[i][j]=math.sin(beta-2*math.pi*kay/int(N)*(i-j))
            B[i][j]=cos1[i][j]*adj_mat[i][j]
            C[i][j]=sin1[i][j]*adj_mat[i][j]
    for k in range(int(N)):
        for i in range(len(One)):
            for j in range(len(One)):
                sum=0
                for n in range(int(N)):
                    if i==0 and j==1:
                        sum+=-adj_mat[k][n]*sin1[k][n]
                    elif i==1 and j==0:
                        sum+=adj_mat[k][n]*sin1[k][n]
                    else:
                        sum+=adj_mat[k][n]*cos1[k][n]
                h[k][k][i][j]=sum


    count=0
    for ro in r:
        J2=np.array([[-2*ro*ro,0],[0,0]])
        v=0
        for lam in lambd:
            for sig in sigma:
                M=[]
                for i in range(int(N)):
                    for j in range(int(N)):
                        for k in range(len(J2)):
                            for l in range(len(J2)):
                                e[i][j][k][l]=I[i][j]*J2[k][l]
                                f[i][j][k][l]=sig*B[i][j]*One[k][l]
                                g[i][j][k][l]=sig*C[i][j]*One1[k][l]
                            
                final=e+f+g-h
                for i in range(int(N)):
                    for k in range(len(J2)):
                        for j in range(int(N)):
                            for l in range(len(J2)):
                                M.append(final[i][j][k][l])
    
                M=np.array(M).reshape((int(N)*2,int(N)*2))
          
           #print c
           #evals_large,evecs = eigsh(c, 1, which='LA', return_eigenvectors='False')
                evals1, evecs1=eig(M)
         
                evals1.sort()
           #print evals
                hello.append(evals1[(len(evals1)-1)])
           #count+=1
                if (hello[count]<0):
                    lambdasol3[kay][v]=lam
                    sigmasol3[kay][v]=sig
                    v+=1
                count+=1
    
    #for i in range(len(sigmasol3[kay])):
       # sigmasol3f.append(sigmasol3[kay][i])
       # lambdasol3f.append(lambdasol3[kay][i])
for i in range(int(N)):
    plt.scatter(sigmasol3[i], lambdasol3[i], color=colors[i], alpha=0.3, label='cluster '+str(i))
       # plt.show()
       # plt.hold('True')

mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',colors)
sm = pyplot.cm.ScalarMappable(cmap=mymap, norm=plt.Normalize(vmin=1, vmax=int(N)))
sm._A = []

pyplot.colorbar(sm)

        #legend()
    
f = open('cluster.txt', 'a')
    
for i in range(int(N)):
    f.write('\n\n\n k is %d \n' %i)
    for sig,lam in zip(sigmasol3[i],lambdasol3[i]):
        f.write("%f" % sig) 
        f.write("%f \n" %lam)
f.close()

#np.savetxt('stability_clust.txt', zip(sigmasol3, lambdasol3), fmt="%f %f")
#np.savetxt('stability_sync.txt', zip(sigmasol, lambdasol), fmt="%f %f")

#print lambdasol2, len(lambdasol2)
plt.scatter(sigmasol1, lambdasol1, color='c', alpha=0.2, label='trivial steady state')
plt.scatter(sigmasol2, lambdasol2, color='b',  alpha=0.2, label='partial amplitude death')
plt.scatter(sigmasol, lambdasol, color='k',  alpha=0.2, label='synchronization')
# cluster plot is within the loop for different cluster states

#plt.clf()
#plt.imshow(heatmap)
#plt.show()
pyplot.xlim([-10, 10])
pyplot.ylim([-10, 10])
plt.ylabel('lambda')
plt.xlabel('sigma')
#plt.hold(True)
#plt.show()
legend()
plt.title('Stability', fontsize=13)
plt.savefig('stability.jpg')
plt.close()


#lambd = raw_input("enter lambda:")
#omega = raw_input("enter frequency:")
#N = raw_input("enter number of elements:")
lambd=-2.0
omega=2.0
#beta = float(raw_input("enter value of beta: "))
#N=5
#sigma = -5.0
sigmatotal = np.linspace(3,10, 8)

#tmax=100
#intstep=0.5


#adj_mat=matr() #adj.func(N) #this takes input matrix depending on motif from file adj.py
#print adj_mat
N=len(adj_mat[0])
#print "adjacency:" , adj_mat

def eiges():
    e_vals,e_vecs = LA.eig(adj_mat)
    return e_vals,e_vecs

def adj(i,j):
    return adj_mat[i][j]

def stulan(Z,t, sigma):
	dzdt = []
	#print N
	for i in range(int(N)):
		sum = 0.0
                if (np.iscomplexobj(Z)):
                    t=t+0j
	#	print "Z, t, i", Z, t, i
                #print type(Z)
		z = Z[i]
                for j in range(int(N)):
                    sum += adj(i,j)*(Z[j]-Z[i])
		c=(complex(lambd,omega)-pow(abs(z),2))*z + sigma*cmath.exp(beta*1j)*sum
               # print "c, type", c, type(c)
                
                '''if isinstance(c,complex):
                    print "zazzz"
                else:
                    print "hellooo"'''
                dzdt.append(c)
         #       print "typedzdt", type(dzdt)          
	#print "return"
	return dzdt


#Z0=eiges()[1][len(eiges()[1])/2]
#Z0=np.array([1.0/math.sqrt(5),cmath.exp((1j*6*math.pi)/5)/math.sqrt(5),cmath.exp((1j*12*math.pi)/5)/math.sqrt(5),cmath.exp((1j*18*math.pi)/5)/math.sqrt(5),cmath.exp((1j*24*math.pi)/5)/math.sqrt(5)])
#Z0 = np.array([1.0/math.sqrt(5),cmath.exp(1j*2*math.pi/5)/math.sqrt(5),cmath.exp(1j*4*math.pi/5)/math.sqrt(5),cmath.exp(1j*6*math.pi/5)/math.sqrt(5),cmath.exp(1j*8*math.pi/5)/math.sqrt(5)], dtype=np.complex128) #initialize with as many terms as N
Z0 = np.array([1.0,-1.0,0.0,1.0,-1.0], dtype=np.complex128)

#print type(Z0)
#tstart = np.zeros(np.ceil(tmax/intstep)+1)
#t=tstart
time = np.linspace(0, 80, 1600) #ideal for 1000 t is 40000 and 1000
t=time[::1]
time1= np.linspace(0, 2000, 80000)
t1=time1[::20]
np.savetxt('time.txt',t , fmt="%f")
np.savetxt('time1.txt',t1 , fmt="%f")

    
#colors=['b','g','r','c','m','y','k'] #add colors if N greater than 7
colors=cm.rainbow(np.linspace(0, 1, N))


def myodeint(stulan, Z0, t):
    Z0 = np.array(Z0, complex)
    func = lambda t, Z: stulan(Z, t)   # odeint has these in reverse
    sol = ode(func).set_integrator('zvode').set_initial_value(Z0, t=t[0])
    Z = [sol.integrate(tp) for tp in t[1:]]
    #print Z, t
   # sys.stdout.flash()
    Z.insert(0, Z0)
  
    my_list=collections.defaultdict(list)
    #l,m,n,o,p,q,r,s,t,u,v,w=[],[],[],[],[],[],[],[],[],[],[],[]
  
    for row in Z:
        for i in range(int(N)):
            my_list[i].append(real(row[i]))
    
    mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',colors)
    for i in range(int(N)):
        plt.plot(t[:400],my_list[i][1200:1600],color=colors[i], label=str(i))
       # plt.show()
       # plt.hold('True')
    sm = pyplot.cm.ScalarMappable(cmap=mymap, norm=plt.Normalize(vmin=1, vmax=5))
    sm._A = []
    pyplot.colorbar(sm)
    
    
    
        #legend()
    
    f = open('ring hierarch'+str(sigma)+'.txt', 'a')
    
    for i in range(int(N)):
        f.write('\n\n\n')
        for j in my_list[i]:
            ho=str(j)
            f.write(ho + "\n")
    f.close()
    
   # figure()
    plt.xlabel('$t$')
    plt.ylabel('$Re(Z)$')
   # plt.ylim((-10,10))
    #plt.hold('True')
    #plt.show()
    #legend()
  
    plt.savefig('stulan' + str(sigma) +'.png')
    plt.close()
    return np.array(Z)

def myodeint1(stulan, Z0, t):
    Z0 = np.array(Z0, complex)
    func = lambda t, Z: stulan(Z, t)   # odeint has these in reverse
    sol = ode(func).set_integrator('zvode').set_initial_value(Z0, t=t1[0])
    Z1 = [sol.integrate(tp) for tp in t1[1:]]
    #print Z, t
   # sys.stdout.flash()
    Z1.insert(0, Z0)
  
    my_list1=collections.defaultdict(list)
    #l,m,n,o,p,q,r,s,t,u,v,w=[],[],[],[],[],[],[],[],[],[],[],[]
  
    for row in Z1:
        for i in range(int(N)):
            my_list1[i].append(real(row[i]))
    
    mymap = mpl.colors.LinearSegmentedColormap.from_list('mycolors',colors)
    for i in range(int(N)):
        plt.plot(t1[3500:4000],my_list1[i][3500:4000],color=colors[i], label=str(i))
       # plt.show()
       # plt.hold('True')
    sm = pyplot.cm.ScalarMappable(cmap=mymap, norm=plt.Normalize(vmin=1, vmax=int(N)))
    sm._A = []
    pyplot.colorbar(sm)
    
    
    
        #legend()
    
    f = open('ring hierarchagain'+str(sigma)+'.txt', 'a')
    
    for i in range(int(N)):
        f.write('\n\n\n')
        for j in my_list1[i]:
            ho=str(j)
            f.write(ho + "\n")
    f.close()
    
   # figure()
    plt.xlabel('$t$')
    plt.ylabel('$Re(Z)$')
    #plt.hold('True')
    #plt.show()
    #legend()
  
    plt.savefig('stulanagain' + str(sigma) +'.png')
    plt.close()
    return np.array(Z1)


'''s.set_integrator('dopri5', nsteps=1000)
count=0
time=[]
while s.successful() and s.t<tmax-intstep:
	count+=1
	s.integrate(s.t+intstep)
	time.append(s.t)
	print("%g %g" % (s.t, s.Z))
sol = odeint(stulan, t, Z)'''

for k in range(len(sigmatotal)):
    sigma=sigmatotal[k]
#z=odeint(stulan,Z0,t)
    #print ("hello")
    Z = myodeint(lambda Z, t: stulan(Z, t, sigma), Z0, t)
   # Z1 = myodeint1(lambda Z, t: stulan(Z, t, sigma), Z0, t)
#print type(Z)
#igure()
#lt.plot(t, Z)
#show()

