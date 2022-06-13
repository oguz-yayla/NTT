import math
import numpy as np
import time
import random


#[0]=n, [1]=p, [2]=w, [3]=r, [4]=r^-1
param1=[3,19,7,4,5]                 # f1=x^3-7,    and f2=x^3-11     , alpha=7    , beta=11                 
param2=[9,109,16,9,97]              # f1=x^9-45    and f2=x^9-63     , alpha=45   , beta=63
param3=[27,163,64,161,81]           # f1=x^27-58   and f2=x^27-64    , alpha=58   , beta=64
param4=[81,487,242,4,122]           # f1=x^81-232  and f2=x^81-254   , alpha=232  , beta=254
param5=[243,1459,729,9, 1297]       # f1=x^243-339 and f2=x^243-1119 , alpha=339  , beta=1119
param6=[729,17497,12013,4, 13123]   # f=x^729-4719 and f2=x^729-12777, alpha=4719 , beta=12777

param811459=[81,1459,547,729, 1457]        # f1=x^81-339 and f2=x^81-1119, alpha=339 , beta=1119
param812917=[81,2917,764,2890,108]         # f1=x^81-247 and f2=x^81-2669, alpha=247 , beta=2669
param2432917=[243, 2917, 1040, 2914, 972]  # f1=x^243-247 and f2=x^243-2669 , alpha=247 , beta=2669 


param=param812917


############################ Parameters ##########################

n=param[0] # n 
p=param[1] # prime
w=param[2] # nth root of unity 
r=param[3] # nth root of alpha
s=param[4] # s=r^-1, also nth root of beta

winv=inverse_mod(w,p)

mu=w**(n/3) % p
mu2=mu**2 % p 

l=log(n,3)

gamma=[[r*w**i%p for i in range(n/3)],[s*w**i%p for i in range(n/3)]]
gamma_inverse=[[r*winv^i%p for i in range(n/3)],[s*winv^i%p for i in range(n/3)]]

inv3=inverse_mod(3,p)

wk=[0,0] 
wk[0]=[0]*l
wk[1]=[0]*l

for level in range(1,l+1):   
    m=3**level
    wk[0][level]=[0]*(m/3)
    wk[1][level]=[0]*(m/3)

    for j in range(0,m/3): 
        wk[0][level][j]=(gamma[0][j]**(n/m))%p
        wk[1][level][j]=(gamma[1][j]**(n/m))%p

# r and w are precomputed


############################## MAIN ##############################

def trit_reverse(i,ntrit):
    rev=np.base_repr(i,3)[::-1]
    rev=rev+(ntrit - len(rev))*'0'
    return int(rev,3)

def madd(a,b):
    return (a+b) % p 

def mmult(a,b):
    return (a*b) % p 

def scramble(a,n):
    b=[0]*n
    for i in range(n):
        b[i]=a[trit_reverse(i,int(math.log(n,3)))]
    return b

# base level multiplication in Z[x]/<x^3-r>
def basemul2(a,b,r):
    c=[0,0,0]
    c[0]= (a[0]*b[0]+(a[1]*b[2]+a[2]*b[1])*r)%p
    c[1]= (a[0]*b[1]+a[2]*b[2]*r+a[1]*b[0])%p
    c[2]= (a[0]*b[2]+a[1]*b[1]+a[2]*b[0])%p
    return c    

def butterfly(a,b,c,w,kind):
    result=[0,0,0] 
    w2=w**2 % p
    if kind == 0:
        wb =mmult(b,w)
        w2c=mmult(c,w2)
        muwb=mmult(mu,wb)
        mu2w2c=mmult(mu2,w2c)
        mu2wb=mmult(mu2,wb)
        muw2c=mmult(mu,w2c)
        result[0]=madd(a, madd(wb,w2c))        #exp: a+w^k b+w^2k c
        result[1]=madd(a, madd(muwb,mu2w2c))   #exp: a+mu w^k b+mu^2 w^2k c
        result[2]=madd(a, madd(mu2wb,muw2c))   #exp: a+mu^2 w^k b+mu w^2k c
    else:
        result[0]= mmult(inv3, madd(a,madd(b,c)))                               #exp: 1/3 * (a+b+c)
        result[1]= mmult(mmult(inv3,w),madd(a,madd(mmult(mu2,b),mmult(mu,c))))  #exp: 1/3 * 1/w^k * (a+mu2*b+mu*c)
        result[2]= mmult(mmult(inv3,w2),madd(a,madd(mmult(mu,b),mmult(mu2,c)))) #exp: 1/3 * 1/w^(2k) * (a+mu*b+mu2*c)
    return result
        
def NTT(a,ring):
    A=scramble(a,n)
    for level in range(1,l+1):   
        m=3**level
        for j in range(0,m/3):   
            wk=(gamma[ring][j]**(n/m))%p
            for k in range(0,n/m):  
                [A[k*m+j], A[k*m+j+m/3],A[k*m+j+2*m/3]]=butterfly(A[k*m+j],A[k*m+j+m/3],A[k*m+j+2*m/3],wk,0)
    return A

def INTT(a,ring):
    for level in range(l,0,-1):
        m=3**level
        for j in range(0,m/3):
            wkinv=(gamma_inverse[ring][j]^(n/m))%p  
            for k in range(0,n/m):
                [a[k*m+j], a[k*m+j+m/3], a[k*m+j+2*m/3]]=butterfly(a[k*m+j],a[k*m+j+m/3],a[k*m+j+2*m/3],wkinv,1)
    A = scramble(a,n)
    return A
       
def pmult(NTTa,NTTb):
    return [mmult(NTTa[i],NTTb[i]) for i in range(n)]

def crt(a,zeta):
    a_zeta=[0]*n
    for i in range(n):
        a_zeta[i]=a_zeta[i] + (a[i]+zeta*a[n+i])%p
    return a_zeta

def icrt(c_alpha,c_beta,alpha,beta):
    
    def comb(A,d,t):
        A += [0]*n
        tAxd = [t*j%p for j in list(np.array(A[n:]+A[:n]) - np.array([i*d%p for i in A]))]
        return tAxd
    
    t1=inverse_mod(alpha-beta,p)
    t2=inverse_mod(beta-alpha,p)
    c=[i%p for i in list(np.array(comb(c_alpha,beta,t1))+np.array(comb(c_beta,alpha,t2)))]
    return c

# -------------------------------------------------------------------

##uncomment the parameters according to the ring.
## for Z_1459
#alpha=339 
#beta=1119

## for Z_2917
alpha=247
beta=2669

runtime=[0]*10000
for i in range(10000):
    a=[random.randrange(0, p) for i in range(162)]
    b=[random.randrange(0, p) for i in range(162)]
    t1=time.time()
    c_alpha=INTT(pmult(NTT(crt(a,alpha),0),NTT(crt(b,alpha),0)),1)
    c_beta =INTT(pmult(NTT(crt(a,beta),1),NTT(crt(b,beta),1)),0)
    icrt(c_alpha,c_beta,alpha,beta);
    t2=time.time()
    runtime[i]=t2-t1

print(sum(runtime)/10000, max(runtime), min(runtime))
