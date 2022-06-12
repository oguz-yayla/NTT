import math, time


def bit_reverse(i,nbit):
    reverse=bin(i)[:1:-1]
    reverse=reverse+(nbit - len(reverse))*'0'
    return int(reverse,2)

def scramble(a,N):
    b=[0]*N
    for i in range(N):
        b[i]=a[bit_reverse(i,int(math.log(N,2)))]
    return b

def moduler_addition(a,b):
    return (a+b) % 3329

def moduler_substraction(a,b):
    return (a-b) % 3329

def moduler_multiplication(a,b):
    return a*b % 3329

def basemul2(a,b,r):
    c=[0,0,0]
    c[0]= (a[0]*b[0]+(a[1]*b[2]+a[2]*b[1])*r)%p
    c[1]= (a[0]*b[1]+a[2]*b[2]*r+a[1]*b[0])%p
    c[2]= (a[0]*b[2]+a[1]*b[1]+a[2]*b[0])%p
    return c    

def butterfly(a,b,w,kind):
    result=[0,0]
    if kind==0:
        result[0]=moduler_addition(a,moduler_multiplication(b,w))
        result[1]=moduler_substraction(a,moduler_multiplication(b,w))
    else:
        result[0]=moduler_multiplication(1665, (moduler_addition(a,b)))
        result[1]=moduler_multiplication(1665, (moduler_multiplication(moduler_substraction(a,b),w)))
    return result

def NTT(a):
    N=256 #predefined
    q=3329 #predefined
    gamma=[(17**i) % q for i in range(N)] #predefined
    A=scramble(a,N)  
    for s in range(1,8):
        m=2**s
        for j in range(0,m//2):
            w=gamma[(2*j+1)*N/m]
            for k in range(0,N/m):
                [A[k*m+j], A[k*m+j+m/2]]=butterfly(A[k*m+j],A[k*m+j+m/2],w,0)
    return A
def INTT(a):
    N=256 #predefined
    q=3329 #predefined
    gamma_inverse=[(1175**i) % q for i in range(N)] #predefined
    for s in range(8,1,-1):
        m=2**s
        for j in range(0,m/2):
            w=gamma_inverse[(2*j+1)*N/m]
            for k in range(0,N/m):
                [a[k*m+j], a[k*m+j+m/2]]=butterfly(a[k*m+j],a[k*m+j+m/2],w,1)
    A=scramble(a,N)
    return A

def pmult(NTTa,NTTb):
    N=256
    return [moduler_multiplication(NTTa[i],NTTb[i]) for i in range(N)]
