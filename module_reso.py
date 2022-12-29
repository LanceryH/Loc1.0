import numpy as np
from tqdm import tqdm

def force(m1,X1,m2,X2) :
    G=4*np.pi**2    
    d=np.linalg.norm(X1-X2)
    f=-(G*m1*m2/(d**3))*(X1-X2)
    return f
    
def f(t,y,M) :
    nb_corps=len(M)
    F=np.zeros((nb_corps*6,1))
    for i in range(nb_corps) :
        S=np.zeros((3,1))
        for j in range(nb_corps) : 
            #print(y[i*3:(i+1)*3])
            if j != i:
                S+=force(M[i], y[i*3:(i+1)*3,0].reshape(-1,1),M[j],y[j*3:(j+1)*3,0].reshape(-1,1))
            
        F[i*3:(i+1)*3,0]=y[nb_corps*3+i*3:(i+1)*3+nb_corps*3,0]  
        F[nb_corps*3+i*3:(i+1)*3+nb_corps*3,0]=(1/M[i]*S).reshape(-1)    
    return F    

def Rk4(y,t0,tf,N,nb_corps,M):
    h=(tf-t0)/N
    t=t0
    Y=np.zeros((3*nb_corps,N))
    for i in tqdm(range(0,N)) : 
        k1=h*f(t,y,M)
        k2=h*f(t+h/2, y+k1/2,M)
        k3=h*f(t+h/2, y+k2/2,M)
        k4=h*f(t+h, y+k3,M)
        y=y+(1/6)*(k1 +2*k2 +2*k3 + k4)
        t+=h    
        for j in range(0,nb_corps) :
            Y[j*3:(j+1)*3,i]=y[j*3:(j+1)*3,0]
    return(Y)   

def adaptive_rkf45(y, t0, tf, N, nb_corps, M):
    h=(tf-t0)/N 
    t=t0
    ErreurAdmis=1e-6
    Y=np.zeros((3*nb_corps,N))
    for i in tqdm(range(0,N)) : 
        k1 = f(t, y, M)
        k2 = f(t+h/5, y+h*k1/5, M)
        k3 = f(t+3*h/10, y+3*h*(k1+3*k2)/40, M)
        k4 = f(t+3*h/5, y+h*((k1*3/10)-(k2*9/10)+(k3*6/5)), M)
        k5 = f(t+h, y+h*((-k1*11/54)+(k2*5/2)-(k3*70/27)+(k4*35/27)), M)
        k6 = f(t+h*7/8, y+h*((k1*1631/55296)+(k2*175/512)+(k3*575/13824)+(k4*44275/110592)+(k5*253/4096)), M)
        t+=h    
        ya=y+np.array(h*(37*k1/378+250*k3/621+125*k4/594+512*k6/1771))
        yb=y+np.array(h*(2825*k1/27648+18575*k3/48384+13525*k4/55296+277*k5/14336+k6/4))
        ErreurActuel=np.linalg.norm(np.abs(ya-yb))
        if (ErreurAdmis<=ErreurActuel):
            h=h*(ErreurAdmis/ErreurActuel)**0.2
        if (ErreurAdmis>ErreurActuel and ErreurActuel>0.0):
            h=h*(ErreurAdmis/ErreurActuel)**0.25
        if (h>(tf-t0)/N):
            h=(tf-t0)/N     
        for j in range(0,nb_corps):
            Y[j*3:(j+1)*3,i]=yb[j*3:(j+1)*3,0] 
        y=yb    
        print("etape:{} et pas:{}".format(i,h))    
    return Y


def Rk2(y,t0,tf,N,nb_corps,M):
    h=(tf-t0)/N
    t=t0
    Y=np.zeros((3*nb_corps,N))
    for i in tqdm(range(0,N)) : 
        k1=h*f(t,y,M)
        k2=h*f(t+h, y+k1,M)
        y=y+(k1+k2)/2
        t+=h    
        for j in range(0,nb_corps) :
            Y[j*3:(j+1)*3,i]=y[j*3:(j+1)*3,0]
    return(Y)  

def Verlet(y,t0,tf,N,nb_corps,M):
    h=(tf-t0)/N
    t=t0
    Y=np.zeros((3*nb_corps,N))
    for i in range(0,N) : 
        k1=h*f(t,y,M)
        k2=h*f(t+h, y+k1,M)
        y=y+(k1+k2)/2
        t+=h    
        for j in range(0,nb_corps) :
            Y[j*3:(j+1)*3,i]=y[j*3:(j+1)*3,0]
    return(Y)  
