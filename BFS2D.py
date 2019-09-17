# -*- coding: utf-8 -*-
"""
Created on Sun May 19 21:53:44 2019

@author: harsh arora
"""

import numpy as np
import matplotlib.pyplot as plt

L=1.0
Lc=1.0
L2=2.0
L3=10.0
dx=0.1
dy=0.1
nu=1.0e-3
dt=0.1
dt_stream = 0.002
epsilon_stream = 1.0e-2
eps_vort=1.0e-2

N=int((L2+L3)/dx)+1
M=int((Lc+L)/dy)+1

u0=4.0

Re=(u0/6.0)/nu

# Step 1 - initial and boundary conditions (that do not change with time)
# (a) initialise u,v, stream, vort at t=0
vort_n = np.zeros((N,M))
vort_np1 = np.zeros((N,M))

stream_k = np.zeros((N,M))
stream_kp1 = np.zeros((N,M))

u= np.zeros((N,M))
v= np.zeros((N,M))

x=np.linspace(0,12.,N)
y=np.linspace(0,2.,M)
xv,yv=np.meshgrid(x,y)

# 1(b) Inlet BC for u, v, stream and vort
i=0
for j in range(10,M):
    y=(j-10)*dy
    u[i][j]=u0*(y/L)*(1.0-y/L)
    v[i][j]=0.0
    stream_k[i][j]=stream_kp1[i][j]=(u0/L)*(y**2/2. - y**3./(3.*L))
    vort_n[i][j]=vort_np1[i][j]=-(u0/L)*(1.-2*(y/L))
    
# (c) All walls BC for u and v
    #leave it as initialised
# (d) Bottom walls BC for stream
    #leave it as initialised
# (e) stream on top wall
j=M-1
for i in range(1,N-1):
    stream_k[i][j]=stream_kp1[i][j]=u0/6.0

#Step 2 - Time loop
for n in range(100):
 # 2(a) Solve vort equation
# 2(a)(i) Find vort_np1 on walls
# Top wall
    j=M-1
    for i in range(1,N-1):
        vort_np1[i][j]=-(2.0/(dy*dy))*(stream_k[i][j-1]-stream_k[i][j])
# Bottom wall
    j=0
    for i in range(20,N-1):   
        vort_np1[i][j]=-(2.0/(dy*dy))*(stream_k[i][j+1]-stream_k[i][j])
# Left wall
    i=20
    for j in range(1,11):
        vort_np1[i][j]=-(2.0/(dx*dx))*(stream_k[i+1][j]-stream_k[i][j])
        
# Left bottom wall
    j=10
    for i in range(1,20):   
        vort_np1[i][j]=-(2.0/(dy*dy))*(stream_k[i][j+1]-stream_k[i][j])        
        
# 2(a)(ii) Find vort_np1 for all interior points
    for i in range(21,N-1):
        for j in range(1,11):
            ue=(u[i][j]+u[i+1][j])/2.0
            uw=(u[i-1][j]+u[i][j])/2.0
            vn=(v[i][j]+v[i][j+1])/2.0
            vs=(v[i][j-1]+v[i][j])/2.0
            
            if(ue>=0):
                we=vort_n[i][j]
            else: 
                we=vort_n[i+1][j]
                
            if(uw>=0):
                ww=vort_n[i-1][j]
            else: 
                ww=vort_n[i][j]
                
            if(vn>=0):
                wn=vort_n[i][j]
            else: 
                wn=vort_n[i][j+1]
                
            if(vs>=0):
                ws=vort_n[i][j-1]
            else: 
                ws=vort_n[i][j]    

            vort_np1[i][j]=((nu*dt)/(dx*dx))*(vort_n[i-1][j]-2*vort_n[i][j]+vort_n[i+1][j])+((nu*dt)/(dy*dy))*(vort_n[i][j-1]-2*vort_n[i][j]+vort_n[i][j+1])\
                            -ue*we*(dt/dx)+uw*ww*(dt/dx)-vn*wn*(dt/dy)+vs*ws*(dt/dy) + vort_n[i][j]

    for i in range(1,N-1):
        for j in range(11,M-1):
            ue=(u[i][j]+u[i+1][j])/2.0
            uw=(u[i-1][j]+u[i][j])/2.0
            vn=(v[i][j]+v[i][j+1])/2.0
            vs=(v[i][j-1]+v[i][j])/2.0
            
            if(ue>=0):
                we=vort_n[i][j]
            else: 
                we=vort_n[i+1][j]
                
            if(uw>=0):
                ww=vort_n[i-1][j]
            else: 
                ww=vort_n[i][j]
                
            if(vn>=0):
                wn=vort_n[i][j]
            else: 
                wn=vort_n[i][j+1]
                
            if(vs>=0):
                ws=vort_n[i][j-1]
            else: 
                ws=vort_n[i][j]    

            vort_np1[i][j]=((nu*dt)/(dx*dx))*(vort_n[i-1][j]-2*vort_n[i][j]+vort_n[i+1][j])+((nu*dt)/(dy*dy))*(vort_n[i][j-1]-2*vort_n[i][j]+vort_n[i][j+1])\
                            -ue*we*(dt/dx)+uw*ww*(dt/dx)-vn*wn*(dt/dy)+vs*ws*(dt/dy) +vort_n[i][j]

# 2(a)(iii) Find vort_np1 for outlet boundary points
    i=N-1
    for j in range(0,M):
        vort_np1[i][j]=vort_np1[i-1][j]
        
# 2(b) Solve stream eqn for al interior and boundary points
    for k in range(1000):
        # Interior part I
        for i in range(21,N-1):
            for j in range(1,11):
                stream_kp1[i][j] = stream_k[i][j] + vort_np1[i][j]*dt_stream + ((stream_k[i-1][j]-2.*stream_k[i][j]+stream_k[i+1][j])*(dt_stream))/(dx*dx)\
                                    + ((stream_k[i][j-1]-2.*stream_k[i][j]+stream_k[i][j+1])*(dt_stream))/(dy*dy)
        # Interior part II
        for i in range(1,N-1):
            for j in range(11,M-1):
                stream_kp1[i][j] = stream_k[i][j] + vort_np1[i][j]*dt_stream + ((stream_k[i-1][j]-2.*stream_k[i][j]+stream_k[i+1][j])*(dt_stream))/(dx*dx)\
                                    + ((stream_k[i][j-1]-2.*stream_k[i][j]+stream_k[i][j+1])*(dt_stream))/(dy*dy)
        i=N-1
        for j in range(0,M):
            stream_kp1[i][j] = stream_kp1[i-1][j]
        # k=k+1    
        error = 0.0
        for i in range(0,N):
            for j in range(0,M):
                error = error + (stream_kp1[i][j] - stream_k[i][j])**2
                stream_k[i][j] = stream_kp1[i][j]
        error = np.sqrt(error/(M*N))
        print("timestep", n, "error=", error/epsilon_stream, "iteration number", k)
        if (error/dt_stream)<epsilon_stream:
            break

# 2(c) Get u and v from stream_kp1
    # Interior part I
    for i in range(21,N-1):
        for j in range(1,11):
            u[i][j]=(stream_kp1[i][j+1]-stream_kp1[i][j-1])/(2.0*dy)
            v[i][j]=-(stream_kp1[i+1][j]-stream_kp1[i-1][j])/(2.0*dx)
    # Interior part II
    for i in range(1,N-1):
        for j in range(11,M-1):
            u[i][j]=(stream_kp1[i][j+1]-stream_kp1[i][j-1])/(2.0*dy)
            v[i][j]=-(stream_kp1[i+1][j]-stream_kp1[i-1][j])/(2.0*dx)
            
    i=N-1
    for j in range(1,M-1):
        u[i][j]=(stream_kp1[i][j+1]-stream_kp1[i][j-1])/(2.0*dy)
        v[i][j]=-(stream_kp1[i][j]-stream_kp1[i-1][j])/(dx)
    #n=n+1
    for i in range(1,N):
        for j in range(0,M):
            vort_n[i][j]=vort_np1[i][j]
            
    if (n%20==0):
        nn=n/20
        f, axarr = plt.subplots(2, 2)
        axarr[0, 0].contour(xv,yv,np.transpose(stream_kp1),30)
        axarr[0, 0].set_title('Streamlines')
        axarr[0, 1].plot(u[:,1])
        axarr[0, 1].set_title('u-vel above bottom wall')
        axarr[0, 1].grid()
        axarr[1, 0].quiver(xv,yv,np.transpose(u),np.transpose(v),scale=15)
        axarr[1, 0].set_title('Velocity vectors')
        axarr[1, 1].contourf(xv,yv,np.transpose(vort_np1),50)
        axarr[1, 1].set_title('Vorticity')
        f.suptitle('t='+str(n*dt)+' s', fontsize=18)
        
        #plt.contour(xv,yv,np.transpose(stream_np1),30)
        axarr[0, 0].set_aspect("equal")
        axarr[1, 0].set_aspect("equal")
        axarr[1, 1].set_aspect("equal")
        plt.savefig('test'+str(nn)+'.png', bbox_inches='tight', dpi=400)
    
    
            
            