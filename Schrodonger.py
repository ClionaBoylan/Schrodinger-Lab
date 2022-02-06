import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import simps

N=1000
gamma2=200
l=1.0/(N-1)
psi=np.zeros(N)
potentials=-1*np.ones(N)
x=np.linspace(0,1,N)

psi[0]=0
psi[1]=1e-4
trial=-0.96

k2=[]
i=0
while i<N:
    k2.append(gamma2*(trial-potentials[i]))
    i+=1

#Correct answers for trial energy
trials=[]
i=0
while i<10:
    trials.append((i**2*np.pi**2)/gamma2-1)
    i+=1

#Function to find wavefunction
def findwavefunction(trial):
    k2=[]
    i=0
    while i<N:
        k2.append(gamma2*(trial-potentials[i]))
        i+=1
    
    psi=np.zeros(N)
    psi[0]=0
    psi[1]=1e-4
    
    i=1
    while i<N-1:
        psi[i+1]=(2.0*(1.0-(5.0/12.0)*(l**2)*k2[i])*psi[i]-(1.0+(1.0/12.0)*(l**2)*(k2[i-1]))*psi[i-1])/(1.0+(1.0/12.0)*(l**2)*(k2[i+1]))
        i+=1
    
    return psi

#First calculation
psi=findwavefunction(trial)

print("Last psi",  psi[N-1], "FIRST TRIAL")

# plt.plot(psi)
# plt.show()

plt.plot(x,findwavefunction(-0.98), label = r"$\epsilon=-0.99$",color='#D62800')
plt.plot(x,findwavefunction(-0.97), label = r"$\epsilon=-0.98$",color='#FF9B56')
plt.plot(x,findwavefunction(-0.96), label = r"$\epsilon=-0.97$",color='#D462A6')
plt.plot(x,findwavefunction(-0.95), label = r"$\epsilon=-0.95$",color='#A40062')
plt.xlabel(r'$\tilde{x}$')
plt.ylabel(r'$\psi(\tilde{x})$')
plt.legend(loc="upper left")
plt.savefig("Driveby.pdf")
plt.show()

accuracy=1.0e-12 

#Function to find eigenvalue
def findeigenvalue(start,accuracy):
    delta=1.0e-2
    psi=findwavefunction(start)
    sign=np.sign(psi[N-1])
    while abs(psi[N-1])>accuracy:
        psi=findwavefunction(start+delta)
        if np.sign(psi[N-1])!=sign:
            delta=-delta/2.0
        else:
            delta=delta
        sign=np.sign(psi[N-1])
        start+=delta
    return start

eigenvalue=findeigenvalue(trial,accuracy)
print("Last psi:",psi[N-1],"AFTER LOOP")
print("Epsilon",trial,"AFTER LOOP")

def findeigenvalues(start,number,accuracy):
    length=number
    i=0
    eigenvalues=[findeigenvalue(start,accuracy)]
    while i<length-1:
        eigenvalues.append(findeigenvalue(eigenvalues[i]+0.001,accuracy))
        i+=1
    return eigenvalues

eigenvalues1=findeigenvalues(-1.0,10,1.0e-12)

psi=findwavefunction(eigenvalues1[0])
A=np.sqrt(1.0/(simps(psi**2,dx=l)))
print("Normalising factor", A)
normalpsi=A*psi

i=0
analytic=[]
while i<N:
    analytic.append(np.sqrt(2)*np.sin(np.pi*x[i]))
    i+=1
plt.plot(x,analytic-normalpsi,color='#D60270')
plt.xlabel(r'$\tilde{x}$')
plt.ylabel(r'Analytic $\psi(\tilde{x})$ - Numerical $\psi(\tilde{x})$')
plt.savefig('Error1.pdf')
plt.show()

psi=findwavefunction(eigenvalues1[5])
A=np.sqrt(1.0/(simps(psi**2,dx=l)))
print("Normalising factor", A)
normalpsi=A*psi

i=0
analytic=[]
while i<N:
    analytic.append(np.sqrt(2)*np.sin(6*np.pi*x[i]))
    i+=1

plt.plot(x,analytic-normalpsi,color='#9B4F96')
plt.xlabel(r'$\tilde{x}$')
plt.ylabel(r'Analytic $\psi(\tilde{x})$ - Numerical $\psi(\tilde{x})$')
plt.savefig('Error2.pdf')
plt.show()

i=0
colors=(np.linspace(216/255,36/255,len(eigenvalues1)),np.linspace(9/255,70/255,len(eigenvalues1)),np.linspace(126/255,142/255,len(eigenvalues1)))
while i<len(eigenvalues1):
    psi=findwavefunction(eigenvalues1[i])
    A=np.sqrt(1.0/(simps(psi**2,dx=l)))
    normalpsi=A*psi
    plt.plot(x,normalpsi,color=(colors[0][i],colors[1][i],colors[2][i]),zorder=2)
    plt.xlabel(r'$\tilde{x}$')
    plt.ylabel(r'$\psi(\tilde{x})$')
    plt.title('Normalised Wavefunction for $n=%s$' %(i+1))
    plt.hlines(0,0,1,color='black',zorder=1,alpha=0.9)
    plt.xlim(0,1)
    plt.savefig('Square Well %s.pdf' %i)
    plt.show()
    i+=1

def uncertainty():
    x=np.linspace(0,1,N)
    firsttermdeltax=simps(x**2*normalpsi**2,dx=l)
    secondtermdeltax=simps(x*normalpsi**2,dx=l)
    deltax=np.sqrt(firsttermdeltax-secondtermdeltax**2)

    psidotdot=np.zeros(N)
    i=1
    while i<N-1:
        psidotdot[i]=(normalpsi[i-1]-2.0*normalpsi[i]+normalpsi[i+1])/(l**2)
        i+=1

    deltap=np.sqrt(-simps(normalpsi*psidotdot,dx=l))
    
    return deltax, deltap, deltax*deltap

uncertainties1=[]
i=0
while i<len(eigenvalues1):
    psi=findwavefunction(eigenvalues1[i])
    A=np.sqrt(1.0/(simps(psi**2,dx=l)))
    normalpsi=A*psi
    uncertainties1.append(uncertainty()[2])
    i+=1
i=0
diff1=[]
while i<len(eigenvalues1)-1:
    diff1.append(eigenvalues1[i+1]-eigenvalues1[i])
    i+=1
    
index=[]
i=0
while i<len(uncertainties1):
    index.append(i+1)
    i+=1
plt.scatter(index,uncertainties1,zorder=2,color='red')
plt.plot(index,uncertainties1,zorder=1,color='black')
plt.xlabel(r'$n$')
plt.ylabel(r'$\Delta \tilde{x} \Delta \tilde{p}$')
plt.hlines(0.5,index[0],index[-1],color='black',linestyle='dashed',zorder=0)
plt.savefig('UncertaintySquare.pdf')
plt.show()


i=0
potentials=np.zeros(N)
while i<N:
    potentials[i]=8.0*(x[i]-0.5)**2-1
    i+=1

i=0
potentials=np.zeros(N)
while i<N:
    potentials[i]=-np.sin(np.pi*x[i])+1
    i+=1
    
i=0
potentials=np.zeros(N)
while i<N:
    potentials[i]=-np.cos(2*np.pi*x[i])+2
    i+=1

i=0
potentials=np.zeros(N)
while i<N:
    potentials[i]=np.tan((0.9/2)*np.pi*x[i])
    i+=1

gamma2=1000
k2=[]
i=0
while i<N:
    k2.append(gamma2*(trial-potentials[i]))
    i+=1

eigenvalues2=findeigenvalues(2,20,1.0e-7)
uncertainties2=[]
i=0
while i<len(eigenvalues2):
    psi=findwavefunction(eigenvalues2[i])
    A=np.sqrt(1.0/(simps(psi**2,dx=l)))
    normalpsi=A*psi
    uncertainties2.append(uncertainty()[2])
    i+=1

i=0
diff2=[]
while i<len(eigenvalues2)-1:
    diff2.append(eigenvalues2[i+1]-eigenvalues2[i])
    i+=1

index=[]
i=0
while i<len(uncertainties2):
    index.append(i+1)
    i+=1
plt.scatter(index,uncertainties2,zorder=2,color='red')
plt.plot(index,uncertainties2,zorder=1,color='black')
plt.xlabel(r'$n$')
plt.ylabel(r'$\Delta \tilde{x} \Delta \tilde{p}$')
plt.hlines(0.5,index[0],index[-1],color='black',linestyle='dashed',zorder=0)
plt.savefig('UncertaintyHarmonic.pdf')
plt.show()

i=0
colors=(np.linspace(216/255,36/255,len(eigenvalues1)),np.linspace(9/255,70/255,len(eigenvalues1)),np.linspace(126/255,142/255,len(eigenvalues1)))
while i<len(eigenvalues1):
    psi=findwavefunction(eigenvalues2[i])
    A=np.sqrt(1.0/(simps(psi**2,dx=l)))
    normalpsi=A*psi
    plt.plot(x,normalpsi,color=(colors[0][i],colors[1][i],colors[2][i]),zorder=2)
    plt.xlabel(r'$\tilde{x}$')
    plt.ylabel(r'$\psi(\tilde{x})$')
    plt.title('Normalised Wavefunction for $n=%s$' %(i+1))
    plt.hlines(0,0,1,color='black',zorder=1,alpha=0.9)
    plt.xlim(0,1)
    plt.savefig('Harmonic %s.pdf' %i)
    plt.show()
    i+=1

index1=[]
i=0
while i<len(diff1):
    index1.append(i+1)
    i+=1

index2=[]
i=0
while i<len(diff2):
    index2.append(i+1)
    i+=1

plt.loglog(index1,diff1,color='black',zorder=1)
plt.scatter(index1,diff1,color='red',zorder=2)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$n$')
plt.ylabel(r'$\Delta \epsilon$')
plt.savefig('UncLogLog1.pdf')
plt.show()

plt.loglog(index2,diff2,color='black',zorder=1)
plt.scatter(index2,diff2,color='red',zorder=2)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$n$')
plt.ylabel(r'$\Delta \epsilon$')
plt.savefig('UncLogLog2.pdf')
plt.show()

# i=0
# while 10*i<N:
#     print(normalpsi[10*i])
#     i+=1
    
# print()

# i=0
# while 10*i<N:
#     print(x[10*i])
#     i+=1


