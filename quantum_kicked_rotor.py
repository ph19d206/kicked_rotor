#importing libraries
import numpy as np
import scipy as sp
import qutip as qt
import matplotlib.pyplot as plt

#--------------------------------------------
# input paramater
kicks=10
N=200
Nb=N*N
k=[9]# 1st one is k1. It is in list so that we can add k2 as next element.
b=2/N
pow_N=15
ϵ=10**(-pow_N)# tollerance. If value is bellow ϵ then it will taken as zero.
#--------------------------------------------

#converting to float
N=float(N)
Nb=float(Nb)
b=float(b)
k[0]=float(k[0])
#--------------------------------------------


# matrix elements are in position basis
# Flouet matrix element(m,n) of each unperturbed subsystem. kj input is k[j].
def U_k_elmt(kj,m,n): 
    part1=(np.pi/N)*pow(n-m,2)
    part2=(2*np.pi/N)*(n)
    part3=(kj*N/(2*np.pi))*np.cos(part2)
    uk_mn=(1/np.sqrt(N))*np.exp(-1j*part3)*np.exp(1j*part1)
    return uk_mn
#--------------------------------------------



# Floquet matrices

U_k1=sp.sparse.lil_matrix((int(N),int(N)),dtype=np.complex_)

# U_k1
for row in range(int(N)):
    for column in range(int(N)):
        U_k1_mn=U_k_elmt(k[0],float(row),float(column))
        if np.abs(U_k1_mn.real) >= ϵ or np.abs(U_k1_mn.imag) >= ϵ:
            U_k1[(row,column)]=np.round(U_k1_mn,pow_N)
            
                        
# As csr sparce or csc are efficient sparse matrix. So we change U_k1,U_k2,U_b to csr
U_k1 = U_k1.tocsr()


## Evolution of quantum state in Schrodinger Picture
#If the system is chaotic we can see the signature of chos as scrambling of hussimi distribution.
#--------------------------------------------
kicks=10
# Initial state.
#psi=qt.rand_ket(int(N)); ; state_name = 'Random state' #qt.rand_ket means any random initial state.
psi = qt.basis(int(N),0); state_name =  'fock basis' # fock basis
#psi = qt.coherent(int(N),1); state_name = 'coherent state' 
#--------------------------------------------

qt_U_k1 = qt.Qobj(U_k1.toarray())

# ploting the initial state
xvec = np.linspace(-5, 5, 500)
qt_Qfun = qt.qfunc(psi, xvec, xvec)

fig, axes = plt.subplots(1, 2, figsize=(12,3))
cont0 = axes[0].contourf(xvec, xvec, qt_Qfun, 100)
lbl0 = axes[0].set_title(f'{state_name}(initial)')

#time evolution of state
for i in range(kicks):
    psi = qt_U_k1*psi

qt_Qfun = qt.qfunc(psi, xvec, xvec)
cont1 = axes[1].contourf(xvec, xvec, qt_Qfun, 100)
lbl1 = axes[1].set_title(f'{state_name} (after {kicks} kicks)')

plt.show()
