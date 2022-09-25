#!/usr/bin/env python
# coding: utf-8

# # Lista 3 - RCO
# i am doing something here

# Autor: Norton Martin Reichert Trennepohl

# ## Questão 1

# In[127]:


import numpy as np
import math
import matplotlib.pyplot as plt


# Constantes dadas

# In[125]:


# Direções
direc = [math.radians(45),math.radians(0),math.radians(0),math.radians(45)]
n = len(direc)
#Vetor de carregamentos (elementos não-nulos estão em unidades de N/mm)
#carreg = [[Nx],[Ny],[Nxy],[Mx],[My],[Mxy]]
carreg = [[1000],[200],[0],[0],[0],[0]]

#Espessura da placa
h = 3E-3 #m

E11 = 19.76E9 #Pa
E22 = 1.97E9 #Pa
nu12 = 0.35
G12 = 0.7E9 #Pa

nu21 = (E22*nu12)/E11


# Matriz de rigidez

# In[83]:


Q11 = E11/(1-nu12*nu21)
Q22 = E22/(1-nu12*nu21)
Q66 = G12
Q12 = (nu12*E22)/(1-nu12*nu21)
Q21 = Q12

Q = np.array([[Q11, Q12, 0],[Q21, Q22, 0], [0, 0, Q66]])

cos = np.cos(direc)
sin = np.sin(direc)

# T = np.zeros((1,4))
# T_inv = np.zeros((1,4))
# Q_dash = np.zeros((1,4))

T = [[0],[0],[0],[0]]
T_inv = [[0],[0],[0],[0]]
Q_dash = [[0],[0],[0],[0]]

for i in range(n):
    T[i] = np.array([[cos[i]**2, sin[i]**2, 2*sin[i]*cos[i]],[sin[i]**2, cos[i]**2, -2*sin[i]*cos[i]],[-sin[i]*cos[i], sin[i]*cos[i], cos[i]**2-sin[i]**2]])
    T_inv[i] = np.linalg.inv(T[i])
    Q_dash[i] = T_inv[i]@Q@T[i]


# In[84]:


z0 = n*h/2
A_local = [[0],[0],[0],[0]]
A_global = 0
for i in range(n):
    A_local[i] = Q_dash[i]*(((((n/2)-(i+1))/n)*-h) - (((n/2 - i)/n)*-h))
    #A_local[i] = Q_dash[i]*(h/n)
    A_global = A_global + A_local[i]
    
#print(A_global)
    
B_local = [[0],[0],[0],[0]]
B_global = 0
for i in range(n):
    B_local[i] = 0.5*Q_dash[i]*(((((n/2)-(i+1))/n)*-h)**2 - (((n/2 - i)/n)*-h)**2)
    #print(B_local[i])
    B_global = B_global + B_local[i]
    #print(B_global)
#print(B_global)


D_local = [0,0,0,0]
D_global = 0
for i in range(n):
    D_local[i] = (1/3)*Q_dash[i]*(((((n/2)-(i+1))/n)*-h)**3 - (((n/2 - i)/n)*-h)**3)
    #print(D_local[i])
    D_global = D_global + D_local[i]
#print(D_global)
    
linha1 = np.vstack((A_global,B_global))
linha2 = np.vstack((B_global,D_global))
ABBD = np.hstack((linha1,linha2))
#print(ABBD)


# In[139]:


# Passo 5
def_curv = [[0],[0],[0],[0],[0],[0]]

def_curv = np.linalg.inv(ABBD)@carreg

epsilon_0_global = np.vstack((def_curv[0],def_curv[1],def_curv[2]))
K_global = np.vstack((def_curv[3],def_curv[4],def_curv[5]))

#print(epsilon_0_global)
#print(K_global)

# Passo 6 (coordenadas z referenciadas no plano médio de cada lâmina)
sigma_global = [[0],[0],[0],[0]]
sigma_local = [[0],[0],[0],[0]]
z = [[0],[0],[0],[0]]
epsilon_global = [[0],[0],[0],[0]]
epsilon_local = [[0],[0],[0],[0]]
y = [[0],[0],[0],[0]]

for i in range(n):
    z[i] = 0.5*(((((n/2)-(i+1))/n)*-h) + (((n/2 - i)/n)*-h))
    #print(z[i])
    sigma_global[i] = (Q_dash[i]@(epsilon_0_global + (z[i]*K_global)))*10**(-3)
    print("==============================")
    print("Resultados lâmina %d:" %(i+1))
    print("Tensão na lâmina no sistema global de coordenadas (kPa): ")
    print(sigma_global[i])
    sigma_local[i] = T[i]@sigma_global[i]
    print("Tensão na lâmina no sistema local de coordenadas (kPa):")
    print(sigma_local[i])
    epsilon_global[i] = epsilon_0_global + z[i]*K_global
    print("Deformação no plano médio da lâmina no sistema global de coordenadas:")
    print(epsilon_global[i])
    epsilon_local[i] = T[i]@epsilon_global[i]
    print("Deformação no plano médio da lâmina no sistema local de coordenadas:")
    print(epsilon_local[i])
    


# In[ ]:





# In[ ]:


### Teste se o resultado era o mesmo pra ABBD

# def_curv_2 = [[0],[0],[0],[0],[0],[0]]

# A_est = np.linalg.inv(A_global)
# B_est = -1*np.linalg.inv(A_global)@B_global
# C_est = B_global@np.linalg.inv(A_global)
# D_est = D_global-B_global@np.linalg.inv(A_global)@B_global

# A_ap = A_est-B_est@np.linalg.inv(D_est)@C_est
# B_ap = B_est@np.linalg.inv(D_est)
# C_ap = B_ap
# D_ap = np.linalg.inv(D_est)

# print(B_ap)

# linha1_2 = np.vstack((A_ap,B_ap))
# linha2_2 = np.vstack((C_ap,D_ap))
# ABBD_2 = np.hstack((linha1_2,linha2_2))

# def_curv_2 = ABBD_2@carreg
# print(def_curv)
# print(def_curv_2)




# y[i] = epsilon_local[i][1]
    
# plt.figure()
# plt.plot(z, y)

