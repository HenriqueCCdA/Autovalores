#!/bin/python
from numpy import linalg as La
import numpy as np
import HccaBlas as MyBlas
import random as rd
import sys
import math
import time as tm

__all__ = ['autovalor','condNumber','norm','matVec']
#**********************************************************************
#* AUTOVALOR: calcula os autovetores da matriz a                      *
#* ------------------------------------------------------------------ *
#* Paramentros de entrada:                                            *
#* ------------------------------------------------------------------ *
#* a       -> matriz no formato a(n,n)                                *
#* flagI   -> calculo os autovalores da matriz inversa                *
#* fileOut -> nome do arquivo de saida                                *
#* ------------------------------------------------------------------ *
#* Paramentros de saida:                                              *
#* ------------------------------------------------------------------ *
#* OBS:                                                               *
#* Calcula os autovalores atraves do pacote numpy                     *
#**********************************************************************
def autovalor(a,iFlag,fileOut):
# ... matriz inversa
  if iFlag == 'true':
    fileAux = fileOut + '_autovalorI.txt'
    f = open(fileAux,"w")
    ein = La.eigvals(a.I)
    n = len(a)
    for nEin in range(n):
      aux  = ein[nEin]
      auxi = 1/aux
      f.write("%9i %16.8e %16.8ei %16.8e inv %16.8e %16.8ei\n"
      %(nEin+1,aux.real,aux.imag,abs(aux),auxi.real,auxi.imag))
# ... matriz
  else:
    fileAux = fileOut + '_autovalor.txt'
    f = open(fileAux,"w")
    ein = La.eigvals(a)
    ein = np.sort(ein)
    n = len(a)
    for nEin in range(n):
      aux  = ein[nEin]
      auxi = 1/aux
      f.write("%9i %16.8e %16.8ei %16.8e inv %16.8e %16.8ei\n"
      %(nEin+1,aux.real,aux.imag,abs(aux),auxi.real,auxi.imag))
  f.close()
#**********************************************************************

#**********************************************************************
#* CONDNUMBER : numero de condicionamento da matriz                   *
#* ------------------------------------------------------------------ *
#* Paramentros de entrada:                                            *
#* ------------------------------------------------------------------ *
#* a       -> matriz no formato a(n,n)                                *
#* fileOUt -> nome do arquivo de saida                                *
#* ------------------------------------------------------------------ *
#* Paramentros de saida:                                              *
#* ------------------------------------------------------------------ *
#* OBS:                                                               *
#* Calcula os autovalores atraves do pacote numpy                     *
#**********************************************************************
def condNumber(a,flag,fileOut):
# ... 2-norm 
  condNumber2= La.cond(a,2)
# ... 1-norm 
  condNumber1= La.cond(a,1)
  fileAux = fileOut + '_condNumber.txt'
  f = open(fileAux,"w")
  f.write('Numero de condicionamento (1-norm)  : %e\n'%condNumber1)
  f.write('Numero de condicionamento (2-norm)  : %e\n'%condNumber2)

# ... 2-norm (definicao ||A||||A-1||)
# ||A|| - max sing valor(AAT)         
  if flag == 'true':
# ||A||
    aaT   = a.dot(a.T)
    normA = math.sqrt(max(abs(La.eigvals(aaT)))) 
# ||A-1||
    Ia     = La.inv(a)
    IaaT   = Ia.dot(Ia.T)
    normAi = math.sqrt(max(abs(La.eigvals(IaaT)))) 
    f.write('Numero de condicionamento (2-norm-d): %e\n'%(normA*normAi))
# .....................................................................

#...
  f.close()
#**********************************************************************

#**********************************************************************
#* NORM : calculo da norma de um matriz                               *
#* ------------------------------------------------------------------ *
#* Paramentros de entrada:                                            *
#* ------------------------------------------------------------------ *
#* a       -> matriz no formato a(n,n)                                *
#* fileOUt -> nome do arquivo de saida                                *
#* ------------------------------------------------------------------ *
#* Paramentros de saida:                                              *
#* ------------------------------------------------------------------ *
#* OBS:                                                               *
#* Calcula os autovalores atraves do pacote numpy                     *
#**********************************************************************
def norm(a,p,nl,flag,fileOut):
# ... 2-norm 
  n = La.norm(a,p)
  fileAux = fileOut + '_norm.txt'
  f = open(fileAux,"w")
  f.write('%i-norm:             %f\n'%(p,n))

# ... 2-norm (definicao ||A||||A-1||)
# ||A|| - max sing valor(AAT)         
  if flag == 'true':
    normA = 0.0
    x = np.ones(nl)
    for sample in range(1000):
# ...
      for i in range(nl):
        x[i] = rd.uniform(-300000,300000)
# .....................................................................    
      y = a.dot(x)
      normAx = math.sqrt(y.dot(y.T))
      normX = math.sqrt(x.dot(x))
      normA = max(normA,normAx/normX)
# .....................................................................    
    f.write('%i-norm (definicao): %f\n'%(p,normA))
#...
  f.close()
#**********************************************************************

#**********************************************************************
#* MATVEC : operacoes de matriz vetor                                 *
#* ------------------------------------------------------------------ *
#* Paramentros de entrada:                                            *
#* ------------------------------------------------------------------ *
#* a       -> matriz no formato a(n,n)                                *
#* x       -> vetor de multiplicacao                                  *
#* ------------------------------------------------------------------ *
#* Paramentros de saida:                                              *
#* ------------------------------------------------------------------ *
#* OBS:                                                               *
#* Calcula os autovalores atraves do pacote numpy                     *
#**********************************************************************
def matVec(a,nl,flag,fileOut):
  x = np.ones(nl)
  MyBlas.intRandom(x)
# ...
  fileAux = fileOut + '_matVec.txt'
  f = open(fileAux,"w")
# .....................................................................  
  timeIn  = tm.time()
# ...
  y1 = a.dot(x)
# .....................................................................
  timeOut      = tm.time()
  print 'numPy:', timeOut - timeIn
# .....................................................................

# ... minha implementacao do matvec  
# ||A|| - max sing valor(AAT)         
  if flag == 'true':
# ...
    timeIn  = tm.time()
# .....................................................................

# ...
    y2 = np.ones(nl)
    MyBlas.matVecFull(a,x,y2)
# .....................................................................
    timeOut      = tm.time()
    print 'MyBlas', timeOut - timeIn

#.....................................................................  

# ...  
  for n in range(nl):
    value=' %i '+str(y1[0,n])+' '+str(y2[n])+'\n'
    f.write(value%(n+1))
#.....................................................................  
  f.close()
#**********************************************************************

#**********************************************************************
#* SVD : decompiscao em valor singulares                              *
#* ------------------------------------------------------------------ *
#* Paramentros de entrada:                                            *
#* ------------------------------------------------------------------ *
#* a       -> matriz no formato a(n,n)                                *
#* flagI   -> calculo os autovalores da matriz inversa                *
#* fileOUt -> nome do arquivo de saida                                *
#* ------------------------------------------------------------------ *
#* Paramentros de saida:                                              *
#* ------------------------------------------------------------------ *
#* OBS:                                                               *
#* Calcula os autovalores atraves do pacote numpy                     *
#**********************************************************************
def svd(a,iFlag,fileOut):
  if iFlag == 'true':
    Ia     = La.inv(a)
# ... svd    
    sing1= La.svd(Ia,compute_uv=0)
# ...
    aaT   = Ia.dot(Ia.T)
    sing2 = La.eigvals(aaT) 
# ... 
    fileAux = fileOut + '_svdI.txt'
    f = open(fileAux,"w")
# ... 
    n = len(sing1)
    f.write("#  n           svd        raiz(ro(AAt)) ro(AAt)\n")
# ...    
    for nEin in range(n):
      aux  = sing1[nEin]
      aux1 = sing2[nEin]
      aux2 = math.sqrt(sing2[nEin])
      f.write("%3i %16.8f %16.8fi %16.8f %16.8f %16.8fi\n"
      %(nEin,aux.real,aux.imag,aux2,aux1.real,aux1.imag))
  else:
# ... svd    
    sing1= La.svd(a,compute_uv=0)
# ...
    aaT   = a.dot(a.T)
    sing2 = La.eigvals(aaT) 
# ... 
    fileAux = fileOut + '_svd.txt'
    f = open(fileAux,"w")
# ... 
    n = len(sing1)
    f.write("#  n           svd        raiz(ro(AAt)) ro(AAt)\n")
    for nEin in range(n):
      aux  = sing1[nEin]
      aux1 = sing2[nEin]
      aux2 = math.sqrt(sing2[nEin])
      f.write("%3i %16.8f %16.8fi %16.8f %16.8f %16.8fi\n"
      %(nEin,aux.real,aux.imag,aux2,aux1.real,aux1.imag))
#...
  f.close()
#**********************************************************************
