#!/bin/python
import random

__all__ = ['matVecFull','intRandom']

#**********************************************************************
#* MATVEC: calcula o produto matriz vetor                             *
#* ------------------------------------------------------------------ *
#* Paramentros de entrada:                                            *
#* ------------------------------------------------------------------ *
#* a       -> matriz no formato a(n,n)                                *
#* x       -> vetor no formato x(n)                                   *
#* ------------------------------------------------------------------ *
#* Paramentros de saida:                                              *
#* ------------------------------------------------------------------ *
#* y       -> vetor da multiplicacao y(n)                             *
#* ------------------------------------------------------------------ *
#* OBS:                                                               *
#* Calcula os autovalores atraves do pacote numpy                     *
#**********************************************************************
def matVecFull(a,x,y):
  n = len(x)
  for nl in range(n):
    y[nl] = 0.0
    for cl in range(n):
      y[nl] += a[nl,cl]*x[cl]
#**********************************************************************

#**********************************************************************
#* INTRANDOM: inicializacao de um vetor random                        *
#* ------------------------------------------------------------------ *
#* Paramentros de entrada:                                            *
#* ------------------------------------------------------------------ *
#* a       -> matriz no formato a(n,n)                                *
#* x       -> vetor no formato x(n)                                   *
#* ------------------------------------------------------------------ *
#* Paramentros de saida:                                              *
#* ------------------------------------------------------------------ *
#* y       -> vetor da multiplicacao y(n)                             *
#* ------------------------------------------------------------------ *
#* OBS:                                                               *
#* Calcula os autovalores atraves do pacote numpy                     *
#**********************************************************************
def intRandom(x):
  nl = len(x)
  for i in range(nl):
    x[i] = random.uniform(-1000,1000)
#**********************************************************************

