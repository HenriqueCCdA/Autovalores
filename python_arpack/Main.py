#!/bin/python3.6
import sys
import time as tm
#import networkx as nx
#import MyMatrix as MyM
from scipy import io
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import eigsh
#from scipy.sparse import coo_matrix

#......................................................................
def main(argv):
#checando os argumentos
    nArgs = len(argv)
    if nArgs < 3:
        sys.stderr.write("Usage: %s\nOpecoes:\n" % argv[0])
        print ("FileIn","FileOut")  
        return 1
#......................................................................
 
    fileIn  = argv[1]
    fileOut = argv[2]
# ... leitura do arquivo
    aCoo = io.mmread(fileIn)
    nl   = int(io.mminfo(fileIn)[0])
    nc   = int(io.mminfo(fileIn)[1])
    
#    aDense = aCoo.todense()
#   aDense = aCoo.todense()  
#   evalsAll, evecsAll = eigh(aDense)
#   print evalsAll
#   return
# convertendo COO para CSR
    timeConv = tm.time()
    aCsr   = aCoo.tocsr()
    timeConv = tm.time() - timeConv
#......................................................................

    maxit        = 3500
    tol          = 1.e-15
    lancozVector = 20
#Maior autovalor
    eighLarge = tm.time()
    evalsLarge,evecsLarge = eigsh(aCsr,1
                                 ,ncv=lancozVector
                                 ,which='LA'
                                 ,tol=tol
                                 ,maxiter=maxit)
    eighLarge = tm.time() - eighLarge 
    print (evalsLarge)
#    erro = abs(evalsLarge[0] - evalsAll[nl-1])/evalsAll[nl-1]
#    print evalsLarge[0],evalsAll[nl-1],erro
#......................................................................

#Menor autovalor                   
    eighSmall = tm.time()
    evalsSmall,evecsSmall = eigsh(aCsr,1 
                                 ,ncv=lancozVector
                                 ,which='LA'
                                 ,sigma=0
                                 ,tol=tol
                                 ,maxiter=maxit)
    eighSmall = tm.time() -  eighSmall
    print (evalsSmall)
#    erro = abs(evalsSmall[0] - evalsAll[0])/evalsAll[0]
#     print evalsSmall[0],evalsAll[0],erro
#......................................................................
    print("NUmero de equacoes %d" % nl)
    print("Maior autovalor %0.6e" % evalsLarge[0])
    print("Menor autovalor %0.6e" % evalsSmall[0])
    print("Condicio number %0.6e" % (evalsLarge[0]/evalsSmall[0]))
    print ("Time(Maior): ",eighLarge)
    print ("Time(Menor): ",eighSmall)
# .....................................................................    
    f = open(fileOut,"w")
    f.write("Numero de equacoes %d\n" % nl)
    f.write("Maior autovalor %.6e\n" % evalsLarge[0])
    f.write("Menor autovalor %.6e\n" % evalsSmall[0])
    f.write("Condicio number %.6e\n" %
                  (evalsLarge[0]/evalsSmall[0]))
    f.close()
#......................................................................
if __name__ == '__main__':
    sys.exit(main(sys.argv))
#......................................................................
