#!/bin/python
import sys
import time as tm
import numpy as np
import networkx as nx
import MyMatrix as MyM
from scipy import io
from scipy.sparse import coo_matrix

#......................................................................
def main(argv):
    myArgv=(['-autovalor','-autovalorI','-condNumber'
            ,'-norm'     ,'-matvec'     ,'-svd'
            ,'-svdI'])
#checando os argumentos
    nArgs = len(argv)
    if nArgs < 4:
        sys.stderr.write("Usage: %s\nOpecoes:\n" % argv[0])
        for arg in myArgv:
            print arg,"FileIn","FileOut"  
        return 1
#......................................................................
 
    fileIn  = argv[2]
    fileOut = argv[3]
# ... leitura do arquivo
    aCoo = io.mmread(fileIn)
    nl = int(io.mminfo(fileIn)[0])
    nc = int(io.mminfo(fileIn)[1])
    a    = coo_matrix(aCoo,shape=(4,4)).todense()
#escolhendo a opcao de execucao

#autovalor
    if argv[1] == '-autovalor':
      timeIn  = tm.time()
# ... 
      MyM.autovalor(a,'false',fileOut)
# ...........................................................
      timeOut      = tm.time()
      print timeOut - timeIn
# ...........................................................

#autovalorI
    elif argv[1] == '-autovalorI':
      timeIn  = tm.time()
# ... 
      MyM.autovalor(a,'true',fileOut)
# ...........................................................
      timeOut      = tm.time()
      print timeOut - timeIn
    
#numero de condicionamento    
    elif  argv[1] == '-condNumber':
      timeIn       = tm.time()
# ... 
      MyM.condNumber(a,'true',fileOut)
# ...........................................................
      timeOut      = tm.time()
      print 'time:', timeOut - timeIn
# ...........................................................

#norma da matriz    
    elif  argv[1] == '-norm':
      timeIn  = tm.time()
# ... 
      MyM.norm(a,2,nl,'false',fileOut)
# ...........................................................
      timeOut      = tm.time()
      print timeOut - timeIn
#......................................................................

#matvec    
    elif  argv[1] == '-matvec':
      timeIn  = tm.time()
# ... 
      MyM.matVec(a,nl,'true',fileOut)
# ...........................................................
      timeOut      = tm.time()
      print timeOut - timeIn
#......................................................................

#svd       
    elif  argv[1] == '-svd':
      timeIn  = tm.time()
# ... 
      MyM.svd(a,'false',fileOut)
# ...........................................................
      timeOut      = tm.time()
      print timeOut - timeIn
#......................................................................

#svdI       
    elif  argv[1] == '-svdI':
      timeIn  = tm.time()
# ... 
      MyM.svd(a,'true',fileOut)
# ...........................................................
      timeOut      = tm.time()
      print 'time:', timeOut - timeIn
#......................................................................

    else:
      sys.stderr.write("Usage: %s\nOpecoes:\n" % argv[0])
      for arg in myArgv:
        print arg,"FileIn","FileOut"  
      return 1
#......................................................................
if __name__ == '__main__':
    sys.exit(main(sys.argv))
#......................................................................
