c *********************************************************************
c * CON_GRAPH : conversao do grafo para o formato csr                 *
c *********************************************************************
c * Parametros de entrada:                                            *
c * ------------------------------------------------------------------*
c * viz(nhsared,*) - vizinhos de primeira ordem                       *
c * xadj           -                                                  *
c * adjncy         -                                                  *
c * numel          - numero de elementos                              *
c * nshared        - numero de aresta por elementos                   *
c * ------------------------------------------------------------------*
c * Parmetros de saida:                                               *
c * ------------------------------------------------------------------*
c * xadj    - ponteiro do para grafo ( formato csr )                  *
c * adjncy - arranjo do grafo (formato csr)                           *  
c * ------------------------------------------------------------------*
c *********************************************************************  
       subroutine conv_graph(viz,numel,nshared,xadj,adjncy)
       implicit none
       integer viz(nshared,*)
       integer xadj(*),adjncy(*)
       integer i,j,aux,numel,nshared,nelviz,ipont
c ...
       xadj(1) = 1   
       do i = 1, numel
         aux = 0
         ipont = xadj(i) 
         do j = 1, nshared
           nelviz  = viz(j,i)
           if(nelviz .ne. -1) then
             adjncy(ipont+aux) = nelviz
             aux = aux + 1
           endif
         enddo
         xadj(i+1) = xadj(i) + aux
       enddo
c .....................................................................
c
c ...
       return
       end
c **********************************************************************
c
c **********************************************************************
      subroutine sortgraphCsr(ia,ja,neq)
c **********************************************************************
c *                                                                    *
c *   SORTGRAPH: ordena o grafo armazenado em ja() em ordem crescente  *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *    ia(neq+1) - ia(i) informa a posicao no vetor a do primeiro      *
c *                      coeficiente nao-nulo da equacao i             *
c *    ja(nad) - ja(k) informa a coluna do coeficiente que ocupa       *
c *              a posicao k no vetor a                                *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *    ja(nad) - grafo ordenado                                        *
c *                                                                    *
c **********************************************************************
      implicit none
      integer ia(*),ja(*),neq,i
c ......................................................................               
      do 100 i = 1,neq
         call bubblesort(ja(ia(i)),ia(i+1)-ia(i))
  100 continue
c ......................................................................           
      return
      end
      subroutine bubblesort(ja,n)
c **********************************************************************
c *                                                                    *
c *   BUBBLESORT: ordena o arranjo ja(n), em ordem crescente.          *
c *                                                                    *
c **********************************************************************
      implicit none
      integer ja(*),n,i,j,itroca
c ......................................................................               
   50 continue
      itroca = 0
      do 100 i = 2, n
         if(ja(i) .lt. ja(i-1)) then
            j = ja(i-1)
            ja(i-1) = ja(i)
            ja(i) = j
            itroca = 1
         endif
  100 continue
      if (itroca .eq. 1) goto 50
c ......................................................................               
      return
      end
