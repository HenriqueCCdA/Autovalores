      real*8 function dot(a,b,n)
c **********************************************************************
c *                                                                    *
c *   DOT: Produto escalar a.b                                         *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8 a(*),b(*),dotOmp
c ......................................................................
      dotOmp = 0.0d0
c$omp parallel do reduction(+:dotOmp)
      do  i = 1, n
        dotOmp = dotOmp + a(i)*b(i)
      enddo
c$omp end parallel do
      dot = dotOmp
      return
      end
      subroutine diag(neq,ad,x)
c **********************************************************************
c *                                                                    *
c *   DIAG                                                             *
c *                                                                    *
c **********************************************************************
      integer neq
      real*8 ad(*),x(*)
      do 100 i = 1, neq
       x(i) = x(i)/ad(i)
  100 continue
      return
      end
      subroutine saxpb(a,b,x,n,c)
c **********************************************************************
c *                                                                    *
c *   SAXPB: escalar . vetor + vetor                                   *
c *   -----                                                            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    b(n) - vetor de dimensao n                                      *
c *    x    - escalar                                                  *
c *    n    - dimensao                                                 *
c *    c(n) - nao definido                                             *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   c - resultado: a.x + b                                           *
c *                                                                    *
c **********************************************************************
      integer n
      real*8 a(*),b(*),c(*),x
c ......................................................................
      do 100 i = 1, n
         c(i) = a(i) * x + b(i)
  100 continue
      return
      end
      subroutine azero(a,j)
c **********************************************************************
c *                                                                    *
c *   AZERO: zera as posicoes de 1 ate j do vetor a.                   *
c *                                                                    *
c **********************************************************************
      real*8 a(*)
c ......................................................................
      do 100 k = 1, j
         a(k) = 0.d0
  100 continue
      return
      end
      subroutine mzero(m,j)
c **********************************************************************
c *                                                                    *
c *   MZERO: zera as posicoes de 1 ate j do vetor m.                   *
c *                                                                    *
c **********************************************************************
      integer m(*)
c ......................................................................
      do 100 k = 1, j
         m(k) = 0
  100 continue
      return
      end
      subroutine acopy(a,i1,i2,i3)
c **********************************************************************
c *                                                                    *
c *   ACOPY: desloca um vetor de reais na memoria.                     *
c *                                                                    *
c **********************************************************************
      real*8 a(*)
c ...........................................
      if (i3 .gt. i1) goto 200
      do 100 i = 1, i2-i1
         a(i+i3-i1) = a(i)
  100 continue
      return
  200 continue
      do 300 i = i2-i1, 1, -1
         a(i+i3-i1) = a(i)
  300 continue
      return
      end
      subroutine mcopy(m,i1,i2,i3)
c **********************************************************************
c *                                                                    *
c *   MCOPY: desloca um vetor de inteiros na memoria.                  *
c *                                                                    *
c **********************************************************************
      integer m(*)
c ...........................................
      if (i3 .gt. i1) goto 200
      do 100 i = 1, i2-i1
         m(i+i3-i1) = m(i)
  100 continue
      return
  200 continue
      do 300 i = i2-i1, 1, -1
         m(i+i3-i1) = m(i)
  300 continue
      return
      end
      subroutine pmove(a,b,nn)
c***********************************************************************
c*                                                                     *
c*     PMOVE: Moves real array a into b                                *
c*                                                                     *
c*     Inputs:                                                         *
c*        a(*)      - Array to move                                    *
c*        nn        - Length of array to move                          *
c*                                                                     *
c*     Outputs:                                                        *
c*        b(*)      - Moved array                                      *
c*                                                                     *
c***********************************************************************
      implicit  none
      integer   nn,n
      real*8    a(nn),b(nn)
c.......................................................................
      do 100 n = 1,nn
        b(n) = a(n)
  100 continue
      return
      end
      subroutine avalue(a,j,value)
c **********************************************************************
c *                                                                    *
c *   AVALUE:                                                          *
c *                                                                    *
c **********************************************************************
      real*8 a(*),x,value
c ......................................................................
      do 100 k = 1, j
         a(k) = value
  100 continue
      return
      end
