      subroutine matvec_csrcsym(neq,ia,ja,al,ad,x,y)
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRCSYM: produto matriz-vetor y = Ax  (A simetrica),      *
c *                   coef. de A no formato CSRC.                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da linha   i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor au                              *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de A, no formato CSR, ou      *
c *            parte triangular superior de A, no formato CSC          *
c *   au(*)  - nao utilizado                                           *
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      implicit none
      integer neq,ia(*),ja(*),i,k,jak
      real*8  ad(*),al(*),x(*),y(*),s,t,xi
c ... Loop nas linhas:
c
      do 110 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         do 100 k = ia(i), ia(i+1)-1
            jak = ja(k)
            s   = al(k)
c
c ...       Produto da linha i pelo vetor x (parte triangular inferior):
c
            t   = t + s*x(jak)
c
c ...       Produto dos coef. da parte triangular superior por x(i):
c
            y(jak) = y(jak) + s*xi
  100    continue
c
c ...    Armazena o resultado em y(i):
c
         y(i) = t
  110 continue
c ......................................................................
      return
      end
      subroutine matvec_csrcsym1(neq,ia,ja,al,ad,x,y)
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRCSYM1: produto matriz-vetor y = Ax  (A simetrica),     *
c *                    coef. de A no formato CSRC.                     *
c *                    (versao com loop interno desenrolado)           *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da linha   i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor au                              *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de A, no formato CSR, ou      *
c *            parte triangular superior de A, no formato CSC          *
c *   au(*)  - nao utilizado                                           *
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      implicit none
      integer neq,ia(*),ja(*),i,k,jak
      real*8  al(*),ad(*),x(*),y(*),t,xi,s,s1
      integer k1,k2,n,jak1
c ... Loop nas linhas:
c
c ----------------------------------------------------------------------
      do 300 i = 1, neq
c
c ...    Produto da diagonal de A por x:
c
         xi = x(i)
         t  = ad(i)*xi
         k1 = ia(i)
         k2 = ia(i+1)
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         if (k2 .eq. k1) goto 120
         n = mod(k2-k1,2)
         if (n .eq. 0) goto 100
         jak1 = ja(k1)
         t = t + al(k1)*x(jak1)
         y(jak1) = y(jak1) + al(k1)*xi
         k1 = k1 + 1
  100    do 110 k = k1, k2-1, 2
            jak  = ja(k)
            jak1 = ja(k+1)
            s    = al(k)
            s1   = al(k+1)
            t    = t + s*x(jak) + s1*x(jak1)
            y(jak)  = y(jak)  +  s*xi
            y(jak1) = y(jak1) + s1*xi
  110    continue
c
c ...    Armazena o resultado em y(i):
c  
  120    y(i) = t
c ......................................................................
  300 continue
      return
      end
c **********************************************************************
      subroutine matvec_csrcsym2(neq,ia,ja,al,ad,x,y)
c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRCSYM2: produto matriz-vetor y = Ax  (A simetrica),     *
c *                    coef. de A no formato CSRC.                     *
c *                    (versao com loops desenrolados)                 *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq   - numero de equacoes                                       *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro      *
c *                     coeficiente nao-nulo da linha   i              *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa      *
c *               a posicao k no vetor au                              *
c *   ad(neq)- diagonal da matriz A                                    *
c *   al(nad)- parte triangular inferior de A, no formato CSR, ou      *
c *            parte triangular superior de A, no formato CSC          *
c *   au(*)  - nao utilizado                                           *
c *   x(neq) - vetor a ser multiplicado                                *
c *   y(neq) - nao definido                                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      implicit none
      integer neq,ia(*),ja(*),i,k,jak
      real*8  ad(*),al(*),x(*),y(*),t,xi,t1,xi1,s,s1
      integer k1,k2,n,jak1,jak2,i0,k3
c ... Loop nas linhas:
c
      i0 = 1
      n = mod(neq,2)
      if (n .eq. 0) goto 10
      y(1) = ad(1)*x(1)
      i0 = 2
c -------------------------------------------------------------      
   10 do 300 i = i0, neq, 2
c
c ...    Produto da diagonal de A por x:
c
         xi  = x(i)
         xi1 = x(i+1)
         t   =   ad(i)*xi
         t1  = ad(i+1)*xi1                    
         k1 = ia(i)
         k2 = ia(i+1)
         k3 = ia(i+2)
c
c ...    Loop nos coeficientes nao nulos da linha i:
c
         if (k2 .eq. k1) goto 120
         n = mod(k2-k1,2)
         if (n .eq. 0) goto 100
         jak = ja(k1)         
         t = t + al(k1)*x(jak)
         y(jak) = y(jak) + al(k1)*xi
         k1 = k1 + 1
         
  100    do 110 k = k1, k2-1, 2
            jak1 = ja(k)
            jak2 = ja(k+1)
            s =  al(k)
            s1 = al(k+1)
            t = t + s*x(jak1) + s1*x(jak2)
            y(jak1) = y(jak1) + s*xi
            y(jak2) = y(jak2) + s1*xi                        
  110    continue
c
c ...    Armazena o resultado em y(i):
c  
  120    y(i) = t
c -------------------------------------------------------------
c
c ...    Loop nos coeficientes nao nulos da linha i+1:
c
         if (k3 .eq. k2) goto 220
         n = mod(k3-k2,2)
         if (n .eq. 0) goto 200
         jak = ja(k2)         
         t1 = t1 + al(k2)*x(jak)
         y(jak) = y(jak) + al(k2)*xi1
         k2 = k2 + 1
  200    do 210 k = k2, k3-1, 2
            jak1 = ja(k)
            jak2 = ja(k+1)
            s = al(k)
            s1 = al(k+1)            
            t1 = t1 + s*x(jak1) + s1*x(jak2)
            y(jak1) = y(jak1) +  s*xi1
            y(jak2) = y(jak2) + s1*xi1                        
  210    continue  
c
c ...    Armazena o resultado em y(i):
c
  220 y(i+1) = t1
c ----------------------------------------------------------------------  
  300 continue
c ......................................................................
      return
      end
      real*8 function dot(a,b,n)
c **********************************************************************
c *                                                                    *
c *   DOT: Produto escalar a.b                                         *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8  a(*),b(*),tmp
c ......................................................................
      tmp = 0.0d+0
      do 100 i = 1, n
         tmp = tmp + a(i)*b(i)
  100 continue
      dot = tmp
c ......................................................................    
      return
      end
      real*8 function ddot(x,y,n)
c **********************************************************************
c *                                                                    *
c *   DOT: Produto escalar x.y                                         *
c *                                                                    *
c **********************************************************************
      implicit none    
      real*8 x(*),y(*),temp
      integer i,m,mp1,n,inc
c ......................................................................
      ddot = 0.d0
      temp = 0.d0
      inc  = 10
c ......................................................................
      m = mod(n,inc)
      if(m .eq. 0) goto 40
      do 30 i = 1,m
         temp = temp + x(i)*y(i)
   30 continue
      if(n .lt. inc) goto 60
   40 mp1 = m + 1
      do 50 i = mp1, n, inc
         temp = temp + x(i  )*y(  i) + x(i+1)*y(i+1) + x(i+2)*y(i+2) +
     .                 x(i+3)*y(i+3) + x(i+4)*y(i+4) + x(i+5)*y(i+5) +
     .                 x(i+6)*y(i+6) + x(i+7)*y(i+7) + x(i+8)*y(i+8) +
     .                 x(i+9)*y(i+9)
   50 continue
   60 ddot = temp
      return
      end
      subroutine aequalb(a,b,n)
c **********************************************************************
c *                                                                    *
c *   AEQUALB:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    b(n) - vetor de dimensao n                                      *
c *    n    - dimensao                                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   a = b                                                            *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8 a(*),b(*)
c ......................................................................
      do 100 i = 1, n
         a(i) = b(i)
  100 continue
      return
      end
      subroutine aminusb(a,b,c,n)
c **********************************************************************
c *                                                                    *
c *   AMINUSB:                                                         *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    b(n) - vetor de dimensao n                                      *
c *    n    - dimensao                                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   c = a - b                                                        *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8 a(*),b(*),c(*)
c ......................................................................
      do 100 i = 1, n
         c(i) = a(i) - b(i)
  100 continue
      return
      end
      subroutine vsum(a,b,n,c)
c **********************************************************************
c *                                                                    *
c *   VSUM:                                                            *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    b(n) - vetor de dimensao n                                      *
c *    n    - dimensao                                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   c = a + b                                                        *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8 a(*),b(*),c(*)
c ......................................................................
      do 100 i = 1, n
         c(i) = a(i) + b(i)
  100 continue
      return
      end      
      subroutine vsmul(a,x,n,c)
c **********************************************************************
c *                                                                    *
c *   VSUM:                                                            *
c *   -------                                                          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *   ---------------------                                            *
c *                                                                    *
c *    a(n) - vetor de dimensao n                                      *
c *    x    - escalar                                                  *
c *    n    - dimensao                                                 *
c *                                                                    *
c *   Parametros de saida:                                             *
c *   -------------------                                              *
c *                                                                    *
c *   c = x*a                                                          *
c *                                                                    *
c **********************************************************************
      implicit none
      integer n,i
      real*8 a(*),x,c(*)
c ......................................................................
      do 100 i = 1, n
         c(i) = a(i)*x
  100 continue
      return
      end        
