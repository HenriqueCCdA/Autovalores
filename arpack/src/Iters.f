c*****************************Svn***************************************      
c*$Date: 2011-10-07 10:36:21 -0300 (Fri, 07 Oct 2011) $                 
c*$Rev: 951 $                                                           
c*$Author: henrique $                                                   
c***********************************************************************      
      subroutine pcg(neq,ia,ja,ad,au,al,m,b,x,z,r,tol,maxit,
     .               matvec,dot,my_id,neqf1i,neqf2i,neq_doti,i_fmapi,
     .               i_xfi,i_rcvsi,i_dspli)
c **********************************************************************
c *                                                                    *
c *   Subroutine PCG                                                   *
c *                                                                    *
c *   Solucao de sistemas de equacoes pelo metodo dos gradientes       *
c *   conjugados com precondicionador diagonal para matrizes           *
c *   simetricas.                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq    - numero de equacoes                                      *
c *   ia(*)  - ponteiro do formato CSR                                 *
c *   ja(*)  - ponteiro das colunas no formato CSR                     *
c *   ad(neq)- diagonal da matriz A                                    *
c *   au(*)  - parte triangular superior de A                          *
c *   al(*)  - parte triangular inferior de A                          *
c *   m(*)   - precondicionador diagonal                               *
c *   b(neq) - vetor de forcas                                         *
c *   x(neq) - chute inicial                                           *
c *   z(neq) - arranjo local de trabalho                               *
c *   r(neq) - arranjo local de trabalho                               *
c *   tol    - tolerancia de convergencia                              *
c *   maxit  - numero maximo de iteracoes                              *
c *   matvec - nome da funcao externa para o produto matrix-vetor      *
c *   dot    - nome da funcao externa para o produto escalar           *
c *   energy - nao definido                                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   x(neq) - vetor solucao                                           *
c *   b(neq) - modificado                                              *
c *   ad(*),al(*),au(*) - inalterados                                  *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,maxit,i,j,nad,k
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),m(*),x(*),r(*),z(*),b(*)
      real*8  dot,ddot,tol,conv,energy,d,alpha,beta
      real*8  time0,time
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
      nad = ia(neq+1)-1
      if(my_id.eq.0)print *, 'nad :',nad
c ......................................................................
c
c ... Chute inicial:
c
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c ----------------------------------------------------------------------
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
      do 100 i = 1, neq
         r(i) = b(i) - z(i)
         z(i) = r(i) / m(i)
         b(i) = z(i)
  100 continue
      d    = dot(r(1),z(1),neq_doti)
      conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
      k = 0
      do 230 j = 1, maxit
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               b,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
         alpha = d / dot(b,z,neq_doti)
         do 210 i = 1, neq
            x(i) = x(i) + alpha * b(i)
            r(i) = r(i) - alpha * z(i)
            z(i) = r(i) / m(i)
  210    continue
         beta = dot(r,z,neq_doti) / d
         do 220 i = 1, neq
            b(i) = z(i) + beta * b(i)
  220    continue
         d = beta * d
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ...
         k = k + 1 
         if(k .eq. 10000) then 
           k = 0 
           if(my_id.eq.0) write(*,1300),j,dsqrt(dabs(d)),conv
         endif
c ----------------------------------------------------------------------
  230 continue
c ----------------------------------------------------------------------
      write(*,1200) maxit
      call stop_mef()
  300 continue
c
c ... Energy norm:
c
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
      energy = dot(x(1),z(1),neq_doti)
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id.eq.0)write(*,1100)tol,neq,j,energy,time
c ......................................................................
c     Controle de flops
      if(my_id.eq.0) write(10,'(a,a,i9,a,d20.10,a,d20.10)')"PCG: ",
     .               "it",j, " energy norm ",energy," tol ",tol
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA PCG:',/,5x,'Coeficiente da diagonal nulo
     .ou negativo - equacao ',i7)
 1100 format(' (PCG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
 1300 format (' PCG it =',i9,' norm =',d20.10,' conv =',d20.10)
      end
      subroutine gmres(neq,ia,ja,ad,au,al,m,b,x,k,g,h,y,c,s,e,
     .              tol,maxit,matvec,dot,neqovlp,my_id,neqf1i,neqf2i,
     .              neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli)
c **********************************************************************
c *                                                                    *
c *   GMRES: Solucao iterativa de sistemas simetricos e nao-simetricos *
c *          pelo metodo GMRES com precondicionador diagonal.          *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq    - numero de equacoes                                      *
c *   ia(*)  - ponteiro do formato CSR                                 *
c *   ja(*)  - ponteiro das colunas no formato CSR                     *
c *   ad(neq)- diagonal da matriz A                                    *
c *   au(*)  - parte triangular superior de A                          *
c *   al(*)  - parte triangular inferior de A                          *
c *   m(*)   - precondicionador diagonal                               *
c *   b(neq) - vetor de forcas                                         *
c *   x(neq) - chute inicial                                           *
c *   z(neq) - arranjo local de trabalho                               *
c *   r(neq) - arranjo local de trabalho                               *
c *   tol    - tolerancia de convergencia                              *
c *   maxit  - numero maximo de iteracoes                              *
c *   matvec - nome da funcao externa para o produto matrix-vetor      *
c *   dot    - nome da funcao externa para o produto escalar           *
c *                                                                    *
c *   Arranjos locais de trabalho:                                     *
c *                                                                    *
c *      g(neq+1,k+1)                                                  *
c *      h(k+1,k)                                                      *
c *      y(k)                                                          *
c *      c(k)                                                          *
c *      s(k)                                                          *
c *      e(k+1)                                                        *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   x(neq) - vetor solucao                                           *
c *   b(neq),ad(*),al(*),au(*) - modificados                           *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................      
      integer neq,k,maxit,ia(*),ja(*),neqovlp,nit,i,j,l,ni,ic,nad,nad1
      real*8  ad(neq),au(*),al(*),m(*),b(*),x(*)
c      real*8  g(0:neqovlp,1:k+1),h(k+1,k),y(k),c(k),s(k),e(k+1),tol
      real*8  g(neqovlp,1:k+1),h(k+1,k),y(k),c(k),s(k),e(k+1),tol
      real*8  energy,econv,norm,dot,ddot,r,aux1,aux2,beta
      real*8  time0,time
      external matvec,dot
      integer my_id
c      integer nii(maxit)
c ......................................................................
      time0 = MPI_Wtime()
c ......................................................................
      nad = ia(neq+1)-1
      if(my_id.eq.0)print*,nad
c ----------------------------------------------------------------------
c
c.... Chute inicial:
c
      do 10 i = 1, neq
         x(i) = 0.d0
c ...    pre-condicionador diagonal:                  
         g(i,1) = b(i)/m(i)
   10 continue
c ----------------------------------------------------------------------
c
c ... Limite de convergencia:
c
      norm  = dsqrt(dot(g(1,1),g(1,1),neq_doti))
      econv = tol*norm
c ----------------------------------------------------------------------      
c
c ... Ciclos GMRES:
c
      nit = 0
      do 1000 l = 1, maxit
c
c ...... Residuo g(1) = b - A x:
c
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               x,g(1,1),neqf1i,neqf2i,i_fmapi,i_xfi,
     .               i_rcvsi,i_dspli)
c
c ...... Residuo com precondicionador diagonal:
c
         do 200 i = 1, neq
            g(i,1) = (b(i) - g(i,1))/m(i)
  200    continue
c
c ...... Norma do residuo:
c
         e(1) = dsqrt(dot(g(1,1),g(1,1),neq_doti))
c
c ...... Normalizacao de g1:
c
         do 210 i = 1, neq
            g(i,1) = g(i,1)/e(1)
  210    continue
c
c ...... Iteracoes GMRES:
c
         ni = 0
         do 400 i = 1, k
            nit = nit + 1
            ni  = ni  + 1
c
c ......... Produto g(i+1) = A.g(i):
c
            call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,
     .                  al(nad+1),g(1,i),g(1,i+1),neqf1i,neqf2i,
     .                  i_fmapi,i_xfi,i_rcvsi,i_dspli)
c
c ......... Precondicionador diagonal:
c
            do 300 j = 1, neq
                g(j,i+1) = g(j,i+1)/m(j)
  300       continue
c
c ......... Ortogonalizacao (Gram-Schmidt modificado):
c
            do 320 j = 1, i
               beta = dot(g(1,i+1),g(1,j),neq_doti)
               do 310 ic = 1, neq
                  g(ic,i+1) = g(ic,i+1) - beta * g(ic,j)
  310          continue
               h(j,i) = beta
  320       continue
c
c ......... Norma de g(i+1):
c
            norm = dsqrt(dot(g(1,i+1),g(1,i+1),neq_doti))
c
            h(i+1,i) = norm
c
c ......... Normalizacao de g(i+1):
c
            do 330 ic = 1, neq
               g(ic,i+1) = g(ic,i+1)/norm
  330       continue
c
            do 340 j = 1, i-1
               aux1 =  c(j) * h(j,i) + s(j) * h(j+1,i)
               aux2 = -s(j) * h(j,i) + c(j) * h(j+1,i)
               h(j,i)   = aux1
               h(j+1,i) = aux2
  340       continue
            r = dsqrt(h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i))
            c(i) = h(i,i)/r
            s(i) = h(i+1,i)/r
            h(i,i)   = r
            h(i+1,i) = 0.d0
            e(i+1) = -s(i) * e(i)
            e(i)   =  c(i) * e(i)
            if (dabs(e(i+1)) .le. econv) goto 500
  400    continue
  500    continue
c
c ...... Resolve o sistema h y = e :
c
         y(ni) = e(ni) / h(ni,ni)
         do 520 i = ni-1, 1, -1
            y(i) = 0.d0
            do 510 j = i+1, ni
               y(i) = y(i) - h(i,j)*y(j)
  510       continue
            y(i) = (y(i) + e(i)) / h(i,i)
  520    continue
c
c ...... Atualizacao de x:
c
         do 610 i = 1, neq
            do 600 j = 1, ni
               x(i) = x(i) + y(j) * g(i,j)
  600       continue
  610    continue
c
c ...... Verifica a convergencia:
c
c         nii(l)=ni
         if (dabs(e(ni+1)) .le. econv) goto 1100
 1000 continue
c ----------------------------------------------------------------------
 1100 continue
c
c ... Norma de energia da solucao:
c
      energy = dot(x(1),b(1),neq_doti)
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id.eq.0)write(*,2000) tol,neq,l,nit,dabs(e(ni+1)),energy
     .                           ,time
      if (dabs(e(ni+1)) .gt. econv) then
         write(*,2100) maxit
         call stop_mef()
      endif
c ......................................................................
c     Controle de flops
c      if(my_id.eq.0)write(10,'(999(i4,1x))') l,nit,(nii(j),j=1,l)
c ......................................................................
      return
c ----------------------------------------------------------------------
 2000 format(' (GMRES) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of cycles     = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Norm                 = ',d20.6/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 2100 format(' *** WARNING: no convergence reached for ',i9,' cycles !',
     . /)
      end      
      subroutine bicgstab(neq,ia,ja,ad,au,al,m,b,x,y,z,p,r,s,tol,maxit,
     .                    matvec,dot,my_id,neqf1i,neqf2i,neq_doti,
     .                    i_fmapi,i_xfi,i_rcvsi,i_dspli)
c **********************************************************************
c *                                                                    *
c *   Subroutine BICGSTAB                                              *
c *                                                                    *
c *   Solucao de sistemas de equacoes pelo metodo dos gradientes       *
c *   conjugados com precondicionador diagonal para matrizes           *
c *   simetricas.                                                      *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq    - numero de equacoes                                      *
c *   ia(*)  - ponteiro do formato CSR                                 *
c *   ja(*)  - ponteiro das colunas no formato CSR                     *
c *   ad(neq)- diagonal da matriz A                                    *
c *   au(*)  - parte triangular superior de A                          *
c *   al(*)  - parte triangular inferior de A                          *
c *   m(*)   - precondicionador diagonal                               *
c *   b(neq) - vetor de forcas                                         *
c *   x(neq) - chute inicial                                           *
c *   z(neq) - arranjo local de trabalho                               *
c *   r(neq) - arranjo local de trabalho                               *
c *   tol    - tolerancia de convergencia                              *
c *   maxit  - numero maximo de iteracoes                              *
c *   matvec - nome da funcao externa para o produto matrix-vetor      *
c *   dot    - nome da funcao externa para o produto escalar           *
c *   energy - nao definido                                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   x(neq) - vetor solucao                                           *
c *   b(neq) - modificado                                              *
c *   ad(*),al(*),au(*) - modificado                                   *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................
      integer neq,maxit,nad,i,j,k
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),m(*),x(*),r(*),p(*),b(*),y(*),z(*),s(*)
      real*8  dot,ddot,tol,conv,energy,d,alpha,beta,rr0,w
      real*8  time0,time
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
      nad = ia(neq+1)-1
      if(my_id.eq.0)print *, 'nad :',nad
c ......................................................................
c
c ... Chute inicial:
c
      do 10 i = 1, neq
         x(i) = 0.d0
c ...... pre-condicionador diagonal:         
         b(i)  = b(i)/m(i)
         ad(i) = ad(i)/m(i)
         do 5 k = ia(i), ia(i+1)-1
            j = ja(k)
            al(k) = al(k) / m(i)
            au(k) = au(k) / m(j)
   5     continue      
  10  continue
c ----------------------------------------------------------------------
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,r,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)      
      do 100 i = 1, neq
         r(i) = b(i) - r(i)
         p(i) = r(i)
         b(i) = r(i)
  100 continue
      d    = dot(r(1),r(1),neq)
      conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
      do 230 j = 1, maxit
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               p,z,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
         rr0 = dot(r,b,neq)
         alpha = rr0/dot(z,b,neq)
          do 210 i = 1, neq
            s(i) = r(i) - alpha * z(i)
  210    continue
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               s,y, neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
         w = dot(y,s,neq) / dot(y,y,neq)
         do 220 i = 1, neq
            x(i) = x(i) + alpha*p(i) + w * s(i)
            r(i) = s(i) - w*y(i)
  220    continue
         beta = (dot(r,b,neq) / rr0)*(alpha/w)
         do 225 i = 1, neq
             p(i) = r(i) + beta*(p(i)-w*z(i))
  225    continue
         d = dot(r,r,neq)  
         if (dsqrt(dabs(d)) .lt. conv) goto 300
  230 continue
c ----------------------------------------------------------------------
      write(*,1200) maxit
      call stop_mef()
  300 continue
c
c ... Energy norm:
c
c      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
c     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
      b(1:neq) = b(1:neq)*m(1:neq)
      energy   = dot(x(1),b(1),neq)
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id.eq.0)write(*,1100) neq,j,energy,time
c ......................................................................
c     Controle de flops
      if(my_id.eq.0)write(10,'(a,a,i9,a,d20.10,a,d20.10)')"BICGSTAB: ",
     .              "it",j, " energy norm ",energy," tol ",tol
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA BICGSTAB:',/,5x,'Coeficiente da diagonal
     . nulo ou negativo - equacao ',i7)
 1100 format(' (BICGSTAB) solver:'/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Energy norm          = ',d20.6/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i4,
     .        ' iterations !',/)
      end
c **********************************************************************
      subroutine pbicgstab(neq,ia,ja,ad,au,al,m,b,x,t,v,r,p,z,tol,
     .                      maxit,matvec,dot,my_id,neqf1i,neqf2i,
     .                      neq_doti,i_fmapi,i_xfi,i_rcvsi,i_dspli)
c **********************************************************************
c *                                                                    *
c *   Subroutine PBICGSTAB                                             *
c *                                                                    *
c *   Solucao de sistemas de equacoes pelo metodo dos gradientes       *
c *   biconjugados com precondicionador diagonal para matrizes         *
c *   nao-simetricas.                                                  *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   neq    - numero de equacoes                                      *
c *   ia(*)  - ponteiro do formato CSR                                 *
c *   ja(*)  - ponteiro das colunas no formato CSR                     *
c *   ad(neq)- diagonal da matriz A                                    *
c *   au(*)  - parte triangular superior de A                          *
c *   al(*)  - parte triangular inferior de A                          *
c *   m(*)   - precondicionador diagonal                               *
c *   b(neq) - vetor de forcas                                         *
c *   x(neq) - chute inicial                                           *
c *   t(neq) - arranjo local de trabalho                               *
c *   v(neq) - arranjo local de trabalho
c *   r(neq) - arranjo local de trabalho                               *
c *   p(neq) - arranjo local de trabalho
c *   z(neq) - arranjo local de trabalho                               *
c *   tol    - tolerancia de convergencia                              *
c *   maxit  - numero maximo de iteracoes                              *
c *   matvec - nome da funcao externa para o produto matrix-vetor      *
c *   dot    - nome da funcao externa para o produto escalar           *
c *   energy - nao definido                                            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   x(neq) - vetor solucao                                           *
c *   b(neq) - modificado                                              *
c *   ad(*),al(*),au(*) - inalterados                                  *
c *                                                                    *
c **********************************************************************
      implicit none
      include 'mpif.h'
      integer neqf1i,neqf2i,neq_doti
c ... ponteiros      
      integer*8 i_fmapi,i_xfi
      integer*8 i_rcvsi,i_dspli
c .....................................................................   
      integer neq,maxit,nad,i,j,k
      integer ia(*),ja(*),my_id
      real*8  ad(*),au(*),al(*),m(*),x(*),r(*),p(*),b(*),t(*),v(*),z(*)
      real*8  dot,ddot,tol,conv,energy,d,alpha,beta,rr0,w
      real*8  time0,time
      external matvec,dot
c ======================================================================
      time0 = MPI_Wtime()
c ......................................................................
      nad = ia(neq+1)-1
      if(my_id.eq.0) print *, 'nad :',nad
c ......................................................................
c
c ... Chute inicial:
c
      do 10 i = 1, neq
         x(i) = 0.d0
   10 continue
c ----------------------------------------------------------------------
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)      
      do 100 i = 1, neq
         r(i) = b(i) - z(i)
         p(i) = r(i)
         b(i) = p(i)
         z(i) = p(i)/m(i) 
  100 continue
      d    = dot(r(1),z(1),neq_doti)
      conv = tol*dsqrt(dabs(d))
c ----------------------------------------------------------------------
      do 230 j = 1, maxit
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               z,v,neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
         rr0 = dot(b,r,neq_doti)
         alpha = rr0/dot(v,r,neq_doti)
         do 210 i = 1, neq
            x(i) = x(i) + alpha * z(i)
            b(i) = b(i) - alpha * v(i)
            z(i) = b(i) / m(i)
  210    continue
         call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),
     .               z,t, neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
         w = dot(t,b,neq_doti) / dot(t,t,neq_doti)
         do 220 i = 1, neq
            x(i) = x(i) + w*z(i)
            b(i) = b(i) - w*t(i)
  220    continue
         d = dot(b,z,neq_doti)
         if (dsqrt(dabs(d)) .lt. conv) goto 300
         beta = (dot(r,b,neq_doti) / rr0)*(alpha/w)
         do 225 i = 1, neq
             p(i) = b(i) + beta*(p(i)-w*v(i))
             z(i) = p(i)/m(i)
  225    continue
  230 continue
c ----------------------------------------------------------------------
      write(*,1200) maxit
      call stop_mef()
  300 continue
c
c ... Energy norm:
c
      call matvec(neq,ia,ja,ia(neq+2),ja(nad+1),ad,al,au,al(nad+1),x,z,
     .            neqf1i,neqf2i,i_fmapi,i_xfi,i_rcvsi,i_dspli)
      energy   = dot(x(1),z(1),neq_doti)
c ......................................................................
      time = MPI_Wtime()
      time = time-time0
c ----------------------------------------------------------------------
      if(my_id.eq.0)write(*,1100)tol,neq,j,energy,time
c ......................................................................
c     Controle de flops
      if(my_id.eq.0) write(10,'(a,a,i9,a,d20.10,a,d20.10)')
     .               "PBICGSTAB: ","it",j, " energy norm ",energy,
     .               " tol ",tol
c ......................................................................
      return
c ======================================================================
 1000 format (//,5x,'SUBROTINA BICGSTAB:',/,5x,'Coeficiente da diagonal
     . nulo ou negativo - equacao ',i7)
 1100 format(' (PBICGSTAB) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING: No convergence reached after ',i9,
     .        ' iterations !',/)
      end
