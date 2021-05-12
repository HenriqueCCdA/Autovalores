c *********************************************************************
c * PCG: Solucao de sistemas de equacoes pelo metodo dos gradientes   *
c * conjugados com precondicionador diagonal para matrizes simetricas.*  
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * neq    - numero de equacoes                                       *
c * nad    - numero de elementos nao nulos fora da diagonal principal *
c * ia     - ponteiro do arranjo do csr da matrix A                   *
c * ja     - arranjo csr da matriz A                                  *
c * ad     - diagonal principal de A                                  *
c * au(nad)- coeficiantes fora da diagonal principal de A             *
c * al(nad)- coeficiantes fora da diagonal principal de A             *
c * m(neq) - precondicianado diagonal                                 *
c * b(neq) - termo indepedente                                        *
c * x(neq) - nao definido                                             *
c * z(neq) - auxiliar                                                 *
c * r(neq) - auxiliar                                                 *
c * tol    - tolerancia do solver                                     *
c * maxit  - numero maximo de iteracao do solver                      *
c * nlit     - numero da iteraçao nao lienar                          *
c * init     - inicializacao do vetor do resultados com nulo          *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * x    - solucao do sistema                                         *
c * ----------------------------------------------------------------- *
c *********************************************************************      
      subroutine pcg(neq,nad,ia,ja,al,ad,au,m,b,x,z,r,tol,maxit
     .              ,nlit,init)
      implicit none
      include 'time.fi'
      real*8 au(*),al(*),ad(*),x(*),b(*)
      real*8 r(*),z(*),m(*)
      real*8 d,conv,tol,alpha,beta,energy
      real*8 dot,dx
      integer neq,ia(*),ja(*)
      integer i,j,maxit,nad,nlit,k
      logical init
      external dot,matVecCooSym
c ...  
      time0 = get_time()
      if(init) then      
        do i = 1, neq
          x(i) = 0.d0
        enddo
      endif
c ..................................................................... 
c
c ...
      call matvec_csrcSym(neq,ia,ja,al,ad,x,z)
c     call matVecCooSym(neq,nad,ia,ja,ad,x,z)
      do i = 1, neq
        r(i) = b(i) - z(i)
        z(i) = r(i) / m(i)
        b(i) = z(i)
      enddo
      d    = dot(r,z,neq)
      conv = tol*dsqrt(dabs(d))
c ..................................................................... 
      k = 0
      do j = 1, maxit
         call matvec_csrcSym(neq,ia,ja,al,ad,b,z)
c        call matVecCooSym(neq,nad,ia,ja,ad,b,z)
         alpha = d/dot(b,z,neq)
         dx = 0.d0
         do i = 1, neq
            x(i) = x(i) + alpha * b(i)
            r(i) = r(i) - alpha * z(i)
            z(i) = r(i) / m(i)
         enddo
         beta = dot(r,z,neq) / d
         do i = 1, neq
            b(i) = z(i) + beta * b(i)
         enddo   
         d = beta*d
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ...
         k = k + 1 
         if(k .eq. 5000) then 
           k = 0 
           write(*,1300),j,dsqrt(dabs(d)),conv
         endif
c ..................................................................... 
      enddo
c ----------------------------------------------------------------------
      write(*,1200) maxit,dsqrt(dabs(d)),conv
      stop
  300 continue
c
c ... Energy norm:
c
      call matvec_csrcSym(neq,ia,ja,al,ad,x,z)
c     call matVecCooSym(neq,nad,ia,ja,ad,x,z)
      energy   = dot(x,z,neq)
c ......................................................................
      time = get_time()
      time = time-time0
c ----------------------------------------------------------------------
      if(nlit .eq. 1) write(*,1100)tol,neq,nad,j,energy,time
c ......................................................................
c     Controle de flops
      write(13,'(a,a,i9,a,d20.10,a,d20.10)')
     .      "PCG:       ","it",j, " energy norm ",energy,
     .      " tol ",tol
c ......................................................................
      return
c ======================================================================
 1100 format(' (PCG) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING:  PCG no convergence reached after ',i9,
     .        ' iterations !',d20.6,' Res ',d20.6,' Conv'/)
 1300 format (' PCG it =',i9,' norm =',d20.10,' conv =',d20.10)
      end
c *********************************************************************

c *********************************************************************
c * PCGOMP: Solucao de sistemas de equacoes pelo metodo dos gradientes*
c * conjugados com precondicionador diagonal para matrizes simetricas.*  
c * ----------------------------------------------------------------- *
c * Parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c * neq    - numero de equacoes                                       *
c * nad    - numero de elementos nao nulos fora da diagonal principal *
c * ia     - ponteiro do arranjo do csr da matrix A                   *
c * ja     - arranjo csr da matriz A                                  *
c * ad     - diagonal principal de A                                  *
c * au(nad)- coeficiantes fora da diagonal principal de A             *
c * al(nad)- coeficiantes fora da diagonal principal de A             *
c * m(neq) - precondicianado diagonal                                 *
c * b(neq) - termo indepedente                                        *
c * x(neq) - nao definido                                             *
c * z(neq) - auxiliar                                                 *
c * r(neq) - auxiliar                                                 *
c * tol    - tolerancia do solver                                     *
c * maxit  - numero maximo de iteracao do solver                      *
c * nlit     - numero da iteraçao nao lienar                          *
c * init     - inicializacao do vetor do resultados com nulo          *
c * ----------------------------------------------------------------- *
c * Parametros de saida :                                             *
c * ----------------------------------------------------------------- *
c * x    - solucao do sistema                                         *
c * ----------------------------------------------------------------- *
c *********************************************************************      
      subroutine pcgOmp(neq,nad,ia,ja,al,ad,au,m,b,x,z,r,tol,maxit
     .                 ,nlit,init,bufferY)
      implicit none
      include 'time.fi'
      include 'openmp.fi'
      real*8 au(*),al(*),ad(*),x(*),b(*)
      real*8 r(*),z(*),m(*)
      real*8 d,conv,tol,alpha,beta,energy
      real*8 dot_par_omp,dx
      real*8 bufferY(*)
      integer neq,ia(*),ja(*)
      integer i,j,maxit,nad,nlit,k
      logical init
c ...  
      time0 = get_time()
c$omp parallel private(i,j,k,d,conv,alpha,beta,energy)
!$    thread_id = omp_get_thread_num()
      if(init) then
c$omp do      
        do i = 1, neq
          x(i) = 0.d0
        enddo
c$omp end do      
      endif
c ..................................................................... 
c
c ...
      call matvec_csrcSymOmp(neq,ia,ja,al,ad,x,z,bufferY)
c$omp do      
      do i = 1, neq
        r(i) = b(i) - z(i)
        z(i) = r(i) / m(i)
        b(i) = z(i)
      enddo
c$omp end do      
      d    = dot_par_omp(r,z,neq)
      conv = tol*dsqrt(dabs(d))
c ..................................................................... 
      k = 0
      do j = 1, maxit
         call matvec_csrcSymOmp(neq,ia,ja,al,ad,b,z,bufferY)
         alpha = d/dot_par_omp(b,z,neq)
         dx = 0.d0
c$omp do      
         do i = 1, neq
            x(i) = x(i) + alpha * b(i)
            r(i) = r(i) - alpha * z(i)
            z(i) = r(i) / m(i)
         enddo
c$omp end do      
         beta = dot_par_omp(r,z,neq) / d
c$omp do      
         do i = 1, neq
            b(i) = z(i) + beta * b(i)
         enddo   
c$omp end do      
         d = beta*d
         if (dsqrt(dabs(d)) .lt. conv) goto 300
c ...
         k = k + 1 
         if(k .eq. 5000) then 
           k = 0 
c$omp single
           write(*,1300),j,dsqrt(dabs(d)),conv
c$omp end single
         endif
c ..................................................................... 
      enddo
c ----------------------------------------------------------------------
c$omp master
      write(*,1200) maxit,dsqrt(dabs(d)),conv
      stop
c$omp end master
  300 continue
c
c ... Energy norm:
c
      call matvec_csrcSymOmp(neq,ia,ja,al,ad,x,z,bufferY)
c     call matVecCooSym(neq,nad,ia,ja,ad,x,z)
      energy   = dot_par_omp(x,z,neq)
c ......................................................................
c$omp single
      time = get_time()
      time = time-time0
c ----------------------------------------------------------------------
      if(nlit .eq. 1) write(*,1100)tol,neq,nad,j,energy,time
c ......................................................................
c     Controle de flops
      write(13,'(a,a,i9,a,d20.10,a,d20.10)')
     .      "PCG:       ","it",j, " energy norm ",energy,
     .      " tol ",tol
c$omp end single
c$omp end parallel 
c ......................................................................
      return
c ======================================================================
 1100 format(' (PCG_OMP) solver:'/
     . 5x,'Solver tol           = ',d20.6/
     . 5x,'Number of equations  = ',i20/
     . 5x,'Nad                  = ',i20/
     . 5x,'Number of iterations = ',i20/
     . 5x,'Energy norm          = ',d20.10/
     . 5x,'CPU time (s)         = ',f20.2/)
 1200 format (' *** WARNING:  PCG no convergence reached after ',i9,
     .        ' iterations !',d20.6,' Res ',d20.6,' Conv'/)
 1300 format (' PCG it =',i9,' norm =',d20.10,' conv =',d20.10)
      end
c *********************************************************************
