      program condNumber
      implicit none
      include 'time.fi' 
      include 'openmp.fi'
c     include 'omp_lib.h'
c     integer          maxn, maxnev, maxncv, ldv
      integer          ldv,maxncv
      parameter       (maxncv=100)
c     parameter       (maxn=256, maxnev=10, maxncv=25, 
c    .                 ldv=maxn )
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
c     real(8) v(ldv,maxncv), workl(maxncv*(maxncv+8)),
c    .        workd(3*maxn), d(maxncv,2), resid(maxn),
c    .        ax(maxn)
c     logical          select(maxncv)
      real(8),allocatable :: v(:,:),workl(:),workd(:),d(:,:),resid(:)
     .                     , ax(:)
      logical,allocatable ::  select(:)
      integer          iparam(11), ipntr(11)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character        bmat*1, which*2
      integer          ido, n, nev, ncv, lworkl, info, ierr, j, 
     .                 nconv, maxitr, mode, ishfts
      logical          rvec
      real(8)          tol, sigma
c
c     %------------%
c     | Parameters |
c     %------------%
c
      real(8) zero
      parameter        (zero = 0.0D+0)
c  
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      real(8) precision,dnrm2
      external         dnrm2, daxpy
c
c     %--------------------%
c     | Intrinsic function |
c     %--------------------%
c
      intrinsic        dabs
c ... IO
      integer nargs
      integer iData
      character arg*80
      character nameIN*80
      data iData /11/
c ... arquivo COO
      character*10 rep   
      character*7  field
      character*19  symm
      integer   rows,cols,nnz
      complex cdum
      integer idum
c ... formato COO
      real(8) ,allocatable :: aCoo(:)
      integer ,allocatable :: iaCoo(:),jaCoo(:)
      integer ier
c ... PCG
      real(8) ,allocatable :: m(:),z(:),r(:)
      external  matVecCooSym
c ... opnemp
      real(8) ,allocatable :: bufferY(:)
c ... CsrC
      real(8) ,allocatable :: adCsr(:),alCsr(:)
      integer ,allocatable :: iaCsr(:),jaCsr(:),aux(:)
      integer nadCsr 
c ...
      real(8) eigenLarge(10),eigenSmall(10)
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------%
c     | Abertura do arquivo COO |
c     %-------------------------%
c
      nargs = iargc()
      if(nargs .lt. 1) then
        print*,'Usage: ./condNumer fileIn'
        stop
      endif
      call getarg(1,arg)
      nameIn = arg
      open(iData,status='old',action= 'read'
     .    ,err= 1,file=nameIn)
      goto  2
    1 continue
      print*,'Arquivo ',trim(nameIn), ' nao encontrado.'
      stop
    2 continue
c
c     %-------------------------%
c     | leitura do arquivo COO  |
c     %-------------------------%
c
      call mminfo(idata,rep,field,symm,rows,cols,nnz)
      if( rows .ne. cols) then
        print*,'Matriz retangular:',rows,'x',cols
        stop
      endif
      if( symm .ne. 'symmetric') then
        print*,'Matriz nao e simetria.'
        stop
      endif
      allocate(iaCoo(nnz),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor iaCoo.'
        stop
      endif 
      allocate(jaCoo(nnz),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor jaCoo.'
        stop
      endif
      allocate(aCoo(nnz),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor aCoo.'
        stop
      endif
c ......................................................................
c
c ...
      call mmread(idata,rep,field,symm,rows,cols,nnz,nnz,iaCoo,jaCoo
     .           ,idum,aCoo,cdum)
c ......................................................................
      close(idata)
c ......................................................................
c
c ... csrC
      nadCsr = nnz - rows
      allocate(adCsr(rows),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor ad.'
        stop
      endif
      allocate(alCsr(nadCsr),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor al.'
        stop
      endif
      allocate(iaCsr(rows+1),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor iaCsr.'
        stop
      endif
      allocate(jaCsr(nadCsr),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor jaCsr.'
        stop
      endif
c .....................................................................
c
c ... COO -> CSRC
      allocate(aux(rows),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor rows.'
        stop
      endif
c ......................................................................
c
c ...
      call cooToCsr(iaCoo,jaCoo,aCoo
     .             ,iaCsr  ,jaCsr
     .             ,alCsr  ,adCsr ,alCsr
     .             ,rows   ,nnz   ,nadCsr   
     .             ,3      ,aux  
     .             ,.false.,.false. ,.true.)
c ......................................................................
c
c ...
      deallocate(aux)
      deallocate(iaCoo)
      deallocate(jaCoo)
      deallocate(aCoo)
c ......................................................................
c
c ......................................................................
c
c ... aloca do memoria para o PCG
      allocate(m(rows),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor m.'
        stop
      endif
      allocate(z(rows),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor z.'
        stop
      endif
      allocate(r(rows),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor r.'
        stop
      endif
c .....................................................................
c
c ... precondicionador
      call aequalb(m,adCsr,rows)
c .....................................................................
c
c ... openmpSolver
      openmpSolver = .true.
      if(openmpSolver) then
        num_threads = omp_get_max_threads()
        allocate(bufferY(num_threads*rows),stat=ier)
        if( ier .gt. 0 ) then 
          print*,'Erro na alocacao da memoria para o vetor bufferY.'
          stop
        endif
c ... particionamento da matriz por threads
        call partition_matrix(iaCsr,jaCsr,rows)
c .....................................................................
      endif
c ......................................................................
c
c     %----------------------------------------------------% 
c     |  A standard eigenvalue problem is solved           |
c     | (BMAT = 'I'). NEV is the number                    |
c     | of eigenvalues to be approximated.  The user can   |
c     | modify NEV, NCV, WHICH to solve problems of        |
c     | different sizes, and to get different parts of the |
c     | spectrum.  However, The following conditions must  |
c     | be satisfied:                                      |
c     |             NEV + 1 <= NCV <= MAXNCV               | 
c     %----------------------------------------------------% 
c
      n   = rows
      ldv = n
      nev = 10
c ... ncv >= nev +1
      if( n .lt. 5) then 
        ncv =  nev+1
      else
        ncv =  max(25,nev+1)
      endif
      if ( ncv .lt. (nev + 1)) then
         print *, ' ERROR: ncv >= nev + 1'
         print *, ' ncv  :',ncv
         print *, ' nev  :',nev + 1
         go to 9000
      end if
      if ( ncv .gt. maxncv) then
         print *, ' ERROR: ncv > maxncv'
         print *, ' Numero da base execedido'
         print *, ' ncv     :',ncv
         print *, ' maxncv  :',maxncv
         go to 9000
      end if
c ... problema padrão de autovalor
      bmat = 'I'
c ... auto valores com maior valor algebrico
      which = 'LA'
c ... alocando vetores para o arpack
c ... autovetores
      allocate(v(n,ncv),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor v.'
        stop
      endif
c 
      allocate(workl(ncv*(ncv+8)),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor workl.'
        stop
      endif
c
      allocate(workd(3*n),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor workd.'
        stop
      endif
c ... autovalores
      allocate(d(ncv,2),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor d.'
        stop
      endif
c ... residuo 
      allocate(resid(n),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor resid.'
        stop
      endif
c ... ax 
      allocate(ax(n),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor ax.'
        stop
      endif
c .....................................................................
c
c ... ax 
      allocate(select(ncv),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor select.'
        stop
      endif
c .....................................................................
c     %--------------------------------------------------%
c     | The work array WORKL is used in DSAUPD as        |
c     | workspace.  Its dimension LWORKL is set as       |
c     | illustrated below.  The parameter TOL determines |
c     | the stopping criterion.  If TOL<=0, machine      |
c     | precision is used.  The variable IDO is used for |
c     | reverse communication and is initially set to 0. |
c     | Setting INFO=0 indicates that a random vector is |
c     | generated in DSAUPD to start the Arnoldi         |
c     | iteration.                                       |
c     %--------------------------------------------------%
c
      lworkl = ncv*(ncv+8)
      tol   = zero 
      info  = 0
      ido   = 0
c
c     %---------------------------------------------------%
c     | This program uses exact shifts with respect to    |
c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of DSAUPD is used     |
c     | (IPARAM(7) = 1).  All these options may be        |
c     | changed by the user. For details, see the         |
c     | documentation in DSAUPD.                          |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = 300
      mode   = 1
c      
      iparam(1) = ishfts 
      iparam(3) = maxitr 
      iparam(7) = mode 
c
c     %-------------------------------------------%
c     | M A I N   L O O P (Reverse communication) |
c     |          5 Large Eigenvalues              |
c     %-------------------------------------------%
c
      timeArnold1 = 0.0d+0
      timeArnold2 = 0.0d+0
      timePcg     = 0.0d+0
      timeCsr     = 0.0d+0
 10   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine DSAUPD and take | 
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c 
         timeArnold1 = get_time() - timeArnold1
         call dsaupd ( ido, bmat, n, which, nev, tol, resid, 
     .                 ncv, v, ldv, iparam, ipntr, workd, workl,
     .                 lworkl, info )
         timeArnold1 = get_time() - timeArnold1
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %--------------------------------------%
c           | Perform matrix vector multiplication |
c           |              y <--- OP*x             |
c           | The user should supply his/her own   |
c           | matrix vector multiplication routine |
c           | here that takes workd(ipntr(1)) as   |
c           | the input, and return the result to  |
c           | workd(ipntr(2)).                     |
c           %--------------------------------------%
c
c           call matVecCooSym (n,nnz,iaCoo,jaCoo,aCoo
c    .                        ,workd(ipntr(1)), workd(ipntr(2)))
            timeCsr    = get_time() - timeCsr   
            call matvec_csrcSym(n,iaCsr,jaCsr,alCsr,adCsr
     .                         ,workd(ipntr(1)),workd(ipntr(2)))
            timeCsr    = get_time() - timeCsr
c
c           %-----------------------------------------%
c           | L O O P   B A C K to call DSAUPD again. |
c           %-----------------------------------------%
c
            go to 10
c
         end if 
c
c     %----------------------------------------%
c     | Either we hmatVecCooSyme convergence or there is |
c     | an error.                              |
c     %----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %--------------------------%
c        | Error message. Check the |
c        | documentation in DSAUPD. |
c        %--------------------------%
c
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '
c
      else 
c
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |
c        | Post-Process using DSEUPD.                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |  
c        |                                           |
c        | Eigenvectors may also be computed now if  |
c        | desired.  (indicated by rvec = .true.)    | 
c        %-------------------------------------------%
c           
         rvec = .true.
c
         call dseupd ( rvec, 'All', select, d, v, ldv, sigma, 
     .        bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     .        iparam, ipntr, workd, workl, lworkl, ierr )

c        %----------------------------------------------%
c        | Eigenvalues are returned in the first column |
c        | of the two dimensional array D and the       |
c        | corresponding eigenvectors are returned in   |
c        | the first NEV columns of the two dimensional |
c        | array V if requested.  Otherwise, an         |
c        | orthogonal basis for the invariant subspace  |
c        | corresponding to the eigenvalues in D is     |
c        | returned in V.                               |
c        %----------------------------------------------%
c
         if ( ierr .ne. 0) then
c
c            %------------------------------------%
c            | Error condition:                   |
c            | Check the documentation of DSEUPD. |
c            %------------------------------------%
c
             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _seupd. '
             print *, ' '
c
         else
c
             nconv =  iparam(5)
             do 20 j=1, nconv
c
c               %---------------------------%
c               | Compute the residual norm |
c               |                           |
c               |   ||  A*x - lambda*x ||   |
c               |                           |
c               | for the NCONV accurately  |
c               | computed eigenvalues and  |
c               | eigenvectors.  (iparam(5) |
c               | indicates how many are    |
c               | accurate to the requested |
c               | tolerance)                |
c               %---------------------------%
c
c               call matVecCooSym(n,nnz,iaCoo,jaCoo,aCoo,v(1,j), ax)
c ...
                if(openmpSolver) then
                  timeCsr    = get_time() - timeCsr   
                  call matvec_csrcSymOmp(n,iaCsr,jaCsr,alCsr,adCsr
     .                              ,v(1,j),ax
     .                              ,bufferY)
                  timeCsr    = get_time() - timeCsr
                else
                  timeCsr    = get_time() - timeCsr   
                  call matvec_csrcSym(n,iaCsr,jaCsr,alCsr,adCsr
     .                               ,v(1,j),ax)
                  timeCsr    = get_time() - timeCsr
                endif
c .....................................................................   
                call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / dabs(d(j,1))
c
 20          continue
c
c            %-------------------------------%
c            | Display computed residuals    |
c            %-------------------------------%
c
             write(*,'(a)')'Autovalores de maiores:' 
             do j = 1, nev
               write(*,'(a,i3,2es16.6e3)')'Autovalor:',j,d(j,1),d(j,2) 
               eigenLarge(j) = d(j,1)
             enddo
c             call dmout(6, nconv, 2, d, maxncv, -6,
c     .            'Ritz values and relative residuals')
         end if
c
c        %------------------------------------------%
c        | Print additional convergence information |
c        %------------------------------------------%
c
         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' ' 
            print *, ' No shifts could be applied during implicit',
     .               ' Arnoldi update, try increasing NCV.'
            print *, ' '
         end if      
c
         print *, ' '
         print *, ' LARGE EIGENVAULES'
         print *, ' ================'
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated',
     .            ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', 
     .              nconv 
         print *, ' The number of Implicit Arnoldi update',
     .            ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
c
      end if
c .....................................................................
c
      lworkl = ncv*(ncv+8)
      tol   = zero 
      info  = 0
      ido   = 0
c
      ishfts = 1
      maxitr = 300
      mode   = 3
      sigma  = zero 
c
      iparam(1)  = ishfts 
      iparam(3)  = maxitr 
      iparam(7)  = mode
c 
c     %-------------------------------------------%
c     | M A I N   L O O P (Reverse communication) |
c     |          5 small Eigenvalues              |
c     %-------------------------------------------%
c
 15   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine DSAUPD and take | 
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c   
         timeArnold2 = get_time() - timeArnold2
         call dsaupd ( ido, bmat, n, which, nev, tol, resid, 
     .                 ncv, v, ldv, iparam, ipntr, workd, workl,
     .                 lworkl, info )
         timeArnold2 = get_time() - timeArnold2
c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %----------------------------------------%
c           | Perform y <-- OP*x = inv[A-sigma*I]*x. |
c           | The user only need the linear system   |
c           | solver here that takes workd(ipntr(1)) |
c           | as input, and returns the result to    |
c           | workd(ipntr(2)).                       |
c           %----------------------------------------%
c
c
c          call diagPrecondCoo(m,nnz,iaCoo,jaCoo,aCoo)
c ..
           timePcg    = get_time() - timePcg   
           if(openmpSolver) then   
             call pcgOmp(n,nnz,iaCsr,jaCsr,alCsr,adCsr,alCsr
     .               ,m
     .               ,workd(ipntr(1))
     .               ,workd(ipntr(2))
     .               ,z,r
     .               ,tol    
     .               ,60000
     .               ,1
     .               ,.false.
     .               ,bufferY)
           else    
             call pcg(n,nnz,iaCsr,jaCsr,alCsr,adCsr,alCsr
     .               ,m
     .               ,workd(ipntr(1))
     .               ,workd(ipntr(2))
     .               ,z,r
     .               ,tol     
     .               ,60000
     .               ,1
     .               ,.false.)
           endif 
           timePcg    = get_time() - timePcg   
c .....................................................................
c:
c           %-----------------------------------------%
c           | L O O P   B A C K to call DSAUPD again. |
c           %-----------------------------------------%
c
            go to 15
c
         end if 
c
c     %----------------------------------------%
c     | Either we hmatVecCooSyme convergence or there is |
c     | an error.                              |
c     %----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %--------------------------%
c        | Error message. Check the |
c        | documentation in DSAUPD. |
c        %--------------------------%
c
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '
c
      else 
c
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |
c        | Post-Process using DSEUPD.                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |  
c        |                                           |
c        | Eigenvectors may also be computed now if  |
c        | desired.  (indicated by rvec = .true.)    | 
c        %-------------------------------------------%
c           
         rvec = .true.
c
         call dseupd ( rvec, 'All', select, d, v, ldv, sigma, 
     .        bmat, n, which, nev, tol, resid, ncv, v, ldv, 
     .        iparam, ipntr, workd, workl, lworkl, ierr )

c        %----------------------------------------------%
c        | Eigenvalues are returned in the first column |
c        | of the two dimensional array D and the       |
c        | corresponding eigenvectors are returned in   |
c        | the first NEV columns of the two dimensional |
c        | array V if requested.  Otherwise, an         |
c        | orthogonal basis for the invariant subspace  |
c        | corresponding to the eigenvalues in D is     |
c        | returned in V.                               |
c        %----------------------------------------------%
c
         if ( ierr .ne. 0) then
c
c            %------------------------------------%
c            | Error condition:                   |
c            | Check the documentation of DSEUPD. |
c            %------------------------------------%
c
             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _seupd. '
             print *, ' '
c
         else
c
             nconv =  iparam(5)
             do 25 j=1, nconv
c
c               %---------------------------%
c               | Compute the residual norm |
c               |                           |
c               |   ||  A*x - lambda*x ||   |
c               |                           |
c               | for the NCONV accurately  |
c               | computed eigenvalues and  |
c               | eigenvectors.  (iparam(5) |
c               | indicates how many are    |
c               | accurate to the requested |
c               | tolerance)                |
c               %---------------------------%
c
c               call matVecCooSym(n,nnz,iaCoo,jaCoo,aCoo,v(1,j), ax)
c ...
                if(openmpSolver) then
                  timeCsr    = get_time() - timeCsr   
                  call matvec_csrcSymOmp(n,iaCsr,jaCsr,alCsr,adCsr
     .                              ,v(1,j),ax
     .                              ,bufferY)
                  timeCsr    = get_time() - timeCsr
                else
                  timeCsr    = get_time() - timeCsr   
                  call matvec_csrcSym(n,iaCsr,jaCsr,alCsr,adCsr
     .                               ,v(1,j),ax)
                  timeCsr    = get_time() - timeCsr
                endif
                timeCsr    = get_time() - timeCsr
c .....................................................................   
                call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / dabs(d(j,1))
c
 25          continue
c
c            %-------------------------------%
c            | Display computed residuals    |
c            %-------------------------------%
c
             write(*,'(a)')'Autovalores menores:' 
             do j = 1, nev
               write(*,'(a,i3,2es16.6e3)')'Autovalor:',j,d(j,1),d(j,2) 
              eigenSmall(j) = d(j,1)

             enddo
c             call dmout(6, nconv, 2, d, maxncv, -6,
c     .            'Ritz values and relative residuals')
         end if
c
c        %------------------------------------------%
c        | Print additional convergence information |
c        %------------------------------------------%
c
         if ( info .eq. 1) then
            print *, ' '
            print *, ' Maximum number of iterations reached.'
            print *, ' '
         else if ( info .eq. 3) then
            print *, ' ' 
            print *, ' No shifts could be applied during implicit',
     .               ' Arnoldi update, try increasing NCV.'
            print *, ' '
         end if      
c
         print *, ' '
         print *, ' SMALL EIGENVAULES'
         print *, ' ================'
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated',
     .            ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', 
     .              nconv 
         print *, ' The number of Implicit Arnoldi update',
     .            ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
c
      end if
c .....................................................................
      print *, ' CONDITION NUMBER'
      print *, ' ================'
      write(*,'(a,es16.6e3)')'Time Arnold 1            :',timeArnold1
      write(*,'(a,es16.6e3)')'Time Arnold 2            :',timeArnold2
      write(*,'(a,es16.6e3)')'Time Pcg                 :',timePcg    
      write(*,'(a,es16.6e3)')'Time Csr                 :',timeCsr    
      do j = 1, nev
        write(*,'(a,es16.6e3)')'Maior Autovalor        :',eigenLarge(j)
      enddo
      do j = 1, nev
        write(*,'(a,es16.6e3)')'Menor Autovalor        :',eigenSmall(j)
      enddo
      write(*,'(a,es16.6e3)')'Numero de condicionamento:'
     .       ,eigenLarge(nev)/eigenSmall(nev)
c .....................................................................
c
c ...
 9000 continue
c
      stop
      end
c ********************************************************************* 
c
c ********************************************************************* 
c * MATVECCOOSYM : matriz vetor para matriz simetrica no formato COO  * 
c ********************************************************************* 
      subroutine matVecCooSym(neq,nnz,ia,ja,a,x,y)
      implicit none
      real(8) a(*),af,x(*),y(*)
      integer neq,nnz,ia(*),ja(*)
      integer i,nLin,nCol

c
      y(1:neq) = 0.0d0
c
c ... 
      do i = 1, nnz
        nLin = ia(i)
        nCol = ja(i)
        af   =  a(i)
        if(nLin .eq. nCol) then
          y(nLin) = y(nLin) + af*x(nLin)
        else
          y(nLin) = y(nLin) + af*x(nCol)
          y(nCol) = y(nCol) + af*x(nLin)
        endif
      enddo     
c .....................................................................
c
c ...
      return  
      end
c ********************************************************************* 
c
c ********************************************************************* 
c * DIAGPRECONDCOO : Precodicionador diagonal de uma matriz armazenada*
c * no formato COO                                                    * 
c ********************************************************************* 
      subroutine diagPrecondCoo(m,nnz,iaCoo,jaCoo,aCoo) 
      integer iaCoo(*),jaCoo(*)
      real(8) m(*),aCoo(*)
      integer i,nLin
      do i = 1, nnz 
        nLin = iaCoo(i)  
        if( nLin .eq. jaCoo(i) ) m(nLin) = aCoo(i)
      enddo
      return
      end
c ********************************************************************* 
