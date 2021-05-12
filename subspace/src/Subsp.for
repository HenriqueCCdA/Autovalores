c***********************************************************************
c*                                                                     *
c*     SUBSP: Subspace iteration to extract lowest nf eigenpairs       *
c*            Solves:  (A - shift*B)*V = B*V*d for V and d             *
c*                                                                     *
c*     Inputs:                                                         *
c*        a(*)      - Coefficient matrix (tangent stiffness)           *
c*        b(*)      - Coefficient matrix (mass or geometric stiffness) *
c*        nf        - Number of pairs to converge                      *
c*        nv        - Size of subspace problem > or = nf               *
c*        neq       - Size of A,B                                      *
c*        imas      - Switch: =1 for consistent B; =2 for diagonal B   *
c*        shift     - Value of shift                                   *
c*        tol       - Tolerance to converge pairs                      *
c*        prt       - Flag, output iteration arrays if true            *
c*        its       - Maximum number of subspace iterations (=25)      *
c*                                                                     *
c*     Scratch:                                                        *
c*        t(neq)     - Working vector                                  *
c*        g(nv*(nv+1)/2) - Projection of A - shift*B matrix            *
c*        h(nv*(nv+1)/2) - Projection of B           matrix            *
c*        dp(nv)     - Previous iteration values                       *
c*        dtol(nv)   - Tolerance of eigenvalue iterations              *
c*        p(nv,nv)   - Eigenvectors of G*P = H*P*d                     *
c*                                                                     *
c*     Outputs:                                                        *
c*        v(neq,nv) - Eigenvectors                                     *
c*        d(nv)     - Eigenvalues                                      *
c*                                                                     *
c***********************************************************************
      implicit none
      logical   conv,prt
      integer   i,j,k,n
      integer   neq,nv,nad
      integer   iow,idata
      integer   nf, imas,its,itlim,it,itt,num,nmas,imtyp
      real*8    shift,tol,dm,dr,epmac,smachn,det,prodE
      real*8    dpp(4)
      data      iow/10/,idata/11/
      integer   ier
      integer nargs
      character arg*80

c ... market matrix file
      character*10 rep   
      character*7  field
      character*19  symm
      integer   rows,cols,nnz
      character*80 nameIN
      complex cdum
      integer idum
      integer*8 mem1,lneq,lnv,lnad
c ... formato COO
      real*8 ,allocatable :: aCoo(:)
      integer,allocatable :: iaCoo(:),jaCoo(:)
c .....................................................................
      integer,allocatable :: jp(:)
      real*8   ,allocatable :: ad(:),au(:),bd(:),bu(:),v(:,:),t(:)
      real*8   ,allocatable :: g(:),h(:),d(:),dp(:),dtol(:),p(:,:)
c .....................................................................
      open(iow,file='autovalres.dat')
c.......................................................................
c
c ... Variaveis de inicializacao:
c
      prt    = .true.
      imtyp = 1
      imas  = 2
      tol   = 1.d-15
      its   = 50
      shift = 0.0d0
c .....................................................................
c
c ... lendo a matrix no formato COO
c
c ... abrindo arquivo
c      nameIn = 'matriz3x3.mtx'
      nargs = iargc()
      call getarg(1,arg)
      nameIn = arg
      open(idata,status='old',action= 'read'
     .    ,err=10,file=nameIn)
      goto 20
   10 continue
      print*,'Arquivo ',trim(nameIn), ' nao encontrado.'
      stop
   20 continue
c .....................................................................
c
c .. alocacao de COO
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
        print*,'Erro na alocacao da memoria para o vetor au.'
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
c ... alocao de variaveis do skyline da matriz A
      neq = rows
      nv  = min(100,neq)  
c ... vetor da altura das colunas
      allocate(jp(neq),stat=ier)
      if( ier .ne. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor jp.'
        stop
      endif 
      call  heightSkylineFromCoo(iaCoo,jaCoo,neq,nnz 
     .                          ,jp,nad,.false.)
      lneq   = neq
      lnad   = nad
      mem1   = (lnad+lneq**2)*8
      print*,'convertendo para o formato skyline.'
      print*,'numero de valores no perfil da matriz:',nad
      print*,'numero de linhas:                     ',neq
      print*,'Skyline:                              ',
     .      (mem1)/(1024**2),'MB'
c .....................................................................
c
c ... coeficientes da matriz A
      allocate(ad(neq),stat=ier)
      if( ier .ne. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor ad.'
        stop
      endif
      allocate(au(nad),stat=ier)
      if( ier .ne. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor au.'
        print*,(nad*8)/(1024**3),'GB'
        stop
      endif
      call cooToSkyline(au,ad,au,jp,nad,aCoo,iaCoo,jaCoo
     .                 ,neq,nnz,.false.)
c .......................................................................
c
c ... dealloc
      deallocate(iaCoo)
      deallocate(jaCoo)
      deallocate( aCoo)
c .......................................................................
c
c .......................................................................
c
c ... alocao de variaveis
      allocate(bd(neq),stat=ier)
      if( ier .ne. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor bd.'
        stop
      endif
      if(imas .eq. 1) then
        allocate(bu(nad),stat=ier)
        if( ier .gt. 0 ) then 
          print*,'Erro na alocacao da memoria para o vetor bu.'
          stop
        endif
      endif
      allocate( v(neq,nv),stat=ier)
c      allocate( v(neq,neq),stat=ier)
      if( ier .ne. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor v.'
        print*,(neq*nv*8)/(1024**3),'GB'
        stop
      endif
      allocate( t(neq),stat=ier)
      if( ier .ne. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor t.'
        stop
      endif
      lnv = nv          
      mem1= (lnv*(lnv+1))/2
      allocate(g(mem1),stat=ier)
      if( ier .ne. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor g.'
        stop
      endif
      allocate(h(mem1),stat=ier)
      if( ier .ne. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor h.'
        stop
      endif
      allocate(d(nv),stat=ier)
      if( ier .gt. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor d.'
        stop
      endif
      allocate(dp(nv),stat=ier)
      if( ier .ne. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor dp.'
        stop
      endif
      allocate(dtol(nv),stat=ier)
      if( ier .ne. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor dtol.'
        stop
      endif
      allocate(p(nv,nv),stat=ier)
      if( ier .ne. 0 ) then 
        print*,'Erro na alocacao da memoria para o vetor p.'
        stop
      endif
      lneq = neq
      mem1 = 2*mem1 + lnv*lnv + lneq*lnv + 2*lneq + 3*lnv
      print*,'SubSpace                              ',
     .      (mem1)/(1024**2),'MB'
c ......................................................................
c
c ... 
      nf    = nv/2 
c .......................................................................
c
c ... Precisao da maquina:
c
      epmac = smachn()
c
c ... Shift:
c
      bd(1:neq) = 1.0d0
      if( imas .eq. 1) then 
        bu(1:nad) = 0.0d0
        j = jp(neq)
        do i = 1, j
          au(i) = au(i) - shift*bu(i)
        enddo
      endif
      do i = 1, neq
        ad(i) = ad(i) - shift*bd(i)
      enddo
c
c ... Fatoracao da matriz A (LDL):
c
      print*,'Fatoracao LDLt'
      call dtri(ad,au,au,jp,neq,.false.)
c
c ... Compute initial iteration vectors:
c
      call azero(v,nv*neq)
      num    = 0
      nmas   = 0
      do 100 n = 1,neq
c
c ...   Count number of values less than current shift:
c
        if(ad(n).lt.0.0d0) num = num + 1
c
c ...   Count number of nonzero masses:
c
        if(bd(n).ne.0.0d0) nmas = nmas + 1
  100 continue
c......................................
      if(prt) then
        write(iow,2002) num,shift
      endif
      nmas = nmas/nv
      i = 0
      j = 1
      do 110 n = 1,neq
        dm = bd(n)
        if(dm.ne.0.0d0) then
          v(n,j) = dm
          i      = i + 1
          if(mod(i,nmas).eq.0) j = min(j+1,nv)
        endif
  110 continue
      k = 1
      do 120 i = 1,nv
        dp(i)   = 0.0
        dtol(i) = 1.0d0
c        call scalev0(v(1,i),k,1,1,neq)
         call scalev(v(1,i),neq,k)
  120 continue
c
c ... Compute new vectors and project 'a' onto 'g':
c
      conv = .false.
      itlim = its
      if(nv.eq.nf) itlim = 1
      do 400 it = 1,itlim
        itt = it
c
c ...   Project 'b/' matrix to form 'h' and compute 'z' vectors:
c
        print*,"Projb"
        call sprojb(bd,bu,v,t,h,jp,neq,nv,imas)
c
c ...   Project 'a' matrix to form 'g'
c
        print*,"Proja"
        call sproja(ad,au,v,t,g,jp,neq,nv)
c
c ...   Solve reduced eigenproblem:
c
        if(imtyp.eq.1) then
c
c ...     Eigenproblem: 'g*p = h*p*d'; e.g., vibration modes:
c
          call geig(g,h,d,p,t,nv,1,epmac)
c
        else
c
c....     Eigenproblem: 'h*p = g*p*1/d'; e.g., structural buckling:
c
          call geig(h,g,d,p,t,nv,0,epmac)
c
          do 200 n = 1,nv
            if(d(n).ne.0.0d0) d(n) = 1.d0/d(n)
  200     continue
c
c....     Resort eigenvalues and vectors to ascending order:
c
          num = nv + 1
          do 220 n = 1,nv/2
            dm       = d(n)
            d(n)     = d(num-n)
            d(num-n) = dm
            do 210 i = 1,nv
              dm         = p(i,n)
              p(i,n)     = p(i,num-n)
              p(i,num-n) = dm
  210       continue
  220     continue
        endif
c
c....   Compute new iteration vector 'u' from 'z':
c
        do 330 i = 1,neq
c
c....     Move row of 'v' into temporary vector 't':
c
          do 300 j = 1,nv
            t(j) = v(i,j)
  300     continue
c
c....     Compute new iiteration vector entries:
c
          do 320 j = 1,nv
            v(i,j) = 0.0d0
            do 310 k = 1,nv
              v(i,j) = v(i,j) + t(k)*p(k,j)
  310       continue
  320     continue
  330   continue
c
c....   Check for convergence:
c
        do 340 n = 1,nv
          if(d(n).ne.0.0d0) dtol(n) = dabs((d(n)-dp(n))/d(n))
          dp(n) = d(n)
  340   continue
c
c....   Compute maximum error for iteration:
c
        conv = .true.
        dm   =  0.0d0
        do 350 n = 1,nf
          if(dtol(n).gt.tol) conv = .false.
          dm = max(dm,dtol(n))
  350   continue
        if(conv) goto 500
  400 continue
c
c.... Scaled vectors mass orthonormal:
c
c.............................................................
c1000   dm = 1.d0/gtan(1)
  500 continue
      dm = 1.d0
      do 510 n = 1,nv
        d(n) = (1.0/d(n) + shift)*dm
        dp(n)= dsqrt(dabs(d(n)))
  510 continue
c      k = 1
c      call scalev(v(1,2),k,1,1,neq)
      call scalev(v(1,1),neq,nf)
c
c.... Output solution values:
c
c ... determinante
      prodE = d(1)
      do n = 2,nv 
        prodE= prodE*d(n)
      enddo
c .....................................................................        
      write(iow,2000) itt
      write(iow,'(a,d17.8)') 'condition number:',d(nv)/d(1)
      write(iow,'(a,d17.8)') 'Prod autovalores:',prodE     
      do n = 1,nv
        write(iow,'(i9,2e20.8e3)')n,d(n)
      enddo
      if (prt) then
c        write(iow,2000) itt,(d(n),n=1,nv)
         if(itt.gt.1) write(iow,2001) itt,(dtol(n),n=1,nv)
         if(imtyp.eq.1) then
           write(iow,2003) (dp(n),n=1,nv)
           dm = 0.5d0/dacos(-1.0d0)
           write(iow,2004) (dp(n)*dm,n=1,nv)
           write(iow,2007)
           dr = 1.d0/dm
           do 530 i = 0,nv-1,4
             do 520 j = 1,min(4,nv-i)
               if(dabs(dp(i+j)).gt.tol*10.d0) then
                  dpp(j) = dr/dp(i+j)
               else
                  dpp(j) = 0.0d0
               endif
  520        continue
             write(iow,2008) (dpp(j),j=1,min(4,nv-i))
  530      continue
         endif
      endif
      stop
c=======================================================================
2000  format(/'  SUBSPACE: Current eigenvalues, iteration',i4/)
2001  format( '  SUBSPACE: Current residuals,   iteration',i4/
     +        (5x,1p,4d17.8))
2002  format(5x,'There are',i4,' eigenvalues less than shift',1p,e12.4)
2003  format(/'  SUBSPACE: Square root of eigenvalues (rad/t)'/
     +        (5x,1p,4d17.8))
2004  format(/'  SUBSPACE: Square root of eigenvalues (Hz.)'/
     +        (5x,1p,4d17.8))
2005  format('  Completed subspace iteration',i4,' Time =',f9.2,
     +       '  Max. tol =',1p,e12.5)
2006  format(/'  SUBSPACE: Inverse eigenvalues, iteration',i4/
     +        (5x,1p,4d17.8))
2007  format(/'  SUBSPACE: Period in units of time      (T)')
2008  format(5x,1p,4d17.8)
      end
      subroutine sproja(ad,au,v,t,g,jp,neq,nv)
c***********************************************************************
c*                                                                     *
c*     SPROJA: Compute subspace projection of 'a' to form 'g'          *
c*                                                                     *
c*     Inputs:                                                         *
c*       ad(*) - Diagonal entries for b matrix                         *
c*       au(*) - Off-diagonal entries for symmetric matrix (upper part)*
c*       v(neq,*) - Set of iteration vectors                           *
c*       neq      - Number of equations in A                           *
c*       nv       - Size of projected matrix                           *
c*                                                                     *
c*     Scratch:                                                        *
c*        t(neq)   - Working vector                                    *
c*                                                                     *
c*     Outputs:                                                        *
c*        g(*)     - Projected matrix V_trans * A * V                  *
c*                                                                     *
c***********************************************************************
      implicit none
      integer jp(*)
      integer neq,nv
      integer  i,j,k
      real*8   ad(*),au(*),v(neq,*),t(neq),g(*)
      real*8   dot,en
c.......................................................................
c
c ... Forward reduce eigenvector estimates
c
      k = 0
      do 200 j = 1,nv
c
c ...   Copy vector 'v' into 'z' and solve equations:
c
        call pmove (v(1,j),t(1),neq)
c
        call dsolv(ad,au,au,v(1,j),jp,neq,en,.false.)
c
c ...   Compute projection of stiffness:
c
        do 100 i = 1,j
          k = k + 1
          g(k) = dot(v(1,i),t(1),neq)
  100   continue
  200 continue
      return
      end
      subroutine sprojb(bd,bu,v,t,h,jp,neq,nv,imas)
c***********************************************************************
c*                                                                     *
c*     Purpose: Compute subspace projection of 'b' to form 'h'         *
c*                                                                     *
c*     Inputs:                                                         *
c*       bd(*) - Diagonal entries for b matrix                         *
c*       bu(*) - Off-diagonal entries for symmetric matrix (upper part)*
c*       v(neq,*) - Set of iteration vectors                           *
c*       neq      - Number of equations in b                           *
c*       nv       - Size of projected matrix                           *
c*       imas     - Mass type: 1 = consistent; 2 = diagonal.           *
c*                                                                     *
c*     Scratch:                                                        *
c*        t(neq)   - Working vector                                    *
c*                                                                     *
c*     Outputs:                                                        *
c*        h(*)     - Projected matrix V_trans * B * V                  *
c*                                                                     *
c***********************************************************************
      implicit none
      integer jp(*)
      integer neq,nv,imas,i,j,k
      real*8  bd(*),bu(*),v(neq,*),t(*),h(*),dot
c ......................................................................
c
c ... Compute 'z' and 'b' projection to form 'h'
c
      do 1000 j = 1,nv
c
c ..... Consistent mass:
c
        if(imas.eq.1) then
           call skyprod(bd,bu,v(1,j),t,jp,neq)
c
c ..... Lumped mass:
c
        else
           do 100 i = 1,neq
              t(i) = v(i,j)*bd(i)
  100      continue
        endif

c       Project 'z' and 'v' vectors to form 'h'

        k = j*(j+1)/2
        do 200 i = j,nv
          h(k) = dot(t,v(1,i),neq)
          k = k + i
  200   continue
        do 300 i = 1,neq
          v(i,j) = t(i)
  300   continue
 1000 continue
      return
      end
      subroutine skyprod(ad,au,p,v,jp,neq)
c***********************************************************************
c*                                                                     *
c*     SKYPROD:  Peform a matrix vector product (A*p) for a symmetric  *
c*               matrix stored in skyline form.                        *
c*                                                                     *
c*     Inputs:                                                         *
c*       ad(*) - Diagonal entries for A array                          *
c*       au(*) - Off-diagonal entries for symmetric matrix (upper part)*
c*       p(*)  - Specified vector for product                          *
c*       jp(*) - Skyline pointer                                       *
c*       neq   - Number of equations                                   *
c*                                                                     *
c*     Outputs:                                                        *
c*       v(*)   - Matrix product of A*p.                               *
c*                                                                     *
c***********************************************************************
      implicit none
      real*8  ad(*),au(*),p(*),v(*),an
      integer neq
      integer jp(*)
      integer ni,jh,k0,k,n0
c ......................................................................
c
c ... Diagonal part
c
      do 100 ni = 1,neq
        v(ni) = ad(ni)*p(ni)
  100 continue
c
c ... Loop over columns:
c
      do 210 ni = 2,neq
c
c ...... Column height:
c         
         jh  = jp(ni)-jp(ni-1)
c
c ...... Starting row/column:
c
         k0  = ni-jh
c
c ...... Pointer to first au in column ni:
c
         n0 = jp(ni-1)+1
c
c ...... Loop over column/row height:
c
         do 200 k = k0, ni-1
           an = au(n0)
c
c ......... lower part:
c
           v(ni) = v(ni) + an*p(k)
c
c ......... upper part:
c
           v(k)  = v(k)  + an*p(ni)
           n0 = n0 + 1
  200    continue
  210 continue
      return
      end
      subroutine scalev(v,numnp,nv)
c***********************************************************************
c*                                                                     *
c*     SCALEV: Scale vector to have maximum element of +1.0            *
c*                                                                     *
c*     Inputs:                                                         *
c*        v(ndf,*) - Vector of values                                  *
c*        pdf(*)   - DOF to scale on                                   *
c*        ndm      - Space dimension of mesh                           *
c*        ndf      - DOF's/node (maximum)                              *
c*        numnp    - Number of nodes                                   *
c*                                                                     *
c*     Outputs:                                                        *
c*        v(ndf,*) - Unit vector of values                             *
c*                                                                     *
c***********************************************************************
      implicit none
      integer  numnp,nv
      integer  i,n
      real*8   v(numnp,*),vmax
c.......................................................................
c
c.... Locate maximum:
c
      do 300 i = 1, nv
        vmax = 0.0d0
        do 100 n = 1,numnp
            vmax = max(vmax,dabs(v(n,i)))
  100   continue
c.... Perform scaling:
        if(vmax.gt.0.0d0) then
           vmax = 1.d0/vmax
           do 200 n = 1,numnp
              v(n,i) = v(n,i)*vmax
  200      continue
        else
           write(*,*) 'Zero length vector in SCALEV'
        endif
  300 continue
      return
      end
      subroutine scalev0(v,pdf,ndm,ndf,numnp)
c***********************************************************************
c*                                                                     *
c*     SCALEV: Scale vector to have maximum element of +1.0            *
c*                                                                     *
c*     Inputs:                                                         *
c*        v(ndf,*) - Vector of values                                  *
c*        pdf(*)   - DOF to scale on                                   *
c*        ndm      - Space dimension of mesh                           *
c*        ndf      - DOF's/node (maximum)                              *
c*        numnp    - Number of nodes                                   *
c*                                                                     *
c*     Outputs:                                                        *
c*        v(ndf,*) - Unit vector of values                             *
c*                                                                     *
c***********************************************************************
      implicit none
      integer  pdf(*),ndm,ndf,numnp
      integer  i,n
      real*8   v(ndf,numnp),vmax
c.......................................................................
c
c.... Locate maximum:
c
      vmax = 0.0d0
      do 110 i = 1,ndm
        if(pdf(i).ge.1.and.pdf(i).le.ndf) then
          do 100 n = 1,numnp
            vmax = max(vmax,dabs(v(pdf(i),n)))
  100     continue
        endif
  110 continue
c
c.... Perform scaling:
c
      if(vmax.gt.0.0d0) then
        vmax = 1.d0/vmax
        do 210 n = 1,numnp
          do 200 i = 1,ndf
            v(i,n) = v(i,n)*vmax
  200     continue
  210   continue
      else
        write(*,*) 'Zero length vector in SCALEV'
      endif
      return
      end
      subroutine geig(g,h,d,p,t,nv,nvs,epmac)
c***********************************************************************
c*                                                                     *
c*     GEIG: Solve general eigenproblem 'g*p = h*p*d'                  *
c*                                                                     *
c*     Inputs:                                                         *
c*        g(*)   - Left hand projected array                           *
c*        h(*)   - Right hand projected array                          *
c*        nv     - Size of problem                                     *
c*        nvs    - Eigenvectors h orthogonal if positive               *
c*        epmac  - smallest posite number in real*8 precision          *
c*                                                                     *
c*     Outputs:                                                        *
c*        d(*)   - Eigenvalues of problem                              *
c*        p(*,*) - Eigenvectors of problem                             *
c*                                                                     *
c***********************************************************************
      implicit none
      integer  nv,nvs,ir, i, j
      real*8   g(*),h(*),d(*),p(nv,*),t(*),epmac
c.......................................................................
c
c...  Compute standard eigenvalue problem matrix 'c'
c
      call chlfwd(h,g,p,nv)
c
c...  Perform eignfunction decomposition of 'c'
c
      call eisql(g,d,t,p,nv,epmac,ir)
c
c...  Compute vectors of original problem
c
      call chlbac(h,p,nv)
c
c...  Divide eigenvectors by eigenvalue to prevent overflows
c
      ir = 0
      do 110 j = 1,nv
        ir = ir + j
        if(nvs.gt.0 .and. d(j).ne.0.0d0) then
          t(1) = 1.0d0/d(j)
        elseif(nvs.le.0 .and. h(ir).ne.0.0d0) then
          if(d(j).ne.0.0d0) then
            t(1) = dabs(d(j))/dsqrt(dabs(h(ir)))
          else
            t(1) = 1.0d0/dsqrt(dabs(h(ir)))
          endif
        else
          t(1) = 1.0d0
        endif
        if(p(j,j).lt.-0.00001d0) t(1) = -t(1)
        do 100 i = 1,nv
          p(i,j) = p(i,j)*t(1)
  100   continue
  110 continue
      return
      end
      subroutine eisql(a,d,e,z,n,epmac,ierr)
c***********************************************************************
c*                                                                     *
c*     EISQL: Compute eigen-pairs for standard eigenproblem.           *
c*                                                                     *
c*     Inputs:                                                         *
c*        a(*)   - Matrix for wanted eigenvalues                       *
c*        n      - size of eigenproblem                                *
c*        epmac  - smallest posite number in real*8 precision          *
c*                                                                     *
c*     Outputs:                                                        *
c*        d(n)   - Eigenvalues                                         *
c*        z(n,n) - Eigenvectors                                        *
c*        ierr   - Error indicator                                     *
c*                                                                     *
c*     Scratch:                                                        *
c*        e(*)   - Working vector                                      *
c*                                                                     *
c***********************************************************************
      implicit none
      integer i,j,k,l,m,n,ierr,jp1,ii,l1,mml,n2
      real*8  b,c,f,g,h,hh,p,r,s,scale,epmac
      real*8  a(*),d(*),e(*),z(n,n)
c.......................................................................
c
c...  Eispac QL algorithm adapted from 'tred2' and 'tql2'
c
c.......................................
      n2 = 0
      do 1100 i = 1,n
        do 1000 j = 1,i
          n2 = n2 + 1
          z(i,j) = a(n2)
 1000   continue
 1100 continue
c.......................................
      if(n.eq.1) goto 3000
      n2 = n + 2
      do 2300 ii = 2,n
         i = n2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if(l.lt.2) then
            e(i) = z(i,l)
         else
            do 2120 k = 1,l
               scale = scale + dabs(z(i,k))
 2120       continue
            if(scale.ne.0.0d0) goto 2130
            e(i) = z(i,l)
         endif
         goto 2200
c..............................................
 2130    do 2140 k = 1,l
            z(i,k) = z(i,k)/scale
            h = h + z(i,k)*z(i,k)
 2140    continue
         f = z(i,l)
         g = -sign(dsqrt(h),f)
         e(i) = scale*g
         h = h - f*g
         z(i,l) = f - g
         f = 0.0d0
         do 2170 j = 1,l
            z(j,i) = z(i,j)/h
            g = 0.0d0
            do 2150 k = 1,j
               g = g + z(j,k)*z(i,k)
 2150       continue
            jp1 = j + 1
            if(l.ge.jp1) then
               do 2160 k = jp1,l
                  g = g + z(k,j)*z(i,k)
 2160          continue
            endif
            e(j) = g/h
            f = f + e(j)*z(i,j)
 2170    continue
         hh = f/(h+h)
         do 2190 j = 1,l
            f = z(i,j)
            g = e(j) - hh*f
            e(j) = g
            do 2180 k = 1,j
               z(j,k) = z(j,k) - f*e(k) - g*z(i,k)
 2180       continue
 2190    continue
c.....................................................
 2200    d(i) = h
 2300 continue
c.....................................................
c
c...  Set transformation array for ql:
c
c....................................................
 3000 d(1) = z(1,1)
      z(1,1) = 1.0d0
      e(1) = 0.0d0
      ierr = 0
      if(n.eq.1) return
      do 3100 i = 2,n
        l = i - 1
        if(d(i).ne.0.0d0) then
          do 3030 j = 1,l
            g = 0.0d0
            do 3010 k = 1,l
              g = g + z(i,k)*z(k,j)
 3010       continue
            do 3020 k = 1,l
              z(k,j) = z(k,j) - g*z(k,i)
 3020       continue
 3030     continue
        endif
        d(i) = z(i,i)
        z(i,i) = 1.0d0
        do 3040 j = 1,l
           z(i,j) = 0.0d0
           z(j,i) = 0.0d0
 3040   continue
 3100 continue
c......................................................................
c
c...  Begin 'QL' algorithm on tridagonal matrix now stored in 'd' and 'e
c
      do 4000 i = 2,n
        e(i-1) = e(i)
 4000 continue
      f = 0.0d0
      b = 0.0d0
      e(n) = 0.0d0
c......................................................................
      do 5300 l = 1,n
        j = 0
        h = epmac*(dabs(d(l)) + dabs(e(l)))
        if(b.lt.h) b = h
        do 5010 m = l,n
          if(dabs(e(m)).le.b) goto 5100
 5010   continue
c.................................................
 5100   if(m.ne.l) then
c......................................................................
 5200     if(j.eq.30) goto 7000
          j = j + 1
          l1 = l + 1
          g = d(l)
          p = (d(l1)-g)/(e(l)+e(l))
          r = dsqrt(p*p+1.0d0)
          d(l) = e(l)/(p+sign(r,p))
          h = g - d(l)
          do 5210 i = l1,n
            d(i) = d(i) - h
 5210     continue
          f = f + h
          p = d(m)
          c = 1.0d0
          s = 0.0d0
          mml = m - l
          do 5230 ii = 1,mml
             i = m - ii
             g = c*e(i)
             h = c*p
             if(dabs(p).ge.dabs(e(i))) then
                c = e(i)/p
                r = dsqrt(c*c+1.0d0)
                e(i+1) = s*p*r
                s = c/r
                c = 1.0d0/r
             else
                c = p/e(i)
                r = dsqrt(c*c+1.0d0)
                e(i+1) = s*e(i)*r
                s = 1.0d0/r
                c = c*s
             endif
             p = c*d(i) - s*g
             d(i+1) = h + s*(c*g + s*d(i))
             do 5220 k = 1,n
                h = z(k,i+1)
                z(k,i+1) = s*z(k,i) + c*h
                z(k,i  ) = c*z(k,i) - s*h
 5220        continue
 5230     continue
          e(l) = s*p
          d(l) = c*p
          if(dabs(e(l)).gt.b) goto 5200
        endif
        d(l) = d(l) + f
 5300 continue
c...................................................................... 
      do 6500 ii = 2,n
        i = ii - 1
        k = i
        p = d(i)
        do 6210 j = ii,n
          if(dabs(d(j)).gt.dabs(p)) then
            k = j
            p = d(j)
          endif
 6210   continue
        if(k.ne.i) then
          d(k) = d(i)
          d(i) = p
          do 6220 j = 1,n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
 6220     continue
        endif
 6500 continue
      return
c......................................................................
 7000 ierr = l
      return
      end
      subroutine chlbac(u,s,nn)
c***********************************************************************
c*                                                                     *
c*      CHLBAC: Back substitution for Cholesky factors in eigen        *
c*               solutions.                                            *
c*                                                                     *
c*     Inputs:                                                         *
c*       u(*)   - Unreduced array                                      *
c*       s(*,*) - Factored array of matrix                             *
c*       nn     - Size of arrays                                       *
c*                                                                     *
c*     Outputs:                                                        *
c*       u(*)   - Solution after back substitution                     *
c*                                                                     *
c***********************************************************************
      implicit none
      integer  i,j,jd, nn
      real*8   u(*),s(nn,nn)
c.......................................................................
c...  Compute eigenvalues of general linear problem by backsubstitution
c
      j  = nn
      jd = nn*(nn+1)/2
      do 100 i = 1,nn
        s(nn,i) = s(nn,i)*u(jd)
  100 continue
c.......................
      do 210 j = nn,2,-1
        jd = jd - j
        do 200 i = 1,nn
          call colbac(u(jd+1),s(1,i),u(jd),j-1)
  200   continue
  210 continue
      return
      end
      subroutine chlfwd(u,g,s,nn)
c***********************************************************************
c*                                                                     *
c*    CHLFWD: Use Cholesky factors to project onto a standard          *
c*            eigenproblem                                             *
c*                                                                     *
c*    Inputs:                                                          *
c*        g(*)  - Symmetric projected matrix                           *
c*        u(*)  - Upper factor for projection                          *
c*        nn    - Size of arrays                                       *
c*                                                                     *
c*    Outputs:                                                         *
c*        s(*,*) - Projected array                                     *
c*                                                                     *
c***********************************************************************
      implicit none
      integer  i,j,id,im,jd,nn
      real*8   u(*),g(*),s(nn,nn),dot
c.......................................................................
c     Choleski factorization of a symmetric, positive definite matrix

      u(1) = 1.d0/dsqrt(dabs(u(1)))
      jd   = 1
      do 200 j = 2,nn
        id = 0
        do 100 i = 1,j-1
          if(i.gt.1) u(jd+i) = u(jd+i) - dot(u(id+1),u(jd+1),i-1)
          id = id + i
          u(jd+i) = u(jd+i)*u(id)
  100   continue
        u(jd+j) = 1.d0/dsqrt(dabs(u(jd+j) - dot(u(jd+1),u(jd+1),j-1)))
        jd = jd + j
  200 continue
c
c...  Perform forward solutions to get projected matrix
c
      s(1,1) = g(1)*u(1)
      id = 1
      do 400 i = 2,nn
        s(1,i) = g(id+1)*u(1)
        im = i - 1
        jd = 0
        do 300 j = 1,im
         s(i,j) = (g(id+j) - dot(u(id+1),s(1,j),im))*u(id+i)
         if(j.gt.1) s(j,i) = (g(id+j)-dot(u(jd+1),s(1,i),j-1))*u(jd+j)
         jd = jd + j
  300   continue
        id = id + i
        s(i,i) = (g(id) - dot(u(id-im),s(1,i),im))*u(id)
  400 continue
c
c...  Complete projection
c
      g(1) = s(1,1)*u(1)
      jd = 2
      do 600 j = 2,nn
        g(jd) = s(j,1)*u(1)
        id = 2
        do 500 i = 2,j
          im = i - 1
          g(jd+im) = (s(j,i) - dot(u(id),g(jd),im))*u(id+im)
          id = id + i
  500   continue
        jd = jd + j
  600 continue
      return
      end
      subroutine colbac(u,s,d,jj)
c***********************************************************************
c*                                                                     *
c*     COLBAC: Backsubstitution macro for eigen solution               *
c*                                                                     *
c*     Inputs:                                                         *
c*        s(*)  - Unreduced column                                     *
c*        u(*)  - Column of upper array already reduced                *
c*        d     - Solution value for 'u' column                        *
c*        jj    - Length to reduce                                     *
c*                                                                     *
c*     Outputs:                                                        *
c*        s(*)  - Reduced column                                       *
c*                                                                     *
c***********************************************************************
      implicit none
      integer  j,jj
      real*8   d,u(*),s(*)
c.......................................................................
      do 100 j = 1,jj
        s(j) = s(j) - u(j)*s(jj+1)
  100 continue
      s(jj) = s(jj)*d
      return
      end
      subroutine colred(au,xj,nn, b)
c***********************************************************************
c*                                                                     *
c*     COLRED: Columnwise reduction for back substitution              *
c*                                                                     *
c*     Inputs:                                                         *
c*        au(*)   - Upper column of reduced array A                    *
c*        xj      - Solution of reduced column                         *
c*        nn      - Length to reduce                                   *
c*        b(*)    - Unreduced column                                   *
c*                                                                     *
c*     Outputs:                                                        *
c*        b(*)    - Reduced column                                     *
c*                                                                     *
c***********************************************************************
      implicit none
      integer  n,nn
      real*8   xj,au(*),b(*)
c.......................................................................
      do 100 n = 1,nn
        b(n) = b(n) - au(n)*xj
  100 continue
      return
      end
      real*8 function smachn()
c **********************************************************************
c *                                                                    *
c *   SMACHN: calcula o menor numero real*8 positivo da maquina        *
c *                                                                    *
c **********************************************************************
      smachn = 1.0d0
100   smachn = smachn*0.5d0
      if(smachn+1.d0 .ne. 1.d0) go to 100
      smachn = 2.0d0*smachn
      return
      end
      subroutine getcon(a,nmax,epmac,dops)
c      * * F E A P * * A Finite Element Analysis Program
c....  Copyright (c) 1984-2000: Robert L. Taylor
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Perform installation parameter computations
c      Inputs:
c      Outputs:
c         epmac   - Smallest number that can be added to 1.0
c         dops    - Estimate on computer speed
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none
      integer   i, nmax
      real*4    tarr1(2),tarr2(2) , etime , tt
      real*8    a(nmax),epmac,dops
c
c     Compute machine epsilon estimate
c
      epmac = 1.0d0
100   epmac = epmac*0.5d0
      if(1.d0 .ne. 1.d0 + epmac) goto 100
      epmac = 2.0d0*epmac

c     Compute estimate of floating point speed

c     nmas = min(500000,maxm/2-1)
      do i = 1,nmax
        a(i) = 0.33333333333333d0
      end do
      dops = 0.0d0
      tt = etime(tarr1)
      do i = 1, nmax
        dops = dops + a(i)*a(i)
      end do
      tt = etime(tarr2)
      if(tarr2(1)-tarr1(1) .gt.0.0) then
        dops = real(nmax)/(tarr2(1) - tarr1(1))*2.d0
      endif
      end
      subroutine numass(d,n,np)
c***********************************************************************
c*                                                                     *
c*     NUMASS: Conta o numero de coeficientes da diagonal diferentes   *
c*             de zero.                                                *
c*     Inputs:                                                         *
c*        d(*)    - coeficientes da diagonal                           *
c*        n       - numero de coeficientes                             *
c*                                                                     *
c*     Outputs:                                                        *
c*        np      - numero de coeficientes diferentes de zero          *
c*                                                                     *
c***********************************************************************
      implicit none
      integer  n,np,i
      real*8   d(*)
c.......................................................................
      np = 0
      do 100 i = 1,n
        if (d(i).ne.0.d0) np = np + 1
  100 continue
      if (np .eq. 0) then
      print*, 'NUMASS: Coeficientes da diagonal iguais a ZERO !'
      stop
      endif
      return
      end
