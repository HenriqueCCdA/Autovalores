c *********************************************************************
c * COOTOSKYLINE: converte a matriz do formato COO para o skyLine     *
c * ------------------------------------------------------------------*
c * Parametros de Entrada :                                           *
c * ------------------------------------------------------------------*
c * al         -> nao definido                                        *
c * ad         -> nao definido                                        *
c * au         -> nao definido                                        *
c * jp         -> altura da coluna do formato skyLine                 *
c * nnzSkyline -> numero total de elementos np SkyLine                *
c * iaCoo      -> arranjo da com as linhas no formato COO             *
c * jaCoo      -> arranjo da com as colunas no formato COO            *
c * neq        -> numero de equacoes                                  *
c * nnzCoo     -> numero total de nao zeros no formato COO            *
c * unsymm     -> matriz nao simetrica                                *
c * ------------------------------------------------------------------*
c * Parametros de Saida:                                              *
c * ------------------------------------------------------------------*
c * al(nnzSky) -> coeficientes da matriz inferior no formato skyline  *
c * ad(neq)    -> coeficientes da matriz diagonal                     *
c * au(nnzSky) -> coeficientes da matriz suoerio  no formato skyline  *
c * ------------------------------------------------------------------*
c *********************************************************************
      subroutine cooToSkyline(al,ad,au,jp,nnzSkyline,aCoo,iaCoo,jaCoo
     .                       ,neq,nnzCoo,unsymm)
      implicit none
      real*8 al(*),ad(*),au(*)
      real*8 aCoo(*)
      integer iaCoo(*),jaCoo(*)
      integer nnzSkyline,jp(*)
      integer nnzCoo,neq
      integer i,linha,coluna
      integer k0,jh
      logical unsymm
c ...
      ad(1:neq)        = 0.0d0
      au(1:nnzSkyline) = 0.0d0
      if(unsymm) al(1:nnzSkyline) = 0.0d0
c .....................................................................
c
c ... diagonal
      do i = 1, nnzCoo
        linha  = iaCoo(i)
        coluna = jaCoo(i)
        if( linha .eq. coluna ) ad(linha) = aCoo(i)
      enddo
c .....................................................................
c
c ...
      do i = 1, nnzCoo
        linha  = iaCoo(i)
        coluna = jaCoo(i)
        if(unsymm) then
          print*,'Nao implementado para matriz nao simetricas.'
          stop
        else
c ... triangular superior
          if( linha .lt. coluna ) then
            jh         = jp(coluna) - jp(coluna-1)
            if( jh .le. coluna) then
              k0         = jp(coluna) + linha - coluna + 1
              au(k0   )  = aCoo(i)
            endif
c ... triangular inferior
          elseif( linha .gt. coluna ) then
            jh         = jp(linha) - jp(linha-1)
            if( jh .le. linha) then
              k0         = jp(linha) + coluna - linha + 1
              au(k0   )  = aCoo(i)
            endif
          endif
        endif 
      enddo
c .....................................................................
      return
      end
c *********************************************************************
c
c *********************************************************************
c * HEIGHTSKYLINEFROMCOO : Calcula a altura da coluna apartir de      *
c * uma matriz no formato COO                                         *
c * ------------------------------------------------------------------*
c * Parametros de Entrada :                                           *
c * ------------------------------------------------------------------*
c * iaCoo      -> arranjo da com as linhas no formato COO             *
c * jaCoo      -> arranjo da com as colunas no formato COO            *
c * neq        -> numero de equacoes                                  *
c * nnzCoo     -> numero total de nao zeros no formato COO            *
c * jp         -> nao definido                                        *
c * nnzSkyline -> nao definido                                        *
c * unsymm     -> matriz nao simetrica                                *
c * ------------------------------------------------------------------*
c * Parametros de Saida:                                              *
c * ------------------------------------------------------------------*
c * jp         -> alura de coluna no formato SkyLine                  *
c * nnzSkyline -> numero total de elementos np SkyLine                *
c * ------------------------------------------------------------------*
c *********************************************************************
      subroutine heightSkylineFromCoo(iaCoo,jaCoo,neq,nnzCoo 
     .                         ,jp,nnzSkyLine,unsymm)
      implicit none
      integer iaCoo(*),jaCoo(*)
      integer nnzSkyline,jp(*)
      integer nnzCoo,neq,jh
      integer i,linha,coluna
      logical unsymm
c ...
      jp(1:neq) = 0
      do i = 1, nnzCoo
        linha  = iaCoo(i)
        coluna = jaCoo(i)
        if(unsymm) then
          print*,'Nao implementado para matriz nao simetricas.'
          stop
        else
c ... triangular superior
          if( linha .lt. coluna ) then
            jh         = coluna - linha
            jp(coluna) = max(jh,jp(coluna))
c ... triangular inferior
          elseif( linha .gt. coluna ) then
            jh         = linha - coluna
            jp(linha)  = max(jh,jp(linha))
          endif
        endif 
      enddo
c .....................................................................
c 
c ...
      jp(1) = 0
      do i = 2, neq
        jp(i) = jp(i) + jp(i-1)
      enddo
      nnzSkyLine = jp(neq)
c .....................................................................
      return
      end
c********************************************************************* 
c
c********************************************************************* 
c* CSRTOCOO : conveter do formato CSR para COO                       * 
c*-------------------------------------------------------------------* 
c* Parametros de entrada:                                            * 
c*-------------------------------------------------------------------* 
c* iaCoo  -> indefinido                                              * 
c* jaCoo  -> indefinido                                              * 
c* aCoo   -> indefinido                                              * 
c* ia     -> vetor CSR                                               * 
c* ja     -> vetor CSR                                               * 
c* au     -> matrix de coeficientes                                  * 
c* ad     -> matrix de coeficientes                                  * 
c* al     -> matrix de coeficientes                                  * 
c* neq    -> numero de equacoes                                      * 
c* bin    -> matriz binaria                                          * 
c* nadCsr -> numero de termos nao nulos                              * 
c*-------------------------------------------------------------------* 
c* Parametros de saida:                                              * 
c*-------------------------------------------------------------------* 
c* iaCoo -> numero da linha                                          * 
c* jaCoo -> numero da coluna                                         * 
c* aCoo  -> valor                                                    * 
c*-------------------------------------------------------------------* 
c* OBS:                                                              * 
c*-------------------------------------------------------------------* 
c********************************************************************* 
      subroutine csrToCoo(iaCoo,jaCoo,aCoo
     .                   ,ia   ,ja
     .                   ,al   ,ad ,au
     .                   ,neq  ,nadCsr
     .                   ,unsym,bin)
      implicit none
      integer iaCoo(*),jaCoo(*)
      real*8  aCoo(*)
      integer ia(*),ja(*)
      real*8 al(*),ad(*),au(*)   
      integer neq ,nadCsr
      logical unsym,bin
      integer nl,nc,kk,n,ipoint
c
      kk = 0
      do nl = 1, neq
        kk     = kk + 1
        iaCoo(kk) = nl
        jaCoo(kk) = nl
        if(bin) then
          aCoo(kk) = 1.0
        else
          aCoo(kk) = ad(nl)
        endif
        n      = ia(nl+1) - ia(nl)
        ipoint = ia(nl)
        do nc = ia(nl), ia(nl+1) - 1
          kk = kk + 1
          iaCoo(kk)  = nl
          jaCoo(kk)  = ja(nc)
          if(bin) then
            aCoo(kk) = 1.0
          else
            if( jaCoo(kk) .lt. nl) then
              aCoo(kk) = al(nc)
            endif
            if( jaCoo(kk) .gt. nl .and. unsym) then
              aCoo(kk) = au(nc)
            endif
          endif
        enddo
      enddo
c
      return
      end
c *********************************************************************
c
c********************************************************************* 
c* COOTOCSR : conveter do formato COO para CSR                       * 
c*-------------------------------------------------------------------* 
c* Parametros de entrada:                                            * 
c*-------------------------------------------------------------------* 
c* iaCoo  -> indefinido                                              * 
c* jaCoo  -> indefinido                                              * 
c* aCoo   -> indefinido                                              * 
c* ia     -> vetor CSR                                               * 
c* ja     -> vetor CSR                                               * 
c* au     -> matrix de coeficientes                                  * 
c* ad     -> matrix de coeficientes                                  * 
c* al     -> matrix de coeficientes                                  * 
c* neq    -> numero de equacoes                                      * 
c* bin    -> matriz binaria                                          * 
c* nadCsr -> numero de termos nao nulos                              * 
c*-------------------------------------------------------------------* 
c* Parametros de saida:                                              * 
c*-------------------------------------------------------------------* 
c* iaCoo -> numero da linha                                          * 
c* jaCoo -> numero da coluna                                         * 
c* aCoo  -> valor                                                    * 
c*-------------------------------------------------------------------* 
c* OBS:                                                              * 
c*-------------------------------------------------------------------* 
c********************************************************************* 
      subroutine cooToCsr(iaCoo,jaCoo,aCoo
     .                   ,ia   ,ja
     .                   ,al   ,ad    ,au
     .                   ,neq  ,nnz   ,nadCsr   
     .                   ,code ,aux  
     .                   ,upper,diag  ,lower)
      implicit none
      integer iaCoo(*),jaCoo(*)
      real*8  aCoo(*)
      integer ia(*),ja(*),aux(*)
      real*8  al(*),ad(*),au(*)   
      integer neq ,nnz,nadCsr
      integer nLin,nCol,ipont,i,j,code
      logical aCooUp,aCooL,transP
      logical upper,diag,lower
c ...
      aCooUp = .false.
      aCooL  = .false.
      transP = .false.
      do i = 1, nnz
        nLin = iaCoo(i)
        nCol = jaCoo(i)
        if( nCol .gt. nLin) then
          aCooUp = .true.
        else if( nCol .lt. nLin) then
          aCooL  = .true.
        endif
      enddo
c .....................................................................
c
c ...
      aux(1:neq) = 0
c .....................................................................
c
c ... vetor ia
      do i = 1, nnz
        nLin = iaCoo(i)
        nCol = jaCoo(i)
        if     ( upper .and. nCol  .gt. nLin) then
          aux(nLin) = aux(nLin) + 1
        else if( diag  .and. nCol  .eq. nLin) then
          aux(nLin) = aux(nLin) + 1
        else if( lower  .and. nCol .lt. nLin) then
          aux(nLin) = aux(nLin) + 1
        endif
      enddo
c
      ia(1) = 1
      do i = 2, neq + 1
        ia(i) = ia(i-1) + aux(i-1)
      enddo
c .....................................................................
c
c ... vetor ja
      aux(1:neq) = 0
c
      do i = 1, nnz
        nLin  = iaCoo(i)
        nCol  = jaCoo(i)
        ipont = ia(nLin) + aux(nLin)
        if     ( upper .and. nCol .gt. nLin) then
          ja(ipont) = nCol
          aux(nLin) = aux(nLin) + 1
        else if( diag  .and. nCol .eq. nLin) then
          ja(ipont) = nCol
          aux(nLin) = aux(nLin) + 1
        else if( lower .and. nCol .lt. nlin) then
          ja(ipont) = nCol
          aux(nLin) = aux(nLin) + 1
        endif
      enddo
c .....................................................................
c
c ...
      call sortGraphCsr(ia,ja,neq)
c .....................................................................
c 
c ... csr padrao - (ia,ja,a) 
      if(code .eq. 1) then
        print*,'Csr padrao nao implementado!!'
        stop
c .....................................................................
c
c ... csrD       - (ia,ja,ad,a)
      else if(code .eq. 2) then
        print*,'CsrD nao implementado!!'
        stop
c .....................................................................
c
c ... csrC       - (ia,ja,al,ad,au)
      else if(code .eq. 3) then
        do i = 1, nnz
          nLin  = iaCoo(i)
          nCol  = jaCoo(i)
c ... diagonal principal 
          if( nLin .eq. nCol) then 
            ad(nLin) = aCoo(i)
c ... triagonal inferior
          else if( nLin .gt. nCol) then
            do j = ia(nLin), ia(nLin+1)
              if( nCol .eq. ja(j) ) then 
                al(j) = aCoo(i)
              endif
            enddo
c ... triagonal superior
          else if( nLin .lt. nCol) then
            do j = ia(nLin), ia(nLin+1)
              if( nLin .eq. ja(j) ) then 
                au(j) = aCoo(i)
              endif
            enddo
          endif
        enddo
      endif
c .....................................................................
c 
c ...
      return
      end
c *********************************************************************
     
