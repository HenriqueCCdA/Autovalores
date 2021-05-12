c*****************************Svn**************************************
c*$Date: 2011-12-02 22:10:15 -0200 (Fri, 02 Dec 2011) $                
c*$Rev: 959 $                                                          
c*$Author: henrique $                                                  
c**********************************************************************
c *********************************************************************
c * PARTITION_MATRIX : dividi o trabalho do matevec entre as threads  *
c * por linhas e inicializa a estruturas do buffer do matvec          *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro     *
c *                     coeficiente nao-nulo da linha   i             *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa     *
c *               a posicao k no vetor au                             *
c *   neq       - numero de equacoes                                  *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c *  thread_begin(1:nthreads) - primeira linha do sistema da thread i *
c *  thread_end(1:nthreads)   - ultima linha do sistema da thread i   *
c *  thread_heigth(1:nthread) - altura onde a thread escreve no vetor *
c *                             y                                     *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine partition_matrix(ia,ja,neq)
      implicit none
      include 'openmp.fi'
      integer ia(*),ja(*),neq
      integer i
c
      call partition_csrc_bynonzeros(ia,neq)
      call compute_effective_work(ia,ja,neq)
c      
c ... checa se a divisao da matriz ocorreu sem problemas
c
      do i = 1 ,num_threads
        if (thread_height(i) .le. 0 ) goto 1000
        if ( thread_begin(i) .le. 0 ) goto 1000
        if ( thread_end(i)   .le. 0 ) goto 1000
      enddo
c .....................................................................      
      return
c ... controle de erro       
1000  continue 
      print*,'***error: divisao da matrix para o openmp falhou!'
      print*,'Diminua o numero de threads usado ou desabilite o openmp.'
      print*,'Log da partition_matrix:'
      do i = 1 ,num_threads
        print*,'thread id',i,'height' , thread_height(i)
        print*,'thread id',i,'begin ' , thread_begin(i)
        print*,'thread id',i,'end   ' , thread_end(i)
      enddo
      stop
      end
c *********************************************************************
c
c *********************************************************************
c * PARTITION_CSRC_BYNOZEROS :dividi o trabalho do matvec entre as    *
c * threads considerando valores nao nulos e inicializa a estruturas  *
c * do buffer do matvec                                               *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro     *
c *                     coeficiente nao-nulo da linha   i             *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa     *
c *               a posicao k no vetor au                             *
c *   neq       - numero de equacoes                                  *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c *  thread_begin(1:nthreads) - primeira linha do sistema da thread i *
c *  thread_end(1:nthreads)   - ultima linha do sistema da thread i   *
c *  thread_size(1:nthread)   - numero de termos no buffer da thread i*
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine partition_csrc_bynonzeros(ia,neq)
      implicit none
      include 'openmp.fi'
      integer ia(*),nnzr,neq
      integer mean_variables,line,thread_size(num_threads),tam,i
      integer*8 nad  
c
      nad = ia(neq+1)-1
      mean_variables = (2*nad + neq)/num_threads + 1
      line = 2
      thread_begin(1) = 1
      do i = 1, num_threads - 1
        thread_size(i) = 0
        tam = 0
c
 100    tam = 2*(ia(line) - ia(line - 1)) + 1
        if (thread_size(i) + tam .le. mean_variables) then
           thread_size(i) = thread_size(i) + tam
           thread_end(i) = line - 1
           thread_begin(i+1) = line
           line = line + 1
           goto 100
        endif
      enddo
      thread_size(num_threads) = 2*(ia(neq+1) -
     .   ia(thread_begin(num_threads))) +
     .   neq + 1 - thread_begin(num_threads)
      thread_end(num_threads) = neq
      return
      end
c *********************************************************************
c
c *********************************************************************
c * PARTITION_CSRC_EVENLY : dividi o trabalho entre as threads        *
c * considerando valores nulos e nao nulos  inicializa a estruturas   *
c * do buffer do matvec                                               *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c *   neq       - numero de equacoes                                  *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c *  thread_begin(1:nthreads) - primeira linha do sistema da thread i *
c *  thread_end(1:nthreads)   - ultima linha do sistema da thread i   *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine partition_csrc_evenly(neq)
      implicit none
      include 'openmp.fi'
      integer neq,split
      integer i,j
c
c$omp parallel private(split)
!$    thread_id = omp_get_thread_num()
      split = mod(neq, num_threads)
      thread_begin(thread_id+1) = thread_id*(neq/num_threads)
c
      if (thread_id .lt. split) then
         thread_begin(thread_id+1) = thread_begin(thread_id+1) +
     .      thread_id + 1
         thread_end(thread_id+1) = thread_begin(thread_id+1) +
     .      neq/num_threads
      else
         thread_begin(thread_id+1) = thread_begin(thread_id+1) +
     .      split + 1
         thread_end(thread_id+1) = thread_begin(thread_id+1) +
     .      neq/num_threads - 1
      endif
c$omp end parallel
      return
      end
c *********************************************************************
c
c *********************************************************************
c * COMPUTE_EFFECTIVE_WORK: Calculo do trabalho efeitivo por thread   *
c * ----------------------------------------------------------------- *
c * parametros de entrada :                                           *
c * ----------------------------------------------------------------- *
c *   ia(neq+1) - ia(i) informa a posicao no vetor au do primeiro     *
c *                     coeficiente nao-nulo da linha   i             *
c *   ja(neq+1) - ja(k) informa a coluna do coeficiente que ocupa     *
c *               a posicao k no vetor au                             *
c *   neq       - numero de equacoes                                  *
c * ----------------------------------------------------------------- *
c * parametros de saida                                               *
c * ----------------------------------------------------------------- *
c *  thread_heigth(1:nthread) - altura onde a thread escreve no vetor *
c *                             y                                     *
c * ----------------------------------------------------------------- *
c *********************************************************************
      subroutine compute_effective_work(ia,ja,neq)
      implicit none
      include 'openmp.fi'
      integer ia(*),ja(*),neq,h,i,k
c
c$omp parallel private(h) num_threads(num_threads)
c$    thread_id = omp_get_thread_num()
      h = thread_begin(thread_id+1)
      do i = thread_begin(thread_id+1), thread_end(thread_id+1) 
         do k = ia(i), ia(i+1) - 1
           h = min(h, ja( k ))
          enddo
      enddo
      thread_height(thread_id+1) = h
c$omp end parallel
      return
      end
c **********************************************************************
c
c **********************************************************************
c * GET_BUFFER_SIZE: quantidade de memoria usada no buffer do matvec   *
c * ------------------------------------------------------------------ *
c * Parametro de entrada :                                             *
c * cod - 1 Bytes, 2 KBytes, 3 MBytes, 4 GBytes                        *
c * ------------------------------------------------------------------ *
c * ------------------------------------------------------------------ *
c **********************************************************************
      real*8 function get_buffer_size(tp,neq)
      implicit none
      include 'openmp.fi'
      integer nbytes
      integer neq
      character*2 tp
c ... 8 byts vector
      nbytes = 8*num_threads*neq
c ......................................................................
c
c ...      
      if(tp .eq.' B') then
        get_buffer_size = nbytes
      else if (tp .eq. 'KB') then
        get_buffer_size = nbytes / 1024
      else if (tp .eq. 'MB') then
        get_buffer_size = nbytes / (1024**2)
      else if (tp .eq. 'GB') then
        get_buffer_size = nbytes / (1024**3)
      endif
c ......................................................................
c
c ...
      return
      end
c ......................................................................
c **********************************************************************
      
