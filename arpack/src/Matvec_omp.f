c **********************************************************************
c *                                                                    *
c *   MATVEC_CSRCSYM_OMP: produto matriz-vetor y = Ax  (A simetrica),  *
c *                    coef. de A no formato CSRC e grafo simetrico.   *
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
c *   thread_y - buffer de equacoes para o vetor y (openmp)            *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   y(neq) - vetor contendo o resultado do produto y = Ax            *
c *                                                                    *
c **********************************************************************
      subroutine matvec_csrcSymOmp(neq,ia,ja,al,ad,x,y,thread_y)
      implicit none
      include 'openmp.fi'
      include 'time.fi'
      integer neq,ia(*),ja(*),i,j,k,jak,inc
      real*8  ad(*),al(*),x(*),y(*),s,t,xi
c ... ponteiros      
      real*8  thread_y(*)
c ......................................................................
c$    thread_id = omp_get_thread_num() + 1
      do i = 1, num_threads
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i)+inc, thread_begin(i)+inc-1
            thread_y(j) = 0.d0
         enddo
c$omp end do
      enddo
      inc = (thread_id - 1)*neq
c$omp barrier
      do 110 i = thread_begin(thread_id), thread_end(thread_id)
         y(i) = 0.d0
         xi = x(i)
         t  = ad(i)*xi
         do 100 k = ia(i), ia(i+1)-1
            jak = ja(k)
            s   = al(k)
            t   = t + s*x(jak)
            jak = jak + inc
            thread_y(jak) = thread_y(jak) + s*xi
  100    continue
         thread_y(i+inc) = t
  110 continue
c$omp barrier
c
c ... Accumulate thread_y(i) into y(i)
c
      do i = 1, num_threads
         inc = (i-1)*neq
c$omp do
         do j = thread_height(i), thread_end(i)
            y(j) = y(j) + thread_y(j+inc)
         end do
c$omp end do
      end do
c ......................................................................
      return
      end
c **********************************************************************
c 
c **********************************************************************
c *                                                                    *
c *   DOT_PAR_OMP: produto escalar omp(versao de diretiva orfã,  funcao*
c *   chamada de uma regiao paralela, portanto nao ha a necessidade de *
c *   abertura de uma nova regia paralela).                            *
c *                                                                    *
c *   Parametros de entrada:                                           *
c *                                                                    *
c *   a        - vetor                                                 *
c *   b        - vetor                                                 *
c *   neq_doti - dimensao dos vetores                                  *
c *                                                                    *
c *   Parametros de saida:                                             *
c *                                                                    *
c *   dot    - valor do produto escalar                                *
c *                                                                    *
c **********************************************************************
      real*8 function dot_par_omp(a,b,neq_doti)
      implicit none
      include 'time.fi'
      include 'openmp.fi'
      integer i,neq_doti
      real*8  a(*),b(*),tmp
c ......................................................................
c$omp single
      omp_dot = 0.d0
c$omp end single
c$omp do reduction(+:omp_dot)
      do 100 i = 1, neq_doti
         omp_dot = omp_dot + a(i)*b(i)
  100 continue
c$omp end do
      dot_par_omp = omp_dot
c$omp barrier 
c ......................................................................
      return
      end
c *********************************************************************
