      subroutine ppttxjcallme2born (M2, ccP, Mt, Mh)

      implicit none

      real*8 :: Mass_E, Mass_M, Mass_L, Mass_U, Mass_D, Mass_S,
     &     Mass_C,  Width_C, Mass_B, Width_B, Mass_T, Width_T,
     &     Mass_W, Width_W, Mass_Z, Width_Z, Mass_H, Width_H,
     &     Coupl_Alpha_QED, Coupl_Alpha_QCD
      integer :: last_switch, amp_switch, amp_switch_rescue,
     &     use_coli_cache, check_Ward_tree, check_Ward_loop,
     &     out_symmetry

      real*8 ,  dimension(0:3,5):: P
      real*8  ccP(20)
      real*8  M2
      real*8  Mt, Mh
      integer i,j

      Mass_E  = 0d0
      Mass_M  = 0d0
      Mass_L  = 0d0
      Mass_U  = 0d0
      Mass_D  = 0d0
      Mass_S  = 0d0
      Mass_C  = 0d0
      Width_C = 0d0
      Mass_B  = 0d0 
      Width_B = 0d0 
c      Mass_T  = 172d0 
      Mass_T  = Mt 
      Width_T = 0d0
      Mass_W  = 80.399d0 
      Width_W = 0d0 
      Mass_Z  = 91.118d0
      Width_Z = 0d0
c      Mass_H  = 120d0
      Mass_H  = Mh
      Width_H = 0d0
      Coupl_Alpha_QED = 1/128d0
      Coupl_Alpha_QCD = 0.1258086856923967d0
      last_switch = 1
      amp_switch = 1
      amp_switch_rescue = 7
      use_coli_cache = 1
      check_Ward_tree = 0
      check_Ward_loop = 0
      out_symmetry = 1

c      do i = 1, 20
c      print *, "=> " , ccP(i)
c      enddo

      do i = 0, 3
      do j = 1, 5
      P(i,j) = ccP(i + (j - 1) * 4 + 1)
c      print *, "P(i,j) ", P(i,j)
      enddo
      enddo
      call parameters_init( Mass_E, Mass_M, Mass_L, Mass_U, Mass_D, 
     & Mass_S,Mass_C, Width_C, Mass_B, Width_B, Mass_T, Width_T,
     & Mass_W, Width_W, Mass_Z, Width_Z, Mass_H, Width_H,
     & Coupl_Alpha_QED, Coupl_Alpha_QCD,
     & last_switch, amp_switch, amp_switch_rescue,
     & use_coli_cache, check_Ward_tree, check_Ward_loop, out_symmetry)
c      print *, "Param init"
      call set_permutation_ppttj_ttxggg_1([4, 5, 3, 2, 1])
c      print *, "set_permut"
      call AMP2tree_ppttj_ttxggg_1(P, M2)
c      print *, "Amp2"
      end 
