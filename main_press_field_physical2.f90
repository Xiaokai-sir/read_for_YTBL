  program main
  implicit none
  integer,parameter::ndx=192,my=191,ndz=3456,nt=60000
  integer,parameter::mx=192,ndy=191,mz=3456,nump=1
  integer,parameter::mzz=2200,mzz2=2200,ndz1=mz-mzz,mz3=ndz-mzz
  real,parameter::pi=3.1415926535
  integer,parameter::itst=10,iz=1,ix=1,numt=NT/ITST,t_win=9,z_win=3
  integer,parameter::mx3=10/ix,i_win_m=5
  INTEGER,PARAMETER::n_twin=2*numt/(t_win+1),n_zwin=2*mz3/(z_win+1)
  REAL,PARAMETER::DZ=0.0072722*IZ,DT=0.0005*ITST
  real,PARAMETER::lt_win=dt*real(n_twin),lz_win=dz*real(n_zwin)
  real,PARAMETER::f0=1.0/(lt_win),k0=1.0/(lz_win)
  real,PARAMETER::w_basic=2.0*pi*f0,k_basic=2.0*pi*k0
  integer,PARAMETER::n_sample=t_win*z_win*mx3,n_sample3=t_win*z_win


  logical::exist,exist1
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) ldim
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) lenr
  character(len=6) cname
  character(len=4) cname4
  character(len=80) FILE_NAME

  integer::iwin,ii,win_c
  integer::i,j,k,it,itbp,kk,IW,ixx,kx,iwn,iww,ik,itt,tt,nit,nik,ikk,jj,izz
  real::xL,ZL,a,r21,gama,r12,eps,arf
  real::w_p,w_n
  real::dqadx(1:n_zwin,1:n_twin),dqbdx(1:n_zwin,1:n_twin)
  real::dqadt(1:n_zwin,1:n_twin),dqbdt(1:n_zwin,1:n_twin)
  complex::dqwdx(1:n_zwin,1:n_twin),dqkdt(1:n_zwin,1:n_twin)
  real::press_mat(1:n_zwin,1:n_twin,1:n_sample),press_mat2(1:2*n_zwin-1,1:n_twin,1:n_sample)
  real::t_self(1:n_twin,1:n_zwin),x_self(1:n_zwin,1:n_twin)
  real::R_xt(1:n_zwin,1:n_twin),phi_kw(1:n_zwin,1:n_twin)
  real::kxc1(1:n_twin),kxc(1:n_twin),B1(1:n_twin),B1_width(1:n_twin)
  real::wxc1(1:n_zwin),wxc(1:n_zwin),B2(1:n_zwin),B2_width(1:n_zwin)
  real::swk(1:n_zwin,1:n_twin),swk2(1:n_zwin,1:n_twin),swk_LMW1(1:n_zwin,1:n_twin),swk_LMW2(1:n_zwin,1:n_twin)
  real::ave1(0:mx,0:my,0:mz3,1:7)

  REAL (KIND=4),allocatable,dimension(:,:) ::BP
  real (kind=4),allocatable,dimension(:,:,:) ::BP_TBL
  real (kind=4),allocatable,dimension(:,:,:) ::BP_TBL2
  real (kind=4),allocatable,dimension(:,:,:) ::temp0
  real (kind=4),allocatable,dimension(:,:,:) ::temp1
  real (kind=4),allocatable,dimension(:,:,:) ::temp11
  real (kind=4),allocatable,dimension(:,:,:) ::temp2
  real (kind=4),allocatable,dimension(:,:,:) ::temp22
  real (kind=4),allocatable,dimension(:,:,:) ::temp3
  real (kind=4),allocatable,dimension(:,:,:) ::temp33
  real (kind=4),allocatable,dimension(:,:,:) ::temp4
  real (kind=4),allocatable,dimension(:,:) ::s_kw
  real (kind=4),allocatable,dimension(:) ::r_tt
  real (kind=4),allocatable,dimension(:) ::r_xx
  real (kind=4),allocatable,dimension(:,:) ::am
  real (kind=4),allocatable,dimension(:,:) ::ph
  real (kind=4),allocatable,dimension(:) ::ts
  real (kind=4),allocatable,dimension(:) ::xs
  real (kind=4),allocatable,dimension(:,:) ::e_wk
  real (kind=4),allocatable,dimension(:,:) ::temp_wk
  real (kind=4),allocatable,dimension(:,:) ::temp_ax
  real (kind=4),allocatable,dimension(:,:) ::temp_bx
  real (kind=4),allocatable,dimension(:,:) ::temp_at
  real (kind=4),allocatable,dimension(:,:) ::temp_bt





  COMPLEX(KIND=4),allocatable,dimension(:):: temp_p
  COMPLEX(KIND=4),allocatable,dimension(:):: temp_fft
  COMPLEX(KIND=4),allocatable,dimension(:):: p_tk
  COMPLEX(KIND=4),allocatable,dimension(:):: p_wz
  COMPLEX(KIND=4),allocatable,dimension(:):: s_ww
  COMPLEX(KIND=4),allocatable,dimension(:):: temp_r
  COMPLEX(KIND=4),allocatable,dimension(:):: temp_rr
  COMPLEX(KIND=4),allocatable,dimension(:,:) ::p_xt
  COMPLEX(KIND=4),allocatable,dimension(:,:) ::temp_xt
  COMPLEX(KIND=4),allocatable,dimension(:,:) ::p_kw
  COMPLEX(KIND=4),allocatable,dimension(:,:) ::pp_kw
  COMPLEX(KIND=4),allocatable,dimension(:,:) ::r_xt1
  COMPLEX(KIND=4),allocatable,dimension(:,:) ::p_xw
  COMPLEX(KIND=4),allocatable,dimension(:,:) ::p_kt
  COMPLEX(KIND=4),allocatable,dimension(:,:) ::temp
  COMPLEX(KIND=4),allocatable,dimension(:,:) ::temp_w
  COMPLEX(KIND=4),allocatable,dimension(:,:) ::temp_k

  real::theta(0:mx),x(0:mx),Z(0:NDZ),y(0:ndy)

  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  real ( kind = 4 ), allocatable, dimension ( : ) :: work
  real ( kind = 4 ), allocatable, dimension ( : ) :: wsave
  
  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'pressure_post'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the pressure-behaviour code.'
  

! basic parameters for the DNS
  r21=50.0
  gama=3.0
  r12=1.0/r21
	EPS=(1.0-R12)/(1.0+R12)
  arf=1.0*eps
  a=2.0/(r21-1.0)
  XL=2.0*PI/arf
  DO I=0,MX
    X(I)=XL*I/MX
    THETA(I)=X(I)*arf
  ENDDO
  DO J=0,MY
    JJ=MY-J
    Y(JJ)=-1.0+2.0*(1+TANH(GAMA*(real(J)/MY-1))/TANH(GAMA))
  ENDDO
  ZL=8.0*PI
  DO K=0,MZ
    Z(k)=ZL*k/MZ
  ENDDO

  open(20,file='cc.dat')
  WRITE(20,*)"t_win	 z_win	n_twin  n_zwin  n_sample"
  WRITE(20,*) t_win,z_win,n_twin,n_zwin,n_sample
  WRITE(20,*)"nt	numt	itst mx3 ix i_win_method"
  WRITE(20,*) nt,numt,itst,mx3,ix,i_win_m
  WRITE(20,*)"dt dz LT LZ"
  WRITE(20,*) dt,dz,nt*0.0005,dz*1256
  WRITE(20,*)"f0  k0  w_basic  k_basic"
  WRITE(20,*) f0,k0,w_basic,k_basic
  WRITE(20,*)"LZ_win LT_win"
  WRITE(20,*) lz_win,lt_win
  close(20)
! basic parameters for the DNS
  OPEN(15,FILE='u_mean_values.out',access='direct',recl=4*(MY+1)*(MX+1)*(mz3+1),form='binary')
  read(15,rec=1) ave1(0:mx,0:my,0:mz3,1)
  read(15,rec=nump+1) ave1(0:mx,0:my,0:mz3,2)
  read(15,rec=2*nump+1) ave1(0:mx,0:my,0:mz3,3)
  read(15,rec=3*nump+1) ave1(0:mx,0:my,0:mz3,4)
  read(15,rec=4*nump+1) ave1(0:mx,0:my,0:mz3,5)
  read(15,rec=5*nump+1) ave1(0:mx,0:my,0:mz3,6)
  read(15,rec=6*nump+1) ave1(0:mx,0:my,0:mz3,7)
CLOSE(15)

! information for windows
 
! information for windows

! read the data
  allocate(BP_TBL(1:numt,1:ndz1,1:ndx))
  allocate(BP(0:ndz,0:ndx))
 
  do it=1,NT,itst
      BP=0.0
      CALL ITOA6(IT,CNAME)
      FILE_NAME='../pressure2/PRESS_BP'//CNAME//'.dat'
      INQUIRE(FILE=FILE_NAME,EXIST=EXIST1)
    IF(EXIST1) THEN
          itbp=(it)/ITST+1
!           if itst=1
!      itbp=it/itst
      OPEN(15,FILE='../pressure2/PRESS_BP'//CNAME//'.dat')
        do k=0,ndz-1
          read(15,*) (BP(k,i),i=1,ndx)
        enddo
      CLOSE(15)
      DO K=MZZ,mz-1
        KK=(K-MZZ)/IZ+1
        BP_TBL(ITBP,KK,1:NDX)=BP(K,1:ndx)-ave1(0:mx-1,my,kk-1,7) !sum(BP(k,1:NDX))/ndx
      ENDDO
      PRINT *,"READING ITbp=",ITbp
    ELSE
      PRINT *,"FILE DO NOT EXIST!"
!             ITBP=(it)/ITST
!             CYCLE
    ENDIF
   
  ENDDO
  PRINT *,"**end reading**"
  deallocate(BP)

  allocate(BP_TBL2(1:numt,1:ndz1,1:mx3))
  DO I=140,149,ix
!    II=I/ix+1
    ! if ix=1
    ii=i-139
    BP_TBL2(1:numt,1:MZ3,II)=BP_TBL(1:numt,1:MZ3,I)
  enddo
  deallocate(BP_TBL)
    PRINT *,"**end BP_TBL2**"
!cc    PRINT *,"**writting pressure**"
!cc  FILE_NAME='pressure_data/press_series004_40T.out'
!cc  open(20,file=FILE_NAME,access='direct',recl=4*numt*mz3*mx3,form='binary')
!cc  write(20,rec=1) BP_TBL2(1:numt,1:MZ3,1:MX3)
!cc  close(20)
!cc      PRINT *,"** end writting pressure**"

! read the data

! divide several samples
        PRINT *,"**divide several samples**"
        win_c=0
  allocate(temp0(1:numt,1:mz3,1:ndx))
  allocate(temp1(1:z_win,1:n_zwin,1:numt))
  allocate(temp11(1:z_win,1:n_zwin,1:numt))
  allocate(temp2(1:n_zwin,1:numt,1:z_win))
  allocate(temp22(1:n_zwin,1:numt,1:z_win))
  allocate(temp3(1:t_win,1:n_twin,1:n_zwin))
  allocate(temp33(1:t_win,1:n_twin,1:n_zwin))
  allocate(temp4(1:n_zwin,1:n_twin,1:n_sample3))

  do i=1,mx3

    do it=1,numt
      temp0(it,1:mz3,i)=BP_TBL2(it,1:mz3,i)
      call windowz(temp0(it,1:mz3,i),mz3,z_win,n_zwin,temp11(1:z_win,1:n_zwin,it),i_win_m)
      
        temp1(1:z_win,1:n_zwin,it)=temp11(1:z_win,1:n_zwin,it)
    enddo

    do ik=1,z_win
      do it=1,numt
        temp2(1:n_zwin,it,ik)=temp1(ik,1:n_zwin,it)
      enddo
    enddo

    do ik=1,z_win
      do izz=1,n_zwin
        temp22(izz,1:numt,ik)=temp2(izz,1:numt,ik)
        call windowt(temp22(izz,1:numt,ik),numt,t_win,n_twin,temp33(1:t_win,1:n_twin,izz),i_win_m)
          temp3(1:t_win,1:n_twin,izz)=temp33(1:t_win,1:n_twin,izz)
      enddo
      do itt=1,t_win
        do izz=1,n_zwin
          temp4(izz,1:n_twin,itt+(ik-1)*t_win)=temp3(itt,1:n_twin,izz)
        enddo
      enddo
    enddo

    do ii=1,n_sample3
      press_mat(1:n_zwin,1:n_twin,ii+(i-1)*n_sample3)=temp4(1:n_zwin,1:n_twin,ii)
    enddo

  ENDDO

  deallocate(temp0)
  deallocate(temp1)
  deallocate(temp11)
  deallocate(temp2)
  deallocate(temp22)
  deallocate(temp3)
  deallocate(temp33)
  deallocate(temp4)
      PRINT *,"**end dividing blocks**"

      PRINT *,"**writting press_mat**"
  FILE_NAME='BP_result_x3/press_mat.out'
  open(20,file=FILE_NAME,access='direct',recl=4*n_twin*n_zwin*10,form='binary')
  write(20,rec=1) press_mat(1:n_zwin,1:n_twin,1:10)
  close(20)
  deallocate(BP_TBL2)
  PRINT *,"**end writting press_mat**"
!symmetric 
!c  do it=1,n_twin
!c    do ii=1,n_sample
!c      press2_mat(1:n_zwin,it,ii)=press_mat(1:n_zwin,it,ii)
!c      do k=1,n_zwin-1
!c        press2_mat(n_zwin+k,it,ii)=press_mat(n_zwin-k,it,ii)
!c      enddo
!c    enddo
!c  enddo

!symmetric 


! divide several samples

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc                Time-space energy investigation code                                    !cc
!cc                                   BY Y.K. Xu                                           !cc
!cc                                                                                        !cc
!cc                     press(1:n_twin,1:n_zwin,1:win)                                     !cc
!cc                                                                                        !cc
!cc                                                                                        !cc
!cc                                                                                        !cc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! time self-correlation
  allocate(TEMP_P(1:n_twin))
  allocate(TEMP_R(1:n_twin))
  allocate(TEMP_RR(1:n_twin))
  allocate(TEMP_fft(1:n_twin))
  allocate(p_wz(1:n_twin))
  allocate(s_ww(1:n_twin))
  allocate(r_tt(1:n_twin))
  t_self=0.0
  do k=1,n_zwin
    do iwin=1,n_sample
      TEMP_P(1:n_twin)=cmplx(PRESS_MAT(K,1:n_twin,IWIN),0.0)
      call FFT1D(temp_p(1:n_twin),n_twin,temp_fft(1:n_twin))
      p_wz(1:n_twin)=temp_fft(1:n_twin)
      s_ww(1:n_twin)=p_wz(1:n_twin)*conjg(p_wz(1:n_twin))
      TEMP_R(1:n_twin)=s_ww(1:n_twin)
      call IFFT1D(temp_r(1:n_twin),n_twin,temp_rr(1:n_twin))
      R_tt(1:n_twin)=real(temp_rr(1:n_twin))
      t_self(1:n_twin,k)=t_self(1:n_twin,k)+R_tt(1:n_twin)
    enddo
    t_self(1:n_twin,k)=t_self(1:n_twin,k)/n_sample
  enddo
  open(10,file='BP_result_x3/R_tt.dat')
  do k=1,n_zwin
    write(10,1) t_self(1:n_twin,k)
  enddo 
  close(10)

  deallocate(temp_p)
  deallocate(temp_r)
  deallocate(temp_rr)
  deallocate(temp_fft)
  deallocate(p_wz)
  deallocate(s_ww)
  deallocate(r_tt)


! time self-correlation

! streamwise-space self-correlation
  allocate(TEMP_P(1:n_zwin))
  allocate(TEMP_R(1:n_zwin))
  allocate(TEMP_RR(1:n_zwin))
  allocate(TEMP_fft(1:n_zwin))
  allocate(p_tk(1:n_zwin))
  allocate(s_ww(1:n_zwin))
  allocate(r_xx(1:n_zwin))
  x_self=0.0
  do it=1,n_twin
  do iwin=1,n_sample
      TEMP_P(1:n_zwin)=cmplx(PRESS_MAT(1:n_zwin,it,IWIN),0.0)
      call FFT1D(temp_p(1:n_zwin),n_zwin,temp_fft(1:n_zwin))
      p_tk(1:n_zwin)=temp_fft(1:n_zwin)
      s_ww(1:n_zwin)=p_tk(1:n_zwin)*conjg(p_tk(1:n_zwin))
      TEMP_R(1:n_zwin)=s_ww(1:n_zwin)
      call IFFT1D(temp_r(1:n_zwin),n_zwin,temp_rr(1:n_zwin))
      R_xx(1:n_zwin)=real(temp_rr(1:n_zwin))
      x_self(1:n_zwin,it)=x_self(1:n_zwin,it)+R_xx(1:n_zwin)
    enddo
    x_self(1:n_zwin,it)=x_self(1:n_zwin,it)/n_sample
  enddo

  open(10,file='BP_result_x3/R_xx.dat')
  do it=1,n_twin
    write(10,1) x_self(1:n_zwin,it)
  enddo 
  close(10)

  deallocate(temp_p)
  deallocate(temp_r)
  deallocate(temp_rr)
  deallocate(temp_fft)
  deallocate(p_tk)
  deallocate(s_ww)
  deallocate(r_xx)
! streamwise-space self-correlation
    PRINT *,"**end self-correlation**"
! SWK calculation
  allocate(P_kw(1:n_zwin,1:n_twin))
  allocate(Pp_kw(1:n_zwin,1:n_twin))
  allocate(P_xt(1:n_zwin,1:n_twin))
  allocate(s_kw(1:n_zwin,1:n_twin))
  allocate(r_xt1(1:n_zwin,1:n_twin))
  phi_kw=0.0
  R_xt=0.0
  do iwin=1,n_sample
    p_xt(1:n_zwin,1:n_twin)=cmplx(press_mat(1:n_zwin,1:n_twin,iwin),0.0)
    call FFT2S(p_xt(1:n_zwin,1:n_twin),n_zwin,n_twin,p_kw(1:n_zwin,1:n_twin))
    pp_kw(1:n_zwin,1:n_twin)=p_kw(1:n_zwin,1:n_twin)
    s_kw(1:n_zwin,1:n_twin)=real( pp_kw(1:n_zwin,1:n_twin) * conjg(pp_kw(1:n_zwin,1:n_twin)) )
    phi_kw=phi_kw+s_kw
    call IFFT2S(cmplx(s_kw(1:n_zwin,1:n_twin),0.0),n_zwin,n_twin,r_xt1(1:n_zwin,1:n_twin))
    R_xt=R_xt+real(r_xt1)
  enddo
  R_xt=R_xt/n_sample
  phi_kw=phi_kw/n_sample

  deallocate(p_kw)
  deallocate(pp_kw)
  deallocate(p_xt)
  deallocate(s_kw)
  deallocate(r_xt1)
    FILE_NAME='BP_result_x3/Swk_DNS.out'
    open(20,file=FILE_NAME,access='direct',recl=4*n_zwin*n_twin,form='binary')
      write(20,rec=1) phi_kw(1:n_zwin,1:n_twin)
    close(20)

    FILE_NAME='BP_result_x3/Rxt_DNS.out'
    open(20,file=FILE_NAME,access='direct',recl=4*n_zwin*n_twin,form='binary')
      write(20,rec=1) R_xt(1:n_zwin,1:n_twin)
    close(20)

    PRINT *,"**end SWK**"
! SWK calculation


!LMW2 reconstruct
    allocate(p_kt(1:n_zwin,1:n_twin))
    allocate(temp(1:n_zwin,1:n_twin))
    allocate(temp_wk(1:n_zwin,1:n_twin))
    allocate(temp_XT(1:n_zwin,1:n_twin))
    allocate(P_xt(1:n_zwin,1:n_twin))
    allocate(am(1:n_zwin,1:n_twin))
    allocate(ph(1:n_zwin,1:n_twin))
    allocate(xs(1:n_zwin))
    allocate(e_wk(1:n_zwin,1:n_twin))
    allocate(temp_at(1:n_zwin,1:n_twin))
    allocate(temp_bt(1:n_zwin,1:n_twin))
    allocate(temp_k(1:n_zwin,1:n_twin))

    xs=0.0
    swk2=0.0
    wxc1=0.0
        PRINT *,"**begain SWK_LMW2**"
    do iwin=1,n_sample
      e_wk=0.0
      temp=(0.0,0.0)
      am=0.0
      ph=0.0
        p_xt(1:n_zwin,1:n_twin)=cmplx(press_mat(1:n_zwin,1:n_twin,iwin),0.0)
        temp_xt(1:n_zwin,1:n_twin)=p_xt(1:n_zwin,1:n_twin)
      do it=1,n_twin
        call FFT1D(temp_xt(1:n_zwin,it),n_zwin,temp(1:n_zwin,it))
        p_kt(1:n_zwin,it)=temp(1:n_zwin,it)
      enddo
      am=abs(p_kt)
!phi=atan2(imag(p_kt),real(p_kt))
      do ik=1,n_zwin
        do it=1,n_twin
        xs(ik)=xs(ik)+am(ik,it)*am(ik,it)
        enddo
      enddo

      do ik=1,n_zwin
        do it=3,n_twin
          dqkdt(ik,it)=(p_kt(ik,it-2)-4.0*p_kt(ik,it-1)+3.0*p_kt(ik,it))/(2.0*dt)
        enddo
        dqkdt(ik,1)=(-3.0*p_kt(ik,1)+4.0*p_kt(ik,2)-p_kt(ik,3))/(2.0*dt)
        dqkdt(ik,2)=(p_kt(ik,3)-p_kt(ik,1))/(2.0*dt)
      enddo
        PRINT *,"**iwin=",iwin

        do ik=1,n_zwin
          do it=1,n_twin
        dqadt(ik,it)=real(conjg(p_kt(ik,it))*dqkdt(ik,it))/am(ik,it)
        dqbdt(ik,it)=imag(conjg(p_kt(ik,it))*dqkdt(ik,it))/am(ik,it)/am(ik,it)
          ENDDO
        ENDDO
        temp_at=dqadt
        temp_bt=dqbdt
        temp_k=p_kt

      call SWK_LMW2_FUN(temp_k,temp_at,temp_bt,LT_win,n_zwin,n_twin,e_wk)
      temp_wk=e_wk
      Swk2=Swk2+temp_wk
!new1106
      do ik=1,n_zwin
        do it=1,n_twin
          wxc1(ik)=wxc1(ik)+am(ik,it)*am(ik,it)*dqbdt(ik,it)
        enddo
      enddo
!new1106
    enddo

    do ik=1,n_zwin
      swk_LMW2(ik,1:n_twin)=0.5*Swk2(ik,1:n_twin)/xs(ik)/w_basic
!new1106
      wxc(ik)=-wxc1(ik)/xs(ik)
!new1106
    enddo

    deallocate(p_kt)
    deallocate(temp)
    deallocate(temp_wk)
    deallocate(p_xt)
    deallocate(temp_xt)
    deallocate(am)
    deallocate(ph)
    deallocate(e_wk)
    deallocate(temp_at)
    deallocate(temp_bt)
    deallocate(temp_k)


    FILE_NAME='BP_result_x3/Swk_LMW22.out'
    open(20,file=FILE_NAME,access='direct',recl=4*n_zwin*n_twin,form='binary')
      write(20,rec=1) swk2(1:n_zwin,1:n_twin)
    close(20)

    FILE_NAME='BP_result_x3/Swk_LMW2.out'
    open(20,file=FILE_NAME,access='direct',recl=4*n_zwin*n_twin,form='binary')
      write(20,rec=1) swk_LMW2(1:n_zwin,1:n_twin)
    close(20)

    FILE_NAME='BP_result_x3/xs.out'
    open(20,file=FILE_NAME,access='direct',recl=4*n_zwin,form='binary')
      write(20,rec=1) xs(1:n_zwin)
    close(20)

    FILE_NAME='BP_result_x3/wxc.out'
    open(20,file=FILE_NAME,access='direct',recl=4*n_zwin,form='binary')
      write(20,rec=1) wxc1(1:n_zwin)
    close(20)

    FILE_NAME='BP_result_x3/B2_width.out'
    open(20,file=FILE_NAME,access='direct',recl=4*n_zwin,form='binary')
      write(20,rec=1) B2_width(1:n_zwin)
    close(20)

    FILE_NAME='BP_result_x3/B2.out'
    open(20,file=FILE_NAME,access='direct',recl=4*n_zwin,form='binary')
      write(20,rec=1) B2(1:n_zwin)
    close(20)

    deallocate(xs)
    PRINT *,"**end SWK_LMW2**"
!LMW2 reconstruct



!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FFTPACK5.1_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

!  stop 0
1 format(1x,5000e15.6)  
2 format(1x,2000e15.6)  
3 format(1x,10e15.6)
4 format(1x,20005e15.6)
    end





    subroutine SWK_LMW2_FUN(q_kt,dadt,dbdt,T_win,n_win1,n_win2,e_wk_lmw)
      implicit none
      integer::n_win1,n_win2,wtid,wt_id,itt,ik
      complex(kind=4)::q_kt(1:n_win1,1:n_win2)
      real(kind=4)::dadt(1:n_win1,1:n_win2),dbdt(1:n_win1,1:n_win2)
      real(kind=4)::T_win,w0_lmw,pi0,e_tmp,wt_p,wt_n,a
      real(kind=4)::e_wk_LMW(1:n_win1,1:n_win2)
      pi0=3.1415926535
      w0_lmw=2*pi0/T_win
      do  itt=1,n_win2
        do ik=1,n_win1
          a=abs(q_kt(ik,itt))
          e_tmp=a*a
          wt_p=dbdt(ik,itt)+dadt(ik,itt)/a
          if(wt_p .ge. 0) then
            wtid=int(abs(wt_p/w0_lmw)+0.5)
          else
            wtid=-(int(abs(wt_p/w0_lmw)+0.5))
          endif

          if((wtid .ge. n_win2/2) .or. (wtid .lt. -n_win2/2)) then
            e_tmp=0.0
            CYCLE
          elseif(wtid .lt. 0) then
            wt_id=n_win2+1-abs(wtid)
          else
            wt_id=wtid+1
          endif
          e_wk_LMW(ik,wt_id)=e_wk_LMW(ik,wt_id)+e_tmp
        enddo
      enddo

      do  itt=1,n_win2
        do ik=1,n_win1
          a=abs(q_kt(ik,itt))
          e_tmp=a*a
          wt_n=dbdt(ik,itt)-dadt(ik,itt)/a
          if(wt_n .ge. 0) then
            wtid=int(abs(wt_n/w0_lmw)+0.5)
          else
            wtid=-(int(abs(wt_n/w0_lmw)+0.5))
          endif
          if((wtid .ge. n_win2/2) .or. (wtid .lt. -n_win2/2)) then
            e_tmp=0.0
            CYCLE
          elseif(wtid .lt. 0) then
            wt_id=n_win2+1-abs(wtid)
          else
            wt_id=wtid+1
          endif
          e_wk_LMW(ik,wt_id)=e_wk_LMW(ik,wt_id)+e_tmp
        enddo
      enddo

      return
    end subroutine




    subroutine SWK_LMW1_FUN(q_xw,dadx,dbdx,L_win,n_win1,n_win2,e_wk_lmw)
      implicit none
      integer::n_win1,n_win2,kxid,kx_id,ixx,iw
      complex(kind=4)::q_xw(1:n_win1,1:n_win2)
      real(kind=4)::dadx(1:n_win1,1:n_win2),dbdx(1:n_win1,1:n_win2)
      real(kind=4)::L_win,k0_lmw,pi0,e_tmp,kx_p,kx_n,a
      real(kind=4)::e_wk_LMW(1:n_win1,1:n_win2)
      pi0=3.1415926535
      k0_lmw=2*pi0/L_win
      do  Iw=1,n_win2
        do ixx=1,n_win1
          a=abs(q_xw(ixx,iw))
          e_tmp=a*a
          kx_p=dbdx(ixx,iw)+dadx(ixx,iw)/a
          if(kx_p .ge. 0) then
            kxid=int(abs(kx_p/k0_lmw)+0.5)
          else
            kxid=-(int(abs(kx_p/k0_lmw)+0.5))
          endif

          if((kxid .ge. n_win1/2) .or. (kxid .lt. -n_win1/2)) then
            e_tmp=0.0
            CYCLE
          elseif(kxid .lt. 0) then
            kx_id=n_win1+1-abs(kxid)
          else
            kx_id=kxid+1
          endif
          e_wk_LMW(kx_id,iw)=e_wk_LMW(kx_id,iw)+e_tmp
        enddo
      enddo

      do  iw=1,n_win2
        do ixx=1,n_win1
          a=abs(q_xw(ixx,iw))
          e_tmp=a*a
          kx_n=dbdx(ixx,iw)-dadx(ixx,iw)/a
          if(kx_n .ge. 0) then
            kxid=int(abs(kx_n/k0_lmw)+0.5)
          else
            kxid=-(int(abs(kx_n/k0_lmw)+0.5))
          endif
          if((kxid .ge. n_win1/2) .or. (kxid .lt. -n_win1/2)) then
            e_tmp=0.0
            CYCLE
          elseif(kxid .lt. 0) then
            kx_id=n_win1+1-abs(kxid)
          else
            kx_id=kxid+1
          endif
          e_wk_LMW(kx_id,iw)=e_wk_LMW(kx_id,iw)+e_tmp
        enddo
      enddo

      return
    end subroutine
    
       Subroutine itoa6(intnum, chars)
        Character(len=1) chd1, chd2, chd3, chd4,chd5,chd6
        Character(len=6) chars
        integer(kind=4) intnum
        id1   = intnum/100000
        id2   = (intnum - id1*100000)/10000
        id3   = (intnum - id1*100000 - id2*10000)/1000
        id4   = (intnum - id1*100000 - id2*10000 - id3*1000)/100
        id5=(intnum-id1*100000-id2*10000-id3*1000-id4*100)/10
        id6 = (intnum-id1*100000-id2*10000-id3*1000-id4*100-id5*10)/1
        chd1  = CHAR(id1+48)
        chd2  = CHAR(id2+48)
        chd3  = CHAR(id3+48)
        chd4  = CHAR(id4+48)
        chd5  = CHAR(id5+48)
        chd6  = CHAR(id6+48)
        chars = chd1//chd2//chd3//chd4//chd5//chd6
        Return
    End

    Subroutine itoa4(intnum, chars)
      Character(len=1) chd1, chd2, chd3, chd4
      Character(len=4) chars
      integer(kind=4) intnum
      id1   = intnum/1000
      id2   = (intnum - id1*1000)/100
      id3   = (intnum - id1*1000 - id2*100)/10
      id4   = (intnum - id1*1000- id2*100 - id3*10)/1
      chd1  = CHAR(id1+48)
      chd2  = CHAR(id2+48)
      chd3  = CHAR(id3+48)
      chd4  = CHAR(id4+48)
      chars = chd1//chd2//chd3//chd4
      Return
  End

    subroutine window1(signal,N_orig,num_win,n_win,signal_win)
    implicit none
    integer::N_orig,iwin,num_win,n_win
    integer::i,ii
    real::pi
    real::signal(1:N_orig),signal_win(1:num_win,1:n_win),win_fun(1:n_win)

    !number of points for each window
    pi=3.1415926535
    do ii=1,n_win
      win_fun(ii)=0.5-0.5*cos(2*pi*(real(ii)-1)/(real(n_win)-1))
    enddo


    do iwin=1,num_win
!      signal_win(iwin,1:N_win)=signal( (iwin-1)*n_win/2+1 : (iwin-1)*n_win/2+n_win )*win_fun(1:n_win)
      signal_win(iwin,1:N_win)=signal( (iwin-1)*n_win/2+1 : (iwin-1)*n_win/2+n_win )
    enddo
        return
    end

    subroutine windowt(signal,N_orig,num_win,n_win,signal_win,i_method)
      implicit none
      integer::N_orig,iwin,num_win,n_win
      integer::i,ii,i_method,k,win_count
      real::pi
      real::signal(1:N_orig),signal_win(1:num_win,1:n_win),win_fun1(1:n_win),win_fun2(1:n_win)
  
      !number of points for each window
      pi=3.1415926535
      if(i_method.eq.4) then
        do ii=1,n_win
          win_fun1(ii)=0.5-0.5*cos(2*pi*(real(ii)-1)/(real(n_win)-1))
          win_fun2(ii)=0.54-0.46*cos(2*pi*(real(ii)-1)/(real(n_win)-1))
        enddo
      elseif(i_method.eq.1) then
        open(10,file='Case4/Blackmanharris_wint.dat')
        do k=1,n_win
          read(10,*) win_fun2(k)
        enddo
        close(10)
      elseif(i_method.eq.2) then
        open(10,file='Case4/Kaiser_wint.dat')
        do k=1,n_win
          read(10,*) win_fun2(k)
        enddo
        close(10)
      elseif(i_method.eq.3) then
        open(10,file='Case4/Flattop_wint.dat')
        do k=1,n_win
          read(10,*) win_fun2(k)
        enddo
        close(10)
      elseif(i_method.eq.5) then
          do k=1,n_win
          win_fun2(k)=1.0
          enddo
      endif

      do iwin=1,num_win
        signal_win(iwin,1:N_win)=signal( (iwin-1)*n_win/2+1 : (iwin-1)*n_win/2+n_win )*win_fun2(1:n_win)
      enddo
          return
    end subroutine

    subroutine windowz(signal,N_orig,num_win,n_win,signal_win,i_method)
      implicit none
      integer::N_orig,iwin,num_win,n_win
      integer::i,ii,i_method,k,win_count
      real::pi
      real::signal(1:N_orig),signal_win(1:num_win,1:n_win),win_fun1(1:n_win),win_fun2(1:n_win)
  
      !number of points for each window
      pi=3.1415926535
      if(i_method.eq.4) then
        do ii=1,n_win
          win_fun1(ii)=0.5-0.5*cos(2*pi*(real(ii)-1)/(real(n_win)-1))
          win_fun2(ii)=0.54-0.46*cos(2*pi*(real(ii)-1)/(real(n_win)-1))
        enddo
      elseif(i_method.eq.1) then
        open(10,file='Case4/Blackmanharris_winz.dat')
        do k=1,n_win
          read(10,*) win_fun2(k)
        enddo
        close(10)
      elseif(i_method.eq.2) then
        open(10,file='Case4/Kaiser_winz.dat')
        do k=1,n_win
          read(10,*) win_fun2(k)
        enddo
        close(10)
      elseif(i_method.eq.3) then
        open(10,file='Case4/Flattop_winz.dat')
        do k=1,n_win
          read(10,*) win_fun2(k)
        enddo
        close(10)
      elseif(i_method.eq.5) then
        do k=1,n_win
        win_fun2(k)=1.0
        enddo
      endif

  
      do iwin=1,num_win
 !       signal_win(iwin,1:N_win)=signal( (iwin-1)*n_win/2+1 : (iwin-1)*n_win/2+n_win )*win_fun1(1:n_win)
        signal_win(iwin,1:N_win)=signal( (iwin-1)*n_win/2+1 : (iwin-1)*n_win/2+n_win )*win_fun2(1:n_win)
  !      signal_win(iwin,1:N_win)=signal( (iwin-1)*n_win/2+1 : (iwin-1)*n_win/2+n_win )
      enddo

!c      if(win_count.eq.1) then
!c      open(10,file='Case2/win_check.dat')
!c      do k=1,n_win
!c        write(10,6) win_fun2(k)
!c      enddo
!c      close(10)
!c      endif

          return
6 format(1x,5e15.6)
    end subroutine
    
    subroutine FFT1D(signal,N,signal_fft)
    
    implicit none
    integer::N
    integer::i,ib,ier,inc,lenc
    complex::signal(1:N)
    complex::signal_fft(1:n)

  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk
  real ( kind = 4 ), allocatable, dimension ( : ) :: work
  real ( kind = 4 ), allocatable, dimension ( : ) :: wsave

  signal_fft=(0.0,0.0)
          lenwrk = 2 * N
          lensav = 2 * N + int ( log ( real ( N, kind = 4 ) ) / log ( 2.0E+00 ) ) + 4
  
     allocate(work(1:lenwrk))
     allocate(wsave(1:lensav))
     call cfft1i(N,wsave,lensav,ier)
     inc=1
     lenc=N
       call cfft1f(N,inc,signal(1:N),lenc,wsave,lensav,work,lenwrk,ier)
       signal_fft(1:N)=signal(1:N)
       
    deallocate ( work )
    deallocate ( wsave )
        return
    end

    subroutine IFFT1D(signal_fft,N,signal)
          
          implicit none
          integer::n_block,n_overlap,N,n_eve,nn
          integer::i,ib,iwin,ier,inc,lenc
          complex::signal(1:N)
          complex::signal_fft(1:N)
      
        integer ( kind = 4 ) lensav
        integer ( kind = 4 ) lenwrk
        real ( kind = 4 ), allocatable, dimension ( : ) :: work
        real ( kind = 4 ), allocatable, dimension ( : ) :: wsave
              
                lenwrk = 2 * N
                lensav = 2 * N + int ( log ( real ( N, kind = 4 ) ) / log ( 2.0E+00 ) ) + 4
        
           allocate(work(1:lenwrk))
           allocate(wsave(1:lensav))
           call cfft1i(N,wsave,lensav,ier)
           inc=1
           lenc=N
             call cfft1b(N,inc,signal_fft(1:N),lenc,wsave,lensav,work,lenwrk,ier)
             signal(1:N)=signal_fft(1:N)
             
          deallocate ( work )
          deallocate ( wsave )
              return
          end
    
          subroutine fft2S(RU,N1,N2,RC)
            implicit none
            INTEGER::N1,N2,INUM,IC
            INTEGER::lenwrk,lensav,ier,ldim
            COMPLEX::RU(1:N1,1:N2),R(1:N1,1:N2),RC(1:N1,1:N2)
            real ( kind = 4 ), allocatable, dimension ( : ) :: work
            real ( kind = 4 ), allocatable, dimension ( : ) :: wsave
            
            R=RU
            
            LENWRK=2*N1*N2
            lensav = 2 * n1 + int ( log ( real ( n1, kind = 4 ) ) / log ( 2.0E+00 ) ) &
            + 2 * n2 + int ( log ( real ( n2, kind = 4 ) ) / log ( 2.0E+00 ) ) &
            + 8 
              allocate ( work(1:lenwrk) )
              allocate ( wsave(1:lensav) )
                call cfft2i ( N1, N2, wsave, lensav, ier )
                LDIM=N1
               call cfft2F ( ldim, N1, N2, R(1:N1,1:N2), wsave, lensav, work, lenwrk, ier )   
        !       RU(0:MX-1,0:MZ-1,IC)=REAL(CU(0:MX-1,0:MZ-1,IC))
                DEALLOCATE(WORK)
                DEALLOCATE(WSAVE)
            RC=R
        
            RETURN
            end

            subroutine Ifft2S(RC,N1,N2,RU)
              implicit none
              INTEGER::N1,N2,INUM,IC
              INTEGER::lenwrk,lensav,ier,ldim
              COMPLEX::RU(1:N1,1:N2),R(1:N1,1:N2),RC(1:N1,1:N2)
              real ( kind = 4 ), allocatable, dimension ( : ) :: work
              real ( kind = 4 ), allocatable, dimension ( : ) :: wsave
              
              R=RC
              
              LENWRK=2*N1*N2
              lensav = 2 * n1 + int ( log ( real ( n1, kind = 4 ) ) / log ( 2.0E+00 ) ) &
              + 2 * n2 + int ( log ( real ( n2, kind = 4 ) ) / log ( 2.0E+00 ) ) &
              + 8 
                allocate ( work(1:lenwrk) )
                allocate ( wsave(1:lensav) )
                  call cfft2i ( N1, N2, wsave, lensav, ier )
                  LDIM=N1
                 call cfft2B ( ldim, N1, N2, R(1:N1,1:N2), wsave, lensav, work, lenwrk, ier )   
          !       RU(0:MX-1,0:MZ-1,IC)=REAL(CU(0:MX-1,0:MZ-1,IC))
                  DEALLOCATE(WORK)
                  DEALLOCATE(WSAVE)
              RU=R
          
              RETURN
              end
    
    
    
    
    
    
    
subroutine c4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! C4MAT_PRINT_SOME prints some of a C4MAT.
!
!  Discussion:
!
!    A C4MAT is a matrix of C4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the matrix.
!
!    Input, complex ( kind = 4 ) A(M,N), the matrix.
!
!    Input, integer ( kind = 4 ) ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 4
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 4 ) a(m,n)
  character ( len = 20 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * )  title
  complex ( kind = 4 ) zero

  zero = cmplx ( 0.0E+00, 0.0E+00, kind = 4 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)' 
    return
  end if
!
!  Print the columns of the matrix, in strips of INCX.
!
  do j2lo = jlo, min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == zero ) then
          ctemp(j2) = '       0.0          '
        else if ( imag ( a(i,j) ) == 0.0E+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( a(i,j), kind = 4 )
        else
          write ( ctemp(j2), '(2g10.3)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a1,4a20)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine c4mat_uniform_01 ( m, n, seed, c )

!*****************************************************************************80
!
!! C4MAT_UNIFORM_01 returns a unit pseudorandom C4MAT.
!
!  Discussion:
!
!    A C4MAT is a matrix of C4's.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 4 ) C(M,N), the pseudorandom complex matrix.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  complex ( kind = 4 ) c(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  integer ( kind = 4 ) seed
  real ( kind = 4 ) theta

  do j = 1, n
    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r = sqrt ( real ( seed, kind = 4 ) * 4.656612875E-10 )

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      theta = 2.0E+00 * r4_pi * ( real ( seed, kind = 4 ) * 4.656612875E-10 )

      c(i,j) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 4 )

    end do

  end do

  return
end
subroutine c4vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! C4VEC_PRINT_PART prints "part" of a C4VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, complex ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) n

  complex   ( kind = 4 ) a(n)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) max_print
  character ( len = * )  title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ........  ..............  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,a)' ) i, ':', a(i), &
      '...more entries...'

  end if

  return
end
subroutine c4vec_uniform_01 ( n, seed, c )

!*****************************************************************************80
!
!! C4VEC_UNIFORM_01 returns a unit pseudorandom C4VEC.
!
!  Discussion:
!
!    A C4VEC is a vector of C4's.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of values to compute.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which should 
!    NOT be 0.  On output, SEED has been updated.
!
!    Output, complex ( kind = 4 ) C(N), the pseudorandom complex vector.
!
  implicit none

  integer ( kind = 4 ) n

  complex ( kind = 4 ) c(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  real ( kind = 4 ), parameter :: r4_pi = 3.141592653589793E+00
  integer ( kind = 4 ) seed
  real ( kind = 4 ) theta

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r = sqrt ( real ( seed, kind = 4 ) * 4.656612875E-10 )

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    theta = 2.0E+00 * r4_pi * ( real ( seed, kind = 4 ) * 4.656612875E-10 )

    c(i) = r * cmplx ( cos ( theta ), sin ( theta ), kind = 4 )

  end do

  return
end
subroutine r4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R4MAT_PRINT prints an R4MAT.
!
!  Discussion:
!
!    An R4MAT is an array of R4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 4 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  character ( len = * )  title

  call r4mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R4MAT_PRINT_SOME prints some of an R4MAT.
!
!  Discussion:
!
!    An R4MAT is an array of R4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * )  title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 4 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r4mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R4MAT_UNIFORM_01 fills an R4MAT with unit pseudorandom numbers.
!
!  Discussion:
!
!    An R4MAT is an array of R4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 May 2011
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns in
!    the array.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 4 ) r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = real ( seed, kind = 4 ) * 4.656612875E-10

    end do
  end do

  return
end
subroutine r4vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! R4VEC_PRINT_PART prints "part" of an R4VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries of the vector.
!
!    Input, real ( kind = 4 ) A(N), the vector to be printed.
!
!    Input, integer ( kind = 4 ) MAX_PRINT, the maximum number of lines
!    to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer   ( kind = 4 ) n

  real      ( kind = 4 ) a(n)
  integer   ( kind = 4 ) i
  integer   ( kind = 4 ) max_print
  character ( len = * )  title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ........  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,a)' ) i, ':', a(i), '...more entries...'

  end if

  return
end
subroutine r4vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R4VEC_UNIFORM_01 returns a unit pseudorandom R4VEC.
!
!  Discussion:
!
!    An R4VEC is an array of R4's.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value,
!    which should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 4 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = 4 ) * 4.656612875E-10

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8  ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9  ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end                          