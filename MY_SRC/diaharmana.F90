MODULE harmonic_analysis
   !!======================================================================
   !!                       ***  MODULE  example  ***
   !! Ocean physics:  On line harmonic analyser
   !!                 
   !!=====================================================================

#if defined key_harm_ana

   !!----------------------------------------------------------------------
   !!   'key_harm_ana'  :                Calculate harmonic analysis
   !!----------------------------------------------------------------------
   !!   harm_ana        :
   !!   harm_ana_init   :
   !!   NB: 2017-12 : add 3D harmonic analysis of velocities
   !!                 integration of Maria Luneva's development
   !!   'key_3Ddiaharm'
   !!----------------------------------------------------------------------

   USE oce             ! ocean dynamics and tracers 
   USE dom_oce         ! ocean space and time domain
   USE iom
   USE in_out_manager  ! I/O units
   USE phycst          ! physical constants
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE bdy_oce         ! ocean open boundary conditions
   USE bdytides        ! tidal bdy forcing
   USE daymod          ! calendar
   USE tideini
   USE restart
   USE ioipsl, ONLY : ju2ymds    ! for calendar
   !
   !
   USE timing          ! preformance summary
# if defined key_3Ddiaharm
   USE zdf_oce
#endif

   IMPLICIT NONE
   PRIVATE

   !! *  Routine accessibility
   PUBLIC harm_ana    ! routine called in step.F90 module

   !! * Module variables
   INTEGER, PARAMETER ::  nharm_max  = jpmax_harmo  ! max number of harmonics to be analysed 
   INTEGER, PARAMETER ::  nhm_max    = 2*nharm_max+1 
   INTEGER, PARAMETER ::  nvab       = 2 ! number of 3D variables
   INTEGER            ::  nharm
   INTEGER            ::  nhm 
   INTEGER ::                 & !!! ** toto namelist (namtoto) **
      nflag  =  1                ! default value of nflag 
   REAL(wp), DIMENSION(nharm_max) ::                & 
      om_tide                     ! tidal frequencies ( rads/sec)
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:)   ::                & 
      bzz,c,x    ! work arrays
   REAL(wp) :: cca,ssa,zm,bt
   REAL(wp) :: ccau,ccav,ssau,ssav
   REAL(wp), PUBLIC, ALLOCATABLE,DIMENSION(:) :: anau, anav, anaf   ! nodel/phase corrections used by diaharmana
!
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:,:,:) ::   &
      bssh, bubar, bvbar        ! work array for ssh anaylsis
!
   REAL(wp), ALLOCATABLE,SAVE, DIMENSION(:,:,:,:) :: cosamp2D, sinamp2D
!
# if defined key_3Ddiaharm
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:,:,:)     :: ana_bfrU, ana_bfrV
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:,:,:,:)   :: ana_D, ana_I, ana_U, ana_V, ana_W
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:,:,:,:,:) :: cosamp3D, sinamp3D  
#endif
   REAL(WP), ALLOCATABLE,SAVE,DIMENSION(:,:) :: cc,a
!
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:,:) ::   &
      gout,hout        ! arrays for output
# if defined key_3Ddiaharm
   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:,:,:) ::   &
      gout3D,hout3D        ! arrays for 3D output
#endif
!
!   REAL(wp), ALLOCATABLE,SAVE,DIMENSION(:,:) ::   &
!      huout,hvout,guout,gvout        ! arrays for output
   REAL(wp), PUBLIC ::   fjulday_startharm       !: Julian Day since start of harmonic analysis

!  NAMELIST
   INTEGER ::   nit000_han    ! First time step used for harmonic analysis
   INTEGER ::   nitend_han    ! Last time step used for harmonic analysis
   INTEGER ::   nstep_han     ! Time step frequency for harmonic analysis
   INTEGER ::   nb_ana        ! Number of harmonics to analyse
   CHARACTER (LEN=4), DIMENSION(jpmax_harmo) ::   tname   ! Names of tidal constituents ('M2', 'K1',...)

   INTEGER , ALLOCATABLE, DIMENSION(:)       ::   ntide_all ! INDEX within the full set of constituents (tide.h90)
   INTEGER , ALLOCATABLE, DIMENSION(:)       ::   ntide_sub ! INDEX within the subset of constituents pass in input

   !! * Substitutions

   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! or LIM 2.0 , UCL-LOCEAN-IPSL (2005)
   !! or  TOP 1.0 , LOCEAN-IPSL (2005)
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/module_example,v 1.3 2005/03/27 18:34:47 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE harm_ana( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE harm_ana  ***
      !!
      !! ** Purpose :   Harmonic analyser
      !!
      !! ** Method  :   
      !!
      !! ** Action  : - first action (share memory array/varible modified
      !!                in this routine
      !!              - second action .....
      !!              - .....
      !!
      !! References :
      !!   Give references if exist otherwise suppress these lines
      !!
      !! History :
      !!   9.0  !  03-08  (Autor Names)  Original code
      !!        !  02-08  (Author names)  brief description of modifications
      !!----------------------------------------------------------------------
      !! * Modules used
      
      !! * arguments
      INTEGER, INTENT( in  ) ::   &  
         kt                          ! describe it!!!

      !! * local declarations
      INTEGER ::   ji, jj, jk,  &        ! dummy loop arguments
                   ih,i1,i2,iv,jgrid

      !!--------------------------------------------------------------------

      IF( nn_timing == 1 )   CALL timing_start('harm_ana')

      IF( kt == nit000  )   CALL harm_ana_init    ! Initialization (first time-step only)

      ! this bit done every time step
      nhm=2*nb_ana+1
      c(1) = 1.0


      IF(lwp) WRITE(numout,*) "ztime NEW", kt, (fjulday-fjulday_startharm)*86400._wp,  sshn(314,264)


      DO ih=1,nb_ana
         c(2*ih  ) = anaf(ih)*cos((fjulday-fjulday_startharm)*86400._wp*om_tide(ih) +anau(ih)+anav(ih) )
         c(2*ih+1) = anaf(ih)*sin((fjulday-fjulday_startharm)*86400._wp*om_tide(ih) +anau(ih)+anav(ih) )
      ENDDO 

      ! CUMULATES
      DO ji=1,jpi
         DO jj=1,jpj
            DO ih=1,nhm
               bssh (ih,ji,jj) = bssh (ih,ji,jj) + c(ih) * sshn(ji,jj) * ssmask (ji,jj)
               bubar(ih,ji,jj) = bubar(ih,ji,jj) + c(ih) * un_b(ji,jj) * ssumask(ji,jj)
               bvbar(ih,ji,jj) = bvbar(ih,ji,jj) + c(ih) * vn_b(ji,jj) * ssvmask(ji,jj)
# if defined key_3Ddiaharm
               ana_bfrU(ih,ji,jj) = ana_bfrU(ih,ji,jj) + c(ih) * bfrua(ji,jj) * un(ji,jj,mbku(ji,jj)) * ssumask(ji,jj)
               ana_bfrV(ih,ji,jj) = ana_bfrV(ih,ji,jj) + c(ih) * bfrva(ji,jj) * vn(ji,jj,mbkv(ji,jj)) * ssvmask(ji,jj)
               DO jk=1,jpkm1     
                  ana_D(ih,ji,jj,jk) = ana_D(ih,ji,jj,jk) + c(ih) *  rhd(ji,jj,jk)               * tmask(ji,jj,jk)
                  ana_U(ih,ji,jj,jk) = ana_U(ih,ji,jj,jk) + c(ih) * ( un(ji,jj,jk)-un_b(ji,jj) ) * umask(ji,jj,jk)
                  ana_V(ih,ji,jj,jk) = ana_V(ih,ji,jj,jk) + c(ih) * ( vn(ji,jj,jk)-vn_b(ji,jj) ) * vmask(ji,jj,jk)
                  ana_W(ih,ji,jj,jk) = ana_W(ih,ji,jj,jk) + c(ih) *   wn(ji,jj,jk)               * wmask(ji,jj,jk)
                  !IF (kk.le.jpkm1) 
                  ana_I(ih,ji,jj,jk) = ana_I(ih,ji,jj,jk) + 0.5*grav*c(ih)*(rhd(ji,jj,jk)+rhd(ji,jj,jk+1) )/max(rn2(ji,jj,jk),1.e-8_wp)* tmask(ji,jj,jk)
!                  IF ( jk <= mbathy(ji,jj) ) THEN
!                     ana_I(ih,ji,jj,jk) = ana_I(ih,ji,jj,jk) + 0.5*grav*c(ih)*(rhd(ji,jj,jk)+rhd(ji,jj,jk+1) )/max(rn2(ji,jj,jk),1.e-8_wp)
!                  ENDIF
               ENDDO
#endif
            ENDDO
         ENDDO
      ENDDO

      ! Compute nodal factor cumulatve cross-product
      DO i1=1,nhm
         DO i2=1,nhm
            cc(i1,i2)=cc(i1,i2)+c(i1)*c(i2)
         ENDDO
      ENDDO

      ! Output RESTART
      IF( kt == nitrst ) THEN
         CALL harm_rst_write(kt) ! Dump out data for a restarted run 
      ENDIF

      ! At End of run
      IF ( kt ==  nitend ) THEN

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'harm_ana : harmonic analysis of tides at end of run'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~'

         IF( ln_harm_ana_compute ) THEN
             IF(lwp) WRITE(numout,*) "Computing harmonics at last step"

             ! INITIALISE TABLE TO 0
             cosamp2D = 0.0_wp
             sinamp2D = 0.0_wp

             ! FIRST OUTPUT 2D VARIABLES
# if defined key_3Ddiaharm
             DO jgrid=1,5 ! elevation, Ubar, Vbar, Surface friction U and V
# else
             DO jgrid=1,3 ! elevation, Ubar, Vbar
# endif
                DO ji=1,jpi
                   DO jj=1,jpj
                      bt = 1.0_wp; bzz(:) = 0.0_wp
                      DO ih=1,nhm
                         IF( jgrid .eq. 1 ) bzz(ih) = bssh (ih,ji,jj)
                         IF( jgrid .eq. 2 ) bzz(ih) = bubar(ih,ji,jj)
                         IF( jgrid .eq. 3 ) bzz(ih) = bvbar(ih,ji,jj)
# if defined key_3Ddiaharm
                         IF( jgrid .eq. 4 ) bzz(ih) = ana_bfrU(ih,ji,jj)
                         IF( jgrid .eq. 5 ) bzz(ih) = ana_bfrV(ih,ji,jj)
# endif
                         bt = bt*bzz(ih)
                      ENDDO

                      ! Copy back original cumulated nodal factor
                      a(:,:) = cc(:,:)

!                     now do gaussian elimination of the system
!                     a * x = b
!                     the matrix x is (a0,a1,b1,a2,b2 ...)
!                     the matrix a and rhs b solved here for x
                      x=0.0_wp
                      IF(bt.ne.0.) THEN
                        CALL gelim( a, bzz, x, nhm )
!                       Backup output in variables
                        DO ih=1,nb_ana
                           cosamp2D(ih,ji,jj,jgrid) = x(ih*2  )
                           sinamp2D(ih,ji,jj,jgrid) = x(ih*2+1)
                        ENDDO
                        cosamp2D( 0,ji,jj,jgrid) = x(1)
                        sinamp2D( 0,ji,jj,jgrid) = 0.0_wp
                      ENDIF     ! bt.ne.0.
                   ENDDO        ! jj
                ENDDO           ! ji
             ENDDO              ! jgrid

# if defined key_3Ddiaharm
             ! INITIALISE TABLE TO 0
             cosamp3D = 0.0_wp
             sinamp3D = 0.0_wp
             ! SECOND OUTPUT 3D VARIABLES
             DO jgrid=1,5 ! dzi, rho, U, V, W
                DO jk=1,jpkm1
                   DO ji=1,jpi
                      DO jj=1,jpj
                         bt = 1.0_wp; bzz(:) = 0.0_wp
                         DO ih=1,nhm
                            IF( jgrid .eq. 1 ) bzz(ih) = ana_I(ih,ji,jj,jk)
                            IF( jgrid .eq. 2 ) bzz(ih) = ana_D(ih,ji,jj,jk)
                            IF( jgrid .eq. 3 ) bzz(ih) = ana_U(ih,ji,jj,jk)
                            IF( jgrid .eq. 4 ) bzz(ih) = ana_V(ih,ji,jj,jk)
                            IF( jgrid .eq. 5 ) bzz(ih) = ana_W(ih,ji,jj,jk)
                            bt = bt*bzz(ih)
                         ENDDO

                         ! Copy back original cumulated nodal factor
                         a(:,:) = cc(:,:)                      

!                        now do gaussian elimination of the system
!                        a * x = b
!                        the matrix x is (a0,a1,b1,a2,b2 ...)
!                        the matrix a and rhs b solved here for x
                         x=0.0_wp
                         IF(bt.ne.0.) THEN
                           CALL gelim( a, bzz, x, nhm )
!                          Backup output in variables
                           DO ih=1,nb_ana
                              cosamp3D(ih,ji,jj,jk,jgrid) = x(ih*2  )
                              sinamp3D(ih,ji,jj,jk,jgrid) = x(ih*2+1)
                           ENDDO
                           cosamp3D   ( 0,ji,jj,jk,jgrid) = x(1)
                           sinamp3D   ( 0,ji,jj,jk,jgrid) = 0.0_wp
                        ENDIF     ! bt.ne.0.
                      ENDDO       ! jj
                   ENDDO          ! ji
                ENDDO             ! jk
             ENDDO                ! jgrid
# endif

             CALL harm_ana_out     ! output analysis (last time step)

         ELSE    ! ln_harmana_compute = False 
             IF(lwp) WRITE(numout,*) " Skipping Computing harmonics at last step"

         ENDIF   ! ln_harmana_compute 
      ENDIF      ! kt ==  nitend

      IF( nn_timing == 1 )   CALL timing_stop('harm_ana')

   END SUBROUTINE harm_ana

   SUBROUTINE harm_ana_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE harm_ana_init  ***
      !!                   
      !! ** Purpose :   initialization of ....
      !!
      !! ** Method  :   blah blah blah ...
      !!
      !! ** input   :   Namlist namexa
      !!
      !! ** Action  :   ...  
      !!
      !! history :
      !!   9.0  !  03-08  (Autor Names)  Original code
      !!----------------------------------------------------------------------
      !! * local declarations
      INTEGER ::   ji, jj, jk, jit,ih   ! dummy loop indices

      INTEGER ::   ios                 ! Local integer output status for namelist read
      NAMELIST/nam_diaharm/ nit000_han, nitend_han, nstep_han, tname

      !!----------------------------------------------------------------------
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'harm_init : initialization of harmonic analysis of tides'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~'

      ! GET NAMELIST DETAILS
      REWIND( numnam_ref )              ! Namelist nam_diaharm in reference namelist : Tidal harmonic analysis
      READ  ( numnam_ref, nam_diaharm, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diaharm in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist nam_diaharm in configuration namelist : Tidal harmonic analysis
      READ  ( numnam_cfg, nam_diaharm, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diaharm in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nam_diaharm )

      ! GET NUMBER OF HARMONIC TO ANALYSE - from diaharm.F90
      nb_ana = 0
      DO jk=1,jpmax_harmo
         DO ji=1,nb_harmo
            IF(TRIM(tname(jk)) == Wave( ntide(ji) )%cname_tide ) THEN
               nb_ana=nb_ana+1
            ENDIF
         END DO
      END DO
      !
      IF(lwp) THEN
         WRITE(numout,*) '        Namelist nam_diaharm'
         WRITE(numout,*) '        nb_ana    = ', nb_ana
         CALL flush(numout)
      ENDIF
      !
      IF (nb_ana > nharm_max) THEN
        IF(lwp) WRITE(numout,*) ' E R R O R harm_ana : nb_ana must be lower than nharm_max, stop'
        IF(lwp) WRITE(numout,*) ' nharm_max = ', nharm_max
        nstop = nstop + 1
      ENDIF

      ALLOCATE(ntide_all(nb_ana))
      ALLOCATE(ntide_sub(nb_ana))

      DO jk=1,nb_ana
       DO ji=1,nb_harmo
          IF (TRIM(tname(jk)) .eq. Wave( ntide(ji) )%cname_tide ) THEN
             ntide_sub(jk) = ji
             ntide_all(jk) = ntide(ji)
             EXIT
          END IF
       END DO
      END DO

      ! DO ALLOCATIONS
      ALLOCATE( bubar(nb_ana*2+1,jpi,jpj) )
      ALLOCATE( bvbar(nb_ana*2+1,jpi,jpj) )
      ALLOCATE( bssh (nb_ana*2+1,jpi,jpj) )
# if defined key_3Ddiaharm
      ALLOCATE( ana_D (nb_ana*2+1,jpi,jpj,jpk) )
      ALLOCATE( ana_I (nb_ana*2+1,jpi,jpj,jpk) )
      ALLOCATE( ana_U (nb_ana*2+1,jpi,jpj,jpk) )
      ALLOCATE( ana_V (nb_ana*2+1,jpi,jpj,jpk) )
      ALLOCATE( ana_W (nb_ana*2+1,jpi,jpj,jpk) )
      ALLOCATE( ana_bfrU (nb_ana*2+1,jpi,jpj) )
      ALLOCATE( ana_bfrV (nb_ana*2+1,jpi,jpj) )

      ALLOCATE( cosamp2D(0:nb_ana*2+1,jpi,jpj,5))
      ALLOCATE( sinamp2D(0:nb_ana*2+1,jpi,jpj,5))
      ALLOCATE( cosamp3D(0:nb_ana*2+1,jpi,jpj,jpk,5))
      ALLOCATE( sinamp3D(0:nb_ana*2+1,jpi,jpj,jpk,5))

      ALLOCATE( gout3D (jpi,jpj,jpk) )
      ALLOCATE( hout3D (jpi,jpj,jpk) )

# else
      ALLOCATE( cosamp2D(0:nb_ana*2+1,jpi,jpj,3))
      ALLOCATE( sinamp2D(0:nb_ana*2+1,jpi,jpj,3))

# endif

      ALLOCATE( cc(nb_ana*2+1,nb_ana*2+1) )
      ALLOCATE( a (nb_ana*2+1,nb_ana*2+1) )

      ALLOCATE( anau(nb_ana) )
      ALLOCATE( anav(nb_ana) )
      ALLOCATE( anaf(nb_ana) )

      ALLOCATE( bzz(nb_ana*2+1) )
      ALLOCATE( x  (nb_ana*2+1) )
      ALLOCATE( c  (nb_ana*2+1) )

      ALLOCATE( gout (jpi,jpj) )
      ALLOCATE( hout (jpi,jpj) )
!      ALLOCATE( guout(jpi,jpj) )
!      ALLOCATE( huout(jpi,jpj) )
!      ALLOCATE( gvout(jpi,jpj) )
!      ALLOCATE( hvout(jpi,jpj) )
      ! END ALLOCATE 

      ! SELECT AND STORE FREQUENCIES
      IF(lwp)    WRITE(numout,*) 'Analysed frequency  : ',nb_ana ,'Frequency '
      DO ih=1,nb_ana
         om_tide(ih) = omega_tide( ntide_sub(ih) ) 
         IF(lwp) WRITE(numout,*) '                    : ',tname(ih),' ',om_tide(ih)
      ENDDO

      ! READ RESTART IF 
      IF ( ln_harmana_read ) THEN
         IF (lwp) WRITE(numout,*) "Reading previous harmonic data from previous run"
         ! Need to read in  bssh bz, cc anau anav and anaf 
         call harm_rst_read  ! This reads in from the previous day
                             ! Currrently the data in in assci format
      ELSE 

         IF (lwp) WRITE(numout,*) "Starting harmonic analysis from Fresh "
         bssh (:,:,:) = 0.0_wp
         bubar(:,:,:) = 0.0_wp
         bvbar(:,:,:) = 0.0_wp
# if defined key_3Ddiaharm
         ana_D (:,:,:,:) = 0.0_wp
         ana_I (:,:,:,:) = 0.0_wp
         ana_U (:,:,:,:) = 0.0_wp
         ana_V (:,:,:,:) = 0.0_wp
         ana_W (:,:,:,:) = 0.0_wp
         ana_bfrU(:,:,:) = 0.0_wp
         ana_bfrV(:,:,:) = 0.0_wp
# endif
         cc           = 0.0_wp
         anau (:)     = 0.0_wp
         anav (:)     = 0.0_wp
         anaf (:)     = 0.0_wp
         bzz  (:)     = 0.0_wp
         x    (:)     = 0.0_wp
         c    (:)     = 0.0_wp
         a    (:,:)   = 0.0_wp ! NB

         DO ih = 1, nb_ana
            anau(ih) = utide ( ntide_sub(ih) )
            anav(ih) = v0tide( ntide_sub(ih) )
            anaf(ih) = ftide ( ntide_sub(ih) )
         END DO

         fjulday_startharm=fjulday !Set this at very start and store

         IF (lwp) THEN
            WRITE(numout,*) '--------------------------'
            WRITE(numout,*) '   - Output anaf for check'
            WRITE(numout,*) 'ANA F', anaf
            WRITE(numout,*) 'ANA U', anau
            WRITE(numout,*) 'ANA V', anav
            WRITE(numout,*) fjulday_startharm
            WRITE(numout,*) '--------------------------'
         ENDIF
      ENDIF

   END SUBROUTINE harm_ana_init
!
   SUBROUTINE gelim (a,b,x,n)
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE harm_ana  ***
      !!
      !! ** Purpose :   Guassian elimination
      !!
      !!
      !! ** Action  : - first action (share memory array/varible modified
      !!                in this routine
      !!              - second action .....
      !!              - .....
      !!
      !! References :
      !!   Give references if exist otherwise suppress these lines
      !!
      !! History :
        implicit none
!
        integer  :: n
        REAL(WP) :: b(nb_ana*2+1), a(nb_ana*2+1,nb_ana*2+1)
        REAL(WP) :: x(nb_ana*2+1)
        INTEGER  :: row,col,prow,pivrow,rrow,ntemp
        REAL(WP) :: atemp
        REAL(WP) :: pivot
        REAL(WP) :: m
        REAL(WP) :: tempsol(nb_ana*2+1)

        do row=1,n-1
           pivrow=row
           pivot=a(row,n-row+1)
           do prow=row+1,n
              if (abs(a(prow,n-row+1)).gt.abs(pivot)  ) then
                 pivot=a(prow,n-row+1)
                 pivrow=prow
              endif
           enddo
!	swap row and prow
           if ( pivrow .ne. row ) then
              atemp=b(pivrow)
              b(pivrow)=b(row)
              b(row)=atemp
              do col=1,n
                 atemp=a(pivrow,col)
                 a(pivrow,col)=a(row,col)
                 a(row,col)=atemp
              enddo
           endif

           do rrow=row+1,n
              if (a(row,row).ne.0) then
   
                 m=-a(rrow,n-row+1)/a(row,n-row+1)
                 do col=1,n
                    a(rrow,col)=m*a(row,col)+a(rrow,col)
                 enddo
                 b(rrow)=m*b(row)+b(rrow)
              endif
           enddo
        enddo
!	back substitution now

        x(1)=b(n)/a(n,1)
        do row=n-1,1,-1
           x(n-row+1)=b(row)
           do col=1,(n-row)
              x(n-row+1)=(x(n-row+1)-a(row,col)*x(col)) 
           enddo

           x(n-row+1)=(x(n-row+1)/a(row,(n-row)+1))
        enddo

        return
   END SUBROUTINE gelim

   SUBROUTINE harm_ana_out
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE harm_ana_init  ***
      !!                   
      !! ** Purpose :   initialization of ....
      !!
      !! ** Method  :   blah blah blah ...
      !!
      !! ** input   :   Namlist namexa
      !!
      !! ** Action  :   ...  
      !!
      !! history :
      !!   9.0  !  03-08  (Autor Names)  Original code
      !!----------------------------------------------------------------------
        USE dianam          ! build name of file (routine)
 
      !! * local declarations
      INTEGER :: ji, jj, jk, jgrid,jit, ih, jjk   ! dummy loop indices
      INTEGER :: nh_T
      INTEGER :: nid_harm
      CHARACTER (len=40) :: clhstnamt, clop1, clop2 ! temporary names 
      CHARACTER (len=40) :: clhstnamu, clhstnamv    ! temporary names 
      CHARACTER (len=40) :: suffix
      REAL(wp) :: zsto1, zsto2, zout, zmax, zjulian, zdt, zmdi  ! temporary scalars

      !IF (lwp) THEN
      !      WRITE(numout,*) '--------------------------'
      !      WRITE(numout,*) 'ANA F NEW', anaf
      !      WRITE(numout,*) 'ANA U NEW', anau
      !      WRITE(numout,*) 'ANA V NEW', anav
      !      WRITE(numout,*) '--------------------------'
      !ENDIF


# if defined key_3Ddiaharm
      do jgrid=1,5
# else        
      do jgrid=1,3
# endif
          do ih=1,nb_ana
             hout = 0.0
             gout = 0.0
             do jj=1,nlcj
                do ji=1,nlci
                   cca=cosamp2D(ih,ji,jj,jgrid)
                   ssa=sinamp2D(ih,ji,jj,jgrid)
                   hout(ji,jj)=sqrt(cca**2+ssa**2)
                   IF (cca.eq.0.0 .and. ssa.eq.0.0) THEN 
                      gout(ji,jj)= 0.0_wp
                   ELSE
                      gout(ji,jj)=(180.0/rpi)*atan2(ssa,cca)       
                   ENDIF 
                   IF (hout(ji,jj).ne.0) THEN
                       hout(ji,jj)=hout(ji,jj)/anaf(ih)
                   ENDIF
                   IF (gout(ji,jj).ne.0) THEN  !Correct and take modulus
                       gout(ji,jj) = gout(ji,jj) + MOD( (anau(ih)+anav(ih))/rad , 360.0)
                       if (gout(ji,jj).gt.360.0) then
                           gout(ji,jj)=gout(ji,jj)-360.0
                       else if (gout(ji,jj).lt.0.0) then
                           gout(ji,jj)=gout(ji,jj)+360.0
                       endif
                   ENDIF
                enddo
             enddo
             !
             ! NETCDF OUTPUT
             if ( jgrid==1 ) suffix = ''
             if ( jgrid==2 ) suffix = '_u2D'
             if ( jgrid==3 ) suffix = '_v2D'
# if defined key_3Ddiaharm
             if ( jgrid==4 ) suffix = '_bfrU'
             if ( jgrid==5 ) suffix = '_bfrV'
# endif
             CALL iom_put( TRIM(Wave(ntide_all(ih))%cname_tide)//'amp'  //TRIM(suffix), hout(:,:) )
             CALL iom_put( TRIM(Wave(ntide_all(ih))%cname_tide)//'phase'//TRIM(suffix), gout(:,:) )
             if ( jgrid==1 ) then
!              CALL iom_put( TRIM(Wave(ntide_all(ih))%cname_tide)//'x_new'//TRIM(suffix), cosamp2D(ih,:,:,jgrid) )
!              CALL iom_put( TRIM(Wave(ntide_all(ih))%cname_tide)//'y_new'//TRIM(suffix), sinamp2D(ih,:,:,jgrid) )
              CALL iom_put( TRIM(Wave(ntide_all(ih))%cname_tide)//'x'//TRIM(suffix), cosamp2D(ih,:,:,jgrid) )
              CALL iom_put( TRIM(Wave(ntide_all(ih))%cname_tide)//'y'//TRIM(suffix), sinamp2D(ih,:,:,jgrid) )
             endif

          enddo
      enddo
# if defined key_3Ddiaharm
!
! DO THE SAME FOR 3D VARIABLES
!
      do jgrid=1,5
          do ih=1,nb_ana
             hout3D = 0.0
             gout3D = 0.0
             DO jk=1,jpkm1
                do jj=1,nlcj
                   do ji=1,nlci
                      cca=cosamp3D(ih,ji,jj,jk,jgrid)
                      ssa=sinamp3D(ih,ji,jj,jk,jgrid)
                      hout3D(ji,jj,jk)=sqrt(cca**2+ssa**2)
                      IF (cca.eq.0.0 .and. ssa.eq.0.0) THEN
                         gout3D(ji,jj,jk) = 0.0_wp
                      ELSE
                         gout3D(ji,jj,jk) = (180.0/rpi)*atan2(ssa,cca)
                      ENDIF
                      IF (hout3D(ji,jj,jk).ne.0) THEN
                          hout3D(ji,jj,jk) = hout3D(ji,jj,jk)
!!/anaf(ih)
                      ENDIF
                      IF (gout3D(ji,jj,jk).ne.0) THEN  !Correct and take modulus
!!                          gout3D(ji,jj,jk) = gout3D(ji,jj,jk) + MOD( (anau(ih)+anav(ih))/rad , 360.0)
                          if      (gout3D(ji,jj,jk).gt.360.0) then
                                   gout3D(ji,jj,jk) = gout3D(ji,jj,jk)-360.0
                          else if (gout3D(ji,jj,jk).lt.0.0) then
                                   gout3D(ji,jj,jk) = gout3D(ji,jj,jk)+360.0
                          endif
                      ENDIF
                   enddo    ! ji
                enddo       ! jj
             ENDDO          ! jk
             !
             ! NETCDF OUTPUT
             if ( jgrid==1 ) suffix = '_zdi'
             if ( jgrid==2 ) suffix = '_rhd'
             if ( jgrid==3 ) suffix = '_u3D'
             if ( jgrid==4 ) suffix = '_v3D'
             if ( jgrid==5 ) suffix = '_w3D'
             CALL iom_put( TRIM(Wave(ntide_all(ih))%cname_tide)//'amp'  //TRIM(suffix), hout3D(:,:,:) )
             CALL iom_put( TRIM(Wave(ntide_all(ih))%cname_tide)//'phase'//TRIM(suffix), gout3D(:,:,:) )
          enddo             ! ih 
      enddo                 ! jgrid
# endif
!
   END SUBROUTINE harm_ana_out
!
!   SUBROUTINE Compute_Amp_Phase
!      !!----------------------------------------------------------------------
!      !!                  ***  ROUTINE Compute_Amp_Phase  ***
!      !!   Compute Amplitude and Phase from cosinus and sin
!   END SUBROUTINE Compute_Amp_Phase



   SUBROUTINE harm_rst_write(kt)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE harm_ana_init  ***
      !!                   
      !! ** Purpose :  To write out cummulated Tidal Harmomnic data to file for
      !!               restarting
      !!
      !! ** Method  :   restart files will be dated by default
      !!
      !! ** input   :   
      !!
      !! ** Action  :   ...  
      !!
      !! history :
      !!   0.0  !  01-16  (Enda O'Dea)  Original code
      !! ASSUMES  dated file for rose  , can change later to be more generic
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ocean time-step
      !!
      INTEGER             ::   ih
      CHARACTER(LEN=20)   ::   clkt     ! ocean time-step define as a character
      CHARACTER(LEN=50)   ::   clname   ! ocean output restart file name
      CHARACTER(LEN=150)  ::   clpath   ! full path to ocean output restart file
      CHARACTER(LEN=250)  ::   clfinal   ! full name

      !These Three which contain the most data can be moved to a regular
      !restart file
      CALL iom_rstput( kt, nitrst, numrow, 'Mean_bssh'     , bssh (1,:,:)       )
      CALL iom_rstput( kt, nitrst, numrow, 'Mean_bubar'    , bubar(1,:,:)       )
      CALL iom_rstput( kt, nitrst, numrow, 'Mean_bvbar'    , bvbar(1,:,:)       )
# if defined key_3Ddiaharm
      CALL iom_rstput( kt, nitrst, numrow, 'Mean_ana_bfrU' , ana_bfrU(1,:,:)    )
      CALL iom_rstput( kt, nitrst, numrow, 'Mean_ana_bfrV' , ana_bfrV(1,:,:)    )
      CALL iom_rstput( kt, nitrst, numrow, 'Mean_ana_D'    , ana_D   (1,:,:,:)  )
      CALL iom_rstput( kt, nitrst, numrow, 'Mean_ana_I'    , ana_I   (1,:,:,:)  )
      CALL iom_rstput( kt, nitrst, numrow, 'Mean_ana_U'    , ana_U   (1,:,:,:)  )
      CALL iom_rstput( kt, nitrst, numrow, 'Mean_ana_V'    , ana_V   (1,:,:,:)  )
      CALL iom_rstput( kt, nitrst, numrow, 'Mean_ana_W'    , ana_W   (1,:,:,:)  )
#endif

      do ih=1,nb_ana
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'bssh_cos'    , bssh (ih*2  , : , : )    )
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'bssh_sin'    , bssh (ih*2+1, : , : )    )
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'bubar_cos'   , bubar(ih*2  , : , : )    )
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'bubar_sin'   , bubar(ih*2+1, : , : )    )
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'bvbar_cos'   , bvbar(ih*2  , : , : )    )
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'bvbar_sin'   , bvbar(ih*2+1, : , : )    )
# if defined key_3Ddiaharm
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_bfrU_cos', ana_bfrU(ih*2  , : , : ) )
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_bfrU_sin', ana_bfrU(ih*2+1, : , : ) )
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_bfrV_cos', ana_bfrV(ih*2  , : , : ) )
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_bfrV_sin', ana_bfrV(ih*2+1, : , : ) )

         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_D_cos'   , ana_D(ih*2  , : , :, : ) )
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_D_sin'   , ana_D(ih*2+1, : , :, : ) )
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_I_sin'   , ana_I(ih*2+1, : , :, : ) )
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_I_sin'   , ana_I(ih*2+1, : , :, : ) )
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_U_sin'   , ana_U(ih*2+1, : , :, : ) )
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_U_sin'   , ana_U(ih*2+1, : , :, : ) )
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_V_sin'   , ana_V(ih*2+1, : , :, : ) )
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_V_sin'   , ana_V(ih*2+1, : , :, : ) )
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_W_sin'   , ana_W(ih*2+1, : , :, : ) )
         CALL iom_rstput( kt, nitrst, numrow, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_W_sin'   , ana_W(ih*2+1, : , :, : ) )
# endif
      enddo

      IF(lwp) THEN
        IF( kt > 999999999 ) THEN ; WRITE(clkt, *       ) kt
        ELSE                      ; WRITE(clkt, '(i8.8)') kt
        ENDIF
        clname = TRIM(cexper)//"_"//TRIM(ADJUSTL(clkt))//"_restart_harm_ana.bin"
        clpath = TRIM(cn_ocerst_outdir)
        IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath) // '/'
        IF (lwp) WRITE(numout,*) 'Open tidal harmonics restart file for writing: ',TRIM(clpath)//clname

        WRITE(clfinal,'(a)') trim(clpath)//trim(clname)
        OPEN( 66, file=TRIM(clfinal), form='unformatted', access="stream" )
        WRITE(66) cc
        WRITE(66) anau
        WRITE(66) anav
        WRITE(66) anaf
        WRITE(66) fjulday_startharm
        CLOSE(66)
        WRITE(numout,*) '----------------------------'
        WRITE(numout,*) '   harm_rst_write: DONE '
        WRITE(numout,*) cc
        WRITE(numout,*) anaf
        WRITE(numout,*) fjulday_startharm
        WRITE(numout,*) MAXVAL(bssh), MAXVAL(bubar), MAXVAL(bvbar)
        WRITE(numout,*) '----------------------------'
      ENDIF
 
   END SUBROUTINE harm_rst_write

   SUBROUTINE harm_rst_read
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE harm_ana_init  ***
      !!                   
      !! ** Purpose :  To read in  cummulated Tidal Harmomnic data to file for
      !!               restarting
      !!
      !! ** Method  :   
      !!
      !! ** input   :   
      !!
      !! ** Action  :   ...  
      !!
      !! history :
      !!   0.0  !  01-16  (Enda O'Dea)  Original code
      !! ASSUMES  dated file for rose  , can change later to be more generic
      !!----------------------------------------------------------------------
      CHARACTER(LEN=20)   ::   clkt     ! ocean time-step define as a character
      CHARACTER(LEN=50)   ::   clname   ! ocean output restart file name
      CHARACTER(LEN=150)  ::   clpath   ! full path to ocean output restart file
      CHARACTER(LEN=250)  ::   clfinal   ! full name
      INTEGER             ::   ih

      IF( nit000 > 999999999 ) THEN ; WRITE(clkt, *       ) nit000-1
      ELSE                      ; WRITE(clkt, '(i8.8)') nit000-1
      ENDIF
      clname = TRIM(cexper)//"_"//TRIM(ADJUSTL(clkt))//"_restart_harm_ana.bin"
      clpath = TRIM(cn_ocerst_outdir)
      IF( clpath(LEN_TRIM(clpath):) /= '/' ) clpath = TRIM(clpath) // '/'

      IF (lwp) WRITE(numout,*) 'Open tidal harmonics restart file for reading: ',TRIM(clpath)//clname

      !2D regular fields can be read from normal restart this saves space and handy to
      !view in netcdf format also.
      CALL iom_get( numror,jpdom_autoglo, 'Mean_bssh'     , bssh (1,:,:)       )
      CALL iom_get( numror,jpdom_autoglo, 'Mean_bubar'    , bubar(1,:,:)       )
      CALL iom_get( numror,jpdom_autoglo, 'Mean_bvbar'    , bvbar(1,:,:)       )
# if defined key_3Ddiaharm
      CALL iom_get( numror,jpdom_autoglo, 'Mean_ana_bfrU' , ana_bfrU(1,:,:)    )
      CALL iom_get( numror,jpdom_autoglo, 'Mean_ana_bfrV' , ana_bfrV(1,:,:)    )
      CALL iom_get( numror,jpdom_autoglo, 'Mean_ana_D'    , ana_D   (1,:,:,:)  )
      CALL iom_get( numror,jpdom_autoglo, 'Mean_ana_I'    , ana_I   (1,:,:,:)  )
      CALL iom_get( numror,jpdom_autoglo, 'Mean_ana_U'    , ana_U   (1,:,:,:)  )
      CALL iom_get( numror,jpdom_autoglo, 'Mean_ana_V'    , ana_V   (1,:,:,:)  )
      CALL iom_get( numror,jpdom_autoglo, 'Mean_ana_W'    , ana_W   (1,:,:,:)  )
#endif          

      do ih=1,nb_ana
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'bssh_cos'    , bssh (ih*2  , : , : )    )
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'bssh_sin'    , bssh (ih*2+1, : , : )    )
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'bubar_cos'   , bubar(ih*2  , : , : )    )
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'bubar_sin'   , bubar(ih*2+1, : , : )    )
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'bvbar_cos'   , bvbar(ih*2  , : , : )    )
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'bvbar_sin'   , bvbar(ih*2+1, : , : )    )
# if defined key_3Ddiaharm
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_bfrU_cos', ana_bfrU(ih*2  , : , : ) )
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_bfrU_sin', ana_bfrU(ih*2+1, : , : ) )
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_bfrV_cos', ana_bfrV(ih*2  , : , : ) )
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_bfrV_sin', ana_bfrV(ih*2+1, : , : ) )

         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_D_cos'   , ana_D(ih*2  , : , :, : ) )
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_D_sin'   , ana_D(ih*2+1, : , :, : ) )
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_I_sin'   , ana_I(ih*2+1, : , :, : ) )
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_I_sin'   , ana_I(ih*2+1, : , :, : ) )
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_U_sin'   , ana_U(ih*2+1, : , :, : ) )
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_U_sin'   , ana_U(ih*2+1, : , :, : ) )
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_V_sin'   , ana_V(ih*2+1, : , :, : ) )
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_V_sin'   , ana_V(ih*2+1, : , :, : ) )
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_W_sin'   , ana_W(ih*2+1, : , :, : ) )
         CALL iom_get( numror,jpdom_autoglo, TRIM(Wave(ntide_all(ih))%cname_tide)//'ana_W_sin'   , ana_W(ih*2+1, : , :, : ) )
# endif
      enddo

      WRITE(clfinal,'(a)') trim(clpath)//trim(clname)
      OPEN( 66, file=TRIM(clfinal), form='unformatted', access="stream" )
      READ(66) cc
      READ(66) anau
      READ(66) anav
      READ(66) anaf
      READ(66) fjulday_startharm
      CLOSE(66)

      IF(lwp) THEN
        WRITE(numout,*) '----------------------------'
        WRITE(numout,*) '   Checking anaf is correct'
        WRITE(numout,*) cc
        WRITE(numout,*) anaf
        WRITE(numout,*) fjulday_startharm
        WRITE(numout,*) MAXVAL(bssh), MAXVAL(bubar), MAXVAL(bvbar)
        WRITE(numout,*) '----------------------------'
      ENDIF
 
   END SUBROUTINE harm_rst_read

   !!======================================================================
#else
!!---------------------------------------------------------------------------------
!!   Dummy module                                   NO harmonic Analysis
!!---------------------------------------------------------------------------------
        CONTAINS
           SUBROUTINE harm_rst_write(kt)     ! Dummy routine
           END SUBROUTINE harm_rst_write
           SUBROUTINE harm_rst_read    ! Dummy routine
           END SUBROUTINE harm_rst_read
           SUBROUTINE harm_ana_out      ! Dummy routine
           END SUBROUTINE harm_ana_out
           SUBROUTINE harm_ana_init
           END SUBROUTINE harm_ana_init
           SUBROUTINE harm_ana( kt )
!--- NB : end call not properly written
           END SUBROUTINE harm_ana
!           END SUBROUTINE harm_ana_init
!--- END NB
           SUBROUTINE gelim (a,b,x,n)
!--- NB : end call not properly written
           END SUBROUTINE gelim
!           END SUBROUTINE gelim (a,b,x,n)
!--- END NB           
#endif

END MODULE harmonic_analysis
