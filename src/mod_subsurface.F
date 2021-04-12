!------------------------------------
! Subsurface module
!------------------------------------

# if defined (SUBSURFACE)

MODULE SUBSURFACE_MODULE
  USE PARAM
  USE GLOBAL,ONLY : Mloc,Nloc,Mloc1,Nloc1,Nghost,  &
                    Mglob,Nglob, &
                    Ibeg,Iend,Jbeg,Jend,Iend1,Jend1,  &
                    DX,DY, DX2,DY2,DXDY,TIME, &
                    H,H0,H_u,H_v,ETA,U,V,P,Q,PQ,  &
                    MASK,MASKu,MASKv,MaxDepth, &
                    MinDepth,MASK,DT,tmp4preview, &
                    Numpixel, SubMainGridRatio, DepSubGrid

  USE INPUT_Util
#if defined (PARALLEL)
  USE GLOBAL,ONLY : myid,ier, npx,npy,PX,PY,n_west,n_east,n_suth,n_nrth
  USE MPI
# endif
  IMPLICIT NONE
  SAVE

  LOGICAL :: SATURATION =.TRUE.
  LOGICAL :: OUT_ETA_GROUND = .FALSE.
  LOGICAL :: OUT_P_GROUND = .FALSE.
  LOGICAL :: OUT_Q_GROUND = .FALSE.
  LOGICAL :: OUT_MASK_GROUND = .FALSE.
  LOGICAL :: OUT_MASKu_GROUND = .FALSE.
  LOGICAL :: OUT_MASKv_GROUND = .FALSE.

  REAL(SP):: PLOT_INTV_ground,PLOT_COUNT_ground
  REAL(SP),DIMENSION(:,:),ALLOCATABLE::PORO_ground,ETA_ground,P_ground,Q_ground, &
                         Du_ground,Su_ground,Dv_ground,Sv_ground,  &
                         KAu_ground,KAv_ground,DepthX,DepthY, &
                         ETAu_ground,ETAv_ground,POROgrn
  INTEGER, DIMENSION(:,:),ALLOCATABLE::MASK_ground,MASKu_ground,MASKv_ground
  REAL(SP) :: S_constant,KA_constant,POR_constant

! statistics

    REAL(SP), DIMENSION(:,:),ALLOCATABLE :: MaxEtaGrn,MinEtaGrn
    REAL(SP) :: StartTimeStatistics

!-------------FOR UNSATURATION-----------------------------------
    REAL(SP) :: Gamma, Parameter_n
    REAL(SP),DIMENSION(:,:),ALLOCATABLE:: Omega, Omega0,Depth,tmp_mos
    REAL(SP) :: gamma_rcp,tmp_coeff
    CHARACTER(LEN=80):: MODEL_TYPE
# if defined (MULTIRELATIONSHIP)
    REAL(SP) :: ddep
# endif
    LOGICAL :: OUT_Omega = .FALSE.

!-------------FOR VEGEROOT--------------------------------------
    LOGICAL :: VEGEROOT_CLAMPED =.FALSE.
    LOGICAL :: HOT_IntialBIOMASS = .FALSE.
    CHARACTER(LEN=80):: IntialBIO1_FILE,IntialBIO2_FILE,&
                        Eta0_FILE
                        

    INTEGER :: num_species

    REAL(SP):: MRD1,MRD2,E11,E12,E13,E14, &
               E21,E22,E23,E24, &
               ET_t_max1,ET_t_max2,ET_e, &
               growth1,growth2,Tidal_min,  &
               Bio_init1, Bio_init2, &
# if defined (MULTIRELATIONSHIP)
               E_dot, &
# endif
               n_vege, &
               Beta1, Beta2, Beta12
   REAL(SP),DIMENSION(:,:),ALLOCATABLE::avgEr1,avgEr2,ET_t1,ET_t2,&
                                   Bio1, Bio2, DeltaBio1, DeltaBio2,&
                                   Bio10, Bio20,&
                                   Coeff_uptake1, Coeff_capct1,&
                                   Coeff_uptake2, Coeff_capct2

    LOGICAL :: OUT_Biomass = .FALSE.
    LOGICAL :: OUT_avgERoot = .FALSE.
    LOGICAL :: OUT_Coefficient = .FALSE.
    LOGICAL :: OUT_Transpiration = .FALSE.
!---------------------------------------------------------------

#if defined (PARALLEL)
    REAL(SP) :: myvar
# endif 

CONTAINS

! ------------------------------  
SUBROUTINE SUBSURFACE_INITIAL
!  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5
                    
  USE Input_Util
  IMPLICIT NONE

  CHARACTER(LEN=80) FILE_NAME
  INTEGER::LINE
  INTEGER :: ierr

! read parameters from input.txt
      FILE_NAME='input.txt'

      CALL GET_LOGICAL_VAL(SATURATION,FILE_NAME,'SATURATION',line,ierr)

      IF(ierr==1)THEN
        SATURATION = .TRUE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'No SATURATION defined, USE defalt: TRUE'
         WRITE(3,'(A80)')'No SATURATION defined, USE defalt: TRUE'
      endif
# else
         WRITE(*,'(A80)')'No SATURATION defined, USE defalt: TRUE'
         WRITE(3,'(A80)')'No SATURATION defined, USE defalt: TRUE'
# endif
       ENDIF

      CALL GET_LOGICAL_VAL(VEGEROOT_CLAMPED,FILE_NAME,'VEGEROOT_CLAMPED',line,ierr)

      IF(ierr==1)THEN
        VEGEROOT_CLAMPED = .FALSE.
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'No VEGEROOT_CLAMPED defined, USE defalt: FALSE'
         WRITE(3,'(A80)')'No VEGEROOT_CLAMPED defined, USE defalt: FALSE'
      endif
# else
         WRITE(*,'(A80)')'No VEGEROOT_CLAMPED defined, USE defalt: FALSE'
         WRITE(3,'(A80)')'No VEGEROOT_CLAMPED defined, USE defalt: FALSE'
# endif
       ENDIF
      

      CALL GET_LOGICAL_VAL(OUT_P_GROUND,FILE_NAME,'P_GROUND',line,ierr)
      CALL GET_LOGICAL_VAL(OUT_Q_GROUND,FILE_NAME,'Q_GROUND',line,ierr)
      CALL GET_LOGICAL_VAL(OUT_ETA_GROUND,FILE_NAME,'ETA_GROUND',line,ierr)
      CALL GET_LOGICAL_VAL(OUT_MASK_GROUND,FILE_NAME,'MASK_GROUND',line,ierr)
      CALL GET_LOGICAL_VAL(OUT_MASKu_GROUND,FILE_NAME,'MASKu_GROUND',line,ierr)
      CALL GET_LOGICAL_VAL(OUT_MASKv_GROUND,FILE_NAME,'MASKv_GROUND',line,ierr)

      CALL GET_Float_VAL(KA_constant,FILE_NAME,'Kappa',line,ierr)
      IF(ierr==1)THEN
        KA_constant = 0.01
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'Kappa: Use default,', '0.01'
         WRITE(3,'(A80)')'Kappa: Use default,', '0.01'
      endif
# else
         WRITE(*,'(A80)')'Kappa: Use default,', '0.01'
         WRITE(3,'(A80)')'Kappa: Use default,', '0.01'
# endif
      ENDIF

! constant porosity
      CALL GET_Float_VAL(POR_constant,FILE_NAME,'POROSITY_GROUND',line,ierr)
      IF(ierr==1)THEN
        POR_constant = 0.3
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'POROgrn: Use default,', '0.3'
         WRITE(3,'(A80)')'POROgrn: Use default,', '0.3'
      endif
# else
         WRITE(*,'(A80)')'POROgrn: Use default,', '0.3'
         WRITE(3,'(A80)')'POROgrn: Use default,', '0.3'
# endif
      ENDIF

      CALL GET_Float_VAL(S_constant,FILE_NAME,'HardBottom',line,ierr)
      IF(ierr==1)THEN
        S_constant = MaxDepth
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'S_constant: Use default,', 'MaxDepth'
         WRITE(3,'(A80)')'S_constant: Use default,', 'MaxDepth'
      endif
# else
         WRITE(*,'(A80)')'S_constant: Use default,', 'MaxDepth'
         WRITE(3,'(A80)')'S_constant: Use default,', 'MaxDepth'
# endif
      ENDIF

      CALL GET_Float_VAL(StartTimeStatistics,FILE_NAME,'StartTimeStatistics',line,ierr)
      IF(ierr==1)THEN
          StartTimeStatistics = 0.0
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'StartTimeStatistics: Use default,', 'ZERO'
         WRITE(3,'(A40,A40)')'StartTimeStatistics: Use default,', 'ZERO'
      endif
# else
         WRITE(*,'(A40,A40)')'StartTimeStatistics: Use default,', 'ZERO'
         WRITE(3,'(A40,A40)')'StartTimeStatistics: Use default,', 'ZERO'
# endif
      ENDIF


!write in LOG.txt

# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A50)') '-----FOR SUBSURFAVE MODULE------'
       if(myid.eq.0) WRITE(3,'(A19,L12)') 'SATURATION =', SATURATION
       if(myid.eq.0) WRITE(3,'(A19,F12.2)') 'Hardbottom =', S_constant
       if(myid.eq.0) WRITE(3,'(A19,F12.2)') 'POROSITY_GROUND =', POR_constant
       if(myid.eq.0) WRITE(3,'(A19,F12.6)') 'Kappa =', KA_constant
# else
       WRITE(3,'(A50)') '-----FOR SUBSURFAVE MODULE------'
       WRITE(3,'(A19,L12)') 'SATURATION =', SATURATION
       WRITE(3,'(A19,F12.2)') 'Hardbottom =', S_constant
       WRITE(3,'(A19,F12.2)') 'POROSITY_GROUND =', POR_constant
       WRITE(3,'(A19,F12.2)') 'Kappa =', KA_constant

# endif

! ALLOCATION

        ALLOCATE(ETA_ground(Mloc,Nloc),P_ground(Mloc1,Nloc),Q_ground(Mloc,Nloc1), &
                 Du_ground(Mloc1,Nloc),Dv_ground(Mloc,Nloc1), &
                 Su_ground(Mloc1,Nloc),Sv_ground(Mloc,Nloc1), &
                 KAu_ground(Mloc1,Nloc),KAv_ground(Mloc,Nloc1), &
                 MASK_ground(Mloc,Nloc),MASKu_ground(Mloc1,Nloc), &
 !                POROgrn(Mloc,Nloc), &
                 MASKv_ground(Mloc,Nloc1), &
                 DepthX(Mloc1,Nloc),DepthY(Mloc,Nloc1), &
                 ETAu_ground(Mloc1,Nloc),ETAv_ground(Mloc,Nloc1))

    ALLOCATE(MaxEtaGrn(Mloc,Nloc),MinEtaGrn(Mloc,Nloc))

    MaxEtaGrn = -1.0*LARGE
    MinEtaGrn = LARGE

! assign variables

    MASKu_ground = 0
    MASKv_ground = 0
    MASK_ground = 0
    ETA_ground = ZERO
    ETAu_ground = ZERO
    ETAv_ground = ZERO
    P_ground = ZERO
    Q_ground = ZERO
    Du_ground = ZERO
    Dv_ground = ZERO

    Su_ground = S_constant
    Sv_ground = S_constant
    KAu_ground = KA_constant 
    KAv_ground = KA_constant
!    POROgrn = POR_constant

    CALL UPDATE_GROUND_MASK_ETA
  
    IF (.NOT.SATURATION)THEN
     CALL UNSATURATION_INITIAL
     CALL CALCULATE_MOISTURE
        IF(VEGEROOT_CLAMPED)THEN
           CALL VEGEROOT_INITIAL
        ENDIF
    ENDIF

END SUBROUTINE SUBSURFACE_INITIAL

!------------------------

SUBROUTINE UPDATE_GROUND_MASK_ETA


! depthx and depthy
! update with time
     DepthX = H_u
     DepthY = H_v

! eta at u v points

   DO J=1,Nloc
   DO I=Ibeg,Iend1
       ETAu_ground(I,J)=0.5_SP*(ETA_ground(I-1,J)+ETA_ground(I,J))
   ENDDO
   ENDDO

   DO J=Jbeg,Jend1
   DO I=1,Mloc
       ETAv_ground(I,J)=0.5_SP*(ETA_ground(I,J-1)+ETA_ground(I,J))
   ENDDO
   ENDDO

# if defined (PARALLEL)
     CALL PHI_EXCH_2(ETAu_ground)
     CALL PHI_EXCH_3(ETAv_ground)
# endif

! thickness

    DO J=1,Nloc
    DO I=1,Mloc1
      IF(MASKu(I,J) == 0)THEN
        Du_ground(I,J) = Su_ground(I,J) +ETAu_ground(I,J)
      ELSE
        Du_ground(I,J) = Su_ground(I,J)-DepthX(I,J)
      ENDIF
    ENDDO
    ENDDO

    DO J=1,Nloc1
    DO I=1,Mloc
      IF(MASKv(I,J) == 0)THEN
        Dv_ground(I,J) = Sv_ground(I,J) +ETAv_ground(I,J)
      ELSE
        Dv_ground(I,J) = Sv_ground(I,J)-DepthY(I,J)
      ENDIF
    ENDDO
    ENDDO


END SUBROUTINE UPDATE_GROUND_MASK_ETA

SUBROUTINE UPDATE_SUBSURFACE
!  USE GLOBAL,ONLY : itmp1,itmp2,itmp3,itmp4,itmp5,SMALL

  IMPLICIT NONE
!  INTEGER::ISTEP,ISTAGE,IVAR

     CALL UPDATE_GROUND_MASK_ETA

     DO J=Jbeg,Jend
     DO I=Ibeg,Iend1
       P_ground(I,J) = - KAu_ground(I,J)*Du_ground(I,J) &
                    *(ETA_ground(I,J)-ETA_ground(I-1,J))/DX
     ENDDO
     ENDDO

     DO J=Jbeg,Jend1
     DO I=Ibeg,Iend
       Q_ground(I,J) = - KAv_ground(I,J)*Dv_ground(I,J) &
                    *(ETA_ground(I,J)-ETA_ground(I,J-1))/DY
     ENDDO
     ENDDO

# if defined (PARALLEL)
     CALL PHI_EXCH_2(P_ground)
     CALL PHI_EXCH_3(Q_ground)
# endif

     IF (SATURATION)THEN

        DO J=Jbeg,Jend
        DO I=Ibeg,Iend
           IF(MASK(I,J)==0)THEN
              ETA_ground(I,J) = ETA_ground(I,J)  &
                  -1.0_SP/POR_constant*((P_ground(I+1,J)-P_ground(I,J))*DT/DX  &
                  +(Q_ground(I,J+1)-Q_ground(I,J))*DT/DY)
           ELSE
              ETA_ground(I,J) = ETA(I,J)
           ENDIF
        ENDDO
        ENDDO
   
 
      ELSE
 
          CALL CALCULATE_MOISTURE
          IF (VEGEROOT_CLAMPED) THEN
             DO J=1,Nloc
             DO I=1,Mloc
                IF(MASK(I,J)==0)THEN
                   CALL CALCULATE_avgEr(MRD1,avgEr1(I,J),ETA_GROUND(I,J),Depth(I,J))
                ELSE
                   avgEr1(I,J)=1.0
                ENDIF
             ENDDO
             ENDDO
               CALL GET_ROOTUPTAKE_COEFF(E11,E12,E13,E14,avgEr1,Coeff_uptake1)
               CALL GET_CARRYCAPACITY_COEFF(E11,E14,avgEr1,Coeff_capct1)

              IF(num_species == 2) THEN
                 DO J=1,Nloc
                 DO I=1,Mloc
                    IF(MASK(I,J)==0)THEN
                       CALL CALCULATE_avgEr(MRD2,avgEr2(I,J),ETA_GROUND(I,J),Depth(I,J))
                    ELSE
                       avgEr2(I,J)=1.0
                    ENDIF
                 ENDDO
                 ENDDO

                 CALL GET_ROOTUPTAKE_COEFF(E21,E22,E23,E24,avgEr2,Coeff_uptake2)
                 CALL GET_CARRYCAPACITY_COEFF(E21,E24,avgEr2,Coeff_capct2)
              ENDIF
             CALL UPDATE_BIOMASS

          ENDIF
          

        DO J=Jbeg,Jend
        DO I=Ibeg,Iend
           IF(MASK(I,J)==0)THEN
              IF (VEGEROOT_CLAMPED) THEN
                IF(num_species==1)THEN
                   ETA_ground(I,J) = ETA_ground(I,J)  &
                     -1.0_SP/POR_constant*((P_ground(I+1,J)-P_ground(I,J))*DT/DX  &
                     +(Q_ground(I,J+1)-Q_ground(I,J))*DT/DY)+Omega0(I,J)-Omega(I,J) &
                     -DT*(ET_e+Bio1(I,J)*ET_t_max1*Coeff_uptake1(I,J))
                 ELSEIF(num_species==2)THEN
                   ETA_ground(I,J) = ETA_ground(I,J)  &
                     -1.0_SP/POR_constant*((P_ground(I+1,J)-P_ground(I,J))*DT/DX  &
                     +(Q_ground(I,J+1)-Q_ground(I,J))*DT/DY)+Omega0(I,J)-Omega(I,J) &
                     -DT*(ET_e+Bio1(I,J)*ET_t_max1*Coeff_uptake1(I,J) &
                     +Bio2(I,J)*ET_t_max2*Coeff_uptake2(I,J))
                ENDIF
              ELSE
               ETA_ground(I,J) = ETA_ground(I,J)  &
                  -1.0_SP/POR_constant*((P_ground(I+1,J)-P_ground(I,J))*DT/DX  &
                  +(Q_ground(I,J+1)-Q_ground(I,J))*DT/DY)+Omega0(I,J)-Omega(I,J)
              ENDIF
           ELSE
              ETA_ground(I,J) = ETA(I,J)
           ENDIF
       ENDDO
       ENDDO

     ENDIF

# if defined (PARALLEL)
     CALL PHI_EXCH_1(ETA_ground)
# endif

    CALL FILLIN_GRN_GHOST

     IF(TIME>StartTimeStatistics)THEN
       DO J=1,Nloc
       DO I=1,Mloc
         IF (MASK(I,J) < 1.0_SP)THEN
           IF(ETA_ground(I,J)>MaxEtaGrn(I,J))MaxEtaGrn(I,J)=ETA_ground(I,J)
           IF(ETA_ground(I,J)<MinEtaGrn(I,J))MinEtaGrn(I,J)=ETA_ground(I,J)
         ENDIF       
       ENDDO
       ENDDO
     ENDIF


END SUBROUTINE UPDATE_SUBSURFACE


SUBROUTINE FILLIN_GRN_GHOST

!	west
#	if defined(PARALLEL)
        if(n_west .eq.MPI_PROC_NULL)then
#	endif

        DO J=1,Nloc
        DO I=1,Nghost
            ETA_ground(I,J) = ETA_ground(Ibeg,J)
            P_ground(I,J)=P_ground(Ibeg,J)
            Du_ground(I,J)=Du_ground(Ibeg,J)
        ENDDO
        ENDDO

        DO J=Jbeg,Jend1
        DO I=1,Nghost
            Q_ground(I,J)=Q_ground(Ibeg,J)
            Dv_ground(I,J)=Dv_ground(Ibeg,J)
        ENDDO
        ENDDO

#	if defined(PARALLEL)
       endif
#	endif

!	east
#	if defined(PARALLEL)
        if(n_east .eq.MPI_PROC_NULL)then
#	endif

        DO J=1,Nloc
        DO I=Iend1+1,Mloc1
            ETA_ground(I-1,J) = ETA_ground(Iend,J)
            P_ground(I,J)=P_ground(Iend1,J)
            Du_ground(I,J)=Du_ground(Iend1,J)
        ENDDO
        ENDDO

        DO J=Jbeg,Jend1
        DO I=Iend1,Mloc
            Q_ground(I,J)=Q_ground(Iend,J)
            Dv_ground(I,J)=Dv_ground(Iend,J)
        ENDDO
        ENDDO

#	if defined(PARALLEL)
        endif
#	endif

!	south
#	if defined(PARALLEL)
       if(n_suth .eq.MPI_PROC_NULL)then
#	endif

        DO J=1,Nghost
        DO I=1,Mloc
            ETA_ground(I,J) = ETA_ground(I,Jbeg)
            Q_ground(I,J)=Q_ground(I,Jbeg)
            Dv_ground(I,J)=Dv_ground(I,Jbeg)
        ENDDO
        ENDDO

        DO J=1,Nghost
        DO I=Ibeg,Iend1
            P_ground(I,J)=P_ground(I,Jbeg)
            Du_ground(I,J)=Du_ground(I,Jbeg)
        ENDDO
        ENDDO

#	if defined(PARALLEL)
       endif
#	endif

!	north
#	if defined(PARALLEL)
        if(n_nrth .eq.MPI_PROC_NULL)then
#	endif

        DO J=Jend1+1,Nloc1
        DO I=1,Mloc
           ETA_ground(I,J-1) = ETA_ground(I,Jend)
           Q_ground(I,J)=Q_ground(I,Jend1)
           Dv_ground(I,J)=Dv_ground(I,Jend1)
        ENDDO
        ENDDO

        DO J=Jend1,Nloc
        DO I=Ibeg,Iend1
           P_ground(I,J)=P_ground(I,Jend)
           Du_ground(I,J)=Du_ground(I,Jend)
        ENDDO
        ENDDO

#	if defined(PARALLEL)
       endif
#	endif

END SUBROUTINE FILLIN_GRN_GHOST

SUBROUTINE UNSATURATION_INITIAL
     USE Input_Util
     IMPLICIT NONE
     CHARACTER(LEN=80) FILE_NAME
     INTEGER::LINE
     INTEGER :: ierr
    
     REAL(SP) :: tmp_dep
 
     FILE_NAME='input.txt'

     CALL GET_LOGICAL_VAL(OUT_Omega,FILE_NAME,'MeanOmega',line,ierr)

     CALL GET_STRING_VAL(MODEL_TYPE,FILE_NAME,'MODEL_TYPE',line,ierr)
!    MODEL_TYPE: DARCY/B-C/vGen
           IF(ierr==1)THEN
             MODEL_TYPE = 'Brooks-Corey'
# if defined (PARALLEL)
           if (myid.eq.0) THEN
              WRITE(*,'(A40,A40)')'MODEL_TYPE: Use default,', 'B-C'
              WRITE(3,'(A40,A40)')'MODEL_TYPE: Use default,', 'B-C'
           endif
# else
              WRITE(*,'(A40,A40)')'MODEL_TYPE: Use default,', 'B-C'
              WRITE(3,'(A40,A40)')'MODEL_TYPE: Use default,', 'B-C'
# endif
           ENDIF

     ALLOCATE(Omega(Mloc,Nloc),Omega0(Mloc,Nloc),Depth(Mloc,Nloc),tmp_mos(I,J))

      CALL GET_Float_VAL(Gamma,FILE_NAME,'Gamma',line,ierr)
      IF(ierr==1)THEN
        Gamma = 29.0
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'Gamma: Use default,', '29.0'
         WRITE(3,'(A80)')'Gamma: Use default,', '29.0'
      endif
# else
         WRITE(*,'(A80)')'Gamma: Use default,', '29.0'
         WRITE(3,'(A80)')'Gamma: Use default,', '29.0'
# endif
      ENDIF

      CALL GET_Float_VAL(Parameter_n,FILE_NAME,'Parameter_n',line,ierr)
      IF(ierr==1)THEN
        Gamma = 4.0
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'Parameter_n: Use default,', '4.0'
         WRITE(3,'(A80)')'Parameter_n: Use default,', '4.0'
      endif
# else
         WRITE(*,'(A80)')'Parameter_n: Use default,', '4.0'
         WRITE(3,'(A80)')'Parameter_n: Use default,', '4.0'
# endif
      ENDIF

!write in LOG.txt
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A50)') '-----for UNSUBSURFACE-----'
       if(myid.eq.0) WRITE(3,'(A18,F12.2)') 'Gamma =', Gamma
       if(myid.eq.0) WRITE(3,'(A18,F12.2)') 'Parameter_n =', Parameter_n
# else
       WRITE(3,'(A50)') '-----FOR SUBSURFACE------'
       WRITE(3,'(A18,F12.2)') 'Gamma =', Gamma
       WRITE(3,'(A18,F12.2)') 'Parameter_n =', Parameter_n
# endif

     Omega0 = ZERO
     Omega = ZERO
     Depth = ZERO
     tmp_mos = ZERO

     Gamma_rcp = 1 / Gamma
     tmp_coeff = Gamma_rcp/(1-Parameter_n)

# if defined (MULTIRELATIONSHIP)
     ddep = MinDepth
# endif


! get depth in coarse grid
     
      DO J=1,Nloc
      DO I=1,Mloc
         tmp_dep = ZERO
         DO JJ=1,SubMainGridRatio
         DO II=1,SubMainGridRatio
            tmp_dep=tmp_dep+DepSubGrid(I,J,II,JJ)
         ENDDO
         ENDDO
            Depth(I,J)=tmp_dep/Numpixel
      ENDDO
      ENDDO
END SUBROUTINE UNSATURATION_INITIAL


# if defined (MULTIRELATIONSHIP)

SUBROUTINE CALCULATE_MOISTURE
     REAL(SP)::SUM
     Omega0=Omega
     DO J=1,Nloc
     DO I=1,Mloc
      IF (MASK(I,J)==0)THEN
         IF(MODEL_TYPE=='B-C')THEN
             IF((-DEPTH(I,J)).GE.(ETA_GROUND(I,J)+Gamma_rcp))THEN
                 CALL INTEGRATION(ETA_GROUND(I,J),ETA_GROUND(I,J)+Gamma_rcp,(-DEPTH(I,J)),ddep,SUM)
                 Omega(I,J)=(SUM+Gamma_rcp)*POR_constant
             ELSE
                 Omega(I,J)=(-DEPTH(I,J)-ETA_GROUND(I,J))*POR_constant
             ENDIF
         ENDIF
      ELSE
         Omega(I,J)=ZERO
      ENDIF
    ENDDO
    ENDDO
END SUBROUTINE CALCULATE_MOISTURE

# else

SUBROUTINE CALCULATE_MOISTURE
     Omega0=Omega
     DO J=1,Nloc
     DO I=1,Mloc
      IF (MASK(I,J)==0)THEN
             IF((-DEPTH(I,J)).GE.(ETA_GROUND(I,J)+Gamma_rcp))THEN
                 tmp_mos(I,J)=(-Gamma*min(-MinDepth,(ETA_GROUND(I,J)+Depth(I,J))))**(1-Parameter_n)
                 Omega(I,J)=max(ZERO,(tmp_coeff*(tmp_mos(I,J)-1)+Gamma_rcp)*POR_constant)
             ELSE
                 Omega(I,J)=(-DEPTH(I,J)-ETA_GROUND(I,J))*POR_constant
             ENDIF
      ELSE
         Omega(I,J)=ZERO
      ENDIF
    ENDDO
    ENDDO
END SUBROUTINE CALCULATE_MOISTURE

# endif

SUBROUTINE VEGEROOT_INITIAL
                    
  USE Input_Util
  IMPLICIT NONE

  CHARACTER(LEN=80) FILE_NAME
  INTEGER::LINE
  INTEGER :: ierr

! read parameters from input.txt
      FILE_NAME='input.txt'

      CALL GET_INTEGER_VAL(num_species,FILE_NAME,'num_species',line,ierr)

      IF(ierr==1)THEN
        num_species = 1
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'No num_species defined, USE defalt: 1'
         WRITE(3,'(A80)')'No num_species defined, USE defalt: 1'
      endif
# else
         WRITE(*,'(A80)')'No num_species defined, USE defalt: 1'
         WRITE(3,'(A80)')'No num_species defined, USE defalt: 1'
# endif
       ENDIF

       CALL GET_Float_VAL(n_vege,FILE_NAME,'VegeDynamics_factor',line,ierr)

             IF(ierr==1)THEN
              n_vege = 1.0
# if defined (PARALLEL)
             if (myid.eq.0) THEN
                WRITE(*,'(A80)')'No VegeDynamics_factor defined, USE defalt: 1.0'
                WRITE(3,'(A80)')'No VegeDynamics_factor defined, USE defalt: 1.0'
             endif
# else
                WRITE(*,'(A80)')'No VegeDynamics_factor defined, USE defalt: 1.0'
                WRITE(3,'(A80)')'No VegeDynamics_factor defined, USE defalt: 1.0'
# endif
              ENDIF

      CALL GET_LOGICAL_VAL(HOT_IntialBIOMASS,FILE_NAME,'HOT_IntialBIOMASS',line,ierr)
      CALL GET_LOGICAL_VAL(OUT_Biomass,FILE_NAME,'Biomass',line,ierr)
      CALL GET_LOGICAL_VAL(OUT_avgERoot,FILE_NAME,'AvgEroot',line,ierr)
      CALL GET_LOGICAL_VAL(OUT_Coefficient,FILE_NAME,'Coefficient',line,ierr)
      CALL GET_LOGICAL_VAL(OUT_Transpiration,FILE_NAME,'Transpiration',line,ierr)

      CALL GET_Float_VAL(ET_e,FILE_NAME,'Evaporation',line,ierr)
      IF(ierr==1)THEN
        ET_e = 2.0
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'Evaporation: Use default,', '2.0 mm/d'
         WRITE(3,'(A80)')'Evaporation: Use default,', '2.0 mm/d'
      endif
# else
         WRITE(*,'(A80)')'Evaporation: Use default,', '2.0 mm/d'
         WRITE(3,'(A80)')'Evaporation: Use default,', '2.0 mm/d'
# endif
      ENDIF
      ET_e = ET_e/3600000

      CALL GET_Float_VAL(Tidal_min,FILE_NAME,'Tidal_min',line,ierr)
      IF(ierr==1)THEN
        Tidal_min = 2.0
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'Tidal_min: Use default,', '2.0 m'
         WRITE(3,'(A80)')'Tidal_min: Use default,', '2.0 m'
      endif
# else
         WRITE(*,'(A80)')'Tidal_min: Use default,', '2.0 m'
         WRITE(3,'(A80)')'Tidal_min: Use default,', '2.0 m'
# endif
      ENDIF

! FOR Species1
      CALL GET_Float_VAL(MRD1,FILE_NAME,'Root_Depth1',line,ierr)
      IF(ierr==1)THEN
        MRD1 = 0.4
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'Root_Depth1: Use default,', '0.4 m'
         WRITE(3,'(A80)')'Root_Depth1: Use default,', '0.4 m'
      endif
# else
         WRITE(*,'(A80)')'Root_Depth1: Use default,', '0.4 m'
         WRITE(3,'(A80)')'Root_Depth1: Use default,', '0.4 m'
# endif
      ENDIF

      IF(.NOT.HOT_IntialBIOMASS)THEN

      CALL GET_Float_VAL(Bio_init1,FILE_NAME,'Initial_Biomass1',line,ierr)
      IF(ierr==1)THEN
          Bio_init1 = 0.2
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'Initial_Biomass1: Use default,', '0.2'
         WRITE(3,'(A80)')'Initial_Biomass1: Use default,', '0.2'
      endif
# else
         WRITE(*,'(A80)')'Initial_Biomass1: Use default,', '0.2'
         WRITE(3,'(A80)')'Initial_Biomass1: Use default,', '0.2'
# endif
      ENDIF

      ELSE

      CALL GET_STRING_VAL(IntialBIO1_FILE,FILE_NAME,'IntialBIO1_FILE',line,ierr)

      CALL GET_STRING_VAL(Eta0_FILE,FILE_NAME,'Eta0_FILE',line,ierr)

     ENDIF

      CALL GET_Float_VAL(E14,FILE_NAME,'Anaerobiosis_point1',line,ierr)
      IF(ierr==1)THEN
        E14 = 0.95
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'Anaerobiosis_point1: Use default,', '0.95'
         WRITE(3,'(A40,A40)')'Anaerobiosis_point1: Use default,', '0.95'
      endif
# else
         WRITE(*,'(A40,A40)')'Anaerobiosis_point1: Use default,', '0.95'
         WRITE(3,'(A40,A40)')'Anaerobiosis_point1: Use default,', '0.95'
# endif
      ENDIF

      CALL GET_Float_VAL(E13,FILE_NAME,'Threshold3_1',line,ierr)
      IF(ierr==1)THEN
        E11 = 0.75
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'Threshold3_1: Use default,', '0.75'
         WRITE(3,'(A40,A40)')'Threshold3_1: Use default,', '0.75'
      endif
# else
         WRITE(*,'(A40,A40)')'Threshold3_1: Use default,', '0.75'
         WRITE(3,'(A40,A40)')'Threshold3_1: Use default,', '0.75'
# endif
      ENDIF

      CALL GET_Float_VAL(E12,FILE_NAME,'Threshold2_1',line,ierr)
      IF(ierr==1)THEN
        E12 = 0.65
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'Threshold2_1: Use default,', '0.65'
         WRITE(3,'(A40,A40)')'Threshold2_1: Use default,', '0.65'
      endif
# else
         WRITE(*,'(A40,A40)')'Threshold2_1: Use default,', '0.65'
         WRITE(3,'(A40,A40)')'Threshold2_1: Use default,', '0.65'
# endif
      ENDIF

      CALL GET_Float_VAL(E11,FILE_NAME,'Wilting_point1',line,ierr)
      IF(ierr==1)THEN
        E11 = 0.45
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'Wilting_point1: Use default,', '0.45'
         WRITE(3,'(A40,A40)')'Wilting_point1: Use default,', '0.45'
      endif
# else
         WRITE(*,'(A40,A40)')'Wilting_point1: Use default,', '0.45'
         WRITE(3,'(A40,A40)')'Wilting_point1: Use default,', '0.45'
# endif
      ENDIF

      CALL GET_Float_VAL(ET_t_max1,FILE_NAME,'Transpiration_max1',line,ierr)
      IF(ierr==1)THEN
        ET_t_max1 = 6.0
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'Transpiration_max1: Use default,', '6 mm/d'
         WRITE(3,'(A40,A40)')'Transpiration_max1: Use default,', '6 mm/d'
      endif
# else
         WRITE(*,'(A40,A40)')'Transpiration_max1: Use default,', '6 mm/d'
         WRITE(3,'(A40,A40)')'Transpiration_max1: Use default,', '6 mm/d'
# endif
      ENDIF
       ET_t_max1 = ET_t_max1/3600000

      CALL GET_Float_VAL(Beta1,FILE_NAME,'Intra_competition1',line,ierr)
      IF(ierr==1)THEN
        Beta1 = 1.0
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'Intra_competition1: Use default,', '1.0'
         WRITE(3,'(A40,A40)')'Intra_competition1: Use default,', '1.0'
      endif
# else
         WRITE(*,'(A40,A40)')'Intra_competition1: Use default,', '1.0'
         WRITE(3,'(A40,A40)')'Intra_competition1: Use default,', '1.0'
# endif
      ENDIF

      CALL GET_Float_VAL(growth1,FILE_NAME,'Intrisic_growth1',line,ierr)
      IF(ierr==1)THEN
        growth1 = 0.0000001
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'Intrisic_growth1: Use default,', '1.0*10^-7'
         WRITE(3,'(A40,A40)')'Intrisic_growth1: Use default,', '1.0*10^-7'
      endif
# else
         WRITE(*,'(A40,A40)')'Intrisic_growth1: Use default,', '1.0*10^-7'
         WRITE(3,'(A40,A40)')'Intrisic_growth1: Use default,', '1.0*10^-7'
# endif
      ENDIF

     IF(E11>E12 .OR. E12>E13 .OR. E13>E14 )THEN
         WRITE(*,'(A40,A40)')'ERROR: WRONG Threshold values!'
         WRITE(3,'(A40,A40)')'ERROR: WRONG Threshold values!'
     ENDIF

!write in LOG.txt
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A50)') '-----for VEGEROOT ----------'
       if(myid.eq.0) WRITE(3,'(A22,F10.2)') 'VegeDynamics_factor =', n_vege
       if(myid.eq.0) WRITE(3,'(A19,F15.10)') 'Evaporation =', ET_e
       if(myid.eq.0) WRITE(3,'(A19,F12.2)') 'Tidal_min =', Tidal_min
       if(myid.eq.0) WRITE(3,'(A50)') '! For species 1 '
       if(myid.eq.0) WRITE(3,'(A19,F12.2)') 'Root_Depth =', MRD1
       if(myid.eq.0) WRITE(3,'(A19,F12.5)') 'Initial_Biomass =', Bio_init1
       if(myid.eq.0) WRITE(3,'(A19,F12.2)') 'E1 =', E11
       if(myid.eq.0) WRITE(3,'(A19,F12.2)') 'E2 =', E12
       if(myid.eq.0) WRITE(3,'(A19,F12.2)') 'E3 =', E13
       if(myid.eq.0) WRITE(3,'(A19,F12.2)') 'E4 =', E14
       if(myid.eq.0) WRITE(3,'(A19,F15.10)') 'Transpiration_max =', ET_t_max1
       if(myid.eq.0) WRITE(3,'(A19,F12.2)') 'Intra_competition =', Beta1
       if(myid.eq.0) WRITE(3,'(A19,F12.10)') 'Intrisic_growth =', growth1
# else
       WRITE(3,'(A50)') '-----for VEGEROOT ----------'
       WRITE(3,'(A22,F10.2)') 'VegeDynamics_factor =', n_vege
       WRITE(3,'(A19,F15.10)') 'Evaporation =', ET_e
       WRITE(3,'(A19,F12.2)') 'Tidal_min =', Tidal_min
       WRITE(3,'(A50)') '! For species 1 '
       WRITE(3,'(A19,F12.2)') 'Root_Depth =', MRD1
       WRITE(3,'(A19,F12.2)') 'Initial_Biomass =', Bio_init1
       WRITE(3,'(A19,F12.2)') 'E1 =', E11
       WRITE(3,'(A19,F12.2)') 'E2 =', E12
       WRITE(3,'(A19,F12.2)') 'E3 =', E13
       WRITE(3,'(A19,F12.2)') 'E4 =', E14
       WRITE(3,'(A19,F15.10)') 'Transpiration_max =', ET_t_max1
       WRITE(3,'(A19,F12.2)') 'Intra_competition =', Beta1
       WRITE(3,'(A19,F12.10)') 'Intrisic_growth =', growth1
# endif

! FOR Species2
      CALL GET_Float_VAL(MRD2,FILE_NAME,'Root_Depth2',line,ierr)
      IF(ierr==1)THEN
        MRD2 = 0.2
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'Root_Depth2: Use default,', '0.2 m'
         WRITE(3,'(A80)')'Root_Depth2: Use default,', '0.2 m'
      endif
# else
         WRITE(*,'(A80)')'Root_Depth2: Use default,', '0.2 m'
         WRITE(3,'(A80)')'Root_Depth2: Use default,', '0.2 m'
# endif
      ENDIF

      IF(.NOT.HOT_IntialBIOMASS)THEN

      CALL GET_Float_VAL(Bio_init2,FILE_NAME,'Initial_Biomass2',line,ierr)
      IF(ierr==1)THEN
          Bio_init2 = 0.2
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A80)')'Initial_Biomass2: Use default,', '0.2'
         WRITE(3,'(A80)')'Initial_Biomass2: Use default,', '0.2'
      endif
# else
         WRITE(*,'(A80)')'Initial_Biomass2: Use default,', '0.2'
         WRITE(3,'(A80)')'Initial_Biomass2: Use default,', '0.2'
# endif
      ENDIF

      ELSE
   
         CALL GET_STRING_VAL(IntialBIO2_FILE,FILE_NAME,'IntialBIO2_FILE',line,ierr)

      ENDIF
      

      CALL GET_Float_VAL(E24,FILE_NAME,'Anaerobiosis_point2',line,ierr)
      IF(ierr==1)THEN
        E24 = 0.75
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'Anaerobiosis_point2: Use default,', '0.75'
         WRITE(3,'(A40,A40)')'Anaerobiosis_point2: Use default,', '0.75'
      endif
# else
         WRITE(*,'(A40,A40)')'Anaerobiosis_point2: Use default,', '0.75'
         WRITE(3,'(A40,A40)')'Anaerobiosis_point2: Use default,', '0.75'
# endif
      ENDIF

      CALL GET_Float_VAL(E23,FILE_NAME,'Threshold3_2',line,ierr)
      IF(ierr==1)THEN
        E23 = 0.55
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'Threshold3_2: Use default,', '0.55'
         WRITE(3,'(A40,A40)')'Threshold3_2: Use default,', '0.55'
      endif
# else
         WRITE(*,'(A40,A40)')'Threshold3_2: Use default,', '0.55'
         WRITE(3,'(A40,A40)')'Threshold3_2: Use default,', '0.55'
# endif
      ENDIF

      CALL GET_Float_VAL(E22,FILE_NAME,'Threshold2_2',line,ierr)
      IF(ierr==1)THEN
        E22 = 0.4
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'Threshold2_2: Use default,', '0.4'
         WRITE(3,'(A40,A40)')'Threshold2_2: Use default,', '0.4'
      endif
# else
         WRITE(*,'(A40,A40)')'Threshold2_2: Use default,', '0.4'
         WRITE(3,'(A40,A40)')'Threshold2_2: Use default,', '0.4'
# endif
      ENDIF

      CALL GET_Float_VAL(E21,FILE_NAME,'Wilting_point2',line,ierr)
      IF(ierr==1)THEN
        E21 = 0.2
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'Wilting_point2: Use default,', '0.2'
         WRITE(3,'(A40,A40)')'Wilting_point2: Use default,', '0.2'
      endif
# else
         WRITE(*,'(A40,A40)')'Wilting_point2: Use default,', '0.2'
         WRITE(3,'(A40,A40)')'Wilting_point2: Use default,', '0.2'
# endif
      ENDIF

      CALL GET_Float_VAL(ET_t_max2,FILE_NAME,'Transpiration_max1',line,ierr)
      IF(ierr==1)THEN
        ET_t_max2 = 6.0
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'Transpiration_max2: Use default,', '6 mm/d'
         WRITE(3,'(A40,A40)')'Transpiration_max2: Use default,', '6 mm/d'
      endif
# else
         WRITE(*,'(A40,A40)')'Transpiration_max2: Use default,', '6 mm/d'
         WRITE(3,'(A40,A40)')'Transpiration_max2: Use default,', '6 mm/d'
# endif
      ENDIF
        ET_t_max2 = ET_t_max2/3600000

      CALL GET_Float_VAL(Beta2,FILE_NAME,'Intra_competition2',line,ierr)
      IF(ierr==1)THEN
        Beta2 = 1.0
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'Intra_competition2: Use default,', '1.0'
         WRITE(3,'(A40,A40)')'Intra_competition2: Use default,', '1.0'
      endif
# else
         WRITE(*,'(A40,A40)')'Intra_competition2: Use default,', '1.0'
         WRITE(3,'(A40,A40)')'Intra_competition2: Use default,', '1.0'
# endif
      ENDIF

      CALL GET_Float_VAL(growth2,FILE_NAME,'Intrisic_growth2',line,ierr)
      IF(ierr==1)THEN
        growth2 = 0.0000001
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'Intrisic_growth2: Use default,', '1.0*10^-7'
         WRITE(3,'(A40,A40)')'Intrisic_growth2: Use default,', '1.0*10^-7'
      endif
# else
         WRITE(*,'(A40,A40)')'Intrisic_growth2: Use default,', '1.0*10^-7'
         WRITE(3,'(A40,A40)')'Intrisic_growth2: Use default,', '1.0*10^-7'
# endif
      ENDIF

      IF(E21>E22 .OR. E22>E23 .OR. E23>E24 )THEN
          WRITE(*,'(A40,A40)')'ERROR: WRONG Threshold values!'
          WRITE(3,'(A40,A40)')'ERROR: WRONG Threshold values!'
      ENDIF

!write in LOG.txt
     IF(num_species==2)THEN
# if defined (PARALLEL)
       if(myid.eq.0) WRITE(3,'(A50)') '! For species 2 '
       if(myid.eq.0) WRITE(3,'(A19,F12.2)') 'Root_Depth =', MRD2
       if(myid.eq.0) WRITE(3,'(A19,F12.2)') 'Initial_Biomass =', Bio_init2
       if(myid.eq.0) WRITE(3,'(A19,F12.2)') 'E1 =', E21
       if(myid.eq.0) WRITE(3,'(A19,F12.2)') 'E2 =', E22
       if(myid.eq.0) WRITE(3,'(A19,F12.2)') 'E3 =', E23
       if(myid.eq.0) WRITE(3,'(A19,F12.2)') 'E4 =', E24
       if(myid.eq.0) WRITE(3,'(A19,F15.10)') 'Transpiration_max =', ET_t_max2
       if(myid.eq.0) WRITE(3,'(A19,F12.2)') 'Intra_competition =', Beta2
       if(myid.eq.0) WRITE(3,'(A19,F12.2)') 'Intrisic_growth =', growth2
# else
       WRITE(3,'(A50)') '! For species 2 '
       WRITE(3,'(A19,F12.2)') 'Root_Depth =', MRD2
       WRITE(3,'(A19,F12.2)') 'Initial_Biomass =', Bio_init2
       WRITE(3,'(A19,F12.2)') 'E1 =', E21
       WRITE(3,'(A19,F12.2)') 'E2 =', E22
       WRITE(3,'(A19,F12.2)') 'E3 =', E23
       WRITE(3,'(A19,F12.2)') 'E4 =', E24
       WRITE(3,'(A19,F15.10)') 'Transpiration_max =', ET_t_max2
       WRITE(3,'(A19,F12.10)') 'Intra_competition =', Beta2
       WRITE(3,'(A19,F12.10)') 'Intrisic_growth =', growth2
# endif
     ENDIF

      CALL GET_Float_VAL(Beta12,FILE_NAME,'Inter_competition',line,ierr)
      IF(ierr==1)THEN
        Beta12 = 0.0
# if defined (PARALLEL)
      if (myid.eq.0) THEN
         WRITE(*,'(A40,A40)')'Inter_competition: Use default,', '0.0'
         WRITE(3,'(A40,A40)')'Inter_competition: Use default,', '0.0'
      endif
# else
         WRITE(*,'(A40,A40)')'Inter_competition: Use default,', '0.0'
         WRITE(3,'(A40,A40)')'Inter_competition: Use default,', '0.0'
# endif
      ENDIF

! ALLOCATION

        ALLOCATE(Bio1(Mloc,Nloc),Bio2(Mloc,Nloc),DeltaBio1(Mloc,Nloc),DeltaBio2(Mloc,Nloc),&
                 Bio10(Mloc,Nloc),Bio20(Mloc,Nloc),&
            avgEr1(Mloc,Nloc),avgEr2(Mloc,Nloc),ET_t1(Mloc,Nloc),ET_t2(Mloc,Nloc), &
            Coeff_uptake1(Mloc,Nloc),Coeff_uptake2(Mloc,Nloc), &
            Coeff_capct1(Mloc,Nloc),Coeff_capct2(Mloc,Nloc))


! assign variables

       IF(.NOT.HOT_IntialBIOMASS)THEN
        DO J = 1, Nloc
        DO I = 1, Mloc
            IF(DEPTH(I,J).le.Tidal_min)THEN
                Bio10(I,J) = Bio_init1
            ELSE
                Bio10(I,J) = ZERO
           ENDIF
        ENDDO
        ENDDO

       ELSE
         CALL GetFile(IntialBIO1_FILE,Bio10)
         CALL GetFile(ETA0_FILE,ETA_ground)
       ENDIF


        avgEr1 = ZERO
        ET_t1 = ZERO
        Coeff_uptake1 = ZERO
        Coeff_capct1 = ZERO
        DeltaBio1 = ZERO

        IF(num_species==2)THEN

           IF(.NOT.HOT_IntialBIOMASS)THEN
            DO J = 1, Nloc
            DO I = 1, Mloc
                IF(DEPTH(I,J).le.Tidal_min)THEN
                    Bio20(I,J) = Bio_init2
                ELSE
                    Bio20(I,J) = ZERO
               ENDIF
            ENDDO
            ENDDO

           ELSE
             CALL GetFile(IntialBIO2_FILE,Bio20)
           ENDIF
            
        avgEr2 = ZERO
        ET_t2 = ZERO
        Coeff_uptake2 = ZERO
        Coeff_capct2 = ZERO
        DeltaBio2 = ZERO
 
         ENDIF


END SUBROUTINE VEGEROOT_INITIAL

# if defined (MULTIRELATIONSHIP)

SUBROUTINE CALCULATE_avgEr(MRD,avgEr,ETA_GROUND,DEPTH)
! Calculate the mean soil saturation over the mean root depth

     
     IMPLICIT NONE
     REAL(SP),INTENT(IN)::MRD,ETA_GROUND,DEPTH
     REAL(SP),INTENT(OUT)::avgEr
     REAL(SP):: tmp_r, SUM
      IF(-DEPTH-MRD .GE. ETA_GROUND+gamma_rcp)THEN
          CALL INTEGRATION(ETA_GROUND,-DEPTH-MRD,-DEPTH,ddep,SUM)
          avgEr=SUM/MRD
      ELSEIF(-DEPTH .GT. ETA_GROUND+gamma_rcp) THEN
          CALL INTEGRATION(ETA_GROUND,ETA_GROUND+gamma_rcp,-DEPTH,ddep,SUM)
          avgEr=(SUM+DEPTH+ETA_GROUND+gamma_rcp+MRD)/MRD
      ELSE
          avgEr=1.0
      ENDIF
          avgEr= min(1.0, avgEr)
          avgEr= max(ZERO, avgEr)
        
END SUBROUTINE CALCULATE_avgEr

# else

SUBROUTINE CALCULATE_avgEr(MRD,avgEr,ETA_GROUND,DEPTH)
! Calculate the mean soil saturation over the mean root depth
     IMPLICIT NONE
     REAL(SP),INTENT(IN)::MRD,ETA_GROUND,DEPTH
     REAL(SP),INTENT(OUT)::avgEr
     REAL(SP):: tmp_r, tmp_mos2
            
            tmp_r=ETA_GROUND+Depth+MRD
            tmp_mos2=(-Gamma*min(-MinDepth,(ETA_GROUND+Depth)))**(1-Parameter_n)
            IF (tmp_r<-gamma_rcp) THEN
                avgEr= tmp_coeff*(tmp_mos2-(-gamma*min(tmp_r,-MinDepth))**(1-parameter_n))/MRD
            ELSEIF(ETA_GROUND+Depth<-gamma_rcp)THEN
                avgEr= (tmp_coeff*(tmp_mos2-1)+tmp_r+gamma_rcp)/MRD
            ELSE
                avgEr= 1.0
            ENDIF
                avgEr= min(1.0, avgEr)
                avgEr= max(ZERO, avgEr)
        
END SUBROUTINE CALCULATE_avgEr

# endif

SUBROUTINE GET_ROOTUPTAKE_COEFF(e1,e2,e3,e4,avgEr,coeff_uptake)
      USE GLOBAL, ONLY : Mloc,Nloc
      REAL(SP),INTENT(IN):: e1,e2,e3,e4
      REAL(SP),INTENT(IN),DIMENSION(Mloc,Nloc)::avgEr
      REAL(SP),INTENT(OUT),DIMENSION(Mloc,Nloc)::coeff_uptake

        coeff_uptake = ZERO

        DO J = 1, Nloc
        DO I = 1, Mloc
          IF(MASK(I,J)==0)THEN
            IF(avgEr(I,J)<=e1)THEN
              coeff_uptake(I,J)=ZERO
            ELSEIF(avgEr(I,J)<=e2)THEN
              coeff_uptake(I,J)=(avgEr(I,J)-e1)/(e2-e1)
            ELSEIF(avgEr(I,J)<=e3)THEN
              coeff_uptake(I,J)=1.0
            ELSEIF(avgEr(I,J)<=e4)THEN
              coeff_uptake(I,J)=(avgEr(I,J)-e4)/(e3-e4)
            ELSE
              coeff_uptake(I,J)=ZERO
            ENDIF
          ENDIF
          coeff_uptake(I,J)=max(ZERO,coeff_uptake(I,J))
          coeff_uptake(I,J)=min(1.0,coeff_uptake(I,J))
        ENDDO
        ENDDO
END SUBROUTINE GET_ROOTUPTAKE_COEFF

SUBROUTINE GET_CARRYCAPACITY_COEFF(e1,e4,avgEr,coeff_capct)
       USE GLOBAL, ONLY : Mloc,Nloc
       REAL(SP),INTENT(IN):: e1,e4
       REAL(SP),INTENT(IN),DIMENSION(Mloc,Nloc)::avgEr
       REAL(SP),INTENT(OUT),DIMENSION(Mloc,Nloc)::coeff_capct

        coeff_capct = ZERO

        DO J = 1, Nloc
        DO I = 1, Mloc
          IF(MASK(I,J)==0)THEN
              IF(avgEr(I,J)<=e1)THEN
                 coeff_capct(I,J)=ZERO
              ELSEIF(avgEr(I,J)<=e4)THEN
                 coeff_capct(I,J)=1.0
              ELSE
                 coeff_capct(I,J)=ZERO
              ENDIF
           ENDIF
         ENDDO
         ENDDO
END SUBROUTINE GET_CARRYCAPACITY_COEFF

SUBROUTINE UPDATE_BIOMASS

  !    REAL(SP),DIMENSION(:,:),ALLOCATABLE::tmpBio1,tmpBio2
      REAL(SP):: tmp1, tmp2
 !     ALLOCATE(tmpBio1(Mloc,Nloc),tmpBio2(Mloc,Nloc))


 !     tmpBio1=DeltaBio1
 !     tmpBio2=DeltaBio2

   DO J = 1, Nloc
   DO I = 1, Mloc
        IF(DEPTH(I,J).le.Tidal_min)THEN
              IF(num_species==1)THEN
                 
                 tmp1=DT*Bio1(I,J)*growth1*(Coeff_capct1(I,J)-Beta1*Bio1(I,J))
                 DeltaBio1(I,J)=DeltaBio1(I,J)+tmp1*n_vege
                 Bio1(I,J)=Bio10(I,J)+DeltaBio1(I,J)
                
!                Bio1(I,J)=Bio10(I,J)+DT*Bio10(I,J)*growth1*(Coeff_capct1(I,J)-Beta1*Bio10(I,J))
                Bio1(I,J)=max(Bio1(I,J),ZERO)

              ELSEIF(num_species==2)THEN

                 tmp1=DT*Bio1(I,J)*growth1*(Coeff_capct1(I,J)-Beta1*Bio1(I,J)-Beta12*Bio2(I,J))
                 DeltaBio1(I,J)=DeltaBio1(I,J)+tmp1*n_vege

                 tmp2=DT*Bio2(I,J)*growth2*(Coeff_capct2(I,J)-Beta2*Bio2(I,J)-Beta12*Bio1(I,J))
                 DeltaBio2(I,J)=DeltaBio2(I,J)+tmp2*n_vege
                 
                 Bio1(I,J)=Bio10(I,J)+DeltaBio1(I,J)
                 Bio2(I,J)=Bio20(I,J)+DeltaBio2(I,J)

                Bio1(I,J)=max(Bio1(I,J),ZERO)
                Bio2(I,J)=max(Bio2(I,J),ZERO)
              ENDIF
        ENDIF
   ENDDO
   ENDDO
END SUBROUTINE UPDATE_BIOMASS

# if defined (MULTIRELATIONSHIP)
SUBROUTINE Model_Func(Z,ETA_GROUND,E_dot)
 REAL,INTENT(IN):: Z,ETA_GROUND
 REAL(SP),INTENT(OUT):: E_dot
! For B-C Model
 IF(MODEL_TYPE == 'B-C')THEN
   E_dot = (Gamma*(Z-ETA_GROUND))**(-parameter_n)
 ENDIF
END SUBROUTINE

SUBROUTINE INTEGRATION(ETA_GROUND,LOWBOUND,UPBOUND,ddep,SUM)
 IMPLICIT NONE
 REAL(SP),INTENT(IN):: ETA_GROUND,LOWBOUND,UPBOUND,ddep
 REAL(SP),INTENT(OUT):: SUM
 REAL(SP):: Z
 INTEGER :: num,i
    IF(LOWBOUND .LT. UPBOUND)THEN
       SUM = ZERO
       Z = LOWBOUND
       DO WHILE(Z .ge. LOWBOUND .and. Z .le. UPBOUND )
          CALL Model_Func(Z,ETA_GROUND,E_dot)
          SUM = SUM + E_dot*ddep
          Z = Z + ddep
       ENDDO
       SUM = max(SUM,ZERO)
    ELSE
       WRITE(*,'(A80)')'UPBOUND should be bigger than LOWBOUND'
       STOP
    ENDIF

END SUBROUTINE INTEGRATION

# endif

END MODULE SUBSURFACE_MODULE

# endif
! end the module

