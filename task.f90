MODULE AFUN_d
  USE CUDAFOR
  IMPLICIT NONE
  CONTAINS
  
  ATTRIBUTES (DEVICE) REAL(8) FUNCTION ACOF1(X1, Y1)
    IMPLICIT NONE
    REAL(8) :: ACOF
    REAL(8), INTENT(IN) :: X1, Y1

    IF (X1  == 0.0_8) THEN 
      ACOF = 0.0_8 
    ELSE
      ACOF = DMAX1(0.0_8, (1.0_8 - 0.1_8 * DABS(Y1 / X1))**5)
    END IF

    ACOF1 = X1*ACOF + DMAX1(0.0_8, Y1)

  END FUNCTION ACOF1


  ATTRIBUTES (DEVICE) REAL(8) FUNCTION Q(FLOW, x, a, b, c, fa, fb, fc)
    IMPLICIT NONE
    REAL(8) :: ka, kc
    REAL(8), INTENT(IN) :: x, a, b, c, fa, fb, fc, FLOW

    !IF (FLOW == 0.0_8) THEN 
    !  Q = 0.0_8 
    !ELSE
    !  ka = (x - b) * (x - c) / ((b - a) * (c - a) ) 
    !  kc = (x - b) * (x - a) / ((c - b) * (c - a) )
    !  Q = DMAX1(0.0_8, FLOW) * (ka * (fa - fb) + kc * (fc - fb)) 
    !END IF

    IF ((FLOW < 0.0_8) .OR. (FLOW > 0.0_8)) THEN  
      ka = (x - b) * (x - c) / ((b - a) * (c - a) ) 
      kc = (x - b) * (x - a) / ((c - b) * (c - a) )
      Q = DMAX1(0.0_8, FLOW) * (ka * (fa - fb) + kc * (fc - fb)) 
    ELSE
      Q = 0.0_8 
    END IF

  END FUNCTION Q


  ATTRIBUTES (DEVICE) REAL(8) FUNCTION Qsc(FLOW, x, a, b, c, fa, fb, fc)
    IMPLICIT NONE
    REAL(8) :: ka, kc
    REAL(8), INTENT(IN) :: x, a, b, c, fa, fb, fc, FLOW

    !IF (FLOW == 0.0_8) THEN 
    !  Qsc = 0.0_8 
    !ELSE
    !  ka = (x - b) * (x - c) / ((b - a) * (c - a) ) 
    !  kc = (x - b) * (x - a) / ((c - b) * (c - a) )
    !  Qsc = DMAX1(0.0_8, FLOW) * (fb + ka * (fa - fb) + kc * (fc - fb)) 
    !END IF

    IF ((FLOW < 0.0_8) .OR. (FLOW > 0.0_8)) THEN
      ka = (x - b) * (x - c) / ((b - a) * (c - a) ) 
      kc = (x - b) * (x - a) / ((c - b) * (c - a) )
      Qsc = DMAX1(0.0_8, FLOW) * (fb + ka * (fa - fb) + kc * (fc - fb)) 
    ELSE
      Qsc = 0.0_8
    END IF

  END FUNCTION Qsc

  ATTRIBUTES (DEVICE) REAL(8) FUNCTION WEIGHT_AVG(a1, a2, coeff)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: a1, a2, coeff
    REAL(8) :: WEIGHT_AVG
    
    WEIGHT_AVG = coeff * a1 + (1.0_8 - coeff) * a2

  END FUNCTION WEIGHT_AVG

END MODULE AFUN_d


MODULE USER
  USE CUDAFOR
  IMPLICIT NONE
  CONTAINS

  !------------------------------------------------------------------
  SUBROUTINE SET_TASK_PARAMS
    USE VAR
    USE VAR_d, ONLY: GRD_IJK, GRD_IJ, GRD_JK, GRD_IK, GRD_I, TBL_IJK, TBL_IJ, TBL_JK, TBL_IK, TBL_I
    USE VAR_PSI_OUT_d, ONLY: PSI_OUT_GRD_IJK => GRD_IJK, PSI_OUT_GRD_IJ => GRD_IJ, &
                             PSI_OUT_GRD_JK => GRD_JK, PSI_OUT_GRD_IK => GRD_IK, &
                             PSI_OUT_TBL_IJK => TBL_IJK, PSI_OUT_TBL_IJ => TBL_IJ, &
                             PSI_OUT_TBL_JK => TBL_JK, PSI_OUT_TBL_IK => TBL_IK, &
                             N1_PSI_OUT => N1 
    USE VAR_PSI_IN_d, ONLY:  PSI_IN_GRD_IJK => GRD_IJK, PSI_IN_GRD_IJ => GRD_IJ, &
                             PSI_IN_GRD_JK => GRD_JK, PSI_IN_GRD_IK => GRD_IK, &
                             PSI_IN_TBL_IJK => TBL_IJK, PSI_IN_TBL_IJ => TBL_IJ, &
                             PSI_IN_TBL_JK => TBL_JK, PSI_IN_TBL_IK => TBL_IK, &
                             N1_PSI_IN => N1 
    IMPLICIT NONE

    TIME = 0.0_8; DT = 1.0D-4
    LAST = 1; IWRITE = 1000
    
    !Геометрические параметры сферического слоя   
    X1_0 = 5.0_8 / 13.0_8; XL1 = 1.0_8 !XL1 - не радиус а, толщина.
    X2_0 = 0.0_8; XL2 = 0.5_8 * PI
    X3_0 = 0.0_8; XL3 = 2.0_8 * PI

    r_i = X1_0; r_o = X1_0 + XL1
   
    !Безразмерные физические параметры
    E = 1.0D-3
    Ranum = 140.0_8
    Pm = 5.0_8
    Pr = 1.0_8

   !физические параметры кода

    viscm = 1.0_8 / Pm
    visk = E
    dens = E
    thc = 1.0_8 / Pr
    rho_c = 1.0_8
    omega_ic = - 2.6595_8
    
    NTIMES_UVW = 3
    NTIMES_P   = 15
    NTIMES_PC  = 30
    NTIMES_T   = 3
    NTIMES_F1  = 20
    NTIMES_F2  = 20
    NTIMES_F3  = 20

    TOL_UVW = 1.0D-5  
    TOL_P   = 1.0D-10 
    TOL_T   = 1.0D-5
    tolPSI  = 1.D-6
    tolPSI_C  = 2.D-9
    
    FirstCallMP = .TRUE.
 
    !--------------------------------

    TBL_IJK = DIM3(4, 8, 8)
    GRD_IJK = DIM3(CEILING(REAL(N1) / TBL_IJK%x), &
                  CEILING(REAL(N2) / TBL_IJK%y), &
                  CEILING(REAL(N3) / TBL_IJK%z))

    TBL_IJ = DIM3(4, 8, 1)
    GRD_IJ = DIM3(CEILING(REAL(N1) / TBL_IJ%x), &
                 CEILING(REAL(N2) / TBL_IJ%y), 1)

    TBL_JK = DIM3(4, 8, 1)
    GRD_JK = DIM3(CEILING(REAL(N2) / TBL_JK%x), &
                 CEILING(REAL(N3) / TBL_JK%y), 1)

    TBL_IK = DIM3(4, 8, 1)
    GRD_IK   = DIM3(CEILING(REAL(N1) / TBL_IK%x), &
                    CEILING(REAL(N3) / TBL_IK%y), 1)

    TBL_I = DIM3(4, 1, 1)
    GRD_I = DIM3(CEILING(REAL(N1) / TBL_I%x), 1, 1)
                    
    !----------------------------------
    
    PSI_OUT_TBL_IJK = DIM3(4, 8, 8)
    PSI_OUT_GRD_IJK = DIM3(CEILING(REAL(N1_PSI_OUT) / PSI_OUT_TBL_IJK%x), &
                           CEILING(REAL(N2)         / PSI_OUT_TBL_IJK%y), &
                           CEILING(REAL(N3)         / PSI_OUT_TBL_IJK%z))

    PSI_OUT_TBL_IJ = DIM3(4, 8, 1)
    PSI_OUT_GRD_IJ = DIM3(CEILING(REAL(N1_PSI_OUT) / PSI_OUT_TBL_IJ%x), &
                          CEILING(REAL(N2)         / PSI_OUT_TBL_IJ%y), 1)

    PSI_OUT_TBL_JK = DIM3(4, 8, 1)
    PSI_OUT_GRD_JK = DIM3(CEILING(REAL(N2)         / PSI_OUT_TBL_JK%x), &
                          CEILING(REAL(N3)         / PSI_OUT_TBL_JK%y), 1)

    PSI_OUT_TBL_IK = DIM3(4, 8, 1)
    PSI_OUT_GRD_IK = DIM3(CEILING(REAL(N1_PSI_OUT) / PSI_OUT_TBL_IK%x), &
                          CEILING(REAL(N3)         / PSI_OUT_TBL_IK%y), 1)
                          
    !--------------------------------------

    PSI_IN_TBL_IJK = DIM3(4, 8, 8)
    PSI_IN_GRD_IJK = DIM3(CEILING(REAL(N1_PSI_IN) / PSI_IN_TBL_IJK%x), &
                           CEILING(REAL(N2)       / PSI_IN_TBL_IJK%y), &
                           CEILING(REAL(N3)       / PSI_IN_TBL_IJK%z))

    PSI_IN_TBL_IJ = DIM3(4, 8, 1)
    PSI_IN_GRD_IJ = DIM3(CEILING(REAL(N1_PSI_IN) / PSI_IN_TBL_IJ%x), &
                          CEILING(REAL(N2)       / PSI_IN_TBL_IJ%y), 1)

    PSI_IN_TBL_JK = DIM3(4, 8, 1)
    PSI_IN_GRD_JK = DIM3(CEILING(REAL(N2)        / PSI_IN_TBL_JK%x), &
                          CEILING(REAL(N3)       / PSI_IN_TBL_JK%y), 1)

    PSI_IN_TBL_IK = DIM3(4, 8, 1)
    PSI_IN_GRD_IK = DIM3(CEILING(REAL(N1_PSI_IN) / PSI_IN_TBL_IK%x), &
                          CEILING(REAL(N3)       / PSI_IN_TBL_IK%y), 1)



!     WRITE(*, *) 'GRD_IJK: ', CEILING(REAL(N1) / TBL_IJK%x), &
!                              CEILING(REAL(N2) / TBL_IJK%y), &
!                              CEILING(REAL(N3) / TBL_IJK%z) 
! 
!     WRITE(*, *) 'GRD_IJ: ', CEILING(REAL(N1) / TBL_IJ%x), &
!                             CEILING(REAL(N2) / TBL_IJ%y)
! 
!     WRITE(*, *) 'GRD_JK: ', CEILING(REAL(N2) / TBL_JK%x), &
!                             CEILING(REAL(N3) / TBL_JK%y)
! 
!     WRITE(*, *) 'GRD_IK: ', CEILING(REAL(N1) / TBL_IK%x), &
!                             CEILING(REAL(N3) / TBL_IK%y)


         
  END SUBROUTINE SET_TASK_PARAMS
 
  !------------------------------------------------------------------

  SUBROUTINE START
    USE VAR
    USE VAR_PSI_IN, ONLY: PSI_IN
    USE VAR_PSI_OUT, ONLY: PSI_OUT
    IMPLICIT NONE
    INTEGER :: I, J, K
    REAL(8) :: H1, H2, H3
    REAL(8) :: AMPLT, XR, XR1, a_Br, a_Btheta, a_Bphi, &
               b1_Br, b2_Br, b3_Br, b4_Br 


    B1 = 0.0_8; F1 = 0.0_8; U1 = 0.0_8 
    B2 = 0.0_8; F2 = 0.0_8; U2 = 0.0_8
    B3 = 0.0_8; F3 = 0.0_8; U3 = 0.0_8

    PSI_IN = 0.0_8; PSI_OUT = 0.0_8
    PSI = 0.0_8; PSI_C = 0.0_8
    P = 0.0_8; T = 1.0_8 
  
    AMPLT = 21.0_8 / DSQRT(17920.0_8 * PI)

    DO K = 1, N3
      DO J = 2, N2
        DO I = N1_core, N1
          XR = 2.0_8 * X1(I) - 2.0_8 * r_i - 1.0_8
          XR1 = 1.0_8 - 3.0_8 * (XR**2) + 3.0_8 * (XR**4) - (XR**6)
          T(I,J,K) = AMPLT * XR1 * (sinX2(J)**4) * DCOS(3.0_8 * X3(K)) - r_i + r_i * r_o / X1(I)
        ENDDO 
      ENDDO
    ENDDO


    DO K = 1, N3
      DO J = 2, N2
        T(N1, J, K) = 0.0_8
      ENDDO
    ENDDO
    
    DO K = 2, N3-1
      DO I = 1, N1
        U2(I, N2, K) = 0.0_8
      ENDDO 
    ENDDO


    !DO K = 2, N3-1 
    !  DO J = 2, N2-1
    !    DO I = 2, N1_core-1
    !      U3(I,J,K) = omega_ic * X1(I) * sinX2(J)
    !    ENDDO
    !  ENDDO
    !ENDDO


    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 3, N1
          H1 = 1.0_8    
          B1(I,J,K) = 5.0_8 * DCOS(X2(J)) * (4.0_8 * r_o - 3.0_8 * Xs1(I)) / (3.0_8 + r_o)
          F1(I,J,K) = H1 * B1(I,J,K) 
        ENDDO
      ENDDO
    ENDDO

    
    DO K = 2, N3-1
      DO I = 3, N1
        B1(I,N2,K) = 0.0_8
        F1(I,N2,K) = 0.0_8
      ENDDO
    ENDDO



    DO K = 2, N3-1
      DO J = 3, N2
        DO I = 2, N1
          H2 = X1(I)
          B2(I,J,K) = 5.0_8 * DSIN(Xs2(J)) * (9.0_8 * X1(I) - 8.0_8 * r_o) / (2.0_8 * r_o + 6.0_8)
          F2(I,J,K) = H2 * B2(I,J,K) 
        ENDDO
      ENDDO
    ENDDO


    DO K = 2, N3-1 
      DO J = 2, N2-1
        DO I = 2, N1
          H3 = X1(I) * sinX2(J)
          B3(I,J,K) = 5.0_8 * DSIN(PI * X1(I) / r_o) * DSIN(2.0_8 * X2(J))
          F3(I,J,K) = H3 * B3(I,J,K)
        ENDDO
      ENDDO
    ENDDO

    
    DO K = 2, N3-1 
      DO I = 2, N1
        B3(I,N2,K) = 0.0_8
        F3(I,N2,K) = 0.0_8
      ENDDO
    ENDDO
    

    OPEN(UNIT = 100, FILE = 'STR3D_V.DAT', STATUS = 'OLD', FORM = 'UNFORMATTED', ERR = 100)
    READ(100) F1, F2, F3, B1, B2, B3, U1, U2, U3, P, T, PSI_C, PSI_OUT, PSI_IN, TIME
    CLOSE(UNIT = 100)

    WRITE(*,*) ' FILE WAS READ. PROCESS CONTINUES.'
    RETURN
    100 WRITE(*,*) ' NO DATA FILE. PROCESS STARTS FROM ZERO.'
 

  END SUBROUTINE START
  
  !------------------------------------------------------------------
   
  
  SUBROUTINE BOUND
  USE VAR
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: H1, H2, H3, R_PI, PI_R, EXP_M

  END SUBROUTINE BOUND

  !------------------------------------------------------------------
  SUBROUTINE SAVE_DATA
  USE VAR
  USE VAR_PSI_OUT, ONLY: PSI_OUT
  USE VAR_PSI_IN, ONLY: PSI_IN
  USE VAR_d, ONLY: F1_d => F1, F2_d => F2, F3_d => F3, &
                   P_d => P, PSI_C_d => PSI_C
  USE VAR_PSI_OUT_d, ONLY: PSI_OUT_d => PSI_OUT
  USE VAR_PSI_IN_d, ONLY: PSI_IN_d => PSI_IN

  
  
    IF(MOD(TIME_ITER,IWRITE).NE.0) RETURN

    F1 = F1_d; F2 = F2_d; F3 = F3_d; P = P_d; PSI_C = PSI_C_d 
    PSI_OUT = PSI_OUT_d; PSI_IN = PSI_IN_d 

    OPEN(UNIT=101,FILE='STR3D_V.DAT',FORM='UNFORMATTED',STATUS='UNKNOWN',ERR=705)
    WRITE(101) F1, F2, F3, B1, B2, B3, U1, U2, U3, P, T, PSI_C, PSI_OUT, PSI_IN, TIME
    CLOSE(UNIT=101)
    WRITE(*,*) ' DATA WAS STORED IN DATA FILE.'
    RETURN

  705 WRITE(*,*) ' ERROR IN THE STORING PROCESS!'

  END SUBROUTINE SAVE_DATA


  !-------------------------------------------------------

  PURE FUNCTION WEIGHT_AVG(a1, a2, coeff)
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: a1, a2, coeff
  REAL(8) WEIGHT_AVG
    
    WEIGHT_AVG = coeff * a1 + (1.0_8 - coeff) * a2

  END FUNCTION WEIGHT_AVG

  !-------------------------------------------------

  SUBROUTINE ComputeSourceU1
    USE VAR_d
    IMPLICIT NONE
    
    CALL ComputeForcesU1_d <<<GRD_IJK, TBL_IJK>>>
    CALL ComputeForcesThetaU1_d <<<GRD_IK, TBL_IK>>>

  END SUBROUTINE ComputeSourceU1

  !-------------------------------------------------

  SUBROUTINE ComputeSourceU2
    USE VAR_d
    IMPLICIT NONE

    CALL ComputeForcesU2_d <<<GRD_IJK, TBL_IJK>>>

  END SUBROUTINE ComputeSourceU2

  !-------------------------------------------------

  SUBROUTINE ComputeSourceU3
    USE VAR_d
    IMPLICIT NONE
    
    CALL ComputeForcesU3_d <<<GRD_IJK, TBL_IJK>>>
    CALL ComputeForcesThetaU3_d <<<GRD_IK, TBL_IK>>>

  END SUBROUTINE ComputeSourceU3

  !-------------------------------------------------

  SUBROUTINE ComputeSourceP
    USE VAR_d
    IMPLICIT NONE
    
    CALL ComputeForcesP_d <<<GRD_IJK, TBL_IJK>>>

  END SUBROUTINE ComputeSourceP

  !-------------------------------------------------

  SUBROUTINE ComputeSourcePC
    USE VAR_d
    IMPLICIT NONE
    
    CALL ComputeForcesP_d <<<GRD_IJK, TBL_IJK>>>

  END SUBROUTINE ComputeSourcePC

  !-------------------------------------------------

  SUBROUTINE ComputeSourceT
    USE VAR_d
    IMPLICIT NONE
    
    CALL ComputeForcesT_d <<<GRD_IJK, TBL_IJK>>>

  END SUBROUTINE ComputeSourceT


 !-------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeForcesU1_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: CBOT, CTOP
  REAL(8) :: WOM1, WOM2, WOM 
  REAL(8) :: J2B3m, J2B3p, J3B2m, J3B2p, S_1

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z


    IF ( (I >= (N1_core + 1)) .AND. (I <= (N1-1))  .AND. &
         (J >= 3) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      !------------Archimedes Force-----------------

      CBOT = T(I-1,J,K) * 0.25_8 * XsDifP4(I)
      CTOP = T(I,J,K) * 0.25_8 * XDifP4(I)

      CON(I,J,K) = Ranum*(cosXs2(J)-cosXs2(J+1))*(CBOT+CTOP)*(Xs3(K+1) - Xs3(K)) / r_o
      AP(I, J, K) = 0.0_8

      !------------Coriolis Force--------------------

      WOM1 = ( U3(I,J,K)*VVX1S(I) + U3(I-1,J,K)*VVX2S(I) )*(X3(K)-Xs3(K))
      WOM2 = ( U3(I,J,K+1)*VVX1S(I) + U3(I-1,J,K+1)*VVX2S(I) )*(Xs3(K+1)-X3(K))

      WOM = (WOM1+WOM2)*( (Xs2(J+1) - Xs2(J)) + sinXs2(J) * cosXs2(J) - &
            sinXs2(J+1) * cosXs2(J+1) ) / 3.0_8

      CON(I,J,K) = CON(I,J,K) + WOM

      !Amperes Force (Для вакуумных г.у. рассмотреть вклад токов на границах I =2 и I=N1)------
           
        J2B3m = J2(I, J, K)     *  WEIGHT_AVG( B3(I - 1, J, K),     B3(I, J, K),     Kx1(I) ) 
        J2B3p = J2(I, J, K + 1) *  WEIGHT_AVG( B3(I - 1, J, K + 1), B3(I, J, K + 1), Kx1(I) )
  
        J3B2m = J3(I, J, K)     *  WEIGHT_AVG( B2(I - 1, J, K),     B2(I, J, K),     Kx1(I) )
        J3B2p = J3(I, J + 1, K) *  WEIGHT_AVG( B2(I - 1, J + 1, K), B2(I, J + 1, K), Kx1(I) )
  
        !S_1 = E1jk Jj Bk = J2B3 - J3B2; Eijk - Levi-Civita symbol
        S_1 = 0.5_8 * (J2B3m + J2B3p) - 0.5_8 * (J3B2m + J3B2p) 
        CON(I,J,K) = CON(I,J,K) + S_1 * VCV1(I, J, K) / Pm
     
    ENDIF

  END SUBROUTINE ComputeForcesU1_d
  
 !-------------------------------------------------
  
  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeForcesThetaU1_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, K
  REAL(8) :: CBOT, CTOP
  REAL(8) :: WOM1, WOM2, WOM 
  REAL(8) :: J2B3m, J2B3p, J3B2m, J3B2p, S_1
  
    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= (N1_core + 1)) .AND. (I <= (N1-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN
         
      !------------Archimedes Force-----------------
      !J = 2
      CBOT = T(I-1,2,K) * 0.25_8 * XsDifP4(I)
      CTOP = T(I,2,K) * 0.25_8 * XDifP4(I)

      CON(I,2,K) = Ranum * (cosXs2(2) - cosXs2(3)) * (CBOT + CTOP) * (Xs3(K+1) - Xs3(K)) / r_o
      AP(I,2,K) = 0.0_8
      

      !------------Coriolis Force--------------------
      !J = 2
      WOM1 = ( U3(I,2,K)*VVX1S(I) + U3(I-1,2,K)*VVX2S(I) )*(X3(K)-Xs3(K))
      WOM2 = ( U3(I,2,K+1)*VVX1S(I) + U3(I-1,2,K+1)*VVX2S(I) )*(Xs3(K+1)-X3(K))

      WOM = (WOM1 + WOM2)*( (Xs2(3) - Xs2(2)) + sinXs2(2) * cosXs2(2) - &
            sinXs2(3) * cosXs2(3) ) / 3.0_8

      CON(I,2,K) = CON(I,2,K) + WOM
      

      !------------Amperes Force--------------------
       !J = 2
       J2B3m = J2(I, 2, K)     *  WEIGHT_AVG( B3(I - 1, 2, K),     B3(I, 2, K),     Kx1(I) ) 
       J2B3p = J2(I, 2, K + 1) *  WEIGHT_AVG( B3(I - 1, 2, K + 1), B3(I, 2, K + 1), Kx1(I) )
       J3B2p = J3(I, 3, K)     *  WEIGHT_AVG( B2(I - 1, 3, K),     B2(I, 3, K),     Kx1(I) )
 
       !S_1 = E1jk Jj Bk = J2B3 - J3B2; Eijk - Levi-Civita symbol
       S_1 = 0.5_8 * (J2B3m + J2B3p) - J3B2p 
       CON(I,2,K) = CON(I,2,K) + S_1 * VCV1(I,2,K) / Pm
       

    ENDIF

  END SUBROUTINE ComputeForcesThetaU1_d  

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeForcesU2_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: WW1, WW2 
  REAL(8) :: J1B3m, J1B3p, J3B1m, J3B1p, S_2 


    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (J >= 3) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      !----------------Coriolis force------------------

      WW1 = (U3(I,J,K) * WWY1S(J) + U3(I,J-1,K) * WWY2S(J)) * (X3(K) - Xs3(K))
      WW2 = (U3(I,J,K+1) * WWY1S(J) + U3(I,J-1,K+1)*WWY2S(J)) * (Xs3(K+1) - X3(K))

      CON(I,J,K) = (WW1 + WW2) * (Xs1(I+1)**3 - Xs1(I)**3) / 3.0_8
      AP(I, J, K) = 0.0_8

      !---------------Ampere's force----------------

      J1B3m = J1(I, J, K)     *  WEIGHT_AVG( B3(I, J - 1, K),     B3(I, J, K),     Kx2(J) ) 
      J1B3p = J1(I, J, K + 1) *  WEIGHT_AVG( B3(I, J - 1, K + 1), B3(I, J, K + 1), Kx2(J) )

      J3B1m = J3(I, J, K)     *  WEIGHT_AVG( B1(I, J - 1, K),     B1(I, J, K),     Kx2(J) )
      J3B1p = J3(I + 1, J, K) *  WEIGHT_AVG( B1(I + 1, J - 1, K), B1(I + 1, J, K), Kx2(J) )

      !S_2 = E2jk Jj Bk = J3B1 - J1B3; Eijk - Levi-Civita symbol
      S_2 = 0.5_8 * (J3B1m + J3B1p) - 0.5_8 * (J1B3m + J1B3p) 
      CON(I, J, K) = CON(I, J, K) + S_2 * VCV2(I, J, K) / Pm

    ENDIF

  END SUBROUTINE ComputeForcesU2_d

  !-------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeForcesU3_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: UOM1, UOM2 
  REAL(8) :: VOM1, VOM2
  REAL(8) :: J1B2m, J1B2p, J2B1m, J2B1p, S_3     


    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (J >= 3) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      !----------------Coriolis force------------------

      UOM1 = ( U1(I,J,K) * 0.5_8 * X3cv(K) + U1(I,J,K-1) * 0.5_8 * X3cv(K-1) )* &
             (X1(I)**3 - Xs1(I)**3)
      UOM2 = ( U1(I+1,J,K) * 0.5_8 * X3cv(K) + U1(I+1,J,K-1) * 0.5_8 * X3cv(K-1) )* &
             (Xs1(I+1)**3 - X1(I)**3)
      CON(I,J,K) = -( UOM1 + UOM2 )*( X2cv(J) + sinXs2(J) * cosXs2(J) - &
                   sinXs2(J+1) * cosXs2(J+1) ) / 3.0_8
      AP(I, J, K) = 0.0_8

      VOM1 = (U2(I,J,K) * 0.5_8 * X3cv(K) + U2(I,J,K-1) * 0.5_8 * X3cv(K-1) ) * VZY1S(J) 
      VOM2 = (U2(I,J+1,K) * 0.5_8 * X3cv(K) + U2(I,J+1,K-1) * &
             0.5_8 * X3cv(K-1)) * VZY2S(J)                  

      CON(I,J,K) = CON(I,J,K) - (VOM1 + VOM2) * (Xs1(I+1)**3 - Xs1(I)**3) / 3.0_8

      !---------------Ampere's force----------------

       J1B2m = J1(I, J, K)     *  WEIGHT_AVG( B2(I, J, K - 1),     B2(I, J, K),     Kx3(K) ) 
       J1B2p = J1(I, J + 1, K) *  WEIGHT_AVG( B2(I, J + 1, K - 1), B2(I, J + 1, K), Kx3(K) ) 
 
       J2B1m = J2(I, J, K)     *  WEIGHT_AVG( B1(I, J, K - 1),     B1(I, J, K),     Kx3(K) )
       J2B1p = J2(I + 1, J, K) *  WEIGHT_AVG( B1(I + 1, J, K - 1), B1(I + 1, J, K), Kx3(K) )
 
       !S_3 = E3jk Jj Bk = J1B2 - J2B1; Eijk - Levi-Civita symbol
       S_3 = 0.5_8 * (J1B2m + J1B2p) - 0.5_8 * (J2B1m + J2B1p) 
       CON(I, J, K) = CON(I, J, K) + S_3 * VCV3(I, J, K) / Pm

    ENDIF

  END SUBROUTINE ComputeForcesU3_d
  
    !-------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeForcesThetaU3_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, K
  REAL(8) :: UOM1, UOM2 
  REAL(8) :: VOM1, VOM2
  REAL(8) :: J1B2m, J1B2p, J2B1m, J2B1p, S_3     


    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN
           
      !----------------Coriolis force------------------
      !J = 2
      UOM1 = ( U1(I,2,K) * 0.5_8 * X3cv(K) + U1(I,2,K-1) * 0.5_8 * X3cv(K-1) )* &
             (X1(I)**3 - Xs1(I)**3)
      UOM2 = ( U1(I+1,2,K) * 0.5_8 * X3cv(K) + U1(I+1,2,K-1) * 0.5_8 * X3cv(K-1) )* &
             (Xs1(I+1)**3 - X1(I)**3)
      CON(I,2,K) = -( UOM1 + UOM2 )*( X2cv(2) + sinXs2(2) * cosXs2(2) - &
                     sinXs2(3) * cosXs2(3) ) / 3.0_8
      AP(I,2,K) = 0.0_8

      VOM1 = (U2(I,2,K) * 0.5_8 * X3cv(K) + U2(I,2,K-1) * 0.5_8 * X3cv(K-1) ) * VZY1S(2) 
      VOM2 = (U2(I,3,K) * 0.5_8 * X3cv(K) + U2(I,3,K-1) * &
             0.5_8 * X3cv(K-1)) * VZY2S(2)                  

      CON(I,2,K) = CON(I,2,K) - (VOM1 + VOM2) * (Xs1(I+1)**3 - Xs1(I)**3) / 3.0_8
      
       
      !---------------Ampere's force----------------
       !J = 2  
       
       J1B2p = 0.5_8 * (J1(I, 2, K) + J1(I, 3, K)) * &
               WEIGHT_AVG( B2(I, 3, K - 1), B2(I, 3, K), Kx3(K) ) 
       !Старый вариант:
       !J1B2p = J1(I, 3, K)     *  WEIGHT_AVG( B2(I, 3, K - 1),     B2(I, 3, K),     Kx3(K) ) 
       J2B1m = J2(I, 2, K)     *  WEIGHT_AVG( B1(I, 2, K - 1),     B1(I, 2, K),     Kx3(K) )
       J2B1p = J2(I + 1, 2, K) *  WEIGHT_AVG( B1(I + 1, 2, K - 1), B1(I + 1, 2, K), Kx3(K) )
 
       !S_3 = E3jk Jj Bk = J1B2 - J2B1; Eijk - Levi-Civita symbol
       S_3 = J1B2p - 0.5_8 * (J2B1m + J2B1p) 
       CON(I, 2, K) = CON(I, 2, K) + S_3 * VCV3(I, 2, K) / Pm

    ENDIF

  END SUBROUTINE ComputeForcesThetaU3_d
  !-------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeForcesP_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 1) .AND. (I <= N1)  .AND. &
         (J >= 1) .AND. (J <= N2)  .AND. &
         (K >= 1) .AND. (K <= N3) ) THEN

      CON(I, J, K) = 0.0_8
      AP(I, J, K) = 0.0_8
      PC(I, J, K) = 0.0_8

    ENDIF

  END SUBROUTINE ComputeForcesP_d

  !-------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeForcesT_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 1) .AND. (I <= N1)  .AND. &
         (J >= 1) .AND. (J <= N2)  .AND. &
         (K >= 1) .AND. (K <= N3) ) THEN

      CON(I, J, K) = 0.0_8
      AP(I, J, K) = 0.0_8

    ENDIF

  END SUBROUTINE ComputeForcesT_d

  !-------------------------------------------------

  SUBROUTINE OUTPUT
    USE VAR
    IMPLICIT NONE
    INTEGER :: I,J,K
    REAL(8) :: delF1, delF2, delF3, delF1_1, delF2_1, delF3_1, &
               mF1, mF2, mF3, H1, H2, H3, B1p, B2p, B3p, VOL, Emag, Ekin, &
               U1p, U2p, U3p, FluxB, FluxB_Sin, FluxB_Sout, FluxB_bot, divB, &
               VOL_B1, VOL_B2, VOL_B3, FluxQ, FluxQ_Sin, FluxQ_Sout, Emag_ic

    Emag = 0.0_8
    Ekin = 0.0_8


    DO K = 2, N3-1 
      DO J = 2, N2-1
        DO I = N1_core, N1-1
          VOL_B1 = (X1(I)**3-X1(I-1)**3) * (DCOS(Xs2(J))-DCOS(Xs2(J+1))) * (Xs3(K+1) - Xs3(K)) / 3.0_8
          !IF (I .EQ. 3) VOL_B1 = (X1(3)**3-X1(1)**3) * (DCOS(Xs2(J))-DCOS(Xs2(J+1))) * &
          !                       (Xs3(K+1) - Xs3(K)) / 3.0_8
          !IF (I .EQ. (N1-1)) VOL_B1 = (X1(N1)**3-X1(N1-2)**3) * (DCOS(Xs2(J))-DCOS(Xs2(J+1))) * &
          !                            (Xs3(K+1) - Xs3(K)) / 3.0_8
          Emag = Emag + B1(I,J,K)**2 * VOL_B1
          Ekin = Ekin + U1(I,J,K)**2 * VOL_B1                                         
        ENDDO 
      ENDDO
    ENDDO

    DO K = 2, N3-1 
      DO J = 3, N2-1
        DO I = N1_core, N1-1
          VOL_B2 = (Xs1(I+1)**3-Xs1(I)**3)*(DCOS(X2(J-1)) - DCOS(X2(J))) * (Xs3(K+1) - Xs3(K)) / 3.0_8
          IF (J .EQ. 3) VOL_B2 = (Xs1(I+1)**3-Xs1(I)**3)*(DCOS(X2(1)) - DCOS(X2(3))) * (Xs3(K+1) - Xs3(K)) / 3.0_8
          IF (J .EQ. (N2 -1)) VOL_B2 = (Xs1(I+1)**3-Xs1(I)**3)*(DCOS(X2(N2-2)) - DCOS(X2(N2))) * &
                                       (Xs3(K+1) - Xs3(K)) / 3.0_8
          Emag = Emag + B2(I,J,K)**2 * VOL_B2
          Ekin = Ekin + U2(I,J,K)**2 * VOL_B2    
        ENDDO 
      ENDDO
    ENDDO

    DO K = 2, N3-1 
      DO J = 2, N2-1
        DO I = N1_core, N1-1
          VOL_B3 = (Xs1(I+1)**3-Xs1(I)**3) * (DCOS(Xs2(J))-DCOS(Xs2(J+1))) * Xdif3(K) / 3.0_8
          Emag = Emag + B3(I,J,K)**2 * VOL_B3
          Ekin = Ekin + U3(I,J,K)**2 * VOL_B3    
        ENDDO 
      ENDDO
    ENDDO

    VOL = 2.0_8 * PI * (r_o**3 - r_i**3) / 3.0_8 

    Emag = Emag / ( 2.0_8 * VOL*E* Pm )
    Ekin = Ekin / ( 2.0_8 * VOL)


   Emag_ic = 0.0_8

    DO K = 2, N3-1 
      DO J = 2, N2-1
        DO I = 3, N1_core-1
          VOL_B1 = (X1(I)**3-X1(I-1)**3) * (DCOS(Xs2(J))-DCOS(Xs2(J+1))) * (Xs3(K+1) - Xs3(K)) / 3.0_8
          Emag_ic = Emag_ic + B1(I,J,K)**2 * VOL_B1                                
        ENDDO 
      ENDDO
    ENDDO

    DO K = 2, N3-1 
      DO J = 3, N2-1
        DO I = 2, N1_core-1
          VOL_B2 = (Xs1(I+1)**3-Xs1(I)**3)*(DCOS(X2(J-1)) - DCOS(X2(J))) * (Xs3(K+1) - Xs3(K)) / 3.0_8

          IF (J .EQ. 3) VOL_B2 = (Xs1(I+1)**3-Xs1(I)**3)*(DCOS(X2(1)) - DCOS(X2(3))) * (Xs3(K+1) - Xs3(K)) / 3.0_8
          IF (J .EQ. (N2 -1)) VOL_B2 = (Xs1(I+1)**3-Xs1(I)**3)*(DCOS(X2(N2-2)) - DCOS(X2(N2))) * &
                                       (Xs3(K+1) - Xs3(K)) / 3.0_8
          Emag_ic = Emag_ic + B2(I,J,K)**2 * VOL_B2
        ENDDO 
      ENDDO
    ENDDO

    DO K = 2, N3-1 
      DO J = 2, N2-1
        DO I = 2, N1_core-1
          VOL_B3 = (Xs1(I+1)**3-Xs1(I)**3) * (DCOS(Xs2(J))-DCOS(Xs2(J+1))) * Xdif3(K) / 3.0_8
          Emag_ic = Emag_ic + B3(I,J,K)**2 * VOL_B3 
        ENDDO 
      ENDDO
    ENDDO


    VOL = 2.0_8 * PI * r_i**3 / 3.0_8 

    Emag_ic = Emag_ic / ( 2.0_8 * VOL * E * Pm )


    divB = 0.0_8
    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 2, N1-1
          divB = divB + S1(I+1,J,K)*B1(I+1,J,K) - S1(I,J,K)*B1(I,J,K) + &
                        S2(I,J+1,K)*B2(I,J+1,K) - S2(I,J,K)*B2(I,J,K) + &
                        S3(I,J,K+1)*B3(I,J,K+1) - S3(I,J,K)*B3(I,J,K)    
        ENDDO          
      ENDDO
    ENDDO
   
    FluxB_Sout = 0.0_8
    FluxB_Sin = 0.0_8
    FluxB_bot = 0.0_8
    
    FluxQ_Sout = 0.0_8
    FluxQ_Sin = 0.0_8

    
    DO K = 2, N3-1
      DO J = 2, N2-1
        FluxB_Sout = FluxB_Sout + S1(N1,J,K) * B1(N1,J,K) 
        FluxB_Sin = FluxB_Sin + S1(2,J,K) * B1(2,J,K)
        
        FluxQ_Sout = FluxQ_Sout - S1(N1,J,K) * (T(N1,J,K) - T(N1-1,J,K)) / Xdif1(N1)
        FluxQ_Sin  = FluxQ_Sin  - S1(2,J,K)  * (T(2,J,K)  - T(1,J,K)   ) / Xdif1(2)
      ENDDO
    ENDDO
    
    DO K = 2, N3-1
      DO I = 2, N1-1
        FluxB_bot = FluxB_bot + S2(I,N2,K) * B2(I,N2,K)
      ENDDO
    ENDDO
  
    
    FluxB = FluxB_Sout - FluxB_Sin + FluxB_bot
    
    FluxQ_Sout = 2.0_8 * FluxQ_Sout
    FluxQ_Sin =  2.0_8 * FluxQ_Sin 

    !----------------------------------------------------------------------------------------------

  71  FORMAT(4X,'ITER',8X,'TIME')
  72  FORMAT(7X,'Smax',11X,'SmaxB',10X,'Emag',10X,'Emag_ic',11X,'Ekin',11X,'divB',11X,'FluxB')
  73  FORMAT(5X,'FluxQ_Sin',6X,'FluxQ_Sout')
  74  FORMAT(5X,'E1_Oz(2)',7X,'E1_oz(3)',7X,'E1_oz(4)',7X,'E1_oz(5)',7X,'E1_oz(6)',7X,'E1_oz(7)')
  75  FORMAT(7X,'F1(2)',10X,'F1(3)',10X,'F1(4)',10X,'F1(5)',10X,'F1(6)',10X,'F1(7)')
  76  FORMAT(7X,'F2(2)',10X,'F2(3)',10X,'F2(4)',10X,'F2(5)',10X,'F2(6)',10X,'F2(7)')
  77  FORMAT(7X,'F3(2)',10X,'F3(3)',10X,'F3(4)',10X,'F3(5)',10X,'F3(6)',10X,'F3(7)')


    WRITE(*, *) '==============================================================================='
    WRITE(*, 71)
    WRITE(*,'(I8,1P1E15.5)') TIME_ITER, TIME
    WRITE(*, *) '-------------------------------------------------------------------------------'

    WRITE(*, 72)
    WRITE(*,'(1P67E15.5)') SMAX, SMAX_B, Emag, Emag_ic, Ekin, divB, FluxB
    
    WRITE(*, 73)
    WRITE(*,'(1P2E15.5)') FluxQ_Sin, FluxQ_Sout  

    WRITE(*, 74)
    WRITE(*,'(1P6E15.5)') E1_Oz(2), E1_oz(3), E1_oz(4), E1_oz(5), E1_oz(6), E1_oz(7)

    WRITE(*, 75)
    WRITE(*,'(1P6E15.5)') F1(2,2,5), F1(3,2,5), F1(4,2,5), F1(5,2,5), F1(6,2,5), F1(7,2,5)

    WRITE(*, 76)
    WRITE(*,'(1P6E15.5)') F2(2,3,5), F2(3,3,5), F2(4,3,5), F2(5,3,5), F2(6,3,5), F2(7,3,5)

    WRITE(*, 77)
    WRITE(*,'(1P6E15.5)') F3(2,2,5), F3(3,2,5), F3(4,2,5), F3(5,2,5), F3(6,2,5), F3(7,2,5)


    WRITE(*, *) '===============================================================================' 

    !------------------------------------------File Q.out----------------------------
    WRITE(10, *) '==============================================================================='
    WRITE(10, 71)
    WRITE(10,'(I8,1P1E15.5)') TIME_ITER, TIME
    WRITE(10, *) '-------------------------------------------------------------------------------'

    WRITE(10, 72)
    WRITE(10,'(1P7E15.5)') SMAX, SMAX_B, Emag, Emag_ic, Ekin, divB, FluxB

    WRITE(10, *) '==============================================================================='  


    WRITE(11,'(1P3E15.5)') TIME, SMAX_B, SMAX
    WRITE(12,'(1P4E15.5)') TIME, Emag, Emag_ic, Ekin

    WRITE(13,'(1P7E15.5)') TIME, B1(22, 28, 40), B2(22, 28, 40), B3(22, 28, 40), &
                                 B1(52, 56, 40), B2(52, 56, 40), B3(52, 56, 40)

    WRITE(14,'(1P2E15.5)') TIME, divB
    WRITE(15,'(1P2E15.5)') TIME, FluxB

    WRITE(16,'(1P9E15.5)') TIME, U1(22, 28, 40), U2(22, 28, 40), U3(22, 28, 40), T(22, 28, 40), &
                                 U1(52, 56, 40), U2(52, 56, 40), U3(52, 56, 40), T(52, 56, 40)
                                 
    WRITE(17,'(1P3E15.5)') TIME, FluxQ_Sin, FluxQ_Sout

    WRITE(*, *) '==============================================================================='


  END SUBROUTINE OUTPUT


  !-------------------------------------------------

  SUBROUTINE OUTPUT_END
    USE VAR
    IMPLICIT NONE
    INTEGER :: I, J, K
    REAL(8) :: xd, yd, zd, U1t, U2t, U3t, Ux, Uy, Uz, &
               B1t, B2t, B3t, Bx, By, Bz, s, Emag, Ekin, &
               S1B1p, S1B1, S2B2p, S2B2, S3B3p, S3B3




    OPEN(UNIT=37,FILE='B1.DAT',STATUS='UNKNOWN')
    WRITE(37,*)'VARIABLES ="X","Y","Z","B1","r"'
    WRITE(37,*)'ZONE I=',N1-1,', J=',N2-1,', K=',N3-2,', F=POINT'

    DO K = 2, N3-1 
      DO J = 2, N2
        DO I = 2, N1
          xd=Xs1(I)*DSIN(X2(J))*DCOS(X3(K))
          yd=Xs1(I)*DSIN(X2(J))*DSIN(X3(K))
          zd=Xs1(I)*DCOS(X2(J))

          WRITE(37,'(1P5E15.5)') xd,yd,zd,B1(I,J,K),Xs1(I)
        ENDDO
      ENDDO
    ENDDO

    ENDFILE 37
    CLOSE(37)

    OPEN(UNIT=38,FILE='B2.DAT',STATUS='UNKNOWN')
    WRITE(38,*)'VARIABLES ="X","Y","Z","B2"'
    WRITE(38,*)'ZONE I=',N1,', J=',N2-2,', K=',N3-2,', F=POINT'

    DO K = 2, N3-1 
      DO J = 3, N2
        DO I = 1, N1
          xd=X1(I)*DSIN(Xs2(J))*DCOS(X3(K))
          yd=X1(I)*DSIN(Xs2(J))*DSIN(X3(K))
          zd=X1(I)*DCOS(Xs2(J))

          WRITE(38,'(1P4E15.5)') xd,yd,zd,B2(I,J,K)
        ENDDO
      ENDDO
    ENDDO

    ENDFILE 38
    CLOSE(38)


    OPEN(UNIT=39,FILE='B3.DAT',STATUS='UNKNOWN')
    WRITE(39,*)'VARIABLES ="X","Y","Z","B3"'
    WRITE(39,*)'ZONE I=',N1,', J=',N2-1,', K=',N3-2,', F=POINT'

    DO K = 2, N3-1 
      DO J = 2, N2
        DO I = 1, N1
          xd=X1(I)*DSIN(X2(J))*DCOS(Xs3(K))
          yd=X1(I)*DSIN(X2(J))*DSIN(Xs3(K))
          zd=X1(I)*DCOS(X2(J))

          WRITE(39,'(1P4E15.5)') xd,yd,zd,B3(I,J,K)
        ENDDO
      ENDDO
    ENDDO

    ENDFILE 39
    CLOSE(39)



    OPEN(UNIT=39,FILE='U3.DAT',STATUS='UNKNOWN')
    WRITE(39,*)'VARIABLES ="X","Y","Z","U3"'
    WRITE(39,*)'ZONE I=',N1,', J=',N2-1,', K=',N3-2,', F=POINT'

    DO K = 2, N3-1 
      DO J = 2, N2
        DO I = 1, N1
          xd=X1(I)*DSIN(X2(J))*DCOS(Xs3(K))
          yd=X1(I)*DSIN(X2(J))*DSIN(Xs3(K))
          zd=X1(I)*DCOS(X2(J))

          WRITE(39,'(1P4E15.5)') xd,yd,zd,U3(I,J,K)
        ENDDO
      ENDDO
    ENDDO

    ENDFILE 39
    CLOSE(39)


    OPEN(UNIT=40,FILE='T.DAT',STATUS='UNKNOWN')
    WRITE(40,*)'VARIABLES ="X","Y","Z","T"'
    WRITE(40,*)'ZONE I=',N1-1,', J=',N2-2,', K=',N3-2,', F=POINT'

    DO K = 2, N3-1 
      DO J = 2, N2-1
        DO I = 2, N1
          xd=X1(I)*DSIN(X2(J))*DCOS(X3(K))
          yd=X1(I)*DSIN(X2(J))*DSIN(X3(K))
          zd=X1(I)*DCOS(X2(J))

          WRITE(40,'(1P4E15.5)') xd,yd,zd,T(I,J,K)
        ENDDO
      ENDDO
    ENDDO

    ENDFILE 40
    CLOSE(40)
    
    OPEN(UNIT=41,FILE='ALL.DAT',STATUS='UNKNOWN')
    WRITE(41,*)'VARIABLES = "X", "Y", "Z", "U", "V", "W", "BU", "BV", "BW", "T", "Emag", "Ekin"'
    WRITE(41,*)'ZONE I=',N1-2,', J=',N2-2,', K=',N3-2,', F=POINT'
                
    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 2, N1-1
        
          xd=X1(I)*DSIN(X2(J))*DCOS(X3(K))
          yd=X1(I)*DSIN(X2(J))*DSIN(X3(K))
          zd=X1(I)*DCOS(X2(J))
           
          U1t = 0.5_8 * (U1(I,J,K) + U1(I+1,J,K))
          U2t = 0.5_8 * (U2(I,J,K) + U2(I,J+1,K))
          U3t = 0.5_8 * (U3(I,J,K) + U3(I,J,K+1))
          
          B1t = 0.5_8 * (B1(I,J,K) + B1(I+1,J,K))
          B2t = 0.5_8 * (B2(I,J,K) + B2(I,J+1,K))
          B3t = 0.5_8 * (B3(I,J,K) + B3(I,J,K+1))

          
          IF (J .EQ. 2) THEN
            U2t = U2(I,3,K)
            B2t = B2(I,3,K)
          ENDIF

          IF (J .EQ. (N2-1)) THEN
            U2t = U2(I,N2-1,K)
            B2t = B2(I,N2-1,K)
          ENDIF

          Ux = DSIN(X2(J)) * DCOS(X3(K)) * U1t + DCOS(X2(J)) * DCOS(X3(K)) * U2t - DSIN(X3(K)) * U3t
          Uy = DSIN(X2(J)) * DSIN(X3(K)) * U1t + DCOS(X2(J)) * DSIN(X3(K)) * U2t + DCOS(X3(K)) * U3t
          Uz = DCOS(X2(J)) * U1t - DSIN(X2(J)) * U2t
          
          Bx = DSIN(X2(J)) * DCOS(X3(K)) * B1t + DCOS(X2(J)) * DCOS(X3(K)) * B2t - DSIN(X3(K)) * B3t
          By = DSIN(X2(J)) * DSIN(X3(K)) * B1t + DCOS(X2(J)) * DSIN(X3(K)) * B2t + DCOS(X3(K)) * B3t
          Bz = DCOS(X2(J)) * B1t - DSIN(X2(J)) * B2t

          Emag = (B1t**2 + B2t**2 + B3t**2) / ( 2.0_8 * E * Pm ) 
          Ekin = 0.5_8 * (U1t**2 + U2t**2 + U3t**2)

          WRITE(41,'(1P12E15.5)') xd, yd, zd, Ux, Uy, Uz, Bx, By, Bz, T(I,J,K), Emag, Ekin
        ENDDO           
      ENDDO
    ENDDO

    ENDFILE 41
    CLOSE(41)
    
    
    OPEN(UNIT=42,FILE='S.DAT',STATUS='UNKNOWN')
    WRITE(42,*)'VARIABLES ="X","Y","Z","SmaxB","I","J","K"'
    WRITE(42,*)'ZONE I=',N1-2,', J=',N2-2,', K=',N3-2,', F=POINT'

    DO K = 2, N3-1 
      DO J = 2, N2-1
        DO I = 2, N1-1
          xd = X1(I) * DSIN(X2(J)) * DCOS(X3(K))
          yd = X1(I) * DSIN(X2(J)) * DSIN(X3(K))
          zd = X1(I) * DCOS(X2(J))

          s = DABS(S1(I+1,J,K)*B1(I+1,J,K) - S1(I,J,K)*B1(I,J,K) + &
                   S2(I,J+1,K)*B2(I,J+1,K) - S2(I,J,K)*B2(I,J,K) + &
                   S3(I,J,K+1)*B3(I,J,K+1) - S3(I,J,K)*B3(I,J,K) ) / VCV(I,J,K)
          WRITE(42,'(1P4E15.5,I4,I4,I4)') xd, yd, zd, s, I, J, K 
        ENDDO
      ENDDO
    ENDDO

    ENDFILE 42
    CLOSE(42)

 
  END SUBROUTINE OUTPUT_END


  SUBROUTINE OUTPUT_GRID
    USE VAR_PSI_OUT
    IMPLICIT NONE
    INTEGER :: I

    OPEN(UNIT=43,FILE='X1.DAT',STATUS='UNKNOWN')
    DO I = 1, N1
      WRITE(43,'(1P2E15.5)') DBLE(I), X1(I) 
    ENDDO

    ENDFILE 43
    CLOSE(43)

  END SUBROUTINE OUTPUT_GRID
!------------------------------------------------------------------
END MODULE USER
