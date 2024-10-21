module solver
use cudafor
implicit none
contains

  !------------------------------------------------

  SUBROUTINE AdiSolverF1(NTIMES)
    USE VAR_d
    IMPLICIT NONE
    INTEGER :: MT, NTIMES

    DO MT = 1, NTIMES

      !--------------Theta Dir---------------------
      CALL TempsThetaDirF1_d <<<GRD_IJK, TBL_IJK>>>
      CALL TempsBoundThetaDirF1_d <<<GRD_IJ, TBL_IJ>>>
      CALL AdiThetaDirF1_d <<<GRD_IK, TBL_IK>>>

      !--------------Phi Dir---------------------
      CALL TempsPhiDirF1_d <<<GRD_IJK, TBL_IJK>>>
      CALL AdiPhiDirF1_d <<<GRD_IJ, TBL_IJ>>>

    ENDDO

    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (F1, N1, N1)
 
  END SUBROUTINE AdiSolverF1

  !------------------------------------------------

  SUBROUTINE AdiSolverF2(NTIMES)
    USE VAR_d
    IMPLICIT NONE
    INTEGER :: MT, NTIMES

    DO MT = 1, NTIMES

      !--------------Rho Dir--------------------
      CALL TempsRhoDirF2_d <<<GRD_IJK, TBL_IJK>>>
      CALL TempsBoundRhoDirF2_d <<<GRD_IJ, TBL_IJ>>>
      CALL AdiRhoDirF2_d <<<GRD_JK, TBL_JK>>>

      !--------------Phi Dir--------------------
      CALL TempsPhiDirF2_d <<<GRD_IJK, TBL_IJK>>>
      CALL AdiPhiDirF2_d <<<GRD_IJ, TBL_IJ>>>

    ENDDO

    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (F2, N1, N1)
 
  END SUBROUTINE AdiSolverF2

  !------------------------------------------------

  SUBROUTINE AdiSolverF3(NTIMES)
    USE VAR_d
    IMPLICIT NONE
    INTEGER :: MT, NTIMES

    DO MT = 1, NTIMES

      !----------------Rho Dir-------------------
      CALL TempsRhoDirF3_d <<<GRD_IJK, TBL_IJK>>>
      CALL AdiRhoDirF3_d <<<GRD_JK, TBL_JK>>>

      !----------------Theta Dir-------------------
      CALL TempsThetaDirF3_d <<<GRD_IJK, TBL_IJK>>>
      CALL AdiThetaDirF3_d <<<GRD_IK, TBL_IK>>>

    ENDDO

    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (F3, N1, N1)
 
  END SUBROUTINE AdiSolverF3

  !-----------------------------------------------------------------------

  FUNCTION AdiSolver(F, ISTbeg, ISTend, JST, KST, NTIMES)
    USE VAR_d
    IMPLICIT NONE
    REAL(8), DEVICE, DIMENSION (N1,N2,N3) :: F
    INTEGER :: MT, NTIMES
    INTEGER :: I, J, K
    REAL(8) :: AdiSolver
    INTEGER :: ISTbeg, ISTend, JST, KST

    DO MT = 1, NTIMES

      !---------------- RHO_dir -------------

      CALL TempsRhoDir_d <<<GRD_IJK, TBL_IJK>>> (F, ISTbeg, ISTend, JST, KST)
      CALL TempsBoundRhoDir_d <<<GRD_IJ, TBL_IJ>>> (F, ISTbeg, ISTend, JST, KST)
      CALL AdiRhoDir_d <<<GRD_JK, TBL_JK>>> (F, ISTbeg, ISTend, JST, KST)

      !---------------- THETA_dir -------------

      CALL TempsThetaDir_d <<<GRD_IJK, TBL_IJK>>> (F, ISTbeg, ISTend, JST, KST)
      CALL TempsBoundThetaDir_d <<<GRD_IJ, TBL_IJ>>> (F, ISTbeg, ISTend, JST, KST)
      CALL AdiThetaDir_d <<<GRD_IK, TBL_IK>>> (F, ISTbeg, ISTend, JST, KST)

      !---------------- PHI_dir -------------

      !CALL TempsPhiDir_d <<<GRD_IJK, TBL_IJK>>> (F, ISTbeg, ISTend, JST, KST)
      !CALL TempsBoundPhiDir_d <<<GRD_IJ, TBL_IJ>>> (F, ISTbeg, ISTend, JST, KST)
      !CALL AdiPhiDir_d <<<GRD_IJ, TBL_IJ>>> (F, ISTbeg, ISTend, JST, KST)

    ENDDO

    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (F, N1, N1)

    AdiSolver = NTIMES
 
  END FUNCTION AdiSolver

  !-----------------------------------------------------------------------

  FUNCTION AdiSolverPSI_C(Iend, TOLERANCE)
    USE VAR_d
    IMPLICIT NONE
    INTEGER :: MT
    INTEGER :: I, J, K, Iend
    REAL(8) :: AdiSolverPSI_C, TOLERANCE, RESIDUAL

    MT = 0
    RESIDUAL = 1.0_8

    DO WHILE (RESIDUAL .GE. TOLERANCE)
      MT = MT + 1

      !---------------- RHO_dir -------------

      CALL TempsRhoDirPSI_C_d <<<GRD_IJK, TBL_IJK>>> (Iend)
      CALL TempsBoundRhoDirPSI_C_d <<<GRD_IJ, TBL_IJ>>> (Iend)
      CALL AdiRhoDirPSI_C_d <<<GRD_JK, TBL_JK>>> (Iend)

      !---------------- THETA_dir -------------

      CALL TempsThetaDirPSI_C_d <<<GRD_IJK, TBL_IJK>>> (Iend)
      CALL TempsBoundThetaDirPSI_C_d <<<GRD_IJ, TBL_IJ>>>  (Iend)
      CALL AdiThetaDirPSI_C_d <<<GRD_IK, TBL_IK>>>  (Iend)

      !---------------- PHI_dir -------------

      CALL TempsPhiDirPSI_C_d <<<GRD_IJK, TBL_IJK>>> (Iend)
      CALL AdiPhiDirPSI_C_d <<<GRD_IJ, TBL_IJ>>> (Iend)

      RESIDUAL = ComputeNormaMaxPSI_C_d(Iend)
      !WRITE(*, *) RESIDUAL
    ENDDO

    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (PSI_C, N1, Iend)

    AdiSolverPSI_C = MT
 
  END FUNCTION AdiSolverPSI_C

  !------------------------------------------------
  
    FUNCTION AdiSolverPSI_OUT(TOLERANCE)
    USE VAR_PSI_OUT_d
    IMPLICIT NONE
    INTEGER :: MT
    INTEGER :: I, J, K, ITER
    REAL(8) :: AdiSolverPSI_OUT, TOLERANCE, RESIDUAL

    MT = 0
    RESIDUAL = 1.0_8

    DO WHILE (RESIDUAL .GE. TOLERANCE)

      MT = MT + 1

      !---------------- RHO_dir -------------

      CALL TempsRhoDirPSI_OUT_d <<<GRD_IJK, TBL_IJK>>>
      CALL TempsBoundRhoDirPSI_OUT_d <<<GRD_IJ, TBL_IJ>>>
      CALL AdiRhoDirPSI_OUT_d <<<GRD_JK, TBL_JK>>>

      !---------------- THETA_dir -------------

      CALL TempsThetaDirPSI_OUT_d <<<GRD_IJK, TBL_IJK>>>
      CALL TempsBoundThetaDirPSI_OUT_d <<<GRD_IJ, TBL_IJ>>>
      CALL AdiThetaDirPSI_OUT_d <<<GRD_IK, TBL_IK>>>

      !---------------- PHI_dir -------------

      CALL TempsPhiDirPSI_OUT_d <<<GRD_IJK, TBL_IJK>>>
      CALL AdiPhiDirPSI_OUT_d <<<GRD_IJ, TBL_IJ>>>

      RESIDUAL = ComputeNorma2PSI_OUT_d
      !IF((MT .EQ. 1) .OR. (MOD(MT, 1000) .EQ. 0)) THEN
      !  WRITE(*,*) MT, ': RESIDUAL PSI_OUT = ', RESIDUAL
      !ENDIF
    ENDDO

    WRITE(*,*) MT, ': RESIDUAL PSI_OUT = ', RESIDUAL
    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (PSI_OUT, N1, N1)

    AdiSolverPSI_OUT = MT
 
  END FUNCTION AdiSolverPSI_OUT

  !------------------------------------------------
  
  FUNCTION AdiSolverPSI_IN(TOLERANCE)
    USE VAR_PSI_IN_d
    IMPLICIT NONE
    INTEGER :: MT
    INTEGER :: I, J, K
    REAL(8) :: AdiSolverPSI_IN, TOLERANCE, RESIDUAL

    MT = 0
    RESIDUAL = 1.0_8

    DO WHILE (RESIDUAL .GE. TOLERANCE)

      MT = MT + 1

      !---------------- RHO_dir -------------

      CALL TempsRhoDirPSI_IN_d <<<GRD_IJK, TBL_IJK>>>
      CALL TempsBoundRhoDirPSI_IN_d <<<GRD_IJ, TBL_IJ>>>
      CALL AdiRhoDirPSI_IN_d <<<GRD_JK, TBL_JK>>>

      !---------------- THETA_dir -------------

      CALL TempsThetaDirPSI_IN_d <<<GRD_IJK, TBL_IJK>>>
      CALL TempsBoundThetaDirPSI_IN_d <<<GRD_IJ, TBL_IJ>>>
      CALL AdiThetaDirPSI_IN_d <<<GRD_IK, TBL_IK>>>

      !---------------- PHI_dir -------------

      CALL TempsPhiDirPSI_IN_d <<<GRD_IJK, TBL_IJK>>>
      CALL TempsBoundPhiDirPSI_IN_d <<<GRD_IJ, TBL_IJ>>>
      CALL AdiPhiDirPSI_IN_d <<<GRD_IJ, TBL_IJ>>>

      RESIDUAL = ComputeNorma2PSI_IN_d
      !IF((MT .EQ. 1) .OR. (MOD(MT, 1000) .EQ. 0)) THEN
      !  WRITE(*,*) MT, ': RESIDUAL PSI_IN = ', RESIDUAL
      !ENDIF

    ENDDO

    WRITE(*,*) MT, ': RESIDUAL PSI_IN = ', RESIDUAL
    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (PSI_IN, N1, N1)

    AdiSolverPSI_IN = MT
 
  END FUNCTION AdiSolverPSI_IN
  
  !================================================

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsBoundThetaDirF1_d
    USE VAR_d
    IMPLICIT NONE
    INTEGER :: I, J
    INTEGER, PARAMETER :: IST = 3, JST = 2, KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= IST) .AND. (I <= N1) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) ) THEN

          TEMPS(I,J,2) = CON(I,J,2) &
                       + AKP(I,J,2)*F1(I,J,3) &
                       + AKM(I,J,2)*F1(I,J,N3-1)

          TEMPS(I,J,N3-1) = CON(I,J,N3-1) &
                          + AKP(I,J,N3-1)*F1(I,J,2) &
                          + AKM(I,J,N3-1)*F1(I,J,N3-2)
    ENDIF

  END SUBROUTINE TempsBoundThetaDirF1_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsThetaDirF1_d
    USE VAR_d
    IMPLICIT NONE
    INTEGER :: I, J, K
    INTEGER, PARAMETER :: IST = 3, JST = 2, KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= IST) .AND. (I <= N1) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) .AND. &
         (K >= 3)   .AND. (K <= (N3-2)) ) THEN
            TEMPS(I,J,K) = CON(I,J,K) &
                         + AKP(I,J,K)*F1(I,J,K+1) &
                         + AKM(I,J,K)*F1(I,J,K-1)
    ENDIF

  END SUBROUTINE TempsThetaDirF1_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE AdiThetaDirF1_d
    USE VAR_d
    IMPLICIT NONE
    INTEGER :: I, J, K
    INTEGER, PARAMETER :: IST = 3, JST = 2, KST = 2
    REAL(8), DIMENSION (Nmax) :: PT, QT
    REAL(8) :: DENOM


    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= IST) .AND. (I <= N1) .AND. &
         (K >= KST) .AND. (K <= (N3-1)) ) THEN

      PT(JST-1) = 0.0_8
      QT(JST-1) = F1(I,JST-1,K)

      DO J = JST, N2-1
        DENOM = AP(I,J,K) - PT(J-1) * AJM(I,J,K)
        PT(J) = AJP(I,J,K) / DENOM
        QT(J) = (TEMPS(I,J,K) + AJM(I,J,K)*QT(J-1)) / DENOM
      ENDDO

      DO J = N2-1, JST, -1
        F1(I,J,K) = F1(I,J+1,K) * PT(J) + QT(J)
      ENDDO

    ENDIF

  END SUBROUTINE AdiThetaDirF1_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsPhiDirF1_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: IST = 3, JST = 2, KST = 2
  
    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= IST) .AND. (I <= N1) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) .AND. &
         (K >= KST) .AND. (K <= (N3-1)) ) THEN
            TEMPS(I,J,K) = CON(I,J,K) &
                         + AJP(I,J,K)*F1(I,J+1,K) &
                         + AJM(I,J,K)*F1(I,J-1,K)
   ENDIF

  END SUBROUTINE TempsPhiDirF1_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE AdiPhiDirF1_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: IST = 3, JST = 2, KST = 2
  REAL(8), DIMENSION (Nmax) :: PT, QT, RT, ALFA, BETA 
  REAL(8) :: DENOM

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= IST) .AND. (I <= N1) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) ) THEN

      PT(2) = AKP(I,J,2)/AP(I,J,2)
      RT(2) = AKM(I,J,2)/AP(I,J,2)
      QT(2) = TEMPS(I,J,2)/AP(I,J,2)
   
      DO K = 3, N3-2
        DENOM = AP(I,J,K)-PT(K-1)*AKM(I,J,K)
        PT(K) = AKP(I,J,K)/DENOM
        RT(K) = AKM(I,J,K)*RT(K-1)/DENOM
        QT(K) = (TEMPS(I,J,K)+AKM(I,J,K)*QT(K-1))/DENOM
      ENDDO
      
      ALFA(N3-2) = PT(N3-2)+RT(N3-2)
      BETA(N3-2) = QT(N3-2)

      DO K = N3-3, 2, -1
        ALFA(K) = PT(K)*ALFA(K+1)+RT(K)
        BETA(K) = PT(K)*BETA(K+1)+QT(K)
      ENDDO

      F1(I,J,N3-1) = (AKP(I,J,N3-1)*BETA(2)+AKM(I,J,N3-1)*BETA(N3-2)+TEMPS(I,J,N3-1))/ &
                     (AP(I,J,N3-1)-AKP(I,J,N3-1)*ALFA(2)-AKM(I,J,N3-1)*ALFA(N3-2))         

      DO K = 2, N3-2
        F1(I,J,K) = ALFA(K)*F1(I,J,N3-1)+BETA(K)
      ENDDO

    ENDIF

  END SUBROUTINE AdiPhiDirF1_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsBoundRhoDirF2_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J
  INTEGER, PARAMETER :: IST = 2, JST = 3, KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= IST) .AND. (I <= (N1-1)) .AND. &
         (J >= JST) .AND. (J <= N2) ) THEN

      TEMPS(I,J,2) = CON(I,J,2) &
                   + AKP(I,J,2)*F2(I,J,3) &
                   + AKM(I,J,2)*F2(I,J,N3-1)

      TEMPS(I,J,N3-1) = CON(I,J,N3-1) &
                      + AKP(I,J,N3-1)*F2(I,J,2) &
                      + AKM(I,J,N3-1)*F2(I,J,N3-2)
    ENDIF


  END SUBROUTINE TempsBoundRhoDirF2_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsRhoDirF2_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: IST = 2, JST = 3, KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= IST) .AND. (I <= (N1-1)) .AND. &
         (J >= JST) .AND. (J <= N2) .AND. &
         (K >= 3)   .AND. (K <= (N3-2)) ) THEN
      TEMPS(I,J,K) = CON(I,J,K) &
                   + AKP(I,J,K)*F2(I,J,K+1) &
                   + AKM(I,J,K)*F2(I,J,K-1)
    ENDIF

  END SUBROUTINE TempsRhoDirF2_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsPhiDirF2_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: IST = 2, JST = 3, KST = 2
  
    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= IST) .AND. (I <= (N1-1)) .AND. &
         (J >= JST) .AND. (J <= N2) .AND. &
         (K >= KST) .AND. (K <= (N3-1)) ) THEN
      TEMPS(I,J,K) = CON(I,J,K) &
                   + AIP(I,J,K)*F2(I+1,J,K) &
                   + AIM(I,J,K)*F2(I-1,J,K)
    ENDIF

  END SUBROUTINE TempsPhiDirF2_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE AdiRhoDirF2_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: IST = 2, JST = 3, KST = 2
  REAL(8), DIMENSION (Nmax) :: PT, QT
  REAL(8) :: DENOM

    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= JST) .AND. (J <= N2) .AND. &
         (K >= KST) .AND. (K <= (N3-1)) ) THEN

      PT(IST-1) = 0.0_8
      QT(IST-1) = F2(IST-1,J,K)

      DO I = IST, N1-1
        DENOM = AP(I,J,K)-PT(I-1)*AIM(I,J,K)
        PT(I) = AIP(I,J,K)/ DENOM
        QT(I) = (TEMPS(I,J,K)+AIM(I,J,K)*QT(I-1))/ DENOM
      ENDDO

      DO I = N1-1, IST, -1
        F2(I,J,K) = F2(I+1,J,K)*PT(I)+QT(I)
      ENDDO

    ENDIF

  END SUBROUTINE AdiRhoDirF2_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE AdiPhiDirF2_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: IST = 2, JST = 3, KST = 2
  REAL(8), DIMENSION (Nmax) :: PT, QT, RT, ALFA, BETA 
  REAL(8) :: DENOM

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= IST) .AND. (I <= (N1-1)) .AND. &
         (J >= JST) .AND. (J <= N2) ) THEN

      PT(2) = AKP(I,J,2)/AP(I,J,2)
      RT(2) = AKM(I,J,2)/AP(I,J,2)
      QT(2) = TEMPS(I,J,2)/AP(I,J,2)
   
      DO K = 3, N3-2
        DENOM = AP(I,J,K)-PT(K-1)*AKM(I,J,K)
        PT(K) = AKP(I,J,K)/DENOM
        RT(K) = AKM(I,J,K)*RT(K-1)/DENOM
        QT(K) = (TEMPS(I,J,K)+AKM(I,J,K)*QT(K-1))/DENOM
      ENDDO
      
      ALFA(N3-2) = PT(N3-2)+RT(N3-2)
      BETA(N3-2) = QT(N3-2)

      DO K = N3-3, 2, -1
        ALFA(K) = PT(K)*ALFA(K+1)+RT(K)
        BETA(K) = PT(K)*BETA(K+1)+QT(K)
      ENDDO

      F2(I,J,N3-1) = (AKP(I,J,N3-1)*BETA(2)+AKM(I,J,N3-1)*BETA(N3-2)+TEMPS(I,J,N3-1))/ &
                     (AP(I,J,N3-1)-AKP(I,J,N3-1)*ALFA(2)-AKM(I,J,N3-1)*ALFA(N3-2))         

      DO K = 2, N3-2
        F2(I,J,K) = ALFA(K)*F2(I,J,N3-1)+BETA(K)
      ENDDO

    ENDIF

  END SUBROUTINE AdiPhiDirF2_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsRhoDirF3_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: IST = 2, JST = 2, KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= IST) .AND. (I <= (N1-1)) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) .AND. &
         (K >= KST) .AND. (K <= (N3-1)) ) THEN
      TEMPS(I,J,K) = CON(I,J,K) &
                   + AJP(I,J,K)*F3(I,J+1,K) &
                   + AJM(I,J,K)*F3(I,J-1,K)
    ENDIF

  END SUBROUTINE TempsRhoDirF3_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE AdiRhoDirF3_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: IST = 2, JST = 2, KST = 2
  REAL(8), DIMENSION (Nmax) :: PT, QT
  REAL(8) :: DENOM


    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= JST) .AND. (J <= (N2-1)) .AND. &
         (K >= KST) .AND. (K <= (N3-1)) ) THEN

      PT(IST-1) = 0.0_8
      QT(IST-1) = F3(IST-1,J,K)

      DO I = IST, N1-1
        DENOM = AP(I,J,K)-PT(I-1)*AIM(I,J,K)
        PT(I) = AIP(I,J,K)/ DENOM
        QT(I) = (TEMPS(I,J,K)+AIM(I,J,K)*QT(I-1))/ DENOM
      ENDDO

      DO I = N1-1, IST, -1
        F3(I,J,K) = F3(I+1,J,K)*PT(I)+QT(I)
      ENDDO

    ENDIF

  END SUBROUTINE AdiRhoDirF3_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsThetaDirF3_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: IST = 2, JST = 2, KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= IST) .AND. (I <= (N1-1)) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) .AND. &
         (K >= KST) .AND. (K <= (N3-1)) ) THEN
            TEMPS(I,J,K) = CON(I,J,K) &
                         + AIP(I,J,K)*F3(I+1,J,K) &
                         + AIM(I,J,K)*F3(I-1,J,K)
    ENDIF

  END SUBROUTINE TempsThetaDirF3_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE AdiThetaDirF3_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: IST = 2, JST = 2, KST = 2
  REAL(8), DIMENSION (Nmax) :: PT, QT
  REAL(8) :: DENOM


    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= IST) .AND. (I <= (N1-1)) .AND. &
         (K >= KST) .AND. (K <= (N3-1)) ) THEN

      PT(JST-1) = 0.0_8
      QT(JST-1) = F3(I,JST-1,K)

      DO J = JST, N2-1
        DENOM = AP(I,J,K) - PT(J-1) * AJM(I,J,K)
        PT(J) = AJP(I,J,K) / DENOM
        QT(J) = (TEMPS(I,J,K) + AJM(I,J,K)*QT(J-1)) / DENOM
      ENDDO

      DO J = N2-1, JST, -1
        F3(I,J,K) = F3(I,J+1,K) * PT(J) + QT(J)
      ENDDO

    ENDIF

  END SUBROUTINE AdiThetaDirF3_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsRhoDir_d(F, ISTbeg, ISTend, JST, KST)
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8), DIMENSION (N1,N2,N3) :: F
  INTEGER, VALUE :: ISTbeg, ISTend, JST, KST


    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= ISTbeg) .AND. (I <= ISTend)  .AND. &
         (J >= JST)    .AND. (J <= (N2-1))  .AND. &
         (K >= 3)      .AND. (K <= (N3-2)) ) THEN
      TEMPS(I,J,K) = CON(I,J,K) &
                   + AJP(I,J,K) * F(I,J+1,K) &
                   + AJM(I,J,K) * F(I,J-1,K) &
                   + AKP(I,J,K) * F(I,J,K+1) &
                   + AKM(I,J,K) * F(I,J,K-1)
    ENDIF

  END SUBROUTINE TempsRhoDir_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsBoundRhoDir_d(F, ISTbeg, ISTend, JST, KST)
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J
  REAL(8), DIMENSION (N1,N2,N3) :: F
  INTEGER, VALUE :: ISTbeg, ISTend, JST, KST

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) ) THEN

      TEMPS(I,J,2) = CON(I,J,2) &
                   + AJP(I,J,2) * F(I,J+1,2) &
                   + AJM(I,J,2) * F(I,J-1,2) &
                   + AKP(I,J,2) * F(I,J,3) &
                   + AKM(I,J,2) * F(I,J,N3-1)

      TEMPS(I,J,N3-1) = CON(I,J,N3-1) &
                      + AJP(I,J,N3-1) * F(I,J+1,N3-1) &
                      + AJM(I,J,N3-1) * F(I,J-1,N3-1) &
                      + AKP(I,J,N3-1) * F(I,J,2) &
                      + AKM(I,J,N3-1) * F(I,J,N3-2)
    ENDIF

  END SUBROUTINE TempsBoundRhoDir_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsBoundThetaDir_d(F, ISTbeg, ISTend, JST, KST)
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J
  REAL(8), DIMENSION (N1,N2,N3) :: F
  INTEGER, VALUE :: ISTbeg, ISTend, JST, KST

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) ) THEN

      TEMPS(I,J,2) = CON(I,J,2) &
                   + AIP(I,J,2) * F(I+1,J,2) &
                   + AIM(I,J,2) * F(I-1,J,2) &
                   + AKP(I,J,2) * F(I,J,3) &
                   + AKM(I,J,2) * F(I,J,N3-1)

      TEMPS(I,J,N3-1) = CON(I,J,N3-1) &
                      + AIP(I,J,N3-1) * F(I+1,J,N3-1) &
                      + AIM(I,J,N3-1) * F(I-1,J,N3-1) &
                      + AKP(I,J,N3-1) * F(I,J,2) &
                      + AKM(I,J,N3-1) * F(I,J,N3-2)
    ENDIF

  END SUBROUTINE TempsBoundThetaDir_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsThetaDir_d(F, ISTbeg, ISTend, JST, KST)
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8), DIMENSION (N1,N2,N3) :: F
  INTEGER, VALUE :: ISTbeg, ISTend, JST, KST

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST)    .AND. (J <= (N2-1)) .AND. &
         (K >= 3)      .AND. (K <= (N3-2)) ) THEN

      TEMPS(I,J,K) = CON(I,J,K) &
                   + AIP(I,J,K) * F(I+1,J,K) &
                   + AIM(I,J,K) * F(I-1,J,K) &
                   + AKP(I,J,K) * F(I,J,K+1) &
                   + AKM(I,J,K) * F(I,J,K-1)
    ENDIF

  END SUBROUTINE TempsThetaDir_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsBoundPhiDir_d(F, ISTbeg, ISTend, JST, KST)
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J
  REAL(8), DIMENSION (N1,N2,N3) :: F
  INTEGER, VALUE :: ISTbeg, ISTend, JST, KST

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) ) THEN

      TEMPS(I,J,2) = CON(I,J,2) &
                   + AIP(I,J,2) * F(I+1,J,2) &
                   + AIM(I,J,2) * F(I-1,J,2) &
                   + AJP(I,J,2) * F(I,J+1,2) &
                   + AJM(I,J,2) * F(I,J-1,2)


      TEMPS(I,J,N3-1) = CON(I,J,N3-1) &
                    + AIP(I,J,N3-1) * F(I+1,J,N3-1) &
                    + AIM(I,J,N3-1) * F(I-1,J,N3-1) &
                    + AJP(I,J,N3-1) * F(I,J+1,N3-1) &
                    + AJM(I,J,N3-1) * F(I,J-1,N3-1)
    ENDIF

  END SUBROUTINE TempsBoundPhiDir_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsPhiDir_d(F, ISTbeg, ISTend, JST, KST)
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8), DIMENSION (N1,N2,N3) :: F
  INTEGER, VALUE :: ISTbeg, ISTend, JST, KST
  
    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) .AND. &
         (K >= 3) .AND. (K <= (N3-2)) ) THEN
      TEMPS(I,J,K) = CON(I,J,K) &
                   + AIP(I,J,K) * F(I+1,J,K) &
                   + AIM(I,J,K) * F(I-1,J,K) &
                   + AJP(I,J,K) * F(I,J+1,K) &
                   + AJM(I,J,K) * F(I,J-1,K)
   ENDIF

  END SUBROUTINE TempsPhiDir_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE AdiRhoDir_d(F, ISTbeg, ISTend, JST, KST)
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8), DIMENSION (N1,N2,N3) :: F
  INTEGER, VALUE :: ISTbeg, ISTend, JST, KST
  REAL(8), DIMENSION (Nmax) :: PT, QT
  REAL(8) :: DENOM


    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= JST) .AND. (J <= (N2-1)) .AND. &
         (K >= KST) .AND. (K <= (N3-1)) ) THEN

      PT(ISTbeg-1) = 0.0_8
      QT(ISTbeg-1) = F(ISTbeg-1,J,K)

      DO I = ISTbeg, ISTend
        DENOM = AP(I,J,K)-PT(I-1)*AIM(I,J,K)
        PT(I) = AIP(I,J,K) / DENOM
        QT(I) = (TEMPS(I,J,K) + AIM(I,J,K) * QT(I-1))/ DENOM
      ENDDO

      DO I = ISTend, ISTbeg, -1
        F(I,J,K) = F(I+1,J,K) * PT(I) + QT(I)
      ENDDO

    ENDIF

  END SUBROUTINE AdiRhoDir_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE AdiThetaDir_d(F, ISTbeg, ISTend, JST, KST)
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8), DIMENSION (N1,N2,N3) :: F
  INTEGER, VALUE :: ISTbeg, ISTend, JST, KST
  REAL(8), DIMENSION (Nmax) :: PT, QT
  REAL(8) :: DENOM

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (K >= KST) .AND. (K <= (N3-1)) ) THEN

      PT(JST-1) = 0.0_8
      QT(JST-1) = F(I,JST-1,K)

      DO J = JST, N2-1
        DENOM = AP(I,J,K) - PT(J-1) * AJM(I,J,K)
        PT(J) = AJP(I,J,K) / DENOM
        QT(J) = (TEMPS(I,J,K) + AJM(I,J,K) * QT(J-1)) / DENOM
      ENDDO

      DO J = N2-1, JST, -1
        F(I,J,K) = F(I,J+1,K) * PT(J) + QT(J)
      ENDDO

    ENDIF

  END SUBROUTINE AdiThetaDir_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE AdiPhiDir_d(F, ISTbeg, ISTend, JST, KST)
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8), DIMENSION (N1,N2,N3) :: F
  INTEGER, VALUE :: ISTbeg, ISTend, JST, KST
  REAL(8), DIMENSION (Nmax) :: PT, QT, RT, ALFA, BETA 
  REAL(8) :: DENOM

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) ) THEN

      PT(2) = AKP(I,J,2)/AP(I,J,2)
      RT(2) = AKM(I,J,2)/AP(I,J,2)
      QT(2) = TEMPS(I,J,2)/AP(I,J,2)
   
      DO K = 3, N3-2
        DENOM = AP(I,J,K)-PT(K-1)*AKM(I,J,K)
        PT(K) = AKP(I,J,K)/DENOM
        RT(K) = AKM(I,J,K)*RT(K-1)/DENOM
        QT(K) = (TEMPS(I,J,K)+AKM(I,J,K)*QT(K-1))/DENOM
      ENDDO
      
      ALFA(N3-2) = PT(N3-2)+RT(N3-2)
      BETA(N3-2) = QT(N3-2)

      DO K = N3-3, 2, -1
        ALFA(K) = PT(K)*ALFA(K+1)+RT(K)
        BETA(K) = PT(K)*BETA(K+1)+QT(K)
      ENDDO

      F(I,J,N3-1) = (AKP(I,J,N3-1)*BETA(2)+AKM(I,J,N3-1)*BETA(N3-2)+TEMPS(I,J,N3-1))/ &
                      (AP(I,J,N3-1)-AKP(I,J,N3-1)*ALFA(2)-AKM(I,J,N3-1)*ALFA(N3-2))         

      DO K = 2, N3-2
        F(I,J,K) = ALFA(K) * F(I,J,N3-1) + BETA(K)
      ENDDO

    ENDIF

  END SUBROUTINE AdiPhiDir_d

  !------------------------------------------------ 

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsRhoDirPSI_C_d(Iend)
  USE VAR_d, ONLY: TEMPS, CON, AJP, AJM, AKP, AKM, PSI_C, N2, N3    
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, VALUE :: Iend
  INTEGER :: ISTbeg, ISTend, JST, KST
  
  ISTbeg = 2; ISTend = Iend-1; JST = 2; KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= ISTbeg) .AND. (I <= ISTend)  .AND. &
         (J >= JST)    .AND. (J <= (N2-1))  .AND. &
         (K >= 3)      .AND. (K <= (N3-2)) ) THEN
      TEMPS(I,J,K) = CON(I,J,K) &
                   + AJP(I,J,K) * PSI_C(I,J+1,K) &
                   + AJM(I,J,K) * PSI_C(I,J-1,K) &
                   + AKP(I,J,K) * PSI_C(I,J,K+1) &
                   + AKM(I,J,K) * PSI_C(I,J,K-1)
    ENDIF

  END SUBROUTINE TempsRhoDirPSI_C_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsRhoDirPSI_OUT_d
  USE VAR_PSI_OUT_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= ISTbeg) .AND. (I <= ISTend)  .AND. &
         (J >= JST)    .AND. (J <= (N2-1))  .AND. &
         (K >= 3)      .AND. (K <= (N3-2)) ) THEN
      TEMPS(I,J,K) = CON(I,J,K) &
                   + AJP(I,J,K) * PSI_OUT(I,J+1,K) &
                   + AJM(I,J,K) * PSI_OUT(I,J-1,K) &
                   + AKP(I,J,K) * PSI_OUT(I,J,K+1) &
                   + AKM(I,J,K) * PSI_OUT(I,J,K-1)
    ENDIF

  END SUBROUTINE TempsRhoDirPSI_OUT_d


  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsBoundRhoDirPSI_C_d(Iend)
  USE VAR_d, ONLY: TEMPS, CON, AJP, AJM, AKP, AKM, PSI_C, N2, N3
  IMPLICIT NONE
  INTEGER :: I, J
  INTEGER, VALUE :: Iend
  INTEGER :: ISTbeg, ISTend, JST, KST

    ISTbeg = 2; ISTend = Iend-1; JST = 2; KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) ) THEN

      TEMPS(I,J,2) = CON(I,J,2) &
                   + AJP(I,J,2) * PSI_C(I,J+1,2) &
                   + AJM(I,J,2) * PSI_C(I,J-1,2) &
                   + AKP(I,J,2) * PSI_C(I,J,3) &
                   + AKM(I,J,2) * PSI_C(I,J,N3-1)

      TEMPS(I,J,N3-1) = CON(I,J,N3-1) &
                      + AJP(I,J,N3-1) * PSI_C(I,J+1,N3-1) &
                      + AJM(I,J,N3-1) * PSI_C(I,J-1,N3-1) &
                      + AKP(I,J,N3-1) * PSI_C(I,J,2) &
                      + AKM(I,J,N3-1) * PSI_C(I,J,N3-2)
    ENDIF

  END SUBROUTINE TempsBoundRhoDirPSI_C_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsBoundRhoDirPSI_OUT_d
  USE VAR_PSI_OUT_d
  IMPLICIT NONE
  INTEGER :: I, J
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) ) THEN

      TEMPS(I,J,2) = CON(I,J,2) &
                   + AJP(I,J,2) * PSI_OUT(I,J+1,2) &
                   + AJM(I,J,2) * PSI_OUT(I,J-1,2) &
                   + AKP(I,J,2) * PSI_OUT(I,J,3) &
                   + AKM(I,J,2) * PSI_OUT(I,J,N3-1)

      TEMPS(I,J,N3-1) = CON(I,J,N3-1) &
                      + AJP(I,J,N3-1) * PSI_OUT(I,J+1,N3-1) &
                      + AJM(I,J,N3-1) * PSI_OUT(I,J-1,N3-1) &
                      + AKP(I,J,N3-1) * PSI_OUT(I,J,2) &
                      + AKM(I,J,N3-1) * PSI_OUT(I,J,N3-2)
    ENDIF

  END SUBROUTINE TempsBoundRhoDirPSI_OUT_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsBoundThetaDirPSI_C_d(Iend)
  USE VAR_d, ONLY: TEMPS, CON, AIP, AIM, AKP, AKM, PSI_C, N2, N3
  IMPLICIT NONE  
  INTEGER :: I, J
  INTEGER, VALUE :: Iend
  INTEGER :: ISTbeg, ISTend, JST, KST

   ISTbeg = 2; ISTend = Iend-1; JST = 2; KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) ) THEN

      TEMPS(I,J,2) = CON(I,J,2) &
                   + AIP(I,J,2) * PSI_C(I+1,J,2) &
                   + AIM(I,J,2) * PSI_C(I-1,J,2) &
                   + AKP(I,J,2) * PSI_C(I,J,3) &
                   + AKM(I,J,2) * PSI_C(I,J,N3-1)

      TEMPS(I,J,N3-1) = CON(I,J,N3-1) &
                      + AIP(I,J,N3-1) * PSI_C(I+1,J,N3-1) &
                      + AIM(I,J,N3-1) * PSI_C(I-1,J,N3-1) &
                      + AKP(I,J,N3-1) * PSI_C(I,J,2) &
                      + AKM(I,J,N3-1) * PSI_C(I,J,N3-2)
    ENDIF

  END SUBROUTINE TempsBoundThetaDirPSI_C_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsBoundThetaDirPSI_OUT_d
  USE VAR_PSI_OUT_d
  IMPLICIT NONE
  INTEGER :: I, J
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) ) THEN

      TEMPS(I,J,2) = CON(I,J,2) &
                   + AIP(I,J,2) * PSI_OUT(I+1,J,2) &
                   + AIM(I,J,2) * PSI_OUT(I-1,J,2) &
                   + AKP(I,J,2) * PSI_OUT(I,J,3) &
                   + AKM(I,J,2) * PSI_OUT(I,J,N3-1)

      TEMPS(I,J,N3-1) = CON(I,J,N3-1) &
                      + AIP(I,J,N3-1) * PSI_OUT(I+1,J,N3-1) &
                      + AIM(I,J,N3-1) * PSI_OUT(I-1,J,N3-1) &
                      + AKP(I,J,N3-1) * PSI_OUT(I,J,2) &
                      + AKM(I,J,N3-1) * PSI_OUT(I,J,N3-2)
    ENDIF

  END SUBROUTINE TempsBoundThetaDirPSI_OUT_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsThetaDirPSI_C_d(Iend)
  USE VAR_d, ONLY: TEMPS, CON, AIP, AIM, AKP, AKM, PSI_C, N2, N3
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, VALUE :: Iend
  INTEGER :: ISTbeg, ISTend, JST, KST

    ISTbeg = 2; ISTend = Iend-1; JST = 2; KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST)    .AND. (J <= (N2-1)) .AND. &
         (K >= 3)      .AND. (K <= (N3-2)) ) THEN

      TEMPS(I,J,K) = CON(I,J,K) &
                   + AIP(I,J,K) * PSI_C(I+1,J,K) &
                   + AIM(I,J,K) * PSI_C(I-1,J,K) &
                   + AKP(I,J,K) * PSI_C(I,J,K+1) &
                   + AKM(I,J,K) * PSI_C(I,J,K-1)
    ENDIF

  END SUBROUTINE TempsThetaDirPSI_C_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsThetaDirPSI_OUT_d
  USE VAR_PSI_OUT_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST)    .AND. (J <= (N2-1)) .AND. &
         (K >= 3)      .AND. (K <= (N3-2)) ) THEN

      TEMPS(I,J,K) = CON(I,J,K) &
                   + AIP(I,J,K) * PSI_OUT(I+1,J,K) &
                   + AIM(I,J,K) * PSI_OUT(I-1,J,K) &
                   + AKP(I,J,K) * PSI_OUT(I,J,K+1) &
                   + AKM(I,J,K) * PSI_OUT(I,J,K-1)
    ENDIF

  END SUBROUTINE TempsThetaDirPSI_OUT_d


  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsPhiDirPSI_C_d(Iend)
  USE VAR_d, ONLY: TEMPS, CON, AIP, AIM, AJP, AJM, PSI_C, N2, N3
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, VALUE :: Iend
  INTEGER :: ISTbeg, ISTend, JST, KST

    ISTbeg = 2; ISTend = Iend-1; JST = 2; KST = 2
  
    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) .AND. &
         (K >= KST) .AND. (K <= (N3-1)) ) THEN
      TEMPS(I,J,K) = CON(I,J,K) &
                   + AIP(I,J,K) * PSI_C(I+1,J,K) &
                   + AIM(I,J,K) * PSI_C(I-1,J,K) &
                   + AJP(I,J,K) * PSI_C(I,J+1,K) &
                   + AJM(I,J,K) * PSI_C(I,J-1,K)
   ENDIF

  END SUBROUTINE TempsPhiDirPSI_C_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsPhiDirPSI_OUT_d
  USE VAR_PSI_OUT_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2
  
    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) .AND. &
         (K >= KST) .AND. (K <= (N3-1)) ) THEN
      TEMPS(I,J,K) = CON(I,J,K) &
                   + AIP(I,J,K) * PSI_OUT(I+1,J,K) &
                   + AIM(I,J,K) * PSI_OUT(I-1,J,K) &
                   + AJP(I,J,K) * PSI_OUT(I,J+1,K) &
                   + AJM(I,J,K) * PSI_OUT(I,J-1,K)
   ENDIF

  END SUBROUTINE TempsPhiDirPSI_OUT_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE AdiRhoDirPSI_C_d(Iend)
  USE VAR_d, ONLY: TEMPS, CON, AP, AIM, AIP, PSI_C, Nmax, N2, N3
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, VALUE :: Iend
  INTEGER :: ISTbeg, ISTend, JST, KST
  REAL(8), DIMENSION (Nmax) :: PT, QT
  REAL(8) :: DENOM

  ISTbeg = 2; ISTend = Iend-1; JST = 2; KST = 2


    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= JST) .AND. (J <= (N2-1)) .AND. &
         (K >= KST) .AND. (K <= (N3-1)) ) THEN

      PT(ISTbeg-1) = 0.0_8
      QT(ISTbeg-1) = PSI_C(ISTbeg-1,J,K)

      DO I = ISTbeg, ISTend
        DENOM = AP(I,J,K)-PT(I-1)*AIM(I,J,K)
        PT(I) = AIP(I,J,K) / DENOM
        QT(I) = (TEMPS(I,J,K) + AIM(I,J,K) * QT(I-1))/ DENOM
      ENDDO

      DO I = ISTend, ISTbeg, -1
        PSI_C(I,J,K) = PSI_C(I+1,J,K) * PT(I) + QT(I)
      ENDDO

    ENDIF

  END SUBROUTINE AdiRhoDirPSI_C_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE AdiRhoDirPSI_OUT_d
  USE VAR_PSI_OUT_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2
  REAL(8), DIMENSION (Nmax) :: PT, QT
  REAL(8) :: DENOM


    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= JST) .AND. (J <= (N2-1)) .AND. &
         (K >= KST) .AND. (K <= (N3-1)) ) THEN

      PT(ISTbeg-1) = 0.0_8
      QT(ISTbeg-1) = PSI_OUT(ISTbeg-1,J,K)

      DO I = ISTbeg, ISTend
        DENOM = AP(I,J,K)-PT(I-1)*AIM(I,J,K)
        PT(I) = AIP(I,J,K) / DENOM
        QT(I) = (TEMPS(I,J,K) + AIM(I,J,K) * QT(I-1))/ DENOM
      ENDDO

      DO I = ISTend, ISTbeg, -1
        PSI_OUT(I,J,K) = PSI_OUT(I+1,J,K) * PT(I) + QT(I)
      ENDDO

    ENDIF

  END SUBROUTINE AdiRhoDirPSI_OUT_d

  !------------------------------------------------
 
  ATTRIBUTES (GLOBAL) SUBROUTINE AdiThetaDirPSI_C_d(Iend)
  USE VAR_d, ONLY: TEMPS, CON, AP, AJM, AJP, PSI_C, N2, N3, Nmax
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, VALUE :: Iend
  INTEGER :: ISTbeg, ISTend, JST, KST
  REAL(8), DIMENSION (Nmax) :: PT, QT
  REAL(8) :: DENOM

    ISTbeg = 2; ISTend = Iend-1; JST = 2; KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (K >= KST) .AND. (K <= (N3-1)) ) THEN

      PT(JST-1) = 0.0_8
      QT(JST-1) = PSI_C(I,JST-1,K)

      DO J = JST, N2-1
        DENOM = AP(I,J,K) - PT(J-1) * AJM(I,J,K)
        PT(J) = AJP(I,J,K) / DENOM
        QT(J) = (TEMPS(I,J,K) + AJM(I,J,K) * QT(J-1)) / DENOM
      ENDDO

      DO J = N2-1, JST, -1
        PSI_C(I,J,K) = PSI_C(I,J+1,K) * PT(J) + QT(J)
      ENDDO

    ENDIF

  END SUBROUTINE AdiThetaDirPSI_C_d


  !------------------------------------------------
 
  ATTRIBUTES (GLOBAL) SUBROUTINE AdiThetaDirPSI_OUT_d
  USE VAR_PSI_OUT_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2
  REAL(8), DIMENSION (Nmax) :: PT, QT
  REAL(8) :: DENOM


    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (K >= KST) .AND. (K <= (N3-1)) ) THEN

      PT(JST-1) = 0.0_8
      QT(JST-1) = PSI_OUT(I,JST-1,K)

      DO J = JST, N2-1
        DENOM = AP(I,J,K) - PT(J-1) * AJM(I,J,K)
        PT(J) = AJP(I,J,K) / DENOM
        QT(J) = (TEMPS(I,J,K) + AJM(I,J,K) * QT(J-1)) / DENOM
      ENDDO

      DO J = N2-1, JST, -1
        PSI_OUT(I,J,K) = PSI_OUT(I,J+1,K) * PT(J) + QT(J)
      ENDDO

    ENDIF

  END SUBROUTINE AdiThetaDirPSI_OUT_d
 
  !-------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE AdiPhiDirPSI_C_d(Iend)
  USE VAR_d, ONLY: TEMPS, CON, AP, AKM, AKP, PSI_C, N2, N3, Nmax
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, VALUE :: Iend
  INTEGER :: ISTbeg, ISTend, JST, KST
  REAL(8), DIMENSION (Nmax) :: PT, QT, RT, ALFA, BETA 
  REAL(8) :: DENOM

    ISTbeg = 2; ISTend = Iend-1; JST = 2; KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) ) THEN

      PT(2) = AKP(I,J,2)/AP(I,J,2)
      RT(2) = AKM(I,J,2)/AP(I,J,2)
      QT(2) = TEMPS(I,J,2)/AP(I,J,2)
   
      DO K = 3, N3-2
        DENOM = AP(I,J,K)-PT(K-1)*AKM(I,J,K)
        PT(K) = AKP(I,J,K)/DENOM
        RT(K) = AKM(I,J,K)*RT(K-1)/DENOM
        QT(K) = (TEMPS(I,J,K)+AKM(I,J,K)*QT(K-1))/DENOM
      ENDDO
      
      ALFA(N3-2) = PT(N3-2)+RT(N3-2)
      BETA(N3-2) = QT(N3-2)

      DO K = N3-3, 2, -1
        ALFA(K) = PT(K)*ALFA(K+1)+RT(K)
        BETA(K) = PT(K)*BETA(K+1)+QT(K)
      ENDDO

      PSI_C(I,J,N3-1) = (AKP(I,J,N3-1)*BETA(2)+AKM(I,J,N3-1)*BETA(N3-2)+TEMPS(I,J,N3-1))/ &
                      (AP(I,J,N3-1)-AKP(I,J,N3-1)*ALFA(2)-AKM(I,J,N3-1)*ALFA(N3-2))         

      DO K = 2, N3-2
        PSI_C(I,J,K) = ALFA(K) * PSI_C(I,J,N3-1) + BETA(K)
      ENDDO

    ENDIF

  END SUBROUTINE AdiPhiDirPSI_C_d


  !-------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE AdiPhiDirPSI_OUT_d
  USE VAR_PSI_OUT_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2
  REAL(8), DIMENSION (Nmax) :: PT, QT, RT, ALFA, BETA 
  REAL(8) :: DENOM

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) ) THEN

      PT(2) = AKP(I,J,2)/AP(I,J,2)
      RT(2) = AKM(I,J,2)/AP(I,J,2)
      QT(2) = TEMPS(I,J,2)/AP(I,J,2)
   
      DO K = 3, N3-2
        DENOM = AP(I,J,K)-PT(K-1)*AKM(I,J,K)
        PT(K) = AKP(I,J,K)/DENOM
        RT(K) = AKM(I,J,K)*RT(K-1)/DENOM
        QT(K) = (TEMPS(I,J,K)+AKM(I,J,K)*QT(K-1))/DENOM
      ENDDO
      
      ALFA(N3-2) = PT(N3-2)+RT(N3-2)
      BETA(N3-2) = QT(N3-2)

      DO K = N3-3, 2, -1
        ALFA(K) = PT(K)*ALFA(K+1)+RT(K)
        BETA(K) = PT(K)*BETA(K+1)+QT(K)
      ENDDO

      PSI_OUT(I,J,N3-1) = (AKP(I,J,N3-1)*BETA(2)+AKM(I,J,N3-1)*BETA(N3-2)+TEMPS(I,J,N3-1))/ &
                      (AP(I,J,N3-1)-AKP(I,J,N3-1)*ALFA(2)-AKM(I,J,N3-1)*ALFA(N3-2))         

      DO K = 2, N3-2
        PSI_OUT(I,J,K) = ALFA(K) * PSI_OUT(I,J,N3-1) + BETA(K)
      ENDDO

    ENDIF

  END SUBROUTINE AdiPhiDirPSI_OUT_d

  !-------------------------------------------------

  FUNCTION ComputeNorma2PSI_C_d
  USE VAR_d
  IMPLICIT NONE
  REAL(8) :: ComputeNorma2PSI_C_d, sum_sq_f_h
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2

   
    sum_sq_f = 0.0_8

    !$cuf kernel do (2) <<<*, *>>>
    DO J = JST, N2-1
      DO I = ISTbeg, ISTend
        sum_sq_f = sum_sq_f +                                                               &
           ( AP(I,J,2)*PSI_C(I,J,2) - ( AIP(I,J,2)*PSI_C(I+1,J,2)  + &
                                        AIM(I,J,2)*PSI_C(I-1,J,2)  + &   
                                        AJP(I,J,2)*PSI_C(I,J+1,2)  + &
                                        AJM(I,J,2)*PSI_C(I,J-1,2)  + &
                                        AKP(I,J,2)*PSI_C(I,J,3)    + &
                                        AKM(I,J,2)*PSI_C(I,J,N3-1) + &
                                        CON(I,J,2) ) )**2
      ENDDO
    ENDDO

    !$cuf kernel do (2) <<<*, *>>>
    DO J = JST, N2-1
      DO I = ISTbeg, ISTend
        sum_sq_f = sum_sq_f +                                                   &
           ( AP(I,J,N3-1)*PSI_C(I,J,N3-1) - ( AIP(I,J,N3-1)*PSI_C(I+1,J,N3-1) + &
                                              AIM(I,J,N3-1)*PSI_C(I-1,J,N3-1) + &     
                                              AJP(I,J,N3-1)*PSI_C(I,J+1,N3-1) + &
                                              AJM(I,J,N3-1)*PSI_C(I,J-1,N3-1) + &
                                              AKP(I,J,N3-1)*PSI_C(I,J,2)      + &
                                              AKM(I,J,N3-1)*PSI_C(I,J,N3-2)   + &
                                              CON(I,J,N3-1) ) )**2
      ENDDO
    ENDDO

    !$cuf kernel do (3) <<<*, *>>>
    DO K = 3, N3-2
      DO J = JST, N2-1
        DO I = ISTbeg, ISTend
          sum_sq_f = sum_sq_f +                                        &
              ( AP(I,J,K)*PSI_C(I,J,K) - ( AIP(I,J,K)*PSI_C(I+1,J,K) + &
                                           AIM(I,J,K)*PSI_C(I-1,J,K) + &
                                           AJP(I,J,K)*PSI_C(I,J+1,K) + &
                                           AJM(I,J,K)*PSI_C(I,J-1,K) + &
                                           AKP(I,J,K)*PSI_C(I,J,K+1) + & 
                                           AKM(I,J,K)*PSI_C(I,J,K-1) + &
                                           CON(I,J,K) ) )**2
        ENDDO
      ENDDO
    ENDDO

    sum_sq_f_h = sum_sq_f
    ComputeNorma2PSI_C_d = DSQRT(sum_sq_f_h)

  END FUNCTION ComputeNorma2PSI_C_d

  !------------------------------------------------


  FUNCTION ComputeNormaMaxPSI_C_d(Iend)
  USE VAR_d, ONLY: PSI_C, CON, AP, AIM, AIP, AJM, AJP, AKM, AKP, VCV, N2, N3, sum_sq_f
  IMPLICIT NONE
  REAL(8) :: ComputeNormaMaxPSI_C_d, sum_sq_f_h
  INTEGER :: I, J, K, Iend
  INTEGER :: ISTbeg, ISTend, JST, KST


    ISTbeg = 2; ISTend = Iend-1; JST = 2; KST = 2

    sum_sq_f = 0.0_8
    !$cuf kernel do (2) <<<*, *>>>
    DO J = JST, N2-1
      DO I = ISTbeg, ISTend
        sum_sq_f = max(sum_sq_f, DABS(AP(I,J,2)*PSI_C(I,J,2) - ( AIP(I,J,2)*PSI_C(I+1,J,2)  + &
                                        AIM(I,J,2)*PSI_C(I-1,J,2)  + &   
                                        AJP(I,J,2)*PSI_C(I,J+1,2)  + &
                                        AJM(I,J,2)*PSI_C(I,J-1,2)  + &
                                        AKP(I,J,2)*PSI_C(I,J,3)    + &
                                        AKM(I,J,2)*PSI_C(I,J,N3-1) + &
                                        CON(I,J,2) ))/VCV(I,J,2))
      ENDDO
    ENDDO

    !$cuf kernel do (2) <<<*, *>>>
    DO J = JST, N2-1
      DO I = ISTbeg, ISTend
        sum_sq_f = max(sum_sq_f, DABS(AP(I,J,N3-1)*PSI_C(I,J,N3-1) - ( AIP(I,J,N3-1)*PSI_C(I+1,J,N3-1) + &
                                              AIM(I,J,N3-1)*PSI_C(I-1,J,N3-1) + &     
                                              AJP(I,J,N3-1)*PSI_C(I,J+1,N3-1) + &
                                              AJM(I,J,N3-1)*PSI_C(I,J-1,N3-1) + &

                                              AKP(I,J,N3-1)*PSI_C(I,J,2)      + &
                                              AKM(I,J,N3-1)*PSI_C(I,J,N3-2)   + &
                                              CON(I,J,N3-1) ) )/VCV(I,J,N3-1))
      ENDDO
    ENDDO

    !$cuf kernel do (3) <<<*, *>>>
    DO K = 3, N3-2
      DO J = JST, N2-1
        DO I = ISTbeg, ISTend
          sum_sq_f = max(sum_sq_f, DABS(AP(I,J,K)*PSI_C(I,J,K) - ( AIP(I,J,K)*PSI_C(I+1,J,K) + &
                                           AIM(I,J,K)*PSI_C(I-1,J,K) + &
                                           AJP(I,J,K)*PSI_C(I,J+1,K) + &
                                           AJM(I,J,K)*PSI_C(I,J-1,K) + &
                                           AKP(I,J,K)*PSI_C(I,J,K+1) + & 
                                           AKM(I,J,K)*PSI_C(I,J,K-1) + &
                                           CON(I,J,K) ) )/VCV(I,J,K))
        ENDDO
      ENDDO
    ENDDO

    sum_sq_f_h = sum_sq_f

    ComputeNormaMaxPSI_C_d = sum_sq_f_h

  END FUNCTION ComputeNormaMaxPSI_C_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE FF1N1_d(F, N1, Iend)
    USE VAR_d, ONLY: N2, N3
    IMPLICIT NONE
    INTEGER :: I, J
    INTEGER, VALUE :: N1, Iend
    REAL(8), DIMENSION (N1, N2, N3) :: F


    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= 1) .AND. (I <= Iend) .AND. &
         (J >= 1) .AND. (J <= N2) ) THEN

      F(I,J,1)  = F(I,J,N3-1)
      F(I,J,N3) = F(I,J,2)  
    ENDIF

  END SUBROUTINE FF1N1_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE F1ext_d
    USE VAR_d
    IMPLICIT NONE
    INTEGER :: I, J

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= 1) .AND. (I <= N1) .AND. &
         (J >= 2) .AND. (J <= N2) ) THEN

      F_1(I,J,0)  = F1(I,J,N3-2)
      F_1(I,J,1)  = F1(I,J,N3-1)
      F_1(I,J,2:N3-1) = F1(I,J,2:N3-1)
      F_1(I,J,N3) = F1(I,J,2) 
      F_1(I,J,N3+1) = F1(I,J,3)  
    ENDIF

  END SUBROUTINE F1ext_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE F2ext_d
    USE VAR_d
    IMPLICIT NONE
    INTEGER :: I, J

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= 1) .AND. (I <= N1) .AND. &
         (J >= 2) .AND. (J <= N2) ) THEN

      F_2(I,J,0)  = F2(I,J,N3-2)
      F_2(I,J,1)  = F2(I,J,N3-1)
      F_2(I,J,2:N3-1) = F2(I,J,2:N3-1)
      F_2(I,J,N3) = F2(I,J,2) 
      F_2(I,J,N3+1) = F2(I,J,3)  
    ENDIF

  END SUBROUTINE F2ext_d


  !-------------------------------------------------

  FUNCTION ComputeNorma2PSI_OUT_d
  USE VAR_PSI_OUT_d
  IMPLICIT NONE
  REAL(8) :: ComputeNorma2PSI_OUT_d, sum_sq_f_h
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2

   
    sum_sq_f = 0.0_8

    !$cuf kernel do (2) <<<*, *>>>
    DO J = JST, N2-1
      DO I = ISTbeg, ISTend
        sum_sq_f = sum_sq_f +                                                               &
           ( AP(I,J,2)*PSI_OUT(I,J,2) - ( AIP(I,J,2)*PSI_OUT(I+1,J,2)  + &
                                        AIM(I,J,2)*PSI_OUT(I-1,J,2)  + &   
                                        AJP(I,J,2)*PSI_OUT(I,J+1,2)  + &
                                        AJM(I,J,2)*PSI_OUT(I,J-1,2)  + &
                                        AKP(I,J,2)*PSI_OUT(I,J,3)    + &
                                        AKM(I,J,2)*PSI_OUT(I,J,N3-1) + &
                                        CON(I,J,2) ) )**2
      ENDDO
    ENDDO

    !$cuf kernel do (2) <<<*, *>>>
    DO J = JST, N2-1
      DO I = ISTbeg, ISTend
        sum_sq_f = sum_sq_f +                                                   &
           ( AP(I,J,N3-1)*PSI_OUT(I,J,N3-1) - ( AIP(I,J,N3-1)*PSI_OUT(I+1,J,N3-1) + &
                                              AIM(I,J,N3-1)*PSI_OUT(I-1,J,N3-1) + &     
                                              AJP(I,J,N3-1)*PSI_OUT(I,J+1,N3-1) + &
                                              AJM(I,J,N3-1)*PSI_OUT(I,J-1,N3-1) + &
                                              AKP(I,J,N3-1)*PSI_OUT(I,J,2)      + &
                                              AKM(I,J,N3-1)*PSI_OUT(I,J,N3-2)   + &
                                              CON(I,J,N3-1) ) )**2
      ENDDO
    ENDDO

    !$cuf kernel do (3) <<<*, *>>>
    DO K = 3, N3-2
      DO J = JST, N2-1
        DO I = ISTbeg, ISTend
          sum_sq_f = sum_sq_f +                                        &
              ( AP(I,J,K)*PSI_OUT(I,J,K) - ( AIP(I,J,K)*PSI_OUT(I+1,J,K) + &
                                           AIM(I,J,K)*PSI_OUT(I-1,J,K) + &
                                           AJP(I,J,K)*PSI_OUT(I,J+1,K) + &
                                           AJM(I,J,K)*PSI_OUT(I,J-1,K) + &
                                           AKP(I,J,K)*PSI_OUT(I,J,K+1) + & 
                                           AKM(I,J,K)*PSI_OUT(I,J,K-1) + &
                                           CON(I,J,K) ) )**2
        ENDDO
      ENDDO
    ENDDO

    sum_sq_f_h = sum_sq_f
    ComputeNorma2PSI_OUT_d = DSQRT(sum_sq_f_h)

  END FUNCTION ComputeNorma2PSI_OUT_d


  !-------------------------------------------------
  
    !------------------------------------------------ 

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsRhoDirPSI_IN_d
  USE VAR_PSI_IN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= ISTbeg) .AND. (I <= ISTend)  .AND. &
         (J >= JST)    .AND. (J <= (N2-1))  .AND. &
         (K >= 3)      .AND. (K <= (N3-2)) ) THEN
      TEMPS(I,J,K) = CON(I,J,K) &
                   + AJP(I,J,K) * PSI_IN(I,J+1,K) &
                   + AJM(I,J,K) * PSI_IN(I,J-1,K) &
                   + AKP(I,J,K) * PSI_IN(I,J,K+1) &
                   + AKM(I,J,K) * PSI_IN(I,J,K-1)
    ENDIF

  END SUBROUTINE TempsRhoDirPSI_IN_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsBoundRhoDirPSI_IN_d
  USE VAR_PSI_IN_d
  IMPLICIT NONE
  INTEGER :: I, J
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) ) THEN

      TEMPS(I,J,2) = CON(I,J,2) &
                   + AJP(I,J,2) * PSI_IN(I,J+1,2) &
                   + AJM(I,J,2) * PSI_IN(I,J-1,2) &
                   + AKP(I,J,2) * PSI_IN(I,J,3) &
                   + AKM(I,J,2) * PSI_IN(I,J,N3-1)

      TEMPS(I,J,N3-1) = CON(I,J,N3-1) &
                      + AJP(I,J,N3-1) * PSI_IN(I,J+1,N3-1) &
                      + AJM(I,J,N3-1) * PSI_IN(I,J-1,N3-1) &
                      + AKP(I,J,N3-1) * PSI_IN(I,J,2) &
                      + AKM(I,J,N3-1) * PSI_IN(I,J,N3-2)
    ENDIF

  END SUBROUTINE TempsBoundRhoDirPSI_IN_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsBoundThetaDirPSI_IN_d
  USE VAR_PSI_IN_d
  IMPLICIT NONE
  INTEGER :: I, J
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) ) THEN

      TEMPS(I,J,2) = CON(I,J,2) &
                   + AIP(I,J,2) * PSI_IN(I+1,J,2) &
                   + AIM(I,J,2) * PSI_IN(I-1,J,2) &
                   + AKP(I,J,2) * PSI_IN(I,J,3) &
                   + AKM(I,J,2) * PSI_IN(I,J,N3-1)

      TEMPS(I,J,N3-1) = CON(I,J,N3-1) &
                      + AIP(I,J,N3-1) * PSI_IN(I+1,J,N3-1) &
                      + AIM(I,J,N3-1) * PSI_IN(I-1,J,N3-1) &
                      + AKP(I,J,N3-1) * PSI_IN(I,J,2) &
                      + AKM(I,J,N3-1) * PSI_IN(I,J,N3-2)
    ENDIF

  END SUBROUTINE TempsBoundThetaDirPSI_IN_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsThetaDirPSI_IN_d
  USE VAR_PSI_IN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST)    .AND. (J <= (N2-1)) .AND. &
         (K >= 3)      .AND. (K <= (N3-2)) ) THEN

      TEMPS(I,J,K) = CON(I,J,K) &
                   + AIP(I,J,K) * PSI_IN(I+1,J,K) &
                   + AIM(I,J,K) * PSI_IN(I-1,J,K) &
                   + AKP(I,J,K) * PSI_IN(I,J,K+1) &
                   + AKM(I,J,K) * PSI_IN(I,J,K-1)
    ENDIF

  END SUBROUTINE TempsThetaDirPSI_IN_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsBoundPhiDirPSI_IN_d
  USE VAR_PSI_IN_d
  IMPLICIT NONE
  INTEGER :: I, J
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) ) THEN

      TEMPS(I,J,2) = CON(I,J,2) &
                   + AIP(I,J,2) * PSI_IN(I+1,J,2) &
                   + AIM(I,J,2) * PSI_IN(I-1,J,2) &
                   + AJP(I,J,2) * PSI_IN(I,J+1,2) &
                   + AJM(I,J,2) * PSI_IN(I,J-1,2)


      TEMPS(I,J,N3-1) = CON(I,J,N3-1) &
                    + AIP(I,J,N3-1) * PSI_IN(I+1,J,N3-1) &
                    + AIM(I,J,N3-1) * PSI_IN(I-1,J,N3-1) &
                    + AJP(I,J,N3-1) * PSI_IN(I,J+1,N3-1) &
                    + AJM(I,J,N3-1) * PSI_IN(I,J-1,N3-1)
    ENDIF

  END SUBROUTINE TempsBoundPhiDirPSI_IN_d
 
  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE TempsPhiDirPSI_IN_d
  USE VAR_PSI_IN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2
  
    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) .AND. &
         (K >= 3) .AND. (K <= (N3-2)) ) THEN
      TEMPS(I,J,K) = CON(I,J,K) &
                   + AIP(I,J,K) * PSI_IN(I+1,J,K) &
                   + AIM(I,J,K) * PSI_IN(I-1,J,K) &
                   + AJP(I,J,K) * PSI_IN(I,J+1,K) &
                   + AJM(I,J,K) * PSI_IN(I,J-1,K)
   ENDIF

  END SUBROUTINE TempsPhiDirPSI_IN_d

  !------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE AdiRhoDirPSI_IN_d
  USE VAR_PSI_IN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2
  REAL(8), DIMENSION (Nmax) :: PT, QT
  REAL(8) :: DENOM


    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= JST) .AND. (J <= (N2-1)) .AND. &
         (K >= KST) .AND. (K <= (N3-1)) ) THEN

      PT(ISTbeg-1) = 0.0_8
      QT(ISTbeg-1) = PSI_IN(ISTbeg-1,J,K)

      DO I = ISTbeg, ISTend
        DENOM = AP(I,J,K)-PT(I-1)*AIM(I,J,K)
        PT(I) = AIP(I,J,K) / DENOM
        QT(I) = (TEMPS(I,J,K) + AIM(I,J,K) * QT(I-1))/ DENOM
      ENDDO

      DO I = ISTend, ISTbeg, -1
        PSI_IN(I,J,K) = PSI_IN(I+1,J,K) * PT(I) + QT(I)
      ENDDO

    ENDIF

  END SUBROUTINE AdiRhoDirPSI_IN_d

  !------------------------------------------------
 
  ATTRIBUTES (GLOBAL) SUBROUTINE AdiThetaDirPSI_IN_d
  USE VAR_PSI_IN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2
  REAL(8), DIMENSION (Nmax) :: PT, QT
  REAL(8) :: DENOM


    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (K >= KST) .AND. (K <= (N3-1)) ) THEN

      PT(JST-1) = 0.0_8
      QT(JST-1) = PSI_IN(I,JST-1,K)

      DO J = JST, N2-1
        DENOM = AP(I,J,K) - PT(J-1) * AJM(I,J,K)
        PT(J) = AJP(I,J,K) / DENOM
        QT(J) = (TEMPS(I,J,K) + AJM(I,J,K) * QT(J-1)) / DENOM
      ENDDO

      DO J = N2-1, JST, -1
        PSI_IN(I,J,K) = PSI_IN(I,J+1,K) * PT(J) + QT(J)
      ENDDO

    ENDIF

  END SUBROUTINE AdiThetaDirPSI_IN_d
 
  !-------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE AdiPhiDirPSI_IN_d
  USE VAR_PSI_IN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2
  REAL(8), DIMENSION (Nmax) :: PT, QT, RT, ALFA, BETA 
  REAL(8) :: DENOM

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= ISTbeg) .AND. (I <= ISTend) .AND. &
         (J >= JST) .AND. (J <= (N2-1)) ) THEN

      PT(2) = AKP(I,J,2)/AP(I,J,2)
      RT(2) = AKM(I,J,2)/AP(I,J,2)
      QT(2) = TEMPS(I,J,2)/AP(I,J,2)
   
      DO K = 3, N3-2
        DENOM = AP(I,J,K)-PT(K-1)*AKM(I,J,K)
        PT(K) = AKP(I,J,K)/DENOM
        RT(K) = AKM(I,J,K)*RT(K-1)/DENOM
        QT(K) = (TEMPS(I,J,K)+AKM(I,J,K)*QT(K-1))/DENOM
      ENDDO
      
      ALFA(N3-2) = PT(N3-2)+RT(N3-2)
      BETA(N3-2) = QT(N3-2)

      DO K = N3-3, 2, -1
        ALFA(K) = PT(K)*ALFA(K+1)+RT(K)
        BETA(K) = PT(K)*BETA(K+1)+QT(K)
      ENDDO

      PSI_IN(I,J,N3-1) = (AKP(I,J,N3-1)*BETA(2)+AKM(I,J,N3-1)*BETA(N3-2)+TEMPS(I,J,N3-1))/ &
                      (AP(I,J,N3-1)-AKP(I,J,N3-1)*ALFA(2)-AKM(I,J,N3-1)*ALFA(N3-2))         

      DO K = 2, N3-2
        PSI_IN(I,J,K) = ALFA(K) * PSI_IN(I,J,N3-1) + BETA(K)
      ENDDO

    ENDIF

  END SUBROUTINE AdiPhiDirPSI_IN_d

  !-------------------------------------------------

  FUNCTION ComputeNorma2PSI_IN_d
  USE VAR_PSI_IN_d
  IMPLICIT NONE
  REAL(8) :: ComputeNorma2PSI_IN_d, sum_sq_f_h
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2

   
    sum_sq_f = 0.0_8

    !$cuf kernel do (2) <<<*, *>>>
    DO J = JST, N2-1
      DO I = ISTbeg, ISTend
        sum_sq_f = sum_sq_f +                                                               &
           ( AP(I,J,2)*PSI_IN(I,J,2) - ( AIP(I,J,2)*PSI_IN(I+1,J,2)  + &
                                        AIM(I,J,2)*PSI_IN(I-1,J,2)  + &   
                                        AJP(I,J,2)*PSI_IN(I,J+1,2)  + &
                                        AJM(I,J,2)*PSI_IN(I,J-1,2)  + &
                                        AKP(I,J,2)*PSI_IN(I,J,3)    + &
                                        AKM(I,J,2)*PSI_IN(I,J,N3-1) + &
                                        CON(I,J,2) ) )**2
      ENDDO
    ENDDO

    !$cuf kernel do (2) <<<*, *>>>
    DO J = JST, N2-1
      DO I = ISTbeg, ISTend
        sum_sq_f = sum_sq_f +                                                   &
           ( AP(I,J,N3-1)*PSI_IN(I,J,N3-1) - ( AIP(I,J,N3-1)*PSI_IN(I+1,J,N3-1) + &
                                              AIM(I,J,N3-1)*PSI_IN(I-1,J,N3-1) + &     
                                              AJP(I,J,N3-1)*PSI_IN(I,J+1,N3-1) + &
                                              AJM(I,J,N3-1)*PSI_IN(I,J-1,N3-1) + &
                                              AKP(I,J,N3-1)*PSI_IN(I,J,2)      + &
                                              AKM(I,J,N3-1)*PSI_IN(I,J,N3-2)   + &
                                              CON(I,J,N3-1) ) )**2
      ENDDO
    ENDDO

    !$cuf kernel do (3) <<<*, *>>>
    DO K = 3, N3-2
      DO J = JST, N2-1
        DO I = ISTbeg, ISTend
          sum_sq_f = sum_sq_f +                                        &
              ( AP(I,J,K)*PSI_IN(I,J,K) - ( AIP(I,J,K)*PSI_IN(I+1,J,K) + &
                                           AIM(I,J,K)*PSI_IN(I-1,J,K) + &
                                           AJP(I,J,K)*PSI_IN(I,J+1,K) + &
                                           AJM(I,J,K)*PSI_IN(I,J-1,K) + &
                                           AKP(I,J,K)*PSI_IN(I,J,K+1) + & 
                                           AKM(I,J,K)*PSI_IN(I,J,K-1) + &
                                           CON(I,J,K) ) )**2
        ENDDO
      ENDDO
    ENDDO

    sum_sq_f_h = sum_sq_f
    ComputeNorma2PSI_IN_d = DSQRT(sum_sq_f_h)

  END FUNCTION ComputeNorma2PSI_IN_d

  !-------------------------------------------------

  FUNCTION ComputeNorma2PSI_IN_Bound_d
  USE VAR_PSI_IN_d
  IMPLICIT NONE
  REAL(8) :: ComputeNorma2PSI_IN_Bound_d, sum_sq_f_h
  INTEGER :: I, J, K
  INTEGER, PARAMETER :: ISTbeg = 2, ISTend = N1-1, JST = 2, KST = 2

   
    sum_sq_f = 0.0_8

    !$cuf kernel do  <<<*, *>>>
    DO J = JST, N2-1
        sum_sq_f = sum_sq_f +                                                               &
           ( AP(N1-1,J,2)*PSI_IN(N1-1,J,2) - ( AIP(N1-1,J,2)*PSI_IN(N1,J,2)  + &
                                               AIM(N1-1,J,2)*PSI_IN(N1-2,J,2)  + &   
                                               AJP(N1-1,J,2)*PSI_IN(N1-1,J+1,2)  + &
                                               AJM(N1-1,J,2)*PSI_IN(N1-1,J-1,2)  + &
                                               AKP(N1-1,J,2)*PSI_IN(N1-1,J,3)    + &
                                               AKM(N1-1,J,2)*PSI_IN(N1-1,J,N3-1) + &
                                               CON(N1-1,J,2) ) )**2
    ENDDO

    !$cuf kernel do <<<*, *>>>
    DO J = JST, N2-1
        sum_sq_f = sum_sq_f +                                                   &
                   ( AP(N1-1,J,N3-1)*PSI_IN(N1-1,J,N3-1) - ( AIP(N1-1,J,N3-1)*PSI_IN(N1,J,N3-1) + &
                                              AIM(N1-1,J,N3-1)*PSI_IN(N1-2,J,N3-1) + &     
                                              AJP(N1-1,J,N3-1)*PSI_IN(N1-1,J+1,N3-1) + &
                                              AJM(N1-1,J,N3-1)*PSI_IN(N1-1,J-1,N3-1) + &
                                              AKP(N1-1,J,N3-1)*PSI_IN(N1-1,J,2)      + &
                                              AKM(N1-1,J,N3-1)*PSI_IN(N1-1,J,N3-2)   + &
                                              CON(N1-1,J,N3-1) ) )**2
    ENDDO

    !$cuf kernel do (2) <<<*, *>>>
    DO K = 3, N3-2
      DO J = JST, N2-1
          sum_sq_f = sum_sq_f +                                        &
              ( AP(N1-1,J,K)*PSI_IN(N1-1,J,K) - ( AIP(N1-1,J,K)*PSI_IN(N1,J,K) + &
                                           AIM(N1-1,J,K)*PSI_IN(N1-2,J,K) + &
                                           AJP(N1-1,J,K)*PSI_IN(N1-1,J+1,K) + &
                                           AJM(N1-1,J,K)*PSI_IN(N1-1,J-1,K) + &
                                           AKP(N1-1,J,K)*PSI_IN(N1-1,J,K+1) + & 
                                           AKM(N1-1,J,K)*PSI_IN(N1-1,J,K-1) + &
                                           CON(N1-1,J,K) ) )**2

      ENDDO
    ENDDO

    sum_sq_f_h = sum_sq_f
    ComputeNorma2PSI_IN_Bound_d = DSQRT(sum_sq_f_h)

  END FUNCTION ComputeNorma2PSI_IN_Bound_d

  !-------------------------------------------------

end module solver
