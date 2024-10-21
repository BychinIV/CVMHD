MODULE FD
IMPLICIT NONE
  CONTAINS
  
  !------------------------------------------------------
  
  SUBROUTINE SIMPLER
    USE VAR, ONLY : SMAX, NTIMES_UVW, NTIMES_P, NTIMES_PC, NTIMES_T, N1_core                       
    USE VAR_d
    USE MD, ONLY: MHD, ComputeRotB
    USE SOLVER, ONLY : AdiSolver
    USE USER, ONLY : ComputeSourceU1, ComputeSourceU2, ComputeSourceU3, &
                     ComputeSourceP, ComputeSourcePC, ComputeSourceT 
    IMPLICIT NONE
    INTEGER :: NUM_ITERS, OUT_ITER

    
    !DO OUT_ITER = 1, 2

      SMAX = 1.0_8
      CALL ComputeRotB
      CALL ComputeFOnBot_d <<<GRD_IK, TBL_IK>>> (U1)
      CALL ComputeFOnBot_d <<<GRD_IK, TBL_IK>>> (U3)
      CALL ComputeFOnBot_d <<<GRD_IK, TBL_IK>>> (T)

        
      DO WHILE (SMAX > 1.D-3) !( (SMAX > 1.D-6) .OR. (Rsum > 1.0_8) )


        !-----------------Step 1------------------

        CALL ComputeSourceU1
        CALL CofU1

        CALL ComputeSourceU2
        CALL CofU2

        CALL ComputeSourceU3
        CALL CofU3

        CALL ComputeSourceP
        CALL CofP
        NUM_ITERS = ADISolver(P, N1_core, N1-1, 2, 2, NTIMES_P)

        !-----------------Step 2---------------

        U1tmp = U1; U2tmp = U2; U3tmp = U3

        CALL ComputeSourceU1
        CALL CofU1
        NUM_ITERS = AdiSolver(U1tmp, N1_core + 1, N1-1, 2, 2, NTIMES_UVW)
        CALL ComputeFOnBot_d <<<GRD_IK, TBL_IK>>> (U1tmp)

        CALL ComputeSourceU2
        CALL CofU2
        NUM_ITERS = ADISolver(U2tmp, N1_core, N1-1, 3, 2, NTIMES_UVW)

        CALL ComputeSourceU3
        CALL CofU3
        NUM_ITERS = ADISolver(U3tmp, N1_core, N1-1, 2, 2, NTIMES_UVW)
        CALL ComputeFOnBot_d <<<GRD_IK, TBL_IK>>> (U3tmp)

        !-----------------Step 3---------------

        CALL ComputeSourcePC
        CALL CofPC
        NUM_ITERS = ADISolver(PC, N1_core, N1-1, 2, 2, NTIMES_PC)


        !-----------------Step 4---------------

        U1 = U1tmp; U2 = U2tmp; U3 = U3tmp
        CALL CorrectU

        CALL ComputeSmax

        !78  FORMAT(10X,'Smax')
        !    WRITE(*, 78)
        !    WRITE(*, '(1P1E15.4)') SMAX

      ENDDO 

      !-----------------Step 5---------------
      CALL MHD
      CALL ComputeRotB

      CALL ComputeSourceT
      CALL CofT
      NUM_ITERS = ADISolver(T, N1_core, N1-1, 2, 2, NTIMES_T)
      CALL ComputeFOnBot_d <<<GRD_IK, TBL_IK>>> (T)

    !END DO

  END SUBROUTINE SIMPLER

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeFOnBot_d(F)
    USE VAR_d
    IMPLICIT NONE
    INTEGER :: I, K
    REAL(8), DIMENSION (N1,N2,N3) :: F

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= 1) .AND. (I <= N1)  .AND. &
         (K >= 1) .AND. (K <= N3) ) THEN

      F(I,N2,K) = F(I,N2-1,K)

    ENDIF

  END SUBROUTINE ComputeFOnBot_d

  !------------------------------------------------------

  SUBROUTINE CofU1
    USE VAR_d
    IMPLICIT NONE

    CALL CofAimAipU1_d <<<GRD_IJK, TBL_IJK>>>
    CALL CofAimAipBoundU1_d  <<<GRD_JK, TBL_JK>>>
 
    CALL CofAjmAjpU1_d <<<GRD_IJK, TBL_IJK>>>
    CALL CofAjmAjpOnSpheresU1_d <<<GRD_JK, TBL_JK>>>
    CALL CofAjpOnBotU1_d <<<GRD_IK, TBL_IK>>>
    
    CALL CofAkmAkpU1_d <<<GRD_IJK, TBL_IJK>>>
    CALL CofAkmAkpBoundU1_d <<<GRD_JK, TBL_JK>>>
    
    CALL CofConApU1hatDU1_d <<<GRD_IJK, TBL_IJK>>>
    CALL CofBoundDU1_d <<<GRD_JK, TBL_JK>>>
    CALL CofConU1_d <<<GRD_IJK, TBL_IJK>>>

  END SUBROUTINE CofU1

  !------------------------------------------------------

  SUBROUTINE CofU2
    USE VAR_d
    IMPLICIT NONE

    CALL CofAimAipU2_d <<<GRD_IJK, TBL_IJK>>>
    CALL CofAimAipBoundU2_d <<<GRD_IK, TBL_IK>>>
    
    CALL CofAjmAjpU2_d <<<GRD_IJK, TBL_IJK>>>
    CALL CofAjmAjpBoundU2_d <<<GRD_IK, TBL_IK>>>
    
    CALL CofAkmAkpU2_d <<<GRD_IJK, TBL_IJK>>>
    CALL CofAkmAkpBoundU2_d <<<GRD_IK, TBL_IK>>>
    
    CALL CofConApU2hatDU2_d <<<GRD_IJK, TBL_IJK>>>
    CALL CofBoundDU2_d <<<GRD_IK, TBL_IK>>>
    CALL CofConU2_d <<<GRD_IJK, TBL_IJK>>>

  END SUBROUTINE CofU2

  !------------------------------------------------------

  SUBROUTINE CofU3
    USE VAR_d
    IMPLICIT NONE

    CALL CofU3_d <<<GRD_IJK, TBL_IJK>>>
    CALL CofAjpOnBotU3_d <<<GRD_IK, TBL_IK>>>

    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (DU3)
    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (U3hat)

  END SUBROUTINE CofU3

  !------------------------------------------------------

  SUBROUTINE CofP
    USE VAR_d
    IMPLICIT NONE

    CALL CofAimAipConP_d <<<GRD_IJK, TBL_IJK>>>
    CALL CofAimAipConBoundP_d <<<GRD_JK, TBL_JK>>>
    CALL CofAjmAjpConP_d <<<GRD_IJK, TBL_IJK>>>
    CALL CofAjmAjpConBoundP_d <<<GRD_IK, TBL_IK>>>
    CALL CofAkmAkpConP_d <<<GRD_IJK, TBL_IJK>>>
    CALL CofApP_d <<<GRD_IJK, TBL_IJK>>>

  END SUBROUTINE CofP

  !------------------------------------------------------

  SUBROUTINE CofPC
    USE VAR_d
    IMPLICIT NONE

    CALL CofAimAipPC_d <<<GRD_IJK, TBL_IJK>>>
    CALL CofAimAipBoundPC_d <<<GRD_JK, TBL_JK>>>
    CALL CofAjmAjpPC_d <<<GRD_IJK, TBL_IJK>>>
    CALL CofAjmAjpBoundPC_d <<<GRD_IK, TBL_IK>>>
    CALL CofAkmAkpPC_d <<<GRD_IJK, TBL_IJK>>>
    CALL CofApConPC_d <<<GRD_IJK, TBL_IJK>>>

  END SUBROUTINE CofPC

  !------------------------------------------------------

  SUBROUTINE CofT
    USE VAR_d
    IMPLICIT NONE

    CALL CofT_d <<<GRD_IJK, TBL_IJK>>>
    CALL CofAjpOnBotT_d <<<GRD_IK, TBL_IK>>>
    
  END SUBROUTINE CofT

  !------------------------------------------------------

  SUBROUTINE CorrectU
    USE VAR_d
    IMPLICIT NONE

    CALL CorrectU1_d <<<GRD_IJK, TBL_IJK>>>
    CALL CorrectU2_d <<<GRD_IJK, TBL_IJK>>>
    CALL CorrectU3_d <<<GRD_IJK, TBL_IJK>>>

    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (U1)
    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (U2)
    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (U3)


    CALL ComputeFOnBot_d <<<GRD_IK, TBL_IK>>> (U1)
    CALL ComputeFOnBot_d <<<GRD_IK, TBL_IK>>> (U3)


  END SUBROUTINE CorrectU

  !------------------------------------------------------
  
  ATTRIBUTES (GLOBAL) SUBROUTINE FF1N1_d(F)
    USE VAR_d
    IMPLICIT NONE
    INTEGER :: I, J
    REAL(8), DIMENSION (N1,N2,N3) :: F

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= 1) .AND. (I <= N1) .AND. &
         (J >= 2) .AND. (J <= N2) ) THEN

      F(I,J,1)  = F(I,J,N3-1)
      F(I,J,N3) = F(I,J,2)  
    ENDIF

  END SUBROUTINE FF1N1_d


  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofT_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: AP0, DEN, FLOW, DIFF 

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (J >= 2) .AND. (J <= (N2-2))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      !----------------- AIM, AIP -----------------
      DIFF = thc * S1(I,J,K) / Xdif1(I)
      FLOW = rho_c * U1(I,J,K) * S1(I,J,K)               
      AIM(I,J,K) = ACOF1(DIFF,FLOW)                               

      DIFF = thc * S1(I+1,J,K) / Xdif1(I+1)
      FLOW = rho_c * U1(I+1,J,K) * S1(I+1,J,K)    
      AIP(I,J,K) = ACOF1(DIFF,-FLOW) 

      !----------------- AJM, AJP -----------------
      DIFF = thc * S2(I,J,K) / (X1(I) * Xdif2(J))
      FLOW = rho_c * U2(I,J,K) * S2(I,J,K)               
      AJM(I,J,K) = ACOF1(DIFF,FLOW)            

      DIFF = thc * S2(I,J+1,K) / (X1(I) * Xdif2(J+1))
      FLOW = rho_c * U2(I,J+1,K) * S2(I,J+1,K)
      AJP(I,J,K) = ACOF1(DIFF,-FLOW)

      !----------------- AKM, AKP -----------------

      DIFF = thc * S3(I,J,K) / (X1(I) * sinX2(J) * Xdif3(K))
      FLOW = rho_c * U3(I,J,K) * S3(I,J,K)        
      AKM(I,J,K) = ACOF1(DIFF,FLOW)

      DIFF = thc * S3(I,J,K+1) / (X1(I) * sinX2(J) * Xdif3(K+1))
      FLOW = rho_c * U3(I,J,K+1) * S3(I,J,K+1)    
      AKP(I,J,K) = ACOF1(DIFF,-FLOW)

      !----------------- AP0, CON, AP -----------------

      AP0 = rho_c * VCV(I,J,K) / DT
      CON(I,J,K) = CON(I,J,K)  + AP0 * T0(I,J,K)
      AP(I,J,K) = - AP(I,J,K)  * VCV(I,J,K) + AP0+ &
                    AIM(I,J,K) + AIP(I,J,K) + &
                    AJM(I,J,K) + AJP(I,J,K) + &
                    AKM(I,J,K) + AKP(I,J,K)
    ENDIF

  END SUBROUTINE CofT_d
  
  !-------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAjpOnBotT_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: AP0, DEN, FLOW, DIFF 

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    J = N2-1

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      !----------------- AIM, AIP -----------------
      DIFF = thc * S1(I,J,K) / Xdif1(I)
      FLOW = rho_c * U1(I,J,K) * S1(I,J,K)               
      AIM(I,J,K) = ACOF1(DIFF,FLOW)                               

      DIFF = thc * S1(I+1,J,K) / Xdif1(I+1)
      FLOW = rho_c * U1(I+1,J,K) * S1(I+1,J,K)    
      AIP(I,J,K) = ACOF1(DIFF,-FLOW) 

      !----------------- AJM, AJP -----------------
      DIFF = thc * S2(I,J,K) / (X1(I) * Xdif2(J))
      FLOW = rho_c * U2(I,J,K) * S2(I,J,K)               
      AJM(I,J,K) = ACOF1(DIFF,FLOW)            

      !DIFF = thc * S2(I,J+1,K) / (X1(I) * Xdif2(J+1))
      !FLOW = rho_c * U2(I,J+1,K) * S2(I,J+1,K)
      AJP(I,J,K) = 0.0_8!ACOF1(DIFF,-FLOW)

      !----------------- AKM, AKP -----------------

      DIFF = thc * S3(I,J,K) / (X1(I) * sinX2(J) * Xdif3(K))
      FLOW = rho_c * U3(I,J,K) * S3(I,J,K)        
      AKM(I,J,K) = ACOF1(DIFF,FLOW)

      DIFF = thc * S3(I,J,K+1) / (X1(I) * sinX2(J) * Xdif3(K+1))
      FLOW = rho_c * U3(I,J,K+1) * S3(I,J,K+1)    
      AKP(I,J,K) = ACOF1(DIFF,-FLOW)

      !----------------- AP0, CON, AP -----------------

      AP0 = rho_c * VCV(I,J,K) / DT
      CON(I,J,K) = CON(I,J,K)  + AP0 * T0(I,J,K)
      AP(I,J,K) = - AP(I,J,K)  * VCV(I,J,K) + AP0+ &
                    AIM(I,J,K) + AIP(I,J,K) + &
                    AJM(I,J,K) + AJP(I,J,K) + &
                    AKM(I,J,K) + AKP(I,J,K)
    ENDIF

  END SUBROUTINE CofAjpOnBotT_d  

  !-------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAimAipU1_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: DIFF, FLOW

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= (N1_core + 1)) .AND. (I <= (N1-2))  .AND. &
         (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      DIFF = visk * S1u1(I-1,J,K) / (Xs1(I) - Xs1(I-1))
      FLOW = dens * 0.5_8 * (U1(I,J,K) * S1(I,J,K) + U1(I-1,J,K) * S1(I-1,J,K))          
      AIM(I,J,K) = ACOF1(DIFF, FLOW)      

      DIFF = visk * S1u1(I,J,K) / (Xs1(I+1) - Xs1(I))
      FLOW = dens * 0.5_8 * (U1(I+1,J,K) * S1(I+1,J,K) + U1(I,J,K) * S1(I,J,K))   
      AIP(I,J,K) = ACOF1(DIFF, -FLOW)
    ENDIF

  END SUBROUTINE CofAimAipU1_d


  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAimAipBoundU1_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: DIFF, FLOW

    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      !I = 3
      !DIFF = visk * S1u1(1,J,K) / (Xs1(3) - Xs1(2))
      !FLOW = dens * U1(2,J,K) * S1u1(1,J,K)          
      !AIM(3,J,K) = ACOF1(DIFF, FLOW)
      
      !DIFF = visk * S1u1(3,J,K) / (Xs1(4) - Xs1(3))   
      !FLOW = dens * 0.5_8 * (U1(4,J,K) * S1(4,J,K) + U1(3,J,K) * S1(3,J,K))          
      !AIP(3,J,K) = ACOF1(DIFF, -FLOW) 

      !I = N1-1
      DIFF = visk * S1u1(N1-2,J,K) / (Xs1(N1-1) - Xs1(N1-2))
      FLOW = dens * 0.5_8 * (U1(N1-1,J,K) * S1(N1-1,J,K) + U1(N1-2,J,K) * S1(N1-2,J,K))          
      AIM(N1-1,J,K) = ACOF1(DIFF, FLOW) 
              
      DIFF = visk * S1u1(N1,J,K) / (Xs1(N1) - Xs1(N1-1))
      FLOW = dens * U1(N1,J,K) * S1u1(N1,J,K)    
      AIP(N1-1,J,K) = ACOF1(DIFF, -FLOW)

    ENDIF

  END SUBROUTINE CofAimAipBoundU1_d

  !-------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAjmAjpU1_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: DIFF, FLOW, FLOWL, FLOWR

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= (N1_core + 1)) .AND. (I <= (N1-2))  .AND. &
         (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      DIFF = visk * S2u1_1(I,J,K) / Xdif2(J)
      FLOWL = U2(I-1,J,K) * (Xs1(I) - X1(I-1)) * (Xs1(I) + X1(I-1))
      FLOWR = U2(I,J,K) * (X1(I) - Xs1(I)) * (X1(I) + Xs1(I))
      FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * sinXs2(J) * (Xs3(K+1) - Xs3(K))
      AJM(I,J,K) = ACOF1(DIFF, FLOW)

      DIFF = visk * S2u1_1(I,J+1,K) / Xdif2(J+1)
      FLOWL = U2(I-1,J+1,K) * (Xs1(I) - X1(I-1)) * (Xs1(I) + X1(I-1))
      FLOWR = U2(I,J+1,K) * (X1(I) - Xs1(I)) * (X1(I) + Xs1(I))
      FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * sinXs2(J+1) * (Xs3(K+1) - Xs3(K))
      AJP(I,J,K) = ACOF1(DIFF, -FLOW) 
      
    ENDIF

  END SUBROUTINE CofAjmAjpU1_d

  !-------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAjmAjpOnSpheresU1_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: DIFF, FLOW, FLOWL, FLOWR

    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      !I = 3
      !DIFF = visk * S2u1_1(3,J,K) / Xdif2(J)
      !FLOWL = U2(2,J,K) * (Xs1(3) - X1(1)) * (Xs1(3) + X1(1))
      !FLOWR = U2(3,J,K) * (X1(3) - Xs1(3)) * (X1(3) + Xs1(3))
      !FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * sinXs2(J) * (Xs3(K+1) - Xs3(K))
      !AJM(3,J,K) = ACOF1(DIFF, FLOW)

      !DIFF = visk * S2u1_1(3,J+1,K) / Xdif2(J+1)
      !FLOWL = U2(2,J+1,K) * (Xs1(3) - X1(1)) * (Xs1(3) + X1(1))
      !FLOWR = U2(3,J+1,K) * (X1(3) - Xs1(3)) * (X1(3) + Xs1(3))
      !FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * sinXs2(J+1) * (Xs3(K+1) - Xs3(K))
      !AJP(3,J,K) = ACOF1(DIFF, -FLOW)            

      !I = N1-1
      DIFF = visk * S2u1_1(N1-1,J,K) / Xdif2(J)
      FLOWL = U2(N1-2,J,K) * (Xs1(N1-1) - X1(N1-2)) * (Xs1(N1-1)+X1(N1-2))
      FLOWR = U2(N1-1,J,K) * (X1(N1) - Xs1(N1-1)) * (X1(N1) + Xs1(N1-1))
      FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * sinXs2(J) * (Xs3(K+1) - Xs3(K))
      AJM(N1-1,J,K) = ACOF1(DIFF, FLOW)

      DIFF = visk * S2u1_1(N1-1,J+1,K) / Xdif2(J+1)
      FLOWL = U2(N1-2,J+1,K) * (Xs1(N1-1) - X1(N1-2)) * (Xs1(N1-1) + X1(N1-2))
      FLOWR = U2(N1-1,J+1,K) * (X1(N1) - Xs1(N1-1)) * (X1(N1) + Xs1(N1-1))
      FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * sinXs2(J+1) * (Xs3(K+1) - Xs3(K))
      AJP(N1-1,J,K) = ACOF1(DIFF, -FLOW)      

    ENDIF

  END SUBROUTINE CofAjmAjpOnSpheresU1_d
  
  !-------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAjpOnBotU1_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= (N1_core + 1)) .AND. (I <= (N1-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      AJP(I,N2-1,K) = 0.0_8
      
    ENDIF

  END SUBROUTINE CofAjpOnBotU1_d  

  !-------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAkmAkpU1_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: DIFF, FLOW, FLOWL, FLOWR

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= (N1_core + 1)) .AND. (I <= (N1-2))  .AND. &
         (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      DIFF = visk * S3u1(I,J,K) / (Xs1(I) * sinX2(J) * Xdif3(K))
      FLOWL = U3(I-1,J,K) * (Xs1(I) - X1(I-1)) * (Xs1(I) + X1(I-1))
      FLOWR = U3(I,J,K) * (X1(I) - Xs1(I)) * (X1(I) + Xs1(I))
      FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * (Xs2(J+1) - Xs2(J))         
      AKM(I,J,K) = ACOF1(DIFF, FLOW)

      DIFF = visk * S3u1(I,J,K+1) / (Xs1(I) * sinX2(J) * Xdif3(K+1))
      FLOWL = U3(I-1,J,K+1) * (Xs1(I) - X1(I-1)) * (Xs1(I) + X1(I-1))
      FLOWR = U3(I,J,K+1) * (X1(I) - Xs1(I)) * (X1(I) + Xs1(I))
      FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * (Xs2(J+1) - Xs2(J))               
      AKP(I,J,K) = ACOF1(DIFF, -FLOW)  

    ENDIF

  END SUBROUTINE CofAkmAkpU1_d

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAkmAkpBoundU1_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: DIFF, FLOW, FLOWL, FLOWR

    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      !I = 3
      !DIFF = visk * S3u1(3,J,K) / (Xs1(3) * sinX2(J) * Xdif3(K))
      !FLOWL = U3(2,J,K) * (Xs1(3) - X1(1)) * (Xs1(3) + X1(1))
      !FLOWR = U3(3,J,K) * (X1(3) - Xs1(3)) * (X1(3) + Xs1(3))
      !FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * (Xs2(J+1) - Xs2(J))        
      !AKM(3,J,K) = ACOF1(DIFF, FLOW)

      !DIFF = visk * S3u1(3,J,K+1) / (Xs1(3) * sinX2(J) * Xdif3(K+1))
      !FLOWL = U3(2,J,K+1) * (Xs1(3) - X1(1)) * (Xs1(3) + X1(1))
      !FLOWR = U3(3,J,K+1) * (X1(3) - Xs1(3)) * (X1(3) + Xs1(3))
      !FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * (Xs2(J+1) - Xs2(J))                
      !AKP(3,J,K) = ACOF1(DIFF, -FLOW)                 

      !I = N1-1
      DIFF = visk * S3u1(N1-1,J,K) / (Xs1(N1-1) * sinX2(J) * Xdif3(K))
      FLOWL = U3(N1-2,J,K) * (Xs1(N1-1) - X1(N1-2)) * (Xs1(N1-1) + X1(N1-2))
      FLOWR = U3(N1-1,J,K) * (X1(N1) - Xs1(N1-1)) * (X1(N1) + Xs1(N1-1))
      FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * (Xs2(J+1) - Xs2(J))       
      AKM(N1-1,J,K) = ACOF1(DIFF, FLOW)

      DIFF = visk * S3u1(N1-1,J,K+1) / (Xs1(N1-1) * sinX2(J) * Xdif3(K+1))
      FLOWL = U3(N1-2,J,K+1) * (Xs1(N1-1) - X1(N1-2)) * (Xs1(N1-1) + X1(N1-2))
      FLOWR = U3(N1-1,J,K+1) * (X1(N1) - Xs1(N1-1)) * (X1(N1) + Xs1(N1-1))
      FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * (Xs2(J+1) - Xs2(J))         
      AKP(N1-1,J,K) = ACOF1(DIFF, -FLOW)  

    ENDIF

  END SUBROUTINE CofAkmAkpBoundU1_d

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofConApU1hatDU1_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: VV11, VV12, VV21, VV22, VV, WW1, WW2, WW3, WW4, &
             WW, VY1, VY2, VY, WZ1, WZ2, WZ, AP0

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= (N1_core + 1)) .AND. (I <= (N1-1))  .AND. &
         (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      VV11 = U2(I,J,K) * U2(I,J,K) * VVX1(I) * VVY1(J)
      VV12 = U2(I,J+1,K) * U2(I,J+1,K) * VVX1(I) * VVY2(J)
      VV21 = U2(I-1,J,K) * U2(I-1,J,K) * VVX2(I) * VVY1(J)
      VV22 = U2(I-1,J+1,K) * U2(I-1,J+1,K) * VVX2(I) * VVY2(J)

      VV = 0.5_8 * (VV11 + VV12 + VV21 + VV22) * (Xs3(K+1) - Xs3(K))

      !-----------------------------------

      WW1 = U3(I,J,K) * U3(I,J,K) * VVX1(I) * (X3(K) - Xs3(K))
      WW2 = U3(I,J,K+1) * U3(I,J,K+1) * VVX1(I) * (Xs3(K+1) - X3(K))
      WW3 = U3(I-1,J,K+1) * U3(I-1,J,K+1) * VVX2(I) * (Xs3(K+1) - X3(K))
      WW4 = U3(I-1,J,K) * U3(I-1,J,K) * VVX2(I) * (X3(K) - Xs3(K))

      WW = 0.5_8 * (WW1 + WW2 + WW3 + WW4) * (cosXs2(J) - cosXs2(J+1))     

      !-----------------------------------

      VY1 = (U2(I,J,K) * VYX1(I) + U2(I-1,J,K) * VYX2(I)) * sinXs2(J)
      VY2 = (U2(I,J+1,K) * VYX1(I) + U2(I-1,J+1,K) * VYX2(I)) * sinXs2(J+1)

      VY = 2.0_8 * (VY2 - VY1) * (Xs3(K+1) - Xs3(K))

      !------------------------------------

      WZ1 = U3(I,J,K) * VYX1(I) + U3(I-1,J,K) * VYX2(I)
      WZ2 = U3(I,J,K+1) * VYX1(I) + U3(I-1,J,K+1) * VYX2(I)

      WZ = 2.0_8 * (WZ2 - WZ1) * (Xs2(J+1) - Xs2(J))

      !------------------------------------

      CON(I,J,K) = CON(I,J,K) + dens * (VV + WW) - visk * (VY + WZ)
      AP(I,J,K) = AP(I,J,K) - 2.0_8 * visk * APX(I) * &
                  (cosXs2(J) - cosXs2(J+1)) * (Xs3(K+1) - Xs3(K))

      !------------------------------------
       
      AP0 = dens * VCV1(I,J,K) / DT
      CON(I,J,K) = CON(I,J,K) + AP0 * U1_0(I,J,K)
      AP(I,J,K) = - AP(I,J,K) + AP0 +         &
                    AIM(I,J,K) + AIP(I,J,K) + &
                    AJM(I,J,K) + AJP(I,J,K) + &
                    AKM(I,J,K) + AKP(I,J,K)
      DU1(I,J,K) = 0.5_8 * (S1u1(I,J,K) + S1u1(I-1,J,K)) / AP(I,J,K)

      U1hat(I,J,K) = (AIP(I,J,K) * U1(I+1,J,K) + AIM(I,J,K) * U1(I-1,J,K) + &
                      AJP(I,J,K) * U1(I,J+1,K) + AJM(I,J,K) * U1(I,J-1,K) + &
                      AKP(I,J,K) * U1(I,J,K+1) + AKM(I,J,K) * U1(I,J,K-1) + &
                      CON(I,J,K)) / AP(I,J,K)   

    ENDIF

  END SUBROUTINE CofConApU1hatDU1_d

  !-------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofBoundDU1_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      !I = 3
      !DU1(3,J,K) = 0.5_8 * (S1u1(3,J,K) + S1u1(1,J,K)) / AP(3,J,K)
     
      !I = N1-1
      DU1(N1-1,J,K) = 0.5_8 * (S1u1(N1,J,K) + S1u1(N1-2,J,K)) / AP(N1-1,J,K)

    ENDIF

  END SUBROUTINE CofBoundDU1_d


  ATTRIBUTES (GLOBAL) SUBROUTINE CofConU1_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= (N1_core + 1)) .AND. (I <= (N1-1))  .AND. &
         (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      CON(I,J,K) = CON(I,J,K) + DU1(I,J,K) * AP(I,J,K) * (P(I-1,J,K) - P(I,J,K))

    ENDIF

  END SUBROUTINE CofConU1_d

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAimAipU2_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: DIFF, FLOWL, FLOWR, FLOW

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (J >= 4) .AND. (J <= (N2-2))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      DIFF = visk * S1u2(I,J,K) / Xdif1(I)
      FLOWL = U1(I,J-1,K) * (cosX2(J-1) - cosXs2(J))
      FLOWR = U1(I,J,K) * (cosXs2(J) - cosX2(J))
      FLOW = dens * (FLOWL + FLOWR) * Xs1(I) * Xs1(I) * (Xs3(K+1) - Xs3(K))               
      AIM(I,J,K) = ACOF1(DIFF, FLOW)

      DIFF = visk * S1u2(I+1,J,K) / Xdif1(I+1)
      FLOWL = U1(I+1,J-1,K) * (cosX2(J-1) - cosXs2(J)) 
      FLOWR = U1(I+1,J,K) * (cosXs2(J) - cosX2(J))
      FLOW = dens * (FLOWL + FLOWR) * Xs1(I+1) * Xs1(I+1) * (Xs3(K+1) - Xs3(K))
      AIP(I,J,K) = ACOF1(DIFF, -FLOW)

    ENDIF

  END SUBROUTINE CofAimAipU2_d

  !------------------------------------------------------
  ATTRIBUTES (GLOBAL) SUBROUTINE CofAimAipBoundU2_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: DIFF, FLOWL, FLOWR, FLOW

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      !J = 3
      DIFF = visk * S1u2(I,3,K) / Xdif1(I)
      FLOWL = U1(I,2,K) * (cosX2(1) - cosXs2(3))
      FLOWR = U1(I,3,K) * (cosXs2(3) - cosX2(3))
      FLOW = dens * (FLOWL + FLOWR) * Xs1(I) * Xs1(I) * (Xs3(K+1) - Xs3(K))             
      AIM(I,3,K) = ACOF1(DIFF,FLOW)
                
      DIFF = visk * S1u2(I+1,3,K) / Xdif1(I+1)
      FLOWL = U1(I+1,2,K) * (cosX2(1) - cosXs2(3))
      FLOWR = U1(I+1,3,K) * (cosXs2(3) - cosX2(3))
      FLOW = dens * (FLOWL + FLOWR) * Xs1(I+1) * Xs1(I+1) * (Xs3(K+1) - Xs3(K))    
      AIP(I,3,K) = ACOF1(DIFF,-FLOW)  

      !J = N2-1
      DIFF = visk * S1u2(I,N2-1,K) / Xdif1(I)
      FLOWL = U1(I,N2-2,K) * (cosX2(N2-2) - cosXs2(N2-1))
      FLOWR = U1(I,N2-1,K) * (cosXs2(N2-1) - cosX2(N2))
      FLOW = dens * (FLOWL + FLOWR) * Xs1(I) * Xs1(I) * (Xs3(K+1) - Xs3(K))               
      AIM(I,N2-1,K) = ACOF1(DIFF,FLOW)
              
      DIFF = visk * S1u2(I+1,N2-1,K) / Xdif1(I+1)
      FLOWL = U1(I+1,N2-2,K) * (cosX2(N2-2) - cosXs2(N2-1))
      FLOWR = U1(I+1,N2-1,K) * (cosXs2(N2-1) - cosX2(N2))
      FLOW = dens * (FLOWL + FLOWR) * Xs1(I+1) * Xs1(I+1) * (Xs3(K+1) - Xs3(K))    
      AIP(I,N2-1,K) = ACOF1(DIFF, -FLOW)  

    ENDIF

  END SUBROUTINE CofAimAipBoundU2_d

  !-----------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAjmAjpU2_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: DIFF, FLOWL, FLOWR, FLOW

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (J >= 4) .AND. (J <= (N2-2))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      DIFF = visk * S2u2(I,J-1,K) / (X1(I) * (Xs2(J) - Xs2(J-1)))
      FLOW = dens * 0.5_8 * (U2(I,J,K) * S2(I,J,K) + U2(I,J-1,K) * S2(I,J-1,K))               
      AJM(I,J,K) = ACOF1(DIFF, FLOW)             

      DIFF = visk * S2u2(I,J,K) / (X1(I) * (Xs2(J+1) - Xs2(J)))
      FLOW = dens * 0.5_8 * (U2(I,J+1,K) * S2(I,J+1,K) + U2(I,J,K) * S2(I,J,K))    
      AJP(I,J,K) = ACOF1(DIFF, -FLOW)

    ENDIF

  END SUBROUTINE CofAjmAjpU2_d

  !-----------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAjmAjpBoundU2_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: DIFF, FLOWL, FLOWR, FLOW

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      !J = 3
      DIFF = visk * S2u2(I,1,K) / (X1(I) * (Xs2(3) - Xs2(2)))
      FLOW = dens * U2(I,2,K) * S2u2(I,1,K)	               
      AJM(I,3,K) = ACOF1(DIFF,FLOW)             

      DIFF = visk * S2u2(I,3,K) / (X1(I) * (Xs2(4) - Xs2(3)))
      FLOW = dens * 0.5_8 * (U2(I,4,K) * S2(I,4,K) + U2(I,3,K) * S2(I,3,K))    
      AJP(I,3,K) = ACOF1(DIFF, -FLOW)

      !J = N2-1
      DIFF = visk * S2u2(I,N2-2,K) / (X1(I) * (Xs2(N2-1) - Xs2(N2-2)))
      FLOW = dens * 0.5_8 * (U2(I,N2-1,K) * S2(I,N2-1,K) + U2(I,N2-2,K) * S2(I,N2-2,K))               
      AJM(I,N2-1,K) = ACOF1(DIFF, FLOW)             

      DIFF = visk * S2u2(I,N2,K) / (X1(I) * (Xs2(N2) - Xs2(N2-1)))
      FLOW = dens * U2(I,N2,K) * S2(I,N2,K)    
      AJP(I,N2-1,K) = ACOF1(DIFF, -FLOW)

    ENDIF

  END SUBROUTINE CofAjmAjpBoundU2_d

  !-----------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAkmAkpU2_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: DIFF, FLOWL, FLOWR, FLOW

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (J >= 4) .AND. (J <= (N2-2))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      DIFF = visk * S3u2(I,J,K) / (X1(I) * sinXs2(J) * Xdif3(K))
      FLOWL = U3(I,J-1,K) * (Xs2(J) - X2(J-1))
      FLOWR = U3(I,J,K) * (X2(J) - Xs2(J))
      FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * (Xs1(I+1)**2 - Xs1(I)**2)         
      AKM(I,J,K) = ACOF1(DIFF, FLOW)

      DIFF = visk * S3u2(I,J,K+1) / (X1(I) * sinXs2(J) * Xdif3(K+1))
      FLOWL = U3(I,J-1,K+1) * (Xs2(J) - X2(J-1))
      FLOWR = U3(I,J,K+1) * (X2(J) - Xs2(J))
      FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * (Xs1(I+1)**2 - Xs1(I)**2)                
      AKP(I,J,K) = ACOF1(DIFF, -FLOW)           

    ENDIF

  END SUBROUTINE CofAkmAkpU2_d

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAkmAkpBoundU2_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: DIFF, FLOWL, FLOWR, FLOW

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      !J = 3
      DIFF = visk * S3u2(I,3,K) / (X1(I) * sinXs2(3) * Xdif3(K))
      FLOWL = U3(I,2,K) * (Xs2(3) - X2(1))
      FLOWR = U3(I,3,K) * (X2(3) - Xs2(3))
      FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * (Xs1(I+1)**2 - Xs1(I)**2)         
      AKM(I,3,K) = ACOF1(DIFF, FLOW)

      DIFF = visk * S3u2(I,3,K+1) / (X1(I) * sinXs2(3) * Xdif3(K+1))
      FLOWL = U3(I,2,K+1) * (Xs2(3) - X2(1))
      FLOWR = U3(I,3,K+1) * (X2(3) - Xs2(3))
      FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * (Xs1(I+1)**2 - Xs1(I)**2)                
      AKP(I,3,K) = ACOF1(DIFF, -FLOW)      

      !J = N2-1
      DIFF = visk * S3u2(I,N2-1,K) / (X1(I) * sinXs2(N2-1) * Xdif3(K))
      FLOWL = U3(I,N2-2,K) * (Xs2(N2-1) - X2(N2-2))
      FLOWR = U3(I,N2-1,K) * (X2(N2) - Xs2(N2-1))
      FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * (Xs1(I+1)**2 - Xs1(I)**2)         
      AKM(I,N2-1,K) = ACOF1(DIFF, FLOW)

      DIFF = visk * S3u2(I,N2-1,K+1) / (X1(I) * sinXs2(N2-1) * Xdif3(K+1))
      FLOWL = U3(I,N2-2,K+1) * (Xs2(N2-1) - X2(N2-2))
      FLOWR = U3(I,N2-1,K+1) * (X2(N2) - Xs2(N2-1))
      FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * (Xs1(I+1)**2 - Xs1(I)**2)                
      AKP(I,N2-1,K) = ACOF1(DIFF,-FLOW)                 

    ENDIF

  END SUBROUTINE CofAkmAkpBoundU2_d

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofConApU2hatDU2_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: WW1, WW2, WW3, WW4, WW, UY1, UY2, UY, &
             WZ1, WZ2, WZ, Spv1, Spv21, Spv22, Spv23, &
             Spv24, Spv2, Spv, Scv, AP0 

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (J >= 3) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      WW1 = U3(I,J,K) * U3(I,J,K) * (X3(K) - Xs3(K)) * WWY1(J)
      WW2 = U3(I,J,K+1) * U3(I,J,K+1) * (Xs3(K+1) - X3(K)) * WWY1(J)
      WW3 = U3(I,J-1,K) * U3(I,J-1,K) * (X3(K) - Xs3(K)) * WWY2(J)
      WW4 = U3(I,J-1,K+1) * U3(I,J-1,K+1) * (Xs3(K+1) - X3(K)) * WWY2(J)

      WW = 0.5_8 * (WW1 + WW2 + WW3 + WW4) * (Xs1(I+1) * Xs1(I+1) - Xs1(I) * Xs1(I))     


      UY1 = (U1(I,J,K) - U1(I,J-1,K)) * (X1(I) - Xs1(I)) * UYY(J)
      UY2 = (U1(I+1,J,K) - U1(I+1,J-1,K)) * (Xs1(I+1) - X1(I)) * UYY(J) 

      UY = 2.0_8 * sinXs2(J) * (Xs3(K+1) - Xs3(K)) * (UY1 + UY2)


      WZ1 = (U3(I,J,K+1) - U3(I,J,K)) * WZY1(J)
      WZ2 = (U3(I,J-1,K+1) - U3(I,J-1,K)) * WZY2(J)

      WZ = 2.0_8 * (WZ2 + WZ1) * (Xs1(I+1) - Xs1(I)) / sinXs2(J)

      Scv = dens * WW + visk * (UY - WZ)

      Spv1 = visk * (Xs3(K+1) - Xs3(K)) * Spv1J(J) * (Xs1(I+1) - Xs1(I))

      Spv21 = U1(I,J,K) * Spvy1(J) * (X1(I)**2 - Xs1(I)**2)
      Spv22 = U1(I+1,J,K) * Spvy1(J) * (Xs1(I+1)**2 - X1(I)**2)
      Spv23 = U1(I,J-1,K) * Spvy2(J) * (X1(I)**2 - Xs1(I)**2)
      Spv24 = U1(I+1,J-1,K) * Spvy2(J) * (Xs1(I+1)**2 - X1(I)**2)

      Spv2 = 0.5_8 * dens * (Xs3(K+1) - Xs3(K)) * (Spv21 + Spv22 + Spv23 + Spv24)

      Spv = Spv1 + Spv2

      IF(Spv .GE. 0.0_8) THEN
        CON(I,J,K) = CON(I,J,K) + Scv
        AP(I,J,K) = AP(I,J,K) - Spv                                     
      ELSE
        CON(I,J,K) = CON(I,J,K) + Scv - Spv2 * U2(I,J,K)
        AP(I,J,K) = AP(I,J,K) - Spv1                               
      END IF

      AP0 = dens * VCV2(I,J,K) / DT
      CON(I,J,K) = CON(I,J,K) + AP0 * U2_0(I,J,K)
      AP(I,J,K) = - AP(I,J,K) + AP0 + &
                    AIM(I,J,K) + AIP(I,J,K) + &
                    AJM(I,J,K) + AJP(I,J,K)+ &
                    AKM(I,J,K) + AKP(I,J,K)
      DU2(I,J,K) = 0.5_8 * (S2u2(I,J,K) + S2u2(I,J-1,K)) / AP(I,J,K)

      U2hat(I,J,K) = (AIP(I,J,K) * U2(I+1,J,K) + AIM(I,J,K) * U2(I-1,J,K) + & 
                      AJP(I,J,K) * U2(I,J+1,K) + AJM(I,J,K) * U2(I,J-1,K)+ &
                      AKP(I,J,K) * U2(I,J,K+1) + AKM(I,J,K) * U2(I,J,K-1)+ &
                      CON(I,J,K)) / AP(I,J,K)

    ENDIF

  END SUBROUTINE CofConApU2hatDU2_d

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofBoundDU2_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      !J = 3
      DU2(I,3,K) = 0.5_8 * (Xs3(K+1) - Xs3(K)) * (cosX2(1) - cosX2(3)) * &
                           (Xs1(I+1)**2 - Xs1(I)**2) / (AP(I,3,K) * Xdif2(3))

      !J = N2-1   
      DU2(I,N2-1,K) = 0.5_8 * (Xs3(K+1) - Xs3(K)) * (cosX2(N2-2) - cosX2(N2)) * &
                              (Xs1(I+1)**2 - Xs1(I)**2) / (AP(I,N2-1,K) * Xdif2(N2-1))              

    ENDIF

  END SUBROUTINE CofBoundDU2_d

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofConU2_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (J >= 3) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      CON(I,J,K) = CON(I,J,K) + DU2(I,J,K) * AP(I,J,K) * (P(I,J-1,K) - P(I,J,K))    

    ENDIF

  END SUBROUTINE CofConU2_d

    
  !-----------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofU3_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: DIFF, FLOW, FLOWL, FLOWR, AP0
  REAL(8) :: UZ, VZ1, VZ2, VZ, Scw, &
             Spw1, Spw21, Spw22, Spw2, &
             Spw31, Spw32, Spw3, Spw

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (J >= 2) .AND. (J <= (N2-2))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      !-------------- AIM, AIP ---------------
      DIFF = visk * S1u3(I,J,K) / Xdif1(I)	  
      FLOWL = U1(I,J,K-1) * 0.5_8 * X3cv(K-1)
      FLOWR = U1(I,J,K) * 0.5_8 * X3cv(K)
      FLOW = dens * (FLOWL + FLOWR) * Xs1(I) * Xs1(I) * (cosXs2(J) - cosXs2(J+1))               
      AIM(I,J,K) = ACOF1(DIFF, FLOW)            

      DIFF = visk * S1u3(I+1,J,K) / Xdif1(I+1)
      FLOWL = U1(I+1,J,K-1) * 0.5_8 * X3cv(K-1)
      FLOWR = U1(I+1,J,K) * 0.5_8 * X3cv(K)
      FLOW = dens * (FLOWL + FLOWR) * Xs1(I+1) * Xs1(I+1) * (cosXs2(J) - cosXs2(J+1))   
      AIP(I,J,K) = ACOF1(DIFF, -FLOW) 

      !-------------- AJM, AJP --------------- 
      DIFF = visk * S2u3(I,J,K) / (X1(I) * Xdif2(J))
      FLOWL = U2(I,J,K-1) * 0.5_8 * X3cv(K-1) 
      FLOWR = U2(I,J,K) * 0.5_8 * X3cv(K)
      FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * sinXs2(J) * (Xs1(I+1)**2 - Xs1(I)**2)
      AJM(I,J,K) = ACOF1(DIFF, FLOW)            

      DIFF = visk * S2u3(I,J+1,K) / (X1(I) * Xdif2(J+1))
      FLOWL = U2(I,J+1,K-1) * 0.5_8 * X3cv(K-1) 
      FLOWR = U2(I,J+1,K) * 0.5_8 * X3cv(K)
      FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * sinXs2(J+1) * (Xs1(I+1)**2 - Xs1(I)**2)	               
      AJP(I,J,K) = ACOF1(DIFF, -FLOW)


      !-------------- AKM, AKP ---------------
      DIFF = visk * S3(I,J,K) / (X1(I) * sinX2(J) * X3cv(K-1))
      FLOW = dens * 0.5_8 * (U3(I,J,K) + U3(I,J,K-1)) * S3(I,J,K)        
      AKM(I,J,K) = ACOF1(DIFF, FLOW)

      DIFF = visk * S3(I,J,K+1) / (X1(I) * sinX2(J) * X3cv(K))
      FLOW = dens * 0.5_8 * (U3(I,J,K+1) + U3(I,J,K)) * S3(I,J,K+1)
      AKP(I,J,K) = ACOF1(DIFF, -FLOW) 

      !---------------- Sc, Sp ------------------

      UZ = 2.0_8 * ( (U1(I,J,K) - U1(I,J,K-1)) * (X1(I) - Xs1(I)) + &
                   (U1(I+1,J,K) - U1(I+1,J,K-1)) * (Xs1(I+1) - X1(I)) ) * X2cv(J)

      VZ1 = (U2(I,J,K) - U2(I,J,K-1)) * VZY1(J)
      VZ2 = (U2(I,J+1,K) - U2(I,J+1,K-1)) * VZY2(J)

      VZ = 2.0_8 * X1cv(I) * (VZ1 + VZ2) / sinX2(J)

      Scw = visk * (UZ + VZ)

      Spw1 = visk * X1cv(I) * X2cv(J) * Xdif3(K) / sinX2(J)

      Spw21 = (U1(I,J,K) * 0.5_8 * X3cv(K) + U1(I,J,K-1) * 0.5_8 * X3cv(K-1)) * &
              (X1(I)**2 - Xs1(I)**2)
      Spw22 = (U1(I+1,J,K) * 0.5_8 * X3cv(K) + U1(I+1,J,K-1) * 0.5_8 * X3cv(K-1)) * &
              (Xs1(I+1)**2 - X1(I)**2)
      Spw2 = 0.5_8 * dens * (cosXs2(J) - cosXs2(J+1)) * (Spw21 + Spw22)

      Spw31 = (U2(I,J,K) * Spwy1(J) + U2(I,J+1,K) * Spwy2(J)) * 0.5_8 * X3cv(K)
      Spw32 = (U2(I,J,K-1) * Spwy1(J) + U2(I,J+1,K-1) * Spwy2(J)) * 0.5_8 * X3cv(K-1)

      Spw3 = 0.5_8 * dens * (Xs1(I+1)**2 - Xs1(I)**2) * (Spw31 + Spw32)

      Spw = Spw1 + Spw2 + Spw3

      IF(Spw .GE. 0.0_8) THEN
        CON(I,J,K) = CON(I,J,K) + Scw
        AP(I,J,K) = AP(I,J,K) - Spw     
      ELSE 
        IF((Spw2 .LT. 0.0_8) .AND. (Spw3 .GE. 0.0_8)) THEN           
          CON(I,J,K) = CON(I,J,K) + Scw - Spw2 * U3(I,J,K)
          AP(I,J,K) = AP(I,J,K) - Spw1 - Spw3
        END IF

        IF((Spw2 .GE. 0.0_8) .AND. (Spw3 .LT. 0.0_8)) THEN           
          CON(I,J,K) = CON(I,J,K) + Scw - Spw3 * U3(I,J,K)
          AP(I,J,K) = AP(I,J,K) - Spw1 - Spw2
        END IF

        IF((Spw2 .LT. 0.0_8) .AND. (Spw3 .LT. 0.0_8)) THEN           
          CON(I,J,K) = CON(I,J,K) + Scw - (Spw2 + Spw3) * U3(I,J,K)
          AP(I,J,K) = AP(I,J,K) - Spw1
        END IF
      END IF

      AP0 = dens * VCV3(I,J,K) / DT
      CON(I,J,K) = CON(I,J,K) + AP0 * U3_0(I,J,K)
      AP(I,J,K) = - AP(I,J,K) + AP0 + &
                    AIM(I,J,K) + AIP(I,J,K) + &
                    AJM(I,J,K) + AJP(I,J,K) + &
                    AKM(I,J,K) + AKP(I,J,K)
      DU3(I,J,K) = S3(I,J,K) / AP(I,J,K)


      U3hat(I,J,K) = (AIP(I,J,K) * U3(I+1,J,K) + AIM(I,J,K) * U3(I-1,J,K) + &
                      AJP(I,J,K) * U3(I,J+1,K) + AJM(I,J,K) * U3(I,J-1,K) + &
                      AKP(I,J,K) * U3(I,J,K+1) + AKM(I,J,K) * U3(I,J,K-1) + &
                      CON(I,J,K)) / AP(I,J,K)

      CON(I,J,K) = CON(I,J,K) + S3(I,J,K) * (P(I,J,K-1) - P(I,J,K))
 
    ENDIF

  END SUBROUTINE CofU3_d

  !------------------------------------------------------
  
  ATTRIBUTES (GLOBAL) SUBROUTINE CofAjpOnBotU3_d
  USE VAR_d
  USE AFUN_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: DIFF, FLOW, FLOWL, FLOWR, AP0
  REAL(8) :: UZ, VZ1, VZ2, VZ, Scw, &
             Spw1, Spw21, Spw22, Spw2, &
             Spw31, Spw32, Spw3, Spw

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
   
    J = N2-1

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      !-------------- AIM, AIP ---------------
      DIFF = visk * S1u3(I,J,K) / Xdif1(I)	  
      FLOWL = U1(I,J,K-1) * 0.5_8 * X3cv(K-1)
      FLOWR = U1(I,J,K) * 0.5_8 * X3cv(K)
      FLOW = dens * (FLOWL + FLOWR) * Xs1(I) * Xs1(I) * (cosXs2(J) - cosXs2(J+1))               
      AIM(I,J,K) = ACOF1(DIFF, FLOW)            

      DIFF = visk * S1u3(I+1,J,K) / Xdif1(I+1)
      FLOWL = U1(I+1,J,K-1) * 0.5_8 * X3cv(K-1)
      FLOWR = U1(I+1,J,K) * 0.5_8 * X3cv(K)
      FLOW = dens * (FLOWL + FLOWR) * Xs1(I+1) * Xs1(I+1) * (cosXs2(J) - cosXs2(J+1))   
      AIP(I,J,K) = ACOF1(DIFF, -FLOW) 

      !-------------- AJM, AJP --------------- 
      DIFF = visk * S2u3(I,J,K) / (X1(I) * Xdif2(J))
      FLOWL = U2(I,J,K-1) * 0.5_8 * X3cv(K-1) 
      FLOWR = U2(I,J,K) * 0.5_8 * X3cv(K)
      FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * sinXs2(J) * (Xs1(I+1)**2 - Xs1(I)**2)
      AJM(I,J,K) = ACOF1(DIFF, FLOW)            

      !DIFF = visk * S2u3(I,J+1,K) / (X1(I) * Xdif2(J+1))
      !FLOWL = U2(I,J+1,K-1) * 0.5_8 * X3cv(K-1) 
      !FLOWR = U2(I,J+1,K) * 0.5_8 * X3cv(K)
      !FLOW = dens * (FLOWL + FLOWR) * 0.5_8 * sinXs2(J+1) * (Xs1(I+1)**2 - Xs1(I)**2)	               
      AJP(I,J,K) = 0.0_8!ACOF1(DIFF, -FLOW)


      !-------------- AKM, AKP ---------------
      DIFF = visk * S3(I,J,K) / (X1(I) * sinX2(J) * X3cv(K-1))
      FLOW = dens * 0.5_8 * (U3(I,J,K) + U3(I,J,K-1)) * S3(I,J,K)        
      AKM(I,J,K) = ACOF1(DIFF, FLOW)

      DIFF = visk * S3(I,J,K+1) / (X1(I) * sinX2(J) * X3cv(K))
      FLOW = dens * 0.5_8 * (U3(I,J,K+1) + U3(I,J,K)) * S3(I,J,K+1)
      AKP(I,J,K) = ACOF1(DIFF, -FLOW) 

      !---------------- Sc, Sp ------------------

      UZ = 2.0_8 * ( (U1(I,J,K) - U1(I,J,K-1)) * (X1(I) - Xs1(I)) + &
                   (U1(I+1,J,K) - U1(I+1,J,K-1)) * (Xs1(I+1) - X1(I)) ) * X2cv(J)

      VZ1 = (U2(I,J,K) - U2(I,J,K-1)) * VZY1(J)
      VZ2 = (U2(I,J+1,K) - U2(I,J+1,K-1)) * VZY2(J)

      VZ = 2.0_8 * X1cv(I) * (VZ1 + VZ2) / sinX2(J)

      Scw = visk * (UZ + VZ)

      Spw1 = visk * X1cv(I) * X2cv(J) * Xdif3(K) / sinX2(J)

      Spw21 = (U1(I,J,K) * 0.5_8 * X3cv(K) + U1(I,J,K-1) * 0.5_8 * X3cv(K-1)) * &
              (X1(I)**2 - Xs1(I)**2)
      Spw22 = (U1(I+1,J,K) * 0.5_8 * X3cv(K) + U1(I+1,J,K-1) * 0.5_8 * X3cv(K-1)) * &
              (Xs1(I+1)**2 - X1(I)**2)
      Spw2 = 0.5_8 * dens * (cosXs2(J) - cosXs2(J+1)) * (Spw21 + Spw22)

      Spw31 = (U2(I,J,K) * Spwy1(J) + U2(I,J+1,K) * Spwy2(J)) * 0.5_8 * X3cv(K)
      Spw32 = (U2(I,J,K-1) * Spwy1(J) + U2(I,J+1,K-1) * Spwy2(J)) * 0.5_8 * X3cv(K-1)

      Spw3 = 0.5_8 * dens * (Xs1(I+1)**2 - Xs1(I)**2) * (Spw31 + Spw32)

      Spw = Spw1 + Spw2 + Spw3

      IF(Spw .GE. 0.0_8) THEN
        CON(I,J,K) = CON(I,J,K) + Scw
        AP(I,J,K) = AP(I,J,K) - Spw     
      ELSE 
        IF((Spw2 .LT. 0.0_8) .AND. (Spw3 .GE. 0.0_8)) THEN           
          CON(I,J,K) = CON(I,J,K) + Scw - Spw2 * U3(I,J,K)
          AP(I,J,K) = AP(I,J,K) - Spw1 - Spw3
        END IF

        IF((Spw2 .GE. 0.0_8) .AND. (Spw3 .LT. 0.0_8)) THEN           
          CON(I,J,K) = CON(I,J,K) + Scw - Spw3 * U3(I,J,K)
          AP(I,J,K) = AP(I,J,K) - Spw1 - Spw2
        END IF

        IF((Spw2 .LT. 0.0_8) .AND. (Spw3 .LT. 0.0_8)) THEN           
          CON(I,J,K) = CON(I,J,K) + Scw - (Spw2 + Spw3) * U3(I,J,K)
          AP(I,J,K) = AP(I,J,K) - Spw1
        END IF
      END IF

      AP0 = dens * VCV3(I,J,K) / DT
      CON(I,J,K) = CON(I,J,K) + AP0 * U3_0(I,J,K)
      AP(I,J,K) = - AP(I,J,K) + AP0 + &
                    AIM(I,J,K) + AIP(I,J,K) + &
                    AJM(I,J,K) + AJP(I,J,K) + &
                    AKM(I,J,K) + AKP(I,J,K)
      DU3(I,J,K) = S3(I,J,K) / AP(I,J,K)


      U3hat(I,J,K) = (AIP(I,J,K) * U3(I+1,J,K) + AIM(I,J,K) * U3(I-1,J,K) + &
                      AJP(I,J,K) * U3(I,J+1,K) + AJM(I,J,K) * U3(I,J-1,K) + &
                      AKP(I,J,K) * U3(I,J,K+1) + AKM(I,J,K) * U3(I,J,K-1) + &
                      CON(I,J,K)) / AP(I,J,K)

      CON(I,J,K) = CON(I,J,K) + S3(I,J,K) * (P(I,J,K-1) - P(I,J,K))
 
    ENDIF

  END SUBROUTINE CofAjpOnBotU3_d

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAimAipConP_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= (N1_core + 1)) .AND. (I <= (N1-2))  .AND. &
         (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN
            
      AIM(I,J,K) = S1(I,J,K) * DU1(I,J,K)
      AIP(I,J,K) = S1(I+1,J,K) * DU1(I+1,J,K)
      CON(I,J,K) = S1(I,J,K) * U1hat(I,J,K) - S1(I+1,J,K) * U1hat(I+1,J,K)

    ENDIF

  END SUBROUTINE CofAimAipConP_d

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAimAipConBoundP_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN
            
        AIM(N1_core,J,K) = 0.0_8
        AIP(N1_core,J,K) = S1(N1_core + 1,J,K) * DU1(N1_core + 1,J,K)

        AIM(N1-1,J,K) = S1(N1-1,J,K) * DU1(N1-1,J,K)
        AIP(N1-1,J,K) = 0.0_8    

        CON(N1_core,J,K) = S1(N1_core,J,K) * U1(N1_core,J,K) - S1(N1_core + 1,J,K) * U1hat(N1_core + 1,J,K)
        CON(N1-1,J,K) = S1(N1-1,J,K) * U1hat(N1-1,J,K) - S1(N1,J,K) * U1(N1,J,K)

    ENDIF

  END SUBROUTINE CofAimAipConBoundP_d

  !----------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAjmAjpConP_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (J >= 3) .AND. (J <= (N2-2))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN
            
      AJM(I,J,K) = S2(I,J,K) * DU2(I,J,K)
      AJP(I,J,K) = S2(I,J+1,K) * DU2(I,J+1,K)
      CON(I,J,K) = CON(I,J,K) + S2(I,J,K) * U2hat(I,J,K) - S2(I,J+1,K) * U2hat(I,J+1,K)

    ENDIF

  END SUBROUTINE CofAjmAjpConP_d

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAjmAjpConBoundP_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN
            
      AJM(I,2,K) = 0.0_8
      AJP(I,2,K) = S2(I,3,K) * DU2(I,3,K)
      CON(I,2,K) = CON(I,2,K) + S2(I,2,K) * U2(I,2,K) - S2(I,3,K) * U2hat(I,3,K)

      AJM(I,N2-1,K) = S2(I,N2-1,K) * DU2(I,N2-1,K)
      AJP(I,N2-1,K) = 0.0_8
      CON(I,N2-1,K) = CON(I,N2-1,K) + S2(I,N2-1,K) * U2hat(I,N2-1,K) - S2(I,J+1,K) * U2(I,N2,K)

    ENDIF

  END SUBROUTINE CofAjmAjpConBoundP_d

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAkmAkpConP_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN
            
      AKM(I,J,K) = S3(I,J,K) * DU3(I,J,K)
      AKP(I,J,K) = S3(I,J,K+1) * DU3(I,J,K+1)
      CON(I,J,K) = CON(I,J,K) + S3(I,J,K) * U3hat(I,J,K) - S3(I,J,K+1) * U3hat(I,J,K+1)

    ENDIF

  END SUBROUTINE CofAkmAkpConP_d

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofApP_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      AP(I,J,K) = AIM(I,J,K) + AIP(I,J,K) + AJM(I,J,K) + AJP(I,J,K) + AKM(I,J,K) + AKP(I,J,K)
            
    ENDIF

  END SUBROUTINE CofApP_d


  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAimAipPC_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= (N1_core + 1)) .AND. (I <= (N1-2))  .AND. &
         (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN
            
      AIM(I,J,K) = S1(I,J,K) * DU1(I,J,K)
      AIP(I,J,K) = S1(I+1,J,K) * DU1(I+1,J,K)

    ENDIF

  END SUBROUTINE CofAimAipPC_d

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAimAipBoundPC_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN
            
        AIM(N1_core,J,K) = 0.0_8
        AIP(N1_core,J,K) = S1(N1_core + 1,J,K) * DU1(N1_core + 1,J,K)

        AIM(N1-1,J,K) = S1(N1-1,J,K) * DU1(N1-1,J,K)
        AIP(N1-1,J,K) = 0.0_8    

    ENDIF

  END SUBROUTINE CofAimAipBoundPC_d

  !----------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAjmAjpPC_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (J >= 3) .AND. (J <= (N2-2))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN
            
      AJM(I,J,K) = S2(I,J,K) * DU2(I,J,K)
      AJP(I,J,K) = S2(I,J+1,K) * DU2(I,J+1,K)

    ENDIF

  END SUBROUTINE CofAjmAjpPC_d

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAjmAjpBoundPC_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN
            
      AJM(I,2,K) = 0.0_8
      AJP(I,2,K) = S2(I,3,K) * DU2(I,3,K)

      AJM(I,N2-1,K) = S2(I,N2-1,K) * DU2(I,N2-1,K)
      AJP(I,N2-1,K) = 0.0_8

    ENDIF

  END SUBROUTINE CofAjmAjpBoundPC_d

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofAkmAkpPC_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN
            
      AKM(I,J,K) = S3(I,J,K) * DU3(I,J,K)
      AKP(I,J,K) = S3(I,J,K+1) * DU3(I,J,K+1)

    ENDIF

  END SUBROUTINE CofAkmAkpPC_d

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CofApConPC_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      AP(I,J,K) = AIM(I,J,K) + AIP(I,J,K) + AJM(I,J,K) + AJP(I,J,K) + AKM(I,J,K) + AKP(I,J,K)
      CON(I,J,K) = S1(I,J,K) * U1tmp(I,J,K) - S1(I+1,J,K) * U1tmp(I+1,J,K) + &
                   S2(I,J,K) * U2tmp(I,J,K) - S2(I,J+1,K) * U2tmp(I,J+1,K) + &
                   S3(I,J,K) * U3tmp(I,J,K) - S3(I,J,K+1) * U3tmp(I,J,K+1)
            
    ENDIF

  END SUBROUTINE CofApConPC_d

  !------------------------------------------------------

  SUBROUTINE ComputeSmax
  USE VAR_d
  USE VAR, ONLY: smax, N1_core
  IMPLICIT NONE
  INTEGER I, J, K

    smax = 0.0_8

   !$cuf kernel do(3) <<<*, *>>>
    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = N1_core, N1-1
          smax = max(smax, DABS(S1(I+1,J,K) * U1(I+1,J,K) - S1(I,J,K) * U1(I,J,K) + &
                                S2(I,J+1,K) * U2(I,J+1,K) - S2(I,J,K) * U2(I,J,K) + &
                                S3(I,J,K+1) * U3(I,J,K+1) - S3(I,J,K) * U3(I,J,K) ) / VCV(I,J,K))
        ENDDO          
      ENDDO
    ENDDO

  END SUBROUTINE ComputeSmax

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CorrectU1_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= (N1_core + 1)) .AND. (I <= (N1-1))  .AND. &
         (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      U1(I,J,K) =  U1(I,J,K) + DU1(I,J,K) * (PC(I-1,J,K) - PC(I,J,K))
            
    ENDIF


  END SUBROUTINE CorrectU1_d

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CorrectU2_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (J >= 3) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      U2(I,J,K) = U2(I,J,K) + DU2(I,J,K) * (PC(I,J-1,K) - PC(I,J,K))
            
    ENDIF

  END SUBROUTINE CorrectU2_d

  !------------------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CorrectU3_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= N1_core) .AND. (I <= (N1-1))  .AND. &
         (J >= 2) .AND. (J <= (N2-1))  .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      U3(I,J,K) = U3(I,J,K) + DU3(I,J,K) * (PC(I,J,K-1) - PC(I,J,K))
            
    ENDIF

  END SUBROUTINE CorrectU3_d

 !------------------------------------------------------

END MODULE FD
