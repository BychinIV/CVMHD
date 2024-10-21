MODULE MD
IMPLICIT NONE
  CONTAINS
  
  !---------------------------------------------
  SUBROUTINE MHD
    USE VAR, ONLY: NTIMES_F1, NTIMES_F2, NTIMES_F3, SMAX_B, &
                   resF1, resF2, resF3 
    USE VAR_d
    USE SOLVER, ONLY: AdiSolverF1, AdiSolverF2, AdiSolverF3, FF1N1_d
    IMPLICIT NONE
    INTEGER :: I, J, K, IN_ITER, OUT_ITER, ITER_F3t, ITER_F2t, ITER_F1t
    
    !DO OUT_ITER = 1, 2  

      !CALL ComputeE1_Oz 

      F1t = F1; F2t = F2; F3t = F3
      resF1 = 1.0_8; resF2 = 1.0_8; resF3 = 1.0_8

      !----------Predictor--------
      DO WHILE ( (resF1 + resF2 + resF3) > 1.0d-5 )

        CALL ComputeCofF3
        CALL AdiSolverF3(NTIMES_F3) 
 
        CALL ComputeCofF2
        CALL AdiSolverF2(NTIMES_F2)

        CALL ComputeCofF1
        CALL AdiSolverF1(NTIMES_F1)


        !CALL FtoB
        !CALL ComputeSmaxB   

        !WRITE(*, *) 'SMAX_B  = ', SMAX_B   
        
        CALL ComputeCofF3
        CALL ResidualF3

        CALL ComputeCofF2
        CALL ResidualF2

        CALL ComputeCofF1
        CALL ResidualF1

      ENDDO


      !-----Corrector-------
      CALL ComputeCofF3
      CALL ComputeF3t

      CALL ComputeCofF2
      CALL ComputeF2t

      CALL ComputeCofF1
      CALL ComputeF1t


      F1 = F1t; F2 = F2t; F3 = F3t
      CALL FtoB
      CALL ComputeSmaxB   

      !WRITE(*, *) 'SMAX_B by Ft  = ', SMAX_B


      !CALL ComputeCofF3
      !CALL ResidualF3

      !CALL ComputeCofF2
      !CALL ResidualF2

      !CALL ComputeCofF1
      !CALL ResidualF1   


      CALL ComputeMagPotential
      CALL ComputeBoundCondsF


  !ENDDO !OUT_ITER  
 
  END SUBROUTINE MHD
  
  !---------------------------------------------
  
  SUBROUTINE ComputeMagPotential
    USE VAR, ONLY: FirstCallMP, tolPSI
    USE SOLVER, ONLY: AdiSolverPSI_OUT, AdiSolverPSI_IN
    IMPLICIT NONE
    INTEGER :: numIter

    IF (FirstCallMP .EQ. .TRUE.) THEN
      !CALL ComputeCofPSI_IN
      CALL ComputeCofPSI_OUT
      FirstCallMP = .FALSE.
    ELSE
      !CALL UpdateSourcePSI_IN
      CALL UpdateSourcePSI_OUT
    ENDIF

    !numIter = AdiSolverPSI_IN(tolPSI)
    !WRITE(*, '(A, I8)') 'Solve PSI_IN. NumIters:', numIter

    numIter = AdiSolverPSI_OUT(tolPSI)
    WRITE(*, '(A, I8)') 'Solve PSI_OUT. NumIters:', numIter    
 
  END SUBROUTINE ComputeMagPotential
    
  !---------------------------------------------
  
  
  SUBROUTINE ComputeBoundCondsF
    USE VAR_d
    USE SOLVER, ONLY : FF1N1_d
    IMPLICIT NONE

    CALL ComputeF2OutBound_d <<<GRD_JK, TBL_JK>>>
    CALL ComputeF3OutBound_d <<<GRD_JK, TBL_JK>>>
    !CALL ComputeF2InBound_d <<<GRD_JK, TBL_JK>>>
    !CALL ComputeF3InBound_d <<<GRD_JK, TBL_JK>>>
    
    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (F2, N1, N1)
    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (F3, N1, N1)

  END SUBROUTINE ComputeBoundCondsF
  
  !---------------------------------------------

  SUBROUTINE ComputeCofF1
    USE VAR_d
    USE SOLVER, ONLY: F1ext_d, FF1N1_d
    IMPLICIT NONE

    !CALL F1ext_d <<<GRD_IJ, TBL_IJ>>>
    CALL ComputeF1_Oz_d <<<GRD_I, TBL_I>>>
    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (F1, N1, N1)
    CALL ComputeCofF1_d <<<GRD_IJK, TBL_IJK>>>

  END SUBROUTINE ComputeCofF1

  !---------------------------------------------

  SUBROUTINE ComputeF1t
    USE VAR_d
    USE SOLVER, ONLY: FF1N1_d
    IMPLICIT NONE

    CALL ComputeF1t_d <<<GRD_IJK, TBL_IJK>>>
    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (F1t, N1, N1)
  END SUBROUTINE ComputeF1t

  !---------------------------------------------

  SUBROUTINE ComputeCofF2
    USE VAR_d
    USE SOLVER, ONLY: F2ext_d, FF1N1_d
    IMPLICIT NONE

    !CALL F2ext_d <<<GRD_IJ, TBL_IJ>>>
    CALL ComputeF1_Oz_d <<<GRD_I, TBL_I>>>
    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (F1, N1, N1)
    CALL ComputeCofF2_d <<<GRD_IJK, TBL_IJK>>>

  END SUBROUTINE ComputeCofF2

  !---------------------------------------------

  SUBROUTINE ComputeF2t
    USE VAR_d
    USE SOLVER, ONLY: FF1N1_d
    IMPLICIT NONE

    CALL ComputeF2t_d <<<GRD_IJK, TBL_IJK>>>
    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (F2t, N1, N1)
  END SUBROUTINE ComputeF2t

  !---------------------------------------------

  SUBROUTINE ComputeF3t
    USE VAR_d
    USE SOLVER, ONLY: FF1N1_d
    IMPLICIT NONE

    CALL ComputeF3t_d <<<GRD_IJK, TBL_IJK>>>
    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (F3t, N1, N1)
  END SUBROUTINE ComputeF3t

  !---------------------------------------------


  SUBROUTINE ResidualF1
    USE VAR_d
    USE VAR, ONLY: resF1
    IMPLICIT NONE
    INTEGER :: I, J, K
    REAL(8) :: sum_sq_f_h

    sum_sq_f = 0.0_8

    !$cuf kernel do (3) <<<*, *>>>
    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 3, N1
          sum_sq_f = sum_sq_f +                                        &
              ( AP(I,J,K)*F1(I,J,K) - ( AJP(I,J,K)*F1(I,J+1,K) + AJM(I,J,K)*F1(I,J-1,K) + &
                                        AKP(I,J,K)*F1(I,J,K+1) + AKM(I,J,K)*F1(I,J,K-1) + &
                                        CON(I,J,K) ) )**2
        ENDDO
      ENDDO
    ENDDO

    sum_sq_f_h = sum_sq_f
    resF1 = DSQRT(sum_sq_f_h)
    !WRITE(*,*) 'ResiduaF1 = ', resF1

  END SUBROUTINE ResidualF1

  !---------------------------------------------

  SUBROUTINE ResidualF2
    USE VAR_d
    USE VAR, ONLY: resF2
    IMPLICIT NONE
    INTEGER :: I, J, K
    REAL(8) :: sum_sq_f_h

    sum_sq_f = 0.0_8

    !$cuf kernel do (3) <<<*, *>>>
    DO K = 2, N3-1
      DO J = 3, N2-1
        DO I = 2, N1-1
          sum_sq_f = sum_sq_f +                                        &
              ( AP(I,J,K)*F2(I,J,K) - ( AIP(I,J,K)*F2(I+1,J,K) + AIM(I,J,K)*F2(I-1,J,K) + &
                                        AKP(I,J,K)*F2(I,J,K+1) + AKM(I,J,K)*F2(I,J,K-1) + &
                                        CON(I,J,K) ) )**2
        ENDDO
      ENDDO
    ENDDO

    sum_sq_f_h = sum_sq_f
    resF2 = DSQRT(sum_sq_f_h)
    !WRITE(*,*) 'ResiduaF2 = ', resF2

  END SUBROUTINE ResidualF2

  !---------------------------------------------

  SUBROUTINE ResidualF3
    USE VAR_d
    USE VAR, ONLY: resF3
    IMPLICIT NONE
    INTEGER :: I, J, K
    REAL(8) :: sum_sq_f_h

    sum_sq_f = 0.0_8

    !$cuf kernel do (3) <<<*, *>>>
    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 2, N1-1
          sum_sq_f = sum_sq_f +                                        &
              ( AP(I,J,K)*F3(I,J,K) - ( AIP(I,J,K)*F3(I+1,J,K) + AIM(I,J,K)*F3(I-1,J,K) + &
                                        AJP(I,J,K)*F3(I,J+1,K) + AJM(I,J,K)*F3(I,J-1,K) + &
                                        CON(I,J,K) ) )**2
        ENDDO
      ENDDO
    ENDDO

    sum_sq_f_h = sum_sq_f
    resF3 = DSQRT(sum_sq_f_h)
    !WRITE(*,*) 'ResiduaF3 = ', resF3

  END SUBROUTINE ResidualF3

  !---------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeF1t_d
    USE VAR_d
    IMPLICIT NONE
    INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 3 .AND. I <= N1) .AND.     &
         (J >= 2 .AND. J <= (N2-1)) .AND. &
         (K >= 2 .AND. K <= (N3-1)) ) THEN

      F1t(I,J,K) = DT * (AJM(I,J,K) * (F1(I,J-1,K) - F1(I,J,K)) + AJP(I,J,K) * (F1(I,J+1,K) - F1(I,J,K)) + &
                         AKM(I,J,K) * (F1(I,J,K-1) - F1(I,J,K)) + AKP(I,J,K) * (F1(I,J,K+1) - F1(I,J,K)) + CON(I,J,K)) / S1(I,J,K) 

    ENDIF

  END SUBROUTINE ComputeF1t_d

  !---------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeF2t_d
    USE VAR_d
    IMPLICIT NONE
    INTEGER I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 2 .AND. I <= (N1-1)) .AND. &
         (J >= 3 .AND. J <= N2) .AND. &
         (K >= 2 .AND. K <= (N3-1)) ) THEN


      F2t(I,J,K) = X1(I) * DT * ( AIM(I,J,K) * (F2(I-1,J,K) - F2(I,J,K)) + &
                                  AIP(I,J,K) * (F2(I+1,J,K) - F2(I,J,K)) + &
                                  AKM(I,J,K) * (F2(I,J,K-1) - F2(I,J,K)) + &
                                  AKP(I,J,K) * (F2(I,J,K+1) - F2(I,J,K)) + CON(I,J,K)) / S2(I,J,K)

    ENDIF

  END SUBROUTINE ComputeF2t_d


  !---------------------------------------------


  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeF3t_d
    USE VAR_d
    IMPLICIT NONE
    INTEGER I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 2 .AND. I <= (N1-1)) .AND. &
         (J >= 2 .AND. J <= (N2-1)) .AND. &
         (K >= 2 .AND. K <= (N3-1)) ) THEN


      F3t(I,J,K) = X1(I) * sinX2(J) * DT * (AIM(I,J,K) * (F3(I-1,J,K) - F3(I,J,K)) + &
                                            AIP(I,J,K) * (F3(I+1,J,K) - F3(I,J,K)) + &
                                            AJM(I,J,K) * (F3(I,J-1,K) - F3(I,J,K)) + &
                                            AJP(I,J,K) * (F3(I,J+1,K) - F3(I,J,K)) + CON(I,J,K)) / S3(I,J,K)

    ENDIF

  END SUBROUTINE ComputeF3t_d


  !---------------------------------------------


  SUBROUTINE ComputeCofF3
    USE VAR_d
    USE SOLVER, ONLY: F1ext_d, F2ext_d
    IMPLICIT NONE

    CALL F1ext_d <<<GRD_IJ, TBL_IJ>>>
    CALL F2ext_d <<<GRD_IJ, TBL_IJ>>>
    CALL ComputeCofF3_d <<<GRD_IJK, TBL_IJK>>>
    CALL ComputeE1_Oz_d <<<GRD_I, TBL_I>>>
    CALL ComputeCofF3_Oz_d <<<GRD_IK, TBL_IK>>>

  END SUBROUTINE ComputeCofF3

  !---------------------------------------------

  SUBROUTINE ComputeE1_Oz
    USE VAR_d
    IMPLICIT NONE

    CALL ComputeE1_Oz_d <<<GRD_I, TBL_I>>>

  END SUBROUTINE ComputeE1_Oz

  !---------------------------------------------

  SUBROUTINE CorrectB
    USE VAR, ONLY: tolPSI_C, SMAX_B
    USE VAR_d
    USE SOLVER, ONLY: AdiSolverPSI_C, FF1N1_d
    IMPLICIT NONE
    INTEGER :: numIter

    CALL ComputeCofPSI_C(N1)
    numIter = AdiSolverPSI_C(N1, tolPSI_C)

    WRITE(*, '(A, I8)') 'Solve PSI_C. NumIters:', numIter
    WRITE(10, '(A, I8)') 'Solve PSI_C. NumIters:', numIter


    CALL CorrectF1_d <<<GRD_IJK, TBL_IJK>>> (N1)
    CALL CorrectF2_d <<<GRD_IJK, TBL_IJK>>> (N1)
    CALL CorrectF3_d <<<GRD_IJK, TBL_IJK>>> (N1)

    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (F1, N1, N1)
    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (F2, N1, N1)
    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (F3, N1, N1)

    CALL FtoB

    CALL ComputeSmaxB


    WRITE(*, *) 'SMAX_B AFTER CORRECT = ', SMAX_B
    WRITE(10, *) 'SMAX_B AFTER CORRECT = ', SMAX_B

  END SUBROUTINE CorrectB

  !-------------------------------------------

  SUBROUTINE ComputeRotB
    USE VAR_d
    USE SOLVER, ONLY: FF1N1_d
    IMPLICIT NONE

    CALL ComputeJ1_d <<<GRD_IJK, TBL_IJK>>>
    CALL ComputeJ1_Oz_d <<<GRD_I, TBL_I>>>    
    CALL ComputeJ2_d <<<GRD_IJK, TBL_IJK>>>    
    CALL ComputeJ3_d <<<GRD_IJK, TBL_IJK>>>

    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (J1, N1, N1)
    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (J2, N1, N1)
    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (J3, N1, N1)

  END SUBROUTINE ComputeRotB

  !-------------------------------------------

  SUBROUTINE FtoB
    USE VAR_d
    USE SOLVER, ONLY: FF1N1_d 
    IMPLICIT NONE

    CALL F1toB1_d <<<GRD_IJK, TBL_IJK>>>
    CALL F2toB2_d <<<GRD_IJK, TBL_IJK>>>
    CALL F3toB3_d <<<GRD_IJK, TBL_IJK>>>

    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (B1, N1, N1)
    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (B2, N1, N1)
    CALL FF1N1_d <<<GRD_IJ, TBL_IJ>>> (B3, N1, N1)

  END SUBROUTINE FtoB


  !---------------------------------------------

  SUBROUTINE ComputeCofPSI_C(N1)
    USE VAR_d, ONLY: GRD_IJK, TBL_IJK, GRD_JK, TBL_JK
    IMPLICIT NONE
    INTEGER :: N1


    CALL ComputeCofPSI_C_d <<<GRD_IJK, TBL_IJK>>> (N1)
    CALL ComputeApPSI_C_d <<<GRD_IJK, TBL_IJK>>> (N1)

  END SUBROUTINE ComputeCofPSI_C

 
  !---------------------------------------------
   
  SUBROUTINE ComputeCofPSI_OUT
    USE VAR_PSI_OUT_d
    IMPLICIT NONE
    
    CALL ComputeCofPSI_OUT_d <<<GRD_IJK, TBL_IJK>>>
    CALL ComputeBoundCofPSI_OUT_d <<<GRD_JK, TBL_JK>>>
    CALL ComputeBoundConPSI_OUT_d <<<GRD_JK, TBL_JK>>>
    CALL ComputeApPSI_OUT_d <<<GRD_IJK, TBL_IJK>>>

  END SUBROUTINE ComputeCofPSI_OUT


  !---------------------------------------------
  
  SUBROUTINE UpdateSourcePSI_OUT
    USE VAR_PSI_OUT_d
    IMPLICIT NONE
    
    CALL ComputeBoundConPSI_OUT_d <<<GRD_JK, TBL_JK>>>

  END SUBROUTINE UpdateSourcePSI_OUT

  !---------------------------------------------
  
  SUBROUTINE ComputeCofPSI_IN
    USE VAR_PSI_IN_d
    IMPLICIT NONE
    
    CALL ComputeCofPSI_IN_d <<<GRD_IJK, TBL_IJK>>>
    CALL ComputeBoundCofPSI_IN_d <<<GRD_JK, TBL_JK>>>
    CALL ComputeBoundConPSI_IN_d <<<GRD_JK, TBL_JK>>>
    CALL ComputeApPSI_IN_d <<<GRD_IJK, TBL_IJK>>>

  END SUBROUTINE ComputeCofPSI_IN
  
  !---------------------------------------------
  
  SUBROUTINE UpdateSourcePSI_IN
    USE VAR_PSI_IN_d
    IMPLICIT NONE
    
    CALL ComputeBoundConPSI_IN_d <<<GRD_JK, TBL_JK>>>

  END SUBROUTINE UpdateSourcePSI_IN
    
  !====================================================

  ATTRIBUTES (GLOBAL) SUBROUTINE CorrectF1_d(Iend)
  USE VAR_d, ONLY: F1, PSI_C, Xdif1, N2, N3
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, VALUE :: Iend
  
    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 3) .AND. (I <= (Iend)) .AND. &
         (J >= 2) .AND. (J <= (N2-1)) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      F1(I,J,K) = F1(I,J,K) - (PSI_C(I,J,K) - PSI_C(I-1,J,K)) / Xdif1(I)
   ENDIF

  END SUBROUTINE CorrectF1_d


  !---------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CorrectF2_d(Iend)
  USE VAR_d, ONLY: F2, PSI_C, Xdif2, N2, N3
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, VALUE :: Iend
  
    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 2) .AND. (I <= (Iend-1)) .AND. &
         (J >= 3) .AND. (J <= N2) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      F2(I,J,K) = F2(I,J,K) - (PSI_C(I,J,K) - PSI_C(I,J-1,K)) / Xdif2(J)

    ENDIF

  END SUBROUTINE CorrectF2_d

  !---------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE CorrectF3_d(Iend)
  USE VAR_d, ONLY: F3, PSI_C, Xdif3, N2, N3 
  IMPLICIT NONE
  INTEGER :: I, J, K
  INTEGER, VALUE :: Iend
  
    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 2) .AND. (I <= (Iend-1)) .AND. &
         (J >= 2) .AND. (J <= (N2-1)) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      F3(I,J,K) = F3(I,J,K) - (PSI_C(I,J,K) - PSI_C(I,J,K-1)) / Xdif3(K)
 
    ENDIF

  END SUBROUTINE CorrectF3_d


  !---------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeCofPSI_C_d(Iend)
    USE VAR_d, ONLY: AIM, AIP, AJM, AJP, AKM, AKP, CON, AP, &
                     S1, S2, S3, B1, B2, B3, X1, &
                     Xdif1, Xdif2, Xdif3, sinX2, &
                     N2, N3   
    IMPLICIT NONE
    INTEGER :: I, J, K
    INTEGER, VALUE :: Iend

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 2) .AND. (I <= (Iend-1)) .AND. &
         (J >= 2) .AND. (J <= (N2-1)) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      AIM(I,J,K) = S1(I,J,K) / Xdif1(I)                                
      AIP(I,J,K) = S1(I+1,J,K) / Xdif1(I+1)
        
      AJM(I,J,K) = S2(I,J,K) / (X1(I)*Xdif2(J))              
      AJP(I,J,K) = S2(I,J+1,K) / (X1(I)*Xdif2(J+1))

      AKM(I,J,K) =  S3(I,J,K)/(X1(I)*sinX2(J)*Xdif3(K))
      AKP(I,J,K) =  S3(I,J,K+1)/(X1(I)*sinX2(J)*Xdif3(K+1))

      CON(I,J,K) = - ( S1(I+1,J,K)*B1(I+1,J,K) - S1(I,J,K)*B1(I,J,K) + &
                       S2(I,J+1,K)*B2(I,J+1,K) - S2(I,J,K)*B2(I,J,K) + &
                       S3(I,J,K+1)*B3(I,J,K+1) - S3(I,J,K)*B3(I,J,K) )
      AP(I,J,K)  = 0.0_8
    ENDIF

  END SUBROUTINE ComputeCofPSI_C_d


  !---------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeBoundRhoCofPSI_C_d(N1)
    USE VAR_d, ONLY: AIP, N2, N3
    IMPLICIT NONE
    INTEGER :: I, J, K
    INTEGER, VALUE :: N1

    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= 2) .AND. (J <= (N2-1)) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      AIP(N1-1,J,K) = 0.0_8
    ENDIF

  END SUBROUTINE ComputeBoundRhoCofPSI_C_d
  

  !---------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeApPSI_C_d(Iend)
    USE VAR_d, ONLY: AP, AIM, AIP, AJM, AJP, AKM, AKP, N2, N3
    IMPLICIT NONE
    INTEGER :: I, J, K
    INTEGER, VALUE :: Iend

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 2) .AND. (I <= (Iend-1)) .AND. &
         (J >= 2) .AND. (J <= (N2-1)) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

     AP(I,J,K) = - AP(I,J,K) + AIM(I,J,K)+AIP(I,J,K)+ &
                               AJM(I,J,K)+AJP(I,J,K)+ &
                               AKM(I,J,K)+AKP(I,J,K)
    ENDIF

  END SUBROUTINE ComputeApPSI_C_d


  !-------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeJ1_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 2) .AND. (I <= N1) .AND. &
         (J >= 3) .AND. (J <= N2) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      J1(I,J,K) = ( ( B3(I,J,K) * sinX2(J) - B3(I,J-1,K) * sinX2(J-1) ) * X1(I) * Xdif3(K) - &
                    ( B2(I,J,K) - B2(I,J,K-1) ) * X1(I) * Xdif2(J) ) / &
                    ( X1(I)**2 * Xdif3(K) * ( cosX2(J-1) - cosX2(J)) )
 
    ENDIF

  END SUBROUTINE ComputeJ1_d
  
  !-------------------------------------------
  
  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeJ1_Oz_d
    USE AFUN_d
    USE VAR_d
    IMPLICIT NONE
    INTEGER :: I, K
    REAL(8) :: rotB_r

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
 
    IF ( (I >= 2 .AND. I <= N1) ) THEN
    
      rotB_r = 0.0_8

      DO K = 2, N3-1
        rotB_r = rotB_r + ( ( F3(I,3,K) - F3(I,2,K) ) * Xdif3(K) - ( F2(I,3,K) - F2(I,3,K-1) ) * Xdif2(3) ) / &
                            ( X1(I)**2 * Xdif3(K) * ( cosX2(2) - cosX2(3)) )
      ENDDO

     rotB_r = rotB_r / DBLE(N3-2)
      
      DO K = 1, N3
        J1(I,2,K) = rotB_r
      ENDDO      
        
    ENDIF

  END SUBROUTINE ComputeJ1_Oz_d

  !-------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeJ2_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 2) .AND. (I <= N1) .AND. &
         (J >= 2) .AND. (J <= (N2-1)) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      J2(I,J,K) = ((B1(I,J,K) - B1(I,J,K-1)) * Xdif1(I) - &
                   (B3(I,J,K)*X1(I) - B3(I-1,J,K)*X1(I-1)) * sinX2(J)*Xdif3(K)) / &
                   (0.5_8 * (X1(I)*X1(I)-X1(I-1)*X1(I-1)) * sinX2(J)*Xdif3(K))
 
    ENDIF

  END SUBROUTINE ComputeJ2_d

  !-------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeJ3_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 2) .AND. (I <= N1) .AND. &
         (J >= 3) .AND. (J <= N2) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      J3(I,J,K) = ((B1(I,J-1,K) - B1(I,J,K))*Xdif1(I) + &
                   (B2(I,J,K)*X1(I) - B2(I-1,J,K)*X1(I-1)) * Xdif2(J)) / &
                   (0.5_8 * (X1(I)*X1(I)-X1(I-1)*X1(I-1))*Xdif2(J))
 
    ENDIF

  END SUBROUTINE ComputeJ3_d

  
  !-------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE F1toB1_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: H1

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 3) .AND. (I <= N1) .AND. &
         (J >= 2) .AND. (J <= N2) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      H1 = 1.0_8
      B1(I,J,K) = F1(I,J,K) / H1 
  
    ENDIF

  END SUBROUTINE F1toB1_d

  !---------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE F2toB2_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: H2

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z


    IF ( (I >= 2) .AND. (I <= N1) .AND. &
         (J >= 3) .AND. (J <= N2) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      H2 = X1(I)
      B2(I,J,K) = F2(I,J,K) / H2

    ENDIF

  END SUBROUTINE F2toB2_d

  !---------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE F3toB3_d
  USE VAR_d
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: H3

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 2) .AND. (I <= N1) .AND. &
         (J >= 2) .AND. (J <= N2) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      H3 = X1(I) * sinX2(J)
      B3(I,J,K) = F3(I,J,K) / H3

    ENDIF

  END SUBROUTINE F3toB3_d

  !---------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeF1_Oz_d
    USE AFUN_d
    USE VAR_d
    IMPLICIT NONE
    INTEGER I, K
    REAL(8) SUM_F1

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
 
    IF ( (I >= 3 .AND. I <= N1) ) THEN
    
      SUM_F1 = 0.0_8
      
      DO K = 2, N3-1
        SUM_F1 = SUM_F1 + F1(I,2,K)
      ENDDO
      
      F1(I, 1, K) = SUM_F1 / DBLE(N3-2)
    ENDIF

  END SUBROUTINE ComputeF1_Oz_d


  !---------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeCofF1_d
    USE AFUN_d
    USE VAR_d
    IMPLICIT NONE
    INTEGER :: I, J, K
    REAL(8) :: AP0, V1, V2, V3, SUM_FLOW, SUM_CON, AIM_F2, AIM_F3, DIFF, FLOW, SUM_CON_Q
    REAL(8) :: Umax, B0, SQRT_E, Ha, Ux, V1a, Bx, By, F2_a, F3_a, zd

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 3 .AND. I <= N1) .AND.     &
         (J >= 2 .AND. J <= (N2-1)) .AND. &
         (K >= 2 .AND. K <= (N3-1)) ) THEN

      AIM(I,J,K) = 0.0_8; AIP(I,J,K) = 0.0_8

      AP0 = S1(I,J,K) / DT

      SUM_FLOW = 0.0_8
      SUM_CON_Q = 0.0_8 
      
      DIFF = J321(I,J,K)  
      V2 = WEIGHT_AVG(U2(I-1,J,K), U2(I,J,K), Kx1(I))
      FLOW = F_321(I, J, K) * V2
      SUM_FLOW = SUM_FLOW - FLOW
      SUM_CON_Q = SUM_CON_Q + Q(FLOW, Xs2(J), X2(J-2), X2(J-1), X2(J), F1(I,J-2,K), F1(I,J-1,K), F1(I,J,K)) &
                            - Q(-FLOW, Xs2(J), X2(J-1), X2(J), X2(J+1), F1(I,J-1,K), F1(I,J,K), F1(I,J+1,K))        
      AJM(I,J,K) =  DIFF + DMAX1(0.0_8, FLOW)


      DIFF = J321(I, J+1, K)
      V2 = WEIGHT_AVG(U2(I-1,J+1,K), U2(I,J+1,K), Kx1(I))
      FLOW = F_321(I, J+1, K) * V2
      SUM_FLOW = SUM_FLOW + FLOW
      SUM_CON_Q = SUM_CON_Q - Q(FLOW, Xs2(J+1), X2(J-1), X2(J), X2(J+1), F1(I,J-1,K), F1(I,J,K), F1(I,J+1,K)) &
                            + Q(-FLOW, Xs2(J+1), X2(J), X2(J+1), X2(J+2), F1(I,J,K), F1(I,J+1,K), F1(I,J+2,K))   
      AJP(I,J,K) =  DIFF + DMAX1(0.0_8, -FLOW)


      DIFF = J231(I,J,K)  
      V3 = WEIGHT_AVG(U3(I-1,J,K), U3(I,J,K), Kx1(I))
      FLOW = F_231(I, J, K) * V3
      SUM_FLOW = SUM_FLOW - FLOW
      SUM_CON_Q = SUM_CON_Q &
                + Q(FLOW, Xs3(K), X_3(K-2), X_3(K-1), X_3(K), F_1(I,J,K-2), F_1(I,J,K-1), F_1(I,J,K)) &
                - Q(-FLOW, Xs3(K), X_3(K-1), X_3(K), X_3(K+1), F_1(I,J,K-1), F_1(I,J,K), F_1(I,J,K+1))   
      AKM(I,J,K) = DIFF + DMAX1(0.0_8, FLOW)

      DIFF = J231(I,J,K+1)  
      V3 = WEIGHT_AVG(U3(I-1,J,K+1), U3(I,J,K+1), Kx1(I))
      FLOW = F_231(I, J, K+1) * V3
      SUM_FLOW = SUM_FLOW + FLOW
      SUM_CON_Q = SUM_CON_Q &
                - Q(FLOW, Xs3(K+1), X_3(K-1), X_3(K), X_3(K+1), F_1(I,J,K-1), F_1(I,J,K), F_1(I,J,K+1)) &
                + Q(-FLOW, Xs3(K+1), X_3(K), X_3(K+1), X_3(K+2), F_1(I,J,K), F_1(I,J,K+1), F_1(I,J,K+2))
      AKP(I,J,K) = DIFF + DMAX1(0.0_8, -FLOW)
      
      AP(I,J,K) = AP0 + AJM(I,J,K) + AJP(I,J,K) + AKM(I,J,K) + AKP(I,J,K) ! + SUM_AP
  
      !-------------------- Sc-------------------------
   
      SUM_CON = 0.0_8 

      DIFF = J312(I,J,K)  
      V1 = WEIGHT_AVG(U1(I,J-1,K), U1(I,J,K), Kx2(J))
      FLOW = F_312(I, J, K) * V1          
      SUM_CON = SUM_CON  - DIFF * (F2(I-1,J,K) - F2(I,J,K)) &
                - Qsc(FLOW, Xs1(I), X1(I-2), X1(I-1), X1(I), F2(I-2,J,K), F2(I-1,J,K), F2(I,J,K)) &
                + Qsc(-FLOW, Xs1(I), X1(I-1), X1(I), X1(I+1), F2(I-1,J,K), F2(I,J,K), F2(I+1,J,K))

      DIFF = J312(I,J+1,K)  
      V1 = WEIGHT_AVG(U1(I,J,K), U1(I,J+1,K), Kx2(J+1))
      FLOW = F_312(I,J+1,K) * V1       
      SUM_CON = SUM_CON  + DIFF * (F2(I-1,J+1,K) - F2(I,J+1,K)) &
                + Qsc(FLOW, Xs1(I), X1(I-2), X1(I-1), X1(I), F2(I-2,J+1,K), F2(I-1,J+1,K), F2(I,J+1,K)) &
                - Qsc(-FLOW, Xs1(I), X1(I-1), X1(I), X1(I+1), F2(I-1,J+1,K), F2(I,J+1,K), F2(I+1,J+1,K))

  
      DIFF = J213(I,J,K)  
      V1 = WEIGHT_AVG(U1(I,J,K-1), U1(I,J,K), Kx3(K))
      FLOW = F_213(I, J, K) * V1        
      SUM_CON = SUM_CON  - DIFF * (F3(I-1,J,K) - F3(I,J,K)) &
                - Qsc(FLOW, Xs1(I), X1(I-2), X1(I-1), X1(I), F3(I-2,J,K), F3(I-1,J,K), F3(I,J,K)) &
                + Qsc(-FLOW, Xs1(I), X1(I-1), X1(I), X1(I+1), F3(I-1,J,K), F3(I,J,K), F3(I+1,J,K))  

      DIFF = J213(I,J,K+1)  
      V1 = WEIGHT_AVG(U1(I,J,K), U1(I,J,K+1), Kx3(K+1))
      FLOW = F_213(I,J,K+1) * V1        
      SUM_CON = SUM_CON  + DIFF * (F3(I-1,J,K+1) - F3(I,J,K+1)) &
                + Qsc(FLOW, Xs1(I), X1(I-2), X1(I-1), X1(I), F3(I-2,J,K+1), F3(I-1,J,K+1), F3(I,J,K+1)) &
                - Qsc(-FLOW, Xs1(I), X1(I-1), X1(I), X1(I+1), F3(I-1,J,K+1), F3(I,J,K+1), F3(I+1,J,K+1))  


      CON(I,J,K) = SUM_CON + AP0 * F1_0(I,J,K) - SUM_FLOW * F1(I,J,K)

      !--------------------------------------QUICK----------------------------------

      CON(I,J,K) = CON(I,J,K) + SUM_CON_Q

    ENDIF

  END SUBROUTINE ComputeCofF1_d


  !---------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeCofF2_d
    USE AFUN_d
    USE VAR_d
    IMPLICIT NONE
    INTEGER I, J, K
    REAL(8) AP0, V1, V2, V3, SUM_FLOW, SUM_CON, AJM_F1, AJM_F3, DIFF, FLOW, SUM_CON_Q 

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 2 .AND. I <= (N1-1)) .AND. &
         (J >= 3 .AND. J <= N2) .AND. &
         (K >= 2 .AND. K <= (N3-1)) ) THEN

      AJM(I,J,K) = 0.0_8; AJP(I,J,K) = 0.0_8

      AP0 = S2(I,J,K) / (X1(I) * DT)

      SUM_FLOW = 0.0_8
      SUM_CON_Q = 0.0_8
       
      DIFF = J312(I,J,K)  
      V1 = WEIGHT_AVG(U1(I,J-1,K), U1(I,J,K), Kx2(J))
      FLOW = F_312(I, J, K) * V1
      SUM_FLOW = SUM_FLOW - FLOW 
      SUM_CON_Q = SUM_CON_Q + Q(FLOW, Xs1(I), X1(I-2), X1(I-1), X1(I), F2(I-2,J,K), F2(I-1,J,K), F2(I,J,K)) &
                            - Q(-FLOW, Xs1(I), X1(I-1), X1(I), X1(I+1), F2(I-1,J,K), F2(I,J,K), F2(I+1,J,K))              
      AIM(I,J,K) = DIFF + DMAX1(0.0_8, FLOW)


      DIFF = J312(I+1, J, K)
      V1 = WEIGHT_AVG(U1(I+1,J-1,K), U1(I+1,J,K), Kx2(J))
      FLOW = F_312(I+1, J, K) * V1
      SUM_FLOW = SUM_FLOW + FLOW
      SUM_CON_Q = SUM_CON_Q - Q(FLOW, Xs1(I+1), X1(I-1), X1(I), X1(I+1), F2(I-1,J,K), F2(I,J,K), F2(I+1,J,K)) &
                            + Q(-FLOW, Xs1(I+1), X1(I), X1(I+1), X1(I+2), F2(I,J,K), F2(I+1,J,K), F2(I+2,J,K))    
      AIP(I,J,K) = DIFF + DMAX1(0.0_8, -FLOW)


      DIFF = J132(I,J,K)  
      V3 = WEIGHT_AVG(U3(I,J-1,K), U3(I,J,K), Kx2(J))
      FLOW = F_132(I, J, K) * V3
      SUM_FLOW = SUM_FLOW - FLOW
      SUM_CON_Q = SUM_CON_Q &
                + Q(FLOW, Xs3(K), X_3(K-2), X_3(K-1), X_3(K), F_2(I,J,K-2), F_2(I,J,K-1), F_2(I,J,K)) &
                - Q(-FLOW, Xs3(K), X_3(K-1), X_3(K), X_3(K+1), F_2(I,J,K-1), F_2(I,J,K), F_2(I,J,K+1))   
      AKM(I,J,K) = DIFF + DMAX1(0.0_8, FLOW)


      DIFF = J132(I,J,K+1)  
      V3 = WEIGHT_AVG(U3(I,J-1,K+1), U3(I,J,K+1), Kx2(J))
      FLOW = F_132(I, J, K+1) * V3
      SUM_FLOW = SUM_FLOW + FLOW  
      SUM_CON_Q = SUM_CON_Q &
                - Q(FLOW, Xs3(K+1), X_3(K-1), X_3(K), X_3(K+1), F_2(I,J,K-1), F_2(I,J,K), F_2(I,J,K+1)) &
                + Q(-FLOW, Xs3(K+1), X_3(K), X_3(K+1), X_3(K+2), F_2(I,J,K), F_2(I,J,K+1), F_2(I,J,K+2))  
      AKP(I,J,K) = DIFF + DMAX1(0.0_8, -FLOW)

      
      AP(I,J,K) = AP0 + AIM(I,J,K) + AIP(I,J,K) + AKM(I,J,K) + AKP(I,J,K) !+ SUM_FLOW
  
  
    !-------------------- Sc-------------------------

      SUM_CON = 0.0_8 

      DIFF = J321(I,J,K)  
      V2 = WEIGHT_AVG(U2(I-1,J,K), U2(I,J,K), Kx1(I))
      FLOW = F_321(I, J, K) * V2          
      SUM_CON = SUM_CON - DIFF * (F1(I,J-1,K) - F1(I,J,K))&
                - Qsc(FLOW, Xs2(J), X2(J-2), X2(J-1), X2(J), F1(I,J-2,K), F1(I,J-1,K), F1(I,J,K)) &
                + Qsc(-FLOW, Xs2(J), X2(J-1), X2(J), X2(J+1), F1(I,J-1,K), F1(I,J,K), F1(I,J+1,K)) 


      DIFF = J321(I+1,J,K)  
      V2 = WEIGHT_AVG(U2(I,J,K), U2(I+1,J,K), Kx1(I+1))
      FLOW = F_321(I+1,J,K) * V2       
      SUM_CON = SUM_CON + DIFF * (F1(I+1,J-1,K) - F1(I+1,J,K))&
                + Qsc(FLOW, Xs2(J), X2(J-2), X2(J-1), X2(J), F1(I+1,J-2,K), F1(I+1,J-1,K), F1(I+1,J,K)) &
                - Qsc(-FLOW, Xs2(J), X2(J-1), X2(J), X2(J+1), F1(I+1,J-1,K), F1(I+1,J,K), F1(I+1,J+1,K)) 



      DIFF = J123(I,J,K)  
      V2 = WEIGHT_AVG(U2(I,J,K-1), U2(I,J,K), Kx3(K))
      FLOW = F_123(I, J, K) * V2        
      SUM_CON = SUM_CON  - DIFF * (F3(I,J-1,K) - F3(I,J,K)) &
                - Qsc(FLOW, Xs2(J), X2(J-2), X2(J-1), X2(J), F3(I,J-2,K), F3(I,J-1,K), F3(I,J,K)) &
                + Qsc(-FLOW, Xs2(J), X2(J-1), X2(J), X2(J+1), F3(I,J-1,K), F3(I,J,K), F3(I,J+1,K)) 


      DIFF = J123(I,J,K+1)  
      V2 = WEIGHT_AVG(U2(I,J,K), U2(I,J,K+1), Kx3(K+1))
      FLOW = F_123(I,J,K+1) * V2        
      SUM_CON = SUM_CON  + DIFF * (F3(I,J-1,K+1) - F3(I,J,K+1)) &
                + Qsc(FLOW, Xs2(J), X2(J-2), X2(J-1), X2(J), F3(I,J-2,K+1), F3(I,J-1,K+1), F3(I,J,K+1)) &
                - Qsc(-FLOW, Xs2(J), X2(J-1), X2(J), X2(J+1), F3(I,J-1,K+1), F3(I,J,K+1), F3(I,J+1,K+1))   

      CON(I,J,K) = SUM_CON + AP0 * F2_0(I,J,K) - SUM_FLOW * F2(I,J,K)

      !--------------------------------------QUICK----------------------------------

      CON(I,J,K) = CON(I,J,K) + SUM_CON_Q

    ENDIF

  END SUBROUTINE ComputeCofF2_d


  !---------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeCofF3_d
    USE AFUN_d
    USE VAR_d
    IMPLICIT NONE
    INTEGER I, J, K
    REAL(8) AP0, V1, V2, V3, SUM_FLOW, SUM_CON, AKM_F1, AKM_F2, DIFF, FLOW, SUM_CON_Q

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 2 .AND. I <= (N1-1)) .AND. &
         (J >= 3 .AND. J <= (N2-1)) .AND. &
         (K >= 2 .AND. K <= (N3-1)) ) THEN

      AKM(I,J,K) = 0.0_8; AKP(I,J,K) = 0.0_8

      AP0 = S3(I,J,K) / (X1(I) * sinX2(J) * DT)

      SUM_FLOW = 0.0_8
      SUM_CON_Q = 0.0_8
      
      DIFF = J213(I,J,K)  
      V1 = WEIGHT_AVG(U1(I,J,K-1), U1(I,J,K), Kx3(K))
      FLOW = F_213(I, J, K) * V1
      SUM_FLOW = SUM_FLOW - FLOW 
      SUM_CON_Q = SUM_CON_Q + Q(FLOW, Xs1(I), X1(I-2), X1(I-1), X1(I), F3(I-2,J,K), F3(I-1,J,K), F3(I,J,K)) &
                            - Q(-FLOW, Xs1(I), X1(I-1), X1(I), X1(I+1), F3(I-1,J,K), F3(I,J,K), F3(I+1,J,K))         
      AIM(I,J,K) = DIFF + DMAX1(0.0_8, FLOW)


      DIFF = J213(I+1, J, K)
      V1 = WEIGHT_AVG(U1(I+1,J,K-1), U1(I+1,J,K), Kx3(K))
      FLOW = F_213(I+1, J, K) * V1
      SUM_FLOW = SUM_FLOW + FLOW 
      SUM_CON_Q = SUM_CON_Q - Q(FLOW, Xs1(I+1), X1(I-1), X1(I), X1(I+1), F3(I-1,J,K), F3(I,J,K), F3(I+1,J,K)) &
                            + Q(-FLOW, Xs1(I+1), X1(I), X1(I+1), X1(I+2), F3(I,J,K), F3(I+1,J,K), F3(I+2,J,K))   
      AIP(I,J,K) = DIFF + DMAX1(0.0_8, -FLOW)


      DIFF = J123(I,J,K)  
      V2 = WEIGHT_AVG(U2(I,J,K-1), U2(I,J,K), Kx3(K)) 
      FLOW = F_123(I, J, K) * V2
      SUM_FLOW = SUM_FLOW - FLOW
      SUM_CON_Q = SUM_CON_Q + Q(FLOW, Xs2(J), X2(J-2), X2(J-1), X2(J), F3(I,J-2,K), F3(I,J-1,K), F3(I,J,K)) &
                            - Q(-FLOW, Xs2(J), X2(J-1), X2(J), X2(J+1), F3(I,J-1,K), F3(I,J,K), F3(I,J+1,K))          
      AJM(I,J,K) =  DIFF + DMAX1(0.0_8, FLOW)


      DIFF = J123(I, J+1, K)
      V2 = WEIGHT_AVG(U2(I,J+1,K-1), U2(I,J+1,K), Kx3(K))
      FLOW = F_123(I, J+1, K) * V2
      SUM_FLOW = SUM_FLOW + FLOW 
      SUM_CON_Q = SUM_CON_Q - Q(FLOW, Xs2(J+1), X2(J-1), X2(J), X2(J+1), F3(I,J-1,K), F3(I,J,K), F3(I,J+1,K)) &
                            + Q(-FLOW, Xs2(J+1), X2(J), X2(J+1), X2(J+2), F3(I,J,K), F3(I,J+1,K), F3(I,J+2,K))   
      AJP(I,J,K) =  DIFF + DMAX1(0.0_8, -FLOW)

      
      AP(I,J,K) = AP0 + AIM(I,J,K) + AIP(I,J,K) + AJM(I,J,K) + AJP(I,J,K)
  
  
      !-------------------- Sc-------------------------

      SUM_CON = 0.0_8 

      DIFF = J231(I,J,K)  
      V3 = WEIGHT_AVG(U3(I-1,J,K), U3(I,J,K), Kx1(I))
      FLOW = F_231(I, J, K) * V3        
      SUM_CON = SUM_CON - DIFF * (F1(I,J,K-1) - F1(I,J,K))&
                - Qsc(FLOW, Xs3(K), X_3(K-2), X_3(K-1), X_3(K), F_1(I,J,K-2), F_1(I,J,K-1), F_1(I,J,K)) &
                + Qsc(-FLOW, Xs3(K), X_3(K-1), X_3(K), X_3(K+1), F_1(I,J,K-1), F_1(I,J,K), F_1(I,J,K+1)) 


      DIFF = J231(I+1,J,K)  
      V3 = WEIGHT_AVG(U3(I,J,K), U3(I+1,J,K), Kx1(I+1))
      FLOW = F_231(I+1,J,K) * V3        
      SUM_CON = SUM_CON + DIFF * (F1(I+1,J,K-1) - F1(I+1,J,K))&
                + Qsc(FLOW, Xs3(K), X_3(K-2), X_3(K-1), X_3(K), F_1(I+1,J,K-2), F_1(I+1,J,K-1), F_1(I+1,J,K)) &
                - Qsc(-FLOW, Xs3(K), X_3(K-1), X_3(K), X_3(K+1), F_1(I+1,J,K-1), F_1(I+1,J,K), F_1(I+1,J,K+1)) 



      DIFF = J132(I,J,K)  
      V3 = WEIGHT_AVG(U3(I,J-1,K), U3(I,J,K), Kx2(J))
      FLOW = F_132(I, J, K) * V3          
      SUM_CON = SUM_CON  - DIFF * (F2(I,J,K-1) - F2(I,J,K))&
                - Qsc(FLOW, Xs3(K), X_3(K-2), X_3(K-1), X_3(K), F_2(I,J,K-2), F_2(I,J,K-1), F_2(I,J,K)) &
                + Qsc(-FLOW, Xs3(K), X_3(K-1), X_3(K), X_3(K+1), F_2(I,J,K-1), F_2(I,J,K), F_2(I,J,K+1)) 


      DIFF = J132(I,J+1,K)  
      V3 = WEIGHT_AVG(U3(I,J,K), U3(I,J+1,K), Kx2(J+1))
      FLOW = F_132(I,J+1,K) * V3       
      SUM_CON = SUM_CON  + DIFF * (F2(I,J+1,K-1) - F2(I,J+1,K)) &
                + Qsc(FLOW, Xs3(K), X_3(K-2), X_3(K-1), X_3(K), F_2(I,J+1,K-2), F_2(I,J+1,K-1), F_2(I,J+1,K)) &
                - Qsc(-FLOW, Xs3(K), X_3(K-1), X_3(K), X_3(K+1), F_2(I,J+1,K-1), F_2(I,J+1,K), F_2(I,J+1,K+1)) 


      CON(I,J,K) = SUM_CON + AP0 * F3_0(I,J,K) - SUM_FLOW * F3(I,J,K)

      !--------------------------------------QUICK----------------------------------

      CON(I,J,K) = CON(I,J,K) + SUM_CON_Q


    ENDIF

  END SUBROUTINE ComputeCofF3_d

  
  !---------------------------------------------
 
  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeE1_Oz_d
    USE AFUN_d
    USE VAR_d
    IMPLICIT NONE
    INTEGER :: I, K
    REAL(8) :: UxB_r, rotB_r, H2, H3

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x

 
    IF ( (I >= 2 .AND. I <= (N1-1)) ) THEN
    
      UxB_r = 0.0_8
      rotB_r = 0.0_8
      H2 = X1(I)
      H3 = X1(I) * sinX2(2)      

      DO K = 2, N3-1
        UxB_r = UxB_r + (U2(I,3,K) * 0.5_8 * (F3(I,2,K) + F3(I,2,K+1)) / H3 - 0.5_8 * (U3(I,2,K) + U3(I,2,K+1)) * F2(I,3,K) / H2) 
        !rotB_r = rotB_r + 0.5_8 * (F3(I,2,K) + F3(I,2,K+1)) * (Xs3(K+1) - Xs3(K)) / &
        !         (2.0_8 * PI * X1(I)**2 * (1.0_8 - cosX2(2)))  
        rotB_r = rotB_r + ( ( F3(I,3,K) - F3(I,2,K) ) * Xdif3(K) - ( F2(I,3,K) - F2(I,3,K-1) ) * Xdif2(3) ) / &
                            ( X1(I)**2 * Xdif3(K) * ( cosX2(2) - cosX2(3)) )
      ENDDO
 
      UxB_r = UxB_r / DBLE(N3-2)      
      rotB_r = rotB_r / DBLE(N3-2)  
      E1_Oz(I) = (UxB_r - viscm * rotB_r) * (Xs1(I+1) - Xs1(I)) 
      
    ENDIF


  END SUBROUTINE ComputeE1_Oz_d


  !---------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeCofF3_Oz_d
    USE AFUN_d
    USE VAR_d
    IMPLICIT NONE
    INTEGER I, K
    REAL(8) AP0, V1, V2, V3, SUM_FLOW, SUM_CON, AKM_F1, AKM_F2, DIFF, FLOW, SUM_CON_Q

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
 
    IF ( (I >= 2 .AND. I <= (N1-1)) .AND. &
         (K >= 2 .AND. K <= (N3-1)) ) THEN

      AKM(I,2,K) = 0.0_8; AKP(I,2,K) = 0.0_8; AJM(I,2,K) = 0.0_8

      AP0 = S3(I,2,K) / (X1(I) * sinX2(2) * DT)

      SUM_FLOW = 0.0_8
      SUM_CON_Q = 0.0_8 
         
      DIFF = J213(I,2,K)  
      V1 = WEIGHT_AVG(U1(I,2,K-1), U1(I,2,K), Kx3(K))
      FLOW = F_213(I, 2, K) * V1
      SUM_FLOW = SUM_FLOW - FLOW
      SUM_CON_Q = SUM_CON_Q + Q(FLOW, Xs1(I), X1(I-2), X1(I-1), X1(I), F3(I-2,2,K), F3(I-1,2,K), F3(I,2,K)) &
                            - Q(-FLOW, Xs1(I), X1(I-1), X1(I), X1(I+1), F3(I-1,2,K), F3(I,2,K), F3(I+1,2,K))              
      AIM(I,2,K) = DIFF + DMAX1(0.0_8, FLOW)



      DIFF = J213(I+1, 2, K)
      V1 = WEIGHT_AVG(U1(I+1,2,K-1), U1(I+1,2,K), Kx3(K))
      FLOW = F_213(I+1, 2, K) * V1
      SUM_FLOW = SUM_FLOW + FLOW
      SUM_CON_Q = SUM_CON_Q - Q(FLOW, Xs1(I+1), X1(I-1), X1(I), X1(I+1), F3(I-1,2,K), F3(I,2,K), F3(I+1,2,K)) &
                            + Q(-FLOW, Xs1(I+1), X1(I), X1(I+1), X1(I+2), F3(I,2,K), F3(I+1,2,K), F3(I+2,2,K))     
      AIP(I,2,K) = DIFF + DMAX1(0.0_8, -FLOW)

      DIFF = J123(I, 3, K)
      V2 = WEIGHT_AVG(U2(I,3,K-1), U2(I,3,K), Kx3(K))
      FLOW = F_123(I, 3, K) * V2
      SUM_FLOW = SUM_FLOW + FLOW
      SUM_CON_Q = SUM_CON_Q - Q(FLOW, Xs2(3), X2(1), X2(2), X2(3), F3(I,1,K), F3(I,2,K), F3(I,3,K)) &
                            + Q(-FLOW, Xs2(3), X2(2), X2(3), X2(4), F3(I,2,K), F3(I,3,K), F3(I,4,K))    
      AJP(I,2,K) =  DIFF + DMAX1(0.0_8, -FLOW)

      AP(I,2,K) = AP0 + AIM(I,2,K) + AIP(I,2,K) + AJP(I,2,K)
  
  
      !-------------------- Sc-------------------------

      SUM_CON = 0.0_8 

      DIFF = J231(I,2,K)  
      V3 = WEIGHT_AVG(U3(I-1,2,K), U3(I,2,K), Kx1(I))
      FLOW = F_231(I, 2, K) * V3        
      SUM_CON = SUM_CON - DIFF * (F1(I,2,K-1) - F1(I,2,K))&
                - Qsc(FLOW, Xs3(K), X_3(K-2), X_3(K-1), X_3(K), F_1(I,2,K-2), F_1(I,2,K-1), F_1(I,2,K)) &
                + Qsc(-FLOW, Xs3(K), X_3(K-1), X_3(K), X_3(K+1), F_1(I,2,K-1), F_1(I,2,K), F_1(I,2,K+1)) 

 

      DIFF = J231(I+1,2,K)  
      V3 = WEIGHT_AVG(U3(I,2,K), U3(I+1,2,K), Kx1(I+1))
      FLOW = F_231(I+1,2,K) * V3        
      SUM_CON = SUM_CON + DIFF * (F1(I+1,2,K-1) - F1(I+1,2,K))&
                + Qsc(FLOW, Xs3(K), X_3(K-2), X_3(K-1), X_3(K), F_1(I+1,2,K-2), F_1(I+1,2,K-1), F_1(I+1,2,K)) &
                - Qsc(-FLOW, Xs3(K), X_3(K-1), X_3(K), X_3(K+1), F_1(I+1,2,K-1), F_1(I+1,2,K), F_1(I+1,2,K+1)) 

      DIFF = J132(I,3,K)  
      V3 = WEIGHT_AVG(U3(I,2,K), U3(I,3,K), Kx2(3))
      FLOW = F_132(I,3,K) * V3       
      SUM_CON = SUM_CON  + DIFF * (F2(I,3,K-1) - F2(I,3,K))&
                + Qsc(FLOW, Xs3(K), X_3(K-2), X_3(K-1), X_3(K), F_2(I,3,K-2), F_2(I,3,K-1), F_2(I,3,K)) &
                - Qsc(-FLOW, Xs3(K), X_3(K-1), X_3(K), X_3(K+1), F_2(I,3,K-1), F_2(I,3,K), F_2(I,3,K+1)) 

      CON(I,2,K) = SUM_CON + AP0 * F3_0(I,2,K) + E1_Oz(I) - SUM_FLOW * F3(I,2,K) 

      !--------------------------------------QUICK----------------------------------

      CON(I,2,K) = CON(I,2,K) + SUM_CON_Q

    ENDIF

  END SUBROUTINE ComputeCofF3_Oz_d


  !---------------------------------------------
  
  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeCofPSI_OUT_d
    USE VAR_d, ONLY: sinX2, Xdif2, Xdif3
    USE VAR_PSI_OUT_d
    IMPLICIT NONE
    INTEGER I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 2) .AND. (I <= (N1-1)) .AND. &
         (J >= 2) .AND. (J <= (N2-1)) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      AIM(I,J,K) = S1(I,J,K) / Xdif1(I)                                
      AIP(I,J,K) = S1(I+1,J,K) / Xdif1(I+1)
        
      AJM(I,J,K) = S2(I,J,K) / (X1(I) * Xdif2(J))              
      AJP(I,J,K) = S2(I,J+1,K) / (X1(I) * Xdif2(J+1))

      AKM(I,J,K) = S3(I,J,K) / (X1(I) * sinX2(J) * Xdif3(K))
      AKP(I,J,K) = S3(I,J,K+1) / (X1(I) * sinX2(J) * Xdif3(K+1))

      CON(I,J,K) = 0.0_8
      AP(I,J,K)  = 0.0_8

         
    ENDIF

  END SUBROUTINE ComputeCofPSI_OUT_d


  !---------------------------------------------


  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeBoundCofPSI_OUT_d
    USE VAR_d, ONLY: F1, r_o, N1_mhd => N1
    USE VAR_PSI_OUT_d
    IMPLICIT NONE
    INTEGER I, J, K

    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= 2) .AND. (J <= (N2-1)) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN
  
      AIP(N1-1,J,K) = 0.0_8
      !CON(N1-1,J,K) = CON(N1-1,J,K) - F1(N1_mhd, J, K) * S1(N1,J,K) / (r_o + Xdif1(N1))
      AP(N1-1,J,K)  = AP(N1-1,J,K) - S1(N1,J,K) / (r_o + Xdif1(N1))  
                       
    ENDIF

  END SUBROUTINE ComputeBoundCofPSI_OUT_d

  !---------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeApPSI_OUT_d
    USE VAR_PSI_OUT_d
    IMPLICIT NONE
    INTEGER I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 2) .AND. (I <= (N1-1)) .AND. &
         (J >= 2) .AND. (J <= (N2-1)) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN
         
      AP(I,J,K) = - AP(I,J,K) + AIM(I,J,K) + AIP(I,J,K) + &
                                AJM(I,J,K) + AJP(I,J,K) + &
                                AKM(I,J,K) + AKP(I,J,K)
         
    ENDIF

  END SUBROUTINE ComputeApPSI_OUT_d

  
  !---------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeBoundConPSI_OUT_d
    USE VAR_d, ONLY: F1, r_o, N1_mhd => N1
    USE VAR_PSI_OUT_d
    IMPLICIT NONE
    INTEGER I, J, K

    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= 2) .AND. (J <= (N2-1)) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN
         
      CON(N1-1,J,K) = - F1(N1_mhd, J, K) * S1(N1,J,K) / (r_o + Xdif1(N1))
                     
    ENDIF

  END SUBROUTINE ComputeBoundConPSI_OUT_d

  
  !---------------------------------------------
  
  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeCofPSI_IN_d
    USE VAR_d, ONLY: sinX2, Xdif2, Xdif3
    USE VAR_PSI_IN_d
    IMPLICIT NONE
    INTEGER I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 2) .AND. (I <= (N1-1)) .AND. &
         (J >= 2) .AND. (J <= (N2-1)) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      AIM(I,J,K) = S1(I,J,K) / Xdif1(I)                                
      AIP(I,J,K) = S1(I+1,J,K) / Xdif1(I+1)
        
      AJM(I,J,K) = S2(I,J,K) / (X1(I) * Xdif2(J))              
      AJP(I,J,K) = S2(I,J+1,K) / (X1(I) * Xdif2(J+1))

      AKM(I,J,K) =  S3(I,J,K) / (X1(I) * sinX2(J) * Xdif3(K))
      AKP(I,J,K) =  S3(I,J,K+1) / (X1(I) * sinX2(J) * Xdif3(K+1))

      CON(I,J,K) = 0.0_8
      AP(I,J,K)  = 0.0_8

    ENDIF

  END SUBROUTINE ComputeCofPSI_IN_d

  !---------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeBoundCofPSI_IN_d
    USE VAR_d, ONLY: F1
    USE VAR_PSI_IN_d
    IMPLICIT NONE
    INTEGER I, J, K

    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= 2) .AND. (J <= (N2-1)) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

        AIP(N1-1,J,K) = 0.0_8
        !CON(N1-1,J,K) = CON(N1-1,J,K) + F1(2, J, K) * S1(N1,J,K)     
                
          
    ENDIF

  END SUBROUTINE ComputeBoundCofPSI_IN_d

  !---------------------------------------------

  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeApPSI_IN_d
    USE VAR_PSI_IN_d
    IMPLICIT NONE
    INTEGER I, J, K

    I = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    J = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y
    K = (BlockIdx%z - 1) * BlockDim%z + threadIdx%z

    IF ( (I >= 2) .AND. (I <= (N1-1)) .AND. &
         (J >= 2) .AND. (J <= (N2-1)) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN
         
      AP(I,J,K) = - AP(I,J,K) + AIM(I,J,K) + AIP(I,J,K) + &
                                AJM(I,J,K) + AJP(I,J,K) + &
                                AKM(I,J,K) + AKP(I,J,K)
         
    ENDIF

  END SUBROUTINE ComputeApPSI_IN_d
  
  !---------------------------------------------
  
  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeBoundConPSI_IN_d
    USE VAR_d, ONLY: F1
    USE VAR_PSI_IN_d
    IMPLICIT NONE
    INTEGER I, J, K

    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= 2) .AND. (J <= (N2-1)) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

        CON(N1-1,J,K) = F1(2, J, K) * S1(N1,J,K)                     
          
    ENDIF

  END SUBROUTINE ComputeBoundConPSI_IN_d
  
 
  !---------------------------------------------
  
  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeF2OutBound_d
    USE VAR_d
    USE VAR_PSI_OUT_d, ONLY: PSI_OUT, N1_PSI_OUT => N1
    IMPLICIT NONE
    INTEGER :: I, J, K

    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= 3) .AND. (J <= N2) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      F2(N1,J,K) = (PSI_OUT(N1_PSI_OUT-1, J, K) - PSI_OUT(N1_PSI_OUT-1, J-1, K)) * r_o**2  / ( Xdif2(J) * (r_o + Xdif1(N1)) ) 
    
    ENDIF

  END SUBROUTINE ComputeF2OutBound_d
  
  !---------------------------------------------
  
  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeF3OutBound_d
    USE VAR_d
    USE VAR_PSI_OUT_d, ONLY: PSI_OUT, N1_PSI_OUT => N1
    IMPLICIT NONE
    INTEGER :: I, J, K

    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= 2) .AND. (J <= N2) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN

      F3(N1,J,K) = (PSI_OUT(N1_PSI_OUT-1, J, K) - PSI_OUT(N1_PSI_OUT-1, J, K-1)) * r_o**2  / ( Xdif3(K) * (r_o + Xdif1(N1)))
      
    ENDIF

  END SUBROUTINE ComputeF3OutBound_d

  !---------------------------------------------
  
  ATTRIBUTES (GLOBAL) SUBROUTINE ComputeF2InBound_d
    USE VAR_d, ONLY: F2, N2, N3, Xdif2
    USE VAR_PSI_IN_d, ONLY: PSI_IN, N1_PSI_IN => N1
    IMPLICIT NONE
    INTEGER :: I, J, K

    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= 3) .AND. (J <= N2) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN
 
      F2(1,J,K) = (PSI_IN(N1_PSI_IN-1,J,K) - PSI_IN(N1_PSI_IN-1,J-1,K)) / Xdif2(J)
    
    ENDIF

  END SUBROUTINE ComputeF2InBound_d
  
  !---------------------------------------------
  
    ATTRIBUTES (GLOBAL) SUBROUTINE ComputeF3InBound_d
    USE VAR_d, ONLY: F3, N2, N3, Xdif3
    USE VAR_PSI_IN_d, ONLY: PSI_IN, N1_PSI_IN => N1
    IMPLICIT NONE
    INTEGER :: I, J, K

    J = (BlockIdx%x - 1) * BlockDim%x + threadIdx%x
    K = (BlockIdx%y - 1) * BlockDim%y + threadIdx%y

    IF ( (J >= 2) .AND. (J <= N2) .AND. &
         (K >= 2) .AND. (K <= (N3-1)) ) THEN
      
      F3(1,J,K) = (PSI_IN(N1_PSI_IN-1,J,K) - PSI_IN(N1_PSI_IN-1,J,K-1)) / Xdif3(K)
    
    ENDIF

  END SUBROUTINE ComputeF3InBound_d
  
  !---------------------------------------------

  SUBROUTINE ComputeSmaxB
  USE VAR_d
  USE VAR, ONLY: smax_b
  IMPLICIT NONE
  INTEGER I, J, K

    smax_b = 0.0_8

   !$cuf kernel do(3) <<<*, *>>>
    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 2, N1-1
          smax_b = max(smax_b, DABS(S1(I+1,J,K)*B1(I+1,J,K) - S1(I,J,K)*B1(I,J,K) + &
                               S2(I,J+1,K)*B2(I,J+1,K) - S2(I,J,K)*B2(I,J,K) + &
                               S3(I,J,K+1)*B3(I,J,K+1) - S3(I,J,K)*B3(I,J,K) ) / VCV(I,J,K))
        ENDDO          
      ENDDO
    ENDDO

  END SUBROUTINE ComputeSmaxB

  
  !---------------------------------------------
  
  
END MODULE MD
