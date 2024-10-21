PROGRAM MAIN
  USE VAR
  USE VAR_PSI_OUT, ONLY: PSI_OUT, N1_PSI_OUT => N1
  USE VAR_PSI_IN, ONLY: PSI_IN,  N1_PSI_IN => N1
  USE VAR_d, ONLY : F1_d => F1, F2_d => F2, F3_d => F3, &
                    B1_d => B1, B2_d => B2, B3_d => B3, &
                    U1_d => U1, U2_d => U2, U3_d => U3, &
                    T_d => T, P_d => P, PSI_C_d => PSI_C, &
                    E1_Oz_d => E1_Oz
  USE VAR_PSI_OUT_d, ONLY: PSI_OUT_d => PSI_OUT
  USE VAR_PSI_IN_d, ONLY: PSI_IN_d => PSI_IN
  USE FD, ONLY: SIMPLER
  USE MD, ONLY: CorrectB, MHD, ComputeSmaxB, ComputeMagPotential, ComputeBoundCondsF, FtoB
  USE GEOM
  USE USER
  USE OMP_LIB
  IMPLICIT NONE
  REAL(8) :: startTime, endTime
  INTEGER :: J, K
 
  !CALL OPEN_LOG_FILES  
  CALL ALLOC_MEMS

  startTime = omp_get_wtime()

  CALL SET_TASK_PARAMS
  CALL BUILD_GRID
!  CALL COMPUTE_DIFF_FLOW_SPH
!  CALL START
!  CALL COPY_CONST_ARRAYS_TO_DEVICE

 
!  CALL FF1N1(PSI_OUT, N1_PSI_OUT)
!  CALL FF1N1(PSI_IN, N1_PSI_IN)
!  CALL FF1N1(T, N1)
!  CALL FF1N1(U1, N1)
!  CALL FF1N1(U2, N1)
!  CALL FF1N1(U3, N1)

!  CALL FF1N1(F1, N1)
!  CALL FF1N1(F2, N1)
!  CALL FF1N1(F3, N1)

!  CALL FF1N1(B1, N1)
!  CALL FF1N1(B2, N1)
!  CALL FF1N1(B3, N1)

!  F1_d = F1; F2_d = F2; F3_d = F3
!  B1_d = B1; B2_d = B2; B3_d = B3
!  U1_d = U1; U2_d = U2; U3_d = U3
!  T_d = T; P_d = P; PSI_C_d = PSI_C
!  PSI_OUT_d = PSI_OUT; PSI_IN_d = PSI_IN 



!  CALL ComputeSmaxB
  
!  WRITE(*, *) 'Init SmaxB = ', SMAX_B
!  WRITE(10, *) 'Init SmaxB = ', SMAX_B
 
!  IF (SMAX_B > 1.0D-11) THEN
!    WRITE(*, *) 'Correction B ...'
!    WRITE(10, *) 'Correction B ...'
!    CALL CorrectB
!  END IF


!  CALL F0_F

!  DO TIME_ITER = 1, LAST
!    TIME = TIME + DT

!    CALL SIMPLER
            
!    IF (SMAX_B >= 1.0D-7) THEN
!      CALL CorrectB 
!    END IF

!    CALL F0_F

!    B1 = B1_d; B2 = B2_d; B3 = B3_d
!    F1 = F1_d; F2 = F2_d; F3 = F3_d
!    U1 = U1_d; U2 = U2_d; U3 = U3_d
!    T = T_d
!    E1_Oz = E1_Oz_d
!    CALL SAVE_DATA
!    CALL OUTPUT
!  ENDDO !TIME_ITER


  endTime = omp_get_wtime()
  WRITE(*,'(A, 1P1E16.4)') 'Elapsed time:', endTime - startTime
  WRITE(10,'(A, 1P1E16.4)') 'Elapsed time:', endTime - startTime

!  CALL OUTPUT_END
  CALL OUTPUT_GRID
  !CALL CLOSE_LOG_FILES
  CALL DEALLOC_MEMS

END PROGRAM MAIN

!------------------------------------------------------

SUBROUTINE F0_F
  USE VAR_d 
  IMPLICIT NONE
  
  F1_0 = F1; F2_0 = F2;  F3_0 = F3
  U1_0 = U1; U2_0 = U2;  U3_0 = U3
  T0 = T

END SUBROUTINE F0_F

!------------------------------------------------------

SUBROUTINE FF1N1(F, N1) 
  USE VAR, ONLY: N2, N3
  IMPLICIT NONE
  INTEGER :: I, J, N1
  REAL(8), DIMENSION (N1,N2,N3) :: F

  !$omp parallel do default(shared) private(I)
  DO J = 1, N2
    DO I = 1, N1
      F(I,J,1)  = F(I,J,N3-1)
      F(I,J,N3) = F(I,J,2)  
    ENDDO
  ENDDO

END SUBROUTINE FF1N1

!------------------------------------------------------

SUBROUTINE OPEN_LOG_FILES
  IMPLICIT NONE
  
  OPEN(UNIT=10,FILE='Q.OUT',STATUS='UNKNOWN', POSITION="APPEND", ACTION="WRITE")
  OPEN(UNIT=11,FILE='Smax.dat',STATUS='UNKNOWN', POSITION="APPEND", ACTION="WRITE")
  OPEN(UNIT=12,FILE='Energy.dat',STATUS='UNKNOWN', POSITION="APPEND", ACTION="WRITE")
  OPEN(UNIT=13,FILE='Bloc.dat',STATUS='UNKNOWN', POSITION="APPEND", ACTION="WRITE")
  OPEN(UNIT=14,FILE='divB.dat',STATUS='UNKNOWN', POSITION="APPEND", ACTION="WRITE")
  OPEN(UNIT=15,FILE='FluxB.dat',STATUS='UNKNOWN', POSITION="APPEND", ACTION="WRITE")
  OPEN(UNIT=16,FILE='UTloc.dat',STATUS='UNKNOWN', POSITION="APPEND", ACTION="WRITE")
  OPEN(UNIT=17,FILE='FluxQ.dat',STATUS='UNKNOWN', POSITION="APPEND", ACTION="WRITE")
  
END SUBROUTINE OPEN_LOG_FILES

!------------------------------------------------------

SUBROUTINE CLOSE_LOG_FILES
  IMPLICIT NONE
  
  ENDFILE 10
  CLOSE(10)

  ENDFILE 11
  CLOSE(11)

  ENDFILE 12
  CLOSE(12)

  ENDFILE 13
  CLOSE(13)   

  ENDFILE 14
  CLOSE(14)

  ENDFILE 15
  CLOSE(15)

  ENDFILE 16
  CLOSE(16)
  
  ENDFILE 17
  CLOSE(17)

  
END SUBROUTINE CLOSE_LOG_FILES

!------------------------------------------------------

SUBROUTINE ALLOC_MEMS
  IMPLICIT NONE

  CALL ALLOC_HOST_MEMS 
  CALL ALLOC_DEVICE_MEMS
  
END SUBROUTINE ALLOC_MEMS

!------------------------------------------------------

SUBROUTINE DEALLOC_MEMS
  IMPLICIT NONE

  CALL DEALLOC_HOST_MEMS
  CALL DEALLOC_DEVICE_MEMS

END SUBROUTINE DEALLOC_MEMS

!------------------------------------------------------

SUBROUTINE ALLOC_HOST_MEMS
  USE VAR
  USE VAR_PSI_OUT, ONLY: N1_PSI_OUT => N1, Xs1_PSI_OUT => Xs1, X1_PSI_OUT => X1, &
                         Xdif1_PSI_OUT => Xdif1, S1_PSI_OUT => S1, &
                         S2_PSI_OUT => S2, S3_PSI_OUT => S3, PSI_OUT
  USE VAR_PSI_IN, ONLY:  N1_PSI_IN => N1, Xs1_PSI_IN => Xs1, X1_PSI_IN => X1, &
                         Xdif1_PSI_IN => Xdif1, S1_PSI_IN => S1, &
                         S2_PSI_IN => S2, S3_PSI_IN => S3, PSI_IN
  IMPLICIT NONE
  INTEGER :: IERR
  
  ALLOCATE(Xs1(N1), X1(N1), Xdif1(N1), Kx1(N1), X_3(0:N3+1), &
           Xs1_PSI_OUT(N1_PSI_OUT), X1_PSI_OUT(N1_PSI_OUT), Xdif1_PSI_OUT(N1_PSI_OUT), &
           Xs1_PSI_IN(N1_PSI_IN), X1_PSI_IN(N1_PSI_IN), Xdif1_PSI_IN(N1_PSI_IN), &
           Xs2(N2), X2(N2), Xdif2(N2), Kx2(N2), &
           Xs3(N3), X3(N3), Xdif3(N3), Kx3(N3), &
          cosXs2(N2), cosX2(N2), sinXs2(N2), sinX2(N2), &
          X1cv(N1), X2cv(N2), X3cv(N3), E1_Oz(N1), &
          S1(N1, N2, N3), S2(N1, N2, N3), S3(N1, N2, N3), &
          S1_PSI_OUT(N1_PSI_OUT, N2, N3), S2_PSI_OUT(N1_PSI_OUT, N2, N3), &
          S3_PSI_OUT(N1_PSI_OUT, N2, N3), &
          S1_PSI_IN(N1_PSI_IN, N2, N3), S2_PSI_IN(N1_PSI_IN, N2, N3), &
          S3_PSI_IN(N1_PSI_IN, N2, N3), &
          S1u1(N1, N2, N3), S2u1(N1, N2, N3), S2u1_1(N1, N2, N3), S3u1(N1, N2, N3), &
          S1u2(N1, N2, N3), S2u2(N1, N2, N3), S3u2(N1, N2, N3), &
          S1u3(N1, N2, N3), S2u3(N1, N2, N3), &
          U1(N1, N2, N3), U2(N1, N2, N3), U3(N1, N2, N3), &
          F1(N1, N2, N3), F2(N1, N2, N3), F3(N1, N2, N3), &
          F1_0(N1, N2, N3), F2_0(N1, N2, N3), F3_0(N1, N2, N3), &
          F1a(N1, N2, N3), F2a(N1, N2, N3), F3a(N1, N2, N3), &
          AIM(N1, N2, N3), AIP(N1, N2, N3), AJM(N1, N2, N3), &
          AJP(N1, N2, N3), AKM(N1, N2, N3), AKP(N1, N2, N3), &
          AP(N1, N2, N3), CON(N1, N2, N3), T(N1, N2, N3), T0(N1, N2, N3), &
          D213(N1,N2,N3),  D231(N1,N2,N3), D312(N1,N2,N3), D321(N1,N2,N3), &
          D123(N1,N2,N3), D132(N1,N2,N3), F_123(N1,N2,N3), F_132(N1,N2,N3), & 
          F_213(N1,N2,N3), F_231(N1,N2,N3), F_312(N1,N2,N3), F_321(N1,N2,N3), &
          J213(N1,N2,N3), J231(N1,N2,N3), J312(N1,N2,N3), J321(N1,N2,N3), &
          J123(N1,N2,N3), J132(N1,N2,N3), & 
          B1(N1,N2,N3), B2(N1,N2,N3), B3(N1,N2,N3), &
          B2_o(N2,N3), B2_i(N2,N3), B3_o(N2,N3), B3_i(N2,N3), &
          PSI(N1,N2,N3), PSI_OUT(N1_PSI_OUT,N2,N3), PSI_IN(N1_PSI_IN,N2,N3), &
          B1a(N1,N2,N3), B2a(N1,N2,N3), B3a(N1,N2,N3), TEMPS(N1,N2,N3), &
          U1a(N1,N2,N3), U2a(N1,N2,N3), U3a(N1,N2,N3), &
          A1(N1,N2,N3), A2(N1,N2,N3), A3(N1,N2,N3), &
          J1(N1,N2,N3), J2(N1,N2,N3), J3(N1,N2,N3), &
          VCV(N1,N2,N3), VCV1(N1,N2,N3), VCV2(N1,N2,N3), VCV3(N1,N2,N3), PSI_C(N1,N2,N3), &
          P(N1, N2, N3), PC(N1, N2, N3), U1_0(N1, N2, N3), U2_0(N1, N2, N3), U3_0(N1, N2, N3), &
          U1hat(N1, N2, N3), U2hat(N1, N2, N3), U3hat(N1, N2, N3), &
          DU1(N1, N2, N3), DU2(N1, N2, N3), DU3(N1, N2, N3), &
          VVY1(N2), VVY2(N2), VVX1(N1), VVX2(N1), VYX1(N1), &
          VYX2(N1), APX(N1), &
          WWY1(N2), WWY2(N2), UYY(N2), WZY1(N2), WZY2(N2), &
          Spv1J(N2), Spvy1(N2), Spvy2(N2), &
          VZY1(N2), VZY2(N2), Spwy1(N2), Spwy2(N2), &
          VVX1S(N1), VVX2S(N1), XsDifP4(N1), XDifP4(N1), &
          WWY1S(N2), WWY2S(N2), VZY1S(N2), VZY2S(N2), STAT = IERR)

END SUBROUTINE ALLOC_HOST_MEMS

!------------------------------------------------------

SUBROUTINE ALLOC_DEVICE_MEMS
  USE VAR_d
  USE VAR_PSI_OUT_d, ONLY: N1_PSI_OUT => N1, Xs1_PSI_OUT => Xs1, X1_PSI_OUT => X1, &
                           Xdif1_PSI_OUT => Xdif1, S1_PSI_OUT => S1, &
                           S2_PSI_OUT => S2, S3_PSI_OUT => S3, PSI_OUT, &
                           TEMPS_PSI_OUT => TEMPS, AIM_PSI_OUT => AIM, &
                           AIP_PSI_OUT => AIP, AJM_PSI_OUT => AJM, &
                           AJP_PSI_OUT => AJP, AKM_PSI_OUT => AKM, &
                           AKP_PSI_OUT => AKP, AP_PSI_OUT => AP, &
                           CON_PSI_OUT => CON
  USE VAR_PSI_IN_d, ONLY:  N1_PSI_IN => N1, Xs1_PSI_IN => Xs1, X1_PSI_IN => X1, &
                           Xdif1_PSI_IN => Xdif1, S1_PSI_IN => S1, &
                           S2_PSI_IN => S2, S3_PSI_IN => S3, PSI_IN, &
                           TEMPS_PSI_IN => TEMPS, AIM_PSI_IN => AIM, &
                           AIP_PSI_IN => AIP, AJM_PSI_IN => AJM, &
                           AJP_PSI_IN => AJP, AKM_PSI_IN => AKM, &
                           AKP_PSI_IN => AKP, AP_PSI_IN => AP, &
                           CON_PSI_IN => CON
  IMPLICIT NONE
  INTEGER :: IERR

  WRITE(*, *) 'N1_PSI_OUT', N1_PSI_OUT
  
  ALLOCATE(Xs1(N1), X1(N1), Xdif1(N1), Kx1(N1), X_3(0:N3+1), &
           Xs1_PSI_OUT(N1_PSI_OUT), X1_PSI_OUT(N1_PSI_OUT), Xdif1_PSI_OUT(N1_PSI_OUT), &
           Xs1_PSI_IN(N1_PSI_IN), X1_PSI_IN(N1_PSI_IN), Xdif1_PSI_IN(N1_PSI_IN), &
          Xs2(N2), X2(N2), Xdif2(N2), Kx2(N2), &
          Xs3(N3), X3(N3), Xdif3(N3), Kx3(N3), &
          cosXs2(N2), cosX2(N2), sinXs2(N2), sinX2(N2), &
          X1cv(N1), X2cv(N2), X3cv(N3), E1_Oz(N1), &
          S1(N1, N2, N3), S2(N1, N2, N3), S3(N1, N2, N3), &
          S1_PSI_OUT(N1_PSI_OUT, N2, N3), S2_PSI_OUT(N1_PSI_OUT, N2, N3), &
          S3_PSI_OUT(N1_PSI_OUT, N2, N3), &
          S1_PSI_IN(N1_PSI_IN, N2, N3), S2_PSI_IN(N1_PSI_IN, N2, N3), &
          S3_PSI_IN(N1_PSI_IN, N2, N3), &
          S1u1(N1, N2, N3), S2u1(N1, N2, N3), S2u1_1(N1, N2, N3), S3u1(N1, N2, N3), &
          S1u2(N1, N2, N3), S2u2(N1, N2, N3), S3u2(N1, N2, N3), &
          S1u3(N1, N2, N3), S2u3(N1, N2, N3), &
          U1(N1, N2, N3), U2(N1, N2, N3), U3(N1, N2, N3), &
          U1tmp(N1, N2, N3), U2tmp(N1, N2, N3), U3tmp(N1, N2, N3), &
          F1(N1, N2, N3), F2(N1, N2, N3), F3(N1, N2, N3), &
          F1t(N1, N2, N3), F2t(N1, N2, N3), F3t(N1, N2, N3), &
          F1_0(N1, N2, N3), F2_0(N1, N2, N3), F3_0(N1, N2, N3), &
          F_1(1:N1, 1:N2, 0:N3+1), F_2(1:N1, 1:N2, 0:N3+1), & 
          AIM(N1, N2, N3), AIP(N1, N2, N3), AJM(N1, N2, N3), &
          AJP(N1, N2, N3), AKM(N1, N2, N3), AKP(N1, N2, N3), &
          AP(N1, N2, N3), CON(N1, N2, N3), T(N1, N2, N3), T0(N1, N2, N3), &
          AIM_PSI_OUT(N1_PSI_OUT, N2, N3), AIP_PSI_OUT(N1_PSI_OUT, N2, N3), &
          AJM_PSI_OUT(N1_PSI_OUT, N2, N3), AJP_PSI_OUT(N1_PSI_OUT, N2, N3), &
          AKM_PSI_OUT(N1_PSI_OUT, N2, N3), AKP_PSI_OUT(N1_PSI_OUT, N2, N3), &
          AP_PSI_OUT(N1_PSI_OUT, N2, N3), CON_PSI_OUT(N1_PSI_OUT, N2, N3), &
          TEMPS_PSI_OUT(N1_PSI_OUT, N2, N3), &
          AIM_PSI_IN(N1_PSI_IN, N2, N3), AIP_PSI_IN(N1_PSI_IN, N2, N3), &
          AJM_PSI_IN(N1_PSI_IN, N2, N3), AJP_PSI_IN(N1_PSI_IN, N2, N3), &
          AKM_PSI_IN(N1_PSI_IN, N2, N3), AKP_PSI_IN(N1_PSI_IN, N2, N3), &
          AP_PSI_IN(N1_PSI_IN, N2, N3), CON_PSI_IN(N1_PSI_IN, N2, N3), &
          TEMPS_PSI_IN(N1_PSI_IN, N2, N3), &
          D213(N1,N2,N3),  D231(N1,N2,N3), D312(N1,N2,N3), D321(N1,N2,N3), &
          D123(N1,N2,N3), D132(N1,N2,N3), F_123(N1,N2,N3), F_132(N1,N2,N3), &
          J213(N1,N2,N3), J231(N1,N2,N3), J312(N1,N2,N3), J321(N1,N2,N3), &
          J123(N1,N2,N3), J132(N1,N2,N3), &  
          F_213(N1,N2,N3), F_231(N1,N2,N3), F_312(N1,N2,N3), F_321(N1,N2,N3), &
          B1(N1,N2,N3), B2(N1,N2,N3), B3(N1,N2,N3), &
          PSI(N1,N2,N3), PSI_OUT(N1_PSI_OUT,N2,N3), PSI_IN(N1_PSI_IN,N2,N3), &
          TEMPS(N1,N2,N3), &
          B2_o(N2,N3), B2_i(N2,N3), B3_o(N2,N3), B3_i(N2,N3), &
          A1(N1,N2,N3), A2(N1,N2,N3), A3(N1,N2,N3), &
          J1(N1,N2,N3), J2(N1,N2,N3), J3(N1,N2,N3), &
          VCV(N1,N2,N3), VCV1(N1,N2,N3), VCV2(N1,N2,N3), VCV3(N1,N2,N3), PSI_C(N1,N2,N3), &
          P(N1, N2, N3), PC(N1, N2, N3), U1_0(N1, N2, N3), U2_0(N1, N2, N3), U3_0(N1, N2, N3), &
          U1hat(N1, N2, N3), U2hat(N1, N2, N3), U3hat(N1, N2, N3), &
          DU1(N1, N2, N3), DU2(N1, N2, N3), DU3(N1, N2, N3), &
          VVY1(N2), VVY2(N2), VVX1(N1), VVX2(N1), VYX1(N1), &
          VYX2(N1), APX(N1), &
          WWY1(N2), WWY2(N2), UYY(N2), WZY1(N2), WZY2(N2), &
          Spv1J(N2), Spvy1(N2), Spvy2(N2), &
          VZY1(N2), VZY2(N2), Spwy1(N2), Spwy2(N2), &
          VVX1S(N1), VVX2S(N1), XsDifP4(N1), XDifP4(N1), &
          WWY1S(N2), WWY2S(N2), VZY1S(N2), VZY2S(N2), STAT = IERR)

END SUBROUTINE ALLOC_DEVICE_MEMS

!------------------------------------------------------

SUBROUTINE DEALLOC_HOST_MEMS
  USE VAR
  USE VAR_PSI_OUT, ONLY: Xs1_PSI_OUT => Xs1, X1_PSI_OUT => X1, &
                         Xdif1_PSI_OUT => Xdif1, S1_PSI_OUT => S1, &
                         S2_PSI_OUT => S2, S3_PSI_OUT => S3, PSI_OUT
  USE VAR_PSI_IN, ONLY:  Xs1_PSI_IN => Xs1, X1_PSI_IN => X1, &
                         Xdif1_PSI_IN => Xdif1, S1_PSI_IN => S1, &
                         S2_PSI_IN => S2, S3_PSI_IN => S3, PSI_IN
  IMPLICIT NONE

  INTEGER :: IERR
  
  DEALLOCATE(Xs1, X1, Xdif1, Kx1, X_3, &
             Xs1_PSI_OUT, X1_PSI_OUT, Xdif1_PSI_OUT, &
             Xs1_PSI_IN, X1_PSI_IN, Xdif1_PSI_IN, &
             Xs2, X2, Xdif2, Kx2, &
             Xs3, X3, Xdif3, Kx3, &
             cosXs2, cosX2, sinXs2, sinX2, &
             X1cv, X2cv, X3cv, E1_Oz, &
             S1, S2, S3, S1u1, S2u1, S2u1_1, S3u1, &
             S1_PSI_OUT, S2_PSI_OUT, S3_PSI_OUT, &
             S1_PSI_IN, S2_PSI_IN, S3_PSI_IN, &
             S1u2, S2u2, S3u2, &
             S1u3, S2u3, &
             U1, U2, U3, F1, F2, F3, T, T0, P, PC, &
             U1a, U2a, U3a, &
             B1a, B2a, B3a, &
             F1a, F2a, F3a, &
             F1_0, F2_0, F3_0, &
             B1, B2, B3, PSI, &
             B2_o, B2_i, B3_o, B3_i, &
             PSI_OUT, PSI_IN, &
             TEMPS, &
             A1, A2, A3, J1, J2, J3, &
             AIM, AIP, AJM, AJP, AKM, AKP, AP, CON, &
             VCV, VCV1, VCV2, VCV3, &
             D213, D231, D312, D321, D123, D132, &
             F_213, F_231, F_312, F_321, F_123, F_132, PSI_C, &
             J213, J231, J312, J321, J123, J132, &
             U1_0, U2_0, U3_0, U1hat, U2hat, U3hat, DU1, DU2, DU3, &
             VVY1, VVY2, VVX1, VVX2, VYX1, VYX2, APX, &
             WWY1, WWY2, UYY, WZY1, WZY2, Spv1J, Spvy1, Spvy2, &
             VZY1, VZY2, Spwy1, Spwy2, &
             VVX1S, VVX2S, XsDifP4, XDifP4, WWY1S, WWY2S, &
             VZY1S, VZY2S, STAT = IERR)

END SUBROUTINE DEALLOC_HOST_MEMS

!------------------------------------------------------

SUBROUTINE DEALLOC_DEVICE_MEMS
  USE VAR_d
  USE VAR_PSI_OUT_d, ONLY: Xs1_PSI_OUT => Xs1, X1_PSI_OUT => X1, &
                           Xdif1_PSI_OUT => Xdif1, S1_PSI_OUT => S1, &
                           S2_PSI_OUT => S2, S3_PSI_OUT => S3, PSI_OUT, &
                           TEMPS_PSI_OUT => TEMPS, AIM_PSI_OUT => AIM, &
                           AIP_PSI_OUT => AIP, AJM_PSI_OUT => AJM, &
                           AJP_PSI_OUT => AJP, AKM_PSI_OUT => AKM, &
                           AKP_PSI_OUT => AKP, AP_PSI_OUT => AP, &
                           CON_PSI_OUT => CON
  USE VAR_PSI_IN_d, ONLY:  Xs1_PSI_IN => Xs1, X1_PSI_IN => X1, &
                           Xdif1_PSI_IN => Xdif1, S1_PSI_IN => S1, &
                           S2_PSI_IN => S2, S3_PSI_IN => S3, PSI_IN, &
                           TEMPS_PSI_IN => TEMPS, AIM_PSI_IN => AIM, &
                           AIP_PSI_IN => AIP, AJM_PSI_IN => AJM, &
                           AJP_PSI_IN => AJP, AKM_PSI_IN => AKM, &
                           AKP_PSI_IN => AKP, AP_PSI_IN => AP, &
                           CON_PSI_IN => CON
  IMPLICIT NONE

  INTEGER :: IERR
  
  DEALLOCATE(Xs1, X1, Xdif1, Kx1, X_3, &
             Xs1_PSI_OUT, X1_PSI_OUT, Xdif1_PSI_OUT, &
             Xs1_PSI_IN, X1_PSI_IN, Xdif1_PSI_IN, &
             Xs2, X2, Xdif2, Kx2, &
             Xs3, X3, Xdif3, Kx3, &
             cosXs2, cosX2, sinXs2, sinX2, &
             X1cv, X2cv, X3cv, E1_Oz, &
             S1, S2, S3, S1u1, S2u1, S2u1_1, S3u1, &
             S1_PSI_OUT, S2_PSI_OUT, S3_PSI_OUT, &
             S1_PSI_IN, S2_PSI_IN, S3_PSI_IN, &
             S1u2, S2u2, S3u2, &
             S1u3, S2u3, &
             U1, U2, U3, U1tmp, U2tmp, U3tmp, F1, F2, F3, T, T0, P, PC, &
             F1t, F2t, F3t, & 
             F1_0, F2_0, F3_0, &
             F_1, &
             B1, B2, B3, PSI, &
             PSI_OUT, PSI_IN, &
             TEMPS, &
             B2_o, B2_i, B3_o, B3_i, &
             TEMPS_PSI_OUT, TEMPS_PSI_IN, &
             A1, A2, A3, J1, J2, J3, &
             AIM, AIP, AJM, AJP, AKM, AKP, AP, CON, &
             AIM_PSI_OUT, AIP_PSI_OUT, AJM_PSI_OUT, AJP_PSI_OUT, &
             AKM_PSI_OUT, AKP_PSI_OUT, AP_PSI_OUT, CON_PSI_OUT, &
             AIM_PSI_IN, AIP_PSI_IN, AJM_PSI_IN, AJP_PSI_IN, &
             AKM_PSI_IN, AKP_PSI_IN, AP_PSI_IN, CON_PSI_IN, &
             VCV, VCV1, VCV2, VCV3, &
             D213, D231, D312, D321, D123, D132, &
             F_213, F_231, F_312, F_321, F_123, F_132, PSI_C, &
             J213, J231, J312, J321, J123, J132, &
             U1_0, U2_0, U3_0, U1hat, U2hat, U3hat, DU1, DU2, DU3, &
             VVY1, VVY2, VVX1, VVX2, VYX1, VYX2, APX, &
             WWY1, WWY2, UYY, WZY1, WZY2, Spv1J, Spvy1, Spvy2, &
             VZY1, VZY2, Spwy1, Spwy2, &
             VVX1S, VVX2S, XsDifP4, XDifP4, WWY1S, WWY2S, &
             VZY1S, VZY2S, STAT = IERR)

END SUBROUTINE DEALLOC_DEVICE_MEMS

!------------------------------------------------------

SUBROUTINE COPY_CONST_ARRAYS_TO_DEVICE
  USE VAR 
  USE VAR_PSI_OUT, ONLY: Xs1_PSI_OUT => Xs1, X1_PSI_OUT => X1, &
                         Xdif1_PSI_OUT => Xdif1, S1_PSI_OUT => S1, &
                         S2_PSI_OUT => S2, S3_PSI_OUT => S3
  USE VAR_PSI_IN, ONLY:  Xs1_PSI_IN => Xs1, X1_PSI_IN => X1, &
                         Xdif1_PSI_IN => Xdif1, S1_PSI_IN => S1, &
                         S2_PSI_IN => S2, S3_PSI_IN => S3
  USE VAR_d, ONLY : Xs1_d => Xs1, X1_d => X1, Xdif1_d => Xdif1, &
                    Kx1_d => Kx1, Xs2_d => Xs2, X2_d => X2, &  
                    Xdif2_d => Xdif2, Kx2_d => Kx2, Xs3_d => Xs3, & 
                    X3_d => X3, Xdif3_d => Xdif3, Kx3_d => Kx3, X_3_d => X_3, & 
                    cosXs2_d => cosXs2, cosX2_d => cosX2, sinXs2_d => sinXs2, &  
                    sinX2_d => sinX2, X1cv_d => X1cv, X2cv_d => X2cv, & 
                    X3cv_d => X3cv, S1_d => S1, S2_d => S2, S3_d => S3, &  
                    S1u1_d => S1u1, S2u1_d => S2u1, S2u1_1_d => S2u1_1, &  
                    S3u1_d => S3u1, S1u2_d => S1u2, S2u2_d => S2u2, &  
                    S3u2_d => S3u2, S1u3_d => S1u3, S2u3_d => S2u3, & 
                    VCV_d => VCV, VCV1_d => VCV1, VCV2_d => VCV2, &  
                    VCV3_d => VCV3, D213_d => D213, D231_d => D231, &  
                    D312_d => D312, D321_d => D321, D123_d => D123, &  
                    D132_d => D132, F_213_d => F_213, F_231_d => F_231, &
                    F_312_d => F_312, F_321_d => F_321, F_123_d => F_123, &  
                    F_132_d => F_132, &
                    J312_d => J312, J321_d => J321, J123_d => J123, &  
                    J132_d => J132, J213_d => J213, J231_d => J231, &
                    DT_d => DT, &
                    viscm_d => viscm, Ranum_d => Ranum, Pm_d => Pm, &
                    thc_d => thc, rho_c_d => rho_c, r_i_d => r_i, &
                    r_o_d => r_o, Pr_d => Pr, &
                    E_d => E, visk_d => visk, dens_d => dens, &
                    VVY1_d => VVY1, VVY2_d => VVY2, VVX1_d => VVX1, &
                    VVX2_d => VVX2, VYX1_d => VYX1, VYX2_d => VYX2, &
                    APX_d => APX, &
                    WWY1_d => WWY1, WWY2_d => WWY2, UYY_d => UYY, &
                    WZY1_d => WZY1, WZY2_d => WZY2, Spv1J_d => Spv1J, &
                    Spvy1_d => Spvy1, Spvy2_d => Spvy2, &
                    VZY1_d => VZY1, VZY2_d => VZY2, Spwy1_d => Spwy1, Spwy2_d => Spwy2, &
                    VVX1S_d => VVX1S, VVX2S_d => VVX2S, XsDifP4_d => XsDifP4, &
                    XDifP4_d => XDifP4, WWY1S_d => WWY1S, WWY2S_d => WWY2S, &
                    VZY1S_d => VZY1S, VZY2S_d => VZY2S                
  USE VAR_PSI_OUT_d, ONLY: Xs1_PSI_OUT_d => Xs1, X1_PSI_OUT_d => X1, &
                           Xdif1_PSI_OUT_d => Xdif1, S1_PSI_OUT_d => S1, &
                           S2_PSI_OUT_d => S2, S3_PSI_OUT_d => S3
  USE VAR_PSI_IN_d, ONLY:  Xs1_PSI_IN_d => Xs1, X1_PSI_IN_d => X1, &
                           Xdif1_PSI_IN_d => Xdif1, S1_PSI_IN_d => S1, &
                           S2_PSI_IN_d => S2, S3_PSI_IN_d => S3
  IMPLICIT NONE

  Xs1_d = Xs1; X1_d = X1; Xdif1_d = Xdif1; Kx1_d = Kx1
  Xs1_PSI_OUT_d = Xs1_PSI_OUT; X1_PSI_OUT_d = X1_PSI_OUT; Xdif1_PSI_OUT_d = Xdif1_PSI_OUT
  Xs1_PSI_IN_d = Xs1_PSI_IN; X1_PSI_IN_d = X1_PSI_IN; Xdif1_PSI_IN_d = Xdif1_PSI_IN
  Xs2_d = Xs2; X2_d = X2; Xdif2_d = Xdif2; Kx2_d = Kx2 
  Xs3_d = Xs3; X3_d = X3; X_3_d = X_3; Xdif3_d = Xdif3; Kx3_d = Kx3
  cosXs2_d = cosXs2; cosX2_d = cosX2; sinXs2_d = sinXs2  
  sinX2_d = sinX2; X1cv_d = X1cv; X2cv_d = X2cv
  X3cv_d = X3cv; S1_d = S1; S2_d = S2; S3_d = S3  
  S1_PSI_OUT_d = S1_PSI_OUT; S2_PSI_OUT_d = S2_PSI_OUT; S3_PSI_OUT_d = S3_PSI_OUT
  S1_PSI_IN_d = S1_PSI_IN; S2_PSI_IN_d = S2_PSI_IN; S3_PSI_IN_d = S3_PSI_IN
  S1u1_d = S1u1; S2u1_d = S2u1; S2u1_1_d = S2u1_1  
  S3u1_d = S3u1; S1u2_d = S1u2; S2u2_d = S2u2  
  S3u2_d = S3u2; S1u3_d = S1u3; S2u3_d = S2u3 
  VCV_d = VCV; VCV1_d = VCV1; VCV2_d = VCV2  
  VCV3_d = VCV3; D213_d = D213; D231_d = D231  
  D312_d = D312; D321_d = D321; D123_d = D123  
  D132_d = D132; F_213_d = F_213; F_231_d = F_231  
  F_312_d = F_312; F_321_d = F_321; F_123_d = F_123  
  F_132_d = F_132; DT_d = DT
  J213_d = J213; J231_d = J231  
  J312_d = J312; J321_d = J321; J123_d = J123  
  J132_d = J132;
  viscm_d = viscm; Ranum_d = Ranum; Pm_d = Pm; thc_d = thc
  rho_c_d = rho_c; r_i_d = r_i; r_o_d = r_o; Pr_d = Pr
  E_d = E; visk_d = visk; dens_d = dens
  VVY1_d = VVY1; VVY2_d = VVY2; VVX1_d = VVX1
  VVX2_d = VVX2; VYX1_d = VYX1; VYX2_d = VYX2
  APX_d = APX
  WWY1_d = WWY1; WWY2_d = WWY2; UYY_d = UYY
  WZY1_d = WZY1; WZY2_d = WZY2; Spv1J_d = Spv1J
  Spvy1_d = Spvy1; Spvy2_d = Spvy2
  VZY1_d = VZY1; VZY2_d = VZY2; Spwy1_d = Spwy1
  Spwy2_d = Spwy2
  VVX1S_d = VVX1S; VVX2S_d = VVX2S; XsDifP4_d = XsDifP4
  XDifP4_d = XDifP4; WWY1S_d = WWY1S; WWY2S_d = WWY2S
  VZY1S_d = VZY1S; VZY2S_d = VZY2S

END SUBROUTINE COPY_CONST_ARRAYS_TO_DEVICE

!------------------------------------------------------
