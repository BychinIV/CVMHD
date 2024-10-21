MODULE VAR
  IMPLICIT NONE 
  REAL(8), PARAMETER :: PI = 4.0_8 * DATAN(1.0_8) 
  INTEGER :: LAST, IWRITE, TIME_ITER
  REAL(8) :: TIME, DT
  INTEGER, PARAMETER :: N1cv_core = 20, N1cv_shell = 50, N1 = N1cv_shell + N1cv_core + 2, N1_core = N1cv_core + 2,  &
                        N2 = 72, N3 = 142, Nmax = 142
  INTEGER            :: NTIMES_P, NTIMES_PC, NTIMES_T, NTIMES_UVW, NTIMES_F1, NTIMES_F2, NTIMES_F3
  LOGICAL            :: FirstCallMP
  REAL(8) :: XL1, XL2, XL3, X1_0, X2_0, X3_0, viscm, tolPSI, tolPSI_C, SMAX_B, sum_sq_f, &
             Pm, Br_out_avg, Br_in_avg, &
             thc, rho_c, r_i, r_o, Ranum, Pr, E, visk, dens, SMAX, TOL_P, TOL_T, TOL_UVW, Ha, B0, Umax, omega_ic, &
             resF1, resF2, resF3  
  REAL(8), DIMENSION (:),  ALLOCATABLE :: Xs1, X1, Xdif1, Kx1, &
                                          Xs2, X2, Xdif2, Kx2, &
                                          Xs3, X3, Xdif3, Kx3, X1_out, &
                                          cosXs2, cosX2, sinXs2, sinX2, &
                                          X1cv, X2cv, X3cv, &
                                          VVY1, VVY2, VVX1, VVX2, VYX1, VYX2, APX, &
                                          WWY1, WWY2, UYY, WZY1, WZY2, Spv1J, Spvy1, Spvy2, &
                                          VZY1, VZY2, Spwy1, Spwy2, VVX1S, VVX2S, XsDifP4, &
                                          XDifP4, WWY1S, WWY2S, VZY1S, VZY2S, E1_Oz 
  REAL(8), DIMENSION (:), ALLOCATABLE :: X_3
  REAL(8), DIMENSION (:,:,:),  ALLOCATABLE :: U1, U2, U3, T, T0, P, PC, F1, F2, F3, S1, S2, S3, &
                                              U1_0, U2_0, U3_0, U1hat, U2hat, U3hat, DU1, DU2, DU3, &
                                              S1u1, S2u1, S2u1_1, S3u1, &
                                              S1u2, S2u2, S3u2, &
                                              S1u3, S2u3, &
                                              B1, B2, B3, PSI, A1, A2, A3, &
                                              VCV, VCV1, VCV2, VCV3, &
                                              F1_0, F2_0, F3_0, F1a, F2a, F3a, &
                                              U1a, U2a, U3a, B1a, B2a, B3a, TEMPS, &
                                              AIM, AIP, AJM, AJP, AKM, AKP, AP, CON, &
                                              D213,  D231, D312, D321, D123, D132,&
                                              F_213, F_231, F_312, F_321, F_123, F_132, &
                                              J213,  J231, J312, J321, J123, J132, &
                                              J1, J2, J3, PSI_C
  REAL(8), DIMENSION (:,:),  ALLOCATABLE :: B2_o(:,:), B2_i(:,:), B3_o(:,:), B3_i(:,:) 

END MODULE VAR

MODULE VAR_PSI_OUT
  IMPLICIT NONE 
  INTEGER, PARAMETER :: N1 = 52, N2 = 72, N3 = 142, Nmax = 142
  REAL(8), DIMENSION (:),  ALLOCATABLE :: Xs1, X1, Xdif1
  REAL(8), DIMENSION (:,:,:), ALLOCATABLE :: PSI_OUT, S1, S2, S3
END MODULE VAR_PSI_OUT


MODULE VAR_PSI_IN
  IMPLICIT NONE 
  INTEGER, PARAMETER :: N1 = 32, N2 = 72, N3 = 142, Nmax = 142
  REAL(8), DIMENSION (:),  ALLOCATABLE :: Xs1, X1, Xdif1
  REAL(8), DIMENSION (:,:,:), ALLOCATABLE :: PSI_IN, S1, S2, S3
END MODULE VAR_PSI_IN

!===============================================================================

MODULE VAR_d
  USE CUDAFOR
  USE VAR, ONLY: N1_h => N1, N2_h => N2, N3_h => N3, Nmax_h => Nmax, N1_core_h => N1_core
  IMPLICIT NONE 
  INTEGER, PARAMETER :: N1 = N1_h, N2 = N2_h, N3 = N3_h, Nmax = Nmax_h, N1_core = N1_core_h
  REAL(8), PARAMETER :: PI = 4.0_8 * DATAN(1.0_8) 
  TYPE(DIM3) :: GRD_IJK, GRD_IJ, GRD_JK, GRD_IK, GRD_I, TBL_IJK, TBL_IJ, TBL_JK, TBL_IK, TBL_I 
  REAL(8), DEVICE :: DT, Br_in_avg, Br_out_avg, Btheta_avg, Bn_avg, sum_sq_f
  REAL(8), DEVICE :: viscm, Ranum, Pm, thc, rho_c, r_i, r_o, Pr, E, visk, dens
  REAL(8), DEVICE, DIMENSION (:),  ALLOCATABLE :: Xs1, X1, Xdif1, Kx1, &
                                                  Xs2, X2, Xdif2, Kx2, &
                                                  Xs3, X3, Xdif3, Kx3, X1_out, &
                                                  cosXs2, cosX2, sinXs2, sinX2, &
                                                  X1cv, X2cv, X3cv, &
                                                  VVY1, VVY2, VVX1, VVX2, &
                                                  VYX1, VYX2, APX, &
                                                  WWY1, WWY2, UYY, WZY1, &
                                                  WZY2, Spv1J, Spvy1, Spvy2, &
                                                  VZY1, VZY2, Spwy1, Spwy2, &
                                                  VVX1S, VVX2S, XsDifP4, XDifP4, &
                                                  WWY1S, WWY2S, VZY1S, VZY2S, E1_Oz

  REAL(8), DEVICE, DIMENSION (:),  ALLOCATABLE :: X_3
  REAL(8), DEVICE, DIMENSION (:,:,:),  ALLOCATABLE :: U1, U2, U3, U1tmp, U2tmp, U3tmp, T, T0, P, PC, &
                                                      F1, F2, F3, S1, S2, S3, &
                                                      F1t, F2t, F3t, & 
                                                      F_1, F_2, &
                                                      U1_0, U2_0, U3_0, U1hat, U2hat, &
                                                      U3hat, DU1, DU2, DU3, &
                                                      S1u1, S2u1, S2u1_1, S3u1, &
                                                      S1u2, S2u2, S3u2, &
                                                      S1u3, S2u3, &
                                                      B1, B2, B3, PSI, PSI_in, A1, A2, A3, &
                                                      VCV, VCV1, VCV2, VCV3, &
                                                      F1_0, F2_0, F3_0, &
                                                      TEMPS, &
                                                      AIM, AIP, AJM, AJP, AKM, AKP, &
                                                      AP, CON, D213, D231, D312, D321, &
                                                      D123, D132, F_213, F_231, F_312, &
                                                      J213,  J231, J312, J321, J123, J132, &
                                                      F_321, F_123, F_132, J1, J2, J3, PSI_C
  REAL(8), DEVICE, DIMENSION (:,:),  ALLOCATABLE ::   B2_o(:,:), B2_i(:,:), B3_o(:,:), B3_i(:,:) 
END MODULE VAR_d


MODULE VAR_PSI_OUT_d
  USE CUDAFOR
  USE VAR_PSI_OUT, ONLY: N1_h => N1, N2_h => N2, N3_h => N3, Nmax_h => Nmax
  IMPLICIT NONE 
  INTEGER, PARAMETER :: N1 = N1_h, N2 = N2_h, N3 = N3_h, Nmax = Nmax_h
  TYPE(DIM3) :: GRD_IJK, GRD_IJ, GRD_JK, GRD_IK, TBL_IJK, TBL_IJ, TBL_JK, TBL_IK 
  REAL(8), DEVICE :: sum_sq_f
  REAL(8), DEVICE, DIMENSION (:),  ALLOCATABLE :: Xs1, X1, Xdif1
  REAL(8), DEVICE, DIMENSION (:,:,:),  ALLOCATABLE :: PSI_OUT, S1, S2, S3, &
                                                      TEMPS, AIM, AIP, AJM, AJP, AKM, AKP, &
                                                      AP, CON
END MODULE VAR_PSI_OUT_d


MODULE VAR_PSI_IN_d
  USE CUDAFOR
  USE VAR_PSI_IN, ONLY: N1_h => N1, N2_h => N2, N3_h => N3, Nmax_h => Nmax
  IMPLICIT NONE 
  INTEGER, PARAMETER :: N1 = N1_h, N2 = N2_h, N3 = N3_h, Nmax = Nmax_h
  TYPE(DIM3) :: GRD_IJK, GRD_IJ, GRD_JK, GRD_IK, TBL_IJK, TBL_IJ, TBL_JK, TBL_IK
  REAL(8), DEVICE :: sum_sq_f
  REAL(8), DEVICE, DIMENSION (:),  ALLOCATABLE :: Xs1, X1, Xdif1
  REAL(8), DEVICE, DIMENSION (:,:,:),  ALLOCATABLE :: PSI_IN, S1, S2, S3, &
                                                      TEMPS, AIM, AIP, AJM, AJP, AKM, AKP, &
                                                      AP, CON
END MODULE VAR_PSI_IN_d
