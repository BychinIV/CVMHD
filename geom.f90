MODULE GEOM
IMPLICIT NONE
  CONTAINS

  !------------------------------------------------------  
  
  FUNCTION AP01(I, J, K)
    USE VAR
    IMPLICIT NONE
    INTEGER I, J, K
    REAL(8) AP01, H1

    !-----------------Cart coords----------------------------------------
    !H1 = 1.0_8
    !-----------------End Sph coords-------------------------------------

    !-----------------Sph coords-----------------------------------------
    H1 = 1.0_8
    !-----------------End Sph coords-------------------------------------

    AP01 = S1(I,J,K) / H1  
  END FUNCTION AP01
  
  !------------------------------------------------------ 
  
  FUNCTION AP02(I, J, K)
    USE VAR
    IMPLICIT NONE
    INTEGER I, J, K
    REAL(8) AP02, H2

    !-----------------Sph coords-----------------------------------------
    !H2 = 1.0_8
    !-----------------End Sph coords-------------------------------------

    !-----------------Sph coords-----------------------------------------
     H2 = X1(I)
    !-----------------End Sph coords-------------------------------------

    AP02 = S2(I,J,K) / H2  
  END FUNCTION AP02

  !------------------------------------------------------

  FUNCTION AP03(I, J, K)
    USE VAR
    IMPLICIT NONE
    INTEGER I, J, K
    REAL(8) AP03, H3

    !-----------------Sph coords-----------------------------------------
    !H3 = 1.0_8
    !-----------------End Sph coords-------------------------------------

    !-----------------Sph coords-----------------------------------------
     H3 = X1(I) * DSIN(X2(J))
    !-----------------End Sph coords-------------------------------------

    AP03 = S3(I,J,K) / H3  
  END FUNCTION AP03

  !------------------------------------------------------

  SUBROUTINE COMPUTE_DIFF_FLOW_SPH
    USE VAR
    IMPLICIT NONE
    INTEGER I, J, K
    REAL(8) H1, H2, H3, dl1, dl2, dl3, dl_hat1, dl_hat2, dl_hat3  


    D123 = 0.0_8;  D132 = 0.0_8; D213 = 0.0_8; D231 = 0.0_8; D312 = 0.0_8; D321 = 0.0_8
    F_123 = 0.0_8;  F_132 = 0.0_8; F_213 = 0.0_8; F_231 = 0.0_8; F_312 = 0.0_8; F_321 = 0.0_8 

    !--------------------------------------------D312, D321, F_312, F_321-----------------------------------------
    DO K = 2, N3-1
      DO J = 3, N2     !For full shell : DO J = 3, N2-1
        DO I = 3, N1   !Для шара: DO I = 3, N1
          H1 = 1.0_8
          H2 = Xs1(I)

          dl3 = Xs1(I) * DSIN(Xs2(J)) * (Xs3(K+1) - Xs3(K))
          dl_hat1 = xdif1(I)
          dl_hat2 = Xs1(I) * xdif2(J)
 
          F_312(I,J,K) = dl3 / H2
          F_321(I,J,K) = dl3 / H1

          J321(I,J,K) = viscm * dl3 * Xdif1(I) / &                    
                       (0.5_8 * (X1(I)*X1(I)-X1(I-1)*X1(I-1))*Xdif2(J))

          J312(I,J,K) = viscm * dl3 / &  
                       (0.5_8 * (X1(I)*X1(I)-X1(I-1)*X1(I-1)))
 
        ENDDO
      ENDDO
    ENDDO


    DO K = 2, N3-1
      DO J = 3, N2     !For full shell : DO J = 3, N2-1
        dl3 = Xs1(N1) * DSIN(Xs2(J)) * (Xs3(K+1) - Xs3(K))
        dl_hat2 = Xs1(N1) * xdif2(J)

        J312(N1,J,K) = viscm * dl3 / ( 2.0_8 * X1(N1) * xdif1(N1) )
        J321(N1,J,K) = viscm * dl3 / dl_hat2
      ENDDO
    ENDDO


    !--------------------------------------------D213, D231, F_213, F_231------------------------------
    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 3, N1  !Для шара: DO I = 3, N1
          H1 = 1.0_8
          H3 = Xs1(I) * DSIN(X2(J))

          dl2 = Xs1(I) * (Xs2(J+1) -Xs2(J))

          dl_hat1 = xdif1(I)
          dl_hat3 = Xs1(I) * DSIN(X2(J)) * xdif3(K)

        
          F_213(I,J,K) = dl2 / H3 
          F_231(I,J,K) = dl2 / H1

          J231(I,J,K) = viscm * dl2 * Xdif1(I) / &
                       (0.5_8 * (X1(I)*X1(I)-X1(I-1)*X1(I-1)) * sinX2(J)*Xdif3(K))
          J213(I,J,K) = viscm * dl2 / &
                       (0.5_8 * (X1(I)*X1(I)-X1(I-1)*X1(I-1)) * sinX2(J))


        ENDDO
      ENDDO
    ENDDO


    DO K = 2, N3-1
      DO J = 2, N2-1
          H3 = Xs1(N1) * DSIN(X2(J))
          dl2 = Xs1(N1) * (Xs2(J+1) - Xs2(J))
          dl_hat1 = xdif1(N1)
          dl_hat3 = Xs1(N1) * DSIN(X2(J)) * xdif3(K)

          J231(N1,J,K) = viscm * dl2 / dl_hat3
          J213(N1,J,K) = viscm * dl2 / (2.0_8 * H3 * dl_hat1)
      ENDDO
    ENDDO


    DO J = 2, N2-1
      DO I = 3, N1  !Для шара: DO I = 3, N1
          J213(I,J,N3) = J213(I,J,2) 
          J231(I,J,N3) = J231(I,J,2)

        
          F_213(I,J,N3) = F_213(I,J,2)
          F_231(I,J,N3) = F_231(I,J,2)
      ENDDO
    ENDDO

   !--------------------------------------------D123, D132, F_123, F_132-----------------------------

    DO K = 2, N3-1
      DO J = 3, N2  !For full shell : DO J = 3, N2-1
        DO I = 2, N1-1
          H2 = X1(I)
          H3 = X1(I) * DSIN(Xs2(J))

          dl1 = Xs1(I+1) - Xs1(I)


          F_123(I,J,K) = dl1 / H3 
          F_132(I,J,K) = dl1 / H2


          J123(I,J,K) = viscm * dl1 / ( X1(I)**2 * ( cosX2(J-1) - cosX2(J)) )
          J132(I,J,K) = viscm * dl1 * Xdif2(J) / ( X1(I)**2 * Xdif3(K) * ( cosX2(J-1) - cosX2(J)) )
        ENDDO
      ENDDO
    ENDDO	


    DO J = 3, N2  !For full shell : DO J = 3, N2-1
      DO I = 2, N1-1
        J123(I,J,N3) = J123(I,J,2)
        J132(I,J,N3) = J132(I,J,2)


        F_123(I,J,N3) = F_123(I,J,2) 
        F_132(I,J,N3) = F_132(I,J,2) 
      ENDDO
    ENDDO
    

  END SUBROUTINE COMPUTE_DIFF_FLOW_SPH


  !------------------------------------------------------
  
  SUBROUTINE UGRID(N1, XL1, X1_0, Xs1, X1, Xdif1)
    IMPLICIT NONE
    INTEGER N1, I
    REAL(8) XL1, X1_0,DX1,Xs1(N1),X1(N1),Xdif1(N1)

    DX1=XL1/DBLE(N1-2)
    DO I=2,N1
     Xs1(I)=X1_0+(I-2)*DX1
    END DO

    X1(1) = Xs1(2)
    DO I = 2, N1-1
     X1(I) = 0.5_8*(Xs1(I+1)+Xs1(I))
    ENDDO
    X1(N1) = Xs1(N1)

    DO I = 2, N1
     Xdif1(I) = X1(I)-X1(I-1)
    ENDDO
  END SUBROUTINE UGRID

  !------------------------------------------------------

  SUBROUTINE NUGRID_RHO_PSI_IN
    USE VAR, ONLY: r_i, Xdif1_shell => Xdif1
    USE VAR_PSI_IN
    IMPLICIT NONE
    INTEGER :: I, LSI
    REAL(8) :: xalfa, DX

    !-----------NUGRID BY RHO-------------
    DX = 2.0_8 * Xdif1_shell(2) 
    LSI = 0
 
    DO I = 1, N1 - 2
      LSI = LSI + I
    ENDDO

    xalfa = (r_i - DBLE(N1 - 1) * DX) / DBLE(LSI)
   
    Xs1(2) = 0.0_8
    Xs1(N1) = r_i
    DO I = N1-1, 3, -1
      Xs1(I) = Xs1(I+1) - DX - DBLE((N1-1) - I) * xalfa
    END DO

    X1(1) = Xs1(2)
    DO I = 2, N1-1
      X1(I) = 0.5_8 * (Xs1(I+1) + Xs1(I))
    ENDDO
    X1(N1) = Xs1(N1)

    DO I = 2, N1
     Xdif1(I) = X1(I)-X1(I-1)
    ENDDO

  END SUBROUTINE NUGRID_RHO_PSI_IN
    
  !------------------------------------------------------
  
  SUBROUTINE NUGRID_RHO_PSI_OUT
    USE VAR, ONLY: r_o, Xdif1_shell => Xdif1, N1_shell => N1 
    USE VAR_PSI_OUT
    IMPLICIT NONE
    INTEGER :: I, LSI
    REAL(8) :: xalfa, DX

    !-----------NUGRID BY RHO-------------
    DX = 2.0_8 * r_o * Xdif1_shell(N1_shell) / (r_o + Xdif1_shell(N1_shell))  
    LSI = 0
 
    DO I = 1, N1 - 2
      LSI = LSI + I
    ENDDO

    xalfa = (r_o - DBLE(N1 - 1) * DX) / DBLE(LSI)
   
    Xs1(2) = 0.0_8
    Xs1(N1) = r_o
    DO I = N1-1, 3, -1
      Xs1(I) = Xs1(I+1) - DX - DBLE((N1-1) - I) * xalfa
    END DO

    X1(1) = Xs1(2)
    DO I = 2, N1-1
      X1(I) = 0.5_8 * (Xs1(I+1) + Xs1(I))
    ENDDO
    X1(N1) = Xs1(N1)

    DO I = 2, N1
     Xdif1(I) = X1(I)-X1(I-1)
    ENDDO

  END SUBROUTINE NUGRID_RHO_PSI_OUT


  
  !------------------------------------------------------

  SUBROUTINE CHEB_GRID_SHELL
    USE VAR
    IMPLICIT NONE
    INTEGER :: I, J

    J = N1
    DO I = 2, N1cv_shell + 2      
      Xs1(J) = 0.5_8 * (r_i + r_o) + 0.5_8 * (r_o - r_i) * DCOS( (2.0_8 * (I - 1) - 1.0_8) * PI / (2.0_8 * (N1cv_shell + 1)))
      J = J - 1
    END DO

    Xs1(N1_core) = r_i
    Xs1(N1) = r_o

  END SUBROUTINE CHEB_GRID_SHELL

  !------------------------------------------------------

  SUBROUTINE NUGRID_CORE
    USE VAR
    IMPLICIT NONE
    INTEGER :: I, LSI
    REAL(8) :: xalfa, DX

    DX = Xs1(N1_core + 1) - Xs1(N1_core) 
    LSI = 0
 
    DO I = 1, N1_core - 2
      LSI = LSI + I
    ENDDO

    xalfa = (r_i - DBLE(N1_core - 1) * DX) / DBLE(LSI)
   
    Xs1(2) = 0.0_8

    DO I = N1_core-1, 3, -1
      Xs1(I) = Xs1(I+1) - DX - DBLE((N1_core-1) - I) * xalfa
    END DO

    X1(1) = Xs1(2)
    DO I = 2, N1-1
      X1(I) = 0.5_8 * (Xs1(I+1) + Xs1(I))
    ENDDO
    X1(N1) = Xs1(N1)

    DO I = 2, N1
      Xdif1(I) = X1(I)-X1(I-1)
    ENDDO

  END SUBROUTINE NUGRID_CORE

  !------------------------------------------------------

  SUBROUTINE NUGRID2_CORE
    USE VAR
    IMPLICIT NONE
    INTEGER :: I, LSI
    REAL(8) :: xalfa, DX, r_eps

    r_eps = 0.05_8
    DX = Xs1(N1_core + 1) - Xs1(N1_core) 
    LSI = 0
 
    DO I = 1, N1_core - 3
      LSI = LSI + I
    ENDDO

    xalfa = ((r_i - r_eps) - DBLE(N1_core - 2) * DX) / DBLE(LSI)
   
    Xs1(2) = 0.0_8
    Xs1(3) = r_eps

    DO I = N1_core-1, 4, -1
      Xs1(I) = Xs1(I+1) - DX - DBLE((N1_core-1) - I) * xalfa
    END DO



    X1(1) = Xs1(2)
    DO I = 2, N1-1
      X1(I) = 0.5_8 * (Xs1(I+1) + Xs1(I))
    ENDDO
    X1(N1) = Xs1(N1)

    DO I = 2, N1
      Xdif1(I) = X1(I)-X1(I-1)
    ENDDO

  END SUBROUTINE NUGRID2_CORE


  !------------------------------------------------------
  
  SUBROUTINE NUGRID_RHO_MHD
    USE VAR
    IMPLICIT NONE
    INTEGER :: I, LK, LSI
    REAL(8) :: xalfa, DX, XCV(N1)


    !-----------NUGRID BY RHO-------------
    DX = 0.005_8
    LK = N1 / 2
    LSI = 0
 
    DO I = 1, LK - 2
      LSI = LSI + I
    ENDDO

    xalfa = (0.5_8 * XL1 - DBLE(LK - 1) * DX) / DBLE(LSI)
      
    Xs1(2) = X1_0
    DO I = 3, LK + 1
      Xs1(I) = Xs1(I-1) + DX + DBLE(I-3) * xalfa
    END DO

    DO I = 2, LK
      XCV(I) = Xs1(I+1) - Xs1(I)
    ENDDO

    DO I = 1, LK-1
      XCV(LK+I) = XCV(LK+1-I)
    ENDDO

    DO I = LK+1, N1
      Xs1(I) = Xs1(I-1) + XCV(I-1)
    END DO

    X1(1) = Xs1(2)
    DO I = 2, N1-1
      X1(I) = 0.5_8 * (Xs1(I+1) + Xs1(I))
    ENDDO
    X1(N1) = Xs1(N1)

    DO I = 2, N1
     Xdif1(I) = X1(I)-X1(I-1)
    ENDDO


  END SUBROUTINE NUGRID_RHO_MHD
  
  !------------------------------------------------------
  
  SUBROUTINE BUILD_GRID
    USE VAR
    IMPLICIT NONE
    INTEGER :: I, LK, LSI
    REAL(8) :: xalfa, DX, XCV(N1)

    !CALL UGRID(N1, XL1, X1_0, Xs1, X1, Xdif1)
    !CALL NUGRID_RHO_MHD
    CALL CHEB_GRID_SHELL
    CALL NUGRID_CORE
    CALL UGRID(N2, XL2, X2_0, Xs2, X2, Xdif2)
    CALL UGRID(N3, XL3, X3_0, Xs3, X3, Xdif3)
    
    CALL NUGRID_RHO_PSI_OUT

    CALL COMPUTE_GEOM

  END SUBROUTINE BUILD_GRID
   
  !------------------------------------------------------
  
  SUBROUTINE COMPUTE_GEOM
    USE VAR
    USE VAR_PSI_OUT, ONLY: N1_PSI_OUT => N1, S1_PSI_OUT => S1, &
                           S2_PSI_OUT => S2, S3_PSI_OUT => S3, &
                           X1_PSI_OUT => X1, Xs1_PSI_OUT => Xs1 
    USE VAR_PSI_IN, ONLY:  N1_PSI_IN => N1, S1_PSI_IN => S1, &
                           S2_PSI_IN => S2, S3_PSI_IN => S3, &
                           X1_PSI_IN => X1, Xs1_PSI_IN => Xs1 
    IMPLICIT NONE
    INTEGER :: I, J, K

    !X3 - циклическая координата =>

    X_3(0) = - X3(3)
    X_3(1) = - X3(2) 
    X_3(2:N3-1) = X3(2:N3-1)
    X_3(N3) = X3(N3-1) + Xdif3(N3-1)
    X_3(N3+1) = X3(N3-1) + 2.0_8 * Xdif3(N3-1)


    Xdif3(2) = Xdif3(2) + Xdif3(N3)    
    Xdif3(N3) = Xdif3(2)

    DO I = 3, N1-1
     Kx1(I) = (X1(I) - Xs1(I)) / Xdif1(I)
    ENDDO
    Kx1(2) = 1.0_8
    Kx1(N1) = 0.0_8

    DO J = 3, N2-1
     Kx2(J) = (X2(J) - Xs2(J)) / Xdif2(J)
    ENDDO
    Kx2(2) = 1.0_8
    Kx2(N2) = 0.0_8

    DO K = 2, N3-1
     Kx3(K) = (X3(K) - Xs3(K)) / Xdif3(K)
    ENDDO
    Kx3(N3) = Kx3(2)

    DO I = 2, N1-1
      X1cv(I) = Xs1(I+1) - Xs1(I)
    ENDDO

    DO J = 2, N2-1
      X2cv(J) = Xs2(J+1) - Xs2(J)
    ENDDO

    DO K = 2, N3-1
      X3cv(K) = Xs3(K+1) - Xs3(K)
    ENDDO

    !X3 - циклическая координата => 
    X3cv(1) = X3cv(N3-1) 
    X3cv(N3) = X3cv(2)

    !---------Trigonometric---------
    DO J = 1, N2
      cosXs2(J) = DCOS(Xs2(J))
      cosX2(J)  = DCOS(X2(J))
      sinXs2(J) = DSIN(Xs2(J))
      sinX2(J)  = DSIN(X2(J))
    END DO

    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 2, N1
         S1(I,J,K) = Xs1(I)*Xs1(I)*(DCOS(Xs2(J)) - DCOS(Xs2(J+1)))*(Xs3(K+1) - Xs3(K))
        ENDDO
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO J = 2, N2
        DO I = 2, N1-1
          S2(I,J,K) = 0.5_8*DSIN(Xs2(J))*(Xs3(K+1) - Xs3(K))*(Xs1(I+1)**2-Xs1(I)**2)
        ENDDO
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO I = 2, N1-1
        S2(I,2,K) = 0.0_8
        !For full shell:
        !S2(I,N2,K) = 0.0_8
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 2, N1-1
          S3(I,J,K) = X1(I)*(Xs1(I+1) - Xs1(I))*(Xs2(J+1) - Xs2(J))
        ENDDO
      ENDDO
    ENDDO

    !X3 - циклическая координата => 
    DO J = 2, N2-1
      DO I = 2, N1-1
        S3(I,J,1) = S3(I,J,N3-1)
        S3(I,J,N3) = S3(I,J,2)
      ENDDO
    ENDDO

    !-------------------- PSI_OUT-GEOM------------------
    
    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 2, N1_PSI_OUT
         S1_PSI_OUT(I,J,K) = Xs1_PSI_OUT(I) * Xs1_PSI_OUT(I) * &
                             (DCOS(Xs2(J)) - DCOS(Xs2(J+1))) * (Xs3(K+1) - Xs3(K))
        ENDDO
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO J = 2, N2
        DO I = 2, N1_PSI_OUT-1
          S2_PSI_OUT(I,J,K) = 0.5_8 * DSIN(Xs2(J)) * (Xs3(K+1) - Xs3(K)) * &
                              (Xs1_PSI_OUT(I+1)**2 - Xs1_PSI_OUT(I)**2)
        ENDDO
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO I = 2, N1_PSI_OUT-1
        S2_PSI_OUT(I,2,K) = 0.0_8
        !For full shell:
        !S2_PSI_OUT(I,N2,K) = 0.0_8
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 2, N1_PSI_OUT-1
          S3_PSI_OUT(I,J,K) = X1_PSI_OUT(I) * (Xs1_PSI_OUT(I+1) - Xs1_PSI_OUT(I)) * &
                              (Xs2(J+1) - Xs2(J))
        ENDDO
      ENDDO
    ENDDO

    !X3 - циклическая координата => 
    DO J = 2, N2-1
      DO I = 2, N1_PSI_OUT-1
        S3_PSI_OUT(I,J,1) = S3_PSI_OUT(I,J,N3-1)
        S3_PSI_OUT(I,J,N3) = S3_PSI_OUT(I,J,2)
      ENDDO
    ENDDO
    
    !-------------------- PSI_IN-GEOM------------------
    
    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 2, N1_PSI_IN
         S1_PSI_IN(I,J,K) = Xs1_PSI_IN(I) * Xs1_PSI_IN(I) * &
                            (DCOS(Xs2(J)) - DCOS(Xs2(J+1))) * (Xs3(K+1) - Xs3(K))
        ENDDO
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO J = 2, N2
        DO I = 2, N1_PSI_IN-1
          S2_PSI_IN(I,J,K) = 0.5_8 * DSIN(Xs2(J)) * (Xs3(K+1) - Xs3(K)) * &
                              (Xs1_PSI_IN(I+1)**2 - Xs1_PSI_IN(I)**2)
        ENDDO
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO I = 2, N1_PSI_IN-1
        S2_PSI_IN(I,2,K) = 0.0_8
        !For full shell:
        !S2_PSI_IN(I,N2,K) = 0.0_8
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 2, N1_PSI_IN-1
          S3_PSI_IN(I,J,K) = X1_PSI_IN(I) * (Xs1_PSI_IN(I+1) - Xs1_PSI_IN(I)) * &
                             (Xs2(J+1) - Xs2(J))
        ENDDO
      ENDDO
    ENDDO

    !X3 - циклическая координата => 
    DO J = 2, N2-1
      DO I = 2, N1_PSI_IN-1
        S3_PSI_IN(I,J,1) = S3_PSI_IN(I,J,N3-1)
        S3_PSI_IN(I,J,N3) = S3_PSI_IN(I,J,2)
      ENDDO
    ENDDO
    
    !-------------------- U1-GEOM ----------------------

    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 1, N1
          S1u1(I,J,K) = X1(I) * X1(I) * (DCOS(Xs2(J)) - DCOS(Xs2(J+1))) * (Xs3(K+1) - Xs3(K)) 
        ENDDO
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO J = 2, N2
        DO I = 4, N1-2
          S2u1(I,J,K) = 0.5_8 * (X1(I) * X1(I) - X1(I-1) * X1(I-1)) * DSIN(Xs2(J)) * (Xs3(K+1) - Xs3(K)) 
          S2u1_1(I,J,K) = (X1(I)-X1(I-1)) * DSIN(Xs2(J)) * (Xs3(K+1) - Xs3(K))
        ENDDO
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO J = 2, N2
        S2u1(3,J,K) = 0.5_8 * (X1(3) * X1(3) - X1(1) * X1(1)) * DSIN(Xs2(J)) * (Xs3(K+1) - Xs3(K)) 
        S2u1_1(3,J,K) = (X1(3) - X1(1)) * DSIN(Xs2(J)) * (Xs3(K+1) - Xs3(K))

        S2u1(N1-1,J,K) = 0.5_8 * (X1(N1) * X1(N1) - X1(N1-2) * X1(N1-2)) * DSIN(Xs2(J)) * (Xs3(K+1) - Xs3(K)) 
        S2u1_1(N1-1,J,K) = (X1(N1) - X1(N1-2)) * DSIN(Xs2(J)) * (Xs3(K+1) - Xs3(K)) 
      ENDDO
    ENDDO

!     For full shell:
!     DO K = 2, N3-1
!       DO I = 3, N1-1
!         S2u1(I, N2, K) = 0.0_8
!         S2u1_1(I, N2, K) = 0.0_8
!       ENDDO
!     ENDDO

    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 4, N1-2
          S3u1(I,J,K) = 0.5_8 * (X1(I) * X1(I) - X1(I-1) * X1(I-1)) * (Xs2(J+1) - Xs2(J))
        ENDDO
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO J = 2, N2-1
        S3u1(3,J,K) = 0.5_8 * (X1(3) * X1(3) - X1(1) * X1(1)) * (Xs2(J+1) - Xs2(J))
        S3u1(N1-1,J,K) = 0.5_8 * (X1(N1) * X1(N1) - X1(N1-2) * X1(N1-2)) * (Xs2(J+1) - Xs2(J))
      ENDDO
    ENDDO

    DO J = 2, N2-1
      DO I = 2, N1-1
        S3u1(I,J,N3) = S3u1(I,J,2)
      ENDDO
    ENDDO

   !--------------------For SourceU1------------------------

    DO I = 4, N1-2
      VVX1S(I) = X1(I)**3 - Xs1(I)**3
      VVX2S(I) = Xs1(I)**3 - X1(I-1)**3
    ENDDO

    VVX1S(3) = X1(3)**3 - Xs1(3)**3
    VVX2S(3) = Xs1(3)**3 - X1(1)**3

    VVX1S(N1-1) = X1(N1)**3 - Xs1(N1-1)**3
    VVX2S(N1-1) = Xs1(N1-1)**3 - X1(N1-2)**3


    DO I = 4, N1-2
      XsDifP4(I) = Xs1(I)**4 - X1(I-1)**4
      XDifP4(I) = X1(I)**4 - Xs1(I)**4
    ENDDO

    XsDifP4(3) = Xs1(3)**4 - X1(1)**4                              
    XDifP4(3) = X1(3)**4 - Xs1(3)**4

    XsDifP4(N1-1) =  Xs1(N1-1)**4 - X1(N1-2)**4                             
    XDifP4(N1-1) = X1(N1)**4 - Xs1(N1-1)**4

  
   !--------------------For CofU1----------------------------

    DO J = 3, N2-1 !For Full Shell: DO J = 3, N2-2
      VVY1(J) = cosXs2(J) - cosX2(J)
      VVY2(J) = cosX2(J) - cosXs2(J+1)
    ENDDO
    VVY1(2) = 0.0_8
    VVY2(2) = 1.0_8 - cosXs2(3)
    !For Full Shell
    !VVY1(N2-1)  = cosXs2(N2-1) + 1.0_8
    !VVY2(N2-1)  = 0.0_8

    DO I = 4, N1-2
      VVX1(I) = X1(I)**2 - Xs1(I)**2
      VVX2(I) = Xs1(I)**2 - X1(I-1)**2
      VYX1(I) = X1(I) - Xs1(I)
      VYX2(I) = Xs1(I) - X1(I-1)
      APX(I)  = X1(I) - X1(I-1)
    ENDDO
    VVX1(3) = X1(3)**2 - Xs1(3)**2
    VVX2(3) = Xs1(3)**2 - X1(1)**2
    VYX1(3) = X1(3) - Xs1(3)
    VYX2(3) = Xs1(3) - X1(1)
    APX(3)  = X1(3) - X1(1)
    VVX1(N1-1) = X1(N1)**2 - Xs1(N1-1)**2
    VVX2(N1-1) = Xs1(N1-1)**2 - X1(N1-2)**2
    VYX1(N1-1) = X1(N1) - Xs1(N1-1)
    VYX2(N1-1) = Xs1(N1-1) - X1(N1-2) 
    APX(N1-1) = X1(N1) - X1(N1-2)
    

    !---------------- U2-GEOM -----------------

    DO K = 2, N3-1
      DO J = 4, N2-2
        DO I = 2, N1
          S1u2(I,J,K) = Xs1(I) * Xs1(I) * (DCOS(X2(J-1)) - DCOS(X2(J))) * (Xs3(K+1) - Xs3(K))
        ENDDO
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO I = 2, N1
        S1u2(I,3,K) = Xs1(I) * Xs1(I) * (DCOS(X2(1)) - DCOS(X2(3))) * (Xs3(K+1) - Xs3(K))
        S1u2(I,N2-1,K) = Xs1(I) * Xs1(I) * (DCOS(X2(N2-2)) - DCOS(X2(N2))) * (Xs3(K+1) - Xs3(K))
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO J = 1, N2
        DO I = 2, N1-1
          S2u2(I,J,K) = 0.5_8 * DSIN(X2(J)) * (Xs3(K+1) - Xs3(K)) * (Xs1(I+1)**2 - Xs1(I)**2)
        ENDDO
      ENDDO
    ENDDO

    !For full shell
!     DO K = 2, N3-1
!       DO I = 2, N1-1 
!         S2u2(I,N2,K) = 0.0_8
!       ENDDO
!     ENDDO

    DO K = 2, N3-1
      DO J = 4, N2-2
        DO I = 2, N1-1
          S3u2(I,J,K) = 0.5_8 * (Xs1(I+1)**2 - Xs1(I)**2) * Xdif2(J)
        ENDDO
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO I = 2, N1-1
        S3u2(I,3,K) = 0.5_8 * (Xs1(I+1)**2 - Xs1(I)**2) * (X2(3) - X2(1))
        S3u2(I,N2-1,K) = 0.5_8 * (Xs1(I+1)**2 - Xs1(I)**2) * (X2(N2) - X2(N2-2))
      ENDDO
    ENDDO

    DO J = 2, N2-1
      DO I = 2, N1-1
        S3u2(I,J,N3) = S3u2(I,J,2)
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 2, N1-1
          VCV(I,J,K) = (Xs1(I+1)**3-Xs1(I)**3)*(DCOS(Xs2(J))-DCOS(Xs2(J+1)))*(Xs3(K+1) - Xs3(K))/3.0_8
        ENDDO
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 4, N1-2
          VCV1(I,J,K) = (X1(I)**3-X1(I-1)**3)*(DCOS(Xs2(J))-DCOS(Xs2(J+1)))*(Xs3(K+1) - Xs3(K))/3.0_8
        ENDDO
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO J = 2, N2-1     
         VCV1(3,J,K) = (X1(3)**3-X1(1)**3)*(DCOS(Xs2(J))-DCOS(Xs2(J+1)))*(Xs3(K+1) - Xs3(K))/3.0_8

         VCV1(N1-1,J,K) = (X1(N1)**3-X1(N1-2)**3)*(DCOS(Xs2(J))-DCOS(Xs2(J+1)))*(Xs3(K+1) - Xs3(K))/3.0_8    
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO J = 4, N2-2
        DO I = 2, N1-1
           VCV2(I,J,K) = (Xs1(I+1)**3-Xs1(I)**3)*(DCOS(X2(J-1))-DCOS(X2(J)))*(Xs3(K+1) - Xs3(K))/3.0_8
        ENDDO
      ENDDO
    ENDDO

    DO K = 2, N3-1
      DO I = 2, N1-1
         VCV2(I,3,K) = (Xs1(I+1)**3-Xs1(I)**3)*(DCOS(X2(1))-DCOS(X2(3)))*(Xs3(K+1) - Xs3(K))/3.0_8
         VCV2(I,N2-1,K) = (Xs1(I+1)**3-Xs1(I)**3)*(DCOS(X2(N2-2))-DCOS(X2(N2)))*(Xs3(K+1) - Xs3(K))/3.0_8
      ENDDO
    ENDDO

    !----------------------For SourceU2------------

    DO J = 4, N2-2
      WWY1S(J) = sinX2(J)**2 - sinXs2(J)**2       
      WWY2S(J) = sinXs2(J)**2 - sinX2(J-1)**2
    ENDDO

    WWY1S(3) = sinX2(3)**2 - sinXs2(3)**2
    WWY2S(3) = sinXs2(3)**2 - sinX2(1)**2

    WWY1S(N2-1) = sinX2(N2)**2 - sinXs2(N2-1)**2
    WWY2S(N2-1) = sinXs2(N2-1)**2 - sinX2(N2-2)**2

    !----------------------For CofU2---------------

    DO J = 4, N2-2
      WWY1(J) = sinX2(J) - sinXs2(J)       
      WWY2(J) = sinXs2(J) - sinX2(J-1)

      UYY(J) = 1.0_8

      WZY1(J) = sinX2(J) - sinXs2(J)
      WZY2(J) = sinXs2(J) - sinX2(J-1)

      Spv1J(J) = Xdif2(J) / sinXs2(J)

      Spvy1(J) = cosXs2(J) - cosX2(J)
      Spvy2(J) = cosX2(J-1) - cosXs2(J) 
    ENDDO

    WWY2(3) = sinXs2(3) - sinX2(1)
    WWY1(N2-1) = sinX2(N2) - sinXs2(N2-1)

    UYY(3) = 1.0_8 + Xdif2(2) / Xdif2(3)
    UYY(N2-1) = 1.0_8 + Xdif2(N2) / Xdif2(N2-1)

    WZY2(3) = sinXs2(3) - sinX2(1)                   
    WZY1(N2-1) = sinX2(N2) - sinXs2(N2-1) 

    Spv1J(3) = (X2(3) - X2(1)) / sinXs2(3)
    Spv1J(N2-1) = (X2(N2) - X2(N2-2)) / sinXs2(N2-1) 

    Spvy1(3) = cosXs2(3) - cosX2(3)
    Spvy2(3) = cosX2(1) - cosXs2(3)

    Spvy1(N2-1) = cosXs2(N2-1) - cosX2(N2)
    Spvy2(N2-1) = cosX2(N2-2) - cosXs2(N2-1)


    !----------------- U3-GEOM --------------------

    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 2, N1
          S1u3(I,J,K) = Xs1(I) * Xs1(I) * (cosXs2(J) - cosXs2(J+1)) * Xdif3(K)
        ENDDO
      ENDDO
    ENDDO


    DO K = 2, N3-1
      DO J = 2, N2
        DO I = 2, N1-1
          S2u3(I,J,K) = 0.5_8 * sinXs2(J) * Xdif3(K) * (Xs1(I+1)**2 - Xs1(I)**2)
        ENDDO
      ENDDO
    ENDDO

    !For full shell
!     DO K = 2, N3-1
!       DO I = 2, N1-1
!         S2u3(I,N2,K) = 0.0_8
!       ENDDO
!     ENDDO

    DO K = 2, N3-1
      DO J = 2, N2-1
        DO I = 2, N1-1
         VCV3(I,J,K) = (Xs1(I+1)**3 - Xs1(I)**3) * (cosXs2(J) - cosXs2(J+1)) * Xdif3(K) / 3.0_8
       ENDDO
     ENDDO
    ENDDO
 
    !--------------For SourceU3-------------

    DO J = 3, N2-2
      VZY1S(J) = sinX2(J)**2 - sinXs2(J)**2
      VZY2S(J) = sinXs2(J+1)**2 - sinX2(J)**2
    ENDDO

    VZY1S(2) = 0.0_8
    VZY2S(2) = sinXs2(3)**2 - sinX2(1)**2

    VZY1S(N2-1) = sinX2(N2)**2 - sinXs2(N2-1)**2
    VZY2S(N2-1) = 0.0_8

    !--------------For CofU3----------------
    DO J = 3, N2-1 !For Full Shell: DO J = 3, N2-2
      VZY1(J) = sinX2(J) - sinXs2(J)
      VZY2(J) = sinXs2(J+1) - sinX2(J)

      Spwy1(J) = sinX2(J) - sinXs2(J)
      Spwy2(J) = sinXs2(J+1) - sinX2(J)
    ENDDO

    VZY1(2) = 0.0_8
    VZY2(2) = sinXs2(3) - sinX2(1)
    
    !For Full Shell:
    !VZY1(N2-1) = sinX2(N2) - sinXs2(N2-1)
    !VZY2(N2-1)  = 0.0_8

    Spwy1(2) = 0.0_8
    Spwy2(2) = sinXs2(3)

    !For Full Shell:
    !Spwy1(N2-1) = - sinXs2(N2-1)
    !Spwy2(N2-1) = 0.0_8

  END SUBROUTINE COMPUTE_GEOM 

END MODULE GEOM
