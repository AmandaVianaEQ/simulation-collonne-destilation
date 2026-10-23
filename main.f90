PROGRAM colonne_dest
    IMPLICIT NONE
    INTEGER i,j,k,p,Nconst,NET,N_ALIM,Nvariables,npairs
    DOUBLE PRECISION Pamb, Qboilleur, Qcondenseur, Alim, sumKx, h_alim, psat_val, T_alim, delta_T,alpha,dx,norme
    DOUBLE PRECISION, ALLOCATABLE :: z_alim(:), x_coll(:,:), y_coll(:,:),omega(:)
    DOUBLE PRECISION, ALLOCATABLE :: A_ant(:), B_ant(:), C_ant(:), D_ant(:), E_ant(:)
    DOUBLE PRECISION, ALLOCATABLE :: CPvap(:), CPliq(:), Delta_Hvap(:), Teb(:)
    DOUBLE PRECISION, ALLOCATABLE :: COEFNRTL(:,:), GAMMA(:)
    DOUBLE PRECISION, ALLOCATABLE :: hliq_melange(:),hvap_melange(:), K_eq(:,:)
    DOUBLE PRECISION, PARAMETER :: Tref      = 150.0D0
    DOUBLE PRECISION, PARAMETER :: R         = 8314.0D0  !8,314 x 10^3
    DOUBLE PRECISION, PARAMETER :: Rnrtl     = 1.98751D0 ! R[cal/mol]

    DOUBLE PRECISION, ALLOCATABLE :: V(:),L(:),Temp(:),Rec_dest(:),Rec_resid(:)
    DOUBLE PRECISION, ALLOCATABLE :: F(:,:),Vecteur_Variables(:,:),F0(:,:),Fp(:,:), Xp(:,:), delta(:,:), T_temp(:)

    DOUBLE PRECISION :: Calc_Psat, Calc_Hliq, Calc_Hvap
    DOUBLE PRECISION :: D_dist, Lref, Vref, W_res, Reflux

!   PARAMETERS METHODE MRSL21
    INTEGER N,M,N1A,NLA,NC,IER,IRAZ,IMET
    DOUBLE PRECISION, ALLOCATABLE :: A(:,:,:), B(:,:,:), C(:,:,:)
    DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-12

!   INTERFACE
    INTERFACE
        SUBROUTINE lecture_donnees(Nconst,NET,N_ALIM,npairs,Pamb,Qboilleur,Qcondenseur,COEFNRTL,Alim,Teb, &
                                   z_alim, A_ant, B_ant, C_ant, D_ant, E_ant, CPvap, CPliq, Delta_Hvap)
            INTEGER Nconst, NET, N_ALIM, npairs
            DOUBLE PRECISION Pamb, Qboilleur, Qcondenseur, Alim
            DOUBLE PRECISION, ALLOCATABLE :: COEFNRTL(:,:), Teb(:), z_alim(:)
            DOUBLE PRECISION, ALLOCATABLE :: A_ant(:), B_ant(:), C_ant(:), D_ant(:), E_ant(:)
            DOUBLE PRECISION, ALLOCATABLE :: CPvap(:), CPliq(:), Delta_Hvap(:)

        END SUBROUTINE
    END INTERFACE

    CALL lecture_donnees(Nconst,NET,N_ALIM,npairs,Pamb,Qboilleur,Qcondenseur,COEFNRTL,Alim,Teb, &
                        z_alim, A_ant, B_ant, C_ant, D_ant, E_ant,CPvap, CPliq, Delta_Hvap)

    Nvariables = 2*Nconst+3 !Nombre de variables

    ALLOCATE(V(NET),L(NET),Temp(NET),K_eq(Nconst,NET),GAMMA(Nconst),Rec_dest(Nconst),Rec_resid(Nconst))
    ALLOCATE(hliq_melange(NET),hvap_melange(NET),T_temp(NET))
    ALLOCATE(x_coll(Nconst,NET), y_coll(Nconst,NET),omega(NET))

!   INITIALISATION
!    CALL initialisation(Nconst,NET,V,L,Temp,x_coll,y_coll)

    M = NET
    N = Nvariables
    NLA = Nconst+2
    N1A = Nconst+2
    NC  = Nconst+2
    IRAZ = 1
    IMET = 1
    alpha = 1.0D0

    ALLOCATE(A(NLA,N1A:N,M), B(N, N, M), C(N, NC, M))
    ALLOCATE(Xp(N,M))
    ALLOCATE(F(N,M),F0(N,M),Fp(N,M),delta(N,M),Vecteur_Variables(Nvariables,NET))

!   REMPLISSAGE DU VECTEUR DE VARIABLES
!    DO j = 1, NET
!        !Vapeur
!        Vecteur_Variables(1, j) = V(j)
!        DO i = 1, Nconst
!            !Composition de vapeur
!            Vecteur_Variables(1 + i, j) = y_coll(i,j)
!            !Composition de Liquide
!            Vecteur_Variables(Nconst + 2 + i, j) = x_coll(i,j)
!        END DO
!        !Temperature
!        Vecteur_Variables(Nconst + 2, j) = Temp(j)
!        !Liquide
!        Vecteur_Variables(Nvariables, j) = L(j)
!    END DO

!   CALCUL TEMPERATURE DE BULLE         !!!!!!!!!!!!!!!
    T_alim = SUM(z_alim*Teb)

    DO i = 1, 100
        CALL NRTL(Nconst, T_alim, z_alim, COEFNRTL, Rnrtl, GAMMA)
        sumKx = 0.0D0
        DO j = 1, Nconst
            Psat_val = Calc_Psat(T_alim, A_ant(j), B_ant(j), C_ant(j), D_ant(j), E_ant(j))
            sumKx = sumKx + GAMMA(j)*z_alim(j)*Psat_val/Pamb
        END DO

        IF (ABS(sumKx - 1.0D0) < 1.0D-8) EXIT

        delta_T = -30.0D0*LOG(sumKx)
        IF (delta_T >  10.0D0) delta_T =  10.0D0
        IF (delta_T < -10.0D0) delta_T = -10.0D0

        T_alim = T_alim + delta_T

    END DO

    !   CALCUL ENTHALPIE ALIMENTATION
    h_alim = 0.0D0
    DO i = 1, Nconst
        h_alim = h_alim + z_alim(i)*Calc_Hliq(T_alim, CPliq(i))
    END DO

    WRITE(*,'(A,F8.3,A,I3,A)') "Temperature de bulle alim : ", T_alim, " K (en ", i, " iterations)"

    Temp(1) = Teb(1)
    Temp(NET) = Teb(2)

    DO j=1,NET
        Temp(j) = Temp(1) + (Temp(NET) - Temp(1))*(j-1)/(NET-1)
    END DO

    Vref   = Qcondenseur/(SUM(z_alim*Delta_Hvap))
    D_dist = 0.5D0*Alim
    W_res  = Alim - D_dist
    Lref   = Vref - D_dist

    DO j = 1, NET
        IF (j == 1) THEN
            V(j) = D_dist
        ELSE
            V(j) = Vref
        END IF

        IF (j < N_ALIM) THEN
            L(j) = Lref
        ELSE IF (j < NET) THEN
            L(j) = Lref + Alim
        ELSE
            L(j) = W_res
        END IF
    END DO

!   REMPLISSAGE DU VECTEUR DE VARIABLES
    DO j = 1, NET
        Vecteur_Variables(1, j) = V(j)                          !Vapeur
        DO i = 1, Nconst
            Vecteur_Variables(1 + i, j) = z_alim(i)             !x_coll
            Vecteur_Variables(Nconst + 2 + i, j) = z_alim(i)    !y_coll
        END DO
        Vecteur_Variables(Nconst + 2, j) = Temp(j)              !Temperature
        Vecteur_Variables(Nvariables, j) = L(j)                 !Liquide
    END DO

!   RESOLUTION DES EQUATIONS
    DO i = 1,1000
        CALL residus(Vecteur_Variables, F0, Nconst, NET, Nvariables, N_ALIM, Alim, z_alim, h_alim,Pamb, Qcondenseur, Qboilleur, &
                    A_ant, B_ant, C_ant, D_ant, E_ant, CPvap, CPliq, Delta_Hvap, Teb, COEFNRTL, Rnrtl)

        ! ARRET
        norme = 0.0D0

        DO j = 1, NET
            DO k = 1, N
                IF (k == Nconst+2) THEN
                    norme = MAX(norme, ABS(F0(k,j))/1.0D6)
                ELSE
                    norme = MAX(norme, ABS(F0(k,j)))
                END IF
            END DO
        END DO

        WRITE(*,'(A,I3,A,ES12.4)') 'Iteration ', i, '   ||F|| = ', norme

        IF (norme < EPS) EXIT

!       REMPLISSAGE JACOBIANA
        DO p = 1, M
            DO k = 1, N
                ! PERTURBATION
                dx = 1.0D-6*MAX(ABS(Vecteur_Variables(k,p)), 1.0D0)
                Xp = Vecteur_Variables
                Xp(k,p) = Vecteur_Variables(k,p) + dx

                ! RESIDUS PERTURBE
                CALL residus(Xp, Fp, Nconst, NET, Nvariables, N_ALIM, Alim, z_alim, h_alim,Pamb, Qcondenseur, &
                            Qboilleur,A_ant, B_ant, C_ant, D_ant, E_ant, CPvap, CPliq, Delta_Hvap, Teb, COEFNRTL, Rnrtl)

                B(:, k, p) = (Fp(:, p) - F0(:, p))/dx

                IF (p < M .AND. k >= N1A) THEN
                    A(:, k, p+1) = (Fp(1:NLA, p+1) - F0(1:NLA, p+1))/dx
                END IF

                IF (p > 1 .AND. k <= NC) THEN
                    C(:, k, p-1) = (Fp(:, p-1) - F0(:, p-1))/dx
                END IF
            END DO
        END DO

        delta = -F0

        CALL MRSL21(A, B, C, delta, N, M, N1A, NLA, NC, IER, IRAZ, IMET)
        Vecteur_Variables = Vecteur_Variables + alpha*delta
    END DO

    Reflux = L(1)/V(1)
    PRINT *, Reflux

    DO j=1,Nconst
        Rec_dest(j) = (Vecteur_Variables(1,1)*Vecteur_Variables(1+j,1))/(Alim*z_alim(j))
        Rec_resid = 1.0 - Rec_dest
    END DO

    Print *,y_coll(1,2)

    Print *, Rec_dest

    OPEN(50, FILE="results.csv", STATUS="REPLACE")
        WRITE(50,'(A)') 'Étage;T[K];x1;x2;x3;y1;y2;y3'
        DO j = 1, NET
            WRITE(50,'(I0,";",F8.3,*(";",F6.5))') j, Vecteur_Variables(Nconst+2, j), Vecteur_Variables(Nconst+3:2*Nconst+2, j), &
                                          Vecteur_Variables(2:Nconst+1, j)
        END DO
    CLOSE(50)



END PROGRAM

! ------------------------------------------------------------
!   FUNCTION CALCUL DE PSAT (ANTOINE)
    DOUBLE PRECISION FUNCTION Calc_Psat(T, A, B, C, D, E)
        IMPLICIT NONE
        DOUBLE PRECISION :: T, A, B, C, D, E
        Calc_Psat = EXP(A + B/T + C*LOG(T) + D*(T**E))
        RETURN
    END FUNCTION Calc_Psat

! ------------------------------------------------------------
!   FUNCTION CALCUL Hliq
    DOUBLE PRECISION FUNCTION Calc_Hliq(T, CPliq)
        DOUBLE PRECISION :: T, CPliq
        DOUBLE PRECISION, PARAMETER :: Tref  = 150.0D0
        Calc_Hliq = CPliq*(T-Tref)
    END FUNCTION Calc_Hliq

! ------------------------------------------------------------
!   FUNCTION CALCUL Hvap
    DOUBLE PRECISION FUNCTION Calc_Hvap(T, CPvap, CPliq, DHvap, Teb)
        DOUBLE PRECISION :: T, CPvap, CPliq, DHvap, Teb
        DOUBLE PRECISION, PARAMETER :: Tref  = 150.0D0
        Calc_Hvap = CPliq*(Teb - Tref) + DHvap + CPvap*(T - Teb)
    END FUNCTION Calc_Hvap
