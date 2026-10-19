PROGRAM colonne_dest
    IMPLICIT NONE
    INTEGER i,j,Nconst,NET,N_ALIM,Nvariables,npairs,deplac
    DOUBLE PRECISION Pamb, Qboilleur, Qcondenseur, Alim
    DOUBLE PRECISION, ALLOCATABLE :: z_alim(:), x_coll(:,:), y_coll(:,:),omega(:)
    DOUBLE PRECISION, ALLOCATABLE :: A_ant(:), B_ant(:), C_ant(:), D_ant(:), E_ant(:)
    DOUBLE PRECISION, ALLOCATABLE :: CPvap(:), CPliq(:), Delta_Hvap(:), Teb(:)
    DOUBLE PRECISION, ALLOCATABLE :: COEFNRTL(:,:)
    DOUBLE PRECISION, PARAMETER :: Tref  = 150.0D0
    DOUBLE PRECISION, PARAMETER :: R     = 8314.0D0  !8,314 x 10^3
    DOUBLE PRECISION, PARAMETER :: Rnrtl = 1.98751D0 ! R[cal/mol]

    DOUBLE PRECISION, ALLOCATABLE :: V(:),L(:), Temp(:)
    DOUBLE PRECISION, ALLOCATABLE :: Vecteur_Eq(:),Vecteur_Variables(:)

!   PARAMETERS METHODE MRSL21
    INTEGER N,M,N1A,NLA,NC,IER,IRAZ,IMET
    DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-12
    DOUBLE PRECISION A(:),B(:),C(:),F(:)

    CALL lecture_donnees(Nconst,NET,N_ALIM,npairs,Pamb,Qboilleur,Qcondenseur,COEFNRTL &
                        x0, A_ant, B_ant, C_ant, D_ant, E_ant,CPvap, CPliq, Delta_Hvap)

    ALLOCATE(V(NET),L(NET),Temp(NET))
    ALLOCATE(z_alim(NC),x_coll(NET,NC), y_coll(NET,NC),omega(NET))
    ALLOCATE(A_ant(NC), B_ant(NC), C_ant(NC), D_ant(NC), E_ant(NC))
    ALLOCATE(CPvap(NC), CPliq(NC), Delta_Hvap(NC), Teb(NC))
    ALLOCATE(Vecteur_Eq(),Vecteur_Variables())

    Nvariables = 2*NC+3 !Nombre de variables

!   REMPLISSAGE DU VECTEUR DE VARIABLES
    DO j=1,NET
        deplac = (j-1)*Nvariables   !Déplacement a chaque loop
        !Vapeur
        Vecteur_Variables(deplac+1) = V(j)
        !Composition de vapeur
        DO i = 1, NC
            Vecteur_Variables(deplac + 1 + i) = y_coll(i,j)
        END DO
        !Temperature
        Vecteur_Variables(deplac + Nconst + 2) = Temp(j)
        !Composition de liquide
        DO i = 1, NC
            Vecteur_Variables(deplac + Nconst + 2 + i) = x_coll(i,j)
        END DO
        !Liquide
        Vecteur_Variables(desloc + Nvariables) = L(j)
    END DO

!   REMPLISSAGE DU VECTEUR D'EQUATIONS
    DO j=1,NET

    END DO



! ------------------------------------------------------------

!   PREENCHENDO O VETOR DE EQUAÇÕES
    DO i=1,NET

    END DO

! ------------------------------------------------------------


!   DONNEES INITIALISATION
    Temp(1) = 326.5
    Temp(N_ALIM) = 327.0






!    CALL MRSL21(A,B,C,F,N,M,N1A,NLA,NC,IER,IRAZ,IMET)


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
