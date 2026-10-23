SUBROUTINE residus(X, F, Nconst, NET, Nvariables, N_ALIM, Alim, z_alim, h_alim,Pamb, Qcondenseur, Qboilleur, A_ant, B_ant, &
                    C_ant, D_ant, E_ant, CPvap, CPliq, Delta_Hvap, Teb, COEFNRTL, Rnrtl)
    IMPLICIT NONE
    INTEGER :: Nconst, NET, Nvariables, N_ALIM
    DOUBLE PRECISION :: X(Nvariables, NET)
    DOUBLE PRECISION :: F(Nvariables, NET)
    DOUBLE PRECISION :: Alim, h_alim, Pamb, Qcondenseur, Qboilleur, Rnrtl
    DOUBLE PRECISION :: z_alim(Nconst), A_ant(Nconst), B_ant(Nconst), C_ant(Nconst)
    DOUBLE PRECISION :: D_ant(Nconst), E_ant(Nconst), CPvap(Nconst), CPliq(Nconst)
    DOUBLE PRECISION :: Delta_Hvap(Nconst), Teb(Nconst)
    DOUBLE PRECISION :: COEFNRTL(3, Nconst*(Nconst-1)/2)

    INTEGER :: i, j
    DOUBLE PRECISION :: V(NET), L(NET), Temp(NET), x_coll(Nconst,NET), y_coll(Nconst,NET)
    DOUBLE PRECISION :: K_eq(Nconst,NET), Hliq_melange(NET), Hvap_melange(NET), GAMMA(Nconst)
    DOUBLE PRECISION :: Psat_val, sumx, sumy
    DOUBLE PRECISION :: Calc_Psat, Calc_Hliq, Calc_Hvap

    DO j = 1, NET
        V(j) = X(1, j)
        DO i = 1, Nconst
            y_coll(i,j) = X(1+i, j)
            x_coll(i,j) = X(Nconst+2+i, j)
        END DO
        Temp(j) = X(Nconst+2, j)
        L(j)    = X(Nvariables, j)
    END DO

    DO j = 1, NET
        CALL NRTL(Nconst, Temp(j), x_coll(:,j), COEFNRTL, Rnrtl, GAMMA)
        Hliq_melange(j) = 0.0D0
        Hvap_melange(j) = 0.0D0
        DO i = 1, Nconst
            Psat_val  = Calc_Psat(Temp(j), A_ant(i), B_ant(i), C_ant(i), D_ant(i), E_ant(i))
            K_eq(i,j) = GAMMA(i)*Psat_val/Pamb
            Hliq_melange(j) = Hliq_melange(j) + x_coll(i,j)*Calc_Hliq(Temp(j),CPliq(i))
            Hvap_melange(j) = Hvap_melange(j) + y_coll(i,j)*Calc_Hvap(Temp(j),CPvap(i),CPliq(i),Delta_Hvap(i),Teb(i))
        END DO
    END DO

    DO j = 1, NET
        ! Bilan matiere total
        F(1,j) = -L(j) - V(j)
        IF (j > 1)       F(1,j) = F(1,j) + L(j-1)
        IF (j < NET)     F(1,j) = F(1,j) + V(j+1)
        IF (j == N_ALIM) F(1,j) = F(1,j) + Alim

        ! Bilans partiels
        DO i = 1, Nconst
            F(1+i,j) = -L(j)*x_coll(i,j) - V(j)*y_coll(i,j)
            IF (j > 1)       F(1+i,j) = F(1+i,j) + L(j-1)*x_coll(i,j-1)
            IF (j < NET)     F(1+i,j) = F(1+i,j) + V(j+1)*y_coll(i,j+1)
            IF (j == N_ALIM) F(1+i,j) = F(1+i,j) + Alim*z_alim(i)
        END DO

        ! Bilan energie
        F(Nconst+2,j) = -L(j)*Hliq_melange(j) - V(j)*Hvap_melange(j)
        IF (j > 1)       F(Nconst+2,j) = F(Nconst+2,j) + L(j-1)*Hliq_melange(j-1)
        IF (j < NET)     F(Nconst+2,j) = F(Nconst+2,j) + V(j+1)*Hvap_melange(j+1)
        IF (j == N_ALIM) F(Nconst+2,j) = F(Nconst+2,j) + Alim*H_Alim
        IF (j == 1)      F(Nconst+2,j) = F(Nconst+2,j) - Qcondenseur
        IF (j == NET)    F(Nconst+2,j) = F(Nconst+2,j) + Qboilleur

        ! Equilibres
        DO i = 1, Nconst
            F(Nconst+2+i,j) = y_coll(i,j) - K_eq(i,j)*x_coll(i,j)
        END DO

        ! Sommation
        sumx = 0.0D0
        sumy = 0.0D0

        DO i = 1, Nconst
            sumx = sumx + x_coll(i,j)
            sumy = sumy + y_coll(i,j)
        END DO

        F(Nvariables,j) = sumx - sumy

    END DO


END SUBROUTINE
