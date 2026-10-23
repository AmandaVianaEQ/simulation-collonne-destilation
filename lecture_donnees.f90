SUBROUTINE lecture_donnees(Nconst,NET,N_ALIM,npairs,Pamb,Qboilleur,Qcondenseur,COEFNRTL,Alim,Teb, &
                              z_alim, A_ant, B_ant, C_ant, D_ant, E_ant,CPvap, CPliq, Delta_Hvap)
    IMPLICIT NONE
    INTEGER i,Nconst,NET,N_ALIM,npairs
    DOUBLE PRECISION Pamb,Qboilleur,Qcondenseur,Alim
    DOUBLE PRECISION, ALLOCATABLE :: z_alim(:)
    DOUBLE PRECISION, ALLOCATABLE :: A_ant(:), B_ant(:), C_ant(:), D_ant(:), E_ant(:)
    DOUBLE PRECISION, ALLOCATABLE :: CPvap(:), CPliq(:), Delta_Hvap(:), Teb(:)
    DOUBLE PRECISION, ALLOCATABLE :: COEFNRTL(:,:)

! ----- Lecture de donnees.txt -----
    OPEN(10,FILE="donnees.txt",STATUS="OLD")
        READ(10,*) Pamb
        READ(10,*) Qboilleur
        READ(10,*) Qcondenseur
        READ(10,*) Nconst
        READ(10,*) NET
        READ(10,*) Alim
        READ(10,*) N_ALIM

        npairs = Nconst*(Nconst-1)/2

        ALLOCATE(A_ant(Nconst), B_ant(Nconst), C_ant(Nconst), D_ant(Nconst), E_ant(Nconst))
        ALLOCATE(z_alim(Nconst), CPvap(Nconst), CPliq(Nconst), Delta_Hvap(Nconst), Teb(Nconst))
        ALLOCATE(COEFNRTL(3, npairs))

        DO i = 1, Nconst
            READ(10,*) z_alim(i), A_ant(i), B_ant(i), C_ant(i), D_ant(i), E_ant(i)
        END DO

        DO i = 1, Nconst
            READ(10,*) CPvap(i), CPliq(i), Delta_Hvap(i), Teb(i)
        END DO
    CLOSE(10)

! ----- Lecture de thermodat.txt -----
    OPEN(20, FILE="thermodat.txt", STATUS="OLD")
        DO i = 1, npairs
            READ(20,*) COEFNRTL(1,i), COEFNRTL(2,i), COEFNRTL(3,i)
        END DO
    CLOSE(20)

END SUBROUTINE
