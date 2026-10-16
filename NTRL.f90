SUBROUTINE NTRL(NC, T, xi, CNRTL, Rnrtl, GAMMA)
    IMPLICIT NONE
        INTEGER :: NC
        DOUBLE PRECISION :: T, Rnrtl
        DOUBLE PRECISION :: xi(NC), CNRTL(3, NC*(NC-1)/2)
        DOUBLE PRECISION :: GAMMA(NC)

        INTEGER :: i, j, k, l
        DOUBLE PRECISION :: TAU(NC,NC), G(NC,NC)
        DOUBLE PRECISION :: prodGX(NC), prodTAUGX(NC), LNGAMMA(NC)

        TAU = 0.0D0
        G   = 1.0D0

        k = 0
        DO i = 1, NC-1
            DO j = i+1, NC
                k = k+1
                TAU(i,j) = CNRTL(1,k) / (Rnrtl*T)
                TAU(j,i) = CNRTL(2,k) / (Rnrtl*T)
                G(i,j)   = EXP(-CNRTL(3,k)*TAU(i,j))
                G(j,i)   = EXP(-CNRTL(3,k)*TAU(j,i))
            END DO
        END DO

        DO i = 1, NC
            TAU(i,i) = 0.0D0
            G(i,i)   = 1.0D0
        END DO

        DO i = 1, NC
            prodGX(i)   = 0.0D0
            prodTAUGX(i) = 0.0D0
            DO l = 1, NC
                prodGX(i)   = prodGX(i)   + G(l,i)*xi(l)
                prodTAUGX(i) = prodTAUGX(i) + TAU(l,i)*G(l,i)*xi(l)
            END DO
        END DO

        DO i = 1, NC
            LNGAMMA(i) = 0.0D0
            IF (prodGX(i) > 0.0D0) LNGAMMA(i) = prodTAUGX(i)/prodGX(i)
            DO j = 1, NC
                LNGAMMA(i) = LNGAMMA(i) + (xi(j)*G(i,j)/prodGX(j))*(TAU(i,j) - prodTAUGX(j)/prodGX(j))
            END DO
            GAMMA(i) = EXP(LNGAMMA(i))
        END DO

    END SUBROUTINE NRTL
