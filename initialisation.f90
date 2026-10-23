!SUBROUTINE initialisation(Nconst,NET,V,L,Temp,x_coll,y_coll)
!    IMPLICIT NONE
!    INTEGER :: Nconst, NET
!    DOUBLE PRECISION :: V(NET), L(NET), Temp(NET), Teb()
!    DOUBLE PRECISION :: x_coll(Nconst,NET), y_coll(Nconst,NET)
!    INTEGER :: i, j, deplac
!
!!    OPEN(40, FILE="donnees_init.txt", STATUS="OLD")
!!        DO j = 1, NET
!!            READ(40,*) Temp(j)
!!        END DO
!!        DO j = 1, NET
!!            READ(40,*) L(j)
!!        END DO
!!        DO j = 1, NET
!!            READ(40,*) V(j)
!!        END DO
!!        DO j = 1, Nconst
!!            DO i = 1, NET
!!                READ(40,*) x_coll(j,i)
!!            END DO
!!        END DO
!!        DO j = 1, Nconst
!!            DO i = 1, NET
!!                READ(40,*) y_coll(j,i)
!!            END DO
!!        END DO
!!    CLOSE(40)
!
!
!END SUBROUTINE
