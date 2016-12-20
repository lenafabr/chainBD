MODULE CHAINUTIL
  ! utilities for defining and dealing with chain object
  IMPLICIT NONE
  

  TYPE CHAIN 
     INTEGER :: NPT ! number of beads
     ! position for each bead
     ! 1st index is coordinate, 2nd index is bead number
     DOUBLE PRECISION, POINTER :: POS(:,:)
     DOUBLE PRECISION :: LS ! segment length
     ! energy parameters 
     ! LP = persistence length
     ! ESTR = stretch modulus
     DOUBLE PRECISION:: LP, ESTR
     ! friction coefficient associated with each bead
     DOUBLE PRECISION :: FRICT
     ! temperature
     DOUBLE PRECISION :: KT
     
     ! external forces applied to each bead
     ! 1st index is coordinate, 2nd index is bead number
     DOUBLE PRECISION, pointer :: FORCEXT(:,:)
     LOGICAL :: HASEXTFORCES = .FALSE.

     ! fixed beads
     INTEGER :: NFIX
     INTEGER, POINTER :: FIXBEADS(:)
     
     ! chain arrays have been allocated
     LOGICAL :: ARRAYSET = .FALSE.
  END type CHAIN

CONTAINS
  SUBROUTINE GETCHAINENERGY(CHAINP,ENERGY,FORCES,GETFORCE)
    ! get the overall energy associated with the current chain configuration,
    ! as well as the overall forces (-dE/dx) on the chain beads
    ! forces(I,J) has force in dimension i on jth bead
    ! if GETFORCE is false, calculate energy only
    IMPLICIT NONE

    TYPE(CHAIN), POINTER :: CHAINP
    DOUBLE PRECISION, INTENT(OUT) :: ENERGY, FORCES(3,CHAINP%NPT)
    LOGICAL, INTENT(IN) :: GETFORCE
    DOUBLE PRECISION :: ESTRETCH, EBEND, EEXT
    DOUBLE PRECISION, DIMENSION(3,CHAINP%NPT) :: FSTRETCH, FBEND, FEXT
    DOUBLE PRECISION :: TMP, FTMP(3), COEFF, DR(3), NDR, DR2(3), NDR2, NDR12
    DOUBLE PRECISION :: POS1(3), POS2(3)
    INTEGER :: S
    
    
    ENERGY = 0D0; ESTRETCH = 0D0; EBEND = 0D0; EEXT = 0D0
    IF (GETFORCE) THEN
       FORCES = 0D0
       FSTRETCH = 0D0
       FBEND = 0D0
       FEXT = 0D0
    ENDIF
    
    DO S = 1,CHAINP%NPT-1 ! go through each SEGMENT
       ! Calculate the stretching energy       
       POS1 = CHAINP%POS(:,S)
       POS2 = CHAINP%POS(:,S+1)       
             
       COEFF = CHAINP%ESTR/(2*CHAINP%LS)
       DR =POS2 - POS1
       NDR = SQRT(DOT_PRODUCT(DR,DR)) ! segment length
       TMP = NDR - CHAINP%LS
       ESTRETCH = ESTRETCH + COEFF*TMP**2

       IF (GETFORCE) THEN
          FTMP = 2*COEFF*TMP/NDR*DR
          FSTRETCH(:,S+1) = FSTRETCH(:,S+1)+FTMP
          FSTRETCH(:,S) = FSTRETCH(:,S) - FTMP
       END IF

       ! Calculate the bending energy
        IF (S.EQ.CHAINP%NPT-1) CYCLE             
        DR2 = CHAINP%POS(:,S+2)-CHAINP%POS(:,S+1)
        NDR2 = SQRT(DOT_PRODUCT(DR2,DR2))
        NDR12 = NDR*NDR2
        
        COEFF = CHAINP%LP/CHAINP%LS
        TMP = DOT_PRODUCT(DR,DR2)

        EBEND = EBEND + COEFF*(1-TMP/NDR12)

        IF (GETFORCE) THEN
           FBEND(:,S+2) = FBEND(:,S+2) - COEFF*(DR/NDR12 - TMP/NDR/NDR2**3*DR2)
           FBEND(:,S+1) = FBEND(:,S+1) - COEFF*(-DR/NDR12 + TMP/NDR/NDR2**3*DR2 + DR2/NDR12 - TMP/NDR2/NDR**3*DR)
           FBEND(:,S) = FBEND(:,S) - COEFF*(-DR2/NDR12 + TMP/NDR2/NDR**3*DR)
        ENDIF
    ENDDO

    ! calculate the energy due to external force
    IF (CHAINP%HASEXTFORCES) THEN
       EEXT = -SUM(CHAINP%FORCEXT*CHAINP%POS)

       IF (GETFORCE) THEN
          FEXT = -CHAINP%FORCEXT
       ENDIF
    END IF
    
    ENERGY = ESTRETCH+EBEND+EEXT
    FORCES = -FSTRETCH-FBEND-FEXT
    
  END SUBROUTINE GETCHAINENERGY

  SUBROUTINE OUTPUTSNAPSHOT(CHAINP,FILENAME,APPEND)
    ! dump out a snapshot of the chain configuration
    ! if append=true, append to an existing file
    TYPE(CHAIN), POINTER :: CHAINP
    CHARACTER(LEN=*), INTENT(IN) :: FILENAME
    LOGICAL, INTENT(IN) :: APPEND
    CHARACTER*100 :: FMT
    INTEGER :: B

    IF (APPEND) THEN
       OPEN(UNIT=99,FILE=FILENAME,POSITION='APPEND')
    ELSE
        OPEN(UNIT=99,FILE=FILENAME,POSITION='REWIND')
    ENDIF

    ! write information line
    WRITE(99,'(A,1X,I12)') 'C',CHAINP%NPT

    FMT = '(A,1X,12G20.10)'

    ! write out positions of chain beads
    DO B = 1,CHAINP%NPT     
       WRITE(99,FMT) 'A', CHAINP%POS(:,B)    
    ENDDO
  END SUBROUTINE OUTPUTSNAPSHOT  
  
  SUBROUTINE INITIALIZECHAIN(CHAINP,STARTEQUIL)
    ! initialize chain configuration
    ! if STARTEQUIL: select configuration from equilibrium
    ! else: straight configuration
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    LOGICAL, INTENT(IN) :: STARTEQUIL
    INTEGER :: S
    
     IF (.NOT.CHAINP%ARRAYSET) THEN
       PRINT*, 'ERROR IN INITIALIZECHAIN: chain arrays have not been allocated'
       STOP 1
    ENDIF

    CHAINP%POS = 0D0
    IF (STARTEQUIL) THEN
       CALL SAMPLEWLC(CHAINP%POS,CHAINP%NPT,CHAINP%LS,CHAINP%LP,CHAINP%ESTR,CHAINP%KT)
    ELSE
       DO S = 1,CHAINP%NPT-1
          CHAINP%POS(3,S+1) = CHAINP%POS(3,S)+CHAINP%LS
       ENDDO
    ENDIF
  END SUBROUTINE INITIALIZECHAIN

  SUBROUTINE SAMPLEWLC(COORDS,NBEAD,LS,LP,ESTR,KT)
    ! Sample from equilibrium distribution of wormlike chain
    USE GENUTIL, ONLY : GRND, RNORM, CROSS_PRODUCT, GETPERP, PI
    
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NBEAD
    DOUBLE PRECISION, INTENT(OUT) :: COORDS(3,NBEAD)
    DOUBLE PRECISION, INTENT(IN) :: LS, LP, ESTR, KT

    DOUBLE PRECISION :: C, D, U, UN, Z, RHO, PHI, ST
    DOUBLE PRECISION :: ZAX(3), XAX(3), YAX(3)
    INTEGER :: SC
    
    ! place initial bead at 0
    COORDS(:,1) = 0D0    
    
    ! parameter for the truncated exponential rho sampling
    C = LP/LS/KT
    ! parameter for the normal z sampling
    D = SQRT(LS*KT/ESTR)
    
    DO SC = 1,NBEAD-1  
       ! sample segment length from a normal distribution
       UN = RNORM()
       Z = UN*D + LS      

       ! sample phi uniformly
       PHI = 2*PI*GRND()
       
       IF (SC.EQ.1) THEN
          ! sample rho uniformly
          RHO = GRND()*2 - 1
          IF (ABS(RHO).GE.1D0) THEN             
             COORDS(:,2) = COORDS(:,1) + (/0D0,0D0,RHO/)*Z
          ELSE
             ST = SQRT(1-RHO*RHO)
             COORDS(:,2) = COORDS(:,1) + (/ST*COS(PHI), ST*SIN(PHI),RHO/)*Z
          ENDIF
       ELSE
          ! sample rho from a truncated exponential
          U = GRND()
          RHO = 1/C*LOG(EXP(-C)+2*U*SINH(C))

          IF (ABS(RHO).GE.1D0) THEN
             ST = 0D0
          ELSE
             ST = SQRT(1-RHO*RHO)
          END IF
          
          ! convert to absolute coords
          ZAX = COORDS(:,SC)-COORDS(:,SC-1)
          CALL GETPERP(ZAX,XAX)
          CALL CROSS_PRODUCT(ZAX,XAX,YAX)

          COORDS(:,SC+1) = COORDS(:,SC) + ST*COS(PHI)*XAX + ST*SIN(PHI)*YAX + RHO*ZAX
       ENDIF

       
    ENDDO
  END SUBROUTINE SAMPLEWLC
  
  SUBROUTINE SETCHAINPARAMS(CHAINP)
    ! set various parameters from the chain based on key arguments
    USE KEYS, ONLY : LS, LP, ESTR, FRICT, NBEAD, kt, FIXBEADS, NEXTFORCE, EXTFORCEIND, EXTFORCE
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER :: i
    
    IF (.NOT.CHAINP%ARRAYSET) THEN
       PRINT*, 'ERROR IN SETCHAINPARAMS: chain arrays have not been allocated'
       STOP 1       
    ELSEIF (CHAINP%NPT.NE.NBEAD) THEN
       PRINT*, 'ERROR IN SETCHAINPARAMS: wrong number of beads',CHAINP%NPT, NBEAD
       STOP 1
    ENDIF

    CHAINP%LS = LS
    CHAINP%LP = LP
    CHAINP%ESTR = ESTR
    CHAINP%FRICT = FRICT
    CHAINP%KT = KT

    ! external forces on specific beads
    CHAINP%FORCEXT = 0D0
    CHAINP%HASEXTFORCES = (NEXTFORCE.GT.0)
    
    IF (NEXTFORCE.GT.0) THEN    
       DO I = 1,NEXTFORCE
          CHAINP%FORCEXT(:,EXTFORCEIND(I)) = EXTFORCE(:,I)
       ENDDO       
    END IF

    ! fixed beads
    CHAINP%FIXBEADS = FIXBEADS(1:CHAINP%NFIX)
  END SUBROUTINE SETCHAINPARAMS

  SUBROUTINE SETUPCHAIN(CHAINP,NPT,NFIX)
    ! set up chain object with the given number of beads and given number of fixed beads
    ! initialize all arrays
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP
    INTEGER, INTENT(IN) :: NPT, NFIX

    CHAINP%NPT = NPT
    CHAINP%NFIX = NFIX
    CHAINP%HASEXTFORCES = .FALSE.

    ALLOCATE(CHAINP%POS(3,NPT), CHAINP%FORCEXT(3,NPT), CHAINP%FIXBEADS(NFIX))
    CHAINP%FORCEXT = 0D0

    CHAINP%ARRAYSET = .TRUE.
  END SUBROUTINE SETUPCHAIN

  SUBROUTINE CLEANUPCHAIN(CHAINP)
    ! deallocate arrays
    IMPLICIT NONE
    TYPE(CHAIN), POINTER :: CHAINP

    DEALLOCATE(CHAINP%POS, CHAINP%FORCEXT)

    CHAINP%ARRAYSET = .FALSE.
  END SUBROUTINE CLEANUPCHAIN
END MODULE CHAINUTIL
