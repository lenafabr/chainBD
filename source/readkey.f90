SUBROUTINE READKEY
  ! this subroutine reads in keywords from a parameter file
  ! it sets the various global variables defined in KEYS module
  ! name of the parameter file is param.* where * is a keyword argument
  ! if no keyword argument is supplied, the default is just a file called param
  ! The EXTRAPARAMFILES keyword will allow extra parameter files to be 
  ! read in as well

  USE KEYS
  USE INPUTPARAMS, ONLY : READLINE, READA, READF, READI, READO
  USE GENUTIL

  IMPLICIT NONE

  ! ---- stuff for inputing the parameter file in free format --------
  CHARACTER*100 :: ARG ! command line argument
  INTEGER :: NUMARG ! number of command line arguments
  INTEGER :: NITEMS ! number of items on the line in the parameter file
  INTEGER :: PF ! input file unit
  LOGICAL :: FILEEND=.FALSE. ! done reading file?
  CHARACTER*100 :: WORD ! keyword
  ! -------------- for reading multiple parameter files --------  
  INTEGER, PARAMETER :: MAXNFILES = 10
  CHARACTER*100 :: PARAMFILES(MAXNFILES)
  INTEGER :: NPARAMFILES, NPARAMREAD
  ! ------ for initializing random number generator
  INTEGER :: TIMEVAL(8), SEED
  ! ---------------- temporary variables ---------------
  INTEGER :: DUMI, I, TMPI, DUMI1, DUMI2, DUMI3
  CHARACTER*100 :: DUMSTR
  LOGICAL :: LDUM

  ! ------------------------
  ! set variable defaults
  ! ------------------------
  ACTION = 'NONE'
  RNGSEED = 0
  VERBOSE = .FALSE.

  ! geometry and energy parameters 
  LS = 1D0; ! segment length
  LP = 1; ! persistence length
  ESTR = 1D3; ! stretching energy

  NBEAD = 10; ! number of beads in the chain

  NFIX = 0 !number of fixed beads
  
  ! external forces
  NEXTforce = 0  

  ! input/output  
  OUTFILE = '*.out'
  DUMPSNAPSHOTS = .FALSE. ! periodically dump chain snapshots
  SNAPSHOTEVERY = 1 ! how often to dump snapshots
  SNAPSHOTFILE = '*.snap.out' ! snapshot file
  APPENDSNAPSHOTS = .FALSE. ! append snapshots to file rather than replacing
  
  STARTEQUIL = .FALSE. ! start with equilibrated configuration
  
  ! brownian dynamics
  DELT = 1D-4 ! time step
  FRICT = 1D0 ! friction coefficient
  KT = 1D0 ! temperature
  BDSTEPS = 1000 ! number of brownian steps to run
  BDPRINTEVERY = 1 ! how often to print output
  DOBROWN = .TRUE. ! include brownian forces
  
  ! -------------------------
  ! Read in all parameter files, starting with the ones specified on command line
  ! --------------------------

  PF = 55 ! i/o unit number to be used for parameter files

  ! get input parameter files from command line
  NPARAMFILES = 0
  NUMARG = COMMAND_ARGUMENT_COUNT()  
  IF (NUMARG==0) THEN
     NPARAMFILES = 1
     PARAMFILES(1) = 'param'
     ARG = ''
  ELSE
     DO I = 1,NUMARG
        CALL GETARG(I, ARG)
        NPARAMFILES = NPARAMFILES + 1
        WRITE(DUMSTR,'(A)') 'param.' //TRIM(ADJUSTL(ARG))
        PARAMFILES(NPARAMFILES) = DUMSTR
     ENDDO
     ! reset arg to its original value
     IF (NUMARG.GT.1) CALL GETARG(1,ARG)
  ENDIF

  NPARAMREAD = 0 ! keep track of how many files have been read
  DO WHILE (NPARAMREAD.LT.NPARAMFILES)
     NPARAMREAD = NPARAMREAD + 1

     PRINT*, 'Reading parameter file: ', PARAMFILES(NPARAMREAD)
     INQUIRE(FILE=PARAMFILES(NPARAMREAD),EXIST=LDUM)
     IF (.NOT.LDUM) THEN
        PRINT*, 'ERROR in READKEY: Parameter file ', TRIM(ADJUSTL(PARAMFILES(NPARAMREAD))), ' does not exist.'
        STOP 1
     ENDIF
     OPEN(UNIT=PF, FILE=PARAMFILES(NPARAMREAD), STATUS='OLD')

     ! read in the keywords one line at a time
     DO 
        CALL READLINE(PF,FILEEND,NITEMS)
        IF (FILEEND.and.nitems.eq.0) EXIT

        ! skip empty lines
        IF (NITEMS.EQ.0) CYCLE

        ! Read in the keyword for this line
        CALL READA(WORD,CASESET=1)

        ! Skip any empty lines or any comment lines
        IF (WORD(1:1).EQ.'#') CYCLE

        SELECT CASE(WORD) ! pick which keyword
        CASE('ACTION')
           CALL READA(ACTION, CASESET=1)
        CASE('BDPRINTEVERY')
           CALL READI(BDPRINTEVERY)   
        CASE('BDSTEPS')
           CALL READI(BDSTEPS)                   
        CASE('DELT')
           CALL READF(DELT)        
        CASE('ESTR')
           CALL READF(ESTR)
        CASE('EXTFORCE')
           NEXTFORCE = NEXTFORCE + 1
           IF (NEXTFORCE.GT.MAXNEXTFORCE) THEN
              PRINT*, 'ERROR: too many external forces specified', NEXTFORCE
              STOP 1
           ENDIF
           CALL READI(EXTFORCEIND(NEXTFORCE))
           DO I = 1,3
              CALL READF(EXTFORCE(I,NEXTFORCE))
           ENDDO
        CASE('FIXBEADS')
           DO I = 1,NITEMS-1
              NFIX = NFIX + 1
              IF (NFIX.GT.MAXNFIX) THEN
                 PRINT*, 'ERROR: exceeded maximum allowed fixed beads', MAXNFIX
                 STOP 1
              ENDIF
              CALL READI(FIXBEADS(NFIX))
           ENDDO
        CASE('FRICT')
           CALL READF(FRICT)
        CASE('KT')
           CALL READF(KT)
        CASE('LP')
           CALL READF(LP)
        CASE('LS')
           CALL READF(LS)     
        CASE('NOBROWN')
           DOBROWN = .FALSE.
        CASE('NBEAD')
           CALL READI(NBEAD)          
        CASE('OUTFILE')
           CALL READA(OUTFILE)
        CASE('RNGSEED')
           CALL READI(RNGSEED)
        CASE('SNAPSHOTFILE')
           CALL READA (SNAPSHOTFILE)
        CASE('SNAPSHOTS')
           DUMPSNAPSHOTS = .TRUE.
           IF (NITEMS.GT.1) CALL READI(SNAPSHOTEVERY)
           IF (NITEMS.GT.2) CALL READA(SNAPSHOTFILE)
           IF (NITEMS.GT.3) CALL READO(APPENDSNAPSHOTS)
        CASE('STARTEQUIL') 
           ! start with properly sampled equilibrium conformations
           STARTEQUIL = .TRUE.
           IF (NITEMS.GT.1) CALL READO(STARTEQUIL)                
        CASE('VERBOSE')
           CALL READO(VERBOSE)
        CASE DEFAULT
           print*, 'ERROR: unidentified keyword ', TRIM(WORD), " Will ignore."
        END SELECT
     ENDDO
     CLOSE(PF)
  ENDDO

  ! -----------------
  ! check validity of some values, raise errors or adjust as necessary
  ! -----------------  
  
  IF (LS.LT.0) THEN
     PRINT*, 'ERROR IN LS VALUE', LS
     STOP 1
  ENDIF
  IF (LP.LT.0) THEN
     PRINT*, 'ERROR IN LP VALUE', LP
     STOP 1
  ENDIF
  IF (ESTR.LT.0) THEN
     PRINT*, 'ERROR IN ESTR VALUE', ESTR
     STOP 1
  ENDIF
  IF (NBEAD.LT.2) THEN
     PRINT*, 'ERROR IN NBEAD VALUE', NBEAD
     STOP 1
  ENDIF

  ! ----------- fix file names -----------
  CALL REPLACESUBSTR(OUTFILE,'*',TRIM(ADJUSTL(ARG)))
  CALL REPLACESUBSTR(SNAPSHOTFILE,'*',TRIM(ADJUSTL(ARG)))
  ! ---------------------------

  ! Initiate random number generator 
  IF (RNGSEED.EQ.0) THEN
     ! use the current time of day in milliseconds
     CALL DATE_AND_TIME(VALUES=TIMEVAL)
     SEED = TIMEVAL(5)*3600*1000 + TIMEVAL(6)*60*1000 + TIMEVAL(7)*1000 + TIMEVAL(8)
  ELSEIF (RNGSEED.EQ.-1) THEN
     ! use the last 5 characters in the command-line argument
     SEED = STRING2NUM(TRIM(ADJUSTL(ARG)))    
  ELSEIF (RNGSEED.EQ.-2) THEN
     ! use the last 4 characters in the command-line argument 
     ! and additionally the millisecond time 
     CALL DATE_AND_TIME(VALUES=TIMEVAL)
     SEED = STRING2NUM(TRIM(ADJUSTL(ARG)),TIMEVAL(8))
  ELSE
     ! use this seed directly
     SEED = RNGSEED
  ENDIF

  print*, 'Initiating Mersenne twister random number generator with seed:', SEED
  CALL SGRND(SEED)

  print*, '------------Parameter values : -------------------'
  print*, 'ACTION: ', TRIM(ADJUSTL(ACTION))
  print*, 'Output file: ', TRIM(OUTFILE)  
  IF (DUMPSNAPSHOTS) THEN
     PRINT*, 'Dumping snapshot every', SNAPSHOTEVERY,'steps. In file:', TRIM(ADJUSTL(SNAPSHOTFILE))
  ENDIF
  
  print*, 'Number of BEADS:', NBEAD
  print*, 'LS, LP, ESTR:', LS, LP, ESTR
  IF (NFIX.GT.0) THEN
     PRINT*, 'FIXED BEADS: ', FIXBEADS(1:NFIX)
  ENd IF
  DO I = 1,NEXTFORCE
     PRINT*, 'EXTERNAL FORCE ON: ', EXTFORCEIND(I), EXTFORCE(:,I)
  ENDDO
  
  PRINT*, 'Friction coefficint:', FRICT  
  IF (STARTEQUIL) PRINT*, 'Starting from equilibrated configurations.'  
  print*, '----------------------------------------------------'


END SUBROUTINE READKEY
