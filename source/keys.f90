MODULE KEYS
  ! keyword parameters that are globally used in many different places in the code
  IMPLICIT NONE

   ! -------- General program control ---------------
  CHARACTER*100 :: ACTION
  INTEGER :: RNGSEED
  LOGICAL :: VERBOSE
  
  ! ----------------------
  ! chain geometry and energetics
  ! ---------------------
  INTEGER :: NBEAD
  DOUBLE PRECISION :: LP, LS, ESTR
  LOGICAL :: DOENDSTRETCH
  LOGICAL :: STARTEQUIL
  ! fixed beads
  INTEGER, PARAMETER :: MAXNFIX = 1000
  INTEGER :: FIXBEADS(MAXNFIX), NFIX
  ! external forces
  INTEGER, PARAMETER :: MAXNEXTFORCE = 1000
  INTEGER :: EXTFORCEIND(MAXNEXTFORCE), NEXTFORCE
  DOUBLE PRECISION :: EXTFORCE(3,MAXNEXTFORCE)
  
  ! ----------------------
  ! Output / input
  ! -----------------------
  CHARACTER*100 :: OUTFILE, SNAPSHOTFILE
  LOGICAL :: DUMPSNAPSHOTS, RESTART, APPENDSNAPSHOTS
  INTEGER :: SNAPSHOTEVERY  

  ! -----------------
  ! brownian dynamics
  ! -----------------
  DOUBLE PRECISION :: FRICT, DELT
  INTEGER :: BDSTEPS
  INTEGER :: BDPRINTEVERY
  LOGICAL :: DOBROWN
  DOUBLE PRECISION :: KT

END MODULE KEYS
