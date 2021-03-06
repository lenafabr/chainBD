PROGRAM MAIN
  ! subroutines for testing parts of the code
  USE KEYS, ONLY : BDSTEPS, DELT, OUTFILE, BDPRINTEVERY, SNAPSHOTFILE, &
       & SNAPSHOTEVERY, APPENDSNAPSHOTS, DOBROWN, STARTEQUIL, NBEAD, NFIX, ACTION
  USE CHAINUTIL, ONLY : CHAIN, SETUPCHAIN, SETCHAINPARAMS, INITIALIZECHAIN, CLEANUPCHAIN
  USE BROWNDYN, ONLY : RUNBROWNDYNSIM

  IMPLICIT NONE
  TYPE(CHAIN), TARGET :: CHAINOBJ
  TYPE(CHAIN), POINTER :: CHAINP

  CHAINP=>CHAINOBJ

  CALL READKEY

  CALL SETUPCHAIN(CHAINP,NBEAD,NFIX)
  CALL SETCHAINPARAMS(CHAINP)
  CALL INITIALIZECHAIN(CHAINP, STARTEQUIL)
  
  IF (ACTION.EQ.'BROWNDYN') THEN
     CALL RUNBROWNDYNSIM(CHAINP, BDSTEPS, DELT,OUTFILE,bdPRINTEVERY,SNAPSHOTFILE,SNAPSHOTEVERY,APPENDSNAPSHOTS,DOBROWN)
  ELSE
     PRINT*, 'ERROR IN MAIN: unidentified action', ACTION
     STOP 1
  ENDIF

  CALL CLEANUPCHAIN(CHAINP)

END PROGRAM MAIN
