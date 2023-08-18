      MODULE umat_data
          INTEGER, PARAMETER :: MAX_ELEMENTS = 1  ! Adjust these as neededstr
          INTEGER, PARAMETER :: MAX_IPTS = 4  
          REAL(8) :: eeArray(MAX_ELEMENTS, MAX_IPTS, 4)  ! 4 is the size of ee
          REAL(8) :: epArray(MAX_ELEMENTS, MAX_IPTS, 4)
          REAL(8) :: eqplsArray(MAX_ELEMENTS, MAX_IPTS)
      END MODULE umat_data
C
      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)
      USE umat_data
C
      INCLUDE "ABA_PARAM.INC"
C
      CHARACTER*8 CMNAME
C
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),
     1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     2 PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), DROT(3, 3),
     3 DFGRD0(3, 3), DFGRD1(3, 3)
C LOCAL ARRAYS
C ----------------------------------------------------------------
C EELAS - ELASTIC STRAINS
C EPLAS - PLASTIC STRAINS
C FLOW - DIRECTION OF PLASTIC FLOW
C ----------------------------------------------------------------
C
      DIMENSION EELAS(NTENS),EPLAS(NTENS),FLOW(NTENS), HARD(NTENS)
      DIMENSION ee(4), ep(4), STRESSINITIAL(4), STRESSINCREMENT(4)
C
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.D0,
     1 ENUMAX=.4999D0, NEWTON=10, TOLER=1.0D-6)
C
C ----------------------------------------------------------------
C UMAT FOR ISOTROPIC ELASTICITY AND ISOTROPIC MISES PLASTICITY
C CANNOT BE USED FOR PLANE STRESS
C ----------------------------------------------------------------
C PROPS(1) - E
C PROPS(2) - NU
C PROPS(3..) - SYIELD AN HARDENING DATA
C CALLS UHARD FOR CURVE OF YIELD STRESS VS. PLASTIC STRAIN
C ----------------------------------------------------------------
C
C ELASTIC PROPERTIES
C
      WRITE(6,*) 'Entering UMAT at TIME=', TIME
      WRITE(6,*) 'NPT=', NPT
      IF (TIME == 0.0D0) THEN
        WRITE(6,*) 'Attempting to open file...'
        OPEN(UNIT=11,
     1     FILE='P:\Abaqus\RestartV4\initial_conditions_empty.txt',
     2     STATUS='OLD', ACTION='READ', IOSTAT=IO_STATUS)
        WRITE(6,*) 'IO_STATUS after OPEN=', IO_STATUS

        DO i = 1, NOEL
            DO j = 1, NPT
                READ(11, *) ee(1:4), ep(1:4), eqpls
                ! Store read values into the arrays
                eeArray(i, j, :) = ee
                epArray(i, j, :) = ep
                eqplsArray(i, j) = eqpls
                WRITE(6,*) 'ee', ee
                WRITE(6,*) 'ep=', ep
                WRITE(6,*) 'eqpls=', eqpls
            ENDDO
        ENDDO
        CLOSE(UNIT=11)
      ELSE
        ! Fetch stored values into local arrays for use
          ee = eeArray(NOEL, NPT, :)
          ep = epArray(NOEL, NPT, :)
          eqpls = eqplsArray(NOEL, NPT)
C
      ENDIF
      STRESSINITIAL=STRESS
      EMOD=PROPS(1)
      ENU=MIN(PROPS(2), ENUMAX)
      EBULK3=EMOD/(ONE-TWO*ENU)
      EG2=EMOD/(ONE+ENU)
      EG=EG2/TWO
      EG3=THREE*EG
      ELAM=(EBULK3-EG2)/THREE
C
C ELASTIC STIFFNESS
C
      DO K1=1, NDI
          DO K2=1, NDI
              DDSDDE(K2, K1)=ELAM
          END DO
          DDSDDE(K1, K1)=EG2+ELAM
      END DO
      DO K1=NDI+1, NTENS
          DDSDDE(K1, K1)=EG
      END DO
C
C CALCULATE PREDICTOR STRESS AND ELASTIC STRAIN
C
      DO K1=1, NTENS
          DO K2=1, NTENS
              STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*DSTRAN(K1)
          END DO
              ee(K1)=ee(K1)+DSTRAN(K1)
      END DO
C
C CALCULATE EQUIVALENT VON MISES STRESS
C
      SMISES=(STRESS(1)-STRESS(2))**2+(STRESS(2)-STRESS(3))**2
     1           +(STRESS(3)-STRESS(1))**2
      DO K1=NDI+1,NTENS
          SMISES=SMISES+SIX*STRESS(K1)**2
      END DO
      SMISES=SQRT(SMISES/TWO)
C
C GET YIELD STRESS FROM THE SPECIFIED HARDENING CURVE
C
      NVALUE=NPROPS/2-1
      CALL UHARD(SYIEL0, HARD, eqpls, EQPLASRT,TIME,DTIME,TEMP,
     1 DTEMP,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,CMNAME,NSTATV,
     2 STATEV,NUMFIELDV,PREDEF,DPRED,NVALUE,PROPS(3))
C
C DETERMINE IF ACTIVELY YIELDING
C
      IF (SMISES.GT.(ONE+TOLER)*SYIEL0) THEN
C
C     ACTIVELY YIELDING
C     SEPARATE THE HYDROSTATIC FROM THE DEVIATORIC STRESS
C     CALCULATE THE FLOW DIRECTION
C
          SHYDRO=(STRESS(1)+STRESS(2)+STRESS(3))/THREE
          DO K1=1,NDI
              FLOW(K1)=(STRESS(K1)-SHYDRO)/SMISES
          END DO
          DO K1=NDI+1, NTENS
              FLOW(K1)=STRESS(K1)/SMISES
          END DO
C
C SOLVE FOR EQUIVALENT VON MISES STRESS
C AND EQUIVALENT PLASTIC STRAIN INCREMENT USING NEWTON ITERATION
C
      SYIELD=SYIEL0
      DEQPL=ZERO
      DO KEWTON=1, NEWTON
          RHS=SMISES-EG3*DEQPL-SYIELD
          DEQPL=DEQPL+RHS/(EG3+HARD(1))
C
      CALL UHARD(SYIELD,HARD,eqpls+DEQPL,EQPLASRT,TIME,DTIME,TEMP,
     1    DTEMP,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,CMNAME,NSTATV,
     2    STATEV,NUMFIELDV,PREDEF,DPRED,NVALUE,PROPS(3))
          IF(ABS(RHS).LT.TOLER*SYIEL0) GOTO 10
      END DO
C
C WRITE WARNING MESSAGE TO THE .MSG FILE
C
      WRITE(7,2) NEWTON
    2   FORMAT(//,30X,'***WARNING - PLASTICITY ALGORITHM DID NOT ',
     1                 'CONVERGE AFTER ",I3," ITERATIONS')
   10 CONTINUE
C
C UPDATE STRESS, ELASTIC AND PLASTIC STRAINS AND
C EQUIVALENT PLASTIC STRAIN
C
      DO K1=1,NDI
          STRESS(K1)=FLOW(K1)*SYIELD+SHYDRO
          STRESSINCREMENT(K1)=STRESS(K1)-STRESSINITIAL(K1)
          ep(K1)=ep(K1)+THREE/TWO*FLOW(K1)*DEQPL
          ee(K1)=ee(K1)-THREE/TWO*FLOW(K1)*DEQPL
      END DO
      DO K1=NDI+1,NTENS
          STRESS(K1)=FLOW(K1)*SYIELD
          ep(K1)=ep(K1)+THREE*FLOW(K1)*DEQPL
          ee(K1)=ee(K1)-THREE*FLOW(K1)*DEQPL
          STRESSINCREMENT(K1)=STRESS(K1)-STRESSINITIAL(K1)
      END DO
      eqpls=eqpls+DEQPL
C
C CALCULATE PLASTIC DISSIPATION
C
      SPD=DEQPL*(SYIEL0+SYIELD)/TWO
C
C FORMULATE THE JACOBIAN (MATERIAL TANGENT)
C FIRST CALCULATE EFFECTIVE MODULI
C
      EFFG=EG*SYIELD/SMISES
      EFFG2=TWO*EFFG
      EFFG3=THREE/TWO*EFFG2
      EFFLAM=(EBULK3-EFFG2)/THREE
      EFFHRD=EG3*HARD(1)/(EG3+HARD(1))-EFFG3
      DO K1=1, NDI
          DO K2=1, NDI
          DDSDDE(K2, K1)=EFFLAM
          END DO
      DDSDDE(K1, K1)=EFFG2+EFFLAM
      END DO
      DO K1=NDI+1, NTENS
          DDSDDE(K1, K1)=EFFG
      END DO
      DO K1=1, NTENS
          DO K2=1, NTENS
              DDSDDE(K2, K1)=DDSDDE(K2, K1)+EFFHRD*FLOW(K2)*FLOW(K1)
          END DO
      END DO
      ENDIF
C
C STORE ELASTIC AND (EQUIVALENT) PLASTIC STRAINS
C IN STATE VARIABLE ARRAY
C
C
      DO i = 1, NOEL
          DO j = 1, NPT
                ! Store read values into the arrays
                eeArray(i, j, :) = ee
                epArray(i, j, :) = ep
                eqplsArray(i, j) = eqpls
          ENDDO
      ENDDO
      WRITE(6,*) 'eeArray', eeArray
      WRITE(6,*) 'epArray', epArray
      WRITE(6,*) 'eqplsArray', eqplsArray
C
      RETURN
      END
      SUBROUTINE UHARD(SYIELD,HARD,eqpls,EQPLASRT,TIME,DTIME,TEMP,
     1 DTEMP,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,
     2 CMNAME,NSTATV,STATEV,NUMFIELDV,
     3 PREDEF,DPRED,NVALUE,TABLE)
      INCLUDE "ABA_PARAM.INC"
      
      CHARACTER*80 CMNAME
      DIMENSION HARD(3),STATEV(NSTATV),TIME(*),
     1 PREDEF(NUMFIELDV),DPRED(*)
C
      DIMENSION TABLE(2, NVALUE)
C
      PARAMETER(ZERO=0.D0)
C
C SET YIELD STRESS TO LAST VALUE OF TABLE, HARDENING TO ZERO
C
      SYIELD=TABLE(1, NVALUE)
      HARD(1)=ZERO
C IF MORE THAN ONE ENTRY, SEARCH TABLE
C
      IF(NVALUE.GT.1) THEN
      DO K1=1, NVALUE-1
          EQPL1=TABLE(2,K1+1)
          IF(eqpls.LT.EQPL1) THEN
              EQPL0=TABLE(2, K1)
          IF(EQPL1.LE.EQPL0) THEN
              WRITE(7, 1)
    1         FORMAT(//, 30X, '***ERROR - PLASTIC STRAIN MUST BE ',
     1                        'ENTERED IN ASCENDING ORDER')
              CALL XIT
          ENDIF
C
C CURRENT YIELD STRESS AND HARDENING
C
              DEQPL=EQPL1-EQPL0
              SYIEL0=TABLE(1, K1)
              SYIEL1=TABLE(1, K1+1)
              DSYIEL=SYIEL1-SYIEL0
              HARD(1)=DSYIEL/DEQPL
              SYIELD=SYIEL0+(EQPLAS-EQPL0)*HARD(1)
              GOTO 10
            ENDIF
          END DO
   10     CONTINUE
      ENDIF
      RETURN
      END