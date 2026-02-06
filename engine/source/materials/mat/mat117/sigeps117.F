Copyright>        OpenRadioss
Copyright>        Copyright (C) 1986-2026 Altair Engineering Inc.
Copyright>
Copyright>        This program is free software: you can redistribute it and/or modify
Copyright>        it under the terms of the GNU Affero General Public License as published by
Copyright>        the Free Software Foundation, either version 3 of the License, or
Copyright>        (at your option) any later version.
Copyright>
Copyright>        This program is distributed in the hope that it will be useful,
Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Copyright>        GNU Affero General Public License for more details.
Copyright>
Copyright>        You should have received a copy of the GNU Affero General Public License
Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
Copyright>
Copyright>
Copyright>        Commercial Alternative: Altair Radioss Software
Copyright>
Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
Copyright>        software under a commercial license.  Contact Altair to discuss further if the
Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
!||====================================================================
!||    sigeps117   ../engine/source/materials/mat/mat117/sigeps117.F
!||--- called by ------------------------------------------------------
!||    suser43     ../engine/source/elements/solid/sconnect/suser43.F
!||--- calls      -----------------------------------------------------
!||    finter      ../engine/source/tools/curve/finter.F
!||====================================================================
      SUBROUTINE SIGEPS117(
     1     NEL     ,NUPARAM ,NUVAR   ,JSMS    ,TIME    ,TIMESTEP,
     2     UPARAM  ,UVAR    ,AREA    ,OFF     ,OFFL    ,
     3     EPSZZ   ,EPSYZ   ,EPSZX   ,DEPSZZ  ,DEPSYZ  ,DEPSZX  ,
     4     SIGNZZ  ,SIGNYZ  ,SIGNZX  ,STIFM   ,DMELS   ,DMG     ,
     5     IPG     ,NFAIL   ,NGL     ,NFUNC   ,IFUNC   ,NPF      ,TF)    
C----------------------------------------------------------
C   I m p l i c i t   T y p e s
C----------------------------------------------------------
#include "implicit_f.inc"
C----------------------------------------------------------
#include "units_c.inc"
#include "comlock.inc"
#include "sms_c.inc"
C----------------------------------------------------------
C   D u m m y   A r g u m e n t s  
C----------------------------------------------------------
      INTEGER ,INTENT(IN)  :: NEL,NUPARAM,NUVAR,JSMS,IPG
      INTEGER ,INTENT(OUT) :: NFAIL
      INTEGER ,DIMENSION(NEL) ,INTENT(IN)  :: NGL
      my_real ,INTENT(IN)  :: TIME,TIMESTEP
      my_real ,DIMENSION(NEL) :: OFF,OFFL,AREA,DMELS,
     .   EPSZZ,EPSYZ,EPSZX,DEPSZZ,DEPSYZ,DEPSZX,SIGNZZ,SIGNYZ,SIGNZX
      my_real ,DIMENSION(NEL)  ,INTENT(OUT)   :: DMG
      my_real ,DIMENSION(NEL)  ,INTENT(INOUT) :: STIFM
      my_real ,DIMENSION(NUPARAM)   ,INTENT(IN)    :: UPARAM
      my_real ,DIMENSION(NEL,NUVAR) ,INTENT(INOUT) :: UVAR
C----------------------------------------------------------
      INTEGER NPF(*), NFUNC, IFUNC(NFUNC)
      my_real FINTER ,TF(*)
      EXTERNAL FINTER
C----------------------------------------------------------
C   L o c a l   V a r i a b l e s
C----------------------------------------------------------
      INTEGER :: I,II, IRUPT,NINDXF
      INTEGER ,DIMENSION(NEL)  :: INDXF
      my_real ::   DM,   DAM,  STF,DTB,
     .   E_ELAS_N,DTINV , T_N,T_S,
     .   EXP_G   ,     EXP_BK,   DELTA0S_CST, DELTA0N_CST, 
     .   GIC     ,     GIIC  ,  E_ELAS_S,  UND_CST,  UTD_CST,GAMA,
     .   IEVOL, ALPHA
c
      my_real, DIMENSION(NEL) :: EPSM,DELTA0M,  BETA,
     .   DELTAFMAX,EPST,FAC1, FAC2,FAC3,EPSMMAX ,EPSN ,
     .   DELTA0N,DELTA0S,UND,UTD,DYDX1,DYDX2,LENGTH,TMAX_N,TMAX_S
c----------------------------------------------------------
c     UVAR(1)  =  
c     UVAR(2)  =  
c     UVAR(3)  =  
c     UVAR(4)  =  
c     UVAR(5)  = 
c     UVAR(6)  =    
c     UVAR(7)  =  
c     UVAR(8)  =  
C=======================================================================
C     INPUT PARAMETERS INITIALIZATION
C=======================================================================
      E_ELAS_N      =  UPARAM(1)         ! Elastic stiffness for pure mode I
      E_ELAS_S      =  UPARAM(2)         ! Elastic stiffness for pure mode II   
      GAMA          =  UPARAM(3)
      T_N           =  UPARAM(4)         ! Peak traction for pure mode I
      T_S           =  UPARAM(5)         ! Peak traction for pure mode II
      IRUPT         =  INT(UPARAM(6))    ! Damage initiation criteria for mixed mode
      GIC           =  UPARAM(7)         ! Fracture energy mode I
      GIIC          =  UPARAM(8)         ! Fracture energy mode II
      NFAIL         =  INT(UPARAM(9))
      EXP_G         =  UPARAM(10)  
      EXP_BK        =  UPARAM(11)  
      DELTA0N_CST   =  UPARAM(13)        ! Separation in which damage starts for pure mode I
      DELTA0S_CST   =  UPARAM(14)        ! Separation in which damage starts for pure mode II
      UND_CST       =  UPARAM(15)        ! Ultimate separation in mode I(UND)
      UTD_CST       =  UPARAM(16)        ! Utimate displacement in tangential direction (UTD)
      IEVOL         = UPARAM(18)         ! Damage evolution type
      ALPHA         = UPARAM(19)         ! Exponential damage evolution exponent
c
      DTINV   = ONE / (MAX(TIMESTEP, EM20))
      STF     = E_ELAS_N + E_ELAS_S
      NINDXF  = 0
c-------------------------
      IF (TIME == ZERO) THEN
        DO I=1,NEL
         EPSMMAX(I)=ZERO      
        ENDDO 
      ELSE
        DO I=1,NEL
         EPSMMAX(I)    = UVAR(I,4)
        ENDDO 
      ENDIF

      DO I=1,NEL             
        EPSN(I) = MAX(EPSZZ(I) , ZERO)      
        EPST(I) = SQRT(EPSYZ(I)**2  + EPSZX(I)**2)
        EPSM(I) = SQRT( EPSN(I)**2  + EPST(I)**2)             
        EPSMMAX(I) = MAX(EPSM(I),EPSMMAX(I))
      ENDDO !I=1,NEL

      !---------------------------------------------------------
      !      Compute damage initiation and max displacement 
      !---------------------------------------------------------

      ! Check if element regularization is used. If yes, computed characteristic length
      IF ( IFUNC(1) /= 0.OR. IFUNC(2) /= 0) THEN 
         LENGTH(1:NEL) = SQRT(AREA(1:NEL))    
      ENDIF

      ! Element regularization for mode I
      IF ( IFUNC(1) /= 0) THEN
         DO I=1,NEL             
           TMAX_N(I) = T_N * FINTER(IFUNC(1),LENGTH(I),NPF,TF,DYDX1(I))
           DELTA0N(I)= TMAX_N(I) /E_ELAS_N
           UND(I)    = TWO*GIC /TMAX_N(I)
         ENDDO
      ELSE      ! No regularization, use input constant value
         DO I=1,NEL           
           DELTA0N(I) =  DELTA0N_CST
           UND(I)     =  UND_CST
         ENDDO  
      ENDIF

      ! Element regularization for mode II
      IF ( IFUNC(2) /= 0) THEN
         DO I=1,NEL             
           TMAX_S(I) = T_S * FINTER(IFUNC(2),LENGTH(I),NPF,TF,DYDX2(I))
           DELTA0S(I)= TMAX_S(I) /E_ELAS_S
           UTD(I)    = TWO*GIIC/TMAX_S(I)
         ENDDO
      ELSE        ! No regularization, use input constant value  
         DO I=1,NEL               
           DELTA0S(I) =  DELTA0S_CST
           UTD(I)     =  UTD_CST
         ENDDO  
      ENDIF

      !---------------------------------------------------------
      !       Update failure initiation strains in mixed mode
      !---------------------------------------------------------
      DO I=1,NEL             
        IF (EPST(I) == ZERO) THEN         ! If ESPT is zero, pure mode I
           DELTA0M(I) =  DELTA0N(I)
        ELSE IF (EPSN(I) == ZERO) THEN    ! If ESPN is zero, pure mode I
           DELTA0M(I) = DELTA0S(I)
        ELSE                              ! Mixed Mode (Quadratic initiation criteria)
           BETA(I)   = ABS(EPST(I) / EPSN(I))
           DELTA0M(I)= DELTA0S(I)* DELTA0N(I)*SQRT((ONE + BETA(I)**2)/           
     .                    ((DELTA0S(I)**2) + (BETA(I)* DELTA0N(I))**2))
        END IF
      ENDDO
!------------------------------------------------    
!       DAMAGE      : Ultimate Separation                       
!------------------------------------------------    
      ! IRUPT =2  = BK method
      IF (IRUPT == 2) THEN
        DO I=1,NEL
         IF (EPST(I) == ZERO) THEN
           DELTAFMAX(I)= UND(I)
         ELSE IF (EPSN(I) == ZERO) THEN
           DELTAFMAX(I)= UTD(I)
         ELSE      ! mixed mode 
            FAC1(I) = (E_ELAS_N**GAMA)/(ONE + BETA(I)**2)
            FAC2(I) = (E_ELAS_S**GAMA)*(BETA(I)**2)/(ONE + BETA(I)**2)
            FAC3(I) = (FAC1(I) + FAC2(I) )**(ONE/GAMA)
            !ultimate mixed mode displacement           
            DELTAFMAX(I)= TWO/(DELTA0M(I)*  FAC3(I)   )*
     .                    (GIC + (GIIC - GIC)*((E_ELAS_S*BETA(I)**2)/
     .                    (E_ELAS_N + E_ELAS_S*BETA(I)**2))**ABS(EXP_BK))
         END IF
        ENDDO
      ELSE
        ! IRUPT =1  = Power Law method (default)
        DO I=1,NEL
         IF (EPST(I) == ZERO) THEN          ! Pure mode I
           DELTAFMAX(I)= UND(I)
         ELSE IF (EPSN(I) == ZERO) THEN     ! Pure mode II
           DELTAFMAX(I)= UTD(I)
         ELSE                               ! Mixed Mode
           FAC1(I) = TWO*((ONE+BETA(I)**2))/DELTA0M(I)
           ! Ultimate mixed mode separation
           DELTAFMAX(I)= FAC1(I)*((E_ELAS_N/GIC)**EXP_G + 
     .                  (E_ELAS_S*(BETA(I)**2)/GIIC)**EXP_G)**(-ONE/EXP_G)
         END IF
        ENDDO
      ENDIF
!---------------------------------------------------------------
!       Damage Evolution      
!---------------------------------------------------------------
!          DELTAFMAX: Ultimate separation (UND, UTD or Mixed) 
!          ESPMMAX: Maximum separation
!          DELTA0M: Separation at damage onset

      ! Exponential Damage Evolution
      IF (IEVOL == 2) THEN

        DO I=1,NEL
          DM = EPSMMAX(I) - DELTA0M(I)
        IF (DM > ZERO.AND.EPSMMAX(I) /= ZERO ) THEN

          ! Exponential decay damage evolution
          DAM = ONE - (DELTA0M(I)/EPSMMAX(I)) *
     .      ( ONE - ( ONE - EXP( -ALPHA *
     .            (EPSMMAX(I)-DELTA0M(I)) /
     .            MAX((DELTAFMAX(I)-DELTA0M(I)), EM20) ) )
     .        / ( ONE - EXP(-ALPHA) ) )

          DMG(I) = MAX(DMG(I), DAM)!
          DMG(I) = MIN(DMG(I), ONE)

          ! OFFL = 1 means element is active (not eroded/failed yet)
          ! If element is still active and maximum separation > ultimate separation:
          IF (OFFL(I) == ONE .AND. EPSMMAX(I) > DELTAFMAX(I) ) THEN
            NINDXF = NINDXF+1       ! Counter for the number of failed elements, used in the output loop
            INDXF(NINDXF) = I       ! Array that stores the index of all failed elements
            OFFL(I) = ZERO          ! Sets OFFL of element I to zero -> Element has failed
          END IF
        END IF
      ENDDO

      ELSE    ! Linear Damage Evolution
        DO I=1,NEL
          DM = EPSMMAX(I) - DELTA0M(I)
          IF (DM > ZERO.AND.EPSMMAX(I) /= ZERO ) THEN
            DAM =  (DELTAFMAX(I) / EPSMMAX(I))*
     .             (EPSMMAX(I)   - DELTA0M(I))/
     .             MAX((DELTAFMAX(I) - DELTA0M(I)), EM20) ! Ensures no division by zero

            DMG(I) = MAX(DMG(I), DAM)!
            DMG(I) = MIN(DMG(I), ONE)
           IF (OFFL(I) == ONE .AND. EPSMMAX(I) > DELTAFMAX(I) ) THEN
             NINDXF = NINDXF+1
             INDXF(NINDXF) = I
             OFFL(I) = ZERO
           END IF
         END IF
        ENDDO

!---------------------------------------------------------------
!     Stress update
!---------------------------------------------------------------
      DO I=1,NEL
         ! Normal Direction
         IF (EPSZZ(I) < ZERO )  THEN                    ! Compression - No damage effect
            SIGNZZ(I) = E_ELAS_N*EPSZZ(I)            
          ELSE
            SIGNZZ(I) = (ONE-DMG(I))*E_ELAS_N*EPSZZ(I)  ! Damage reduces stress
         ENDIF

          ! Tangential direction - No sign check for shear, damage always have influence  
          SIGNYZ(I) = (ONE-DMG(I))*E_ELAS_S*(EPSYZ(I))  !  Sigma shear Y 
          SIGNZX(I) = (ONE-DMG(I))*E_ELAS_S*(EPSZX(I))  !  Sigma shear X 
      ENDDO
     
      DO I=1,NEL      
        UVAR(I,1)  = EPSZZ(I)
        UVAR(I,2)  = EPSZX(I)
        UVAR(I,3)  = EPSYZ(I)
        UVAR(I,4)  = EPSMMAX(I)

        UVAR(I,6)  = SIGNZZ(I)
        UVAR(I,7)  = SIGNZX(I)
        UVAR(I,8)  = SIGNYZ(I)

        UVAR(I,9)  = DELTA0M(I)
        UVAR(I,10) = BETA(I)
        UVAR(I,12) = DELTAFMAX(I)
      ENDDO
c-----------------------------------------------------      
c-----------------------------------------------------      
c     omega = sqrt(2k/2*dmels), dt=2/omega, 2*dmels=dt**2 * 2k / 4
      IF (IDTMINS==2 .AND. JSMS/=0) THEN
        DTB = (DTMINS/DTFACS)**2
        DO I=1,NEL                                                 
          DMELS(I)=MAX(DMELS(I),HALF*DTB*STF*AREA(I)*OFF(I))
        ENDDO                                                        
      END IF
      STIFM(1:NEL) = STIFM(1:NEL) + STF*AREA(1:NEL)*OFF(1:NEL)    
c-----------------------------------------------------      
      ! If number of failed elements higher than 0:
      IF (NINDXF > 0) THEN
        ! Loop over all failed elements
        DO II=1,NINDXF
          I = INDXF(II) ! Retrieves element index
#include "lockon.inc"
          ! NGL: Element Number; IPG: Integration Point; ESPM: Mixed Mode Strain value
          WRITE(IOUT ,1000) NGL(I),IPG,EPSM(I)
          WRITE(ISTDO,1100) NGL(I),IPG,EPSM(I),TIME
#include "lockoff.inc"
        END DO
      END IF
c-----------------------------------------------------      
 1000 FORMAT(5X,'FAILURE COHESIVE ELEMENT ',I10,
     .          ' INTEGRATION POINT',I2,', MIXED MODE STRAIN=',1PE16.9)
 1100 FORMAT(5X,'FAILURE COHESIVE ELEMENT ',I10,
     .          ' INTEGRATION POINT',I2,', MIXED MODE STRAIN=',1PE16.9,
     .          ' AT TIME ',1PE16.9)
c
      RETURN
      END
