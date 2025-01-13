             SUBROUTINE UHYPER(BI1,BI2,AJ,U,UI1,UI2,UI3,TEMP,NOEL,
     1 CMNAME,INCMPFLAG,NUMSTATEV,STATEV,NUMFIELDV,FIELDV,
     2 FIELDVINC,NUMPROPS,PROPS)
C===============================================================
C     User defined hyperelastic material subroutine
C       for gel with Flory-Rehner free-energy function
C       to be used in Abaqus Standard
C     Fomulated and written by Aaron, April 20, 2022
C---------------------------------------------------------------
C     Material properties to be passed to the subroutine:
C       PROPS(1) - Nv
C       PROPS(2) - mu/kT
C       PROPS(3) - lambda_0   initial swelling
C     State variable:
C       TEMP     - temperature      
C       The initial value of chemical potential should agree
C           with the initial swelling, given as follows:
C       mu_0/kT = Nv*(1/lambda_0-1/detF0) + ln(1-1/detF0)
C                    + 1/detF0+ chi0/detF0^2+2*chi1/detF0^3
C       where detF0 = lambda_0^3
C     Output
C       Free-energy function U(I, J) and its derivatives
C       All free-energy density and stress given by the
C           calculation are normalized by kT/v
C         chi=A0+B0*temperature+(A1+B1*temperature)/AJ/detF0
C===============================================================
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION U(2),UI1(3),UI2(6),UI3(6),STATEV(*),FIELDV(*),
     1 FIELDVINC(*),PROPS(*)
      REAL(8) Nv, chi, lambda0, detF0, mu_kT, temperature, chi1
      Nv = PROPS(1)
      mu_kT = PROPS(2)
      lambda0 = PROPS(3)
      detF0 = lambda0**3
      A0=-12.947
      A1=17.92
      B0=0.0449
      B1=-0.0569
      temperature = TEMP ! TEMP is used to represent chemical potential
      chi1=A1+B1*temperature
      chi0=A0+B0*temperature
     

	U(1) = Nv/2 * ( 1/lambda0*BI1*AJ**(2.0/3.0) - 3/detF0
     & - 2/detF0*(3*LOG(lambda0) + LOG(AJ)) ) 
     & - chi0/detF0**2/AJ-chi1/detF0**3/AJ**2
     & - (AJ-1/detF0)*LOG(AJ/(AJ-1/detF0)) - mu_kT*(AJ-1/detF0) 
     & + chi0/detF0+chi1/detF0**2/AJ
	U(2) = 0
	UI1(1) = Nv/2/lambda0*AJ**(2.0/3.0)
	UI1(2) = 0
	UI1(3) = Nv/3/lambda0*BI1*AJ**(-1.0/3.0)
     & + (1-Nv)/AJ/detF0 - LOG(AJ/(AJ-1/detF0)) - mu_kT
     & + chi0/(detF0*AJ)**2+2*chi1/(detF0*AJ)**3
     & - chi1/detF0**2/AJ**2
	UI2 = 0
	UI2(3) = -Nv/9/lambda0*BI1*AJ**(-4.0/3.0) -(1-Nv)/AJ**2/detF0
     & + 1/AJ/(AJ-1/detF0)/detF0 - 2*chi0/detF0**2/AJ**3
     & - 6*chi1/detF0**3/AJ**4
     & + 2*chi1/detF0**2/AJ**3
      UI2(5) = Nv/3/lambda0 * AJ**(-1.0/3.0)
	UI3 = 0
	UI3(4) = -Nv/9/lambda0 * AJ**(-4.0/3.0)
	UI3(6) = 4*Nv/27/lambda0*BI1*AJ**(-7.0/3.0) + 2*(1-Nv)/AJ**3/detF0
     & - (2*AJ-1/detF0)/(AJ*(AJ-1/detF0))**2/detF0
     & + 6*chi0/detF0**2/AJ**4+24*chi1/detF0**3/AJ**5
     & - 6*chi1/detF0**2/AJ**4
      RETURN
      END