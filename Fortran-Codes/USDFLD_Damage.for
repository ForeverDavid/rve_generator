c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C	
C	
C	- strainARRAY und stressARRAY eingefügt
C	- getvrm-Methodenaufruf?
C	- do-Schleife rausgenommen
C	- statOld - SDV (ARRAY)
C	- Fehlermeldungen rausgenommen, neue schreiben bei 
C	Aktualisierung der user-defined state variables?
C
C
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~









      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,
     3 LACCFLA)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),
     1 T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),
     1 COORD(*)
      
c define the array strainARRAY() and stressARRAY(), the output of fuction
c getvrm will be save in these two arrays
      DIMENSION STRAIN(15), STRESS(15)
C      user coding to define FIELD and, if necessary, STATEV and PNEWDT

CC	START

c define the variables
c
c tot_croft:total value of Cockcroft-Latham parameter
      double precision tot_croft
c tot_freud:total value of Freudenthal parameter
      double precision tot_freud
c tot_oyane:total value of Oyane parameter      
      double precision tot_oyane
c pi: 3.141593
      double precision pi
c I1,I2,I3: stress invariants
      double precision I1
      double precision I2
      double precision I3
c phi: the angle between principal plane and current plane[rad]      
      double precision phi
c princ1,princ2,princ3: principal stresses      
      double precision princ1
      double precision princ2
      double precision princ3
c max_princ:maximum principal stress      
      double precision max_princ
c delpeeq: delta peeq (peeq_new-peeq_old)
      double precision delpeeq
c peeq: equivalent strain      
      double precision peeq
c p,q: variables during caculation, no special meaning      
      double precision p
      double precision q
c s_mises:von mises stress
      double precision s_mises
c dam_croft:Cockcroft-Latham parameter increment in one frame
      double precision dam_croft
c dam_freud:Freudenthal parameter increment in one frame      
      double precision dam_freud
c dam_oyane:Oyane parameter increment in one frame      
      double precision dam_oyane
c ALEX: variable initialization
      double precision triaxiality
      double precision J2
      double precision J3
      double precision LodeParameter
      double precision LodeAngle
c define a logical variable "first" and initialize its value as "true"
      LOGICAL FIRST
      SAVE    FIRST
      DATA    FIRST /.TRUE./
c variable assignment 
      pi=3.141593
c
c use the fuction "getvrm" to get stress tensor (S) and equivalent 
c strain (PEEQ) from the material points.
c if jStatus is not equal to 0, it means the output is failed. And then
c an error message will be written in the message file.
      jrcd = 1

      call getvrm( 'PE', STRAIN, JARRAY, FLGRAY, JRCD, JMAC, 
     &		    	JMATYP, MATLAYO, LACCFLA)
C      if (jrcd .ne. 0 ) then
C	  call stdb_abqerr(-2,'Utility routine GETVRM '//
C     *      'failed to get strain.',0,zero,' ')
C	  call xit
C      end if
c
      call getvrm( 'S', STRESS, JARRAY, FLGRAY, JRCD, JMAC, 
     &	            	JMATYP, MATLAYO, LACCFLA)
C      if( jrcd .ne. 0 ) then
C         call stlb_abqerr(-2,'Utility routine GETVRM '//
C     *      'failed to get stress.',0,zero,' ')
C         call xit
C      end if

C use function getvrm to get values up to this point in time
      call getvrm( 'SDV', ARRAY, JARRAY, FLGRAY, JRCD, JMAC, 
     &			JMATYP, MATLAYO, LACCFLA)
C      if (jrcd .ne. 0 ) then
C	  call stdb_abqerr(-2,'Utility routine GETVRM '//
C     *      'failed to get SDV.',0,zero,' ')
C	  call xit
C      end if

C      WRITE(*,*) 'STRAIN(7) = ', STRAIN(7)
c      WRITE(*,*) 'STRESS(1) = ', STRESS(1)
c      WRITE(*,*) 'STRESS(2) = ', STRESS(2)
c      WRITE(*,*) 'STRESS(3) = ', STRESS(3)
c      WRITE(*,*) 'STRESS(4) = ', STRESS(4)
c      WRITE(*,*) 'SDV(1) = ', ARRAY(1)
c      WRITE(*,*) 'SDV(2) = ', ARRAY(2)
c      WRITE(*,*) 'SDV(3) = ', ARRAY(3)
c      WRITE(*,*) 'SDV(4) = ', ARRAY(4)
c      WRITE(*,*) 'SDV(5) = ', ARRAY(5)
c      WRITE(*,*) 'SDV(6) = ', ARRAY(6)
c      WRITE(*,*) 'SDV(7) = ', ARRAY(7)

      FIELD(1) = STRESS(1)
      FIELD(2) = STRESS(2)
      FIELD(3) = STRESS(3)
      FIELD(4) = STRESS(4)
      FIELD(5) = STRAIN(7)
 


c     
c check whether it is the first time to call this subroutine. If it
c is the first time, initialize the variables delpeeq, dam, totdam
      IF (FIRST) THEN
         FIRST = .FALSE.
         data delpeeq /0.d0/,dam_croft /0.d0/
         data tot_croft/0.d0/
         data tot_freud/0.d0/,dam_freud/0.d0/
         data tot_oyane/0.d0/,dam_oyane/0.d0/
         STATEV(2)=tot_croft
         STATEV(3)=tot_freud
         STATEV(4)=tot_oyane
      END IF
c the main loop to calculate damage parameter in each material point 


C do-Schleife löschen?    do k = 1, nblock  


c extract the equivalent strain at material point k from array strain()
         peeq = STRAIN(7)
c save peeq in the solution dependent variable sdv1
         STATEV(1)=peeq
c calculate the delta peeq

C 

         delpeeq=STATEV(1)-ARRAY(1)

c calculate the stress invariants
c the components of array stress() is:
c stress(k,1)=sigma11,stress(k,2)=sigma22,stress(k,3)=sigma33
c stress(k,4)=sigma12,stress(k,5)=sigma31,stress(k,6)=sigma23
         I1=STRESS(1)+STRESS(2)+STRESS(3)
         I2=STRESS(1)*STRESS(2)+STRESS(2)*STRESS(3)+
     &      STRESS(1)*STRESS(3)-STRESS(4)**2.0
     &      -STRESS(5)**2.0-STRESS(6)**2.0
         I3=STRESS(1)*STRESS(2)*STRESS(3)
     &      -STRESS(1)*(STRESS(6)**2)-STRESS(2)*(STRESS(5)**2) 
     &      -STRESS(3)*(STRESS(4)**2)
     &      +2.0*STRESS(4)*STRESS(5)*STRESS(6)
c
         p=(3.0*I2-I1**2.0)/3.0
         q=(9.0*I1*I2-2.0*(I1**3.0)-27.0*I3)/27.0
c calculate the princial stresses
c in order to avoid the occurrence of imaginary number during calculation,
c two constraints are added: firstly, radicand larger than 0.
c secondly, -1<cos(phi)<1 
         if ((p .LE. 0.0) .and. 
     & ((q*(-0.50)*((-(p**3)/27.0)**(-0.50))) .gt. -1.0) .and.
     & ((q*(-0.50)*((-(p**3)/27.0)**(-0.50))) .LE. 1.0)) then
c if thesse 2 conditions are satisfied, then calculation continues
         phi=ACOS(q*(-0.50)*((-(p**3)/27.0)**(-0.50)))
c calculate principal stresses:
c
         princ1=I1/3.0+2.0*SQRT(p/(-3.0))*COS(phi/3.0)
         princ2=I1/3.0-SQRT(p/(-3.0))*(COS(phi/3.0)-SQRT(3.0)
     &          *sin(phi/3.0))
         princ3=I1/3.0-SQRT(p/(-3.0))*(COS(phi/3.0)+SQRT(3.0)
     &          *sin(phi/3.0))
c find the maximum princial stress
         max_princ=MAX(princ1,princ2,princ3)
c calculate von mises equivalent stresses and hydrostatic stress
         s_mises=SQRT(((princ1-princ2)**2+(princ2-princ3)**2
     &            +(princ3-princ1)**2)/2.d0)
          s_h=(princ1+princ2+princ3)/3.d0
c calculate damage increment of Cockcroft-Latham, Freudenthal and Oyane
         dam_croft=max_princ*delpeeq
         dam_freud=s_mises*delpeeq
         dam_oyane=(1.d0+s_h/s_mises)*delpeeq
c calculate the total damage parameter 

C stateOld unbekannt (siehe oben: getvrm('SDV',...)

         tot_croft=ARRAY(2)
         tot_freud=ARRAY(3)
         tot_oyane=ARRAY(4)
c Cockcroft damage increment will be added to the total damage parameter 
c only when it is larger than 0.
      if (dam_croft .gt. 0.0) then   
         STATEV(2)=dam_croft+tot_croft
      else
         STATEV(2)=tot_croft
      end if
c calculate total frendthal parameter and Oyane parameter
         STATEV(3)=tot_freud+dam_freud
         STATEV(4)=tot_oyane+dam_oyane
         STATEV(5)=max_princ
         STATEV(6)=s_h
         STATEV(7)=delpeeq
      else
         STATEV(2)=ARRAY(2)
         STATEV(3)=ARRAY(3)
         STATEV(4)=ARRAY(4)
         STATEV(5)=ARRAY(5)
         STATEV(6)=ARRAY(6)
         STATEV(7)=ARRAY(7)
c
      end if


C      end do

c-----------------------------------------------------------------------
c Enhancement by Alexander Schowtjak:
c Compute triaxiality and Lode Parameter/Angle
c-----------------------------------------------------------------------

c Triaxiality (1.0D-12 is added to the denominator in order to avoid
c division by zero in case of zero stresses)
      triaxiality = s_h/(s_mises+1.0D-12)

c 2nd invariant of the deviatoric stress tensor
      J2 = s_mises**2/3.0

c 3rd invariant of the deviatoric stress tensor
      J3 = (princ1-s_h)*(princ2-s_h)*(princ3-s_h)

c Lode parameter
      LodeParameter = -27.0/2.0*J3/s_mises**3.0

c Lode Angle
      LodeAngle = acos((27.0/2.0*J3/s_mises**3.0))/3

      STATEV(8)  = triaxiality
      STATEV(9)  = LodeParameter
      STATEV(10) = LodeAngle
c

CC	ENDE


      RETURN
      END