!+
PROGRAM Profile
! ---------------------------------------------------------------------------
! PURPOSE -  Analysis and design of Low Speed Airfoils

! AUTHORS - Richard Eppler, Univ. of Stuttgart
!           Dan Somers, NASA Langley Research Center
!           Ralph Carmichael, Public Domain Aeronautical Software

! NOTE - John Roncz should perhaps be listed as an author, although
!   this particular source code is not his compressibility update.

IMPLICIT NONE


!... G L O B A L   C O N S T A N T S
  INTEGER,PARAMETER:: ILES=1,IDRU=2,DBG=3,FIG=4,ISTA=7  ! unit numbers
  REAL,PARAMETER:: PI=3.14159265, BOGEN=PI/180

!... G L O B A L   V A R I A B L E S
  REAL:: abfa = 1.0   ! scaling factor for TRA1 records; set by ABSZ command
  CHARACTER(LEN=12):: airfoilID   ! read from cols 1:4 of the TRA1 record
                                  ! or from next line following FXPR  
  REAL:: aln   ! zero lift angle of attack
  REAL:: cm
  REAL:: darf
  REAL:: darg
  REAL:: dlt,dltu   ! in degrees
  REAL:: dlv     ! was in /TRIT/  
  REAL:: eta
  INTEGER:: itit1,itit2   ! used in ALFA
  INTEGER:: itp   ! Trapro sets it to 1; Panel sets it to 2
  INTEGER:: izz   ! =0 for trap integration in Cint; >0 for 3rd order
  INTEGER:: jab
  INTEGER:: jst
!  INTEGER:: kfu     ! never used ???
!  INTEGER:: mgc = 0   ! never used ???
!  INTEGER:: mpl = 0   ! never used ???
  INTEGER:: mspli
!  INTEGER:: mtr = 0   ! never used ???
  INTEGER:: nal  ! number of angles of attack
  INTEGER:: nd     ! was in /TRIT/
  INTEGER:: nq   ! the number of points on the complete airfoil
  INTEGER:: nu      ! was in /TRIT/
!  REAL:: sma    ! never used
  REAL:: sump     ! was in /TRIT/
!  REAL:: xtf   ! never used
!  REAL:: xzeh,yzeh       ! /PLTM/   never used


!!      COMMON/TRIT/DLV,SUMP,XTRI(4),NU,ND,XSTPI(4),STHI(4)

!... G L O B A L   A R R A Y S
  REAL,DIMENSION(14):: agam = (/0.,1.,1.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,0. /)
  REAL,DIMENSION(14):: alca
  REAL,DIMENSION(29):: alfa   ! loaded by TRA1
  REAL,DIMENSION(29):: alfr
  REAL,DIMENSION(14):: alv    ! angles of attack
  REAL,DIMENSION(120,120):: amat   ! big matrix used in Panel for solution
  REAL,DIMENSION(28):: ani    ! loaded by TRA1
  REAL,DIMENSION(2):: cae   ! may not be used outside ProcessCommands
  REAL,DIMENSION(14):: cml,crl   ! not sure about these ????
  REAL,DIMENSION(5,2,14):: cw,su,sa
  REAL,DIMENSION(122):: ds
!  REAL,DIMENSION(60,7):: fuw  ! never used??
  REAL,DIMENSION(121,2):: gamma
  REAL,DIMENSION(13):: pures
  REAL,DIMENSION(60):: rs    ! ???
!  REAL,DIMENSION(4):: sthi   ! was in /TRIT/   ! never used
  REAL,DIMENSION(121):: x,y,arg,vf
  REAL,DIMENSION(121):: p,p1,xp,yp
  REAL,DIMENSION(121):: xf,yf,betaf    ! not sure these should be global
!  REAL,DIMENSION(4):: xstpi   ! was in /TRIT/               ! never used
  REAL,DIMENSION(4):: xtri   ! was in /TRIT/
!----------------------------------------------------------------------------
  CALL Welcome()

  CALL ProcessCommands()

  CALL Farewell()

  STOP

CONTAINS

!+                                                                        
SUBROUTINE CDCF(RET,H,CD,CF,H12)
! ---------------------------------------------------------------------------
! PURPOSE - Get the CD and CF coefficients
! called by GRUP
  REAL,INTENT(IN):: ret
  REAL,INTENT(IN):: h    ! boundary layer shape factor, delta-sub-3/delta-sub-2
  REAL,INTENT(OUT):: cd  ! boundary layer dissipation coefficient
  REAL,INTENT(OUT):: cf  ! boundary layer skin friction coefficient
  REAL,INTENT(OUT):: h12 ! boundary layer shape factor, delta-sub-1/delta-sub-2

!      COMMON/GRZK/CDK,AA(7),BB(7)
  REAL,PARAMETER:: CDK=0.01   ! was originally in /GRZK/
  REAL:: epst
  REAL:: rz
!----------------------------------------------------------------------------
  CALL H12B(H,H12,EPST)
  IF (h==0.0) RETURN

  IF (h < 0.0) THEN
    RZ = LOG((H12-1.)*RET) 
    CD = CDK*EXP(-.166692*RZ)                        ! eq 57
    CF = .045716*EXP(-1.26*H12 -.232*RZ)             ! eq 56
  ELSE
    CD = ((H*6.8377961- 20.521103)*H + 15.707952)/RET 
    CF = EPST/RET 
  END IF

  RETURN 
END Subroutine CDCF   ! -----------------------------------------------------

!+                                                                        
SUBROUTINE CINT(P,Z,NQ)
! ---------------------------------------------------------------------------
! PURPOSE - Integrate p over NQ ordinates and store the result in array Z
!   That is, z(k) is the integral up to p(k). z(nq)  is the total integral.
!   All points are assumed to be equidistant, the distance being 1.
!   You have to adjust the answer to the actual dimensions of the
!   independent variable.
!   If IZZ=0, the trapezoidal formula is used; otherwise a 3rd order formula
! called by TraPro
IMPLICIT NONE
  REAL,INTENT(IN),DIMENSION(:):: p   ! function to be integrated
  REAL,INTENT(OUT),DIMENSION(:):: z  ! resulting integral 
  INTEGER,INTENT(IN):: nq         ! size of p and z

!      DIMENSION P(NQ),Z(NQ)
  REAL:: fcint
  INTEGER:: n,nz
  REAL:: pl,pv,pvv,pz
!----------------------------------------------------------------------------
  IF (izz==0) THEN
    FCINT=0.
  ELSE
    FCINT=.08333333333
  END IF

  Z(1)=0. 
  PVV=P(NQ-1)    ! assumes p is periodic
  PV=P(1) 
  PL=P(2)

  DO N=2,NQ 
    NZ=N+1-(NQ-1)*(N/NQ) 
    PZ=P(NZ) 
    Z(N)=Z(N-1)+PL+PV+(PL+PV-PZ-PVV)*FCINT 
    PVV=PV 
    PV=PL 
    PL=PZ
  END DO

  RETURN 
END Subroutine Cint   ! -----------------------------------------------------

!+                                                                        
FUNCTION COSG(a) RESULT(f)
! ---------------------------------------------------------------------------
! PURPOSE - Cosine of an angle in degrees
IMPLICIT NONE
  REAL,INTENT(IN):: a   ! angle, in degrees
  REAL:: f              ! resulting cosine
! ---------------------------------------------------------------------------
  f = COS(A*BOGEN) 
  RETURN 
END Function CosG   ! -------------------------------------------------------                                          
                                                                        
!+                                                                        
FUNCTION CSLG(A,B) RESULT(f)
! ---------------------------------------------------------------------------
! PURPOSE - Evaluate a function that occurs frequently in airfoil theory.
IMPLICIT NONE
  REAL,INTENT(IN):: a,b   ! two angles, in degrees
  REAL:: f
  REAL:: xx
! ---------------------------------------------------------------------------
  xx=ABS(SIN((A-B)*BOGEN))
  IF (xx > 0.0) THEN        
    f = LOG(xx)
  ELSE
    f=0.0                            ! should never get here in this program
  END IF
  RETURN 
END Function CslG   ! -------------------------------------------------------

!+
SUBROUTINE DIA(X,Y,NP,D)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the maximum thickness of the airfoil.
! called several places in ProcessCommands and in TraPro
IMPLICIT NONE

  REAL,INTENT(IN),DIMENSION(:):: x,y
  INTEGER,INTENT(IN):: np
  REAL,INTENT(OUT):: d

  !!!    DIMENSION U(3),V(3),W(3),Z(3),A(3),B(3),X(121),Y(121)

  REAL,DIMENSION(3):: a,b
  REAL:: dnn
  REAL:: f
  INTEGER:: i
  INTEGER:: iit
  INTEGER:: n,nr
  INTEGER:: nm,nrm
  INTEGER:: nmz,nrmz
  INTEGER:: no,nu
  REAL:: ta
  REAL,DIMENSION(3):: u,v,w
  REAL:: xs
  REAL:: xst
  REAL:: xv
  REAL:: yd, yf, ys
  REAL,DIMENSION(3):: z
!----------------------------------------------------------------------------
      D = 0. 
      XV=X(1)     ! never used ???????

  DO N=2,NP 
    NR = NP-1 
    XS=X(N) 
  1 IF(X(NR)-XS)3,2,2 
  2 NR = NR - 1 
    GO TO 1 
  3 IF(NR-N.LE.3) EXIT
    YD=(Y(NR+1)-Y(NR))/(X(NR+1)-X(NR)) 
    DNN=(Y(N)-Y(NR)-(X(N)-X(NR))*YD)*(1.-.5*YD*YD) 
    IF(DNN.LE.D)GO TO 4 
    NM=N 
    NRM=NR 
    YF=YD 
    D=DNN 
 4  XV=XS
  END DO

  IIT=0 
  GO TO 8 
7 IIT=IIT+1 
  IF (IIT.GT.20) GO TO 11 
  DO I=1,3 
    NO=NM+2-I 
    NU=NRM-2+I 
    U(I)=TA*(X(NO)+YF*Y(NO)) 
    V(I)=TA*(Y(NO)-YF*X(NO)) 
    W(I)=TA*(X(NU)+YF*Y(NU)) 
    Z(I)=TA*(Y(NU)-YF*X(NU))
  END DO

  CALL QIP(U,V,A) 
  CALL QIP(W,Z,B) 
  IF(ABS(A(3)-B(3)).LT..0001) GO TO 11

  XST = (B(2)-A(2))*.5/(A(3)-B(3)) 
  YS = A(2)+2.*A(3)*XST 
  NMZ=NM 
  NRMZ=NRM 
  F=1.5 
  IF (F*(U(3)-XST).LT.XST-U(2)) NM=NM-1 
  IF (F*(XST-U(1)).LT.U(2)-XST) NM=NM+1 
  IF (F*(W(3)-XST).LT.XST-W(2)) NRM=NRM+1 
  IF (F*(XST-W(1)).LT.W(2)-XST) NRM=NRM-1 
  IF (NMZ.NE.NM .OR. NRMZ.NE.NRM) GO TO 7 
  IF (ABS(YS) .LT. .0001) GO TO 9 
  YF=YF+YS/TA 
8 TA=1./SQRT(1.+YF*YF) 
  GO TO 7 
9 D = A(1)-B(1)+(A(2)-B(2)+(A(3)-B(3))*XST)*XST 
11 RETURN 
END Subroutine Dia   ! ------------------------------------------------------


!+
SUBROUTINE Diagram(istim,menv,alphaMin,alphaMax)
! ---------------------------------------------------------------------------
! PURPOSE - Make a plot of the airfoil plus the velocity distribution or
!   the pressure distribution
! called by ProcessCommands

  INTEGER,INTENT(IN):: istim   ! plot mode    (needs work)
  INTEGER,INTENT(IN):: menv
  REAL,INTENT(IN):: alphaMin,alphaMax

  INTEGER:: i,m,n
  REAL:: vmx,xmx
  REAL:: v,vq
!----------------------------------------------------------------------------
  WRITE(FIG,*) "PORTRAIT"
  WRITE(FIG,*) "VIEWPORT 0.1 0.7  0.3 0.8"
  IF (menv==1) THEN
    WRITE(FIG, '(A,2F12.4)' ) "WINDOW  -0.3 -4.3 ", alphaMin,alphaMax

  ELSE
    WRITE(FIG,*) "WINDOW  0 1  0.4 2.0"
    WRITE(FIG,*) "LGRID  5 8  1 1  x/c  V/Uinf"
  END IF

  DO M=1,nal
    vmx=0.0
    DO n=1,nq
      v=ABS(Vpr(n,m))
      vq=v*v
      IF (vmx < vq) vmx=vq
      IF (xmx < v) xmx=v
      p(n)=v
    END DO
    IF (menv/=1) THEN
      WRITE(FIG,*) "DATA ", nq, "   for angle# ", m
      p(1:nq) = MAX(0.4, p(1:nq))   ! don't look at values < 0.4
      WRITE(FIG,'(2ES15.5)') (x(i),p(i),i=1,nq)
    END IF
    p1(m)=vmx-1.0
  END DO

  IF (menv==1) THEN
    WRITE(FIG,*) "DATA ", nal
    WRITE(FIG,'(2ES15.5)') (p1(i),alv(i),i=1,nal)
  ELSE
    WRITE(FIG,*) "VIEWPORT 0.1 0.7   0 0.3"
    WRITE(FIG,*) "WINDOW  0 1  -0.25 0.25"
    WRITE(FIG,*) "DATA ",nq
    WRITE(FIG,'(2ES15.5)') (x(i),y(i), i=1,nq)
  END IF

  WRITE(FIG,*) "SCMOVE 0.375 0.96"
  WRITE(FIG,*) "CTEXT DESIGN & ANALYSIS OF LOW SPEED AIRFOILS"
  WRITE(FIG,*) "ENDFRAME ---------------------------------------------------"

  RETURN
END Subroutine Diagram   ! --------------------------------------------------

!+
SUBROUTINE DRAW(WC,WS,WL,DK,DM,FL,AG,MA)
! ---------------------------------------------------------------------------
! PURPOSE -
IMPLICIT NONE
  REAL,INTENT(OUT):: wc
  REAL,INTENT(OUT):: ws
  REAL,INTENT(OUT):: wl
  REAL,INTENT(IN):: dk   !   WHEN MKP = 0, DK = KAPPA; OTHERWISE, DK = K.
  REAL,INTENT(IN):: dm
  REAL,INTENT(IN):: fl
  REAL,INTENT(IN):: ag
  INTEGER,INTENT(IN):: ma

  REAL:: beta
  REAL:: bqm1
  REAL:: cosp,sinp
  REAL:: fk
  REAL:: u
  REAL:: wf
  REAL:: wubeq
!----------------------------------------------------------------------------
      COSP = COSG(AG*FL) 
      FK= DK 
      IF(MA)2,1,2 
    1 FK= FK*(1.-COSP)/(1.+COSP) 
    2 SINP = SING(AG*FL) 
      WL = -DM*LOG(1.+FK) 
      IF(FL.EQ.0.) FK = 1. 
      BETA =(COSP-1.)/FK + COSP 
      BQM1 = BETA**2 - 1. 
      WUBEQ= SQRT(ABS(BQM1)) 
      U= (1.+BETA)*SINP/(1.+COSP) 
      IF(BQM1)3,5,4 
    3 WF = LOG(ABS((WUBEQ+U)/(WUBEQ-U))) 
      GO TO 6 
    4 WF = 2.* ATAN(U/WUBEQ) 
      GO TO 6 
    5 WF = 0. 
    6 WC =(WUBEQ*WF - SINP - BETA*FL*AG*1.745329E-2)*DM 
      WS = (COSP-1.)*(1.-(1./FK + 1.)*LOG(1.+FK))*DM 
  RETURN 
END SUBROUTINE Draw   ! -----------------------------------------------------

!+
SUBROUTINE Farewell()
! ---------------------------------------------------------------------------
! PURPOSE - Perform any final operations; then write a final message and
!   return to the main program and STOP
! called by main program

  CHARACTER(LEN=*),PARAMETER:: FINAL_MESS = "Normal termination of Profile"
!----------------------------------------------------------------------------
  CLOSE(UNIT=IDRU)
  CLOSE(UNIT=DBG)
  WRITE(*,*) FINAL_MESS

  RETURN
END Subroutine Farewell   ! --------------------------------------------------

!+
SUBROUTINE FIXLES()
! ---------------------------------------------------------------------------
! PURPOSE - Read the coordinates of a user-supplied airfoil. This routine
!   just reads the data; it does not make any calcs - see Panel for that.
! called by ProcessCommands

  INTEGER:: mlow   ! number of points on lower surface
  INTEGER:: mup    ! number of points on upper surface
  INTEGER:: nt
!----------------------------------------------------------------------------
  READ(ILES,'(A)') airfoilID 
  READ(ILES,'(2I5)') MUP,MLOW 
  NQ=MUP+MLOW-1   ! one less because the nose is there twice. nq is global

  READ(ILES,'(8F10.5)') x(mup:1:-1)     !      CALL RDC(X,MUP,-1,1)
  READ(ILES,'(8F10.5)') y(mup:1:-1)     !      CALL RDC(Y,MUP,-1,1) 
  READ(ILES,'(8F10.5)') x(mup:nq)       !      CALL RDC(X,MUP,1,NQ) 
  READ(ILES,'(8F10.5)') y(mup:nq)       !      CALL RDC(Y,MUP,1,NQ)

  nt=nq/12
  dlt=y(nt)/(BOGEN*(1.0-x(nt)))           ! dlt and dltu are global
  nt=nq-nt+1                              ! used by __________
  dltu=-y(nt)/(BOGEN*(1.0-x(nt)))
  RETURN
END Subroutine FixLes   ! ---------------------------------------------------

!+
SUBROUTINE FLAP(X,Y,BETA,NQ,XDA,YDA,ARCL,ARCLU,DEFL,XF,YF,BETAF,NQF)
! ---------------------------------------------------------------------------
! PURPOSE - Generate a new set of airfoil coordinates that incorporates the
!   deflection of a simple flap.
! called by ProcessCommands
  REAL,INTENT(IN),DIMENSION(:):: X,Y,BETA    ! original coordinates
  INTEGER,INTENT(IN):: NQ                    ! original number of points
  REAL,INTENT(IN):: XDA,YDA
  REAL,INTENT(IN):: ARCL,ARCLU
  REAL,INTENT(IN):: DEFL
  REAL,INTENT(OUT),DIMENSION(:):: XF,YF,BETAF  ! new coordinates
  INTEGER,INTENT(OUT):: NQF                    ! new number of points

  REAL:: anq
  REAL:: aa,al
  REAL:: alr
  REAL:: aq
  REAL:: cdl,sdl
  REAL:: cdf,sdf
  REAL:: dbt
  REAL:: deit
  REAL:: dl
  REAL:: dtx,dty
  REAL:: dx,dy

  REAL:: dxs,dls
  REAL:: dxlv

  REAL:: etada
  REAL:: ga,gas

  INTEGER:: l
  INTEGER:: nd  ! not the global variable

  INTEGER:: m
  INTEGER:: n
  INTEGER,DIMENSION(2,2):: nab
  INTEGER:: ndp
  INTEGER:: nou
  INTEGER:: np
  INTEGER,DIMENSION(2):: npa
  INTEGER:: nr,nv
  INTEGER:: nsuch

  REAL:: tg,tgs
  REAL:: u,v
  REAL:: xbr
  REAL:: xida
  REAL:: xin
  REAL:: xit
  REAL:: xkp
  REAL:: xst
  REAL:: XVDP,YVDP
!----------------------------------------------------------------------------
      DO 2 N=1,NQ 
      XF(N)=X(N) 
      YF(N)=Y(N) 
    2 BETAF(N)=BETA(N)

      NQF=NQ 
      NSUCH=NQ/2

  DO 22 M=1,2
    AQ=1.E10
    DO N=1,NSUCH
      NP=N+(M-1)*NSUCH 
      DTX=XDA-X(NP) 
      DTY=YDA-Y(NP) 
      ANQ=DTX*DTX+DTY*DTY 
      IF (ANQ.GE.AQ) CYCLE
      NDP=NP 
      AQ=ANQ 
      DBT=ATAN2(DTY,DTX)-BETA(NP) 
    END DO

    IF(SIN(DBT).GE.0.)NDP=NDP-1
    XVDP=X(NDP)
    YVDP=Y(NDP)
    DX=X(NDP+1)-XVDP
    DY=Y(NDP+1)-YVDP
    DLT=ATAN2(DX,-DY)
    GA=BETA(NDP)-DLT
    GAS=DLT-BETA(NDP+1)
    TG=SIN(GA)/COS(GA)
    TGS=SIN(GAS)/COS(GAS)
    U=XDA-XVDP
    V=YDA-YVDP
    SDL=SIN(DLT)
    CDL=COS(DLT)
    DL=SQRT(DX*DX+DY*DY)
    XIDA=(U*SDL-V*CDL)/DL
    ETADA=(U*CDL+V*SDL)/DL
    XIT=XIDA
 12 XKP=1.-XIT
    DEIT=ETADA-XIT*XKP*(XKP*TG+XIT*TGS)
    XIN=XIDA+DEIT*(XKP*(1.-3.*XIT)*TG+XIT*(3.*XKP-1.)*TGS)
    IF(ABS(XIN-XIT).LT..00001)GO TO 14
    XIT=XIN
    GO TO 12
 14 AA=DL*SQRT((XIT-XIDA)*(XIT-XIDA)+DEIT*DEIT)
    AL=ARCL+SIN(.5*DEFL)*SIGN(AA,DX*DEIT)
    IF(ARCLU.GT..00001.AND.DX.GT.0.)AL=AL-ARCL+ARCLU
    IF(AL.LT..0001*X(1))AL=.0001*X(1)
    XST=XKP
    DO 20 L=1,2
      ND=3-2*L 
      ALR=-XST*DL-AL 
      NR=NDP+L-1 
   16 NV=NR+ND 
      DXS=X(NV)-X(NR) 
      DLS=SQRT((Y(NV)-Y(NR))**2+DXS*DXS) 
      IF(ALR.GE.0.)GO TO 18 
      ALR=ALR+DLS 
      NR=NR-ND 
      DXLV=DXS/DLS 
      GO TO 16 
   18 XBR=X(NV)+ALR*DXLV 
      NOU=2*M-3 
      CALL PADD(XF,YF,BETAF,NQF,0,NOU,XBR) 
      NAB(M,L)=NV+2*M-1 
      XST=XIT
   20 END DO
22 END DO


  CDF=COS(DEFL)
  SDF=SIN(DEFL)
  M=1
  N=1
  DO 30 L=1,NQF
    XF(N)=XF(L)
    YF(N)=YF(L)
    BETAF(N)=BETAF(L)
    IF(L.EQ.NAB(M,1))NPA(M)=N
    IF(L.GT.NAB(1,1).AND.L.LT.NAB(2,2))GO TO 26
    BETAF(N)=BETAF(N)-DEFL
    DX=XF(N)-XDA
    DY=YF(N)-YDA
    XF(N)=XDA+DX*CDF+DY*SDF
    YF(N)=YDA-DX*SDF+DY*CDF
 26 IF(L.EQ.NSUCH)M=2
    IF((L-NAB(M,1))*(L-NAB(M,2)).GE.0)N=N+1
30 END DO

  NQF=N-1
  CALL PADD(XF,YF,BETAF,NQF,NPA(2),3,XF(1))
  CALL PADD(XF,YF,BETAF,NQF,NPA(1),3,XF(1))
  RETURN 
END Subroutine Flap   ! -----------------------------------------------------

!+
SUBROUTINE GAUSS(A,B,N,M,EPS)
! ---------------------------------------------------------------------------
! PURPOSE - Solve a system of linear equations by Gaussian elimination
! called by Gauss
  REAL,INTENT(IN OUT),DIMENSION(:,:):: a,b
  INTEGER,INTENT(IN):: n,m
  REAL,INTENT(IN):: eps

  REAL:: epsLocal
  INTEGER:: i,j,k
  REAL:: h,z
!----------------------------------------------------------------------------
  epsLocal = ABS(EPS) 
  DO J=1,N 
    H = .0 
    DO I=J,N             ! look for largest element of a(j:n,j)
      Z = ABS(A(I,J)) 
      IF (Z.LT.H) CYCLE
      K = I 
      H = Z 
    END DO 
    IF (H < epsLocal) THEN
      WRITE(IDRU,*) " SINGULARITY ALARM IN SUBROUTINE GAUSS"
      RETURN
    END IF

    IF (K /= J) THEN
      DO I=J,N 
        H = A(J,I) 
        A(J,I) = A(K,I) 
        A(K,I) = H
      END DO

      DO I=1,M 
        H = B(J,I) 
        B(J,I) = B(K,I) 
        B(K,I) = H
      END DO
    END IF

    DO I=1,N 
      IF (I.EQ.J) CYCLE
      H = A(I,J)/A(J,J) 
      DO K=J,N 
        A(I,K) = A(I,K)-A(J,K)*H
      END DO 
      DO K=1,M 
        B(I,K) = B(I,K)-H*B(J,K)
      END DO
    END DO 
    H = 1.0/A(J,J) 
    DO I=J,N 
      A(J,I) = A(J,I)*H
    END DO
    DO I=1,M 
      B(J,I) = B(J,I)*H
    END DO
  END DO 
  RETURN 
END Subroutine Gauss   ! ----------------------------------------------------                                          

!+
FUNCTION GetDateTimeStr() RESULT(s)
! ---------------------------------------------------------------------------
! PURPOSE - Return a string with the current date and time
!   This is now standard Fortran. It should work on ALL compilers!
!   You can change the first I2.2 below to I2 if you don't like a leading
!   zero on early morning times, i.e., 6:45 instead of 06:45
  IMPLICIT NONE
  CHARACTER(LEN=*),PARAMETER:: MONTH='JanFebMarAprMayJunJulAugSepOctNovDec'
  CHARACTER(LEN=*),PARAMETER:: FMT = '(I2.2,A1,I2.2,I3,A3,I4)'
  CHARACTER(LEN=15):: s
  INTEGER,DIMENSION(8):: v
!----------------------------------------------------------------------------
  CALL DATE_AND_TIME(VALUES=v)

  WRITE(s,FMT) v(5), ':', v(6), v(3), MONTH(3*v(2)-2:3*v(2)), v(1)
  RETURN
END FUNCTION GetDateTimeStr   ! ---------------------------------------------

!+
SUBROUTINE GRP(NAX,RE,MU,JR,ISTIFT)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the boundary-layer development for any given velocity
!   distribution and Reynolds number.
! called by ProcessCommands

  INTEGER,INTENT(IN):: nax   ! the number of angles of attack
  REAL,INTENT(IN),DIMENSION(:):: re  ! the aray of Reynolds numbers
  INTEGER,INTENT(IN),DIMENSION(:):: mu  ! the array of transition modes
  INTEGER,INTENT(IN):: jr     ! the number of Reynolds numbers for each alpha
  INTEGER,INTENT(IN):: istift  ! the pen number (not used)


!      COMMON CW(5,2,14),SA(5,2,14),SU(5,2,14),HP(32),RP(32),PUFF(22),   &
!     &AGAM(6),X(121),Y(121),DS(122),GAP(299),IDU,KDU,NQ,LDU(3),         &
!     &GAPS(382),BETA(121),GAPQ(8),GAMMA(242),AMAT(102,102)              

! inputs/outputs to GRS are in /BL/. Watch out for different spelling!
  real:: HVGL   ! H32 at the beginning of the step
  REAL:: D2i    ! delta-sub-2 at the beginning of the step
  REAL:: UK     ! U at the beginning of the step
  REAL:: UK1    ! U at the end of the step
  REAL:: dsz    ! length of the time step divided by the chord
  REAL:: WR     !
  INTEGER:: mt,ma   ! suction mode
  REAL:: V,V1   !
  REAL:: HR     ! H32 at the end of the step
  REAL:: D2R    ! delta-sub-2 at the end of the step
  REAL:: XA     ! the length within the step that is not separated
  REAL:: XU     ! the length within the step that is laminar
  REAL:: DCQ    ! the contribution of the step to the suction coefficient
  REAL:: USTR   !
  real:: VSTR   !      
  COMMON/BL/HVGL,D2i,Ur,UK1,Dsz,WR,Mt,MA,V,V1,HR,D2R,XA,XU,DCQ,USTR,VSTR                                                


!      COMMON/PRAL/DLT,DLTU,ALN,ALV(14),NAL,ITP,CML(14),CRL(14),         &
!     &DARF,ITIT1,ITIT2,D1I(121,5),AVI(5,5),BVI(5,5),                    &
!     & FMACH,BEMACH,IVI,NAIT(5),XDA,YDA,DEFLG,MOMAG,NAI

!      COMMON/TRIT/DLV,SUMP,XTRI(4),NU,ND,XSTPI(4),STHI(4)


!      DIMENSION RE(5),MU(5),
!      dimension XIT(121),YIT(121),BETI(121),        &
!     &CAE(2),WRE(5),H(5),D2(5),CAD(5),SAR(5),SUR(5),CWD(5),DD(5),CMV(5),&
!     &CPTR(5),BBL(5),CD(5),LD(5),LR(5),ACHI(5),CDLTX(2),&
!     &REVG(5),REDR(5)


  CHARACTER(LEN=*),PARAMETER,DIMENSION(2):: ALTX = (/ &
    "ZERO-LIFT LINE", "CHORD LINE    " /)
  CHARACTER(LEN=*),PARAMETER,DIMENSION(3):: D2TX = (/ &
    " DELTA2", "RDELTA2", " DELTA1" /)
!  CHARACTER,PARAMETER:: RETX = "R = ", STXT = " S"
  CHARACTER(LEN=*),PARAMETER,DIMENSION(2):: SFTX = (/"UPPER","LOWER"/)

  REAL,PARAMETER:: BBLI=0.03, HBBL=-1.6

!  REAL,DIMENSION(5):: cae
  REAL,DIMENSION(5):: wre,h,d2,cad,sar,sur,cwd,dd,cmv
  REAL,DIMENSION(5):: cptr,bbl,cd
!  REAL,DIMENSION(5):: achi,cdlt,xrevg,redr
  CHARACTER(LEN=1),DIMENSION(5):: ld
  INTEGER,DIMENSION(5):: lr



  REAL:: alwr

  REAL:: bbf
!  REAL:: beta

!  REAL:: cadz
!  REAL:: calok
!  REAL:: cdr
!  REAL:: cfr
!  REAL:: clit
!  REAL:: cmb  
!  REAL:: csl,cso,csp
!  REAL:: cwdz
  REAL:: cwpm


!  REAL:: d1a
!  REAL:: d2tr
!  REAL:: d1p,d1pp,d1t,d1u
!  REAL:: dca
!  REAL:: dre
!  REAL:: dsp,dsr,dsrdsz
  REAL:: dsr


  REAL:: fkt=1.0   ! factor for plotting H32 on vertical axis
!  REAL:: gamma

  REAL:: h12r
!  REAL:: htr
  REAL:: hablc
  REAL,DIMENSION(32):: hp,rp   ! OR, should this be GLOBAL
!  REAL:: htx

  INTEGER:: i,j,k
  INTEGER:: ia
!  INTEGER:: idc
  INTEGER:: ip
!  INTEGER:: is
!  INTEGER:: itz,izp,izt


!  INTEGER:: jsl
  INTEGER:: ju

  INTEGER:: kdtx
!  INTEGER:: kst
!  INTEGER:: lj


  INTEGER:: mag

  INTEGER:: nala,nale
  INTEGER:: ndsd
  INTEGER:: nend
  INTEGER:: nkr
!  INTEGER:: nub
  INTEGER:: nundsd
!  INTEGER:: nup,nupp
  
!  REAL:: red
  REAL:: rpl
!  REAL:: rth
  REAL:: qts
  REAL:: s
!  REAL:: slm
!  REAL:: stg

  REAL:: ur
  REAL:: utr

!----------------------------------------------------------------------------
      MA=0 
      V=0. 
      V1=0. 
      KDTX=INT(AGAM(6))-1    ! print code, set by ProcessCommands (RE)
                             ! -1 no print; 0 summary print;
                             ! >0 print that RE
      MAG=INT(AGAM(8))       ! plot code, set by ProcessCommands (RE)
                             ! =0 no plots
      NKR=NQ-1 
      HABLC = 2.5 

  DO J=1,JR
    h(j)=1.0
    LR(J)=NINT(RE(J)*1E-3)    ! Reynolds # in thousands
    WRE(J)= SQRT(RE(J))
  END DO

  NALA = 1
  NALE = ABS(NAX) 
  IF (NAX.LT.0) NALA = NALE 

  DO 38 IA = NALA, NALE    ! the really BIG loop
    IF (MAG.EQ.0.OR.IA.NE.1) GO TO 3 
!      PL=ISTIFT.GE.5 
!      CALL DEFINE(PL,390.,200.,1.43,1.95,.4,4.4) 
!      CALL ACHS(4,1.46,1.457,1.67,10.,10000.,10.,1.,TRUE) 
!      CALL ACHS(1,1.,.96,.92,1.46,1.67,.01,.05,TRUE) 
!      CALL ACHS(4,1.71,1.707,1.92,10.,10000.,10.,1.,TRUE) 
!      CALL ACHS(1,1.,.96,.92,1.71,1.92,.01,.05,TRUE) 
!      HP(1)=1.51509 
!      HP(2)=1.51509 
!      HP(3)=1.672
!      RP(1)=1. 
!      RP(2)=2.6656 
!      RP(3)=3.5039
!      CALL POLZUG(1,HP,RP,3,.TRUE.,TRUE) 
!  hp(1:3)=hp(1:3)+0.2
!      CALL POLZUG(1,HP,RP,3,.TRUE.,TRUE)

 3 continue
    DO 39 ju=1,2      ! upper and lower surfaces  (also a big loop)
      IP=0
      ND = 2*JU-3         !  -1 +1
      NDSD=JU-ITP         !   0  1      or -1 0
      NEND=1+NKR*NDSD     !   1  nkr+1

      DO K=1,NKR
        IF (VPR(K,IA).LE. 0.0) EXIT      ! find first negative vpr
      END DO 
      nu=k+ju-2
      IF (ABS(VPR(NU,IA)) .LT. 0.1) nu=nu+nd
      nundsd=nu-ndsd
      DSR=DS(NUNDSD)*VPR(NU,IA)/(VPR(NU,IA)-VPR(NU-ND,IA))  ! ds is global
      S=0. 
      UK=0. 

      bbl(1:jr)=0.0

      IF (KDTX > 0) THEN
        ALWR=ALV(IA)-DARF*ALN 
        WRITE(IDRU,*) " BOUNDARY LAYER    AIRFOIL "//airfoilID
        WRITE(IDRU,'(A,F6.2,A)' ) "   ALPHA=", ALWR, &
          " deg. relative to the "//ALTX(itit2)
        WRITE(IDRU,14) SFTX(JU),(LR(J),MU(J),J=1,JR) 
     14 FORMAT(1X,A," SURFACE   ",  5("    R=",I6,"000 MU =",I2)) 
        WRITE(IDRU,22) (J,D2TX(KDTX),J=1,JR) 
     22 FORMAT("       S        U ",5(I5,"    H32",4X,A)) 
        GO TO 26
      END IF

   24 NUNDSD=NU-NDSD      ! comes back here frequently
      DSR = DS(NUNDSD)

   26 UK1=ABS(VPR(NU,IA)) 
      S = S + DSR 

      DO 27 J=1,JR   ! the number of Reynolds numbers
        DSZ = DSR 
        UR = UK 
        HVGL=H(J) 
        MT=MU(J) 
        D2I=D2(J) 
        WR=WRE(J) 
        CALL GRS()   ! advance the solution one step
        IF (KDTX==1) dd(j)=d2r
        IF (KDTX==2) DD(J)=uk1*re(j)*d2r*1.E-6 
        IF (KDTX==3) THEN
          CALL H12B (HR,H12R,QTS) 
          DD(J)=D2R*H12R
        END IF
        IF (UK==0.0) THEN
          SAR(J)=0.0
          SUR(J)=0.0
        END IF
        SAR(J)=SAR(J)+XA 
        IF (HVGL*HR < 0.0) THEN
          utr=uk+(uk1-uk)*xu/dsr
          CPTR(J)=1.0-UTR*UTR 
          HVGL=-1.51509
        END IF
        SUR(J)=SUR(J)+XU 
        IF (HVGL.GT.HBBL .AND. HR.LE.HBBL) &
              BBL(J)=s-sur(j)-(dsr-xu)*(hbbl-hr)/(hvgl-hr)
        H(J)=HR 
        D2(J)=D2R
   27 END DO

      IF (mag <= 0) GO TO 31
      IF (ABS(h(mag))>1.62) GO TO 31
      rpl=fkt*LOG(uk1*re(mag)*d2(mag))   !!!!! fkt is never set ?????????
      IF (rpl<1.0 .OR. rpl>4.0) GO TO 31
      ip=ip+1
   29 rp(ip)=rpl
      hp(ip)=ABS(h(mag)) + 0.2*REAL(ju-1)   ! ju=1 upper; =2 lower
      !!! writing omitted




   31 IF (KDTX > 0) WRITE(IDRU,32) S,UK1,(H(K),DD(K),K=1,JR)
   32 FORMAT (F9.5,F8.3,5(F13.4,F10.6))    ! wow, 132 wide
   35 IF (NU==NEND) GO TO 36   ! this is the way out of the loop
      uk=uk1
      nu=nu+nd
      GO TO 24     ! go back for another step


   36 continue
      IF (ip > 0) THEN
        WRITE(FIG,*) "DATA ", ip
        WRITE(FIG,'(2F12.6)' ) (hp(i),rp(i), i=1,ip)
      END IF

      IF (KDTX==3) WRITE(IDRU,'(17X,5(6X,"CP(TRANS) =",F6.2))') cptr(1:jr)

      DO  J=1,jr
        CALL H12B(H(J),H12R,QTS) 
        IF(H12R.GT.HABLC)H12R = HABLC 
        IF(H(J).LT.0. .AND. BBL(J).EQ.0.0) BBL(J)=sur(j)
        BBF=2. 
        IF(bbl(j) > bbli) BBF=-2. 
        CW(J,JU,IA)=(UK1**(2.5+.5*H12R))*D2(J)*BBF 
        SA(J,JU,IA)= S - SAR(J) 
        SU(J,JU,IA)= S - SUR(J)
      END DO

 39 end do
38 END DO   ! end of the big loop
  i=0   ! not sure why.... 


! might want to skip this title if nax==0 OR test on KDTX
  WRITE(IDRU,'(///A)') " SUMMARY       AIRFOIL "//airfoilID
  WRITE(IDRU,*) "  ANGLE OF ATTACK RELATIVE TO THE "//altx(itit2)
  WRITE(IDRU, '("   ALPHA0=" F7.3, " DEGREES")' ) aln
  WRITE(IDRU,'(AF7.3)') " * indicates bubble analog longer than", bbli
  WRITE(IDRU,49) (LR(j),mu(j),j=1,jr)
49 FORMAT(6X, 5("   R=",I6,"000   mu=",I2))

  DO 60 i=nala,nale
    ALWR=ALV(I)-DARF*ALN 
    WRITE(IDRU, '(" ALPHA =",F6.2," DEGREES")' ) alwr
    WRITE(IDRU, '(5X,5(I5," S TURB S SEP    CD"))' )  (j,j=1,jr)

    DO K=1,2
      DO J=1,JR
        CWPM=CW(J,K,I) 
        IF (cwpm < 0.0) THEN
          LD(j)="*"         ! bubble indicator
        ELSE
          LD(j)=" "         ! no bubble
        END IF
        CD(J)=ABS(CWPM)
      END DO
      WRITE(IDRU,54) SFTX(K),(SU(J,K,I),SA(J,K,I),CD(J),LD(J),J=1,JR)
   54 FORMAT(2X,A,5(F9.4,F7.4,F7.4,A))
    END DO

    DO J = 1,JR
      CALL VISC(I,J,CAD(J),CWD(J),CMV(J))
    END DO  ! was loop 56

    WRITE(IDRU,58) (CAD(J),CWD(J),J=1,JR) 
  58 FORMAT("  TOTAL", 5("   CL=",F6.3," CD=",F7.4))
    WRITE(IDRU,102) CMV(1:jr)
  102 FORMAT (6X, 5(10X,"CM=",F7.4)) 
60 END DO 

  RETURN
END Subroutine GRP ! --------------------------------------------------------

!+
SUBROUTINE GRS()
! ---------------------------------------------------------------------------
! PURPOSE - Compute one integration of the boundary layer equations (42)
!    and (43). GRS can also compute the effect of suction on the boundary
!    layer
! called by GRP

! NOTE - the original coding used x as a variable. I renamed it xx
! to avoid confusion with the global array x, even though this is
! not really necessary
IMPLICIT NONE

! inputs/outputs to GRS are in /BL/  ! watch the spelling!
      real:: HVGL   ! H32 at the beginning of the step
      REAL:: D2     ! delta-sub-2 at the beginning of the step
      REAL:: UK     ! U at the beginning of the step
      REAL:: UK1    ! U at the end of the step
      REAL:: DL     ! length of the time step divided by the chord
      REAL:: WRE
      integer:: mu   ! transition mode
      INTEGER:: ma   ! suction mode
      REAL:: V,V1
      REAL:: HR     ! H32 at the end of the step
      REAL:: D2R    ! delta-sub-2 at the end of the step
      REAL:: XA     ! the length within the step that is not separated
      REAL:: XU     ! the length within the step that is laminar
      REAL:: DCQ    ! the contribution of the step to the suction coefficient
      REAL:: USTR
      real:: VSTR
  COMMON/BL/HVGL,D2,UK,UK1,DL,WRE,MU,MA,V,V1,HR,D2R,XA,XU,DCQ,USTR,VSTR                                                

!      integer:: nu,nd
!      real:: DLV,SUMP,XTRI(4),XSTPI(4),STHI(4)
!      COMMON/TRIT/DLV,SUMP,XTRI,NU,ND,XSTPI,STHI


  REAL:: bit
!  REAL:: cd
!  REAL:: cf
  REAL:: d2e,d2m,d2v
  REAL:: dd2,dd2m
!  REAL:: dd2s
  REAL:: dz,dzm
  REAL,PARAMETER:: EAP=0.000005
  REAL:: esp,epst
  REAL,PARAMETER:: EUM=0.0001
  REAL:: fum,fume
  REAL:: h
  REAL:: hap
  REAL:: h12
  REAL:: hdiff
  REAL:: hv
  REAL:: HE=-1.0
  REAL:: hm
  REAL:: hs

  INTEGER:: isg,isgv
!  INTEGER:: ist
  INTEGER:: istab
!  INTEGER:: kast
!  INTEGER:: nxtr


!  REAL:: ret
!  REAL:: sth
  REAL:: stp   ! set but never used
  REAL:: ue
!  REAL:: uh
  REAL:: ukv,uk1v
  REAL:: um
  REAL:: vg
  REAL:: vm
  REAL:: vv
!  REAL:: wcf
  REAL:: wrev
!  REAL:: xan
!  REAL:: xen
  REAL:: xs
!  REAL:: xsf,xstf
  REAL:: xx
  REAL:: z,zm
!----------------------------------------------------------------------------
      BIT=0. 
      H=HVGL 
      IF(MA)3,101,3 
  101 V1 = 0. 
      V = 0. 
      IF(HE.LE.0..OR.UK.LE.0.)GO TO 3
    1 IF(ABS(UK-UKV)+ABS(UK1-UK1V)+ABS(DL-DLV)+ABS(H-HV)-.1E-6)2,3,3 
    2 IF(ABS(D2*WRE-D2V*WREV).GE..1E-6) GO TO 3 
      BIT = 1. 
      D2E = D2E*WREV/WRE

    3 D2V=D2 
      HV = H 
      UKV= UK 
      UK1V = UK1 
      DLV = DL 
      WREV = WRE 
      DCQ = 0. 
      XA = DL 
      XU = DL 
      VV = V 
      IF(DL)4,4,5 
    4 HE=H 
      D2E= D2 
      GO TO 50 
    5 USTR =(UK1-UK)/DL 
      IF(MA.LT.4)VSTR =(V1-V)/DL 
      IF(UK)6,7,6 
    6 IF((UK1-UK)/UK - 1.)8,8,7 
    7 D2E = .290043/(WRE*SQRT(USTR)) 
      HE = 1.619977 
      GO TO 50 
    8 IF(D2)9,9,10 
    9 D2E =(.664108/WRE)*SQRT(2.*DL/(UK1+UK)) 
      HE = 1.572584 
      GO TO 50 
   10 xx = 0.                 ! was x
      SUMP=xx                 ! was x
      ISTAB = 0 
      IF(H.LT.0.)XU = 0. 
      IF(BIT)13,13,11 
   11 CALL UMP(UK,D2,H,FUM) 
      GO TO 40 
   13 IF (xx +DL -DLV) 12,12,19    ! was x
   12 XS = xx                      ! was x
      UK = UKV +  XS*USTR 
      VG= VV + VSTR*XS 
      IF(H)14,14,16 
   14 HAP = 1.46 
!      IF(HAP+EAP.LE.-H.OR.UK1.GE.UK)GO TO 20  ! uk1>=uk from John Roncz
      IF (hap+eap+h <= 0.0) GO TO 20
      XA=xx                         ! was x
      CALL H12B(H,H12,EPST) 
      D2E = D2*(UK/UK1)**((5.+H12)*.5) 
      HE = H 
      GO TO 50 
   16 HAP = 1.515095 
      SUMP=xx                   ! was x
      IF(HAP+EAP-H)17,18,18 
   17 CALL UMP(UK,D2,H,FUM) 
      IF(FUM+EUM)20,18,18 
   18 H = -H 
      XU = xx       ! was x 
      ISTAB = 0 
   19 DL = DLV - xx    ! was x 
      GO TO 12 
   20 CALL GRUP(H,D2,UK,Z,DZ,DD2,VG) 
      ESP = .002 
      IF(MA-3)22,22,21 
   21 ESP = .001 
      IF(xx .EQ.0.)V = VG      ! was x
   22 XS = XS + .5*DL 
      D2M = D2 + .5*DD2 
      ZM = Z + .5*DZ 
      UM = UKV +USTR*XS 
      ISG = 1 
      IF(D2M)25,25,23 
   23 ISG = 2 
      HM = ZM/D2M 
      IF(HM - 2.)30,25,25 
   25 ISGV = ISG 
      IF(ISTAB.GE.9)GO TO 49
   26 ISTAB = ISTAB +1 
      DL = .5*DL 
      GO TO 13 
   30 IF(HM - HAP+EAP)31,32,32 
   31 DL = .5*DL*(ABS(H)-HAP)/(ABS(H)-HM) 
      GO TO 13 
   32 IF(H.LT.0.)HM = -HM 
      IF(MA.LT.4)VM= VV+ XS*VSTR 
      CALL GRUP(HM,D2M,UM,ZM,DZM,DD2M,VM)  ! returns zm,dzm,dd2m,vm
      ISG = 3 
      D2E = D2 +DD2M 
      IF(D2E)25,25,33 
   33 HE= (Z+DZM)/D2E 
      ISG = 4 
      IF(HE-2.)34,25,25 
   34 HDIFF = HE-ABS(H) 
      ISG = 5 
      IF(ABS(HDIFF)-.01)35,35,25 
   35 HS = ABS(ABS(H)-2.*ABS(HM)+HE) 
      ISG = 6 
      IF(HS-ESP)36,36,25 
   36 IF(HE-HAP+EAP)37,38,38 
   37 DL = DL*(ABS(H)-HAP)/(ABS(H)-HE) 
      GO TO 13 
   38 IF(H)39,39,40 
   39 HE = -HE 
      GO TO 42 
   40 UE = UK + USTR*DL 
      SUMP=xx+DL      ! was x
      CALL UMP(UE,D2E,HE,FUME) 
      IF(FUME-EUM)42,42,41 
   41 DL = DL*FUM/(FUM-FUME) 
      GO TO 13 
   42 DCQ = DCQ + DL*VM 
      xx = xx + DL     ! was x
      IF(xx -DLV)47,50,50    ! was x
   47 H = HE 
      D2 = D2E 
      GO TO(43,13,43,45,13,43),ISGV 
   43 IF(HS-.1*ESP)44,13,13 
   44 DL = 2.*DL 
      ISTAB =ISTAB-1 
      GO TO 13 
   45 IF(HDIFF)44,13,13
   49 stp=1/(istab-9)    ! set but never used
   50 HR = HE 
      D2R = D2E 
  RETURN
END Subroutine GRS   ! ------------------------------------------------------

!+
SUBROUTINE GRUP(H,D2,U,Z,DZ,DD2,V)
! ---------------------------------------------------------------------------
! PURPOSE - Compute right hand sides of equations (42) and (43).
!   Equations written as delta-prime-sub-2=RHS and delta-prime-sub-3=RHS
!   Calls CDCF to get CD and CF
! called by GRS (two places)

IMPLICIT NONE
  REAL,INTENT(IN):: h     !
  REAL,INTENT(IN):: d2    !
  REAL,INTENT(IN):: u     !
  REAL,INTENT(OUT):: z    !
  REAL,INTENT(OUT):: dz   ! right-hand-side of (43)
  REAL,INTENT(OUT):: dd2  ! right-hand-side of (42)
  REAL,INTENT(OUT):: v    !
  
!  real:: cdk
  real,dimension(7):: aa=0.0,bb=0.0  ! some sort of suction parameters
                                     ! an undocumented feature
                                     ! never used unless ma /= 0
!!!      COMMON/GRZK/CDK,AA,BB   ! where are aa,bb set ???

  real:: dum(4),dl,wre,du(7),us,dumy
  integer:: md,ma
      COMMON/BL/DUM,DL,WRE,MD,MA,DU,US,DUMY



  REAL:: b
  REAL:: cd
  REAL:: cf
  REAL:: cds,cfs
  REAL:: h12,h12s
  REAL:: hz
  REAL:: psi
  REAL:: ret
  REAL:: usu
  REAL:: vdu
!----------------------------------------------------------------------------
      Z = D2*ABS(H) 
      RET = WRE*WRE*U*D2 
      CALL CDCF(RET,H,CD,CF,H12)   ! returns cd,cf,h12
      IF (MA) 1,1,2                ! ma is 'always' zero
    1 V = 0. 
      GO TO 13

    2 GO TO(13,13,13,6,6,5,5,3),MA     ! ma=suction mode
    3 IF (H) 4,1,1 
    4 V = 25.*(H+1.98)*D2*US 
      GO TO 12 
    5 IF(H)6,1,1 
    6 B = BB(MA) 
      PSI = AA(MA) 
      IF (B.NE.0.0) PSI=B*LOG(RET)+PSI 
      IF (PSI.GE.1.52  .AND.  PSI.LE.1.99) GO TO 11 
      B=0. 
      IF (PSI.LT.1.52) PSI=1.52 
      IF (PSI.GT.1.99) PSI=1.99 
   11 HZ = SIGN(PSI,H) 
      CALL CDCF(RET,HZ,CDS,CFS,H12S) 
      V=(U*(CDS-(PSI+B)*CFS)+D2*US*(B-PSI+H12*(B+PSI)))/(B+ABS(H)-1.) 
   12 IF(V)13,1,1 

   13 VDU =V/U
      USU =US/U 
      DD2 = (CF-(2.+H12)*USU*D2+VDU)*DL      ! RHS of (42)
      DZ  = (CD - 3.*Z*USU + VDU)*DL         ! RHS of (43)
   14 RETURN 
END Subroutine Grup   ! -----------------------------------------------------
                                                                        
!+
SUBROUTINE H12B(H,H12,EPST)
! ---------------------------------------------------------------------------
! PURPOSE - Compute H12 from H32
! called by CDCF GRP GRS
IMPLICIT NONE
  REAL,INTENT(IN):: h
  REAL,INTENT(OUT):: h12
  REAL,INTENT(OUT):: epst
!----------------------------------------------------------------------------
  IF (h==0.0) RETURN   ! no changes to h12 or epst

  IF (h < 0.0) THEN
    H12 = (H-1.36364)/(H*4.36364 + 5.36364) ! don't set epst (?)
  ELSE
    IF (h < 1.57258) THEN
      H12 =(SQRT(H-1.515090))*((-227.18220*H+724.55916)*H-583.60182)    &
        +4.0292200                                                        
      EPST = ((-.03172850655*H12+.3915405523)*H12-1.686094798)*H12      &
        +2.512588652                                                      
    ELSE
      H12 = (25.71578574*H-89.58214201)*H + 79.87084472 
      EPST = (H*2.221687229 - 4.226252829)*H+1.372390703
    END IF
  END IF

  RETURN
END Subroutine H12B   ! -----------------------------------------------------

!+
SUBROUTINE MOMENT (X,Y,NQ, xda,yda, deflg, momag)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the potential flow pitching moment coefficient about the
! quarter chord point and the section hinge coefficient about the flap hinge
! point.
! called by ProcessCommands

  REAL,INTENT(IN),DIMENSION(:):: x,y
  INTEGER,INTENT(IN):: nq
  REAL,INTENT(IN):: xda,yda
  REAL,INTENT(IN):: deflg
  INTEGER,INTENT(IN):: momag   ! a print code. set by ALFA card.

!      COMMON/PRAL/DLT,DLTU,ALN,ALV(14),NAL,ITP,CML(14),CRL(14),         &
!     &DARF,ITIT1,ITIT2,D1I(121,5),AVI(5,5),BVI(5,5),                    &
!     & FMACH,BEMACH,IVI,NAIT(5),XDA,YDA,DEFLG,MOMAG,NAI                 

  CHARACTER(LEN=*),PARAMETER:: FMT1 = &
    '("AIRFOIL ",A,"  HINGE POINT AT  X/C=",F6.3, "   Y/C=", F6.3)'
  CHARACTER(LEN=*),PARAMETER:: FMT2 = &
    '(" ALPHA      CM         CH    DELTA=",F6.2," DEG.")'
  CHARACTER(LEN=*),PARAMETER:: FMT32 = "(F6.2,2F11.6)"

  REAL:: cm      ! not the global variable
  REAL:: cr
  REAL:: drq,drqv
  REAL:: dx,dy
  INTEGER:: i,k,m
  REAL:: rq,rqv
  REAL:: switch
  REAL:: vq,vqv
  REAL:: xd,yd
!----------------------------------------------------------------------------
  IF (MOMAG==1) THEN
    WRITE(IDRU,FMT1) airfoilID, XDA,YDA 
    WRITE(IDRU,FMT2) DEFLG
  END IF

  DO 31 M=1,nal
    XD=0.25
    YD=0.0
    DO 20 K=1,2 
      SWITCH=1. 
      DRQV=-1. 
      DX=X(1)-XD                          ! x,y,vpr are global
      DY=Y(1)-YD 
      RQV=DX*DX+DY*DY 
      VQV=VPR(1,M)**2 
      CR=0. 
      DO 10 I=2,NQ 
        DX=X(I)-XD 
        DY=Y(I)-YD 
        RQ=DX*DX+DY*DY 
        DRQ=RQ-RQV 
        VQ=VPR(I,M)**2 
        IF(XD.EQ.0. .OR. k==1)GO TO 5 
        IF(DRQV.GE.0. .OR. DRQ.LT.0.)GO TO 5 
        CR=CR+SWITCH*RQV*VQV 
        SWITCH=-SWITCH 
    5   IF (SWITCH.GT.0.) CR=CR+RQ*VQV-RQV*VQ 
        RQV=RQ 
        DRQV=DRQ 
        VQV=VQ
   10 END DO
      CR=.25*CR 
      IF (K.EQ.1) CM=CR 
      XD=XDA 
      YD=YDA
   20 END DO

      CML(M)=CM                     !   cml,cmr are global
      crl(m)=cr           ! crl is set but never used
      IF (MOMAG==1) WRITE(IDRU,FMT32) ALV(M),CM,CR 
 31 END DO

  RETURN 
END Subroutine Moment   ! ---------------------------------------------------                                          

!+                                                                        
SUBROUTINE PADD(X,Y,BETA,NQ,MEI,KEI,XST)
! ---------------------------------------------------------------------------
! PURPOSE - Insert additional points along the airfoil surface using a spline.
! called by Flap and ProcessCommands

  REAL,INTENT(OUT),DIMENSION(:):: x,y,beta
  INTEGER,INTENT(IN OUT):: nq
  INTEGER,INTENT(IN):: mei
  INTEGER,INTENT(IN):: kei
  REAL,INTENT(IN):: xst
!      DIMENSION X(65),Y(65),BETA(65) 
!      COMMON/PLTM/MPL,MGC,XZEH,YZEH,MSPLI

  REAL:: dlta
  REAL:: dx,dy
  REAL:: fkp
  REAL:: ga,gas
  INTEGER:: i,k,ke,kr,ks
  INTEGER:: m,me
  REAL:: ph1,ph2
  REAL:: tg,tgs
  REAL:: xi
  REAL:: xm,xmp
  REAL:: x0
  REAL:: xv
!----------------------------------------------------------------------------
      IF(MSPLI.EQ.0)CALL SPLITZ(X,Y,NQ,BETA) 
      M=MEI 
      K=KEI 
      IF(M.NE.0)GO TO 14 
      XV=X(1) 
      DO 10 I=2,NQ 
      XM=X(I) 
      DX=XM-XV 
      IF(K.GT.0.AND.DX.LT.0.)GO TO 10 
      IF((XM-XST)*DX.GE.0.)GO TO 12 
   10 XV=XM

      I=NQ 
   12 M=I-1 
      K=1 
   14 X0=X(1) 
      XM=X(M) 
      XMP=X(M+1) 
      DX=XMP-XM 
      PH1=ATAN2(SQRT(XM*(X0-XM)),XM-.5*X0) 
      PH2=ATAN2(SQRT(XMP*(X0-XMP)),XMP-.5*X0) 
      DY=Y(M+1)-Y(M) 
      DLTA=ATAN2(DX,-DY) 
      GA=BETA(M)-DLTA 
      GAS=DLTA-BETA(M+1) 
      TG=SIN(GA)/COS(GA) 
      TGS=SIN(GAS)/COS(GAS) 
      KS=NQ 
    2 KR=KS+K 
      X(KR)=X(KS) 
      Y(KR)=Y(KS) 
      BETA(KR)=BETA(KS) 
      KS=KS-1 
      IF(KS.GT.M)GO TO 2 
      FKP=REAL(K+1) 

      DO 4 KE=1,K
      XI=REAL(KE)/FKP 
      IF (XM.GE.XST .OR. XMP.GE.XST) &
        XI=(.5*X0*(1.+COS(PH1+XI*(PH2-PH1)))-XM)/DX                                                         
      IF(MEI.EQ.0)XI=(XST-XM)/(XMP-XM) 
      ETA=XI*(1.-XI)*(TG*(1.-XI)+TGS*XI) 
      ME=M+KE 
      X(ME)=X(M)+XI*DX-ETA*DY 
      Y(ME)=Y(M)+ETA*DX+XI*DY 
      BETA(ME)=ATAN(TG*(1.-XI)*(1.-3.*XI)+TGS*XI*(2.-3.*XI))+DLTA
    4 END DO

      NQ=NQ+K
  RETURN 
END Subroutine PADD   ! -----------------------------------------------------



!+
SUBROUTINE PANEL(NK,A,GAMMA,CAE)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the flow characteristics about a given airfoil
! called by ProcessCommands

  INTEGER,INTENT(IN):: nk
  REAL,INTENT(OUT),DIMENSION(:,:):: a
  REAL,INTENT(OUT),DIMENSION(:,:):: gamma
  REAL,INTENT(OUT),DIMENSION(:):: cae

!      COMMON GAP1(506),AGAM(6),GAP2(242),DS(122),VF(121),ARG(121),      &
!     & GAP3(57),IGAP(2),NQ,IGAQ(3),CM,ETA,ABFA,PI,BOGEN                 

!      COMMON/PRAL/DLT,DLTU,ALN,ALV(14),NAL,ITP,CML(14),CRL(14),         &
!     &DARF,ITIT1,ITIT2,D1I(121,5),AVI(5,5),BVI(5,5),                    &
!     & FMACH,BEMACH,IVI,NAIT(5),XDA,YDA,DEFLG,MOMAG,NAI                 

!      COMMON/CHARS/NAMP(12),CPV(2),ALTX(4,2),NNESE,NUFF(80)

!      DIMENSION GAMMA(NK,2),A(NK,NK)
!      REAL BETA(121)
!       REAL DELTA(120),G(121,2),     &
!           ZZL(120),DA(4),DSP(121)    ! X(121),Y(121)

  REAL,DIMENSION(121):: delta
  EQUIVALENCE(delta,p1)
  REAL,DIMENSION(121,2):: g
!  EQUIVALENCE(g(1,1),vf(1))
  REAL,DIMENSION(121):: zzl
  EQUIVALENCE(zzl,p)

  REAL,DIMENSION(4):: da

  REAL:: bem
  REAL:: bw
  REAL:: ca
  REAL:: cbm,sbm
  REAL:: cosab,sinab,cosaq,sinaq
  REAL:: cp
  REAL:: dah,das,dat,dau,dax,day,daz
  REAL:: den,den1
  REAL:: dfx,dfy
  REAL:: dlr
  REAL:: du,dv
  REAL:: duz,dvz
  REAL:: dwu,dwv
  REAL:: dx,dy,dl
  REAL:: edge1,edge2
  REAL:: emax
  REAL:: eps
  REAL:: fab
  REAL,PARAMETER:: FED=0.00001
  REAL:: flom
  REAL:: fm
  REAL:: fnmx,fnmy
  REAL:: fnmxs,fnmys
  REAL:: fp
  REAL:: ga,gas
  REAL:: gr
  INTEGER:: k
  INTEGER:: l
  INTEGER:: m,n
  INTEGER:: naga
  INTEGER:: ne
  INTEGER:: nkor
  REAL:: p,p1
  REAL:: pk,pkl
  REAL:: rab
  REAL:: sp
  REAL:: tad,tal
  REAL:: tg,tgs
  REAL:: u,v
  REAL:: ulm
  REAL:: uq,vq
  REAL:: w
  LOGICAL:: WAKE
  REAL:: wu,wv
  REAL:: wus,wvs
  REAL:: xaa,yaa
  REAL:: xhm,yhm
  REAL:: xi
  REAL:: xmre,ymre
  REAL:: xn,yn
  REAL:: xnp,ynp
  REAL:: xnr,ynr
  REAL:: xz,yz
!----------------------------------------------------------------------------
  WRITE(3,*) "Entering Panel, nk=", nk, "   nq=", nq
  g(:,1)=vf(:)
  g(:,2)=arg(:)
  itp=2
      NAGA=INT(AGAM(10))
!      MMA=0 
      DELTA(1)=0.5*(xp(1)+xp(NQ))-PI 
      XHM=X(NQ) 
      YHM=Y(NQ)

      DO 2 N=1,NQ 
        XN=X(N) 
        YN=Y(N) 
        DX=XN-XHM 
        DY=YN-YHM 
        DL=SQRT(DX*DX+DY*DY) 
        DS(N)=DL 
        IF(N.NE.1)GO TO 22 
        WAKE=DL.GT..0001*XN
        IF(.NOT.WAKE)GO TO 24 
   22   DELTA(N)=ATAN2(DX,-DY) 
   24   XHM=XN 
        YHM=YN
        yp(n)=delta(n)
    2 END DO

    write(3,*) 'x'
    write(3,'(5es15.5)') x(1:nq)
    write(3,*) 'y'
    write(3,'(5es15.5)') y(1:nq)
    write(3,*) 'delta'
    write(3,'(5es15.5)') delta(1:nq)
    write(3,*) 'dd'
    write(3,'(5es15.5)') ds(1:nq)



      IF(DELTA(1).GT.0.)DELTA(1)=DELTA(1)-2.*PI
      EDGE1=xp(1)-DELTA(1) 
      EDGE2=xp(NQ)-DELTA(1)-2.*PI 
      COSAB=COS(EDGE1) 
      SINAB=SIN(EDGE1) 
      COSAQ=COS(EDGE2) 
      SINAQ=SIN(EDGE2) 
      ZZL(1)=PI*.5*(SINAB-SINAQ)+COSAB+COSAQ 



      DO 20 M=1,NQ                                 ! BIG loop
      IF (.NOT.WAKE.AND.M.EQ.NQ) CYCLE   !!! GO TO 20 
      BEM=xp(M) 
      XMRE=X(M) 
      YMRE=Y(M) 
      IF(M.GT.1.AND.M.LT.NQ)GO TO 300 
      XMRE=.5*(X(1)+X(NQ)) 
      YMRE=.5*(Y(1)+Y(NQ)) 
      BEM=DELTA(1)+REAL(M/NK)*.5*PI            ! only cuts in when m >= nk


  300 CBM=SIN(BEM) 
      SBM=-COS(BEM) 
      IF(M.LE.NK)A(M,1)=0. 
      G(M,1)=-CBM*6.283185 
      G(M,2)=-SBM*6.283185 
      XNR=XN 
      YNR=YN 
      DO 9 NKOR=1,NQ               ! inner BIG loop
      XNP=X(NKOR) 
      YNP=Y(NKOR) 
      DAT=0. 
      DAS=0. 
      DAY=0. 
      IF(.NOT.WAKE.AND.(NKOR.EQ.1.OR.NKOR.EQ.M+NK))GO TO 8 
      IF (NKOR.EQ.1 .AND. M.EQ.NQ) GO TO 8
      N=NKOR-1 
      IF (N.EQ.0) N=1 
      DU=XNP-XNR 
      DV=YNP-YNR 
      DL=DS(NKOR) 
      CP=DU/DL 
      SP=DV/DL 
      GA=xp(N)-DELTA(NKOR) 
      GAS=DELTA(NKOR)-xp(NKOR) 
      XAA=XNR 
      YAA=YNR 
      DLR=DL 
      IF(M.NE.NKOR.AND.M.NE.NKOR-1)GO TO 200 
      IF(M.NE.1.AND.M.NE.NQ)GO TO 4 
      IF(.NOT.WAKE)GO TO 6 
      IF(NKOR.EQ.1)GO TO 7 
  200 DX=XMRE-XNR 
      DY=YMRE-YNR 
      U=(DX*CP+DY*SP)/DL 
      V=(DY*CP-DX*SP)/DL 
      PK=1. 
      PKL=1. 
      IF(NKOR.EQ.1)GO TO 204 
      EMAX=ABS(GA+GAS)/8. 
      RAB=ABS(V) 
      IF(U.LT.0..OR.U.GT.1.)RAB=MAX(RAB,ABS(U),ABS(U-1.)) 
      FAB=EMAX/(FED*RAB) 
      IF(FAB.LT.1.)GO TO 204 
      PK = REAL(INT(FAB**.333333)+1) 
      TG=SIN(GA)/COS(GA) 
      TGS=SIN(GAS)/COS(GAS) 
  202 XI=PKL/PK 
      ETA=XI*(1.-XI)*(TG*(1.-XI)+TGS*XI) 
      XZ=XAA+XI*DU-ETA*DV 
      YZ=YAA+ETA*DU+XI*DV 
      DUZ=XZ-XNR 
      DVZ=YZ-YNR 
      DL =SQRT(DUZ*DUZ+DVZ*DVZ) 
      CP=DUZ/DL 
      SP=DVZ/DL 
      DX=XMRE-XNR 
      DY=YMRE-YNR 
      U=(DX*CP+DY*SP)/DL 
      V=(DY*CP-DX*SP)/DL 
  204 VQ=V*V 
      UQ=U*U 
      DEN1=U-.333333 
      DEN = DEN1*DEN1+VQ 
      IF(DEN.LT.10000.)GO TO 205 
      WU=.5*V/DEN 
      WV=-.5*DEN1/DEN 
      DEN1=U-.6666667 
      DEN=DEN1*DEN1+VQ 
      WUS=.5*V/DEN 
      WVS=-.5*DEN1/DEN 
      DWU=0. 
      DWV=0. 
      GO TO 206 
  205 ULM=U-1. 
      TAD=ATAN2(V,U*ULM+VQ) 
      TAL=.5*LOG((VQ+UQ)/(VQ+ULM*ULM)) 
      WU=(1.-U)*TAD+V*TAL 
      WV=V*TAD-(1.-U)*TAL-1. 
      WUS=U*TAD-V*TAL 
      WVS=1.-U*TAL-V*TAD 
      DWU=(U-UQ+VQ)*TAD+V*(2.*U-1.)*TAL-V 
      DWV=(UQ-VQ-U)*TAL+(2.*U-1.)*(V*TAD-.5) 
  206 FNMX=WU*CP-WV*SP 
      FNMY=WV*CP+WU*SP 
      FNMXS=WUS*CP-WVS*SP 
      FNMYS=WVS*CP+WUS*SP 
      DFX=DWU*CP-DWV*SP 
      DFY=DWV*CP+DWU*SP 
      DAX=(FNMX*CBM+FNMY*SBM)/PK 
      DAZ=(FNMXS*CBM+FNMYS*SBM)/PK 
      DAU=DFX*CBM+DFY*SBM 
      IF (NKOR==1) THEN
        DAX=DAX*COSAB+(FNMX*SBM-FNMY*CBM)*SINAB 
        DAZ=-DAZ*COSAQ-(FNMXS*SBM-FNMYS*CBM)*SINAQ
      END IF
      DAT=DAT+DAX*(PK+1.-PKL)+DAZ*(PK-PKL) 
      DAS=DAS+DAX*(PKL-1.)+DAZ*PKL 
      DAY=DAY+(DAU/PK+DAX*(PKL-1.)*(PK+1.-PKL)+DAZ*PKL*(PK-PKL))/PK 
      IF(PKL.GT.PK-.01)GO TO 8 
      PKL=PKL+1. 
      XNR=XZ 
      YNR=YZ 
      GO TO 202 
    4 FLOM=REAL(NKOR-M) 
      DAS=.5*GAS+(FLOM-.3333333)*(GA-GAS) 
      DAT=.5*GA+(FLOM-.6666667)*(GA-GAS)+FLOM*PI 
      DAY=(GA+GAS)/12. 
      GO TO 8 
    6 DAT=LOG(DS(NQ)/DL) 
      DAS=-1. 
      GO TO 8 
    7 DAS=PI*.5*(COSAB-COSAQ)-SINAB-SINAQ

    8 da(1:4)=0.0

      IF(NKOR.EQ.1)GO TO 106
      IF(NKOR.GE.3)FM=DS(NKOR-1)/DLR 
      IF(NKOR.LE.NK)FP=DS(NKOR+1)/DLR 
      IF(NKOR.EQ.2)GO TO 102 
      IF(NKOR.EQ.NQ)GO TO 104 
      DAH=DAY/2. 
      DA(1)=-DAH/(FM*(FM+1.)) 
      DA(2)=DAH*(1./FM-1./(FP+1.)) + DAT 
      DA(3)=DAH*(1./FP-1./(FM+1.)) + DAS 
      DA(4)=-DAH/(FP*(FP+1.)) 
      GO TO 108 
  102 DA(2)=DAT-DAY/(FP+1.) 
      DA(3)=DAS+DAY/FP 
      DA(4)=-DAY/(FP*(FP+1.)) 
      GO TO 108 
  104 DA(1)=-DAY/(FM*(FM+1.)) 
      DA(2)=DAT+DAY/FM 
      DA(3)=DAS-DAY/(FM+1.) 
      GO TO 108 
  106 DA(3)=DAS+DAT
  108 continue

     write(3,'(2i4,10f7.3)') m,nkor, da,fm,fp,dah,das,dat,day
      DO 116 L=1,4
      NE=NKOR+L-3 
      IF (NE.GT.NQ.OR.NE.LT.1) CYCLE   !!! GO TO 116 
      IF(NE.LT.NQ)GO TO 110 
      IF(M.LT.NQ)A(M,1)=A(M,1)-DA(L) 
      IF(M.EQ.NQ)ZZL(1)=ZZL(1)-DA(L) 
      CYCLE   !!! GO TO 116 
  110 IF(M.EQ.NQ)GO TO 112 
      IF(L.EQ.4)A(M,NE)=DA(L) 
      IF(L.LT.4)A(M,NE)=A(M,NE)+DA(L) 
      CYCLE   !!! GO TO 116 
  112 IF(L.EQ.4)ZZL(NE)=DA(L) 
      IF(L.LT.4)ZZL(NE)=ZZL(NE)+DA(L) 
  116 END DO 

   10 XNR=XNP 
    9 YNR=YNP 
   20 END DO     ! end of BIG loop
      write(3,*) "At end of 20 loop, nq=", nq, "   nk=",nk
           write(3,*) a(1:nq,1:nq)

      IF(.NOT.WAKE)A(1,NK)=A(1,NK)-1.
      IF(WAKE)GO TO 139 
      zzl(2:nk)=0.0
      ZZL(1)=1.
      W=(1.+DS(3)/DS(2))**.66667 
      ZZL(3)=.5/(W-1.) 
      ZZL(2)=-W*ZZL(3) 
      W=(1.+DS(NK)/DS(NQ))**.66667 
      ZZL(NK-1)=-.5/(W-1.) 
      ZZL(NK)=-W*ZZL(NK-1)

  139 DO 144 N=1,NK 
        DO 141 M=1,2
          GR=0.
          DO L=1,NK
            GR=GR+G(L,M)*A(L,N)               ! dot product
          END DO
          IF(WAKE)GR=GR+G(NQ,M)*ZZL(N)
          GAMMA(N,M)=GR
  141   END DO
        DO 143 M=N,NK
          BW=0.
          DO L=1,NK
            BW=BW+A(L,N)*A(L,M)               ! dot product
          END DO
          BW=BW+ZZL(M)*ZZL(N)
          DELTA(M)=BW
  143   END DO

        DO M=N,NK 
          A(M,N)=DELTA(M)
        END DO
  144 END DO

      DO N=1,NK
        DO M=N,NK 
          A(N,M)=A(M,N)
        END DO
      END DO

!  write(3,*) "Just before gauss"
!  write(3,*) a(1:nk,1:nk)
      eps=0.0
      CALL GAUSS(A,GAMMA,NK,2,EPS)

!      WRITE(3,*) "After Gauss"
!      WRITE(3,*) "gamma(:,1)=", gamma(1:nk,1)
!      WRITE(3,*) "gamma(:,2)=", gamma(1:nk,2)

      DO 152 M=1,2 
        CA=0. 
        DO K=2,NK 
          CA=CA+(DS(K)+DS(K+1))*GAMMA(K,M)
        END DO
        CAE(M)=CA/X(1)
  152 END DO

      ALN=ATAN2(CAE(1),CAE(2))/BOGEN 
      IF(NAGA.EQ.0)GO TO 48 
      WRITE(IDRU,42)airfoilID,CAE(:),ALN
   42 FORMAT(" PANEL METHOD AIRFOIL ",A," CL=",2F10.5,"   ALPHA0=",F5.2," DEG.")
      IF (NAGA.LT.3) GO TO 48
      WRITE(IDRU,*) &
        "  N       X         Y        BETA      GAMMA1     GAMMA2     DS"
      DO K=1,NK
        WRITE(IDRU,46)K-1,X(K),Y(K),xp(K),GAMMA(K,1),GAMMA(K,2),GA,GAS,DS(K+1)
     END DO
   46 FORMAT (I3,2F10.5,3F11.5,F9.5)
      WRITE(IDRU,46)NK,X(NQ),Y(NQ),xp(NQ),GAMMA(1,1),GAMMA(1,2)

  48  vf(1:nk)=gamma(1:nk,1)
      arg(1:nk)=gamma(1:nk,2)
      vf(nq)=-gamma(1,1)
      arg(nq)=-gamma(1,2)

  RETURN 
END Subroutine Panel   ! ----------------------------------------------------

!+
SUBROUTINE ProcessCommands()
! ---------------------------------------------------------------------------
IMPLICIT NONE
  CHARACTER(LEN=*),PARAMETER,DIMENSION(20):: MARKEN = (/  &
   'TRA1', 'TRA2', 'ALFA', 'AGAM', 'ABSZ',                &
   'STRK', 'ENDE', 'DIAG', 'RE  ', 'STRD',                &
   'FLZW', 'PLWA', 'PLW ', 'TRF ', 'APPR',                &
   'CDCL', 'PAN ', 'FXPR', 'FLAP', 'PUXY' /) 
  CHARACTER(LEN=*),PARAMETER,DIMENSION(2):: CPV= (/ "VELOCITY ","PRESSURE "/)
  CHARACTER(LEN=*),PARAMETER,DIMENSION(2):: ALTX = (/ &
    "ZERO-LIFT LINE", "CHORD LINE    " /)


  REAL:: alphaMax,alphaMin
  REAL,DIMENSION(5):: als
  REAL:: anri
  REAL:: arcl,arclu
  REAL,DIMENSION(5):: bant
  REAL:: batf
  REAL:: bf
  REAL:: bl1,bl2
  CHARACTER(LEN=80):: buffer
  REAL:: ca
  REAL:: cant
  real:: chord
  REAL:: cmdu
  REAL:: cwges
  REAL:: cwp
  REAL:: cwnt
  REAL:: cws
  REAL:: cwsf,cwsfu
  REAL:: d
  REAL:: defl
  REAL:: deflg
  REAL:: delfg
  REAL:: dgf
  REAL:: dichte = 0.12533

  REAL:: dltr,dltur
  REAL:: dst
  INTEGER:: errCode ! tests all I/O for errors
  REAL:: fau
  REAL:: ff
  REAL:: flch
  REAL:: fst
  REAL:: gau
  REAL:: gdf
  REAL:: gew
  REAL:: gltz
  REAL:: gst
  INTEGER:: i,j,k,m,n
  INTEGER:: ipu
  INTEGER:: istift
  INTEGER:: ivmax
  INTEGER:: ivs
  INTEGER:: jp
  INTEGER:: jr ! the number of Reynolds numbers
  INTEGER:: jt
  INTEGER:: jz
  INTEGER:: kei
  INTEGER,DIMENSION(5):: ma,mu
  INTEGER:: mei,meig
  INTEGER:: momag
  INTEGER:: mtr=0   ! counts angle entries in TRA1
  INTEGER,DIMENSION(5):: mur
  INTEGER:: mxz
  INTEGER:: nan
  INTEGER:: nf
  INTEGER:: nkr
!  INTEGER,DIMENSION(5):: nline
!  INTEGER,DIMENSION(5):: npar
  INTEGER:: nqrs
  INTEGER:: nt
  INTEGER:: nupa,nupe,nupi,nupu
  REAL:: pa
  REAL,DIMENSION(14):: puff  ! temporary loc. for 14 reals on most records
  real,dimension(14):: p,v   !don't really belong here
  REAL,DIMENSION(5):: re
  REAL,DIMENSION(5):: rer
  REAL:: rerx
  REAL:: rua
  REAL,DIMENSION(42):: t
  REAL:: thk,thkp
  REAL,DIMENSION(5):: tm
  REAL,DIMENSION(5):: tst
  REAL:: v1
  REAL:: valf
  REAL:: vkmh
  REAL:: vmax
  REAL:: vs
  REAL:: xda,yda
  REAL:: xstx
  REAL:: ybl
  REAL:: zaeh = 13.6E-6
!----------------------------------------------------------------------------
  DO
    READ(ILES,'(A)',IOSTAT=errCode) buffer
    WRITE(DBG,*) "Data:"//Trim(buffer)
    IF (errCode > 0) THEN
      WRITE(*,*) "IO error in ProcessCommands"
      STOP
    END IF
    IF (errCode < 0) THEN
      WRITE(IDRU,*) "End of input file without an ENDE record."
      EXIT
    END IF
    k=TestKeyWord(buffer(1:4), MARKEN)
    SELECT CASE(k)
      CASE(0)

      CASE(1)                                                          ! TRA1
        airfoilID=buffer(7:10)
        READ(buffer(11:80), '(14F5.2)', IOSTAT=errCode) puff
        IF (errCode /= 0) STOP "buffer parsing error in ProcessCommands(TRA1)"
        WRITE(DBG,*) "TRA1:"//airfoilID
        IF (MTR .EQ. 0) JST=0 
        i=0
        DO                                 ! load global ani and alfa arrays
          i=i+1
          ANRI=RUND(PUFF(I)*ABFA, 1000.) 
          IF (ANRI==0.) THEN
            IF (JST .NE. 0) EXIT
            JST=MTR+1    ! global 
          END IF
          MTR=MTR+1 
          ANI(MTR)=ANRI         ! 1,3,5,... of puff in ani
          I=I+1 
          ALFA(MTR)=PUFF(I)     ! 2,4,6,... of puff in alfa
          IF (I>=14) EXIT
        END DO
        JAB=MTR     ! jab is global

      CASE(2)                                                        !   TRA2
        READ(buffer(11:80),'(14F5.2)',IOSTAT=errCode) puff
        IF (errCode /= 0) STOP "buffer parsing error in ProcessCommands(TRA2)"
        DO I=1,13
          PURES(I)=RUND(PUFF(I),1000.)
        END DO
        MSPLI=0 
        ITP=0 
        IZZ=INT(PUFF(14)) 
  alfr(1:jab)=alfa(1:jab)   ! seems to be some confusion
        CALL TRAPRO()
!!!        IF (NUPA.NE.0) READ (efu1,'(A)') NAMP 
        XDA=0. 
        YDA=0. 
        DEFLG=0.  ! flap deflection
        FLCH=0.   ! flap chord   
        mtr=0     ! set up for another set of tra1/tra2 records

      CASE(3)                                                          ! ALFA
        READ(buffer(5:10),'(3I1,I3)') nupa,nupe,nupi,nupu
        READ(buffer(11:80),'(14F5.2)') puff(1:14)
        IF(NUPA/=0) THEN
          MOMAG=NUPA 
          AGAM(2)=REAL(NUPE) 
          IF (NUPA==1) AGAM(10)=REAL(NUPE+1)
        END IF

        IF (NUPU/=0) THEN
          alca(1:14)=puff(1:14)
          ITIT1=NUPI/2+1            !  1,1,2,2 ; nupi=0,1,2,3
          ITIT2=NUPI-2*ITIT1+3      !  1,2,1,2 ; nupi=0,1,2,3
          DARF=0. 
          IF (ITIT2.NE.1) DARF=1. 
          NAL=ABS(NUPU) 
          IF (NAL > 14) NAL=4
        END IF

        DO I=1,NAL 
          PA=ALCA(I) 
          IF (PA.LE.-99.) THEN
            PA=RS(30+I)    ! this needs work. rs is never set ??????
          ELSE
            PA=PA+DARF*DARG
          END IF
          ALV(I)=PA
        END DO

        CALL MOMENT(X,Y,NQ,XDA,YDA,DEFLG,MOMAG)
        p(1:nal)=alv(1:nal)-darf*darg   ! right spot?????

        IF (AGAM(2) .EQ. 0.0) CYCLE
        CALL DIA(X,Y,NQ,THK) 
        THKP=100.*THK 

!        WRITE(IDRU,*) ACHAR(12)
        WRITE(IDRU,'(A,F8.2,"% THICKNESS")') " AIRFOIL "//airfoilID, thkp
        SELECT CASE(ITP)
          CASE(1)
 !           WRITE(IDRU,337)NZT//"AIRFOIL "//airfoilID,THKP,(P(M),M=1,NAL)
 !         337 FORMAT(A,F8.2,"%",F11.2,13F7.2) 
          CASE(2)
            WRITE(IDRU,'(F9.2,"% FLAP",F9.2," DEGREES DEFLECTION")') &
              flch, deflg
          END SELECT
        WRITE(IDRU,*) CPV(itit1)//"-DISTRIBUTION"
        WRITE(IDRU,*) "ANGLES OF ATTACK RELATIVE TO THE "//ALTX(itit2)
        WRITE(IDRU,'(A,14F7.2)') "  N     X        Y   ",p(1:nal)

        DO N=1,NQ
          DO M=1,NAL
            V(m)=ABS(VPR(N,M))                 ! velocity
            IF (itit1==2) v(m)=1.0-v(m)*v(m)   ! pressure
          END DO
          WRITE(IDRU,'(I3,2F9.5,14F7.3)') n-1,x(n),y(n),v(1:nal)
        END DO

        IF (ITP.NE.2) THEN
          WRITE(IDRU,344) DARG,CM,ETA 
  344 FORMAT ("  ALPHA0 =",F5.2," DEG.   CM0 =",F7.4,"    ETA =",F6.3)
        END IF
      CASE(4)                                ! AGAM   (obsolete??? NO such!!)
      CASE(5)                                                          ! ABSZ
        READ(buffer(5:6),'(2I1)') nupa,nupe
        IF (nupa > 0) agam(3)=REAL(nupe)   ! referred to as mtr on p.44
        READ(buffer(16:20),'(F5.2)') puff(2)
        IF (puff(2) /= 0.0) abfa=puff(2)

      CASE(6)                                                          ! STRK
        WRITE(IDRU,*) "STRK currently disabled"
        READ(buffer(5:7),'(2I1)') nupa,nupe,nupi
        READ(buffer(8:10),'(I3)') nupu
        READ(buffer(11:80),'(14F5.2)') puff
        IF (nupu /= 0.0) THEN
          IF (nupu > 0) nt=0
          DO j=1,14
            IF (puff(j)==0.0) GO TO 100
            nt=nt+1
            t(nt)=10.0*puff(j)      ! never used
          END DO
          IF (ABS(nupu) > 14) CYCLE  ! go to next input record
        END IF
    100 continue
!        CALL STRDR(t,nt)
        IF (nupi /= 0) CYCLE   ! go to next input record
        DO i=1,nt
        ! CALL Straak(t(i),rua,ybl,mxz,istift)
        END DO

      CASE(7)                                                          ! ENDE
!        IF (mgc /= 0) CALL Gclose()
!        IF (mpl /= 0) CALL Finish()
        RETURN

      CASE(8)                                                          ! DIAG
        READ(buffer(7:7),'(I1)') nupi
        READ(buffer(8:10),'(I3)') nupu
!        iplot=1
        IF (nupi==1) THEN
          READ(buffer(11:20),'(2F5.2)') alphaMin,alphaMax
        ELSE
          alphaMin=0
          alphaMax=0
        END IF
        CALL Diagram(nupu, nupi, alphaMin,alphaMax) ! args changed - RLC

      CASE(9)                                                            ! RE
        READ(buffer(5:7),'(3I1)') nupa,nupe,nupi
        IF (nupa/=0) THEN
          agam(6)=REAL(nupe)  ! print code
          agam(8)=REAL(nupi)  ! plot code
        END IF
        READ(buffer(11:80),'(14F5.2)') puff

        IF (PUFF(2) /= 0.0) THEN
          DO J=1,5 
            RERX=PUFF(2*J) 
            IF (RERX.EQ.0.0) EXIT
            RE(J)=1.E5*RERX 
            IPU = INT(PUFF(2*J-1)) 
            MA(J) = IPU/100
            MU(J) = IPU/10 - 10*MA(J) 
            JR  = J     ! number of Reynolds numbers
          END DO
          xtri(1:4)=0.01*puff(11:14)
        END IF
        CALL GRP(NAL,RE,MU,JR,ISTIFT)  ! nal is set by ALFA
        MSPLI=0 
        JP=JR   ! jp set but never used

      CASE(10)                                                         ! STRD
        READ(buffer(8:10),'(I3)') mxz
        READ(buffer(11:20),'(2F5.2)') puff(1:2)
        IF (puff(1) /= 0.0) ybl=100.0*puff(1)
        IF (puff(2) /= 0.0) rua=100.0*puff(1)  ! rua,ybl,mxz used by Straak

      CASE(11)                                                         ! FLZW
        READ(buffer(5:7),'(3I1)') nupa,nupe,nupi
        IF (nupa==1) THEN
          agam(6)=REAL(nupe)
          agam(8)=REAL(nupi)
        END IF
        READ(buffer(8:10),'(I3)') nupu
        READ(buffer(11:80),'(14F5.2)') puff
        IF (puff(2) /= 0.0) THEN
          gdf=puff(1)
          vmax=puff(2)
          IF (puff(3) /= 0.0) dichte=0.1*puff(3)/9.806
          IF (puff(4) /= 0.0) zaeh=puff(4)*1E-6
          IF (puff(5)==0.0) CYCLE
          d=0.0
          DO 2499 j=1,5
            jz=2*j+3
            IF (puff(jz)==0.0) EXIT
            tm(j)=puff(jz)
            als(j)=puff(jz+1)
            mur(j)=nupu
            jt=j
2499          END DO
        END IF

     36 jp=jt
        WRITE(IDRU,*) "AIRCRAFT ORIENTED SUMMARY     AIRFOIL "//airfoilID
        WRITE(IDRU,*) "ANGLE OF ATTACK RELATIVE TO "//ALTX(itit2)
        !!! more gets added
        IVMAX=NINT(VMAX*3.6)
        WRITE(IDRU,38) GDF,IVMAX,DICHTE,ZAEH
        38 FORMAT (" W/S =",F6.2," KG/SQ.M    V MAX =",I4," KM/H  RHO =", &
          F6.3," KG*S E2/M E4   NU =",F11.8," SQ.M/S")
        WRITE(IDRU,40) (TM(J),ALS(J),J=1,JT)
        40 FORMAT (6X,5(5X,"C =",F5.2," THETA =",F5.2))
        V1 = SQRT(2.0*GDF/DICHTE)
        DO 48 I=1,NAL
          IVS = -I
          DO J = 1,JT
            IF (ALV(I)-ALS(J))42,44,42
            42 VALF = V1/SQRT(.11*X(1)*ABS(ALV(I)-ALS(J)))
            IF (VALF - VMAX)46,46,44
            44 VALF = VMAX
            46 RER(J) = VALF*TM(J)/ZAEH
          END DO  ! was 46
          CALL GRP(IVS,RER,MUR,JT,ISTIFT)
48        END DO  ! was loop 48
        MSPLI=0
!!!        IF (D /= 0.0) GO TO 72   ! you can't do this. In a do loop

      CASE(12)                                                         ! PLWA
        READ(buffer(8:10),'(I3)') nupu
        NAN = NUPU
        CWSFU=CWSF
        FAU = FF
        GAU = GEW
        DO 88 nf=1,nan
          FF = FF+PUFF(1)
          GEW = GEW+PUFF(2)
          CWSF=CWSF+0.001*PUFF(3)
          V1 = SQRT(2.*GEW/(DICHTE*FF))
          WRITE(IDRU,*) "AIRFOIL POLAR      AIRFOIL "//airfoilID
          CWS = CWSF/FF
          WRITE(IDRU,*) &
            "  B(M)  S(SQ.M) S*(SQ.M) W(KG)   W*(KG)  T/C  (T/C)*AP(SQ.M)  CDP"
          WRITE(IDRU,78)BF,FF,FST,GEW,GST,D,DST,CWSF,CWS
          78 FORMAT (3F8.2,2F8.0,4F8.4)
          WRITE(IDRU,*) &
            "  ALPHA     CL     CDP     CDT   V(KM/H)  VS(M/S)   L/D"
          DO I = 1,NAL
            CA = 0.
            CWP = 0.
            DO J = 1,JT
              BATF=BANT(J)*TST(J)/FST
              CALL VISC(I,J,CANT,CWNT,CMDU)
              CA=CA+CANT*BATF
              CWP=CWP+CWNT*BATF
            END DO
            CWGES = CWP + CWS + 1.03*CA*CA*FF/(PI*BF*BF)
            IF (ABS(CA).LT..01) CA=.01
            VKMH = 3.6*V1/SQRT(ABS(CA))
            VS = VKMH*CWGES/(3.6*ABS(CA))
            GLTZ = CA/CWGES
            WRITE(IDRU,86) ALV(I),CA,CWP,CWGES,VKMH,VS,GLTZ
          END DO
          86 FORMAT (F8.2,F8.3,2F8.4,F8.1,F9.3,F8.2)
          IF (nf /= 0) EXIT
     88 END DO
        IF (NF /= 0) THEN
          CWSF = CWSFU
          FF=FAU
          GEW = GAU
        END IF


      CASE(13)                                                          ! PLW
        READ(buffer(8:10),'(I3)') nupu
        READ(buffer(11:80),'(14F5.2)') puff(:)
        IF (puff(1) /= 0.0) THEN
          DST = PUFF(1)*0.01
          GST = PUFF(2)
          DGF = PUFF(3)
          CWSF = PUFF(4)*0.001
          DO J = 1,5
            JZ = 2*J+3
            IF (PUFF(JZ)==0.0) EXIT
            TST(J) = PUFF(JZ)
            BANT(J) = PUFF(JZ+1)
            MUR(J) = NUPU
            JT = J
          END DO
        END IF
        bf=0.0 
        fst=0.0
        nf=0.0
        CALL Dia(x,y,nq,d)
        IF (DST < 0.0) dst=d
        DO j=1,jt
          tm(j)=tst(j)*dst/d
          als(j)=0.0
          bf=bf+bant(j)
          fst=fst+bant(j)*tst(j)
        END DO
        ff=fst*dst/d
        gew=gst+(ff-fst)*dgf
        gdf=gew/ff
!!!        GO TO 36

      CASE(14)  ! TRF
        !!! nothing - old
      CASE(15)  ! APPR
        !!! nothing - old
      CASE(16)                                                         ! CDCL
        READ(buffer(5:5),'(I1)') nupa
        READ(buffer(8:10),'(I3)') nupu
        IF (nupa==0) THEN
!         CALL CDCL(nupu, jp, istift)
        ELSE
          READ(buffer(11:80),'(14F5.2)') puff
!          bl1=puff(1)+0.005
!          bl2=puff(2)+0.005
          DO k=1,5
!            nline(k)=INT(BL1*(10.**(K-3)))-10*INT(BL1*(10.**(K-4))) 
!            NPAR(K) =INT(BL2*(10.**(K-3)))-10*INT(BL2*(10.**(K-4))) 
          END DO
        END IF

      CASE(17:18)                                     ! PAN (17) or FXPR (18)
        READ(buffer(5:6),'(2I1)') nupa,nupe
        
        READ(buffer(11:80),'(14F5.2)') puff(1:14)
        IF (k==18) THEN
          READ(buffer(8:10),'(I3)') itp        ! only FXPR reads data
          CALL FixLes()
          mspli=0
        END IF
        IF (mspli==0) CALL SPLITZ(X,Y,NQ,XP) 
        IF (NUPA.NE.0) AGAM(10)=REAL(NUPE) 
        IF (NUPA.EQ.9) CYCLE 
        DO I=1,14                          ! was DO 172
          IF (PUFF(I).EQ.0.) CYCLE
          MEIG=INT(PUFF(I)) 
          MEI=MEIG/10 
          KEI =MEIG-10*MEI 
          XSTX=ABS(PUFF(I)-REAL(MEIG))   ! fraction of chord
          CALL PADD(X,Y,XP,NQ,MEI,KEI,XSTX) 
        END DO 

        xf(1:nq)=x(1:nq)      ! was DO 174
        yf(1:nq)=y(1:nq)
        betaf(1:nq)=xp(1:nq)

        DLTR=DLT
        DLTUR=DLTU 
        XDA=0. 
        YDA=0. 
        FLCH=0. 
        DEFLG=0. 
        NQRS=NQ 
        NKR=NQ-1 
        CALL PANEL(NKR,AMAT,GAMMA,CAE) 
        DARG = ALN 

      CASE(19)                                                         ! FLAP
        READ(buffer(11:30),'(4F5.2)') puff(1:4)
        chord=xf(1)
        flch=puff(1)
        xda=(1.0 - 0.01*flch)*chord
        yda=0.01*puff(2)*chord
        arcl=0.01*puff(3)*chord
        delfg=puff(4)
        dlt=dltr+delfg
        dltu=dltur-deflg
        defl=deflg*BOGEN
        arclu=0.01*puff(5)*chord
        mspli=1
        CALL Flap(xf,yf,betaf,nqrs,xda,yda,arcl,arclu,defl,x,y,xp,nq)
!        nkr=nq-1
        CALL Panel(nq-1,amat,gamma,cae)
        darg=aln

      CASE(20)                                                         ! PUXY
        CALL PuDeck()
    END SELECT
  END DO

  RETURN
END Subroutine ProcessCommands   ! ------------------------------------------

!+
SUBROUTINE PUDECK()
! ---------------------------------------------------------------------------
! PURPOSE - Write the computed airfoil coordinates to a file.
!   The original function was to "punch a deck" of cards.
! called by ProcessCommands

  INTEGER:: errCode
  INTEGER,DIMENSION(1):: iloc
!----------------------------------------------------------------------------

  OPEN(UNIT=ISTA,FILE="puxy.dat",STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
  IF (errCode /=0) THEN
    WRITE(*,*)    "Unable to open puxy.dat"
    WRITE(IDRU,*) "Unable to open puxy.dat"
    WRITE(DBG,*)  "Unable to open puxy.dat"
    RETURN
  END IF

  WRITE(ISTA,*) "OUTPUT FROM SUBROUTINE PUDECK"
  WRITE(ISTA,*) airfoilID
  iloc=MINLOC(x)

  WRITE(ISTA,'(2I5)') iloc(1),nq-iloc(1)+1
  WRITE(ISTA,'(5F12.7)' ) x(iloc(1):1:-1)     ! x-upper
  WRITE(ISTA,'(5F12.7)' ) y(iloc(1):1:-1)     ! y-upper
  WRITE(ISTA,'(5F12.7)' ) x(iloc(1):nq)       ! x-lower
  WRITE(ISTA,'(5F12.7)' ) y(iloc(1):nq)       ! y-lower
  CLOSE(UNIT=ISTA)
  RETURN 
END Subroutine PUDECK   ! ---------------------------------------------------

!+
SUBROUTINE QIP(X,Y,A)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the coefficients of a parabola thru 3 points
! called by DIA

  REAL,INTENT(IN),DIMENSION(3):: x,y
  REAL,INTENT(OUT),DIMENSION(3):: a
  REAL:: c1
!----------------------------------------------------------------------------
  C1 = (Y(2)-Y(1))/(X(2)-X(1)) 
  A(3)=(Y(3)-Y(1)-C1*(X(3)-X(1)))/((X(3)-X(1))*(X(2)-X(1))) 
  A(1)=Y(1)-C1*X(1)+A(3)*X(1)*X(2) 
  A(2)=C1-A(3)*(X(1)+X(2)) 
  RETURN 
END Subroutine QIP                                          

!+
FUNCTION Rund(A,B) RESULT(c)
! ---------------------------------------------------------------------------
! PURPOSE - A rounding function. Example: Rund(4376.2,0.01) returns 4400.
  REAL,INTENT(IN):: a,b
  REAL:: c
! ---------------------------------------------------------------------------
  c = (AINT(A*B+SIGN(.5,A)))/B 
  RETURN 
END Function Rund   ! -------------------------------------------------------                                          

!+
FUNCTION SING(A) RESULT(f)                  ! sine of angle in degrees
! ---------------------------------------------------------------------------
  REAL,INTENT(IN):: a     ! angle in DEGREES
  REAL:: f
! ---------------------------------------------------------------------------
  f = SIN(A*BOGEN) 
  RETURN 
END Function SinG   ! -------------------------------------------------------                                          
                                                                        
!+
SUBROUTINE SPLITZ(X,Y,NQ,BETA) 
! ---------------------------------------------------------------------------
! PURPOSE - Perform a spline fit of the nq points contained in arrays x and y.
!   The result is the array of angles BETA normal to the curve at every point
! called by PADD and ProcessCommands

  REAL,INTENT(IN),DIMENSION(:):: x,y
  INTEGER,INTENT(IN):: nq
  REAL,INTENT(OUT),DIMENSION(:):: beta


!      COMMON/PLTM/MGC,MPL,XZEH,YZEH,MSPLI

  REAL:: an
  REAL:: av
  REAL:: dbm
  REAL,PARAMETER:: DBSOLL = 1.E-6
  REAL:: dd,de
  REAL:: dx,dy
  REAL:: dxv,dyv
  REAL:: e,es
  REAL:: gp,gpq,gps,gsq
  REAL:: gn,gs
  INTEGER:: m,n,nk
  INTEGER:: nit
  REAL,DIMENSION(nq):: u,v,ds
  REAL:: un,vn
  REAL:: xn,yn
  REAL:: xv,yv
!----------------------------------------------------------------------------
      MSPLI=1 
      NK=NQ-1 

      DO 10 N=1,NQ
      XN=X(N) 
      YN=Y(N) 
      IF(N.EQ.1)GO TO 8 
      DX=XN-XV 
      DY=YN-YV 
      M=N-1 
      AN=SQRT(DX*DX+DY*DY) 
      DS(M)=AN 
      IF(M.EQ.1)GO TO 6 
      U(M)=DXV*DX+DYV*DY 
      V(M)=DXV*DY-DYV*DX 
      DD=-(AV/(AV+AN))*ATAN2(V(M),U(M)) 
      BETA(M)=SIN(DD)/COS(DD) 
    6 AV=AN 
      DXV=DX 
      DYV=DY 
    8 XV=XN 
      YV=YN
   10 END DO

      BETA(NQ)=0.
      NIT=0 
      DBM=0. 
   12 NIT=NIT+1 
      DBM=0. 
      AV=DS(1) 
      GS=BETA(2) 
      GN=.5*GS 

      DO 20 N=2,NK
      UN=U(N) 
      VN=V(N) 
      GP=(UN*GS+VN)/(VN*GS-UN) 
      AN=DS(N) 
      GPS=BETA(N+1) 
      IF(N.EQ.NK)GPS=.5*GP 
      GSQ=GS*GS 
      GPQ=GP*GP 
      E=AV*(2.*GP-GPS)*(1.+GPQ)-AN*(2.*GS-GN)*(1.+GSQ) 
      ES=AN+AV 
      IF(GSQ+GPQ.LT..09)GO TO 18 
      ES=AV*(1.+3.*GPQ-GP*GPS)*(UN-GP*VN)/(UN-GS*VN)+AN*(1.+3.*GSQ      &
     & -GS*GN)                                                          
   18 DE= .5*E/ES 
      IF(DBM.LT.ABS(DE))DBM=ABS(DE) 
      GS=GS+DE 
      BETA(N)=GS 
      AV=AN 
      GN=(UN*GS+VN)/(VN*GS-UN) 
      GS=GPS 
   20 END DO


      IF(DBM.GT.DBSOLL.AND.NIT.LE.15)GO TO 12
      DO 30 N=1,NQ 
      XN=X(N) 
      YN=Y(N) 
      IF(N.EQ.1)GO TO 28 
      DX=XN-XV 
      DY=YN-YV 
      DLT=ATAN2(DX,-DY) 
      IF(N.EQ.2)BETA(1)=DLT+ATAN(.5*BETA(2)) 
      BETA(N)=DLT-ATAN(BETA(N)) 
   28 XV=XN 
      YV=YN
   30 END DO

      BETA(NQ)=DLT-ATAN(.5*GN)
  RETURN 
END Subroutine Splitz   ! ---------------------------------------------------

!+
SUBROUTINE Straak()
  RETURN
END SUBROUTINE Straak

SUBROUTINE Strdr()
  RETURN
END SUBROUTINE Strdr



!+
FUNCTION TestKeyWord(test,dict) RESULT(k)
! ---------------------------------------------------------------------------
! PURPOSE - Find the position of a word in a dictionary. Return 0 if the
!   word is not in dictionary. Convert test word to uppercase (locally).
! called by ProcessCommands

  CHARACTER(LEN=*),INTENT(IN):: test
  CHARACTER(LEN=*),INTENT(IN),DIMENSION(:):: dict
  INTEGER:: k

  CHARACTER(LEN=LEN(dict(1))):: testWord
  INTEGER:: j
!----------------------------------------------------------------------------
  k=0
  testWord=test         ! makes all the tests between strings of equal length
  CALL UpCase(testWord) ! assumes dict entries are all upper case

  DO j=1,SIZE(dict)
    IF (testWord == dict(j)) THEN
      k=j
      RETURN
    END IF
  END DO
  RETURN
END Function TestKeyWord   ! --------------------------------------------------

!+
SUBROUTINE TRAPRO()
! ---------------------------------------------------------------------------
! PURPOSE -
! called by ProcessCommands

  REAL,DIMENSION(4):: a
  REAL:: aa
  REAL:: abgr
  REAL:: absz=0.0
  REAL:: abzt
  REAL,DIMENSION(4,3):: ac
  REAL:: ak,akp
  REAL:: alis,alisp
  REAL:: aliv
  REAL:: anu
  REAL:: argn
  REAL:: ari
  REAL:: at
  REAL:: b
  REAL:: b2
  REAL:: bi
  REAL:: cosai
  REAL:: csaip,snaip
  REAL:: csli,cslip
  REAL:: csp
  REAL,DIMENSION(3):: d
  REAL:: dal,dald
  REAL:: dalv
!  REAL:: deflg   ! never used
  REAL,DIMENSION(2):: drak,dram
  REAL:: dras
  REAL:: ff1,ff2
  REAL:: f
  REAL:: fg1,fg3
  REAL:: fii,fiip
  REAL,DIMENSION(3):: fint
  REAL,DIMENSION(30):: fkern
  REAL,DIMENSION(2):: fla,fls
  REAL:: fm,fmit
  REAL:: fni
  REAL:: fp,fv
  REAL:: g
  REAL:: habgr
  REAL,DIMENSION(2):: hk
  REAL:: hks,hkst,hksv
  INTEGER:: i,j
  INTEGER:: ib
  INTEGER:: ieppl
!  INTEGER:: iri
  INTEGER:: itmod
  INTEGER:: itmr
!  INTEGER:: ivi
  INTEGER:: jh,jn
  INTEGER:: l
  INTEGER:: m
  INTEGER:: magam
  INTEGER:: mer
  INTEGER:: mit
  INTEGER:: mm,mn,mq
!  INTEGER:: momag
  INTEGER:: n
!  INTEGER:: nai,nait
  INTEGER:: nhkw
  INTEGER:: nkr
  INTEGER:: nu
  REAL:: pb
  REAL:: pdif
  REAL:: phim
  REAL:: phish
  REAL:: ps
  REAL:: q
  REAL,DIMENSION(3):: r
  REAL:: rq,rqv
  REAL:: ruf
  REAL:: shks
!  REAL:: ski
  REAL:: sinai
  REAL:: sq
  REAL:: stref
  REAL:: sx,sy
  REAL:: sxi
  REAL:: tau
  REAL:: v1
  REAL:: vi
  REAL:: wc
  REAL:: wli
  REAL,DIMENSION(2):: wci,wsi
  REAL:: whk
  REAL:: wi
  REAL:: wiln
  REAL:: wl
  REAL:: wq
  REAL:: ws
  REAL:: wstr
  REAL:: wv
  REAL:: x1
!  REAL:: xda,yda
  REAL:: xnas,ynas
  REAL:: xpk,ypk
  REAL:: xr
!  REAL:: xws
  REAL:: zl,zp
!----------------------------------------------------------------------------
  WRITE(DBG,*) "Entering Tranpo, jab=",jab
!!!      CALL WANDEL(NUPRO,NAMP,12,5)

      ALFR(JAB+1)=0.
      ABZT=ANI(JAB) 
      IF(ABS(ABZT-ABSZ).LT..1) GO TO 14 
      IB=INT(.25*ABZT+.1) 
      MQ=2*IB 
      NKR=2*MQ 
      ABSZ=REAL(NKR) 
      ABGR=360./ABSZ 
      HABGR=.5*ABGR 

      DO M=1,IB
        ARI=REAL(MQ+1-2*M)*HABGR 
        fkern(m)=abgr/(PI*TAN(BOGEN*ari))
!        FKERN(M)=ABGR*COSG(ARI)/(SING(ARI)*PI)
      END DO

   14 MAGAM=INT(AGAM(2))
      NQ=NKR+1 
      IF(MAGAM.EQ.0) GO TO 22 
   22 alfa(1:29)=alfr(1:29)
!   22 DO 23 I=1,29
!        ALFA(I)=ALFR(I)
!   23 END DO

      I=1
      J=1 
   24 FLS(J)= PURES(I)*ABFA 
      CALL DRAW(WC,WS,WL,.6,-1.,FLS(J),ABGR,1) 
      CALL DRAW(WCI(J),WSI(J),WLI,-.6,-1.,FLS(J),ABGR,1) 
      WCI(J)= WCI(J)+WC 
      WSI(J)= WSI(J)+WS 
      WLI = WLI+WL 
      FLA(J)= PURES(I+1) * ABFA 

      IF(FLA(J).GT.0.)GO TO 26
      DRAK(J)=0. 
      DRAM(J)=1. 
      GOTO 34

   26 WI = COSG(ABGR*FLA(J)) 
      IF(PURES(I+2)-1.0) 27,30,29 
   27 drak(j)=0.1*pures(I+3)
   28 dram(j)=0.1*pures(i+4)
      go to 34

   29 DRAK(j)=((0.1*pures(i+4))**(-10./PURES(I+3))-1.0)*(1.0+WI)/(1.0-WI)
      drak(j)=Rund(drak(j), 1000.0)
      dram(j)=0.1*pures(i+3)
      GO TO 34

   30 AA=0.5*(1.0-WI)*pures(i+3)
      WILN=LOG(0.1*pures(i+4)) 
      fmit=.5
      mit=0

      DO
        fm= -wiln/LOG(1.0+aa/fmit)
        mit=mit+1
        IF (ABS(fm-fmit) < 1E-6) GO TO 33
        fmit=fm
      END DO

   33 DRAM(J)=RUND(FM,1000.) 
      DRAS=.05*PURES(I+3)*(WI+1.0)/FM
      drak(j)=Rund(dras, 1000.0)

   34 I=I+5 
      J= J+1 
      IF(J-3)24,38,38 

   38 MER = 0
      WSI(2)= -WSI(2) 
      D(1) = WLI *(WSI(2)+WSI(1)) 
      D(2) =-WLI *(WCI(2)+WCI(1)) 
      D(3) = WCI(1)*WSI(2)-WCI(2)*WSI(1) 
      ITMOD=INT(PURES(11)) 
      RUF=100. 
      IF(ITMOD.GE.4.AND.ITMOD.LE.6)RUF=1000. 
      ITMR=ITMOD 
      SHKS = 0.1*PURES(12) 
      HKST=0.1*ABS(PURES(13)) 

   35 ac(1,1:3)=0.0

      ALIV = 0.
      SINAI = 0. 
      COSAI = 1. 
      FNI = 0. 
      J=1 

   37 CSAIP = COSG(2.*ALFA(J))
      SNAIP = SING(2.*ALFA(J)) 
      IF(J-JST-1)40,39,40 
   39 AC(2,1)= SINAI 
      AC(2,2)= -1.-COSAI 
      AC(2,3)= -1. 
      AC(3,1)= -SNAIP 
      AC(3,2)= 1.+CSAIP 
      AC(3,3)= 1. 
      AC(4,1)= COSAI-CSAIP 
      AC(4,2)= SINAI-SNAIP 
      AC(4,3)= 0. 
      ALIS = ALIV 
      ALISP = ALFA(J) 
      GOTO 41 

   40 FII = CSLG(HABGR*FNI-90.,ALIV)
      FIIP= CSLG(HABGR*FNI-90.,ALFA(J)) 
      PB=FNI*HABGR*BOGEN 
      AC(1,1)=-FIIP*SNAIP+FII*SINAI+(COSAI-CSAIP)*PB +AC(1,1) 
      AC(1,2)=-FII*(1.+COSAI)+FIIP*(1.+CSAIP)+(SINAI-SNAIP)*PB+AC(1,2) 
      AC(1,3)= FIIP - FII + AC(1,3) 


   41 IF(J-JAB-1)42,43,43
   42 ALIV =ALFA(J) 
      SINAI=SNAIP 
      COSAI=CSAIP 
      FNI = ANI(J) 
      J=J+1 
      GO TO 37 

   43 DO 47 J=1,2
      IF(FLA(J)) 47,47,49 
   49 CALL DRAW(WC,WS,WL,DRAK(J),DRAM(J),FLA(J),ABGR,0) 
      AC(1,1)= WC+ AC(1,1) 
      IF(J-2)45,44,45 
   44 WS = -WS 
      WL =- WL 
   45 AC(1,2) = WS +AC(1,2) 
      AC(1,3) =-WL +AC(1,3) 
   47 END DO 


  DO J=1,4
    A(J) = 0. 
    DO I=1,3 
      A(J) = A(J)+D(I)*AC(J,I)
    END DO
  END DO

!   SOLUTION OF TRANSCENDENTAL EQUATION
  53 I=0
     fv=9E9
      PHISH = 0.5 *(ALIS+ALISP)

   60 CSLI = CSLG(PHISH,ALIS) 
      CSLIP= CSLG(PHISH,ALISP) 
      FP=A(1)+A(2)*CSLI+A(3)*CSLIP+A(4)*BOGEN*(90.0+PHISH)
      WRITE(DBG,*) "Trapro iteration", fp,fv,phish
      IF(i >= 20) GO TO 66
      IF (ABS(fp)-ABS(fv) < -0.5E-9) GO TO 62
      i=20
      PHISH=PHISH-PDIF
      GO TO 60
      
   62 pdif=-fp / (a(2)/(phish-alis)  +  a(3)/(phish-alisp))
      i=i+1
   65 fv=fp
      phish=phish+pdif
      IF (PHISH.LT.ALIS .AND. PHISH.GT.ALISP)GO TO 60

      WRITE(IDRU,*) " TRANCENDENTAL EQUATION HAS DIVERGED"
      WRITE(IDRU,*) " CHECK TRA1 and TRA2 input records"
      WRITE(IDRU,*) " Iteration",mer, "      mode=",itmod
      WRITE(IDRU,*) " alis=",alis, "   alisp=",alisp, "   phish=",phish
      STOP "Divergence"


   66 ANI(JST) = (PHISH+90.)/HABGR
      DO 71 I=1,3 
        FINT(I)=AC(1,I)+AC(2,I)*CSLI+AC(3,I)*CSLIP+AC(4,I)*BOGEN*(PHISH+90.)
   71 END DO   

   69 HK(1)=(FINT(1)*WLI-FINT(3)*WCI(2))/D(2)
      HK(2)=(FINT(1)*WLI+FINT(3)*WCI(1))/D(2) 
      HKS = HK(1)+HK(2) 
      IF(ITMOD.EQ.0 .OR. ABS(HKS-SHKS).LT.HKST) GO TO 74 
      IF(MAGAM.LT.2 .AND. MAGAM-MER.NE.1) GO TO 100 
      GO TO 76 


   74 ITMOD =0    ! jump back here from below
      IF(MAGAM .EQ. 0) GO TO 300 

   76 continue  
      WRITE(IDRU,77) airfoilID,MER,ITMR 
   77 FORMAT ("TRANSCENDENTAL EQUATION RESULTS FOR AIRFOIL ",A,      &
     &"    ITERATION",I2,"   MODE",I2)                                 
      WRITE(IDRU,*) & 
      "   NU   ALPHA*  OMEGA'  OMEGA     K       MU     K H   LAMBDA LAMBDA*"
      JH= 1 

      DO 85 JN=1,JAB
!        NZT=NZPZ(1,0) 
        IF (JN.NE.1 .AND. JN.NE.JAB) GO TO 83 
        X1 = 0.5*(1.0+ COSG(FLA(JH)*ABGR)) 
        WHK=(1.0+DRAK(JH)*(1.0-X1)/X1)**(-DRAM(JH)) 
        WSTR= DRAM(JH)*DRAK(JH)/X1 
        WRITE(IDRU,82) ANI(JN),ALFA(JN),WSTR,WHK, &
        DRAK(JH),DRAM(JH),HK(JH),FLA(JH),FLS(JH)
   82 FORMAT(F8.3,F7.2,F8.3,F7.3,2F8.3,F9.6,2F7.2) 
        JH=2 
        GO TO 85
   83 WRITE(IDRU,82) ANI(JN),ALFA(JN)
   85 END DO


!   86 FORMAT(1H+,72X,13HWARNING, NUE(,I2,16H) NOT INCREASING)


      IF(ITMOD==0) GO TO 300

  100 IF(MER)103,102,103 
  102 DAL = 0.1 
      GO TO 104 

  103 IF(HKS-HKSV.EQ.0.)GO TO 74

      DAL = (SHKS-HKS)*DAL/(HKS-HKSV) 
      DALD=DAL 
      IF (ITP.EQ.0) DAL=RUND(DAL,RUF) 
      IF (MAGAM.EQ.0.) GO TO 1004 
!      NZT=NZPZ(2,0) 
      WRITE(IDRU,1003) MER,HKS,DALD,DAL 
 1003 FORMAT (" ITERATION",I2,"   K S =",F9.6,"   DELTA =",F12.8,     &
       "    ROUNDED =",F7.3)                                              
 1004 IF (DAL.EQ.0.) GO TO 74 
      IF (MER.GE.3 .AND. ABS(DALV).LE.ABS(DAL)) GO TO 74 

  104 DALV=DAL
      IF (ITMOD.GE.4) GO TO 113 
      DO 111 J=1,JAB 
        IF (ITMOD.NE.2 .AND. J.LE.JST) ALFA(J)=ALFA(J)+DAL 
        IF (ITMOD.NE.1 .AND. J.GT.JST) ALFA(J)=ALFA(J)-DAL 
  111 END DO 
      GO TO 112 

  113 IF (ITMOD.GE.7) GO TO 114
      IF (ITMOD.NE.5) DRAK(1)=DRAK(1)+DAL 
      IF (ITMOD.NE.4) DRAK(2)=DRAK(2)+DAL 
      GO TO 112 

  114 IF (ITMOD.NE.8) ALFA(JST)=ALFA(JST)+DAL
      IF (ITMOD.NE.7) ALFA(JST+1)=ALFA(JST+1)-DAL 

  112 HKSV=HKS
      MER =MER+1 
      GOTO 35         ! big jump

  300 AK=0.5*(COSG(PHISH-ALFA(JST+1))/SING(PHISH-ALFA(JST+1))           &
     &       -COSG(PHISH-ALFA(JST))/  SING(PHISH-ALFA(JST)))              
      AKP  =AK*180./9.8696044         ! magic number??
      PHIM = 0. 
      NU=1 
      I= 1 
      ANU =0. 
      JH=0 
      VI=0.0

  302 JH=JH+1 
      FF1 = COSG(ABGR*FLA(JH)) 
      FF2 = DRAK(JH)/(1.+FF1) 
      FG1 = COSG(ABGR*FLS(JH)) 
      FG3 = 0.6/(FG1-1.0)

  304 VI= VI - CSLG(PHIM-90.0, ALFA(I)) 
      GO TO 310 

  306 ARGN = ANU
      IF (ANU .GT. 0.5* ABSZ) ARGN=ABSZ-ANU 
      CSP = COSG(ARGN*ABGR) 
      F=0. 
      IF (ARGN.LT.FLA(JH)) F=DRAM(JH)*LOG((CSP-FF1)*FF2+1.0) 
      G=0. 
      IF (ARGN.LT.FLS(JH)) G=-HK(JH)*LOG(1.-((CSP-FG1)*FG3)**2) 
      P(NU)= F+G+CSLG(ANU*HABGR-90., ALFA(I)) + VI 
      P1(NU)=P(NU)-AK*ABS(SING((ANU*HABGR - 90.0) - PHISH)) 
      NU = NU + 1 
      ANU= ANU+ 1.0 
  310 IF(ANU-ANI(I))306,306,312 
  312 IF(ANU- ABSZ)314,320,320 
  314 PHIM = ANI(I)*HABGR 
      VI=VI+CSLG(PHIM-90.,ALFA(I)) 
      I = I+1 
      IF(I-1-JST)304,302,304 
  320 PS=0. 
      B2=0. 
      DO 324 I=1,NKR 
        PS=PS+P(I) 
        BI = 2*(I-1) 
        B2 = B2 + SING(BI*ABGR)*P(I)
  324 END DO

      V1 = 2.*EXP(PS/ABSZ) 
      SXI = .00000000     ! ???
      SY=0. 

      DO 328 N=1,NQ
      Q=0. 
      DO 326 M=1,IB 
        MN = N + 1 + MQ - 2*M 
        MM = 2*N - MN 
        IF(MN.GT.NKR) MN = MN - NKR 
        IF(MM.LT.1) MM = MM + NKR 
        Q = Q+ FKERN(M)*(P1(MN)-P1(MM))
  326 END DO

      ANU= N-1
      ZP = ANU*HABGR - 90. 
      ZL = COSG(ZP - PHISH) 
      ZL = ABS((1.-ZL)/(1.+ZL)) 
      IF (ZL.NE.0.) ZL=LOG(ZL) 
      ARG(N) = Q - AKP*SING(ZP-PHISH)*ZL + ZP 
      VF(N) = V1*EXP(-P(N)) 
      WV = COSG(ZP)/VF(N) 
      XP(N)= WV*SING(ARG(N)) 
      YP(N)=-WV*COSG(ARG(N)) 
      SXI= SXI+ XP(N) 
      SY = SY + YP(N)
  328 END DO

      SX = SXI
      XPK = SX/(ABSZ  - 1.) 
      YPK  = SY/(ABSZ -1.) 
      DO 329 N=2,NKR 
        XP(N)=XP(N)-XPK 
        YP(N)=YP(N)-YPK
  329 END DO

      CALL CINT(XP,X,NQ)   ! did have izz as argument 
      CALL CINT(YP,Y,NQ)   ! did have izz as argument 
      RQV = 0. 
      DO 330 N=2,NKR 
        RQ=X(N)*X(N)+Y(N)*Y(N) 
        IF(RQ.GT.RQV)L=N 
        RQV = RQ
  330 END DO

      DO 327 I = 1,3
        IEPPL = L-2+I 
        R(I)=SQRT(X(IEPPL)*X(IEPPL)+Y(IEPPL)*Y(IEPPL))
  327 END DO

  333 TAU = (R(3)-R(1))/(4.*(R(2)+R(2)-R(1)-R(3)))
      XNAS = X(L)+TAU*(X(L+1)-X(L-1)+2.*TAU*(X(L+1)+X(L-1)-X(L)-X(L))) 
      YNAS = Y(L)+TAU*(Y(L+1)-Y(L-1)+2.*TAU*(Y(L+1)+Y(L-1)-Y(L)-Y(L))) 
      SQ = XNAS*XNAS + YNAS*YNAS 
      AT=XNAS/SQ 
      B= YNAS/SQ 
      STREF = 1./SQRT(SQ) 
      ETA =  ABSZ*STREF/PI 
      CM = .5*ETA*STREF*B2 
      DARG = 19.09859 *(3.*YNAS/XNAS - (YNAS/XNAS)**3) 
      IF(ABS(SX)+ABS(SY).LT..0001*ABSZ) GO TO 335 
      SX=STREF*SX*200. 
      SY=STREF*SY*200. 
!      NZT=NZPZ(2,0) 
      WRITE(IDRU,334) SX,SY 
  334 FORMAT (" WARNING - SX =",F6.3,"   SY =",F6.3) 

  335 DO 331 N=2,NQ
      XR=X(N) 
      X(N)= 1.-B*Y(N)-AT*XR 
      Y(N)= B*XR -AT*Y(N) 
      ARG(N) = ARG(N) - DARG 
      WQ = (XP(N)+XP(N-1)-XPK-XPK)**2 + (YP(N)+YP(N-1)-YPK-YPK)**2 
      DS(N-1) = STREF*SQRT(WQ)*(1.+.6666667*((XP(N)*YP(N-1)             &
     &-XP(N-1)*YP(N))/WQ)**2)
  331 END DO

      NHKW=NQ/12
      DLT = Y(NHKW)/(BOGEN*(1.-X(NHKW))) 
      NHKW=NQ-NHKW+1 
      DLTU=-Y(NHKW)/(BOGEN*(1.-X(NHKW))) 
      X(1) = 1. 
  332 ARG(1)=ARG(1)-DARG 
  346 ITP=1 
      ALN=DARG 

      IF (PURES(13).GE.0.) GO TO 11
      PURES(12)=10.0*HKS
      PURES(13)=0.00001 
   11 RETURN 
END Subroutine Trapro   ! ---------------------------------------------------

!+
SUBROUTINE UMP(U,D2,H,FUM)
! ---------------------------------------------------------------------------
! PURPOSE -
! called by GRS
IMPLICIT NONE

  REAL,INTENT(IN):: u
  REAL,INTENT(IN):: d2
  REAL,INTENT(IN):: h
  REAL,INTENT(OUT):: fum

  real:: DUM,WRE,DU,US,DUMY
  integer:: m,md
  COMMON/BL/DUM(5),WRE,M,MD,DU(7),US,DUMY

  INTEGER:: mxt
  INTEGER:: n
!----------------------------------------------------------------------------
  IF (m <= 0.0) THEN
    FUM=-18.4
  ELSE IF (m >= 3.0) THEN
    FUM=LOG(WRE*WRE*D2*U)+21.74+0.36*REAL(M-3)-18.4*H
  ELSE
    MXT=2*M-1+(ND+1)/2 
    N=NU-ND 
    FUM=X(N)+SUMP*(X(NU)-X(N))/DLV -XTRI(MXT)
  END IF
  RETURN 
END Subroutine Ump   ! ------------------------------------------------------

!+
SUBROUTINE UpCase(a)
! ---------------------------------------------------------------------------
! PURPOSE - Convert lower case elements of a character variable to upper
! called by TestKeyWord
  CHARACTER(LEN=*),INTENT(IN OUT):: a
  INTEGER:: j,k
!----------------------------------------------------------------------------
  DO k=1,LEN(a)
    j=IACHAR(a(k:k))
    IF (j>96 .AND. j<123) a(k:k)=ACHAR(j-32)    ! assumes ASCII character set
  END DO
  RETURN
END Subroutine UpCase   ! ---------------------------------------------------

!+
SUBROUTINE VISC (I,J,CL,CD,CM)
! ---------------------------------------------------------------------------
! PURPOSE -
! called by GRP and ProcessCommands(PLWA)
  INTEGER,INTENT(IN):: i,j
  REAL,INTENT(OUT):: cl,cd,cm

  REAL:: alrc
  REAL:: d1,d2
  REAL:: dcl
  REAL:: sa1,sa2
!----------------------------------------------------------------------------
  ALRC=ALV(I)                            ! alv is global
  D1=MAX(ALRC-ALN+DLT,  0.0)             ! dlt,dltu are global
  D2=MIN(ALRC-ALN-DLTU, 0.0)
  SA1=SA(J,1,I)                           ! sa is global
  SA2=SA(J,2,I)
  DCL=0.5*(D1*SA1+D2*SA2)
  CL=0.11*(ALRC-DCL)                            ! cl,cd,cm are arguments
  CD=ABS(CW(J,1,I))+ABS(CW(J,2,I))                 ! cw is global
  CM=CML(I)+.0275*DCL*(ABS(1.0-SA1-SA2)**1.5)      ! cml is global
  RETURN
END Subroutine Visc   ! -----------------------------------------------------

!+
FUNCTION VPR(N,I) RESULT(velocity)
! ---------------------------------------------------------------------------
! PURPOSE -
! called by Diagram,, GRP, Moment

  INTEGER,INTENT(IN):: i,n
  REAL:: velocity

!      COMMON GAP(876),VF(121),ARG(121),GAPS(57),IZZ,KFU,NQ

!      COMMON/PRAL/DLT,DLTU,ALN,ALV(14),NAL,ITP,CML(14),CRL(14),         &
!     &DARF,ITIT1,ITIT2,D1I(121,5),AVI(5,5),BVI(5,5),                    &
!     & FMACH,BEMACH,IVI,NAIT(5),XDA,YDA,DEFLG,MOMAG,NAI                 

!      DATA DALR/99./
!  REAL,SAVE:: alr=99.0, anr=99.0
  REAL:: alr
  REAL:: cast,sast
!----------------------------------------------------------------------------
  IF (itp <= 1) THEN
    velocity=VF(N)*COSG(180.0*REAL(N-1)/REAL(NQ-1)-ALV(I)) 
  ELSE
    alr=alv(i)
    cast=COSG(alr-aln)
    sast=SING(alr-aln)
!    anr=aln
    velocity=cast*vf(n)+sast*arg(n)
  END IF
  RETURN 
END Function Vpr   ! --------------------------------------------------------

!+
SUBROUTINE Welcome()
! ---------------------------------------------------------------------------
! PURPOSE -
! called by main program
  CHARACTER(LEN=*),PARAMETER:: GREETING = &
    "PROFILE - Design and Analysis of Low Speed Airfoils"
  CHARACTER(LEN=15):: dateTime
  INTEGER:: errCode
  CHARACTER(LEN=133):: fileName
!----------------------------------------------------------------------------
  WRITE(*,*) GREETING
  dateTime=GetDateTimeStr()

  DO
    WRITE(*,*) "Enter the name of the input file:"
    READ(*,'(A)') fileName
    IF (Len_Trim(fileName)==0) STOP
    OPEN(UNIT=ILES, FILE=fileName, STATUS='OLD', &
      IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
    IF (errCode==0) EXIT
    OPEN(UNIT=ILES, FILE=Trim(fileName)//".dat", STATUS='OLD', &
      IOSTAT=errCode, ACTION='READ', POSITION='REWIND')
    IF (errCode==0) EXIT
    WRITE(*,*) "Unable to open this file. Try again."
  END DO
  INQUIRE(UNIT=ILES,NAME=fileName)
  WRITE(*,*) "Reading from "//Trim(fileName)

  OPEN(UNIT=DBG, FILE='profile.log', STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
  IF (errCode /= 0) THEN
    WRITE(*,*) "Unable to open profile.log"
    STOP
  END IF

  OPEN(UNIT=IDRU, FILE='profile.out', STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
  IF (errCode==0) THEN
    WRITE(IDRU,*) "EPPLER/SOMERS PROFILE PROGRAM                  "//dateTime
    WRITE(IDRU,*) "INPUT FILE: "//Trim(fileName)
  ELSE
    WRITE(*,*) "Unable to open profile.out"
    STOP
  END IF

  OPEN(UNIT=FIG, FILE='profile.fig', STATUS='REPLACE', &
    IOSTAT=errCode, ACTION='WRITE', POSITION='REWIND')
  IF (errCode /= 0) THEN
    WRITE(*,*) "Unable to open profile.fig"
    STOP
  END IF

  RETURN
END Subroutine Welcome   ! --------------------------------------------------

END Program Profile   ! =================================================

