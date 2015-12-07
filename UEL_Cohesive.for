C
C     ABAQUS user subroutine for a 2D cohesive zone model, based on the paper by
C		Y.F. Gao and A.F. Bower, "A simple technique for avoiding convergence
C		problems in finite element simulations of crack nucleation and growth on
C		cohesive interfaces," Modelling Simul. Mater. Sci. Eng. 12, 453-463, 2004.
C
C	The Fortran code and the user's manual are last modified on December 9, 2004.
C
C		Y.F. Gao & A.F. Bower, Divsion of Engineering, Brown University
C
C=========================== SUBROUTINE UEL ===================
C
      SUBROUTINE UEL(RHS,STIF,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     &       PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     &       KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     &       LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
 	  INCLUDE 'ABA_PARAM.INC'
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      DIMENSION RHS(MLVARX,*),STIF(NDOFEL,NDOFEL),PROPS(*)
      DIMENSION SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL)
      DIMENSION DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*)
      DIMENSION JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*)
      DIMENSION PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)
C
C     Only relevant variables are commented. Refer to ABAQUS manual for others.
C
C     Variables with intent(out)
C     RHS         Residual Vector
C     STIF        Stiffness matrix
C
C     Variables with intent(in)
C     PROPS       Element property array for the Xu-Needleman model
C                 PROPS(1)    SIGMA_max
C                 PROPS(2)    Delta_n
C                 PROPS(3)    Delta_t
C                 PROPS(4)    Q
C                 PROPS(5)    R 
C     COORDS      Nodal coordinate array: COORD(J,N) is jth coord of nth node
C     U           Total accumulated DOF array.  Contains accumulated displacements,
C                 ordered as (u_i^1, u_i^2)
C     DU          Incremental displacements, ordered as
C                 DU(2*(N-1)+I,1) = u_i^n
C
C     NNODE       No. nodes on element.
C     NDOFEL      No. degrees of freedom for element
C     NPROPS      No. real valued element properties
C     MLVARX      Dimensioning variable.
C     NRHS        No. RHS vectors.
C     MCRD        Largest of max value of COORDINATES parameter or active DOF <3.
C
C	NNODE=4		Four-point element (4,3 should collapse onto 1,2 respectively)
C
C				4----3
C				|	 |
C				1----2
C
C	NDOFEL=8	Two-dimensional, d.o.f=4*2
C	NPROPS=5	Five parameters passed in
C	MLVARX=1
C	NRHS=1
C     MCRD=2		(X,Y)
C
      DIMENSION W(4),PQ(4)
      DIMENSION F(10),DF(10)
      DIMENSION GAP(2),TRACT(2),DTDG(2,2)
	DIMENSION RNM(2),TANGENT(2),UREL(2)
	DIMENSION VGAP(2),VREL(2)

	NINTP=2
C
C	The above arrays are used inside the subroutine
C	NINTP		No. integration points
C	PQ(NINTP)	Local coordinates of integration points
C	W(NINTP)	Weights of integration points
C	F(NNODE)	Shape function
C	DF(NNODE)	Derivative of shape function w.r.t local coordinate
C

       DO  I = 1,NDOFEL
         RHS(I,1)=0.D0
         DO  J = 1,NDOFEL
           STIF(J,I) = 0.D0
         END DO
       END DO
C
C      Set up integration points and weights
       CALL KINTPT(PQ,W,NNODE,NINTP)
C
       DO I=1,NINTP
C
C        Shape functions and derivatives
         CALL KSHAPE(NNODE,PQ(I),F)
	   CALL DSHAPE(NNODE,PQ(I),DF)
C
C        Compute Normal and tangent vectors to boundary at int pt.
C        Plane of interface defined as average of upper and lower
C        surfaces (hopefully coincident, but you never know...)

	   TANGENT(1)=.5D0*(COORDS(1,2)-COORDS(1,1)
     &					+COORDS(1,3)-COORDS(1,4))
	   TANGENT(2)=.5D0*(COORDS(2,2)-COORDS(2,1)
     &					+COORDS(2,3)-COORDS(2,4))
	
	   CALL KUNITV(TANGENT,DET)
	   DET=.5D0*DET

	   RNM(1)=-TANGENT(2)
	   RNM(2)=TANGENT(1)
C
C        Relative displacement, in global coords
         UREL(1) = 0.D0
         UREL(2) = 0.D0
	   VREL(1)=0.D0
	   VREL(2)=0.D0
         DO N = 1,NNODE/2
           N2 = N+NNODE/2
	     DO J = 1,2
           UREL(J) = UREL(J) + F(N2)*(U(2*(N2-1)+J)) 
     &                       - F(N)*(U(2*(N-1)+J)) 
		 VREL(J) = VREL(J) + F(N2)*(V(2*(N2-1)+J)) 
     &                       - F(N)*(V(2*(N-1)+J)) 
	     END DO
         END DO

C
C        GAP(1) is normal separation, GAP(2) is separation in tan dirn
         GAP(1) =  RNM(1)*UREL(1)+RNM(2)*UREL(2)
         GAP(2) =  TANGENT(1)*UREL(1)+TANGENT(2)*UREL(2)
	   VGAP(1)=RNM(1)*VREL(1)+RNM(2)*VREL(2)
	   VGAP(2)=TANGENT(1)*VREL(1)+TANGENT(2)*VREL(2)
C
         CALL SEPLAW(PROPS,GAP,VGAP,TRACT,DTDG,DTIME)
C
         SIG = 1.D0
         DO N = 1,NNODE
	     IF (N.GT.NNODE/2) SIG = -1.D0
	     DO K = 1,2
             RHS(2*(N-1)+K,1) = RHS(2*(N-1)+K,1) 
     &          + SIG*F(N)*( TRACT(1)*RNM(K)
     &          +            TRACT(2)*TANGENT(K))*W(I)*DET
	     END DO
	   END DO


         SIGN = 1.D0
         DO N = 1,NNODE
	     IF (N.GT.NNODE/2) SIGN = -1.D0
	     SIGM = 1.D0
	     DO M=1,NNODE
	       IF (M.GT.NNODE/2) SIGM = -1.D0
	       DO KN = 1,2
	       DO KM = 1,2
	       ICOL = 2*(M-1) + KM
	       IROW = 2*(N-1) + KN
C
             STIF(ICOL,IROW) = STIF(ICOL,IROW)
     &         +( (DTDG(1,1)*RNM(KN)+DTDG(2,1)*TANGENT(KN))*RNM(KM) 
     &         +(DTDG(1,2)*RNM(KN)+DTDG(2,2)*TANGENT(KN))*TANGENT(KM) )
     &              *SIGM*SIGN*F(M)*F(N)*W(I)*DET
             END DO
	       END DO

C
           END DO
         END DO

C
       END DO

C
       RETURN
       END
C
C===================== SUBROUTINE SEPLAW =====================
C
      SUBROUTINE SEPLAW(PROPS,GAP,VGAP,TRACT,DTDGAP,DTIME)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION PROPS(5)
      DIMENSION GAP(2),VGAP(2),TRACT(2),DTDGAP(2,2)
C
C     Subroutine to specify traction-separation law for 
C     a debonding interface.
C
C     Currently coded for Xu-Needleman constitutive law, i.e.
C		X.P. Xu and A. Needleman, J. Mech. Phys. Solids, 42, 1397-1434.
C
C     PROPS(1)    SIGMA_max
C     PROPS(2)    Delta_n
C     PROPS(3)    Delta_t
C     PROPS(4)    Q
C     PROPS(5)    R
C
C     GAP(1)      Normal separation
C     GAP(2)      Tangential separation
C     TRACT(1)    Normal traction
C     TRACT(2)    Tangential traction
C     
C
C --  Work of separation
      SEPWRK = DEXP(1.D0)*PROPS(1)*PROPS(2)
      DN = PROPS(2)
      DT = PROPS(3)
      Q = PROPS(4)
      R = PROPS(5)
C
C
      C1 = ( 1.D0-DEXP(-GAP(2)*GAP(2)/(DT*DT)) )
      C1 = C1 * (1.D0-Q)/(R-1.D0) * (R - GAP(1)/DN)
      C2 = (GAP(1)/DN)*DEXP(-GAP(2)*GAP(2)/(DT*DT))
      TRACT(1) = (SEPWRK/DN)*DEXP(-GAP(1)/DN)*(C2+C1)
C
      C1 = Q + (R-Q)/(R-1.D0)*(GAP(1)/DN)
      C1 = C1 * DEXP(-GAP(1)/DN) * DEXP(-GAP(2)*GAP(2)/(DT*DT))
      C1 = C1 * 2.D0*(DN/DT)*(SEPWRK/DN)
	TRACT(2) = C1*GAP(2)/DT
C
      C1 = (1.D0-Q)/(R-1)*(1.D0-DEXP(-GAP(2)*GAP(2)/(DT*DT)))
      C1 = C1 * (R+1.D0-GAP(1)/DN)
      C1 = (1.D0-GAP(1)/DN)*DEXP(-GAP(2)*GAP(2)/(DT*DT)) - C1
      DTDGAP(1,1) = (SEPWRK/(DN*DN))*DEXP(-GAP(1)/DN) * C1
C
      C1 = Q + (GAP(1)/DN) * (R-Q)/(R-1.D0)
      C1 = C1 * DEXP(-GAP(1)/DN)*DEXP(-GAP(2)*GAP(2)/(DT*DT))
      C1 = 2.D0*(SEPWRK/(DT*DT)) * C1
	DTDGAP(2,2) = C1*(1.D0-2.D0*(GAP(2)*GAP(2)/(DT*DT)))
C
      C1 = -GAP(1)/DN + (1.D0-Q)/(R-1.D0)*(R-GAP(1)/DN)
      C1 = C1 * DEXP(-GAP(1)/DN)*DEXP(-GAP(2)*GAP(2)/(DT*DT))
	  C1 = 2.D0*(SEPWRK/(DT*DN))*C1
      DTDGAP(1,2) = (GAP(2)/DT)*C1
      DTDGAP(2,1) = DTDGAP(1,2)
C
C	ZETA is the fictitious viscosity used to regularize the instability
C		problem. In our papere, \zeta_n is ZETA*PROPS(1) here. Simply
C		setting ZETA=0.D0 (or delete the following three lines
C		gives the standard Xu-Needleman model.
      if (dtime>0.d0) then
          ZETA=0.0001D0
	  TRACT(1)=TRACT(1)+ZETA*PROPS(1)*VGAP(1)/DN
	  DTDGAP(1,1)=DTDGAP(1,1)+ZETA*PROPS(1)/DN/DTIME
      endif
C
      RETURN
      END
C
C
C========================== SUBROUTINE INTPTS =================================
C
       SUBROUTINE KINTPT(PQ,W,KNODE,NINTP)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
       DIMENSION PQ(*),W(*)
C
C      Subroutine to initialize integration points PQ and weights W
C      NINTP = No. integration points (DEFAULT=2)
C      KNODE = No. of nodes on element (DEFAULT=4)
C
C      Integration points for 2D interface elements
C
C				4----3
C				|	 |
C				|	 |
C				1----2
c
	IF (NINTP.EQ.2) THEN
		CN=0.5773502691896260D0
		PQ(1)=-CN
		PQ(2)=CN
		W(1)=1.D0
		W(2)=1.D0
	ENDIF

C
       RETURN
       END
C
C=========================== SUBROUTINE SHAPE ===============================
C
       SUBROUTINE KSHAPE(KNODE,PQ,F)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
       DIMENSION F(*)
C
C  
	IF(KNODE.EQ.4) THEN
		F(1)=.5D0*(1-PQ)
		F(2)=.5D0*(1+PQ)
		F(3)=F(2)
		F(4)=F(1)
	ENDIF

       RETURN
       END
C
C============================== SUBROUTINE DSHAPE ====================
C
       SUBROUTINE DSHAPE(KNODE,PQ,DF)
       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
       DIMENSION DF(*)
C
C
	IF(KNODE.EQ.4) THEN
		DF(1)=-.5D0
		DF(2)=.5D0
		DF(3)=DF(2)
		DF(4)=DF(1)
	ENDIF

       RETURN
       END
C
C========================== SUBROUTINE KUNITV ===================
C
       SUBROUTINE KUNITV(A,AMAG)
	DOUBLE PRECISION A(2), AMAG
C
C     Normalize vector A and return its magnitude as AMAG
C

      AMAG = DSQRT(A(1)*A(1)+A(2)*A(2))
	A(1) = A(1)/AMAG
	A(2) = A(2)/AMAG
C
       RETURN
	 END
