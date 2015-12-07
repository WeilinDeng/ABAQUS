C   RUN in the ABAQUS command:
C     > abaqus double user=uel_bodyforce.for 
C       input=user_element_bodyforce.inp job=job1 -inter
C
C   ABAQUS subroutine for rotational body force

C=========================== SUBROUTINE UEL ===================
C
      SUBROUTINE UEL(RHS,STIF,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     &       PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     &       KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     &       LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C

 	INCLUDE 'ABA_PARAM.INC'
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)      

      double precision RHS(MLVARX,*),STIF(NDOFEL,NDOFEL),PROPS(*)
      double precision SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL)
      double precision DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2)      
      double precision ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*)
      double precision PREDEF(2,NPREDF,NNODE),PARAMS(*)

      double precision DTIME
      double precision PNEWDT 
      double precision PERIOD

      integer NDOFEL
      integer NRHS
      integer NSVARS
      integer NPROPS
      integer MCRD
      integer NNODE
      integer JTYPE
      integer KSTEP
      integer KINC
      integer JELEM
      integer NDLOAD
      integer NPREDF
      integer MLVARX
      integer MDLOAD
      integer NJPROP

      integer LFLAGS(*),JPROPS(*)
      integer JDLTYP(MDLOAD,*)
C
C     Only relevant variables are commented. Refer to ABAQUS manual for others.
C
C     Variables with intent(out)
C     RHS         Residual Vector
C     STIF        Stiffness matrix
C
C     Variables with intent(in)
C     PROPS       Element property array for the Xu-Needleman model
C                 PROPS(1)    Mass density
C                 PROPS(2)    Rate of change of Angular speed
C                 PROPS(3:5)    n1:n3
C                 PROPS(6:8)    x01:x03
C     COORDS      Nodal coordinate array: COORD(J,N) is jth coord of nth node
C     U           Total accumulated DOF array.  Contains accumulated displacements,
C                 ordered as (u_i^1, u_i^2)
C     DU          Incremental displacements, ordered as
C                 DU(2*(N-1)+I,1) = u_i^n
C
C     NNODE       No. nodes on element.
C     NDOFEL      No. degrees of freedom for element (total)
C     NPROPS      No. real valued element properties
C     MLVARX      Dimensioning variable.
C     NRHS        No. RHS vectors.
C     MCRD        Largest of max value of COORDINATES parameter or active DOF <3.
C

C     Local variables
        double precision nodal_u(MCRD,NDOFEL/MCRD) 
        double precision xi(MCRD,27)         ! Integration points
        double precision w(27)               ! Integration weights
        double precision N(20)               ! Shape functions
        double precision dNdxi(20,MCRD)      ! Shape function derivatives
        double precision dxdxi(MCRD,MCRD)    ! Jacobian
        double precision dxidx(MCRD,MCRD)    ! Jacobian inverse
        double precision determinant         ! Jacobian determinant

C     element properties
        double precision rho                 
        double precision Omega              
        double precision axis(3)            
        double precision x0(3)   
                    
        double precision temp,xk(MCRD),uk(MCRD),yk(MCRD),delta(3,3)
        integer n_points,kint,ka,ki,kb,kk
        
        nodal_u = reshape(U,(/MCRD,NDOFEL/MCRD/))
        delta = reshape((/ 1.D0, 0.D0, 0.D0, 
     +                     0.D0, 1.D0, 0.D0, 
     +                     0.D0, 0.D0, 1.D0 /), shape(delta))	
        if (NPROPS<8) then
            write(6,*) ' Not enough properties were provided for UEL '
            write(6,*) ' NPROPS ',NPROPS
            write(6,*) ' PROPS ',PROPS(1:nprops) 
            stop
        endif

        rho = PROPS(1)
        Omega = PROPS(2)*TIME(2)
        axis(1:3) = PROPS(3:5)
        x0(1:3) = PROPS(6:8)

        RHS(1:NDOFEL,1) = 0.d0
        STIF(1:NDOFEL,1:NDOFEL) = 0.d0

        if (MCRD==2) then
            if (nnode == 3) n_points = 4
            if (nnode == 4) n_points = 4
            if (nnode == 6) n_points = 4
            if (nnode == 8) n_points = 4
            if (nnode == 9) n_points = 9
        else if (MCRD==3) then
            if (nnode == 4) n_points = 4
            if (nnode == 10) n_points = 5
            if (nnode == 8) n_points = 27
            if (nnode == 20) n_points = 27
        else
            write(6,*) ' Received strange number of coordinates in UEL '
            stop
        endif

C
C       Set up integration points and weights
        call initialize_integration_points(n_points,NNODE,MCRD, 
     +                   xi(1:MCRD,1:n_points), w(1:n_points))
     
C
        DO kint=1,n_points
            call calculate_shapefunctions(xi(1:MCRD,kint),
     +                   nnode,MCRD,N,dNdxi)
            dxdxi = matmul(COORDS(1:MCRD,1:nnode),dNdxi(1:nnode,1:MCRD))
            call invert_small(MCRD,dxdxi,dxidx,determinant)
  
            ! find the ref. coords, displcements, 
            ! and current coords of integration points
            xk(1:MCRD) = matmul(COORDS(1:MCRD,1:NNODE),N(1:NNODE))
            uk(1:MCRD) = matmul(nodal_u(1:MCRD,1:NNODE),N(1:NNODE))
            yk = xk + uk
      
            ! find the residul vector
            do ka = 1, NNODE
                temp = 0.d0
                do ki = 1, MCRD
                   temp = yk(ki)-x0(ki)-dot_product(yk-x0,axis)*axis(ki)
                   RHS(ki+(ka-1)*MCRD,1) = RHS(ki+(ka-1)*MCRD,1) +
     +                     rho*Omega**2*temp*N(ka)*w(kint)*determinant
C                  temp = xk(ki)-x0(ki)-dot_product(xk-x0,axis)*axis(ki)+
C     +                    uk(ki)-x0(ki)-dot_product(uk-x0,axis)*axis(ki)
C                  RHS(ki+(ka-1)*MCRD,1) = RHS(ki+(ka-1)*MCRD,1) +
C     +                    rho*(Omega**2)*temp*N(ka)*w(kint)*determinant
                enddo
            enddo 
            ! find the stiffness matrix
            do ka = 1, NNODE
                do ki = 1, MCRD
                    do kb = 1, NNODE
                        do kk = 1, MCRD
                            STIF(ki+(ka-1)*MCRD,kk+(kb-1)*MCRD) =
     +                       STIF(ki+(ka-1)*MCRD,kk+(kb-1)*MCRD)-
     +                       rho*(Omega**2)*N(ka)*N(kb)*( delta(ki,kk)-
     +                       axis(ki)*axis(kk) )*w(kint)*determinant
                        enddo  
                    enddo
                enddo
            enddo
      
        END DO
       
      RETURN
      END
      
      
       
!==============================================================================
       subroutine initialize_integration_points(n_points, n_nodes,
     +                      n_coord, xi, w)
       
        implicit none

        integer, intent( in )       :: n_nodes             
        integer, intent( in )       :: n_points            
        integer, intent( in )       :: n_coord             
        double precision, intent( out ) :: xi(n_coord,*)            
        double precision, intent( out ) :: w(*)               

        double precision :: x1D(4),w1D(4), cn, w1, w2, w11, w12, w22
        integer :: i,j,k,n
  
        if (n_coord==1) then                       
        ! Integration points for 1d elements  
            select case ( n_points )
                case (2)
                    xi(1,1) = .5773502691896257D+00
                    xi(2,1) = -.5773502691896257D+00
                    w(1) = .1000000000000000D+01
                    w(2) = .1000000000000000D+01
                    return
                case (3)
                    xi(3,1) = 0.7745966692414834D+00
                    xi(2,1) = .0000000000000000D+00
                    xi(3,1) = -.7745966692414834D+00
                    w(1) = .5555555555555556D+00
                    w(2) = .8888888888888888D+00
                    w(3) = .5555555555555556D+00
                    return
                case (4)
                    xi(1,1) = .8611363115940526D+00
                    xi(2,1) = .3399810435848563D+00
                    xi(3,1) = -.3399810435848563D+00
                    xi(4,1) = -.8611363115940526D+00
                    w(1) = .3478548451374538D+00
                    w(2) = .6521451548625461D+00
                    w(3) = .6521451548625461D+00
                    w(4) = .3478548451374538D+00
                    return
                case (5)
                    xi(1,1) = .9061798459386640D+00
                    xi(2,1) = .5384693101056831D+00
                    xi(3,1) = .0000000000000000D+00
                    xi(4,1) = -.5384693101056831D+00
                    xi(5,1) = -.9061798459386640D+00
                    w(1) = .2369268850561891D+00
                    w(2) = .4786286704993665D+00
                    w(3) = .5688888888888889D+00
                    w(4) = .4786286704993665D+00
                    w(5) = .2369268850561891D+00
                    return
                case (6)
                    xi(1,1) = .9324695142031521D+00
                    xi(2,1) = .6612093864662645D+00
                    xi(3,1) = .2386191860831969D+00
                    xi(4,1) = -.2386191860831969D+00
                    xi(5,1) = -.6612093864662645D+00
                    xi(6,1) = -.9324695142031521D+00
                    w(1) = .1713244923791703D+00
                    w(2) = .3607615730481386D+00
                    w(3) = .4679139345726910D+00
                    w(4) = .4679139345726910D+00
                    w(5) = .3607615730481386D+00
                    w(6) = .1713244923791703D+00
                    return
                case DEFAULT
                    write(6,*) ' Error in subroutine 
     +                           initialize_integration_points'
                    write(6,*) ' Invalide number of integration points
     +                           for a 1D integration scheme'
                    write(6,*) ' n_points must be between 1 and 6'
                    stop
             end select
            
          else if (n_coord==2) then                       
          ! Integration points for 2d elements
            if ( n_nodes>9) then
                write (6, 99001) n_nodes
                stop
            end if
            if ( n_points>9 ) then
                write (6, 99002) n_points
                stop
            end if
            if ( n_points==1 ) then
                if ( n_nodes==4 .or. n_nodes==9 ) then    
                !     ---   4 or 9 noded quad
                    xi(1, 1) = 0.D0
                    xi(2, 1) = 0.D0
                    w(1) = 4.D0
                else if ( n_nodes==3 .or. n_nodes==6 ) then
                 !     ---   3 or 6 noded triangle
                    xi(1, 1) = 1.D0/3.D0
                    xi(2, 1) = 1.D0/3.D0
                    w(1) = 1.D0/2.D0
                end if
            else if ( n_points==3 ) then
                xi(1, 1) = 0.5D0
                xi(2, 1) = 0.5D0
                w(1) = 1.D0/6.D0
                xi(1, 2) = 0.D0
                xi(2, 2) = 0.5D0
                w(2) = w(1)
                xi(1, 3) = 0.5D0
                xi(2, 3) = 0.D0
                w(3) = w(1)
            else if ( n_points==4 ) then
                if ( n_nodes==4 .or. n_nodes==8 .or. n_nodes==9 ) then
                    !     2X2 GAUSS INTEGRATION POINTS FOR QUADRILATERAL
                    !     43
                    !     12
                    cn = 0.5773502691896260D0
                    xi(1, 1) = -cn
                    xi(1, 2) = cn
                    xi(1, 3) = cn
                    xi(1, 4) = -cn
                    xi(2, 1) = -cn
                    xi(2, 2) = -cn
                    xi(2, 3) = cn
                    xi(2, 4) = cn
                    w(1) = 1.D0
                    w(2) = 1.D0
                    w(3) = 1.D0
                    w(4) = 1.D0
                else if ( n_nodes==3 .or. n_nodes==6 ) then
                    !     xi integration points for triangle
                    xi(1, 1) = 1.D0/3.D0
                    xi(2, 1) = xi(1, 1)
                    w(1) = -27.D0/96.D0
                    xi(1, 2) = 0.6D0
                    xi(2, 2) = 0.2D0
                    w(2) = 25.D0/96.D0
                    xi(1, 3) = 0.2D0
                    xi(2, 3) = 0.6D0
                    w(3) = w(2)
                    xi(1, 4) = 0.2D0
                    xi(2, 4) = 0.2D0
                    w(4) = w(2)
                end if

            else if ( n_points==7 ) then
                ! Quintic integration for triangle
                xi(1,1) = 1.d0/3.d0
                xi(2,1) = xi(1,1)
                w(1) = 0.1125d0
                xi(1,2) = 0.0597158717d0
                xi(2,2) = 0.4701420641d0
                w(2) = 0.0661970763d0
                xi(1,3) = xi(2,2)
                xi(2,3) = xi(1,2)
                w(3) = w(2)
                xi(1,4) = xi(2,2)
                xi(2,4) = xi(2,2)
                w(4) = w(2)
                xi(1,5) = 0.7974269853d0
                xi(2,5) = 0.1012865073d0
                w(5) = 0.0629695902d0
                xi(1,6) = xi(2,5)
                xi(2,6) = xi(1,5)
                w(6) = w(5)
                xi(1,7) = xi(2,5)
                xi(2,7) = xi(2,5)
                w(7) = w(5)
            else if ( n_points==9 ) then
                !     3X3 GAUSS INTEGRATION POINTS IN PSI-ETA COORDINATES
                !     789
                !     456
                !     123
                cn = 0.7745966692414830D0
                xi(1, 1) = -cn
                xi(1, 2) = 0.D0
                xi(1, 3) = cn
                xi(1, 4) = -cn
                xi(1, 5) = 0.D0
                xi(1, 6) = cn
                xi(1, 7) = -cn
                xi(1, 8) = 0.D0
                xi(1, 9) = cn
                xi(2, 1) = -cn
                xi(2, 2) = -cn
                xi(2, 3) = -cn
                xi(2, 4) = 0.D0
                xi(2, 5) = 0.D0
                xi(2, 6) = 0.D0
                xi(2, 7) = cn
                xi(2, 8) = cn
                xi(2, 9) = cn
                w1 = 0.5555555555555560D0
                w2 = 0.8888888888888890D0
                w11 = w1*w1
                w12 = w1*w2
                w22 = w2*w2
                w(1) = w11
                w(2) = w12
                w(3) = w11
                w(4) = w12
                w(5) = w22
                w(6) = w12
                w(7) = w11
                w(8) = w12
                w(9) = w11
            end if

99001     format ( // ' *** ERROR ***'/  
     +           '  Number of nodes on 2d element is greater ',  
     +           'than 9 '/'  No. nodes = ', I5)
99002     format ( // ' *** ERROR ***'/,  
     +           '  Number of int. pts on element is greater ',  
     +           'than 9 '/'  No. int pts. = ', I5)
  
   
        else if (n_coord==3) then                  
        ! Integration points for 3d elements
            if (n_nodes  == 4.or.n_nodes ==10 ) then
                if (n_points == 1) then
                    xi(1,1) = 0.25
                    xi(2,1) = 0.25
                    xi(3,1) = 0.25
                    w(1) = 1.D0/6.D0
                else if (n_points == 4) then
                    xi(1,1) = 0.58541020
                    xi(2,1) = 0.13819660
                    xi(3,1) = xi(2,1)
                    xi(1,2) = xi(2,1)
                    xi(2,2) = xi(1,1)
                    xi(3,2) = xi(2,1)
                    xi(1,3) = xi(2,1)
                    xi(2,3) = xi(2,1)
                    xi(3,3) = xi(1,1)
                    xi(1,4) = xi(2,1)
                    xi(2,4) = xi(2,1)
                    xi(3,4) = xi(2,1)
                    w(1:4) = 1.D0/24.D0
                else if (n_points == 5) then
                    xi(1,1) = 0.25d0
                    xi(2,1) = 0.25d0
                    xi(3,1) = 0.25d0
                    xi(1,2) = 0.5d0
                    xi(2,2) = 1.d0/6.d0
                    xi(3,2) = 1.d0/6.d0
                    xi(1,3) = 1.d0/6.d0
                    xi(2,3) = 0.5d0
                    xi(3,3) = 1.d0/6.d0
                    xi(1,4) = 1.d0/6.d0
                    xi(2,4) = 1.d0/6.d0
                    xi(3,4) = 0.5d0
                    xi(1,5) = 1.d0/6.d0
                    xi(2,5) = 1.d0/6.d0
                    xi(3,5) = 1.d0/6.d0
                    w(1) = -4.d0/30.d0
                    w(2:5) = 3.d0/40.d0
                else
                    write(6,*) ' Incorrect number of integration points
     +                            for tetrahedral element '
                    write(6, *) ' called with ',n_points
                    stop
                endif
            else if ( n_nodes == 8 .or. n_nodes == 20 ) then
                if (n_points == 1) then
                    xi(1,1) = 0.D0
                    xi(2,1) = 0.D0
                    xi(3,1) = 0.D0
                    w(1) = 8.D0
                else if (n_points == 8) then
                    x1D(1) = -0.5773502692
                    x1D(2) =  0.5773502692
                    do k = 1,2
                        do j = 1,2
                            do i = 1,2
                                n = 4*(k-1) + 2*(j-1) + i
                                xi(1,n) = x1D(i)
                                xi(2,n) = x1D(j)
                                xi(3,n) = x1D(k)
                            end do
                        end do
                    end do
                    w(1:8) = 1.D0
                else if (n_points == 27) then
                    x1D(1) = -0.7745966692
                    x1D(2) = 0.
                    x1D(3) = 0.7745966692
                    w1D(1) = 0.5555555555D0
                    w1D(2) = 0.888888888D0
                    w1D(3) = 0.55555555555D0
                    do k = 1,3
                        do j = 1,3
                            do i = 1,3
                                n = 9*(k-1) + 3*(j-1) + i
                                xi(1,n) = x1D(i)
                                xi(2,n) = x1D(j)
                                xi(3,n) = x1D(k)
                                w(n) = w1D(i)*w1D(j)*w1D(k)
                            end do
                        end do
                    end do
                else if (n_points == 64) then
                    x1D(1) = .8611363115940526D+00
                    x1D(2) = .3399810435848563D+00
                    x1D(3) = -.3399810435848563D+00
                    x1D(4) = -.8611363115940526D+00
                    w1D(1) = .3478548451374538D+00
                    w1D(2) = .6521451548625461D+00
                    w1D(3) = .6521451548625461D+00
                    w1D(4) = .3478548451374538D+00
                    do k = 1,4
                        do j = 1,4
                            do i = 1,4
                                n = 16*(k-1) + 4*(j-1) + i
                                xi(1,n) = x1D(i)
                                xi(2,n) = x1D(j)
                                xi(3,n) = x1D(k)
                                w(n) = w1D(i)*w1D(j)*w1D(k)
                            end do
                        end do
                    end do
                endif
            endif
        endif
       end subroutine initialize_integration_points
  

!==============================================================================
       subroutine calculate_shapefunctions(xi,n_nodes,n_coord,f,df)   
        implicit none
    
        integer, intent(in)   :: n_nodes
        integer, intent(in)   :: n_coord
        double precision, intent(in)  :: xi(n_coord)         
        double precision, intent(out) :: f(20)                
        double precision, intent(out) :: df(20,n_coord)      

        double precision :: z, g1, g2, g3, h1, h2, h3, dg1, dg2, dh1, 
     +                      dh2, dzdp,dzdq, xi4
    
        if (n_coord==1) then    ! 1D shape functions
    
            if (n_nodes==2) then
                f(1) = 0.5d0*(1.d0-xi(1))
                f(2) = 0.5d0*(1.d0+xi(1))
                df(1,1) = -0.5d0
                df(2,1) =  0.5d0
            else if (n_nodes==3) then
                f(1) = -0.5*xi(1)*(1.-xi(1))
                f(2) =  0.5*xi(1)*(1.+xi(1))
                f(3) = (1.-xi(1))*(1.+xi(1))
                df(1,1) = -0.5+xi(1)
                df(2,1) =  0.5+xi(1)
                df(3,1) = -2.d0*xi(1)
            endif
        else if (n_coord==2) then  
        !2D shape functions
            if ( n_nodes==3 ) then      
              !     SHAPE FUNCTIONS FOR 3 NODED TRIANGLE
                f(1) = xi(1)
                f(2) = xi(2)
                f(3) = 1.D0 - xi(1) - xi(2)
                df(1, 1) = 1.D0
                df(1, 2) = 0.D0
                df(2, 1) = 0.D0
                df(2, 2) = 1.D0
                df(3, 1) = -1.D0
                df(3, 2) = -1.D0
            else if ( n_nodes==4 ) then
                !     SHAPE FUNCTIONS FOR 4 NODED QUADRILATERAL
                !     43
                !     12
                g1 = 0.5D0*(1.D0 - xi(1))
                g2 = 0.5D0*(1.D0 + xi(1))
                h1 = 0.5D0*(1.D0 - xi(2))
                h2 = 0.5D0*(1.D0 + xi(2))
                f(1) = g1*h1
                f(2) = g2*h1
                f(3) = g2*h2
                f(4) = g1*h2
                dg1 = -0.5D0
                dg2 = 0.5D0
                dh1 = -0.5D0
                dh2 = 0.5D0
                df(1, 1) = dg1*h1
                df(2, 1) = dg2*h1
                df(3, 1) = dg2*h2
                df(4, 1) = dg1*h2
                df(1, 2) = g1*dh1
                df(2, 2) = g2*dh1
                df(3, 2) = g2*dh2
                df(4, 2) = g1*dh2
       
            else if ( n_nodes==6 ) then

                !     SHAPE FUNCTIONS FOR 6 NODED TRIANGLE
                !          3

                !       6      5

                !     1    4     2

                !     P = L1
                !     Q = L2
                !     Z = 1 - P - Q = L3

                z = 1.D0 - xi(1) - xi(2)
                f(1) = (2.D0*xi(1) - 1.D0)*xi(1)
                f(2) = (2.D0*xi(2) - 1.D0)*xi(2)
                f(3) = (2.D0*z - 1.D0)*z
                f(4) = 4.D0*xi(1)*xi(2)
                f(5) = 4.D0*xi(2)*z
                f(6) = 4.D0*xi(1)*z
                dzdp = -1.D0
                dzdq = -1.D0
                df(1, 1) = 4.D0*xi(1) - 1.D0
                df(2, 1) = 0.D0
                df(3, 1) = 4.D0*z*dzdp - dzdp
                df(4, 1) = 4.D0*xi(2)
                df(5, 1) = 4.D0*xi(2)*dzdp
                df(6, 1) = 4.D0*z + 4.D0*xi(1)*dzdp
                df(1, 2) = 0.D0
                df(2, 2) = 4.D0*xi(2) - 1.D0
                df(3, 2) = 4.D0*z*dzdq - dzdq
                df(4, 2) = 4.D0*xi(1)
                df(5, 2) = 4.D0*z + 4.D0*xi(2)*dzdq
                df(6, 2) = 4.D0*xi(1)*dzdq

            else if ( n_nodes==8 ) then
                !     SHAPE FUNCTIONS FOR 8 NODED SERENDIPITY ELEMENT
                 f(1) = -0.25*(1.-xi(1))*(1.-xi(2))*(1.+xi(1)+xi(2));
                 f(2) = 0.25*(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-1.);
                 f(3) = 0.25*(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-1.);
                 f(4) = 0.25*(1.-xi(1))*(1.+xi(2))*(xi(2)-xi(1)-1.);
                 f(5) = 0.5*(1.-xi(1)*xi(1))*(1.-xi(2));
                 f(6) = 0.5*(1.+xi(1))*(1.-xi(2)*xi(2));
                 f(7) = 0.5*(1.-xi(1)*xi(1))*(1.+xi(2));
                 f(8) = 0.5*(1.-xi(1))*(1.-xi(2)*xi(2));
                 df(1,1) = 0.25*(1.-xi(2))*(2.*xi(1)+xi(2));
                 df(1,2) = 0.25*(1.-xi(1))*(xi(1)+2.*xi(2));
                 df(2,1) = 0.25*(1.-xi(2))*(2.*xi(1)-xi(2));
                 df(2,2) = 0.25*(1.+xi(1))*(2.*xi(2)-xi(1));
                 df(3,1) = 0.25*(1.+xi(2))*(2.*xi(1)+xi(2));
                 df(3,2) = 0.25*(1.+xi(1))*(2.*xi(2)+xi(1));
                 df(4,1) = 0.25*(1.+xi(2))*(2.*xi(1)-xi(2));
                 df(4,2) = 0.25*(1.-xi(1))*(2.*xi(2)-xi(1));
                 df(5,1) = -xi(1)*(1.-xi(2));
                 df(5,2) = -0.5*(1.-xi(1)*xi(1));
                 df(6,1) = 0.5*(1.-xi(2)*xi(2));
                 df(6,2) = -(1.+xi(1))*xi(2);
                 df(7,1) = -xi(1)*(1.+xi(2));
                 df(7,2) = 0.5*(1.-xi(1)*xi(1));
                 df(8,1) = -0.5*(1.-xi(2)*xi(2));
                 df(8,2) = -(1.-xi(1))*xi(2);
            else if ( n_nodes==9 ) then
                !     SHAPE FUNCTIONS FOR 9 NODED LAGRANGIAN ELEMENT
                !     789
                !     456
                !     123
                g1 = -.5D0*xi(1)*(1.D0 - xi(1))
                g2 = (1.D0 - xi(1))*(1.D0 + xi(1))
                g3 = .5D0*xi(1)*(1.D0 + xi(1))
                h1 = -.5D0*xi(2)*(1.D0 - xi(2))
                h2 = (1.D0 - xi(2))*(1.D0 + xi(2))
                h3 = .5D0*xi(2)*(1.D0 + xi(2))
                f(1) = g1*h1
                f(2) = g2*h1
                f(3) = g3*h1
                f(4) = g1*h2
                f(5) = g2*h2
                f(6) = g3*h2
                f(7) = g1*h3
                f(8) = g2*h3
                f(9) = g3*h3
            end if
        else if (n_coord==3) then  !3D shape functions
            if (n_nodes == 4) then
                f(1) = xi(1)
                f(2) = xi(2)
                f(3) = xi(3)
                f(4) = 1.-xi(1)-xi(2)-xi(3)
                df(1,1) = 1.
                df(2,2) = 1.
                df(3,3) = 1.
                df(4,1) = -1.
                df(4,2) = -1.
                df(4,3) = -1.
            else if (n_nodes == 10) then
                xi4 = 1.-xi(1)-xi(2)-xi(3)
                f(1) = (2.*xi(1)-1.)*xi(1)
                f(2) = (2.*xi(2)-1.)*xi(2)
                f(3) = (2.*xi(3)-1.)*xi(3)
                f(4) = (2.*xi4-1.)*xi4
                f(5) = 4.*xi(1)*xi(2)
                f(6) = 4.*xi(2)*xi(3)
                f(7) = 4.*xi(3)*xi(1)
                f(8) = 4.*xi(1)*xi4
                f(9) = 4.*xi(2)*xi4
                f(10) = 4.*xi(3)*xi4
                df(1,1) = (4.*xi(1)-1.)
                df(2,2) = (4.*xi(2)-1.)
                df(3,3) = (4.*xi(3)-1.)
                df(4,1) = -(4.*xi4-1.)
                df(4,2) = -(4.*xi4-1.)
                df(4,3) = -(4.*xi4-1.)
                df(5,1) = 4.*xi(2)
                df(5,2) = 4.*xi(1)
                df(6,2) = 4.*xi(3)
                df(6,3) = 4.*xi(2)
                df(7,1) = 4.*xi(3)
                df(7,3) = 4.*xi(1)
                df(8,1) = 4.*(xi4-xi(1))
                df(8,2) = -4.*xi(1)
                df(8,3) = -4.*xi(1)
                df(9,1) = -4.*xi(2)
                df(9,2) = 4.*(xi4-xi(2))
                df(9,3) = -4.*xi(2)
                df(10,1) = -4.*xi(3)*xi4
                df(10,2) = -4.*xi(3)
                df(10,3) = 4.*(xi4-xi(3))
            else if (n_nodes == 8) then
                f(1) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3))/8.
                f(2) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3))/8.
                f(3) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3))/8.
                f(4) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3))/8.
                f(5) = (1.-xi(1))*(1.-xi(2))*(1.+xi(3))/8.
                f(6) = (1.+xi(1))*(1.-xi(2))*(1.+xi(3))/8.
                f(7) = (1.+xi(1))*(1.+xi(2))*(1.+xi(3))/8.
                f(8) = (1.-xi(1))*(1.+xi(2))*(1.+xi(3))/8.
                df(1,1) = -(1.-xi(2))*(1.-xi(3))/8.
                df(1,2) = -(1.-xi(1))*(1.-xi(3))/8.
                df(1,3) = -(1.-xi(1))*(1.-xi(2))/8.
                df(2,1) = (1.-xi(2))*(1.-xi(3))/8.
                df(2,2) = -(1.+xi(1))*(1.-xi(3))/8.
                df(2,3) = -(1.+xi(1))*(1.-xi(2))/8.
                df(3,1) = (1.+xi(2))*(1.-xi(3))/8.
                df(3,2) = (1.+xi(1))*(1.-xi(3))/8.
                df(3,3) = -(1.+xi(1))*(1.+xi(2))/8.
                df(4,1) = -(1.+xi(2))*(1.-xi(3))/8.
                df(4,2) = (1.-xi(1))*(1.-xi(3))/8.
                df(4,3) = -(1.-xi(1))*(1.+xi(2))/8.
                df(5,1) = -(1.-xi(2))*(1.+xi(3))/8.
                df(5,2) = -(1.-xi(1))*(1.+xi(3))/8.
                df(5,3) = (1.-xi(1))*(1.-xi(2))/8.
                df(6,1) = (1.-xi(2))*(1.+xi(3))/8.
                df(6,2) = -(1.+xi(1))*(1.+xi(3))/8.
                df(6,3) = (1.+xi(1))*(1.-xi(2))/8.
                df(7,1) = (1.+xi(2))*(1.+xi(3))/8.
                df(7,2) = (1.+xi(1))*(1.+xi(3))/8.
                df(7,3) = (1.+xi(1))*(1.+xi(2))/8.
                df(8,1) = -(1.+xi(2))*(1.+xi(3))/8.
                df(8,2) = (1.-xi(1))*(1.+xi(3))/8.
                df(8,3) = (1.-xi(1))*(1.+xi(2))/8.
            else if (n_nodes == 20) then
                f(1) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3))
     +                  *(-xi(1)-xi(2)-xi(3)-2.)/8.
                f(2) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3))
     +                  *(xi(1)-xi(2)-xi(3)-2.)/8.
                f(3) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3))
     +                  *(xi(1)+xi(2)-xi(3)-2.)/8.
                f(4) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3))
     +                  *(-xi(1)+xi(2)-xi(3)-2.)/8.
                f(5) = (1.-xi(1))*(1.-xi(2))*(1.+xi(3))
     +                  *(-xi(1)-xi(2)+xi(3)-2.)/8.
                f(6) = (1.+xi(1))*(1.-xi(2))*(1.+xi(3))
     +                  *(xi(1)-xi(2) +xi(3)-2.)/8.
                f(7) = (1.+xi(1))*(1.+xi(2))*(1.+xi(3))
     +                  *(xi(1)+xi(2)+xi(3)-2.)/8.
                f(8) = (1.-xi(1))*(1.+xi(2))*(1.+xi(3))
     +                  *(-xi(1)+xi(2)+xi(3)-2.)/8.
                f(9)  = (1.-xi(1)**2.)*(1.-xi(2))*(1.-xi(3))/4.
                f(10) = (1.+xi(1))*(1.-xi(2)**2.)*(1.-xi(3))/4.
                f(11) = (1.-xi(1)**2.)*(1.+xi(2))*(1.-xi(3))/4.
                f(12) = (1.-xi(1))*(1.-xi(2)**2.)*(1.-xi(3))/4.
                f(13) = (1.-xi(1)**2.)*(1.-xi(2))*(1.+xi(3))/4.
                f(14) = (1.+xi(1))*(1.-xi(2)**2.)*(1.+xi(3))/4.
                f(15) = (1.-xi(1)**2.)*(1.+xi(2))*(1.+xi(3))/4.
                f(16) = (1.-xi(1))*(1.-xi(2)**2.)*(1.+xi(3))/4.
                f(17) = (1.-xi(1))*(1.-xi(2))*(1.-xi(3)**2.)/4.
                f(18) = (1.+xi(1))*(1.-xi(2))*(1.-xi(3)**2.)/4.
                f(19) = (1.+xi(1))*(1.+xi(2))*(1.-xi(3)**2.)/4.
                f(20) = (1.-xi(1))*(1.+xi(2))*(1.-xi(3)**2.)/4.
                
               df(1,1) = (-(1.-xi(2))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)
     +                   -(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
               df(1,2) = (-(1.-xi(1))*(1.-xi(3))*(-xi(1)-xi(2)-xi(3)-2.)
     +                   -(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
               df(1,3) = (-(1.-xi(1))*(1.-xi(2))*(-xi(1)-xi(2)-xi(3)-2.)
     +                   -(1.-xi(1))*(1.-xi(2))*(1.-xi(3)))/8. 
               df(2,1) = ((1.-xi(2))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)
     +                   +(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
               df(2,2) = (-(1.+xi(1))*(1.-xi(3))*(xi(1)-xi(2)-xi(3)-2.)
     +                   -(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8.
               df(2,3) = (-(1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)-xi(3)-2.)
     +                   -(1.+xi(1))*(1.-xi(2))*(1.-xi(3)))/8. 
               df(3,1) = ((1.+xi(2))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)
     +                   +(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
               df(3,2) = ((1.+xi(1))*(1.-xi(3))*(xi(1)+xi(2)-xi(3)-2.)
     +                   +(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
               df(3,3) = (-(1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)-xi(3)-2.)
     +                   -(1.+xi(1))*(1.+xi(2))*(1.-xi(3)))/8. 
               df(4,1) = (-(1.+xi(2))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)
     +                   -(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
               df(4,2) = ((1.-xi(1))*(1.-xi(3))*(-xi(1)+xi(2)-xi(3)-2.)
     +                   +(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.
               df(4,3) = (-(1.-xi(1))*(1.+xi(2))*(-xi(1)+xi(2)-xi(3)-2.)
     +                   -(1.-xi(1))*(1.+xi(2))*(1.-xi(3)))/8.     
               df(5,1) = (-(1.-xi(2))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)
     +                   -(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
               df(5,2) = (-(1.-xi(1))*(1.+xi(3))*(-xi(1)-xi(2)+xi(3)-2.)
     +                   -(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
               df(5,3) = ((1.-xi(1))*(1.-xi(2))*(-xi(1)-xi(2)+xi(3)-2.)
     +                   +(1.-xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
               df(6,1) = ((1.-xi(2))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)
     +                   +(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
               df(6,2) = (-(1.+xi(1))*(1.+xi(3))*(xi(1)-xi(2)+xi(3)-2.)
     +                   -(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
               df(6,3) = ((1.+xi(1))*(1.-xi(2))*(xi(1)-xi(2)+xi(3)-2.)
     +                   +(1.+xi(1))*(1.-xi(2))*(1.+xi(3)))/8.
               df(7,1) = ((1.+xi(2))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)
     +                   +(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
               df(7,2) = ((1.+xi(1))*(1.+xi(3))*(xi(1)+xi(2)+xi(3)-2.)
     +                   +(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
               df(7,3) = ((1.+xi(1))*(1.+xi(2))*(xi(1)+xi(2)+xi(3)-2.)
     +                   +(1.+xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
               df(8,1) = (-(1.+xi(2))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)
     +                   -(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
               df(8,2) = ((1.-xi(1))*(1.+xi(3))*(-xi(1)+xi(2)+xi(3)-2.)
     +                   +(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
               df(8,3) = ((1.-xi(1))*(1.+xi(2))*(-xi(1)+xi(2)+xi(3) -2.)
     +                   +(1.-xi(1))*(1.+xi(2))*(1.+xi(3)))/8.
                
                df(9,1)  = -2.*xi(1)*(1.-xi(2))*(1.-xi(3))/4.
                df(9,2)  = -(1.-xi(1)**2.)*(1.-xi(3))/4.
                df(9,3)  = -(1.-xi(1)**2.)*(1.-xi(2))/4.
                df(10,1)  = (1.-xi(2)**2.)*(1.-xi(3))/4.
                df(10,2)  = -2.*xi(2)*(1.+xi(1))*(1.-xi(3))/4.
                df(10,3)  = -(1.-xi(2)**2.)*(1.+xi(1))/4.
                df(11,1)  = -2.*xi(1)*(1.-xi(2))*(1.-xi(3))/4.
                df(11,2)  = -(1.-xi(1)**2.)*(1.-xi(3))/4.
                df(11,3)  = -(1.-xi(1)**2.)*(1.-xi(2))/4.
                df(12,1)  = -(1.-xi(2)**2.)*(1.-xi(3))/4.
                df(12,2)  = -2.*xi(2)*(1.-xi(1))*(1.-xi(3))/4.
                df(12,3)  = -(1.-xi(2)**2.)*(1.-xi(1))/4.
                df(13,1)  = -2.*xi(1)*(1.-xi(2))*(1.+xi(3))/4.
                df(13,2)  = -(1.-xi(1)**2.)*(1.+xi(3))/4.
                df(13,3)  = (1.-xi(1)**2.)*(1.-xi(2))/4.
                df(14,1)  = (1.-xi(2)**2.)*(1.+xi(3))/4.
                df(14,2)  = -2.*xi(2)*(1.+xi(1))*(1.+xi(3))/4.
                df(14,3)  = (1.-xi(2)**2.)*(1.+xi(1))/4.
                df(15,1)  = 2.*xi(1)*(1.+xi(2))*(1.+xi(3))/4.
                df(15,2)  = (1.-xi(1)**2.)*(1.+xi(3))/4.
                df(15,3)  = (1.-xi(1)**2.)*(1.+xi(2))/4.
                df(16,1)  = -(1.-xi(2)**2.)*(1.+xi(3))/4.
                df(16,2)  = -2.*xi(2)*(1.-xi(1))*(1.+xi(3))/4.
                df(16,3)  = (1.-xi(2)**2.)*(1.-xi(1))/4.
                df(17,1) = -(1.-xi(2))*(1.-xi(3)**2.)/4.
                df(17,2) = -(1.-xi(1))*(1.-xi(3)**2.)/4.
                df(17,3) = -xi(3)*(1.-xi(1))*(1.-xi(2))/2.
                df(18,1) = (1.-xi(2))*(1.-xi(3)**2.)/4.
                df(18,2) = -(1.+xi(1))*(1.-xi(3)**2.)/4.
                df(18,3) = -xi(3)*(1.+xi(1))*(1.-xi(2))/2.
                df(19,1) = (1.+xi(2))*(1.-xi(3)**2.)/4.
                df(19,2) = (1.+xi(1))*(1.-xi(3)**2.)/4.
                df(19,3) = -xi(3)*(1.+xi(1))*(1.+xi(2))/2.
                df(20,1) = -(1.+xi(2))*(1.-xi(3)**2.)/4.
                df(20,2) = (1.-xi(1))*(1.-xi(3)**2.)/4.
                df(20,3) = -xi(3)*(1.-xi(1))*(1.+xi(2))/2.
            endif
   
        endif
    
      end subroutine calculate_shapefunctions


!==============================================================================
      subroutine invert_small(n_coord,A,A_inverse,determinant)
        
        implicit none

        integer, intent (in)       :: n_coord
        double precision, intent(in)    :: A(n_coord,n_coord)
        double precision, intent(out)   :: determinant
        double precision, intent(out)   :: A_inverse(n_coord,n_coord)

        double precision :: cofactor(3,3)
!
!       A (3x3) or (2x2) input matrix
!       A_inverse (3x3) or (2x2) inverse
!       determinant - determinant of matrix
!
        if (size(A)==4) then
            determinant = A(1,1)*A(2,2)-A(2,1)*A(1,2)
            A_inverse(1,1) = A(2,2)
            A_inverse(2,2) = A(1,1)
            A_inverse(1,2) = -A(1,2)
            A_inverse(2,1) = -A(2,1)
            IF (determinant==0.d0) THEN
                write(6,*) ' Error in element utility invert'
                write(6,*) ' A 2x2 matrix has a zero determinant'
                stop
            endif
            A_inverse = A_inverse/determinant
            
        else if (size(A)==9) then
            determinant = A(1,1)*A(2,2)*A(3,3)  
     +           - A(1,1)*A(2,3)*A(3,2)  
     +           - A(1,2)*A(2,1)*A(3,3)  
     +           + A(1,2)*A(2,3)*A(3,1)  
     +           + A(1,3)*A(2,1)*A(3,2)  
     +           - A(1,3)*A(2,2)*A(3,1)

            IF (determinant==0.d0) THEN
                write(6,*) ' Error in element utility invert'
                write(6,*) ' A 3x3 matrix has a zero determinant'
                stop
            endif
            COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
            COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
            COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
            COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
            COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
            COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
            COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
            COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
            COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

            A_inverse = transpose(COFACTOR) / determinant
        endif

      end subroutine invert_small
