      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
      INCLUDE 'ABA_PARAM.INC'
!     WARNING - the aba_param.inc file declares
!        Implicit real*8(a-h,o-z)
!     This means that, by default, any variables with
!     first letter between a-h or o-z are double precision.
!     The rest are integers.
!     Note that this also means that if you type a variable
!     name incorrectly, the compiler won't catch your typo.
!
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)


!
!      DDSDDE(NTENS,NTENS)
!         Jacobian matrix of the constitutive model.
!         DDSDDE(I,J) defines the change in the Ith stress component
!         at the end of the time increment caused by an infinitesimal
!         perturbation of the Jth component of the strain increment array.
!         Unless you invoke the unsymmetric equation solution capability
!         for the user-defined material, ABAQUS/Standard will use only
!         the symmetric part of DDSDDE. The symmetric part of the matrix
!         is calculated by taking one half the sum of the matrix and its transpose.


!      STRESS(NTENS)
!         This array is passed in as the stress tensor at the beginning
!         of the increment and must be updated in this routine to be the
!         stress tensor at the end of the increment. If you specified
!         initial stresses (“Initial conditions,” Section 19.2.1), this
!         array will contain the initial stresses at the start of the
!         analysis. The size of this array depends on the value of NTENS
!         as defined below. In finite-strain problems the stress tensor
!         has already been rotated to account for rigid body motion in
!         the increment before UMAT is called, so that only the corotational
!         part of the stress integration should be done in UMAT. The
!         measure of stress used is “true” (Cauchy) stress.
!
!      STATEV(NSTATV)
!         An array containing the solution-dependent state variables.
!         These are passed in as the values at the beginning of the
!         increment unless they are updated in user subroutines USDFLD
!        (“USDFLD,” Section 25.2.39) or UEXPAN (“UEXPAN,” Section 25.2.20),
!        in which case the updated values are passed in. In all cases
!         STATEV must be returned as the values at the end of the increment.
!         The size of the array is defined as described in 
!        “Allocating space” in “User subroutines: overview,” Section 25.1.1.
!
!         In finite-strain problems any vector-valued or tensor-valued
!         state variables must be rotated to account for rigid body
!         motion of the material, in addition to any update in the
!         values associated with constitutive behavior. The rotation
!         increment matrix, DROT, is provided for this purpose.
!
!      SSE, SPD, SCD
!         Specific elastic strain energy, plastic dissipation, and
!         “creep” dissipation, respectively. These are passed in as
!         the values at the start of the increment and should be
!         updated to the corresponding specific energy values at
!         the end of the increment. They have no effect on the solution,
!         except that they are used for energy output.
!
!     Only in a fully coupled thermal-stress analysis
!      RPL
!         Volumetric heat generation per unit time at the end of the increment
!         caused by mechanical working of the material.
!
!     DDSDDT(NTENS)
!          Variation of the stress increments with respect to the temperature.
!
!     DRPLDE(NTENS)
!           Variation of RPL with respect to the strain increments.
!
!     DRPLDT
!           Variation of RPL with respect to the temperature.
!
!     Variables that can be updated
!
!     PNEWDT
!        Ratio of suggested new time increment to the time increment being
!        used (DTIME, see discussion later in this section). This variable
!        allows you to provide input to the automatic time incrementation
!        algorithms in ABAQUS/Standard (if automatic time incrementation is chosen).
!        For a quasi-static procedure the automatic time stepping that ABAQUS/Standard
!        uses, which is based on techniques for integrating standard creep laws
!        (see “Quasi-static analysis,” Section 6.2.5), cannot be controlled from within
!        the UMAT subroutine.
!        PNEWDT is set to a large value before each call to UMAT.
!        If PNEWDT is redefined to be less than 1.0, ABAQUS/Standard must abandon the
!        time increment and attempt it again with a smaller time increment. The
!        suggested new time increment provided to the automatic time integration
!        algorithms is PNEWDT × DTIME, where the PNEWDT used is the minimum value
!        for all calls to user subroutines that allow redefinition of PNEWDT for this
!        iteration.
!        If PNEWDT is given a value that is greater than 1.0 for all calls to user
!        subroutines for this iteration and the increment converges in this iteration,
!        ABAQUS/Standard may increase the time increment. The suggested new time increment
!        provided to the automatic time integration algorithms is PNEWDT × DTIME, where
!        the PNEWDT used is the minimum value for all calls to user subroutines for
!        this iteration.
!        If automatic time incrementation is not selected in the analysis procedure,
!        values of PNEWDT that are greater than 1.0 will be ignored and values of
!        PNEWDT that are less than 1.0 will cause the job to terminate.
!
!    Variables passed in for information
!
!     STRAN(NTENS)
!         An array containing the total strains at the beginning of the increment.
!         If thermal expansion is included in the same material definition, the
!         strains passed into UMAT are the mechanical strains only (that is, the
!         thermal strains computed based upon the thermal expansion coefficient have
!         been subtracted from the total strains). These strains are available for output
!         as the “elastic” strains.
!
!         In finite-strain problems the strain components have been rotated to account for
!         rigid body motion in the increment before UMAT is called and are approximations
!         to logarithmic strain.

!     DSTRAN(NTENS)
!         Array of strain increments. If thermal expansion is included in the same
!         material definition, these are the mechanical strain increments (the total
!         strain increments minus the thermal strain increments).
!
!     TIME(1)
!         Value of step time at the beginning of the current increment.
!
!     TIME(2)
!          Value of total time at the beginning of the current increment.
!
!     DTIME
!        Time increment.
!
!     TEMP
!         Temperature at the start of the increment.
!
!     DTEMP
!         Increment of temperature.
!
!     PREDEF
!        Array of interpolated values of predefined field variables at this point
!        at the start of the increment, based on the values read in at the nodes.
!
!      DPRED
!        Array of increments of predefined field variables.
!
!      CMNAME
!        User-defined material name, left justified. Some internal material models are given names starting with the “ABQ_” character string. To avoid conflict, you should not use “ABQ_” as the leading string for CMNAME.
!
!      NDI
!        Number of direct stress components at this point.
!
!      NSHR
!        Number of engineering shear stress components at this point.
!
!      NTENS
!        Size of the stress or strain component array (NDI + NSHR).
!
!      NSTATV
!         Number of solution-dependent state variables that are associated with
!         this material type (defined as described in “Allocating space” in “User
!         subroutines: overview,” Section 25.1.1).
!
!      PROPS(NPROPS)
!         User-specified array of material constants associated with this user material.
!
!      NPROPS
!         User-defined number of material constants associated with this user material.
!
!      COORDS
!         An array containing the coordinates of this point. These are the current
!         coordinates if geometric nonlinearity is accounted for during the step
!         (see “Procedures: overview,” Section 6.1.1); otherwise, the array contains
!         the original coordinates of the point.
!
!     DROT(3,3)
!          Rotation increment matrix. This matrix represents the increment of rigid
!          body rotation of the basis system in which the components of stress
!          (STRESS) and strain (STRAN) are stored. It is provided so that vector- or
!          tensor-valued state variables can be rotated appropriately in this subroutine:
!          stress and strain components are already rotated by this amount before UMAT
!          is called. This matrix is passed in as a unit matrix for small-displacement
!          analysis and for large-displacement analysis if the basis system for the
!          material point rotates with the material (as in a shell element or when a
!          local orientation is used).
!
!      CELENT
!          Characteristic element length, which is a typical length of a line across
!          an element for a first-order element; it is half of the same typical length
!          for a second-order element. For beams and trusses it is a characteristic length
!          along the element axis. For membranes and shells it is a characteristic length
!          in the reference surface. For axisymmetric elements it is a characteristic length
!          in the  plane only. For cohesive elements it is equal to the constitutive
!          thickness.
!
!      DFGRD0(3,3)
!          Array containing the deformation gradient at the beginning of the increment.
!          See the discussion regarding the availability of the deformation gradient for
!          various element types.
!
!     DFGRD1(3,3)
!            Array containing the deformation gradient at the end of the increment.
!           The components of this array are set to zero if nonlinear geometric effects
!           are not included in the step definition associated with this increment. See
!           the discussion regarding the availability of the deformation gradient for
!           various element types.
!
!      NOEL
!           Element number.
!
!      NPT
!           Integration point number.
!
!      LAYER
!          Layer number (for composite shells and layered solids).
!
!      KSPT
!          Section point number within the current layer.
!
!      KSTEP
!         Step number.
!
!     KINC
!         Increment number.

!      user coding to define DDSDDE, STRESS, STATEV, SSE, SPD, SCD
!      and, if necessary, RPL, DDSDDT, DRPLDE, DRPLDT, PNEWDT
!
!     Local variables

  double precision E,xnu
	integer i,j

  ddsdde = 0.d0

  E = props(1)
	xnu = props(2)

!    for debugging, you can use
!      write(6,*) ' Hello '
!    Output is then written to the .dat file
	
	If (ndi==3 .and. nshr==1) then    ! Plane strain or axisymmetry
	  ddsdde(1,1) = 1.d0-xnu
	  ddsdde(1,2) = xnu
	  ddsdde(1,3) = xnu
	  ddsdde(2,1) = xnu
	  ddsdde(2,2) = 1.d0-xnu
	  ddsdde(2,3) = xnu
	  ddsdde(3,1) = xnu
	  ddsdde(3,2) = xnu
	  ddsdde(3,3) = 1.d0-xnu
	  ddsdde(4,4) = 0.5d0*(1.d0-2.d0*xnu)
	  ddsdde = ddsdde*E/( (1.d0+xnu)*(1.d0-2.d0*xnu) )
	else if (ndi==2 .and. nshr==1) then   ! Plane stress
	  ddsdde(1,1) = 1.d0
	  ddsdde(1,2) = xnu
	  ddsdde(2,1) = xnu
	  ddsdde(2,2) = 1.d0
	  ddsdde(3,3) = 0.5d0*(1.d0-xnu)
	  ddsdde = ddsdde*E/( (1.d0+xnu*xnu) )
  else ! 3D
	  ddsdde(1,1) = 1.d0-xnu
	  ddsdde(1,2) = xnu
	  ddsdde(1,3) = xnu
	  ddsdde(2,1) = xnu
	  ddsdde(2,2) = 1.d0-xnu
	  ddsdde(2,3) = xnu
	  ddsdde(3,1) = xnu
	  ddsdde(3,2) = xnu
	  ddsdde(3,3) = 1.d0-xnu
	  ddsdde(4,4) = 0.5d0*(1.d0-2.d0*xnu)
	  ddsdde(5,5) = ddsdde(4,4)
	  ddsdde(6,6) = ddsdde(4,4)
	  ddsdde = ddsdde*E/( (1.d0+xnu)*(1.d0-2.d0*xnu) )
  endif
!
!     NOTE: ABAQUS uses engineering shear strains,
!     i.e. stran(ndi+1) = 2*e_12, etc...
  do i = 1,ntens
	  do j = 1,ntens
	     stress(i) = stress(i) + ddsdde(i,j)*dstran(j)
	  end do
	end do

      RETURN
      END
