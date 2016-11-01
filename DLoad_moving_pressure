! http://www.eng-tips.com/viewthread.cfm?qid=390569
!-------------------------------------BOF---------------------------------------

      subroutine dload(f,kstep,kinc,time,noel,npt,layer,kspt,
     *                 coords,jltyp,sname)
!
      include 'aba_param.inc'
!
      dimension time(2), coords (3)
      character*80 sname

      !---------------------------------------------------------------72--------

      ! user parameters
      parameter(rPressure = 1.d0, ! pressure value
     *          rRadius   = 5.d0) ! ball radius

      ! index
      parameter(iX    = 1, ! x-coord for coords array
     *          iY    = 2, ! y-coord for coords array
     *          iZ    = 3, ! z-coord for coords array
     *          iStepTime = 1, ! step time
     *          iTotTime  = 2) ! total time

      rX0   =   0.d0 ! x-coord ball start point
      rY0   =  50.d0 ! y-coord ball start point
      rXVel = 100.d0 ! ball velocity in x direction
      rYVel =   0.d0 ! ball velocity in y direction

      ! ball current position (s=x0+V*t)
      rXBall = rX0 + (rXVel * time(iTotTime))
      rYBall = rY0 + (rYVel * time(iTotTime))

      ! current element position respect to center of the ball
      ! circle equation: (x-a)^2+(x-b)^2=r^2
      rPosition = sqrt((coords(iX)-rXBall)**2 + (coords(iY)-rYBall)**2)

      ! are you under ball?
      ! ... yes I am
      if (rPosition <= rRadius) then
        f = 1.d0 ! apply pressure
      ! ... no I am not
      else
        f = 0.d0 ! do not apply pressure
      endif

      !---------------------------------------------------------------72--------
      return
      end

!-------------------------------------EOF---------------------------------------

