module HDS
    USE type_HDS

    contains

    subroutine init_pond_Area_Volume(depArea, depVol, totEvap, volFrac, p, pondVol, pondArea)
        ! used to initialize pond volume and area
        ! initialize at capacity minus depth of evaporation
        implicit none
        !subroutine arguments
        real(rkind),  intent(in)  :: depArea                    ! depression area [m2]
        real(rkind),  intent(in)  :: depVol                     ! depression volume [m3]
        real(rkind),  intent(in)  :: totEvap                    ! total evaporation [m] to be subtracted from the ponded water
        real(rkind),  intent(in)  :: volFrac                    ! volume fraction [-]. Used to fill the depressions using percentage (i.e., 50% full)
        real(rkind),  intent(in)  :: p                          ! shape of the slope profile [-]. Exponent for calculating the fractional wet area
        real(rkind),  intent(out) :: pondVol, pondArea          ! pond volume [m3], pond area [m2]
        ! local variables
        real(rkind)               :: depHeight  ! depression height [m]
        real(rkind)               :: pondHeight ! height of the water level in depression [m]

        ! estimate the maximum depth of the depression
        depHeight = depVol*(one + two/p)/depArea  ! depression height [m]

        !subtract ET from ponded water if ET > 0
        if(totEvap > zero) then

            ! estimate the height of the water level after evaporation
            pondHeight = max(depHeight - totEvap, zero)

            ! get the volume of water in the depression
            pondVol = depVol*((pondHeight/depHeight)**(one  + two/p))

        else !use volume frac to get initial pondVol
            pondVol = depVol*volFrac
        end if

        ! calculate pondArea
        pondArea = depArea*((pondVol/depVol)**(two/(p + two)))

    end subroutine init_pond_Area_Volume

    !=============================================================
    !=============================================================
    subroutine runDepression(pondVol,                  &       ! input/output:  state variable = pond volume [m3]
                             qSeas, pRate, etPond,     &       ! input:         forcing data = runoff, precipitation, ET [mm/day]
                             depArea, depVol, upsArea, &       ! input:         spatial attributes = depression area [m2], depression volume [m3], upstream area [m2]
                             p, tau,                   &       ! input:         model parameters = p [-] shape of the slope profile; tau [day-1] time constant linear reservoir
                             b, vmin,                  &       ! input:         model parameters = b [-] shape of contributing fraction curve; vmin [m3] minimum volume
                             dt,                       &       ! input:         model time step [days]
                             Q_det_adj, Q_dix_adj,     &       ! output:        adjusted evapotranspiration & infiltration fluxes [L3 T-1] for mass balance closure (i.e., when losses > pondVol)
                             fVol, fArea,              &       ! output:        fractional volume [-], fractional contributing area [-]
                             pondArea, pondOutflow)            ! output:        pond area at the end of the time step [m2], pond outflow [m3]


        ! Estimate the volumetric storage at the end of the time interval using the numerical solution
        implicit none
        ! subroutine arguments
        real(rkind),  intent(inout) :: pondVol, vmin               ! state variable = pond volume [m3], vmin [m3] minimum volume
        real(rkind),  intent(in)    :: qSeas, pRate, etPond        ! input:         forcing data = runoff, precipitation, ET [mm/day]
        real(rkind),  intent(in)    :: depArea, depVol, upsArea    ! input:         spatial attributes = depression area [m2], depression volume [m3], upstream area [m2]
        real(rkind),  intent(in)    :: p, tau                      ! input:         model parameters = p [-] shape of the slope profile; tau [day-1] time constant linear reservoir
        real(rkind),  intent(in)    :: b                           ! input:         model parameters = b [-] shape of contributing fraction curve; vminold [m3] minimum volume
        real(rkind),  intent(in)    :: dt                          ! input:         model time step [days]
        real(rkind),  intent(out)   :: Q_det_adj, Q_dix_adj        ! output:        adjusted evapotranspiration & infiltration fluxes [L3 T-1] for mass balance closure (i.e., when losses > pondVol). Zero values mean no adjustment needed.
        real(rkind),  intent(out)   :: fVol, fArea                 ! output:        fractional volume [-], fractional contributing area [-]
        real(rkind),  intent(out)   :: pondArea, pondOutflow       ! output:        pond area at the end of the time step [m2]
        ! local variables -- model decisions
        integer(i4b), parameter     :: implicitEuler=1001          ! named variable for the implicit Euler solution
        integer(i4b), parameter     :: shortSubsteps=1002          ! named variable for the short substeps solution
        integer(i4b)                :: solution                    ! numerical solution method
        ! local variables --- implicit Euler solution
        integer(i4b)                :: iter                        ! iteration counter
        integer(i4b) , parameter    :: nIter=100                   ! number of iterations (algorithm control parameter)
        real(rkind)  , parameter    :: xConv=1.e-6_rkind           ! convergence criteria (algorithm control parameter)
        real(rkind)                 :: xMin, xMax                  ! min and max storage
        real(rkind)                 :: xRes                        ! residual
        logical(lgt)                :: failure                     ! failure flag
        ! local variables -- short substeps solution
        integer(i4b) , parameter    :: nSub=100                    ! number of substeps (algorithm control parameter)
        integer(i4b)                :: iSub                        ! counter for shortSubsteps solution
        real(rkind)                 :: dtSub                       ! dt for shortSubsteps solution
        ! local variables -- diagnostic variables (model physics)
        real(rkind)                 ::  Q_di, Q_det, Q_dix, Q_do   ! fluxes [L3 T-1]: sum of water inputs to the pond, evapotranspiration, infiltration, pond outflow
        real(rkind)                 ::  cFrac, g, dgdv             ! contributing fraction, net fluxes and derivative
        real(rkind)                 ::  xVol                       ! pond volume at the end of the time step [m3]

        ! define numerical solution
        solution      = implicitEuler
        !solution      = shortSubsteps

        ! initialize pond volume
        xVol = pondVol
        ! initialize corrected evaporation and infiltration for mass balance closure
        Q_det_adj = zero
        Q_dix_adj = zero

        ! ---------- option 1: implicit Euler ----------

        ! check if implicit Euler is desired
        if(solution == implicitEuler)then

            ! initialize brackets
            xMin = zero
            xMax = depVol

            ! iterate (can start with 1 because the index is not used)
            do iter=1,nIter

                ! compute fluxes and derivatives
                call computFlux(iter, xVol, qSeas, pRate, etPond, depArea, depVol, upsArea, p, tau, b, vMin, Q_di, Q_det, Q_dix, Q_do, cFrac, g, dgdv)

                ! compute residual
                xRes = (pondVol + g*dt) - xVol

                ! check convergence
                if(iter > 1 .and. abs(xRes) < xConv)then
                    failure=.false.
                    exit
                endif

                ! update constraints
                if(xRes > zero) xMin=xVol
                if(xRes < zero) xMax=xVol

                ! special case where xMax is too small
                if(xRes > zero  .and. xVol > 0.99*xMax) xMax = min(xMax*10.0_rkind, depVol)

                ! update state (pondVol)
                xVol = xVol + xRes / (one  - dgdv*dt)

                ! use bi-section if violated constraints
                if(xVol < xMin .or. xVol > xMax) xVol=(xMin+xMax)/2.0_rkind

                ! assign failure
                failure = (iter == nIter)

            enddo ! iterating

        endif  ! if implicit Euler

        ! ---------- option 2: short substeps ----------

        ! check if short substeps are desired
        if(solution == shortSubsteps .or. failure)then

            ! define length of the substeps
            dtSub = one  / real(nSub, kind(rkind))

            ! loop through substeps
            do iSub=1,nSub
                call computFlux(iter, xVol, qSeas, pRate, etPond, depArea, depVol, upsArea, p, tau, b, vMin, Q_di, Q_det, Q_dix, Q_do, cFrac, g, dgdv)
                xVol = xVol + g*dtSub
                ! if xVol is -ve, the code produces NaNs
                ! prevent -ve ponVol values to avoid nans
                if(xVol < zero)then
                    xVol = zero
                    ! adjust evaporation and infiltration fluxes proprtionally (using weights) to close mass balance (i.e., losses = PondVol)
                    Q_det_adj = pondVol * (Q_det/(Q_det+Q_dix))
                    Q_dix_adj = pondVol * (Q_dix/(Q_det+Q_dix))
                    !break the loop as there's no need to continue xvol<zero
                    exit
                end if
            enddo  ! looping through substeps

        endif  ! if short substeps

        ! ---------------------------------------------------------------------------------
        ! ---------------------------------------------------------------------------------
        ! ---------- fill and spill process for the meta depression model -----------------
        fArea = cFrac
        ! ---------------------------------------------------------------------------------
        
        ! save variables
        pondVol = xVol
        ! compute fractional volume and fractional area
        fVol  = pondVol/depVol
        ! compute the pond area at the end of the time step
        pondArea = depArea*((pondVol/depVol)**(two/(p + two)))
        ! pond outflow [m3]
        pondOutflow = Q_do
    end subroutine runDepression

    !=============================================================
    !=============================================================
    subroutine computFlux(iter, pondVol,             & ! input:  iteration index, state variable = pond volume [m3]
                          qSeas, pRate, etPond,      & ! input:  forcing data = runoff, precipitation, ET [mm/day]
                          depArea, depVol, upsArea,  & ! input:  spatial attributes = depression area [m2], depression volume [m3], upstream area [m2]
                          p, tau,                    & ! input:  model parameters = p [-] shape of the slope profile; tau [day-1] time constant linear reservoir
                          b, vmin,                   & ! input:  model parameters = b [-] shape of contributing fraction curve; vmin [m3] minimum volume
                          Q_di, Q_det, Q_dix, Q_do,  & ! output: individual model fluxes
                          cFrac, g, dgdv)              ! output: contributing fraction, net fluxes and derivative
        implicit none
        ! subroutine arguments
        integer, intent(in)         :: iter                        ! iteration index
        real(rkind),  intent(in)    :: pondVol                     ! state variable = pond volume [m3]
        real(rkind),  intent(in)    :: qSeas, pRate, etPond        ! input:  forcing data = runoff, precipitation, ET [mm/day]
        real(rkind),  intent(in)    :: depArea, depVol, upsArea    ! input:  spatial attributes = depression area [m2], depression volume [m3], upstream area [m2]
        real(rkind),  intent(in)    :: p, tau                      ! input:  model parameters = p [-] shape of the slope profile; tau [day-1] time constant linear reservoir
        real(rkind),  intent(in)    :: b                           ! input:  model parameters = b [-] shape of contributing fraction curve; vmin [m3] minimum volume
        real(rkind),  intent(inout) :: vmin                        ! in/out: vmin [m3] minimum volume
        real(rkind),  intent(out)   :: Q_di, Q_det, Q_dix, Q_do    ! output: individual model fluxes
        real(rkind),  intent(out)   :: cFrac, g, dgdv              ! output: contributing fraction, net fluxes and derivative, pond area
        ! local variables
        real(rkind), parameter      ::  ms=0.0001_rkind            ! smoothing parameter (algorithm control)
        real(rkind), parameter      ::  rCoef=zero                 ! runoff coefficient [-] !HDS_standalone 0.050_rkind
        real(rkind)                 ::  pondArea                   ! calculated pond area
        real(rkind)                 ::  pInput                     ! precipitation (or rain+melt) [m/day]
        real(rkind)                 ::  qInput                     ! surface runoff [m/day]
        real(rkind)                 ::  eLosses                    ! evaporation losses [m/day]
        real(rkind)                 ::  vPrime                     ! smoothed pondVol value for Euler solution
        real(rkind)                 ::  vTry                       ! adjusted pondVol for derivatives calculations
        real(rkind)                 ::  dadv, dpdv, dfdv, didv     ! derivatives

        ! compute the pond area
        pondArea = depArea*((pondVol/depVol)**(two/(p + two)))

        ! get the forcing
        pInput  = pRate /rho_w  ! kg m-2 -> m s-1
        qInput  = qSeas /rho_w  + rCoef*pRate/rho_w  ! surface runoff ! kg m-2 s-1 -> m s-1
        eLosses = etPond/rho_w  ! evaporation losses kg m-2 s-1 -> m s-1

        ! get volume fluxes from the host land model
        ! sum of water input to the depression (eq 11)
        Q_di  = upsArea*qInput + (depArea - pondArea)*qInput + pondArea*pInput

        ! evapotranspiration losses (eq 12)
        Q_det = pondArea*eLosses

        ! compute infiltration from the bottom of the pond (eq 13)
        Q_dix = tau*pondVol

        ! update the minimum parameter (force hysteresis)
        if(iter == 1) then
            if(Q_di < Q_det+Q_dix) vmin = pondVol
        endif

        ! compute the outflow from the meta depression
        ! smoothing pondVol calculation ! (eq 24)
        vPrime = half*(vMin + pondVol + sqrt((pondVol-vMin)**two  + ms)) ! (eq 24)
        vPrime = min(vPrime, depVol) !MIA -> limit vPrime to be <= depVol
        ! calculate contributing fraction !(eq 25)
        cFrac  = one  - ((depVol - vPrime)/(depVol - vMin))**b !(eq 25)
        cFrac = max(cFrac, zero ) !prevent -ve values MIA
        ! calculate outlfow (eq 16)
        Q_do   = cFrac*Q_di

        ! compute the net flux (eq 10)
        g = Q_di - Q_det - Q_dix - Q_do

        ! compute derivate in pond area w.r.t. pond volume
        vTry = max(verySmall, pondVol)  ! avoid divide by zero
        dadv = ((two *depArea)/((p + two )*depVol)) * ((vTry/depVol)**(-p/(p + two )))

        ! compute derivatives in volume fluxes w.r.t. pond volume (meta depression)
        dpdv = half*((pondVol-vMin)/sqrt((pondVol-vMin)**two  + ms) + one )  ! max smoother
        dfdv = dpdv*(b*(depVol - pondVol)**(b - one ))/((depVol - vMin)**b)  ! contributing fraction (note chain rule)
        didv = dadv*(pInput - qInput)
        dgdv = didv*(one  - cFrac) - dfdv*Q_di - dadv*eLosses - tau

    end subroutine computFlux

end module HDS
