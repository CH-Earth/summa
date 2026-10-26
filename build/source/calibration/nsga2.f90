! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.

! **************************************************************************************************
! NSGA-II operators (Deb, Pratap, Agarwal and Meyarivan, 2002).
!
! Model-agnostic, like parameter_search. Objectives arrive as f(nObjective,nMember) oriented so that
! smaller is better. Decision vectors are physical parameter values; the variation operators work in
! the transformed search space of the parameter specification, scaled to [0,1] per variable.
! **************************************************************************************************
module nsga2

  USE nr_type, only: i4b, rkind, lgt
  USE, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf

  USE parameter_search, only: parameter_search_info
  USE parameter_search, only: transform_parameter
  USE parameter_search, only: inverse_transform_parameter
  USE parameter_search, only: check_ordered_constraints

  implicit none
  private

  public :: nondominated_sort
  public :: nondominated_set
  public :: crowding_distance
  public :: select_survivors
  public :: make_offspring

contains

  ! **************************************************************************************************
  ! Fast non-dominated sort: rank(i) is the front member i belongs to, 1 for the non-dominated set.
  ! **************************************************************************************************
  subroutine nondominated_sort(f,rank)
    implicit none
    real(rkind),  intent(in)  :: f(:,:)          ! objectives (nObjective,nMember), smaller is better
    integer(i4b), intent(out) :: rank(:)         ! front of each member
    logical(lgt), allocatable :: dom(:,:)        ! dom(i,j) is true when i dominates j
    integer(i4b), allocatable :: nDominating(:)  ! members not yet ranked that dominate each member
    logical(lgt), allocatable :: current(:)      ! members of the front being ranked
    integer(i4b) :: i,j,n,iFront

    n=size(f,2)
    allocate(dom(n,n),nDominating(n),current(n))
    do j=1,n
      do i=1,n
        dom(i,j)=dominates(f(:,i),f(:,j))
      enddo
    enddo
    nDominating=count(dom,dim=1)

    rank=0
    iFront=0
    do while(any(rank==0))
      iFront=iFront+1
      current=(rank==0 .and. nDominating==0)
      where(current) rank=iFront
      do i=1,n
        if(current(i)) where(dom(i,:)) nDominating=nDominating-1
      enddo
    enddo

  end subroutine nondominated_sort

  ! **************************************************************************************************
  ! Flag the members no other member dominates, without the n-by-n table a full sort keeps.
  ! **************************************************************************************************
  function nondominated_set(f) result(flag)
    implicit none
    real(rkind),  intent(in)  :: f(:,:)          ! objectives (nObjective,nMember), smaller is better
    logical(lgt)              :: flag(size(f,2))
    integer(i4b) :: i,j

    flag=.true.
    do i=1,size(f,2)
      do j=1,size(f,2)
        if(dominates(f(:,j),f(:,i)))then
          flag(i)=.false.
          exit
        endif
      enddo
    enddo

  end function nondominated_set

  ! **************************************************************************************************
  ! Crowding distance of each member within its front: the sum over objectives of the normalized gap
  ! between its neighbours. A front's extremes, and every member of a front of one or two, are infinite.
  ! **************************************************************************************************
  subroutine crowding_distance(f,rank,distance)
    implicit none
    real(rkind),  intent(in)  :: f(:,:)          ! objectives (nObjective,nMember), smaller is better
    integer(i4b), intent(in)  :: rank(:)         ! front of each member
    real(rkind),  intent(out) :: distance(:)     ! crowding distance of each member
    integer(i4b), allocatable :: front(:)        ! members of one front
    integer(i4b), allocatable :: order(:)        ! those members sorted on one objective
    real(rkind)  :: span                         ! range of one objective over the front
    real(rkind)  :: infinity
    integer(i4b) :: i,j,k,m,iFront

    infinity=ieee_value(1._rkind,ieee_positive_inf)
    distance=0._rkind
    do iFront=1,maxval(rank)
      front=pack([(i,i=1,size(rank))],rank==iFront)
      m=size(front)
      if(m <= 2)then
        distance(front)=infinity
        cycle
      endif
      do k=1,size(f,1)
        order=front(sort_index(f(k,front)))
        distance(order(1))=infinity
        distance(order(m))=infinity
        span=f(k,order(m))-f(k,order(1))
        if(span <= 0._rkind) cycle
        do j=2,m-1
          distance(order(j))=distance(order(j))+(f(k,order(j+1))-f(k,order(j-1)))/span
        enddo
      enddo
    enddo

  end subroutine crowding_distance

  ! **************************************************************************************************
  ! Elitist survivor selection: fill nKeep places front by front, taking the least crowded members of
  ! the front that does not fit. Returns the survivors' positions in f with their rank and distance.
  ! **************************************************************************************************
  subroutine select_survivors(f,nKeep,keep,rank,distance)
    implicit none
    real(rkind),  intent(in)  :: f(:,:)          ! objectives of parents and offspring together
    integer(i4b), intent(in)  :: nKeep           ! population size
    integer(i4b), intent(out) :: keep(:)         ! positions in f of the survivors (nKeep)
    integer(i4b), intent(out) :: rank(:)         ! their fronts (nKeep)
    real(rkind),  intent(out) :: distance(:)     ! their crowding distances (nKeep)
    integer(i4b), allocatable :: rankAll(:)
    real(rkind),  allocatable :: distanceAll(:)
    integer(i4b), allocatable :: front(:)
    integer(i4b), allocatable :: order(:)
    integer(i4b) :: i,nKept,nTake,iFront

    allocate(rankAll(size(f,2)),distanceAll(size(f,2)))
    call nondominated_sort(f,rankAll)
    call crowding_distance(f,rankAll,distanceAll)

    nKept=0
    iFront=0
    do while(nKept < nKeep)
      iFront=iFront+1
      front=pack([(i,i=1,size(rankAll))],rankAll==iFront)
      nTake=min(size(front),nKeep-nKept)
      if(nTake < size(front))then
        order=sort_index(-distanceAll(front))
        front=front(order(1:nTake))
      endif
      keep(nKept+1:nKept+nTake)=front
      nKept=nKept+nTake
    enddo
    rank=rankAll(keep)
    distance=distanceAll(keep)

  end subroutine select_survivors

  ! **************************************************************************************************
  ! Breed one generation of offspring from the parents: binary crowded tournaments pick each pair,
  ! simulated binary crossover (SBX) recombines it, and polynomial mutation perturbs each child. A
  ! child violating an ordered constraint is discarded and bred again.
  ! **************************************************************************************************
  subroutine make_offspring(search,parents,rank,distance,pc,eta_c,pm,eta_m,children,err,message)
    implicit none
    type(parameter_search_info), intent(in)  :: search          ! parameter-search information
    real(rkind),                 intent(in)  :: parents(:,:)    ! physical parameter values (nParam,nParent)
    integer(i4b),                intent(in)  :: rank(:)         ! front of each parent
    real(rkind),                 intent(in)  :: distance(:)     ! crowding distance of each parent
    real(rkind),                 intent(in)  :: pc              ! crossover probability per pair
    real(rkind),                 intent(in)  :: eta_c           ! SBX distribution index
    real(rkind),                 intent(in)  :: pm              ! mutation probability per variable
    real(rkind),                 intent(in)  :: eta_m           ! mutation distribution index
    real(rkind),                 intent(out) :: children(:,:)   ! physical parameter values (nParam,nChild)
    integer(i4b),                intent(out) :: err
    character(*),                intent(out) :: message
    real(rkind),  allocatable :: z(:,:)          ! parents in scaled search coordinates
    real(rkind),  allocatable :: width(:)        ! search-space width of each variable
    logical(lgt), allocatable :: active(:)       ! variables with room to vary
    real(rkind)  :: c(size(parents,1),2)         ! one pair of children, scaled
    real(rkind)  :: x(size(parents,1))           ! one child, physical
    integer(i4b) :: nParam,nChild,iChild,iPair,iTry,i,d,p1,p2
    integer(i4b), parameter :: maxtry=10000      ! pairs bred per child before giving up
    real(rkind)  :: s
    character(len=256) :: cmessage

    err=0
    message='make_offspring/'
    nParam=size(parents,1)
    nChild=size(children,2)
    if(nParam /= size(search%param_names) .or. size(children,1) /= nParam)then
      message=trim(message)//'incorrect decision-variable vector size'
      err=20; return
    endif

    ! parents in the transformed search space, scaled to [0,1]
    allocate(z(nParam,size(parents,2)),width(nParam),active(nParam))
    width=search%search_upper-search%search_lower
    active=(width > 0._rkind)
    do i=1,size(parents,2)
      do d=1,nParam
        call transform_parameter(parents(d,i),search%transformation(d),s,err,cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
        z(d,i)=0._rkind
        if(active(d)) z(d,i)=min(max((s-search%search_lower(d))/width(d),0._rkind),1._rkind)
      enddo
    enddo

    iChild=0
    iTry=0
    do while(iChild < nChild)
      iTry=iTry+1
      if(iTry > maxtry*nChild)then
        message=trim(message)//'unable to breed offspring satisfying the ordered parameter constraints'
        err=20; return
      endif
      p1=tournament(rank,distance)
      p2=tournament(rank,distance)
      call sbx(z(:,p1),z(:,p2),active,pc,eta_c,c(:,1),c(:,2))
      do iPair=1,2
        if(iChild == nChild) exit
        call mutate(c(:,iPair),active,pm,eta_m)
        do d=1,nParam
          call inverse_transform_parameter(search%search_lower(d)+c(d,iPair)*width(d), &
                                           search%transformation(d),x(d),err,cmessage)
          if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
          x(d)=min(max(x(d),search%lower(d)),search%upper(d))
        enddo
        if(.not.check_ordered_constraints(search,x)) cycle
        iChild=iChild+1
        children(:,iChild)=x
      enddo
    enddo

  end subroutine make_offspring

  ! --------------------------------------------------------------------------------------------------
  ! --- PRIVATE HELPER ROUTINES ----------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------

  ! **************************************************************************************************
  ! True when a is no worse than b in every objective and better in at least one.
  ! **************************************************************************************************
  pure logical(lgt) function dominates(a,b)
    implicit none
    real(rkind), intent(in) :: a(:),b(:)
    dominates=all(a <= b) .and. any(a < b)
  end function dominates

  ! **************************************************************************************************
  ! Binary crowded tournament: the lower front wins, then the larger crowding distance.
  ! **************************************************************************************************
  integer(i4b) function tournament(rank,distance)
    implicit none
    integer(i4b), intent(in) :: rank(:)
    real(rkind),  intent(in) :: distance(:)
    integer(i4b) :: n,i,j
    real(rkind)  :: u

    n=size(rank)
    call random_number(u)
    i=min(1+int(u*n),n)
    tournament=i
    if(n < 2) return
    call random_number(u)
    j=min(1+int(u*(n-1)),n-1)
    if(j >= i) j=j+1
    if(rank(j) < rank(i) .or. (rank(j) == rank(i) .and. distance(j) > distance(i))) tournament=j
  end function tournament

  ! **************************************************************************************************
  ! Bounded simulated binary crossover on [0,1] (Deb and Agrawal, 1995), as in Deb's NSGA-II code:
  ! each variable is recombined with probability one half, its spread bounded on both sides.
  ! **************************************************************************************************
  subroutine sbx(y1,y2,active,pc,eta,c1,c2)
    implicit none
    real(rkind),  intent(in)  :: y1(:),y2(:)     ! parents, scaled
    logical(lgt), intent(in)  :: active(:)       ! variables with room to vary
    real(rkind),  intent(in)  :: pc              ! crossover probability for the pair
    real(rkind),  intent(in)  :: eta             ! distribution index
    real(rkind),  intent(out) :: c1(:),c2(:)     ! children, scaled
    real(rkind), parameter :: eps=1.e-14_rkind
    real(rkind)  :: u,r,lo,hi,a,b
    integer(i4b) :: d

    c1=y1
    c2=y2
    call random_number(u)
    if(u > pc) return
    do d=1,size(y1)
      if(.not.active(d)) cycle
      call random_number(u)
      if(u > 0.5_rkind) cycle
      if(abs(y1(d)-y2(d)) <= eps) cycle
      lo=min(y1(d),y2(d))
      hi=max(y1(d),y2(d))
      call random_number(r)
      a=0.5_rkind*((lo+hi)-sbx_spread(r,1._rkind+2._rkind*lo/(hi-lo),eta)*(hi-lo))
      b=0.5_rkind*((lo+hi)+sbx_spread(r,1._rkind+2._rkind*(1._rkind-hi)/(hi-lo),eta)*(hi-lo))
      a=min(max(a,0._rkind),1._rkind)
      b=min(max(b,0._rkind),1._rkind)
      call random_number(u)
      if(u <= 0.5_rkind)then
        c1(d)=b; c2(d)=a
      else
        c1(d)=a; c2(d)=b
      endif
    enddo
  end subroutine sbx

  ! **************************************************************************************************
  ! SBX spread factor for a uniform draw r, given beta for the distance to the nearer bound.
  ! **************************************************************************************************
  pure real(rkind) function sbx_spread(r,beta,eta)
    implicit none
    real(rkind), intent(in) :: r,beta,eta
    real(rkind) :: alpha

    alpha=2._rkind-beta**(-(eta+1._rkind))
    if(r <= 1._rkind/alpha)then
      sbx_spread=(r*alpha)**(1._rkind/(eta+1._rkind))
    else
      sbx_spread=(1._rkind/(2._rkind-r*alpha))**(1._rkind/(eta+1._rkind))
    endif
  end function sbx_spread

  ! **************************************************************************************************
  ! Bounded polynomial mutation on [0,1] (Deb and Goyal, 1996), each variable with probability pm.
  ! **************************************************************************************************
  subroutine mutate(y,active,pm,eta)
    implicit none
    real(rkind),  intent(inout) :: y(:)          ! child, scaled
    logical(lgt), intent(in)    :: active(:)     ! variables with room to vary
    real(rkind),  intent(in)    :: pm            ! mutation probability per variable
    real(rkind),  intent(in)    :: eta           ! distribution index
    real(rkind)  :: u,r,val,deltaq,mut_pow
    integer(i4b) :: d

    mut_pow=1._rkind/(eta+1._rkind)
    do d=1,size(y)
      if(.not.active(d)) cycle
      call random_number(u)
      if(u > pm) cycle
      call random_number(r)
      if(r <= 0.5_rkind)then
        val=2._rkind*r+(1._rkind-2._rkind*r)*(1._rkind-y(d))**(eta+1._rkind)
        deltaq=val**mut_pow-1._rkind
      else
        val=2._rkind*(1._rkind-r)+2._rkind*(r-0.5_rkind)*y(d)**(eta+1._rkind)
        deltaq=1._rkind-val**mut_pow
      endif
      y(d)=min(max(y(d)+deltaq,0._rkind),1._rkind)
    enddo
  end subroutine mutate

  ! **************************************************************************************************
  ! Stable ascending sort order of a vector (insertion sort; populations are small).
  ! **************************************************************************************************
  function sort_index(v) result(idx)
    implicit none
    real(rkind), intent(in) :: v(:)
    integer(i4b)            :: idx(size(v))
    integer(i4b) :: i,j,t

    idx=[(i,i=1,size(v))]
    do i=2,size(v)
      t=idx(i)
      j=i-1
      do while(j >= 1)
        if(v(idx(j)) <= v(t)) exit
        idx(j+1)=idx(j)
        j=j-1
      enddo
      idx(j+1)=t
    enddo
  end function sort_index

end module nsga2
