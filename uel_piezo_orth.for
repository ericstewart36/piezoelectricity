!************************************************************************
! User element (UEL) for coupled large-deformation elasticity and  
!  piezoelectric behavior in three dimensions.
!************************************************************************
! Element details:
!************************************************************************
!
! Solution variables (or nodal variables) are the displacements (DOFs 1-3)
!  and the electric potential (DOF 11).
!
! Material behavior is linear elasticity with piezoelectricity.
! 
! This subroutine is for a three-dimensional 8-node isoparametric
!  brick element as shown below with 8pt (full) or 1pt (reduced) 
!  integration.
!
! In order to avoid locking for the fully-integrated element, we
!  use the F-bar method of de Souza Neto (1996).
!
!
!  8-node     8-----------7
!  brick     /|          /|       zeta
!           / |         / |       
!          5-----------6  |       |     eta
!          |  |        |  |       |   /
!          |  |        |  |       |  /
!          |  4--------|--3       | /
!          | /         | /        |/
!          |/          |/         O--------- xi
!          1-----------2        origin at cube center
!
! Eric Stewart, Spring 2021
!
!************************************************************************
! Usage:
!************************************************************************
!
! User element statement in the input file:
!  *User Element,Nodes=8,Type=U1,Iproperties=0,Properties=4,Coordinates=3,Variables=1,Unsymm
!  1,2,3,11
!
! Note: No local state variables are used in this element, thus we may set the above 
!  parameter 'Variables' to any non-zero integer.
!
!
!
!     Material Properties Vector
!     --------------------------------------------------------------
!		C1111 = props(1)
!		C1122 = props(1)
!		C2222 = props(1)
!		C1133 = props(1)
!		C2233 = props(1)
!		C3333 = props(1)
!		C1212 = props(1)
!		C1313 = props(1)
!		C2323 = props(1)
!		E31 = props(1)
!		E32 = props(1)
!		E33 = props(1)
!		E24 = props(1)
!		E15 = props(1)
!		e11 = props(1)
!		e22 = props(1)
!		e33 = props(1)
!
!************************************************************************
 
      module global

      ! This module is used to transfer SDV's from the UEL
      !  to the UVARM so that SDV's can be visualized on a
      !  dummy mesh
      !
      !  globalSdv(X,Y,Z)
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      !
      !  numElem
      !   Total number of elements in the real mesh, the dummy
      !   mesh needs to have the same number of elements, and 
      !   the dummy mesh needs to have the same number of integ
      !   points.  You must set that parameter value here.
      !
      !  ElemOffset
      !   Offset between element numbers on the real mesh and
      !    dummy mesh.  That is set in the input file, and 
      !    that value must be set here the same.

      integer numElem,ElemOffset,err

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the number of UEL elements used here
      parameter(numElem=2200)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set the offset here for UVARM plotting, must match input file!
      parameter(ElemOffset=2200)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8, allocatable :: globalSdv(:,:,:)

      integer test

      end module global

***********************************************************************

      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1     LAYER,KSPT)

      use global

      INCLUDE 'ABA_PARAM.INC'

      DIMENSION STATEV(NSTATV),COORDS(NCRDS)


      statev = 0.999
      test = max(test,noel)


      RETURN
      END

***********************************************************************

      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1 NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2 JMAC,JMATYP,MATLAYO,LACCFLA)

      ! This subroutine is used to transfer SDV's from the UEL
      !  onto the dummy mesh for viewing.  Note that an offset of
      !  ElemOffset is used between the real mesh and the dummy mesh.
      !  If your model has more than ElemOffset UEL elements, then
      !  this will need to be modified.
     
      use global
     
      include 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3 FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)

C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.
      do i=1,nuvarm
         uvar(i) = globalSdv(noel-ElemOffset,npt,i)
      enddo
c      for example
c      uvar(2) = globalSdv(noel-ElemOffset,npt,2)
c      uvar(3) = globalSdv(noel-ElemOffset,npt,3)
c      uvar(4) = globalSdv(noel-ElemOffset,npt,4)

      return
      end subroutine uvarm

****************************************************************************
	  
      subroutine UEL(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,
     +     props,nprops,coords,mcrd,nnode,uall,duall,vel,accn,jtype,
     +     time,dtime,kstep,kinc,jelem,params,ndload,jdltyp,adlmag,
     +     predef,npredf,lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,
     +     njprop,period)

      
	  use global
	  !
	  implicit none
      !
	  !
      ! variables defined in uel, passed back to Abaqus
      !
      real*8 rhs(mlvarx,*),amatrx(ndofel,ndofel),svars(*),energy(8),
     +  pnewdt
      !
      ! variables passed into UEL
      !
      integer ndofel,nrhs,nsvars,nprops,mcrd,nnode,jtype,kstep,kinc,
     +  jelem,ndload,jdltyp(mdload,*),npredf,lflags(*),mlvarx,mdload,
     +  jprops(*),njprop
      !
      real*8 props(*),coords(mcrd,nnode),uall(ndofel),duall(mlvarx,*),
     +  vel(ndofel),accn(ndofel),time(2),dtime,params(*),
     +  adlmag(mdload,*),predef(2,npredf,nnode),ddlmag(mdload,*),period
      !
      ! variables defined and used in the UEL
      !
      integer i,j,k,A11,B11,A12,B12,nInt,nIntPt,intpt,nDim,stat
      integer nlSdv, ngSdv, jj
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      parameter(nDim=3) ! number of spatial dimensions, do not change
      parameter(nInt=8) ! number of integration points
      !
      ! nInt=8: fully-integrated
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
      real*8 u(nNode,3),du(nNode,ndofel),phi(nNode),coordsC(mcrd,nNode),
     +  Ru(3*nNode,1),Rphi(nNode,1),Kuu(3*nNode,3*nNode),
     +  Kuphi(3*nNode,nNode),Kphiu(nNode,3*nNode),Kphiphi(nNode,nNode),
     +  Iden(3,3),xi(nInt,3),w(nInt),sh0(nNode),sh(nNode),dsh0(nNode,3),
     +  dshC0(nNode,3),dsh(nNode,3),dshC(nNode,3),dshxi(nNode,3),
     +  detMapJ0,detMapJ0C,detMapJ,detMapJC,Fc_tau(3,3),detFc_tau,
     +  F_tau(3,3),detF_tau,eR_tau(3,1),T_tau(3,3),D_tau(3,1),
     +  dR_tau(3,1),SR_tau(3,3),
     +  SpUUMod(3,3,3,3),SpUPhiMod(3,3,3),SpPhiUMod(3,3,3),
     +  SpPhiPhiMod(3,3),Smat(6,1),Bmat(6,3*nNode),Gmat(9,3*nNode),
     +  G0mat(9,3*nNode),Amat(9,9),Qmat(9,9),AmatPhiU(3,9),
     +  AmatUPhi(9,3),Gmatphi(3,nNode),body(3),BodyForceRes(3*nNode,1),
     +  bodyCharge,BodyChargeRes(nNode,1),Le, cbar_t
      !
      real*8 zero,one,two,half,Pi,three,third
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,Pi=3.141592653d0,
     +     three=3.d0,third=1.d0/3.d0)
      
      
      ! Check the procedure type, this should be a coupled
      !  temperature displacement, which is either 72 or 73
      !
      if((lflags(1).eq.72).or.(lflags(1).eq.73)) then
         !
         ! correct procedure specified
         !
      else
         write(*,*) 'Abaqus does not have the right procedure'
         write(*,*) 'go back and check the procedure type'
         write(*,*) 'lflags(1)=',lflags(1)
         call xit
      endif


      ! Make sure Abaqus knows you are doing a large
      !  deformation problem
      !
      if(lflags(2).eq.0) then
         !
         ! lflags(2)=0 -> small disp.
         ! lflags(2)=1 -> large disp.
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a small displacement analysis'
         write(*,*) 'go in and set nlgeom=yes'
         call xit
      endif


      ! Check to see if you are doing a general
      !  step or a linear purturbation step
      !
      if(lflags(4).eq.1) then
         !
         ! lflags(4)=0 -> general step
         ! lflags(4)=1 -> linear purturbation step
         !
         write(*,*) 'Abaqus thinks you are doing'
         write(*,*) 'a linear purturbation step'
         call xit         
      endif

      ! Get element parameters
      !
      nlSdv = jprops(1) ! number of local sdv's per integ point
      ngSdv = jprops(2) ! number of global sdv's per integ point


      ! Allocate memory for the globalSdv's used for viewing
      !  results on the dummy mesh
      !
      if(.not.allocated(globalSdv)) then
         !
         ! allocate memory for the globalSdv's
         !
         ! numElem needs to be set in the MODULE
         ! nInt needs to be set in the UEL
         !
         stat=0
c         allocate(globalSdv(numElem,nInt,ngSdv))
c         deallocate(globalSdv)
         allocate(globalSdv(numElem,nInt,ngSdv),stat=err)
         if(stat.ne.0) then
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) 'error when allocating globalSdv'
            write(*,*) '//////////////////////////////////////////////'
            write(*,*) '   stat=',stat
            write(*,*) '  ngSdv=',ngSdv
            write(*,*) '   nInt=',nInt
            write(*,*) 'numElem=',numElem
            write(*,*) '  nNode=',nNode
            write(*,*) 'lbound(globalSdv)',lbound(globalSdv)
            write(*,*) 'ubound(globalSdv)',ubound(globalSdv)
            write(*,*) '//////////////////////////////////////////////'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) 'error when allocating globalSdv'
            write(80,*) '//////////////////////////////////////////////'
            write(80,*) '   stat=',stat
            write(80,*) '  ngSdv=',ngSdv
            write(80,*) '   nInt=',nInt
            write(80,*) 'numElem=',numElem
            write(80,*) '  nNode=',nNode
            write(80,*) 'lbound(globalSdv)=',lbound(globalSdv)
            write(80,*) 'ubound(globalSdv)=',ubound(globalSdv)
            write(80,*) '//////////////////////////////////////////////'
            call xit
         endif
         write(*,*) '-------------------------------------------------'
         write(*,*) '----------- globalSDV ALLOCATED -----------------'
         write(*,*) '-------------------------------------------------'
         write(*,*) '------------- U3D8 ELEMENTS ---------------------'
         write(*,*) '-------------------------------------------------'
      endif
	  
      ! Do nothing if a ``dummy'' step
      !
      if(dtime.eq.zero) return


      ! Identity tensor
      !
      call onem(Iden)


      ! Initialize the residual and tangent matrices to zero.
      !
      Ru  = zero
      Rphi  = zero
      Kuu = zero
      Kuphi = zero
      Kphiphi = zero
      Kphiu = zero
      Energy = zero


      ! Obtain nodal displacements and potentials
      !
      k = 0
      do i=1,nNode
         do j=1,nDim
            k = k + 1
            u(i,j) = Uall(k)
            du(i,j) = DUall(k,1)
         enddo
         k = k + 1
         phi(i) = Uall(k)
      enddo


      ! Obtain current nodal coordinates
      !
      do i=1,nNode
         do j=1,nDim
            coordsC(j,i) = coords(j,i) + u(i,j)
         enddo
      enddo


      ! Impose the heuristic restriction that the displacement 
      !  increment is less than ten times the element diagonal,
      !  which is taken as an approximate element size.
      !
      Le = dsqrt(((coordsC(1,1)-coordsC(1,7))**two) + 
     +           ((coordsC(2,1)-coordsC(2,7))**two) +
     +           ((coordsC(3,1)-coordsC(3,7))**two))
      !
      do i=1,nNode
         do j=1,nDim
            if(dabs(du(i,j)).gt.10.d0*Le) then
               pnewdt = 0.5
               return
            endif
         enddo
      enddo


      !----------------------------------------------------------------
      ! 
      ! Take this opportunity to perform calculations at the element
      !  centroid.  Get the deformation gradient for use in the `F-bar' method.
      !
      ! Reference for the F-bar method:
      !  de Souza Neto, E.A., Peric, D., Dutko, M., Owen, D.R.J., 1996.
      !  Design of simple low order finite elements for large strain
      !  analysis of nearly incompressible solids. International Journal
      !  of Solids and Structures 33, 3277-3296.
      !
      ! Obtain shape functions and their local gradients at the element
      !  centroid, that means xi=eta=zeta=0.0, and nIntPt=1
      !
      if(nNode.eq.8) then
         call calcShape3DLinear(1,zero,1,sh0,dshxi)
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.8'
         call xit
      endif


      ! Map shape functions from local to global reference coordinate system
      !
      call mapShape3D(nNode,dshxi,coords,dsh0,detMapJ0,stat)
      if(stat.eq.0) then
         pnewdt = 0.5
         return
      endif


      ! Map shape functions from local to global current coordinate system
      !
      call mapShape3D(nNode,dshxi,coordsC,dshC0,detMapJ0C,stat)
      if(stat.eq.0) then
         pnewdt = 0.5
         return
      endif


      ! Calculate the deformation gradient at the element centroid
      !  at the end of the increment for use in the `F-bar' method
      !  The subscript tau denotes the time at the end of the increment.
      !
      Fc_tau = Iden
      do i=1,nDim
         do j=1,nDim
            do k=1,nNode
               Fc_tau(i,j) = Fc_tau(i,j) + dsh0(k,j)*u(k,i)
            enddo
         enddo
      enddo
      call mdet(Fc_tau,detFc_tau)
      !
      ! With the deformation gradient known at the element centroid
      !  we are now able to implement the `F-bar' method later
      !
      !----------------------------------------------------------------




      !----------------------------------------------------------------
      ! Begin the loop over integration points
      !
      ! Obtain integration point local coordinates and weights
      !
      if(nNode.eq.8) then
         if(nInt.eq.1) then
            call xint3D1pt(xi,w,nIntPt) ! 1-pt integration, nInt=1 above
         elseif(nInt.eq.8) then
            call xint3D8pt(xi,w,nIntPt) ! 8-pt integration, nInt=8 above
         else
            write(*,*) 'Invalid number of int points, nInt=',nInt
            call xit
         endif
      else
         write(*,*) 'Incorrect number of nodes: nNode.ne.8'
         call xit
      endif


      ! Loop over integration points
      !
	  jj=0
      do intpt=1,nIntPt

         ! Obtain state variables from previous increment
         !
         if((kinc.le.1).and.(kstep.eq.1)) then
            !
            ! this is the first increment, of the first step
            !  give initial conditions
            !
            cbar_t  = one
            !
         else
            !
            ! this is not the first increment, read old values
            !
            cbar_t  = svars(1+jj)
            !
         endif
		 
         ! Obtain shape functions and their local gradients
         !
         if(nNode.eq.8) then
            call calcShape3DLinear(nIntPt,xi,intpt,sh,dshxi)
         else
            write(*,*) 'Incorrect number of nodes: nNode.ne.8'
            call xit
         endif


         ! Map shape functions from local to global reference coordinate system
         !
         call mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif


         ! Map shape functions from local to global current coordinate system
         !
         call mapShape3D(nNode,dshxi,coordsC,dshC,detMapJC,stat)
         if(stat.eq.0) then
            pnewdt = 0.5
            return
         endif


         ! Obtain the referential electric field at this integration point
         !
         eR_tau = zero
         do k=1,nNode
            do i=1,nDim
               eR_tau(i,1) = eR_tau(i,1) + phi(k)*dsh(k,i)
            enddo
         enddo
         eR_tau = -eR_tau

         ! Obtain the deformation gradient at this integration point.
         !  The subscript tau denotes the time at the end of the increment.
         !
         F_tau = Iden
         do i=1,nDim
            do j=1,nDim
               do k=1,nNode
                  F_tau(i,j) = F_tau(i,j) + dsh(k,j)*u(k,i)
               enddo
            enddo
         enddo
         
         
         ! Modify the deformation gradient for the `F-bar' method
         !  only when using the 8 node fully integrated linear
         !  element, do not use the `F-bar' method for reduced element
         !
         if((nNode.eq.8).and.(nInt.eq.8)) then
            call mdet(F_tau,detF_tau)
            F_tau = ((detFc_tau/detF_tau)**third)*F_tau
         endif


         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         !
         ! Perform the constitutive update at this integ. point
         !
         !
		 call integ(
!             material prop. inputs:
     +        props,nprops,
!             constitutive inputs:
     +        F_tau,eR_tau,
!             constitutive outputs:
     +        T_tau, SR_tau, d_tau, dR_tau,
!             tangent outputs:
     +        SpUUMod,SpUPhiMod,SpPhiUMod,SpPhiPhiMod
     +        )
	 
	     ! Save the state variables at this integ point in the
         !  global array used for plotting field output
         !
         globalSdv(jelem,intPt,1) = SR_tau(1,1) !-\
         globalSdv(jelem,intPt,2) = SR_tau(2,2) !  |
         globalSdv(jelem,intPt,3) = SR_tau(3,3) !  |-->   Stress
         globalSdv(jelem,intPt,4) = SR_tau(2,3) !  |-->  Components
         globalSdv(jelem,intPt,5) = SR_tau(1,3) !  |
         globalSdv(jelem,intPt,6) = SR_tau(1,2) !-/
		 globalSdv(jelem,intPt,7) = dR_tau(1,1) !-\
         globalSdv(jelem,intPt,8) = dR_tau(2,1) !  |--> Electric disp
         globalSdv(jelem,intPt,9) = dR_tau(3,1) !-/
         !
         !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

         ! Save the state variables at this integ point
         !  at the end of the increment
         !
         svars(1+jj) = one
         jj = jj + nlSdv ! setup for the next intPt

         ! Compute/update the displacement residual vector
         !
         Smat(1,1) = T_tau(1,1)
         Smat(2,1) = T_tau(2,2)
         Smat(3,1) = T_tau(3,3)
         Smat(4,1) = T_tau(1,2)
         Smat(5,1) = T_tau(2,3)
         Smat(6,1) = T_tau(1,3)
         !
         Bmat = zero
         do k=1,nNode
            Bmat(1,1+nDim*(k-1)) = dshC(k,1)
            Bmat(2,2+nDim*(k-1)) = dshC(k,2)
            Bmat(3,3+nDim*(k-1)) = dshC(k,3)
            Bmat(4,1+nDim*(k-1)) = dshC(k,2)
            Bmat(4,2+nDim*(k-1)) = dshC(k,1)
            Bmat(5,2+nDim*(k-1)) = dshC(k,3)
            Bmat(5,3+nDim*(k-1)) = dshC(k,2)
            Bmat(6,1+nDim*(k-1)) = dshC(k,3)
            Bmat(6,3+nDim*(k-1)) = dshC(k,1)
         enddo
         !
         body = zero ! The body force vector may be specified here
         !
         BodyForceRes = zero
         do k=1,nNode
            BodyForceRes(1+nDim*(k-1),1) = sh(k)*body(1)
            BodyForceRes(2+nDim*(k-1),1) = sh(k)*body(2)
            BodyForceRes(3+nDim*(k-1),1) = sh(k)*body(3)
         enddo
         !
         Ru = Ru + matmul(transpose(Bmat),Smat)*detmapJC*w(intpt)
     +        - BodyForceRes*detmapJC*w(intpt)


         ! Compute/update the electric potential residual vector
         !
         Gmatphi = zero
         do i = 1,nDim
            do k = 1,nNode
               Gmatphi(i,k) = dshC(k,i)
            end do
         end do
         !
         bodyCharge = zero ! the free charge density may be specified here
         !
         BodyChargeRes = zero
         do k=1,nNode
            BodyChargeRes(k,1) = sh(k)*bodyCharge
         end do
         !
         Rphi = Rphi + matmul(transpose(Gmatphi),D_tau)
     +                 *detmapJC*w(intpt) 
     +               + BodyChargeRes*detmapJC*w(intpt)


         ! Compute/update the displacement tangent matrix
         !
         Gmat = zero
         do k=1,nNode
            Gmat(1,1+nDim*(k-1)) = dshC(k,1)
            Gmat(2,2+nDim*(k-1)) = dshC(k,1)
            Gmat(3,3+nDim*(k-1)) = dshC(k,1)
            Gmat(4,1+nDim*(k-1)) = dshC(k,2)
            Gmat(5,2+nDim*(k-1)) = dshC(k,2)
            Gmat(6,3+nDim*(k-1)) = dshC(k,2)
            Gmat(7,1+nDim*(k-1)) = dshC(k,3)
            Gmat(8,2+nDim*(k-1)) = dshC(k,3)
            Gmat(9,3+nDim*(k-1)) = dshC(k,3)
         enddo
         !
         G0mat = zero
         do k=1,nNode
            G0mat(1,1+nDim*(k-1)) = dshC0(k,1)
            G0mat(2,2+nDim*(k-1)) = dshC0(k,1)
            G0mat(3,3+nDim*(k-1)) = dshC0(k,1)
            G0mat(4,1+nDim*(k-1)) = dshC0(k,2)
            G0mat(5,2+nDim*(k-1)) = dshC0(k,2)
            G0mat(6,3+nDim*(k-1)) = dshC0(k,2)
            G0mat(7,1+nDim*(k-1)) = dshC0(k,3)
            G0mat(8,2+nDim*(k-1)) = dshC0(k,3)
            G0mat(9,3+nDim*(k-1)) = dshC0(k,3)
         enddo
         !
         Amat = zero
         Amat(1,1) = SpUUMod(1,1,1,1)
         Amat(1,2) = SpUUMod(1,1,2,1)
         Amat(1,3) = SpUUMod(1,1,3,1)
         Amat(1,4) = SpUUMod(1,1,1,2)
         Amat(1,5) = SpUUMod(1,1,2,2)
         Amat(1,6) = SpUUMod(1,1,3,2)
         Amat(1,7) = SpUUMod(1,1,1,3)
         Amat(1,8) = SpUUMod(1,1,2,3)
         Amat(1,9) = SpUUMod(1,1,3,3)
         Amat(2,1) = SpUUMod(2,1,1,1)
         Amat(2,2) = SpUUMod(2,1,2,1)
         Amat(2,3) = SpUUMod(2,1,3,1)
         Amat(2,4) = SpUUMod(2,1,1,2)
         Amat(2,5) = SpUUMod(2,1,2,2)
         Amat(2,6) = SpUUMod(2,1,3,2)
         Amat(2,7) = SpUUMod(2,1,1,3)
         Amat(2,8) = SpUUMod(2,1,2,3)
         Amat(2,9) = SpUUMod(2,1,3,3)
         Amat(3,1) = SpUUMod(3,1,1,1)
         Amat(3,2) = SpUUMod(3,1,2,1)
         Amat(3,3) = SpUUMod(3,1,3,1)
         Amat(3,4) = SpUUMod(3,1,1,2)
         Amat(3,5) = SpUUMod(3,1,2,2)
         Amat(3,6) = SpUUMod(3,1,3,2)
         Amat(3,7) = SpUUMod(3,1,1,3)
         Amat(3,8) = SpUUMod(3,1,2,3)
         Amat(3,9) = SpUUMod(3,1,3,3)
         Amat(4,1) = SpUUMod(1,2,1,1)
         Amat(4,2) = SpUUMod(1,2,2,1)
         Amat(4,3) = SpUUMod(1,2,3,1)
         Amat(4,4) = SpUUMod(1,2,1,2)
         Amat(4,5) = SpUUMod(1,2,2,2)
         Amat(4,6) = SpUUMod(1,2,3,2)
         Amat(4,7) = SpUUMod(1,2,1,3)
         Amat(4,8) = SpUUMod(1,2,2,3)
         Amat(4,9) = SpUUMod(1,2,3,3)
         Amat(5,1) = SpUUMod(2,2,1,1)
         Amat(5,2) = SpUUMod(2,2,2,1)
         Amat(5,3) = SpUUMod(2,2,3,1)
         Amat(5,4) = SpUUMod(2,2,1,2)
         Amat(5,5) = SpUUMod(2,2,2,2)
         Amat(5,6) = SpUUMod(2,2,3,2)
         Amat(5,7) = SpUUMod(2,2,1,3)
         Amat(5,8) = SpUUMod(2,2,2,3)
         Amat(5,9) = SpUUMod(2,2,3,3)
         Amat(6,1) = SpUUMod(3,2,1,1)
         Amat(6,2) = SpUUMod(3,2,2,1)
         Amat(6,3) = SpUUMod(3,2,3,1)
         Amat(6,4) = SpUUMod(3,2,1,2)
         Amat(6,5) = SpUUMod(3,2,2,2)
         Amat(6,6) = SpUUMod(3,2,3,2)
         Amat(6,7) = SpUUMod(3,2,1,3)
         Amat(6,8) = SpUUMod(3,2,2,3)
         Amat(6,9) = SpUUMod(3,2,3,3)
         Amat(7,1) = SpUUMod(1,3,1,1)
         Amat(7,2) = SpUUMod(1,3,2,1)
         Amat(7,3) = SpUUMod(1,3,3,1)
         Amat(7,4) = SpUUMod(1,3,1,2)
         Amat(7,5) = SpUUMod(1,3,2,2)
         Amat(7,6) = SpUUMod(1,3,3,2)
         Amat(7,7) = SpUUMod(1,3,1,3)
         Amat(7,8) = SpUUMod(1,3,2,3)
         Amat(7,9) = SpUUMod(1,3,3,3)
         Amat(8,1) = SpUUMod(2,3,1,1)
         Amat(8,2) = SpUUMod(2,3,2,1)
         Amat(8,3) = SpUUMod(2,3,3,1)
         Amat(8,4) = SpUUMod(2,3,1,2)
         Amat(8,5) = SpUUMod(2,3,2,2)
         Amat(8,6) = SpUUMod(2,3,3,2)
         Amat(8,7) = SpUUMod(2,3,1,3)
         Amat(8,8) = SpUUMod(2,3,2,3)
         Amat(8,9) = SpUUMod(2,3,3,3)
         Amat(9,1) = SpUUMod(3,3,1,1)
         Amat(9,2) = SpUUMod(3,3,2,1)
         Amat(9,3) = SpUUMod(3,3,3,1)
         Amat(9,4) = SpUUMod(3,3,1,2)
         Amat(9,5) = SpUUMod(3,3,2,2)
         Amat(9,6) = SpUUMod(3,3,3,2)
         Amat(9,7) = SpUUMod(3,3,1,3)
         Amat(9,8) = SpUUMod(3,3,2,3)
         Amat(9,9) = SpUUMod(3,3,3,3)
         !
         Qmat = zero
         Qmat(1,1) = third*(Amat(1,1)+Amat(1,5)+Amat(1,9)) 
     +        - (two/three)*T_tau(1,1)
         Qmat(2,1) = third*(Amat(2,1)+Amat(2,5)+Amat(2,9))
     +        - (two/three)*T_tau(2,1)
         Qmat(3,1) = third*(Amat(3,1)+Amat(3,5)+Amat(3,9))
     +        - (two/three)*T_tau(3,1)
         Qmat(4,1) = third*(Amat(4,1)+Amat(4,5)+Amat(4,9))
     +        - (two/three)*T_tau(1,2)
         Qmat(5,1) = third*(Amat(5,1)+Amat(5,5)+Amat(5,9))
     +        - (two/three)*T_tau(2,2)
         Qmat(6,1) = third*(Amat(6,1)+Amat(6,5)+Amat(6,9))
     +        - (two/three)*T_tau(3,2)
         Qmat(7,1) = third*(Amat(7,1)+Amat(7,5)+Amat(7,9))
     +        - (two/three)*T_tau(1,3)
         Qmat(8,1) = third*(Amat(8,1)+Amat(8,5)+Amat(8,9))
     +        - (two/three)*T_tau(2,3)
         Qmat(9,1) = third*(Amat(9,1)+Amat(9,5)+Amat(9,9))
     +        - (two/three)*T_tau(3,3)
         Qmat(1,5) = Qmat(1,1)
         Qmat(2,5) = Qmat(2,1)
         Qmat(3,5) = Qmat(3,1)
         Qmat(4,5) = Qmat(4,1)
         Qmat(5,5) = Qmat(5,1)
         Qmat(6,5) = Qmat(6,1)
         Qmat(7,5) = Qmat(7,1)
         Qmat(8,5) = Qmat(8,1)
         Qmat(9,5) = Qmat(9,1)
         Qmat(1,9) = Qmat(1,1)
         Qmat(2,9) = Qmat(2,1)
         Qmat(3,9) = Qmat(3,1)
         Qmat(4,9) = Qmat(4,1)
         Qmat(5,9) = Qmat(5,1)
         Qmat(6,9) = Qmat(6,1)
         Qmat(7,9) = Qmat(7,1)
         Qmat(8,9) = Qmat(8,1)
         Qmat(9,9) = Qmat(9,1)
         !
         if((nNode.eq.8).and.(nInt.eq.8)) then
            !
            ! This is the tangent using the F-bar method with the
            !  8 node fully-integrated element
            !
            Kuu = Kuu
     +           - matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           *detMapJC*w(intpt)
     +           - matmul(transpose(Gmat),matmul(Qmat,
     +           (G0mat-Gmat)))*detMapJC*w(intpt)
         else
            !
            ! This is the tangent NOT using the F-bar method with all
            !  other elements
            !
            Kuu = Kuu -
     +           matmul(matmul(transpose(Gmat),Amat),Gmat)
     +           *detMapJC*w(intpt)
         endif


         ! Compute/update the electric potential tangent matrix
         !
         Kphiphi = Kphiphi - 
     +             matmul(matmul(transpose(Gmatphi),SpPhiPhiMod),
     +                    Gmatphi)*detMapJC*w(intpt)


         ! Compute/update the electric potential/displacement tangent matrix
         !   Strictly speaking, use of 'F-bar' will affect this tangent;
         !   however, the effect is expected to be small and is neglected here.
         !
         AmatPhiU = zero
         AmatPhiU(1,1) = SpPhiUMod(1,1,1)
         AmatPhiU(1,2) = SpPhiUMod(1,2,1)
         AmatPhiU(1,3) = SpPhiUMod(1,3,1)
         AmatPhiU(1,4) = SpPhiUMod(1,1,2)
         AmatPhiU(1,5) = SpPhiUMod(1,2,2)
         AmatPhiU(1,6) = SpPhiUMod(1,3,2)
         AmatPhiU(1,7) = SpPhiUMod(1,1,3)
         AmatPhiU(1,8) = SpPhiUMod(1,2,3)
         AmatPhiU(1,9) = SpPhiUMod(1,3,3)
         AmatPhiU(2,1) = SpPhiUMod(2,1,1)
         AmatPhiU(2,2) = SpPhiUMod(2,2,1)
         AmatPhiU(2,3) = SpPhiUMod(2,3,1)
         AmatPhiU(2,4) = SpPhiUMod(2,1,2)
         AmatPhiU(2,5) = SpPhiUMod(2,2,2)
         AmatPhiU(2,6) = SpPhiUMod(2,3,2)
         AmatPhiU(2,7) = SpPhiUMod(2,1,3)
         AmatPhiU(2,8) = SpPhiUMod(2,2,3)
         AmatPhiU(2,9) = SpPhiUMod(2,3,3)
         AmatPhiU(3,1) = SpPhiUMod(3,1,1)
         AmatPhiU(3,2) = SpPhiUMod(3,2,1)
         AmatPhiU(3,3) = SpPhiUMod(3,3,1)
         AmatPhiU(3,4) = SpPhiUMod(3,1,2)
         AmatPhiU(3,5) = SpPhiUMod(3,2,2)
         AmatPhiU(3,6) = SpPhiUMod(3,3,2)
         AmatPhiU(3,7) = SpPhiUMod(3,1,3)
         AmatPhiU(3,8) = SpPhiUMod(3,2,3)
         AmatPhiU(3,9) = SpPhiUMod(3,3,3)
         !
         Kphiu = Kphiu - matmul(matmul(transpose(Gmatphi),AmatPhiU),
     +                          Gmat)*detMapJC*w(intpt)


         ! Compute/update the displacement/electric potential tangent matrix
         !
         AmatUPhi = zero
         AmatUPhi(1,1) = SpUPhiMod(1,1,1)
         AmatUPhi(2,1) = SpUPhiMod(2,1,1)
         AmatUPhi(3,1) = SpUPhiMod(3,1,1)
         AmatUPhi(4,1) = SpUPhiMod(1,2,1)
         AmatUPhi(5,1) = SpUPhiMod(2,2,1)
         AmatUPhi(6,1) = SpUPhiMod(3,2,1)
         AmatUPhi(7,1) = SpUPhiMod(1,3,1)
         AmatUPhi(8,1) = SpUPhiMod(2,3,1)
         AmatUPhi(9,1) = SpUPhiMod(3,3,1)
         AmatUPhi(1,2) = SpUPhiMod(1,1,2)
         AmatUPhi(2,2) = SpUPhiMod(2,1,2)
         AmatUPhi(3,2) = SpUPhiMod(3,1,2)
         AmatUPhi(4,2) = SpUPhiMod(1,2,2)
         AmatUPhi(5,2) = SpUPhiMod(2,2,2)
         AmatUPhi(6,2) = SpUPhiMod(3,2,2)
         AmatUPhi(7,2) = SpUPhiMod(1,3,2)
         AmatUPhi(8,2) = SpUPhiMod(2,3,2)
         AmatUPhi(9,2) = SpUPhiMod(3,3,2)
         AmatUPhi(1,3) = SpUPhiMod(1,1,3)
         AmatUPhi(2,3) = SpUPhiMod(2,1,3)
         AmatUPhi(3,3) = SpUPhiMod(3,1,3)
         AmatUPhi(4,3) = SpUPhiMod(1,2,3)
         AmatUPhi(5,3) = SpUPhiMod(2,2,3)
         AmatUPhi(6,3) = SpUPhiMod(3,2,3)
         AmatUPhi(7,3) = SpUPhiMod(1,3,3)
         AmatUPhi(8,3) = SpUPhiMod(2,3,3)
         AmatUPhi(9,3) = SpUPhiMod(3,3,3)
         !
         Kuphi = Kuphi - matmul(matmul(transpose(Gmat),AmatUPhi),
     +                          Gmatphi)*detMapJC*w(intpt)

      enddo
      !
      ! End the loop over integration points
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Traction terms and surface charge terms have not been
      !  implemented in the three-dimensional code at this time.
      ! Mechanical, traction- and pressure-type boundary conditions 
      !  may be applied to the dummy mesh using the Abaqus built-in 
      !  commands *Dload or *Dsload.
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      ! Return Abaqus the RHS vector and the Stiffness matrix.  This
      !  is essentially giving Abaqus the residual and the tangent matrix.
      !
      ! Return Abaqus the right hand side vector
      !
      do i=1,nNode
         A11 = (nDim+1)*(i-1)+1
         A12 = nDim*(i-1)+1
         !
         ! displacement
         !
         rhs(A11,1) = Ru(A12,1)
         rhs(A11+1,1) = Ru(A12+1,1)
         rhs(A11+2,1) = Ru(A12+2,1)
         !
         ! electric potential
         !
         rhs(A11+3,1) = Rphi(i,1)
      enddo
      !
      ! Return Abaqus the tangent matrix
      !
      amatrx = zero
      do i=1,nNode
         do j=1,nNode
            A11 = (nDim+1)*(i-1)+1
            A12 = nDim*(i-1)+1
            B11 = (nDim+1)*(j-1)+1
            B12 = nDim*(j-1)+1
            !
            ! displacement
            !
            amatrx(A11,B11)     = Kuu(A12,B12)
            amatrx(A11,B11+1)   = Kuu(A12,B12+1)
            amatrx(A11,B11+2)   = Kuu(A12,B12+2)
            amatrx(A11+1,B11)   = Kuu(A12+1,B12)
            amatrx(A11+1,B11+1) = Kuu(A12+1,B12+1)
            amatrx(A11+1,B11+2) = Kuu(A12+1,B12+2)
            amatrx(A11+2,B11)   = Kuu(A12+2,B12)
            amatrx(A11+2,B11+1) = Kuu(A12+2,B12+1)
            amatrx(A11+2,B11+2) = Kuu(A12+2,B12+2)
            !
            ! electric potential
            !
            amatrx(A11+3,B11+3) = Kphiphi(i,j)
            !
            ! displacement/electric potential
            !
            amatrx(A11,B11+3) = Kuphi(A12,j)
            amatrx(A11+1,B11+3) = Kuphi(A12+1,j)
            amatrx(A11+2,B11+3) = Kuphi(A12+2,j)
            !
            ! electric potential/displacement
            !
            amatrx(A11+3,B11) = Kphiu(i,B12)
            amatrx(A11+3,B11+1) = Kphiu(i,B12+1)
            amatrx(A11+3,B11+2) = Kphiu(i,B12+2)
            !
         enddo
      enddo
      !
      ! End return of RHS and AMATRX
      !----------------------------------------------------------------


      return
      end subroutine uel

!************************************************************************
!     Material subroutine
!************************************************************************

      subroutine integ(
!             material prop. inputs:
     +        props,nprops,
!             constitutive inputs:
     +        F_tau,eR_tau,
!             constitutive outputs:
     +        T_tau, SR_tau, d_tau, dR_tau,
!             tangent outputs:
     +        SpUUMod,SpUPhiMod,SpPhiUMod,SpPhiPhiMod
     +        )


      ! This subroutine computes everything required for the time integration
      ! of the problem.
      !
      ! Inputs:
      !  1) material parameters, props(nprops)
      !  2) deformation gradient, F_tau(3,3)
      !  3) referential electric field,  eR_tau(3,1)
      !
      ! Outputs:
      !  1) Cauchy stress, T_tau(3,3)
      !  2) 2nd Piola stress, SR_tau(3,3)
      !  3) electric displacement, d_tau(3,1)
	  !  4) spatial tangent moduli

      implicit none

      integer i,j,k,l,m,n,nprops,nargs,stat, nDofN
      parameter(nargs=13)

      ! Inputs
      real*8 props(nprops),dtime
      real*8 F_tau(3,3),phi_tau,e_tau(3,1)
      real*8 F_t(3,3),phi_t
      ! Outputs
      real*8 T_tau(3,3),d_tau(3,1)
	  real*8 SpUUMod(3,3,3,3),SpUPhiMod(3,3,3),SpPhiUMod(3,3,3)
      real*8 SpPhiPhiMod(3,3)


      ! Properties
      real*8 Eyoung,anu,Gshear,Kbulk,Lambda, poisson
      real*8 E31,E32,E33,E24,E15,eps11,eps22,eps33
      real*8 CmatVoigt(6,6), EmatVoigt(3,6), EpsMat(3,3)
      real*8 Cmat(3,3,3,3), Emat(3,3,3)
	  real*8 C1111, C1122, C2222, C1133, C2233, C3333,
     +                    C1212, C1313, C2323


      ! In-subroutine use
      real*8 Iden(3,3)
      real*8 Finv(3,3),detF
      real*8 eR_tau(3,1), dR_tau(3,1)
      real*8 Etot_tau(3,3), Etot_tau_vec(6,1), EeR_tau_vec(6,1)
      real*8 SR_tau(3,3), SR_tau_vec(6,1)
      real*8 args(nargs)
	  real*8 dTRdF(3,3,3,3),dTRde(3,3,3),dddF(3,3,3)
      real*8 ddde(3,3)

      real*8 zero,one,two,three,third,half,Rgas,Farad
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0,
     +     half=1.d0/2.d0)


      ! Identity tensor
      !
      call onem(Iden)


      ! Obtain material properties
		C1111 = props(1)
		C1122 = props(2)
		C2222 = props(3)
		C1133 = props(4)
		C2233 = props(5)
		C3333 = props(6)
		C1212 = props(7)
		C1313 = props(8)
		C2323 = props(9)
		E31 = props(10)
		E32 = props(11)
		E33 = props(12)
		E24 = props(13)
		E15 = props(14)
		eps11 = props(15)
		eps22 = props(16)
		eps33 = props(17)
      
      
      ! Construct the elastic stiffness "C" matrix
	  ! (6x6 Voigt representation of a fourth-order tensor)
	  !
	  call buildCmatVoigt(C1111, C1122, C2222, C1133, C2233, C3333,
     +                    C1212, C1313, C2323, CmatVoigt)
	  
	  ! Construct the piezoelectric "E" matrix
	  ! (3x6 Voigt representation of a third-order tensor)
      !	  
	  call buildEmatVoigt(E31, E32, E33, E24, E15, EmatVoigt)
	  
	  ! Construct the permittivity "Eps" matrix
	  ! (3x3 Voigt representation of a second-order tensor)
      !	  
	  call buildEpsMat(eps11, eps22, eps33, EpsMat)


      ! Green Strain 
	  !
	  Etot_tau = half*(matmul(transpose(F_tau), F_tau) - Iden)
	  ! 
	  ! push to Voigt notation (6x1)
	  !
	  call pushsvStrain(Etot_tau, Etot_tau_vec, 2, 3, 6) 

C 	  write(*,*) 'E_tensor = ', Etot_tau
C 	  write(*,*) 'E_vector = ', Etot_tau_vec

	  
	  ! Second Piola stress, S^R
	  !
	  SR_tau_vec = matmul(CmatVoigt, Etot_tau_vec) - matmul(transpose(EmatVoigt), eR_tau)
	  ! 
	  ! push from Voigt back to matrix notation
	  !
	  call pushsv(SR_tau, SR_tau_vec, 1, 3, 6)
	  
!  	  write(*,*) 'EmatVoigt = ', EmatVoigt

      ! referential electric displacement
	  !
	  dR_tau = matmul(EmatVoigt, Etot_tau_vec) + matmul(EpsMat, eR_tau)

	  
      ! Convert to Cauchy stress and spatial electric displacement
      !
      ! Compute the determinant of F, detF = J.
      !
      call mdet(F_tau,detF)
	  !
      T_tau = matmul(F_tau,matmul(SR_tau,transpose(F_tau)))/detF
      !
      d_tau = matmul(F_tau, dR_tau)/detF

!*****************************************************************
!   Spatial tangent calculations
!*****************************************************************
      
	  ! Build full third- and fourth- order representations of 
	  ! elasticity and piezoelectricity tensors for ease of 
	  ! indicial calculations.
	  !
	  ! (The permittivity tensor is the same in Voigt and 
	  !  standard notation)
	  
	  ! Construct the elastic stiffness "C" tensor
	  !
	  call buildCmat(C1111, C1122, C2222, C1133, C2233, C3333,
     +                    C1212, C1313, C2323, Cmat)
	  
	  ! Construct the piezoelectric "E" tensor
      !	  
	  call buildEmat(E31, E32, E33, E24, E15, Emat)

!*****************************************************************
	  
      ! First, we calculate the material tangent moduli
      !
	  ! d(stress)/d(deformation gradient)
	  !
      dTRdF = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
			      do m=1,3
				     do n=1,3
                  dTRdF(i,j,k,l) = dTRdF(i,j,k,l) + 
     +                (
     +                 iden(i,k)*SR_tau(l,j) + 
     +                  F_tau(i,m)*F_tau(k,n)*Cmat(m,j,l,n)
     +                 )
                     enddo
			      enddo
			   enddo
            enddo
         enddo
      enddo
	  !
	  ! d(Stress)/d(electric field)
	  !
	  dTRde = zero
      do i=1,3
         do j=1,3
            do k=1,3
			   do m=1,3
               dTRde(i,j,k) = dTRde(i,j,k) -
     +              (-F_tau(i,m)*Emat(k,j,m))
	           end do 
            end do
         end do
      end do
	  !
	  ! d(electric displacement)/d(deformation gradient)
	  !
	  dddF = zero
      do i=1,3
         do j=1,3
            do k=1,3
			   do m=1,3
               dddF(i,j,k) = dddF(i,j,k) +
     +              F_tau(j,m)*Emat(i,k,m)
	           end do 
            end do
         end do
      end do
	  !
	  ! d(electric displacement)/d(electric field)
	  !
	  ddde = zero
      do i=1,3
         do j=1,3
               ddde(i,j) = ddde(i,j) -
     +              EpsMat(i,j)
         end do
      end do

!*****************************************************************

      ! Next, convert to spatial tangent moduli
      !
	  ! Spatial stress/deformation gradient 
	  !
      SpUUMod = zero
      do i=1,3
         do j=1,3
            do k=1,3
               do l=1,3
                  do m=1,3
                     do n=1,3
                        SpUUMod(i,j,k,l) = SpUUMod(i,j,k,l) +
     +                      (dTRdF(i,m,k,n)*F_tau(j,m)*F_tau(l,n))/detF
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
	  !
      ! Spatial stress/electric potential modulus
      !
      SpUPhiMod = zero
      do i=1,3
         do j=1,3
            do l=1,3
               do m=1,3
                  do n=1,3
                        SpUPhiMod(i,j,l) = SpUPhiMod(i,j,l) +
     +                      (dTRde(i,m,n)*F_tau(j,m)*F_tau(l,n))/detF
                  enddo
               enddo
            enddo
         enddo
      enddo
	  !
      ! Spatial electric displacement/strain modulus
      !
      SpPhiUMod = zero
      do j=1,3
         do l=1,3
            do k=1,3
               do m=1,3
                  do n=1,3
                        SpPhiUMod(j,l,k) = SpPhiUMod(j,l,k) +
     +                      (dddF(m,k,n)*F_tau(j,m)*F_tau(l,n))/detF
                  enddo
               enddo
            enddo
         enddo
      enddo
      ! 
      ! Spatial electric displacement/electric potential modulus
      !
      SpPhiPhiMod = zero
	  do j=1,3
         do l=1,3
            do m=1,3
               do n=1,3
                        SpPhiPhiMod(j,l) = SpPhiPhiMod(j,l) +
     +                      (ddde(m,n)*F_tau(j,m)*F_tau(l,n))/detF
               enddo
            enddo
         enddo
      enddo
	  
      return
      end subroutine integ
	  
	  
****************************************************************************

	  subroutine buildCmatVoigt(C1111, C1122, C2222, C1133, C2233, C3333,
     +                    C1212, C1313, C2323, Cmat)

      ! This subroutine builds the 6x6 isotropic stiffness matrix C
	  ! (in Voigt notation) from Young's modulus and Poisson ratio
	  
	  
      implicit none
	  
	  real*8 Cmat(6,6), C1111, C1122, C2222, C1133, C2233, C3333,
     +                    C1212, C1313, C2323
	  
	  integer i
	  
	  real*8 zero,one,two,three,third,half
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0,
     +     half=1.d0/2.d0)

	  Cmat = zero
	  
      ! This one is just entry-by-entry construction
	  Cmat(1,1) = C1111
	  !
	  Cmat(1,2) = C1122
	  Cmat(2,1) = C1122
	  !
	  Cmat(2,2) = C2222
	  !
	  Cmat(1,3) = C1133
	  Cmat(3,1) = C1133
	  !
	  Cmat(2,3) = C2233
	  Cmat(3,2) = C2233
	  !
	  Cmat(3,3) = C3333
	  !
	  Cmat(4,4) = C2323
	  Cmat(5,5) = C1313
	  Cmat(6,6) = C1212
	  
      return
      end subroutine buildCmatVoigt
	  
****************************************************************************

	  subroutine buildEmatVoigt(E31, E32, E33, E24, E15, Emat)

      ! This subroutine builds the 3x6 piezoelectric matrix E
	  ! (in Voigt notation) from five piezoelectric coefficients
	  
	  
      implicit none
	  
	  real*8 E31, E32, E33, E24, E15, Emat(3,6)

      ! Initialize Dmat
	  !
	  Emat = 0.d0
	  !
	  ! populate non-zero entries
	  !
	  Emat(3,1) = E31
	  Emat(3,2) = E32
	  Emat(3,3) = E33
	  Emat(2,4) = E24
	  Emat(1,5) = E15
	  
!	  write(*,*) 'Emat =', Emat
	  
      return
      end subroutine buildEmatVoigt
	  
****************************************************************************

	  subroutine buildEpsMat(eps11, eps22, eps33, EpsMat)

      ! This subroutine builds the 3x3 permittivity matrix EpsMat
	  ! from three permittivity values.
	  
	  
      implicit none
	  
	  real*8 eps11, eps22, eps33, EpsMat(3,3)

      ! Initialize EpsMat
	  !
	  EpsMat = 0.d0
	  !
	  ! populate non-zero entries
	  !
	  EpsMat(1,1) = eps11
	  EpsMat(2,2) = eps22
	  EpsMat(3,3) = eps33
	  
      return
      end subroutine buildEpsMat

****************************************************************************
****************************************************************************

	  subroutine buildCmat(C1111, C1122, C2222, C1133, C2233, C3333,
     +                    C1212, C1313, C2323, Cmat)

      ! This subroutine builds the 3x3x3x3 isotropic stiffness matrix C
	  ! from Young's modulus and Poisson ratio.
	  
	  
      implicit none
	  
	  real*8 Cmat(3,3,3,3), C1111, C1122, C2222, C1133, C2233, C3333,
     +                    C1212, C1313, C2323
	  
	  integer i, j, k, l
	  
	  real*8 zero,one,two,three,third,half
      parameter(zero=0.d0,one=1.d0,two=2.d0,three=3.d0,third=1.d0/3.d0,
     +     half=1.d0/2.d0)
	  
	  ! Construct elements of C
	  !
	  cmat = zero
	  !
	  cmat(1,1,1,1) = C1111
	  cmat(2,2,2,2) = C2222
	  cmat(3,3,3,3) = C3333
	  !
	  cmat(1,1,2,2) = C1122
	  cmat(2,2,1,1) = C1122
	  !
	  cmat(1,1,3,3) = C1133
	  cmat(3,3,1,1) = C1133
	  !
	  cmat(2,2,3,3) = C2233
	  cmat(3,3,2,2) = C2233
	  !
	  cmat(1,2,1,2) = C1212
	  cmat(1,2,2,1) = C1212
	  cmat(2,1,1,2) = C1212
	  cmat(2,1,2,1) = C1212
	  !
	  cmat(1,3,1,3) = C1313
	  cmat(1,3,3,1) = C1313
	  cmat(3,1,1,3) = C1313
	  cmat(3,1,3,1) = C1313
      !
	  cmat(2,3,2,3) = C2323
	  cmat(2,3,3,2) = C2323
	  cmat(3,2,2,3) = C2323
	  cmat(3,2,3,2) = C2323

	  
      return
      end subroutine buildCmat
	  
****************************************************************************

	  subroutine buildEmat(E31, E32, E33, E24, E15, Emat)

      ! This subroutine builds the 3x3x3 piezoelectric tensor E
	  ! from five piezoelectric coefficients.
	  
	  
      implicit none
	  
	  real*8 E31, E32, E33, E24, E15, Emat(3,3,3)

      ! Initialize Dmat
	  !
	  Emat = 0.d0
	  !
	  ! populate non-zero entries
	  !
      Emat(3,1,1) = E31;
      Emat(3,2,2) = E32;
      Emat(3,3,3) = E33;
      !
      Emat(2,2,3) = E24;
      Emat(2,3,2) = E24;
      !
      Emat(1,1,3) = E15;
      Emat(1,3,1) = E15;
	  
      return
      end subroutine buildEmat
	  
****************************************************************************
C**********************************************************************
          SUBROUTINE PUSHSV(SYM,VECT,IFLAG,NDI,NTENS)
          IMPLICIT REAL*8 (A-H,O-Z)

C        IFLAG=1   CONVERTS A SYMMETRIC MATRIX WRITTEN AS A
C                  VECTOR VECT(6) TO A MATRIX SYM(3,3)
C        IFLAG=2   CONVERTS A SYMMETRIC MATRIX SYM(3,3) TO 
C                  THE CORRESPONDING VECTOR VECT(6)
C ----------------------------------------------------------------------
         DIMENSION SYM(3,3),VECT(NTENS,1)
         NSHR=NTENS-NDI
         IF (IFLAG.EQ.1) THEN 
               CALL ZEROM(SYM)
               DO 15 I=1,NDI
15             SYM(I,I)=VECT(I,1)
               IF (NSHR.NE.0)    THEN
                   SYM(2,3)=VECT(NDI+1,1)
                   SYM(3,2)=VECT(NDI+1,1)
                        IF  (NSHR.NE.1)  THEN
                        SYM(1,3)=VECT(NDI+2,1)
                        SYM(3,1)=VECT(NDI+2,1)
                               IF  (NSHR.NE.2)  THEN
                               SYM(1,2)=VECT(NDI+3,1)
                               SYM(2,1)=VECT(NDI+3,1)
                               ENDIF
                        ENDIF
               ENDIF
        ELSE 
             IF (IFLAG.EQ.2) THEN 
                 DO 24 I=1,NTENS
24               VECT(I,1)=0.D0
                 DO 25 I=1,NDI
25               VECT(I,1)=SYM(I,I)
                 IF (NSHR.NE.0) THEN
                       VECT(NDI+1,1)=SYM(2,3)
                       IF (NSHR.NE.1) THEN
                             VECT(NDI+2,1)=SYM(1,3)
                             IF (NSHR.NE.2) VECT(NDI+3,1)=SYM(1,2)
                       ENDIF
                 ENDIF
             ELSE
             WRITE(6,*) '** ERROR IN PUSHSV WRONG IFLAG  '
             ENDIF
         ENDIF

         RETURN
         END
C**********************************************************************
          SUBROUTINE PUSHSVSTRAIN(SYM,VECT,IFLAG,NDI,NTENS)
          IMPLICIT REAL*8 (A-H,O-Z)

C 
C         This subroutine is the same as pushsv(), but accounts for
C         factors of two and half in converting between shear strains
C         in Voigt and tensor notation.
C 
C        IFLAG=1   CONVERTS A SYMMETRIC MATRIX WRITTEN AS A
C                  VECTOR VECT(6) TO A MATRIX SYM(3,3)
C        IFLAG=2   CONVERTS A SYMMETRIC MATRIX SYM(3,3) TO 
C                  THE CORRESPONDING VECTOR VECT(6)
C ----------------------------------------------------------------------
         DIMENSION SYM(3,3),VECT(NTENS,1)
		 
         real*8 two,half
         parameter(two=2.d0,half=0.5d0)
		 
         NSHR=NTENS-NDI
         IF (IFLAG.EQ.1) THEN 
               CALL ZEROM(SYM)
               DO 15 I=1,NDI
15             SYM(I,I)=VECT(I,1)
               IF (NSHR.NE.0)    THEN
                   SYM(2,3)=half*VECT(NDI+1,1)
                   SYM(3,2)=half*VECT(NDI+1,1)
                        IF  (NSHR.NE.1)  THEN
                        SYM(1,3)=half*VECT(NDI+2,1)
                        SYM(3,1)=half*VECT(NDI+2,1)
                               IF  (NSHR.NE.2)  THEN
                               SYM(1,2)=half*VECT(NDI+3,1)
                               SYM(2,1)=half*VECT(NDI+3,1)
                               ENDIF
                        ENDIF
               ENDIF
        ELSE 
             IF (IFLAG.EQ.2) THEN 
                 DO 24 I=1,NTENS
24               VECT(I,1)=0.D0
                 DO 25 I=1,NDI
25               VECT(I,1)=SYM(I,I)
                 IF (NSHR.NE.0) THEN
                       VECT(NDI+1,1)=two*SYM(2,3)
                       IF (NSHR.NE.1) THEN
                             VECT(NDI+2,1)=two*SYM(1,3)
                             IF (NSHR.NE.2) VECT(NDI+3,1)=two*SYM(1,2)
                       ENDIF
                 ENDIF
             ELSE
             WRITE(6,*) '** ERROR IN PUSHSVSTRAIN WRONG IFLAG  '
             ENDIF
         ENDIF

         RETURN
         END
C**********************************************************************  
C**********************************************************************
      	SUBROUTINE ZEROM(A)
C
C	THIS SUBROUTINE SETS ALL ENTRIES OF A 3 BY 3 MATRIX TO 0.D0.
C**********************************************************************

        REAL*8 A(3,3)

	DO 1 I=1,3
	  DO 1 J=1,3
	    A(I,J) = 0.D0
1	CONTINUE
C	
	RETURN
	END

C**********************************************************************
!****************************************************************************
!     Element subroutines
!****************************************************************************

      subroutine xint3D1pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using a 1 gauss point for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !      
      integer nIntPt,nDim
      !
      real*8 xi(1,3),w(1)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 1


      ! Gauss weights
      !
      w(1) = 8.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = 0.d0
      xi(1,2) = 0.d0
      xi(1,3) = 0.d0


      return
      end subroutine xint3D1pt
      
!************************************************************************

      subroutine xint3D8pt(xi,w,nIntPt)
      !
      ! This subroutine will get the integration point locations
      !  and corresponding gauss quadrature weights for 3D elements
      !  using 8 gauss points for integration
      !
      !  xi(nIntPt,3): xi,eta,zeta coordinates for the integration pts
      !  w(nIntPt):    corresponding integration weights
      !
      implicit none
      !
      integer nIntPt,nDim
      !
      real*8 xi(8,3),w(8)


      ! Initialize
      !
      w = 0.d0
      xi = 0.d0


      ! Number of Gauss points
      !
      nIntPt = 8


      ! Gauss weights
      !
      w(1) = 1.d0
      w(2) = 1.d0
      w(3) = 1.d0
      w(4) = 1.d0
      w(5) = 1.d0
      w(6) = 1.d0
      w(7) = 1.d0
      w(8) = 1.d0
      

      ! Gauss pt locations in master element
      !
      xi(1,1) = -dsqrt(1.d0/3.d0)
      xi(1,2) = -dsqrt(1.d0/3.d0)
      xi(1,3) = -dsqrt(1.d0/3.d0)
      xi(2,1) = dsqrt(1.d0/3.d0)
      xi(2,2) = -dsqrt(1.d0/3.d0)
      xi(2,3) = -dsqrt(1.d0/3.d0)
      xi(3,1) = -dsqrt(1.d0/3.d0)
      xi(3,2) = dsqrt(1.d0/3.d0)
      xi(3,3) = -dsqrt(1.d0/3.d0)
      xi(4,1) = dsqrt(1.d0/3.d0)
      xi(4,2) = dsqrt(1.d0/3.d0)
      xi(4,3) = -dsqrt(1.d0/3.d0)
      xi(5,1) = -dsqrt(1.d0/3.d0)
      xi(5,2) = -dsqrt(1.d0/3.d0)
      xi(5,3) = dsqrt(1.d0/3.d0)
      xi(6,1) = dsqrt(1.d0/3.d0)
      xi(6,2) = -dsqrt(1.d0/3.d0)
      xi(6,3) = dsqrt(1.d0/3.d0)
      xi(7,1) = -dsqrt(1.d0/3.d0)
      xi(7,2) = dsqrt(1.d0/3.d0)
      xi(7,3) = dsqrt(1.d0/3.d0)
      xi(8,1) = dsqrt(1.d0/3.d0)
      xi(8,2) = dsqrt(1.d0/3.d0)
      xi(8,3) = dsqrt(1.d0/3.d0)


      return
      end subroutine xint3D8pt

!************************************************************************

      subroutine calcShape3DLinear(nIntPt,xi_int,intpt,sh,dshxi)
      !
      ! Calculate the shape functions and their derivatives at the
      ! given integration point in the master element
      !
      ! sh(i) = shape function of node i at the intpt.
      ! dshxi(i,j) = derivative wrt j direction of shape fn of node i
      !
      implicit none
      !
      integer intpt,nDim,nIntPt,i,j
      !
      real*8 xi_int(nIntPt,3),sh(8),dshxi(8,3),xi,eta,zeta
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Location in the master element
      !
      xi = xi_int(intpt,1)
      eta = xi_int(intpt,2)
      zeta = xi_int(intpt,3)
      
      
      ! The shape functions
      !
      sh(1) = eighth*(one - xi)*(one - eta)*(one - zeta)
      sh(2) = eighth*(one + xi)*(one - eta)*(one - zeta)
      sh(3) = eighth*(one + xi)*(one + eta)*(one - zeta)
      sh(4) = eighth*(one - xi)*(one + eta)*(one - zeta)
      sh(5) = eighth*(one - xi)*(one - eta)*(one + zeta)
      sh(6) = eighth*(one + xi)*(one - eta)*(one + zeta)
      sh(7) = eighth*(one + xi)*(one + eta)*(one + zeta)
      sh(8) = eighth*(one - xi)*(one + eta)*(one + zeta)
      
      
      ! The first derivatives
      !
      dshxi(1,1) = -eighth*(one - eta)*(one - zeta)
      dshxi(1,2) = -eighth*(one - xi)*(one - zeta)
      dshxi(1,3) = -eighth*(one - xi)*(one - eta)
      dshxi(2,1) = eighth*(one - eta)*(one - zeta)
      dshxi(2,2) = -eighth*(one + xi)*(one - zeta)
      dshxi(2,3) = -eighth*(one + xi)*(one - eta)
      dshxi(3,1) = eighth*(one + eta)*(one - zeta)
      dshxi(3,2) = eighth*(one + xi)*(one - zeta)
      dshxi(3,3) = -eighth*(one + xi)*(one + eta)
      dshxi(4,1) = -eighth*(one + eta)*(one - zeta)
      dshxi(4,2) = eighth*(one - xi)*(one - zeta)
      dshxi(4,3) = -eighth*(one - xi)*(one + eta)
      dshxi(5,1) = -eighth*(one - eta)*(one + zeta)
      dshxi(5,2) = -eighth*(one - xi)*(one + zeta)
      dshxi(5,3) = eighth*(one - xi)*(one - eta)
      dshxi(6,1) = eighth*(one - eta)*(one + zeta)
      dshxi(6,2) = -eighth*(one + xi)*(one + zeta)
      dshxi(6,3) = eighth*(one + xi)*(one - eta)
      dshxi(7,1) = eighth*(one + eta)*(one + zeta)
      dshxi(7,2) = eighth*(one + xi)*(one + zeta)
      dshxi(7,3) = eighth*(one + xi)*(one + eta)
      dshxi(8,1) = -eighth*(one + eta)*(one + zeta)
      dshxi(8,2) = eighth*(one - xi)*(one + zeta)
      dshxi(8,3) = eighth*(one - xi)*(one + eta)
      
      
      return
      end subroutine calcShape3DLinear

!*************************************************************************

      subroutine mapShape3D(nNode,dshxi,coords,dsh,detMapJ,stat)
      !
      ! Map derivatives of shape fns from xi-eta-zeta domain
      !  to x-y-z domain.
      !
      implicit none
      !
      integer i,j,k,nNode,ieror,stat
      !
      real*8 dshxi(nNode,3),dsh(nNode,3),coords(3,nNode),mapJ(3,3),
     +  mapJ_inv(3,3),detmapJ
      !
      real*8 zero,one,two,half,fourth,eighth
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,fourth=0.25d0,
     +     eighth=1.d0/8.d0)
      

      ! Calculate the mapping Jacobian matrix:
      !
      mapJ = zero
      do i=1,3
        do j=1,3
          do k=1,nNode
              mapJ(i,j) = mapJ(i,j) + dshxi(k,i)*coords(j,k)
          end do
        end do
      end do


      ! Calculate the inverse and the determinant of Jacobian
      !
      call matInv3D(mapJ,mapJ_inv,detMapJ,stat)


      ! Calculate first derivatives wrt x, y, z
      !
      dsh = transpose(matmul(mapJ_inv,transpose(dshxi)))


      return
      end subroutine mapShape3D

!****************************************************************************
!     Utility subroutines
!****************************************************************************

      subroutine matInv3D(A,A_inv,det_A,istat)
      !
      ! Returns A_inv, the inverse and det_A, the determinant
      ! Note that the det is of the original matrix, not the
      ! inverse
      !
      implicit none
      !
      integer istat
      !
      real*8 A(3,3),A_inv(3,3),det_A,det_A_inv


      istat = 1
      
      det_A = A(1,1)*(A(2,2)*A(3,3) - A(3,2)*A(2,3)) -
     +        A(2,1)*(A(1,2)*A(3,3) - A(3,2)*A(1,3)) +
     +        A(3,1)*(A(1,2)*A(2,3) - A(2,2)*A(1,3))
      
      if (det_A .le. 0.d0) then
        write(*,*) 'WARNING: SUBROUTINE matInv3:'
        write(*,*) 'WARNING: DET of MAT=',DET_A
        istat = 0
        return
      end if
          
      det_A_inv = 1.d0/det_A
        
      A_inv(1,1) = det_A_inv*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
      A_inv(1,2) = det_A_inv*(A(3,2)*A(1,3)-A(1,2)*A(3,3))
      A_inv(1,3) = det_A_inv*(A(1,2)*A(2,3)-A(2,2)*A(1,3))
      A_inv(2,1) = det_A_inv*(A(3,1)*A(2,3)-A(2,1)*A(3,3))
      A_inv(2,2) = det_A_inv*(A(1,1)*A(3,3)-A(3,1)*A(1,3))
      A_inv(2,3) = det_A_inv*(A(2,1)*A(1,3)-A(1,1)*A(2,3))
      A_inv(3,1) = det_A_inv*(A(2,1)*A(3,2)-A(3,1)*A(2,2))
      A_inv(3,2) = det_A_inv*(A(3,1)*A(1,2)-A(1,1)*A(3,2))
      A_inv(3,3) = det_A_inv*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
      

      return
      end subroutine matInv3D

!****************************************************************************

      subroutine mdet(A,det)
      !
      ! This subroutine calculates the determinant
      ! of a 3 by 3 matrix [A]
      !
      implicit none
      !
      real*8  A(3,3),det


      det = A(1,1)*A(2,2)*A(3,3) 
     +	  + A(1,2)*A(2,3)*A(3,1)
     +	  + A(1,3)*A(2,1)*A(3,2)
     +	  - A(3,1)*A(2,2)*A(1,3)
     +	  - A(3,2)*A(2,3)*A(1,1)
     +	  - A(3,3)*A(2,1)*A(1,2)


      return
      end subroutine mdet
	
!****************************************************************************

      subroutine onem(A)
      !
      ! This subroutine stores the identity matrix in the
      ! 3 by 3 matrix [A]
      !
      implicit none
      !
      integer i,j
      !
      real*8 A(3,3)


      do i=1,3
         do J=1,3
	    if (i .eq. j) then
              A(i,j) = 1.0
            else
              A(i,j) = 0.0
            end if
         end do
      end do


      return
      end subroutine onem

!****************************************************************************
