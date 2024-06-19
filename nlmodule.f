      MODULE nlmodule
      implicit none				
      type  ::  T_GIP  
         integer                               ::  num,mat     
         integer,allocatable,dimension(:,:)    ::  nb 		 
         real(kind=8),allocatable,dimension(:) ::  ws 
         real(kind=8)                          ::  l_c 
         real(kind=8),dimension(3)             ::  coord
         real(kind=8),allocatable,dimension(:) ::  SDV,ASDV
         real(kind=8)                          ::  wtot
         logical                               ::  in_set   
         character(len=5)                      ::  behav              
      end type T_GIP

      type(T_GIP),allocatable,dimension(:,:),public :: TAB_IP
	  
      double precision,allocatable,dimension(:,:),public :: l_c 
      integer,allocatable,dimension(:,:),public :: materials
	  
      logical,public :: allocation = .FALSE. 
	  
	  contains
	  
      SUBROUTINE get_nl(ASDV,PROPS,NSTATEV,NPROPS,NOEL,NPT)
! 	  
      IMPLICIT NONE
!
!--------------------------------------------------------------- 
!     Declarating UMAT variables and constants 
      INTEGER NSTATEV,NPROPS,NOEL,NPT
      DOUBLE PRECISION PROPS(NPROPS),ASDV(NSTATEV)
! ---------------------------------------------------------------      
!     what get_nl does, depends on where it is called
! ---------------------------------------------------------------  
!     Nonlocal called when TAB_IP has not been allocated
      if ((.NOT.allocated(TAB_IP))) then
	     materials(NOEL,NPT) = INT(PROPS(13))
		 ASDV(:) = 0.0 
!         goto 10 
      else	     
	     ASDV(:) = TAB_IP(NOEL,NPT)%ASDV(:)
      endif	 
! ---------------------------------------------------------------
! 10   continue
! --------------------------------------------------------------- 
      RETURN
      END SUBROUTINE get_nl
!	  
!	  
      SUBROUTINE give_nl(STATEV,NSTATEV,NOEL,NPT)
! 	  
      IMPLICIT NONE
!
!--------------------------------------------------------------- 
!     Declarating UMAT variables and constants 
      INTEGER NSTATEV,NOEL,NPT,J
      DOUBLE PRECISION STATEV(NSTATEV)
! ---------------------------------------------------------------      
!     what get_nl does, depends on where it is called
! ---------------------------------------------------------------  
!     Nonlocal called when TAB_IP has not been allocated does nothing
      if (( allocated(TAB_IP))) then
	     
	     TAB_IP(NOEL,NPT)%SDV(:) = STATEV(:)
		 
!		 STATEV(25) = TAB_IP(NOEL,NPT)%ASDV(1)		 
!		 STATEV(15) = TAB_IP(NOEL,NPT)%ASDV(14) 

      endif	 
! ---------------------------------------------------------------
! 20   continue
! --------------------------------------------------------------- 
      END SUBROUTINE give_nl
	  
      SUBROUTINE WEIGHTS(ra,rb,rc)
 	  
      IMPLICIT NONE
	  
	  REAL(8) :: ra,rb,rc

!--------------------------------------------------------------- 
      ra = 15.0/(16.0*rc)*(1.0-(rb/rc)**2)**2
! --------------------------------------------------------------- 
      END SUBROUTINE WEIGHTS	  
	  
      END MODULE nlmodule