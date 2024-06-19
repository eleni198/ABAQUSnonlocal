      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
!
      USE nlmodule	 

      IMPLICIT NONE   

      REAL(8) :: TIME(2),DTIME
      INTEGER :: LOP,LRESTART,KSTEP,KINC
!

! user coding to set up the FORTRAN environment, open files, close files, 
! calculate user-defined model-independent history information,
! write history information to external files,
! recover history information during restart analyses, etc.
! do not include calls to utility routine XIT
      
!----------------------------------------------------------------------     
!  The use of ABA_PARAM.INC eliminates the need to have different
!  versions of the code for single and double precision.  
!  ABA_PARAM.INC defines an appropriate IMPLICIT REAL statement and 
!  sets the value of NPRECD to 1 or 2, depending on whether the 
!  machine uses single or double precision.              
!----------------------------------------------------------------------     
!
!  Variables:  (all variables are passed in for information. No 
!               variables need to be defined.)
!  LOP
!
!    LOP=0 indicates that the subroutine is being called at the start of 
!          the analysis.
!    LOP=1 indicates that the subroutine is being called at the start of 
!          the current analysis increment. 
!    LOP=2 indicates that the subroutine is being called at the end of 
!          the current analysis increment. When LOP=2, all information 
!          that you need to restart the analysis should be written to 
!          external files.
!    LOP=3 indicates that the subroutine is being called at the end of 
!          the analysis.
!    LOP=4 indicates that the subroutine is being called at the beginning 
!          of a restart analysis. When LOP=4, all necessary external 
!          files should be opened and properly positioned and all 
!          information required for the restart should be read from 
!          the external files.
!  
!  LRESTART
!
!    LRESTART=0 indicates that an analysis restart file is not being 
!               written for this increment.
!    LRESTART=1 indicates that an analysis restart file is being 
!               written for this increment.
!    LRESTART=2 indicates that an analysis restart file is being 
!               written for this increment and that only one increment 
!               is being retained per step so that the current increment 
!               overwrites the previous increment in the restart file.
!
!  TIME =    Current (step,total) time.!    
!  DTIME =   Time increment    
!  KSTEP =   Current step number. When LOP=4, KSTEP gives the restart 
!                 step number.
!  KINC =    Current increment number. When LOP=4, KINC gives the 
!              restart increment number.
!  FNAME =   Root file (without extension) name of input file(s).
!  NRU =     Number of results files (.fil) to be read.
!  LRUNIT =  Array containing unit number and format of results files:
!              LRUNIT(1,*) --> Unit number of input file.
!              LRUNIT(2,*) --> Format of input file. (1-ascii, 2-binary)  
!              lrunit(1,*) == 8 --> result file (.fil) of the current
!                                   analysis
!  JUNIT =   Unit number of file to be opened.
!  JRCD =    Error check return code:
!              .EQ. 0 --> No errors.
!              .NE. 0 --> Errors detected.
!  KEY =     Current record key identifier.
!  OUTDIR =  Path of the current analysis directory

      character (108)    :: fModelSize,outdir
      integer (8)        :: nel,nip,nod,lenoutdir,k1,key,jrcd
      integer (8)        :: i,j,k,i1,j1
      integer :: INTV(500)   ! variable for the subroutine STDB_ABQERR
      real (8) :: REALV(500) ! variable for the subroutine STDB_ABQERR
      CHARACTER (8) :: CHARV(500) ! variable for the function STDB_ABQERR
	  double precision :: tmp2


!     starting with the size of the problem

!      CALL STDB_ABQERR(1,'uexternaldb called at '
!     &         // 'step=%I, increment=%I with lop=%I and lrestart=%I.',
!     &           (/kstep,kinc,lop,lrestart/),REALV,CHARV)
 
      if ((lop.eq.0)) then	  
         fModelSize = 'modelSize.txt'
         call getoutdir( outdir, lenoutdir )
         fModelSize = trim(outdir) // '\'  // 'modelSize.txt'   
         open(100, FILE=fmodelSize, ACCESS='SEQUENTIAL', 
     &          STATUS='OLD', FORM='FORMATTED')

         read(100,*) nel 
         read(100,*) nod
         read(100,*) k1
		 allocate(l_c(k1,2))
		 do i=1,k1
		     read(100,*) l_c(i,1),l_c(i,2)
		 enddo
		 allocate(materials(nel,8))		 
		 materials(:,:) = 0			 
         close(100)
		 allocation = .FALSE.
		 write(6,*) 'read materials'
      endif  
	  
	  if (allocated(TAB_IP).AND.(lop.eq.2)) then  
!     getting the size of the problem
      nel=size(TAB_IP,1)
      nip=size(TAB_IP,2)	  
!     updating the non local plastic strains
      do i=1,nel
       do j=1,nip	  
		if (TAB_IP(i,j)%in_set) then
		 TAB_IP(i,j)%ASDV = 0.0	
!        circling over the neighbors			 
         do k=1,TAB_IP(i,j)%num
		  i1 = TAB_IP(i,j)%nb(k,1)
		  j1 = TAB_IP(i,j)%nb(k,2)	
          tmp2 = TAB_IP(i,j)%ws(k)			 
		  TAB_IP(i,j)%ASDV(:)=TAB_IP(i,j)%ASDV(:)
     1 				 + tmp2*TAB_IP(i1,j1)%SDV(:)		  
         enddo	
         TAB_IP(i,j)%ASDV(:)=	
     1 		TAB_IP(i,j)%ASDV(:)/TAB_IP(i,j)%wtot		 		 
		endif
	   enddo
	  enddo
 	  endif	

      RETURN
      END