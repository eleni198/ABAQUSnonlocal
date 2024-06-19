      SUBROUTINE URDFIL(LSTOP,LOVRWRT,KSTEP,KINC,DTIME,TIME)
!
      USE nlmodule	 
!
      IMPLICIT NONE  
!

      REAL(8) :: TIME(2),DTIME
      INTEGER :: LSTOP,LOVRWRT,JRCD,KSTEP,KINC,NPRECD
      parameter (nprecd=2)	  
      integer, parameter :: blocksize = 513
      real (8) :: ARRAY(blockSize),tmp1(3),tmp2,lc
      character (8) :: crray(blockSize)
      INTEGER :: JRRAY(NPRECD,513)
      EQUIVALENCE (ARRAY(1),JRRAY(1,1)), (ARRAY(1), CRRAY(1))

      character (108)    :: fModelSize,outdir
      integer        :: nel,nip,nod,lenoutdir,k1,key,i1,j1,i,j,k,m1

      integer,allocatable,dimension(:)   :: ipxel
      integer,allocatable,dimension(:,:) :: conn
      REAL(8),allocatable,dimension(:,:) :: rconn

	  
      if ((.NOT.allocated(TAB_IP))) then

         fModelSize = 'modelSize.txt'
         call getoutdir( outdir, lenoutdir )
         fModelSize = trim(outdir) // '\'  // 'modelSize.txt'   
         open(100, FILE=fmodelSize, ACCESS='SEQUENTIAL', 
     &          STATUS='OLD', FORM='FORMATTED')

         read(100,*) nel 
         read(100,*) nod

         close(100)
!        to begin with, we look for the number of integration points of each element
         allocate(ipxel(nel))

		 ipxel(:) = 0
         
         DO K1=1,999999
             CALL DBFILE(0,JRRAY,JRCD)
             IF (JRCD .NE. 0) EXIT
             KEY=JRRAY(1,2) 
			 if (KEY.eq.1) then
			     ipxel(JRRAY(1,3))=max0(ipxel(JRRAY(1,3)),JRRAY(1,4))
             endif        
         END DO
		 


!        for simplicity we just allocate the maximum length
         nip = MAXVAL(ipxel, 1)
!        allocating
         allocate(TAB_IP(nel,nip))
!        allocating the connectivity matrix
         allocate(conn(nel*nip,2))
         allocate(rconn(nel,nip))
		 conn  = 0
		 rconn = 0.0
		 TAB_IP(:,:)%in_set = .FALSE.
!		 
!        rewinding and fishing for the coordinates		 
         CALL DBFILE(2,ARRAY,JRCD)

         DO K1=1,999999
             CALL DBFILE(0,ARRAY,JRCD)
             IF (JRCD .NE. 0) EXIT
             KEY=JRRAY(1,2)         
             IF (KEY.EQ.1) THEN
			     i=JRRAY(1,3) 
			     j=JRRAY(1,4)

             END IF
             IF (KEY.EQ.8) THEN	 
                 TAB_IP(i,j)%coord(1:3) = ARRAY(3:5)
                 TAB_IP(i,j)%mat = materials(i,j)	
                 TAB_IP(i,j)%l_c = l_c(TAB_IP(i,j)%mat,1)	
                 j1=int(l_c(TAB_IP(i,j)%mat,2)+1e-6) 
				 allocate(TAB_IP(i,j)%SDV(J1))
				 TAB_IP(i,j)%SDV=0.0
				 allocate(TAB_IP(i,j)%ASDV(J1))	
				 TAB_IP(i,j)%ASDV=0.0				 
				 TAB_IP(i,j)%in_set = .TRUE.
             END IF			 
         END DO
      
 
 
 
!        it's time to look for the neighbors	
         do i=1,nel
         do j=1,nip
!        only existing nodes can have neighbors
	     if (TAB_IP(i,j)%in_set) then
		     conn = 0
			 rconn = 0.0
			 k1=0
			 lc=TAB_IP(i,j)%l_c
			 m1=TAB_IP(i,j)%mat			 
		     do i1=1,nel 
		     do j1=1,nip 
!            only existing nodes can have neighbors
	         if ((TAB_IP(i1,j1)%in_set).AND.(TAB_IP(i1,j1)%mat.eq.m1)) then
			     tmp1 =TAB_IP(i1,j1)%coord - TAB_IP(i,j)%coord 
				 tmp2 =sqrt(tmp1(1)**2 + tmp1(2)**2 + tmp1(3)**2)
				 if (tmp2.le.lc) then
				    k1=k1+1
				    conn(k1,1) = i1
				    conn(k1,2) = j1		
                    rconn(i1,j1) = tmp2					
				 endif
             endif			
		     enddo					 
		     enddo
!            number of neighbors is k1
             TAB_IP(i,j)%num = k1
!            initializing sum of weights	
             tmp2 = 0.0		 
             allocate(TAB_IP(i,j)%nb(k1,2))	
             allocate(TAB_IP(i,j)%ws(k1))				 
!            saving the array of neighbors
             do i1=1,k1
			     TAB_IP(i,j)%nb(i1,1)=conn(i1,1)
			     TAB_IP(i,j)%nb(i1,2)=conn(i1,2)
				 call weights(TAB_IP(i,j)%ws(i1),rconn(conn(i1,1),conn(i1,2)),lc)		
                 tmp2=tmp2+TAB_IP(i,j)%ws(i1)
             enddo	
			 TAB_IP(i,j)%wtot=tmp2
!            having established the connectivity, we move on to evaluate some quantities	
!            initializing the state variables
!            local
			 TAB_IP(i,j)%SDV  = 0.0	
!            nonlocal			 
			 TAB_IP(i,j)%ASDV = 0.0		
!			 
		 endif
         enddo
         enddo
		 
         deallocate(conn)
         deallocate(rconn)
         deallocate(materials)	 

		 return
      endif 

	  
 200  RETURN
      END
	  
	  
	  
	  
	  
