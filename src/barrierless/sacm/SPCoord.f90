program SPCoord

implicit none

integer							::	NAt,N,i,j,N_Fragx,N_Fragy,N_Fragz,k,l
double precision,allocatable	::	q(:,:),m(:),q_x(:,:),m_x(:),q_z(:,:)
double precision,allocatable	::	q_y(:,:),m_y(:)
double precision,allocatable	::  m_z(:)
double precision				::	x(0:3),y(0:3),z(0:3),mu
double precision				::	mx_ges,my_ges,mz_ges,q2e,q2(0:3),q1(0:3),qe

open(unit=21,file='SPCoord.inp')
open(unit=22,file='SPCoord.out')
!open(unit=7,file='SPCoord2.out')

write(*,*) 'triatomic (3) or diatomic (2) ?'
read(*,*) NAt




!------------------------------------- Einlesen von SPCoord.inp --------------------------------------------------------------------------------------

read(21,*) N

allocate(q(N,3))
allocate(m(N))

do i=1,N
	read(21,*) m(i),q(i,1),q(i,2),q(i,3)	
end do

write(22,*) '---------------- Coordinates of the fragments -------------'

write(22,"(/,X,A11)") 'Fragment X:'	

read(21,*) N_Fragx
allocate(q_x(N_Fragx,3))
allocate(m_x(N_Fragx))

do i=1,N_Fragx
	read(21,*) k
	m_x(i)=m(k)
		do j=1,3
			q_x(i,j)=q(k,j)
		end do
	write(22,"(F10.1,3F14.8)") m_x(i),q_x(i,1),q_x(i,2),q_x(i,3)
end do

write(22,"(/,X,A11)") 'Fragment Y:'

read(21,*) N_Fragy
allocate(q_y(N_Fragy,3))
allocate(m_y(N_Fragy))

do i=1,N_Fragy
	read(21,*) k
	m_y(i)=m(k)
		do j=1,3
			q_y(i,j)=q(k,j)
		end do
	write(22,"(F10.1,3F14.8)") m_y(i),q_y(i,1),q_y(i,2),q_y(i,3)		
end do

if(NAt .gt. 2) Then

	write(22,"(/,X,A11)") 'Fragment Z:'

	read(21,*) N_Fragz
	allocate(q_z(N_Fragz,3))
	allocate(m_z(N_Fragz))

	do i=1,N_Fragz
		read(21,*) k
		m_z(i)=m(k)
			do j=1,3
				q_z(i,j)=q(k,j)
			end do
		write(22,"(F10.1,3F14.8)") m_z(i),q_z(i,1),q_z(i,2),q_z(i,3)
	end do
	
end if


!-----------------------------------------------Berechnung der SP-Koordinaten----------------------------------------------------------------------------------------

write(22,"(A4,/,X)")
write(22,*) '--------------- CMS-coordinates of the fragments ----------'

!Fragment X:

mx_ges=0.0
do i=1,N_Fragx
	mx_ges=mx_ges+m_x(i)
end do

do j=1,3
	x(j)=0.0
	do i=1,N_Fragx
		x(j)=x(j)+m_x(i)*q_x(i,j)
	end do
	x(j)=x(j)/mx_ges
end do

write(22,"(/,X,A11,3F14.8)") "Fragment X:",x(1),x(2),x(3)	

!Fragment Y:

my_ges=0.0
do i=1,N_Fragy
	my_ges=my_ges+m_y(i)
end do

do j=1,3
	y(j)=0.0
	do i=1,N_Fragy
		y(j)=y(j)+m_y(i)*q_y(i,j)
	end do
	y(j)=y(j)/my_ges
end do

write(22,"(X,A11,3F14.8)") "Fragment Y:",y(1),y(2),y(3)



!Fragment Z:

!allocate(q_z(1000,1000))

if(NAt .gt. 2)Then

mz_ges=0.0
do i=1,N_Fragz
	mz_ges=mz_ges+m_z(i)
end do

do j=1,3
	z(j)=0.0
	do i=1,N_Fragz
		z(j)=z(j)+m_z(i)*q_z(i,j)
	end do
	z(j)=z(j)/mz_ges
end do

!write(7,*) mz_ges

write(22,"(X,A11,3F14.8)") "Fragment Z:",z(1),z(2),z(3)

Else

 continue

End if

do i=1,3
	q2(i)=x(i)-y(i)
end do

q2e=0.0
do i=1,3
	q2e=q2e+q2(i)**2.0
end do

q2e=dsqrt(q2e)





if(NAt .gt. 2)Then

	do i=1,3
		q1(i)=z(i)-y(i)
	end do
	
	qe=0.0
	do i=1,3
		qe=qe+q1(i)**2.0
		mu=mu+q1(i)*q2(i)
	end do
	qe=dsqrt(qe)
	mu=mu/(qe*q2e)
	mu=acos(mu)/3.14159*180.0
	
	write(22,"(/,X,A20,3F14.8)") "qe  [in Angstroem] =",qe
	write(22,"(X,A20,3F14.8)") "q2e [in Angstroem] =",q2e
	write(22,"(X,A20,3F14.8)") "v       [in °]     =",mu
	
Else
	
	write(22,"(/,X,A20,3F14.8)") "qe  [in Angstroem] =",q2e

End if

Pause
	
end program SPCoord
	
