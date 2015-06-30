
module initinfo

!integer, save             :: nbases =  10
integer, save             :: nbases
integer                   :: nsteps
integer                   :: natoms
integer,allocatable       :: nheavyatoms(:)
character*40              :: mdcrdfile
integer, allocatable      :: heavyatoms(:,:)                       !atom numbers for the C5 and C6 of each thymine which is part of a ttstep
integer, allocatable      :: check(:,:)                            !atom numbers for the atoms to make sure principal axis have correct sign
character*40              :: twistfile                             !twist output file
character*40              :: rollfile                              !roll output file
character*40              :: tiltfile                              !tilt output file
!character*40, save :: twistfile = 'twist.dat'      !dihedral output file
!character*40, save :: rollfile = 'roll.dat'      !dihedral output file
!character*40, save :: tiltfile = 'tilt.dat'      !dihedral output file
character*40, save :: cfgfile = '.compute_helical_angles.cfg'      !
character*40, save :: heavyatomfile = '.heavyatoms'                !
character*40, save :: checkatomfile = '.checkatoms'                !

endmodule initinfo

module runtime
integer step                   !current step being analyzed
!real, allocatable :: d(:)                 !distance between C5-C6 midpoints 
!real, allocatable :: dih(:)               !dihedral C5-C6-C5*-C6*
real, allocatable :: pos(:,:)             !positions of all atoms at current step
!integer atoms(2,2)             !Stores atom numbers for C5-C6-C5*-C6*
integer counter
integer trial

endmodule runtime

program ttinfo
use initinfo
use runtime
implicit none
real avgd, avgdih, stdevd, stdevdih
integer i

call seconstants ()


!open(30, file=dihfile,form='binary')
open(30, file=twistfile)
open(40, file=rollfile)
open(50, file=tiltfile)
counter = 0

   open(10,file=mdcrdfile)
   read(10,*)
   allocate(pos(3,natoms))
   do step = 1, nsteps
      counter = counter + 1

      call get_amber_coord()

      call compute_dihedral()

   enddo
   close(10)

close(30)
close(40)
close(50)

endprogram ttinfo

subroutine seconstants()
   use initinfo
   use runtime
   implicit none
   character numb1*1, numb2*2
   integer i,j

   open(13,file=cfgfile)
      read(13,'(i10)') natoms
      read(13,'(i10)') nsteps
      read(13,'(i10)') nbases
      read(13,'(a30)') mdcrdfile
      read(13,'(a30)') twistfile
      read(13,'(a30)') rollfile
      read(13,'(a30)') tiltfile
   close(13)
 
   allocate(check(nbases,3),nheavyatoms(nbases))

   open(14,file=heavyatomfile)
   open(15,file=checkatomfile)

      do i=1,nbases
         read(14,'(i5)') nheavyatoms(i)
         read(15,'(3i5)') (check(i,j),j=1,3)
      enddo
 
   close(14)
   close(15)

   allocate(heavyatoms(nbases,maxval(nheavyatoms)))

   open(16,file=heavyatomfile,position='rewind')
     do i=1,nbases
         read(16,'(5x,25i5)') (heavyatoms(i,j),j=1,nheavyatoms(i))
      enddo
   close(16)

!   atoms(1,1) =  11!N1 on cyt 1
!   atoms(1,2) =  22!O2 on cyt 1
!   atoms(2,1) =  41!N1 on cyt 2
!   atoms(2,2) =  52!O2 on cyt 2
!   atoms(3,1) =  71!N1 on cyt 3
!   atoms(3,2) =  82!O2 on cyt 3
!   atoms(4,1) = 107!N1 on cp 4
!   atoms(4,2) = 123!O2 on cp 4
!   atoms(5,1) = 142!N1 on cp 5
!   atoms(5,2) = 158!O2 on cp 5
!   atoms(6,1) = 171!N1 on cyt 6
!   atoms(6,2) = 182!02 on cyt 6
!   atoms(7,1) = 201!N1 on cyt 7
!   atoms(7,2) = 212!O2 on cyt 7
!   atoms(8,1) = 231!N1 on cyt 8
!   atoms(8,2) = 242!O2 on cyt 8
!   atoms(9,1) = 261!N1 on cyt 9
!   atoms(9,2) = 272!O2 on cyt 9
!   atoms(10,1)= 291!N1 on cyt 10
!   atoms(10,2)= 302!O2 on cyt 10

!   check(1,2) =  20!N3  on cyt 1
!   check(1,3) =  17!N4 on cyt 1
!   check(1,1) = atoms(2,1)!N1  on base above cyt 1
!   check(2,2) =  50!N3  on cyt 2
!   check(2,3) =  47!N4 on cyt 2
!   check(2,1) = atoms(3,1)!N1  on base above cyt 2
!   check(3,2) =  80!N3  on cyt 1
!   check(3,3) =  77!N4 on cyt 1
!   check(3,1) = atoms(4,1)!N1  on base above cyt 1
!   check(4,2) = 119!H9  on cp 4 
!   check(4,3) = 114!C8M on cp 4
!   check(4,1) = atoms(5,1)!N1  on base above cp 4
!   check(5,2) = 154!H9  on cp 5
!   check(5,3) = 149!C8M on cp 5
!   check(5,1) = atoms(6,1)!N1 on base above cp 5
!   check(6,2) = 180!N3  on cyt 6
!   check(6,3) = 177!N4 on cyt 6
!   check(6,1) = atoms(7,1)!N1  on base above cyt 6
!   check(7,2) = 210!N3  on cyt 7
!   check(7,3) = 207!N4 on cyt 7
!   check(7,1) = atoms(8,1)!N1  on base above cyt 7
!   check(8,2) = 240!N3  on cyt 8
!   check(8,3) = 237!N4 on cyt 8
!   check(8,1) = atoms(9,1)!N1  on base above cyt 8
!   check(9,2) = 270!N3  on cyt 9
!   check(9,3) = 267!N4 on cyt 9 
!   check(9,1) = atoms(10,1)!N1  on base above cyt 9
!   check(10,2)= 300!N3  on cyt 10
!   check(10,3)= 297!N4 on cyt 10
!   check(10,1)= 305!C2'  on base above cyt 10
!
!   natoms = 10121 
!   mdcrdfile =  'md2.mdcrd'
!   nsteps = 1500
!
endsubroutine seconstants


subroutine open_namd_traj (mdcrdfile, nsteps)
   implicit none
   integer junk, i, nsteps
   character*4 hdr 
   character*30 mdcrdfile


   open(10,file=mdcrdfile,form='binary')
   read(10) junk
   read(10) hdr
   read(10) nsteps

   do i=1,66
      read(10) junk
   enddo
 

endsubroutine open_namd_traj

subroutine get_namd_coord()
   use initinfo, only : natoms
   use runtime,  only : pos, trial
   implicit none
   integer i, j, current
   real temp
   integer junk
   real (kind=8) junkd, axis(3)
   real junkr

   read(10) junk
   read(10) axis(1)
   read(10) junkd
   read(10) axis(2)
   read(10) junkd
   read(10) junkd
   read(10) axis(3)
   read(10) junk
   do j=1,3
      read(10) junkr
      do i=1, natoms
         read(10) pos(j,i)
      enddo
      read(10) junkr
   enddo

endsubroutine get_namd_coord


subroutine get_amber_coord()
use initinfo, only : natoms
use runtime,  only : pos, trial
implicit none
integer i, j, counter

counter = 0
do i=1, natoms

   do j=1,3
      counter = counter + 1
      if (counter.gt.1.and.mod(counter,10).eq.1) read(10,*)
      read(10,'(F8.3)', advance='no') pos(j,i)
   enddo

enddo
read(10,*)
read(10,*)

endsubroutine get_amber_coord

subroutine compute_distance()
use initinfo, only : heavyatoms, nbases
use runtime, only : pos, step, counter
implicit none
integer i, j, base1, base2
real midpoint(3,2), tempdist

do base1=1, nbases-1
!   do base2 = base1+1, nbases
      base2 = base1+1
      do j=1,3
         midpoint(j,1) = (pos(j,heavyatoms(base1,1)) + pos(j,heavyatoms(base1,2)))/2.0
      enddo
      do j=1,3
         midpoint(j,2) = (pos(j,heavyatoms(base2,1)) + pos(j,heavyatoms(base2,2)))/2.0
      enddo

      tempdist = 0

      do j=1,3
         tempdist = tempdist + (midpoint(j,1)-midpoint(j,2))*(midpoint(j,1)-midpoint(j,2))
      enddo

   !d(counter) = sqrt(tempdist)
      tempdist = sqrt(tempdist)
      write(20) tempdist
!   write(20,'(F10.5)',advance='no') tempdist
!   enddo
enddo
!write(20,*)

endsubroutine compute_distance

subroutine compute_dihedral()
use initinfo, only : heavyatoms, nbases, check, nheavyatoms
use runtime, only : pos, step, counter
implicit none
integer, parameter :: lwork =88
real eigenvalues(3)
real work(lwork)
real, dimension(3,3) :: basis1, basis2
integer j, base1, base2
real dih
integer ipiv(3)
integer info
integer base1_atoms
integer base2_atoms
real, allocatable :: base1_pos(:,:), base1_pos_T(:,:)
real, allocatable :: base2_pos(:,:), base2_pos_T(:,:)
real theta
real base1_avg(3)
real base2_avg(3)
real, dimension(3) :: x_hat, y_hat, z_hat
real mag_x_hat, mag_y_hat, mag_z_hat
real sintheta
real costheta
real mag
integer atomcount, i
real, dimension(2) :: vec_1, vec_2, vec_3

do base1=1,nbases-1
   base2 = base1+1

   !determine the number of atoms in this base and allocate the array to store the atom positions
   base1_atoms = nheavyatoms(base1)
   allocate(base1_pos(3,base1_atoms),base1_pos_T(base1_atoms,3))
   base2_atoms = nheavyatoms(base2)
   allocate(base2_pos(3,base2_atoms),base2_pos_T(base2_atoms,3))

   !now get the base1 atom positions
   base1_avg=0
   do i=1,base1_atoms
      atomcount=heavyatoms(base1,i)
      base1_pos(:,i)=pos(:,atomcount)
      base1_avg(:) = base1_avg(:)+base1_pos(:,i)
   enddo
   base1_avg = base1_avg/real(base1_atoms)
   do i=1,base1_atoms
      base1_pos(:,i)=base1_pos(:,i)-base1_avg(:)
      base1_pos_T(i,:)=base1_pos(:,i)
   enddo
   !now get the base2 atom positions
   base2_avg=0
   do i=1,base2_atoms
      atomcount=heavyatoms(base2,i)
      base2_pos(:,i)=pos(:,atomcount)
      base2_avg(:) = base2_avg(:)+base2_pos(:,i)
   enddo
   base2_avg = base2_avg/real(base2_atoms)
   do i=1,base2_atoms
      base2_pos(:,i)=base2_pos(:,i)-base2_avg(:)
      base2_pos_T(i,:)=base2_pos(:,i)
   enddo


   !determine the axis for the base by computing the eigenvectors of the position*position_T matrix
   basis1 = matmul(base1_pos,base1_pos_T)
   call ssyev('V', 'L', 3, basis1, 3, eigenvalues, work, lwork, info )
   basis2 = matmul(base2_pos,base2_pos_T)
   call ssyev('V', 'L', 3, basis2, 3, eigenvalues, work, lwork, info )

   !make sure the basis vectors are pointed in the correct direction
   do i=1,3
      mag = ( pos(1,check(base1,i))-base1_avg(1) )*basis1(1,i) + ( pos(2,check(base1,i))-base1_avg(2) )*basis1(2,i) + ( pos(3,check(base1,i))-base1_avg(3) )*basis1(3,i) 
      if (mag.lt.0) then
         basis1(:,i)=-basis1(:,i)
      endif
   enddo
   do i=1,3
      mag = ( pos(1,check(base2,i))-base2_avg(1) )*basis2(1,i) + ( pos(2,check(base2,i))-base2_avg(2) )*basis2(2,i) + ( pos(3,check(base2,i))-base2_avg(3) )*basis2(3,i) 
      if (mag.lt.0) then
         basis2(:,i)=-basis2(:,i)
      endif
   enddo

   !twist:
   x_hat=dot_product(basis1(:,2),basis2(:,2))*basis1(:,2) + dot_product(basis1(:,3),basis2(:,2))*basis1(:,3)
   mag_x_hat=sqrt(dot_product(x_hat,x_hat)) 
   costheta=dot_product(x_hat,basis1(:,2))/mag_x_hat
   sintheta=dot_product(x_hat,basis1(:,3))/mag_x_hat
   theta=atan2(sintheta,costheta)*180.0/3.1415926535 
   write(30,'(f8.3)', advance='no') theta
!   theta=acos(costheta)*180.0/3.1415926535
!   write(30,'(f8.3)', advance='no') theta

   !tilt:
   y_hat=dot_product(basis1(:,3),basis2(:,3))*basis1(:,3) + dot_product(basis1(:,1),basis2(:,3))*basis1(:,1)
   mag_y_hat=sqrt(dot_product(y_hat,y_hat)) 
   costheta=dot_product(y_hat,basis1(:,3))/mag_y_hat
   sintheta=dot_product(y_hat,basis1(:,1))/mag_y_hat
   theta=atan2(sintheta,costheta)*180.0/3.1415926535 
   write(50,'(f8.3)', advance='no') theta

   !roll:
   z_hat=dot_product(basis1(:,2),basis2(:,1))*basis1(:,2) + dot_product(basis1(:,1),basis2(:,1))*basis1(:,1)
   mag_z_hat=sqrt(dot_product(z_hat,z_hat)) 
   costheta=dot_product(z_hat,basis1(:,1))/mag_z_hat
   sintheta=dot_product(z_hat,basis1(:,2))/mag_z_hat
   theta=-atan2(sintheta,costheta)*180.0/3.1415926535 
   write(40,'(f8.3)', advance='no') theta


   deallocate(base1_pos, base1_pos_T, base2_pos,base2_pos_T)
enddo
write(30,*)
write(40,*)
write(50,*)
endsubroutine compute_dihedral

subroutine stdev (array,arraysize,arrayavg,stdevout)
implicit none
integer arraysize, i, j
real array(arraysize), arrayavg, stdevout

stdevout =0
do i=1, arraysize
   stdevout = stdevout + (array(i)-arrayavg)*(array(i)-arrayavg)
enddo

stdevout = sqrt(stdevout/real(arraysize))

endsubroutine stdev


Subroutine compdih (theta, r1, r2, r3, r4)
IMPLICIT NONE
REAL, PARAMETER :: pi=3.1415926535
REAL theta, cos_theta, sin_theta
REAL, DIMENSION(1:3) :: r12, r23, r34, A, B, C, r1, r2, r3, r4
INTEGER i
do i=1,3
   r12(i)=r1(i)-r2(i)
   r23(i)=r2(i)-r3(i)
   r34(i)=r3(i)-r4(i)
enddo

CALL crossproduct(r12, r23, A)  !A=r12xr23/|r12xr23|
CALL crossproduct(r23, r34, B)  !B=r23xr34/|r23xr34|
CALL crossproduct(r23, A, C)
cos_theta=DOT_PRODUCT(A,B)
sin_theta=DOT_PRODUCT(C,B)
theta=-atan2(sin_theta, cos_theta)

theta = theta*180.0/pi

return
ENDSubroutine compdih

Subroutine crossproduct (ab, bc, cross)
!Computes and normalizes the cross product between vectos ab and bc
IMPLICIT NONE
REAL, DIMENSION(1:3) :: ab, bc, cross
REAL mag
INTEGER i
cross(1)=ab(2)*bc(3)-ab(3)*bc(2)
cross(2)=-ab(1)*bc(3)+ab(3)*bc(1)
cross(3)=ab(1)*bc(2)-ab(2)*bc(1)
mag=SQRT(cross(1)**2+cross(2)**2+cross(3)**2)
DO i=1,3
   cross(i)=cross(i)/mag
ENDDO
RETURN
ENDSubroutine crossproduct
