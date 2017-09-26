module matr_mod
contains
subroutine rand_int(matrix)
 integer, dimension(:,:),intent(inout) :: matrix
 real :: u
integer :: n,m

n = size(matrix,1)

do i = 1,n
do j=1,n
call RANDOM_NUMBER(u)
matrix(i,j) = FLOOR((100)*u)
enddo
enddo
end subroutine rand_int

function mul_int(ai,bi) result (ri)
 integer :: length
 integer, dimension(:,:),allocatable:: ri
 integer, dimension(:,:),intent(in)::ai,bi
length = size(ai,1)
allocate(ri(length,length))
ri = 0
 do i=1,length
 do j=1,length
do k=1,length
  ri(j,i) = ri(j,i) + ai(k,i)*bi(j,k)
enddo
enddo
enddo
end function mul_int

function mul_real(arr,brr) result(rr)
 integer :: length
real, dimension(:,:),allocatable:: rr
real, dimension(:,:),intent(in)::arr,brr
length = size(arr,1)
allocate(rr(length,length))
rr = 0
do i=1,length
do j=1,length
do k=1,length
 rr(i,j) = rr(i,j) + arr(i,k)*brr(k,j)
enddo
enddo
enddo
end function mul_real

function matr_add(ai,bi) result(ri)
  integer :: length
integer, dimension(:,:),allocatable:: ri
integer, dimension(:,:),intent(in)::ai,bi
length = size(ai,1)
allocate(ri(length,length))
ri = 0
do j=1,length
do i=1,length
 ri(i,j) =  ai(i,j)+bi(i,j)
enddo
enddo
end function matr_add

function matr_sub(ai,bi) result(ri)
  integer :: length
integer, dimension(:,:),allocatable:: ri
integer, dimension(:,:),intent(in)::ai,bi
length = size(ai,1)
allocate(ri(length,length))
ri = 0
do j=1,length
do i=1,length
 ri(i,j) =  ai(i,j)-bi(i,j)
enddo
enddo
end function matr_sub

function matr_transpose(ai) result(ri)
  integer :: length
integer, dimension(:,:),allocatable:: ri
integer, dimension(:,:),intent(in)::ai
length = size(ai,1)
allocate(ri(length,length))
ri = 0
do j=1,length
do i=1,length
 ri(i,j) =  ai(j,i)
enddo
enddo
end function matr_transpose

subroutine print_matrix(ai)
 integer :: length
integer, dimension(:,:),intent(in)::ai
length = size(ai,1)
do i=1,length
print *,(ai(j,i), j=1,length)
enddo
end subroutine print_matrix

subroutine print_real_matrix(ar)
 integer :: length
real, dimension(:,:),intent(in)::ar
length = size(ar,1)
do i=1,length
print *,(ar(j,i), j=1,length)
enddo
end subroutine print_real_matrix

function matr_minor(ai,col,row) result(minor)
  integer :: length,curRow,curCol
integer, dimension(:,:),allocatable:: minor
integer, dimension(:,:),intent(in)::ai
integer,intent(in) ::col,row
length = size(ai,1)
allocate(minor(length-1,length-1))
minor  =0;
curRow= 0;
do j=1,length
curCol=1
if (j==row) then
cycle
else
curRow=curRow+1
endif
do i=1,length
if (.NOT. i==col) then
 minor(curCol,curRow) = ai(i,j)
 curCol= curCol + 1
endif
enddo
enddo
end function matr_minor

recursive function matr_det(ai) result(det)
integer, dimension(:,:),intent(in)::ai
integer, allocatable,dimension(:,:) :: minor
integer :: det
length = size(ai,1)
if (length == 1) then
	det = ai(1,1)
else
det = 0
do i=1,length
allocate(minor(length-1,length-1))
minor = matr_minor(ai,i,1);
det = det+((-1)**(i+1))*ai(i,1)*matr_det(minor)
deallocate(minor)
enddo
endif
end function matr_det 

function matr_divider(ai,div) result (ri)
integer, dimension(:,:),intent(in)::ai
real ,intent(in):: div
real, allocatable,dimension(:,:) :: ri
length = size(ai,1)
allocate(ri(length,length))
do i=1,length
do j=1,length
ri(j,i) = ai(j,i)/div
enddo
enddo

end function matr_divider

function div_matr(ai,bi) result (ri)
integer, dimension(:,:),intent(in)::ai
integer, dimension(:,:),intent(in)::bi
real, allocatable,dimension(:,:) :: ri
length = size(ai,1)
allocate(ri(length,length))
ri = mul_check_inversion(ai,matr_inversion(bi))
end function div_matr

function matr_product(ai,prod) result (ri)
integer, dimension(:,:),intent(in)::ai
real ,intent(in):: prod
real, allocatable,dimension(:,:) :: ri
length = size(ai,1)
allocate(ri(length,length))
do i=1,length
do j=1,length
ri(j,i) = ai(j,i)*prod
enddo
enddo

end function matr_product

function matr_inversion(ai) result(ri)
integer, dimension(:,:),intent(in)::ai
real, allocatable,dimension(:,:) :: ri
integer, allocatable,dimension(:,:) ::M,MT
integer :: det
length = size(ai,1)
allocate(ri(length,length),M(length,length),MT(length,length))
det = matr_det(ai)
if (.NOT. det==0) then
do i=1,length
do j=1,length
M(j,i) =((-1)**(i+j))*matr_det(matr_minor(ai,j,i))
enddo
enddo
MT = matr_transpose(M)
ri = matr_product(MT,1./det)
deallocate(M,MT)
else
print *,"Ooops, there is no inverted matrix because det = 0"
endif
end function matr_inversion 

function mul_check_inversion(ai,bi) result (ri)
 integer :: length
 real, dimension(:,:),allocatable:: ri
 integer, dimension(:,:),intent(in)::ai
 real, dimension(:,:),intent(in)::bi
length = size(ai,1)
allocate(ri(length,length))
ri = 0
 do i=1,length
 do j=1,length
do k=1,length
  ri(j,i) = ri(j,i) + ai(k,i)*bi(j,k)
enddo
enddo
enddo
end function mul_check_inversion


end module matr_mod

program matr
use matr_mod
INTEGER,DIMENSION(:,:),  ALLOCATABLE :: a,b
 real,dimension(:,:), ALLOCATABLE :: ar,br
 integer :: lengths(8),n,detLengths(4),divLengths(9)
 integer, dimension (:,:),allocatable ::resi
real, dimension (:,:),allocatable ::resR
integer,allocatable :: seed(:)
real :: startTime,endTime,time
 lengths = (/2,4,10,100,200,400,800,1000/)
detLengths = (/2,4,10,11/)
divLengths=((/2,3,4,5,6,7,8,9,10/))
 do l=1,8
 allocate(a(lengths(l),lengths(l)),b(lengths(l),lengths(l)),resi(lengths(l),lengths(l)))
  call rand_int(a)
  call rand_int(b)
call cpu_time(startTime)
 resi= mul_int(a,b)
call cpu_time(endTime)
time= endTime-startTime
print '(i4,"x",i4" matrix multiplication takes ",f9.4," seconds")',lengths(l),lengths(l),time
deallocate(a,b,resi)

enddo

 do l=1,8
 allocate(a(lengths(l),lengths(l)),b(lengths(l),lengths(l)),resi(lengths(l),lengths(l)))
  call rand_int(a)
  call rand_int(b)
call cpu_time(startTime)
 resi= matr_add(a,b)
call cpu_time(endTime)
time= endTime-startTime
print '(i4,"x",i4" matrix summation takes ",f9.4," seconds")',lengths(l),lengths(l),time
deallocate(a,b,resi)

enddo

do l=1,8
 allocate(a(lengths(l),lengths(l)),b(lengths(l),lengths(l)),resi(lengths(l),lengths(l)))
  call rand_int(a)
  call rand_int(b)
call cpu_time(startTime)
 resi= matr_sub(a,b)
call cpu_time(endTime)
time= endTime-startTime
print '(i4,"x",i4" matrix substraction takes ",f9.4," seconds")',lengths(l),lengths(l),time
deallocate(a,b,resi)

enddo

do l=1,8
 allocate(a(lengths(l),lengths(l)),resi(lengths(l),lengths(l)))
  call rand_int(a)
call cpu_time(startTime)
 resi= matr_transpose(a)
call cpu_time(endTime)
time= endTime-startTime
print '(i4,"x",i4" matrix transposing takes ",f9.4," seconds")',lengths(l),lengths(l),time
deallocate(a,resi)

enddo

do l=1,4
 allocate(a(detLengths(l),detLengths(l)))
  call rand_int(a)
call cpu_time(startTime)
 det =matr_det(a)
call cpu_time(endTime)
time= endTime-startTime
print '(i4,"x",i4" matrix determinant takes ",f9.4," seconds")',detLengths(l),detLengths(l),time
deallocate(a)
enddo

do l=1,4
 allocate(a(detLengths(l),detLengths(l)),resi(lengths(l),lengths(l)))
  call rand_int(a)
call cpu_time(startTime)
 resi =matr_det(a)
call cpu_time(endTime)
time= endTime-startTime
print '(i4,"x",i4" matrix inversion takes ",f9.4," seconds")',detLengths(l),detLengths(l),time
deallocate(a,resi)
enddo

do l=1,9
 allocate(a(divLengths(l),divLengths(l)),b(divLengths(l),divLengths(l)),resR(divLengths(l),divLengths(l)))
  call rand_int(a)
call rand_int(b)
call cpu_time(startTime)
 resR =div_matr(a,b)
call cpu_time(endTime)
time= endTime-startTime
print '(i4,"x",i4" dividing matrix takes ",f9.4," seconds")',divLengths(l),divLengths(l),time
deallocate(a,b,resR)
enddo

end program matr
