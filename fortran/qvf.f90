! Fortran convenience writers over QVF ABI 1. License: Apache-2.0.
module qvf
  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use qvf_c
  implicit none
  private :: pack2, pack3, valid_text, status_message, bounded
  interface
    function qvf_strlen(text) bind(C, name='strlen') result(n)
      import
      type(c_ptr), value :: text
      integer(c_size_t) :: n
    end function
  end interface
  private :: qvf_strlen
contains
  ! This utility is also useful with raw c_qvf_* entry points. The returned
  ! array must have TARGET and remain alive through the C call.
  function qvf_c_string(text) result(bytes)
    character(len=*), intent(in) :: text
    character(c_char), allocatable :: bytes(:)
    integer :: i, n
    n = len_trim(text)
    allocate(bytes(n+1))
    do i=1,n
      bytes(i) = achar(iachar(text(i:i)), kind=c_char)
    end do
    bytes(n+1) = c_null_char
  end function

  function qvf_error_message(handle) result(message)
    type(c_ptr), intent(in), optional :: handle
    character(len=:), allocatable :: message
    type(c_ptr) :: ptr
    character(c_char), pointer :: bytes(:)
    integer :: i, n
    ptr = c_qvf_global_error()
    if (present(handle)) then
      if (c_associated(handle)) ptr = c_qvf_last_error(handle)
    end if
    message = ''
    if (.not. c_associated(ptr)) return
    n = int(qvf_strlen(ptr))
    call c_f_pointer(ptr, bytes, [n])
    do i=1,n
      message = message // achar(iachar(bytes(i)))
    end do
  end function

  logical function valid_text(text)
    character(len=*), intent(in) :: text
    valid_text = index(text, achar(0)) == 0
  end function

  logical function bounded(n)
    integer(c_int64_t), intent(in) :: n
    bounded = n > 0 .and. n <= int(huge(0_c_int), c_int64_t)
  end function

  ! Callers assign the result to their own optional message. gfortran 11 and
  ! 12 lose the length of a deferred-length optional dummy forwarded to a
  ! second procedure, so the message must not be passed through.
  function status_message(status, handle) result(message)
    integer(c_int), intent(in) :: status
    type(c_ptr), intent(in) :: handle
    character(len=:), allocatable :: message
    message = ''
    if (status /= QVF_OK) message = qvf_error_message(handle)
  end function

  function qvf_open(handle, path, program, version, calculation, message) result(status)
    type(c_ptr), intent(inout) :: handle
    character(len=*), intent(in) :: path, program, version, calculation
    character(len=:), allocatable, intent(out), optional :: message
    integer(c_int) :: status
    type(c_ptr), target :: created
    character(c_char), allocatable, target :: p(:), n(:), v(:), c(:)
    status = QVF_ERROR
    if (present(message)) message = 'open requires a null handle and nonempty path/program without NUL bytes'
    if (c_associated(handle)) return
    if (.not. valid_text(path) .or. .not. valid_text(program)) return
    if (.not. valid_text(version) .or. .not. valid_text(calculation)) return
    if (len_trim(path)==0 .or. len_trim(program)==0) return
    p=qvf_c_string(path); n=qvf_c_string(program); v=qvf_c_string(version); c=qvf_c_string(calculation)
    created=c_null_ptr
    status=c_qvf_writer_create(c_loc(n), c_loc(v), c_loc(c), c_loc(p), c_loc(created))
    if (status==QVF_OK) handle=created
    if (present(message)) message=status_message(status,handle)
  end function

  function qvf_finish(handle, message) result(status)
    type(c_ptr), intent(in) :: handle
    character(len=:), allocatable, intent(out), optional :: message
    integer(c_int) :: status
    status=c_qvf_writer_finish(handle)
    if (present(message)) message=status_message(status,handle)
  end function

  function qvf_close(handle) result(status)
    type(c_ptr), intent(inout) :: handle
    integer(c_int) :: status
    status=QVF_OK
    if (c_associated(handle)) status=c_qvf_writer_destroy(handle)
    if (status==QVF_OK) handle=c_null_ptr
  end function

  function pack2(values) result(flat)
    real(c_double), intent(in) :: values(:,:)
    real(c_double), allocatable :: flat(:)
    integer :: i,j,k
    allocate(flat(size(values))); k=0
    do i=1,size(values,1)
      do j=1,size(values,2)
        k=k+1; flat(k)=values(i,j)
      end do
    end do
  end function

  function pack3(values) result(flat)
    real(c_double), intent(in) :: values(:,:,:)
    real(c_double), allocatable :: flat(:)
    integer :: i,j,k,n
    allocate(flat(size(values))); n=0
    do i=1,size(values,1)
      do j=1,size(values,2)
        do k=1,size(values,3)
          n=n+1; flat(n)=values(i,j,k)
        end do
      end do
    end do
  end function

  function qvf_structure(handle, atomic_numbers, positions, message) result(status)
    type(c_ptr), intent(in) :: handle
    integer(c_int), intent(in) :: atomic_numbers(:)
    real(c_double), intent(in) :: positions(:,:) ! (atom, xyz), angstrom
    character(len=:), allocatable, intent(out), optional :: message
    integer(c_int) :: status
    integer(c_int), allocatable, target :: numbers(:)
    real(c_double), allocatable, target :: xyz(:)
    status=QVF_ERROR
    if (present(message)) message='structure needs finite positions(n_atoms,3) and atomic numbers 1..118'
    if (.not. bounded(size(positions,kind=c_int64_t))) return
    if (size(positions,1)/=size(atomic_numbers) .or. size(positions,2)/=3) return
    if (any(atomic_numbers<1) .or. any(atomic_numbers>118)) return
    if (.not. all(ieee_is_finite(positions))) return
    numbers=atomic_numbers; xyz=pack2(positions)
    status=c_qvf_add_structure(handle,c_null_ptr,c_loc(numbers),c_loc(xyz),int(size(numbers),c_int))
    if (present(message)) message=status_message(status,handle)
  end function

  function qvf_volume(handle, kind, grid_json, values, message) result(status)
    type(c_ptr), intent(in) :: handle
    character(len=*), intent(in) :: kind, grid_json
    real(c_double), intent(in) :: values(:,:,:) ! (x,y,z), Z fastest in archive
    character(len=:), allocatable, intent(out), optional :: message
    integer(c_int) :: status
    integer(c_int64_t), target :: dims(3)
    character(c_char), allocatable, target :: k(:), g(:), dtype(:)
    real(c_double), allocatable, target :: data(:)
    status=QVF_ERROR
    if (present(message)) message='volume needs nonempty finite scalar values and text without NUL bytes'
    if (.not. bounded(size(values,kind=c_int64_t))) return
    if (.not. all(ieee_is_finite(values))) return
    if (.not. valid_text(kind) .or. .not. valid_text(grid_json)) return
    dims=int(shape(values),c_int64_t); data=pack3(values)
    k=qvf_c_string(kind); g=qvf_c_string(grid_json); dtype=qvf_c_string('float64')
    status=c_qvf_add_volume(handle,c_loc(k),c_loc(g),c_loc(dtype),3_c_int,c_loc(dims),c_loc(data),c_null_ptr)
    if (present(message)) message=status_message(status,handle)
  end function

  function qvf_trajectory(handle, metadata_json, coords, message, lattice, velocities, gradients, times, energies) result(status)
    type(c_ptr), intent(in) :: handle
    character(len=*), intent(in) :: metadata_json
    real(c_double), intent(in) :: coords(:,:,:) ! (frame,atom,xyz), angstrom
    character(len=:), allocatable, intent(out), optional :: message
    real(c_double), intent(in), optional :: lattice(:,:,:), velocities(:,:,:), gradients(:,:,:), times(:), energies(:)
    integer(c_int) :: status, per_frame
    character(c_char), allocatable, target :: meta(:)
    real(c_double), allocatable, target :: xyz(:), cell(:), vel(:), grad(:), t(:), e(:)
    type(c_ptr) :: pc,pv,pg,pt,pe
    integer :: nf
    status=QVF_ERROR
    if (present(message)) message='trajectory arrays have incompatible shapes or non-finite values, or text contains NUL'
    if (.not. bounded(size(coords,kind=c_int64_t))) return
    if (size(coords,3)/=3 .or. .not. valid_text(metadata_json)) return
    if (.not. all(ieee_is_finite(coords))) return
    nf=size(coords,1); per_frame=0
    pc=c_null_ptr; pv=c_null_ptr; pg=c_null_ptr; pt=c_null_ptr; pe=c_null_ptr
    if (present(lattice)) then
      if (size(lattice,1)/=1 .and. size(lattice,1)/=nf) return
      if (size(lattice,2)/=3 .or. size(lattice,3)/=3) return
      if (.not. bounded(size(lattice,kind=c_int64_t))) return
      if (.not. all(ieee_is_finite(lattice))) return
      cell=pack3(lattice); pc=c_loc(cell)
      if (size(lattice,1)==nf .and. nf>1) per_frame=1
    end if
    if (present(velocities)) then
      if (any(shape(velocities)/=shape(coords))) return
      if (.not. all(ieee_is_finite(velocities))) return
      vel=pack3(velocities); pv=c_loc(vel)
    end if
    if (present(gradients)) then
      if (any(shape(gradients)/=shape(coords))) return
      if (.not. all(ieee_is_finite(gradients))) return
      grad=pack3(gradients); pg=c_loc(grad)
    end if
    if (present(times)) then
      if (size(times)/=nf) return
      if (.not. all(ieee_is_finite(times))) return
      t=times; pt=c_loc(t)
    end if
    if (present(energies)) then
      if (size(energies)/=nf) return
      if (.not. all(ieee_is_finite(energies))) return
      e=energies; pe=c_loc(e)
    end if
    meta=qvf_c_string(metadata_json); xyz=pack3(coords)
    status=c_qvf_add_trajectory(handle,c_loc(meta),int(nf,c_int),int(size(coords,2),c_int), &
      c_loc(xyz),pc,per_frame,pv,pg,pt,pe)
    if (present(message)) message=status_message(status,handle)
  end function

  function qvf_gto(handle, centers, angular_momenta, nprim, exponents, coefficients, mo, &
                   message, pure, libint_normalized, energies, occupations) result(status)
    type(c_ptr), intent(in) :: handle
    integer(c_int), intent(in) :: centers(:), angular_momenta(:), nprim(:) ! centers are ZERO-based
    real(c_double), intent(in) :: exponents(:), coefficients(:), mo(:,:) ! (MO, AO)
    character(len=:), allocatable, intent(out), optional :: message
    logical, intent(in), optional :: pure, libint_normalized
    real(c_double), intent(in), optional :: energies(:), occupations(:)
    integer(c_int) :: status, p, norm
    integer(c_int), allocatable, target :: c(:), l(:), n(:)
    real(c_double), allocatable, target :: exp(:), coeff(:), orbitals(:), e(:), o(:)
    type(c_ptr) :: pe,po
    integer :: ns,nm
    status=QVF_ERROR
    if (present(message)) message='GTO shell, primitive or MO arrays have invalid shapes or values'
    if (.not. bounded(size(mo,kind=c_int64_t)) .or. .not. bounded(size(centers,kind=c_int64_t))) return
    ns=size(centers); nm=size(mo,1)
    if (size(angular_momenta)/=ns .or. size(nprim)/=ns) return
    if (any(centers<0) .or. any(angular_momenta<0) .or. any(angular_momenta>6) .or. any(nprim<1)) return
    if (sum(int(nprim,c_int64_t))/=size(exponents,kind=c_int64_t)) return
    if (size(coefficients)/=size(exponents)) return
    ! Separate checks: .or. need not short-circuit, and an ordered NaN
    ! comparison signals IEEE invalid, which traps under halting mode.
    if (.not. all(ieee_is_finite(exponents))) return
    if (any(exponents<=0)) return
    if (.not. all(ieee_is_finite(coefficients)) .or. .not. all(ieee_is_finite(mo))) return
    p=1; norm=0; pe=c_null_ptr; po=c_null_ptr
    if (present(pure)) p=merge(1_c_int,0_c_int,pure)
    if (present(libint_normalized)) norm=merge(1_c_int,0_c_int,libint_normalized)
    if (present(energies)) then
      if (size(energies)/=nm) return
      if (.not. all(ieee_is_finite(energies))) return
      e=energies; pe=c_loc(e)
    end if
    if (present(occupations)) then
      if (size(occupations)/=nm) return
      if (.not. all(ieee_is_finite(occupations))) return
      o=occupations; po=c_loc(o)
    end if
    c=centers; l=angular_momenta; n=nprim; exp=exponents; coeff=coefficients; orbitals=pack2(mo)
    status=c_qvf_add_wavefunction_gto(handle,int(ns,c_int),c_loc(c),c_loc(l),c_loc(n),c_loc(exp),c_loc(coeff), &
                                    p,norm,int(nm,c_int),int(size(mo,2),c_int),c_loc(orbitals),pe,po)
    if (present(message)) message=status_message(status,handle)
  end function
end module qvf
