! Focused ownership, refusal and raw ABI checks. License: Apache-2.0.
program test_binding
  use qvf
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_invalid, &
    ieee_support_halting, ieee_get_halting_mode, ieee_set_halting_mode
  implicit none
  interface
    function c_size(index) bind(C,name='qvf_test_c_size') result(n)
      import
      integer(c_int), value :: index
      integer(c_size_t) :: n
    end function
  end interface
  type(c_ptr) :: writer=c_null_ptr, original
  type(c_ptr), target :: reader, report, bytes
  integer(c_size_t), target :: count, nbytes
  integer(c_int), target :: ok
  character(len=1024) :: path
  character(len=:), allocatable :: message
  character(c_char), allocatable, target :: text(:)
  type(qvf_limits_t), target :: limits
  type(qvf_tensor_view_t) :: tensor
  type(qvf_json_input_t) :: ji
  type(qvf_binary_input_t) :: bi
  type(qvf_opaque_input_t) :: oi
  type(qvf_section_input_t) :: si
  type(qvf_member_info_t), target :: mi
  type(qvf_finding_t) :: finding
  real(c_double) :: coords(2,2,3), mo(2,2), bad(2,2)
  integer(c_int) :: status
  logical :: halting
  real(c_double) :: nan
  call get_command_argument(1,path)
  if (len_trim(path)==0) error stop 'supply output path'
  call expect(c_size(0_c_int)==c_sizeof(limits),'limits ABI')
  call expect(c_size(1_c_int)==c_sizeof(tensor),'tensor ABI')
  call expect(c_size(2_c_int)==c_sizeof(ji),'json ABI')
  call expect(c_size(3_c_int)==c_sizeof(bi),'binary ABI')
  call expect(c_size(4_c_int)==c_sizeof(oi),'opaque ABI')
  call expect(c_size(5_c_int)==c_sizeof(si),'section ABI')
  call expect(c_size(6_c_int)==c_sizeof(mi),'member ABI')
  call expect(c_size(7_c_int)==c_sizeof(finding),'finding ABI')
  call expect(c_qvf_abi_version()==QVF_ABI_VERSION,'ABI version')
  call expect(c_qvf_limits_default(c_loc(limits))==QVF_OK,'default limits')
  call expect(limits%max_member_bytes>0,'limits output')
  status=qvf_open(writer,trim(path),'test','1','synthetic',message)
  call expect(status==QVF_OK,'open')
  original=writer
  status=qvf_open(writer,trim(path),'test','1','synthetic',message)
  call expect(status/=QVF_OK .and. c_associated(writer,original),'double open preserves handle')
  bad=0
  call expect(qvf_structure(writer,[1_c_int,1_c_int],bad)/=QVF_OK,'wrong xyz shape')
  coords=0; mo=0
  coords(1,1,1)=ieee_value(0.0_c_double,ieee_quiet_nan)
  call expect(qvf_structure(writer,[1_c_int,1_c_int],coords(1,:,:))/=QVF_OK,'non-finite positions')
  coords=0
  call expect(qvf_gto(writer,[0_c_int],[0_c_int],[2_c_int],[1.0_c_double], &
    [1.0_c_double],mo)/=QVF_OK,'primitive count bounds')
  call expect(qvf_gto(writer,[0_c_int],[0_c_int],[1_c_int],[1.0_c_double], &
    [1.0_c_double],mo,energies=[0.0_c_double])/=QVF_OK,'energy bounds')
  ! NaN inputs must be refused, not trapped, when the caller halts on IEEE invalid.
  nan=ieee_value(0.0_c_double,ieee_quiet_nan)
  if (ieee_support_halting(ieee_invalid)) then
    call ieee_get_halting_mode(ieee_invalid,halting)
    call ieee_set_halting_mode(ieee_invalid,.true.)
  end if
  call expect(qvf_gto(writer,[0_c_int],[0_c_int],[1_c_int],[nan], &
    [1.0_c_double],mo)/=QVF_OK,'NaN exponent under halting')
  call expect(qvf_gto(writer,[0_c_int],[0_c_int],[1_c_int],[1.0_c_double], &
    [nan],mo)/=QVF_OK,'NaN coefficient under halting')
  if (ieee_support_halting(ieee_invalid)) call ieee_set_halting_mode(ieee_invalid,halting)
  call expect(qvf_trajectory(writer,'{"symbols":["H","H"]}',coords,times=[0.0_c_double])/=QVF_OK,'time bounds')
  call expect(qvf_trajectory(writer,'{"symbols":["H","H"]}',coords,lattice=coords)/=QVF_OK,'cell bounds')
  call expect(qvf_volume(writer,'volume.generic','{}'//achar(0),coords)/=QVF_OK,'embedded NUL')
  call expect(qvf_volume(writer,'volume.generic','{}',coords,message)/=QVF_OK,'C refusal')
  call expect(len(message)>0,'C error copied')
  call expect(exact_message(message,qvf_error_message(writer)),'C error length')
  call expect(qvf_structure(writer,[1_c_int,1_c_int],coords(1,:,:))==QVF_OK,'valid strided structure')
  call expect(qvf_finish(writer)==QVF_OK,'finish')
  status=qvf_close(writer)
  call expect(status==QVF_OK .and. .not. c_associated(writer),'close clears ownership')
  call expect(qvf_close(writer)==QVF_OK,'repeat close')
  call expect(qvf_finish(writer)/=QVF_OK,'null handle')
  ! Each optional message must carry exactly the C error text and length.
  call expect(qvf_finish(writer,message)/=QVF_OK,'null handle message')
  call expect(exact_message(message,qvf_error_message()),'null handle message length')
  call expect(qvf_open(writer,trim(path)//'.missing/x.qvf','test','1','synthetic',message)/=QVF_OK, &
    'open refusal')
  call expect(exact_message(message,qvf_error_message()),'open refusal message length')
  text=qvf_c_string(trim(path)); reader=c_null_ptr
  call expect(c_qvf_reader_open(c_loc(text),c_loc(limits),c_loc(reader))==QVF_OK,'raw reader')
  call expect(c_qvf_reader_section_count(reader,c_loc(count))==QVF_OK,'raw count')
  call expect(count==1,'refusals wrote no sections')
  call expect(c_qvf_reader_manifest_json(reader,c_loc(bytes),c_loc(nbytes))==QVF_OK,'raw text output')
  call expect(c_associated(bytes) .and. nbytes>0,'borrowed text')
  report=c_null_ptr
  call expect(c_qvf_reader_validate(reader,QVF_VALIDATE_ALL,1_c_int,c_loc(report))==QVF_OK,'raw validate')
  call expect(c_qvf_report_ok(report,c_loc(ok))==QVF_OK,'report status')
  call expect(ok==1,'refusal archive valid')
  call expect(c_qvf_report_destroy(report)==QVF_OK,'report destroy')
  call expect(c_qvf_reader_destroy(reader)==QVF_OK,'reader destroy')
  print *, 'Fortran ownership, raw ABI and refusal checks passed'
contains
  logical function exact_message(actual,expected)
    character(len=*), intent(in) :: actual, expected
    exact_message=len(expected)>0 .and. len(actual)==len(expected) .and. actual==expected
  end function

  subroutine expect(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print *,label
      error stop 'Fortran binding check failed'
    end if
  end subroutine
end program
