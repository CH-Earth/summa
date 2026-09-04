module mizuroute_config

  USE nr_type
  USE summa_type, only:summa1_type_dec

  implicit none
  private

  public :: parse_mizuroute_config

contains

  subroutine parse_mizuroute_config(subtable, section, key, summaStruct, ierr, message)
  
  use tomlf_all, only: toml_table, toml_array, toml_error, toml_key, toml_value ! data types
  use tomlf_all, only: toml_load, get_value, len                                ! procedures

  type(toml_table), pointer, intent(in)    :: subtable
  character(*),              intent(in)    :: section
  character(*),              intent(in)    :: key
  type(summa1_type_dec),     intent(inout) :: summaStruct
  integer,                   intent(out)   :: ierr
  character(*),              intent(out)   :: message

  integer(i4b)       :: istat

  associate(info => summaStruct%mizu_info)

  ierr    = 0
  message = 'parse_mizuroute_config/'

  ! extract configuration values and populate the mizuRoute information structure

  select case(trim(section)//'.'//trim(key))
  
    ! ---- mizuRoute: namelist path/filenames ----
    case ("mizuRoute.namelist_path"      ); call get_value(subtable, trim(key), info%mrout%namelist_path    , stat=istat)
    case ("mizuRoute.namelist_file"      ); call get_value(subtable, trim(key), info%mrout%namelist_file    , stat=istat)
    
    ! ---- mizuRoute: runtime ----
    case ("mizuRoute.dt"                 ); call get_value(subtable, trim(key), info%mrout%dt               , stat=istat)
    case ("mizuRoute.methods"            ); call get_value(subtable, trim(key), info%mrout%methods          , stat=istat)
    
    ! ---- hydrofabric: path/filenames ----
    case ("hydrofabric.hfabric_path"     ); call get_value(subtable, trim(key), info%ntopo%hfabric_path     , stat=istat)
    case ("hydrofabric.hfabric_file"     ); call get_value(subtable, trim(key), info%ntopo%hfabric_file     , stat=istat)
    case ("hydrofabric.hfabric_newfile"  ); call get_value(subtable, trim(key), info%ntopo%hfabric_newfile  , stat=istat)
    
    ! ---- hydrofabric: dimensions ----
    case ("hydrofabric.dname_hru"        ); call get_value(subtable, trim(key), info%ntopo%dname_hru        , stat=istat)
    case ("hydrofabric.dname_seg"        ); call get_value(subtable, trim(key), info%ntopo%dname_seg        , stat=istat)
    
    ! ---- hydrofabric: variable names ----
    case ("hydrofabric.varname_HRUid"    ); call get_value(subtable, trim(key), info%ntopo%varname_HRUid    , stat=istat)
    case ("hydrofabric.varname_segId"    ); call get_value(subtable, trim(key), info%ntopo%varname_segId    , stat=istat)
    case ("hydrofabric.varname_hruSegId" ); call get_value(subtable, trim(key), info%ntopo%varname_hruSegId , stat=istat)
    case ("hydrofabric.varname_downSegId"); call get_value(subtable, trim(key), info%ntopo%varname_downSegId, stat=istat)
    case ("hydrofabric.varname_area"     ); call get_value(subtable, trim(key), info%ntopo%varname_area     , stat=istat)
    case ("hydrofabric.varname_slope"    ); call get_value(subtable, trim(key), info%ntopo%varname_slope    , stat=istat)
    case ("hydrofabric.varname_length"   ); call get_value(subtable, trim(key), info%ntopo%varname_length   , stat=istat)
    
    ! ---- hydrofabric: network topology ----
    case ("hydrofabric.seg_outlet"       ); call get_value(subtable, trim(key), info%ntopo%idSegOut         , stat=istat)
    
    ! ---- remapping: filename ----
    case ("remapping.remap_path"         ); call get_value(subtable, trim(key), info%remap%remap_path       , stat=istat)
    case ("remapping.remap_file"         ); call get_value(subtable, trim(key), info%remap%remap_file       , stat=istat)
    
    ! ---- remapping: dimension names ----
    case ("remapping.dname_hru"          ); call get_value(subtable, trim(key), info%remap%dname_hru        , stat=istat)
    case ("remapping.dname_data"         ); call get_value(subtable, trim(key), info%remap%dname_data       , stat=istat)
    
    ! ---- remapping: variable names ----
    case ("remapping.vname_hruid"        ); call get_value(subtable, trim(key), info%remap%vname_hruid      , stat=istat)
    case ("remapping.vname_weight"       ); call get_value(subtable, trim(key), info%remap%vname_weight     , stat=istat)
    case ("remapping.vname_num_qhru"     ); call get_value(subtable, trim(key), info%remap%vname_num_qhru   , stat=istat)
    case ("remapping.vname_i_index"      ); call get_value(subtable, trim(key), info%remap%vname_i_index    , stat=istat)
    case ("remapping.vname_j_index"      ); call get_value(subtable, trim(key), info%remap%vname_j_index    , stat=istat)
    case ("remapping.vname_qhruid"       ); call get_value(subtable, trim(key), info%remap%vname_qhruid     , stat=istat)
   
    ! ---- default case (something in the table that is not specified above) -----
    case default
      message = trim(message)// "unexpected entry: section = "//trim(section)//"; sub-section = "//trim(key)
      ierr=20; return
   
  end select ! (select key/value pair based on lookup)
     
  ! ---- error checking -----
  if(istat /= 0)then
    message=trim(message)// "get_value error: section = "//trim(section)//"; sub-section = "//trim(key)
    ierr=20; return
  endif
 
  end associate

  end subroutine parse_mizuroute_config

end module mizuroute_config
