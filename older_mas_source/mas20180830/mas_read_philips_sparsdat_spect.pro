;+
;
; This file contains a very simple lookup-table implementation, mainly
; to be used to hold the paramaters read from various data files.
; 
; :Author: wtriplett
;-



function lookuptable::init

    self.num_elements = 0L
    return, 1
    
end

pro lookuptable::cleanup

    ptr_free, self.keys__
    ptr_free, self.data__

end

pro lookuptable::set, key, value

    if (size(key, /type) ne 7) then begin
        message, 'Key must be a string type.'
    endif
        
    strkey = string(self.num_elements, format='(%"STRKEY%05d")')
    
    if (self->keyExists(key, match_id=id)) then begin
        
        (*self.data__).(id) = value
        return
        
    endif else if (not ptr_valid(self.data__)) then begin
    
        self.num_elements = 0L
        self.data__ = ptr_new(create_struct(strkey, value))
        self.keys__ = ptr_new([key])
    
    endif else begin

        *self.data__ = create_struct(*self.data__, strkey, value)
        *self.keys__ = [ *self.keys__, key ]
    
    endelse
    
    self.num_elements++

end

function lookuptable::count

    return, self.num_elements

end

function lookuptable::isEmpty

    if (self.num_elements eq 0) then return, 1
    
    return, 0

end

function lookuptable::keyExists, key, match_id=match_id

    if (not ptr_valid(self.keys__)) then return, 0
    
    matches  = strcmp(key, *self.keys__)
    match_id = where(matches ne 0, n_matches)
    
    return, (n_matches eq 0) ? 0 : 1

end

function lookuptable::getValueForKey, key, found=found

    if (size(key, /type) ne 7) then begin
        message, 'Key must be a string type.'
    endif

    if (not self->keyExists(key, match_id=match_id)) then begin
        found = 0
        message, /info, 'Key ('+key+') not found'
        return, ''
    endif

    found = 1
    return, (*self.data__).(match_id[0])

end

pro lookuptable__define


    struct = { LOOKUPTABLE, $
               data__: ptr_new(), $ ;; holds the values
               keys__: ptr_new(), $ ;; holds the keys or names
               num_elements: 0 $
             }

end


;+
; :Description:
;    Reads a Philips "SDAT" file into a FID array.
;
; :Params:
;    infile - path to the sdat file
; 
; :Keywords:
;    error - this will be 1 if there was a problem opening the file.
;
; :Author: wtriplett
;-
function mas_read_philips_sdat, infile, error=error 

    if (n_elements(infile) eq 0) then begin
        infile = dialog_pickfile(filter='*.SDAT',get_path=working_dir,/must_exist,title='Please Select Philips .SDAT file')
    endif else begin
        working_dir = file_dirname(infile)
    endelse
   
   openr, lun, infile, /get_lun, error=error
   if (error ne 0) then begin
       message, /info, 'Unable to read file ('+infile+').'
       return, -1
   endif
   
   fid = read_binary(lun,data_type = 6,data_start = 8*0) 
   
   byteorder, fid, /lswap 
   byteorder, fid 
   
   close, lun & free_lun, lun 
   
   return, fid 
  
end 

;+
; :Description:
;    Read a Philips "SPAR" file into a lookuptable (see lookuptable__define.pro).
;
; :Params:
;    infile - path to spar file
;
; :Keywords:
;    error - this will be set to 1 if there was a problem opening the file.
;
; :Author: wtriplett
;-
function mas_read_philips_spar, infile, error=error

  openr, lun, infile, /get_lun, error=error
  if (error ne 0) then begin
    message, /info, 'Unable to read file ('+infile+').'
    return, -1
  endif
  
  line = ''
  
  spar = obj_new('lookuptable')
  
  while (not eof(lun)) do begin
  
    readf, lun, line
    
    data = strsplit(line, ':', /extract)
    
    ;; scan_date has : in it already (time) so it causes a bit of trouble.
    if (strtrim(data[0],2) eq 'scan_date') then begin
        temp = strjoin(data[1:*], ':')
        data = [data[0], temp]
    endif
    
    if (n_elements(data) eq 2) then begin
    
      key = strtrim(data[0],2)
      
      if (stregex(key, '^[0-9a-zA-Z_]+$', /boolean) eq 0) then begin
        message, /info, 'Throwing away key: '+key+'.', noprint=1
        continue
      endif
      
      temp = strtrim(data[1], 2)
      
      is_number = stregex(temp, '^[0-9\.\e\+\-]+$', /boolean)
      
      ;; this is way more complicated than it needs to be.
      if (key eq 'patient_name') then begin
        is_number = 0
      endif
      
      if (is_number ne 0) then begin
      
        spar->set, key, float(temp)
        
      endif else begin
        
        spar->set, key, temp
        
      endelse
      
    endif
    
  endwhile

  close, lun & free_lun, lun
  
  return, spar
  
end

;+
; :Description:
;    Read a philips SPAR/SDAT combo..
;
; :Params:
;    file  - path to either .SPAR or .SDAT file
;
; :Keywords:
;    params - on output, this keyword will have the spectroscopy information.
;             if there are most than one FID in the file, then params will be
;             an array of structres of type SPECTROSCOPY_DATA (see end of this file)
;    fid    - on output this keyword will contain the FID data (num_points x num_fids)
;    plot   - set this keyword to have a plot shown, basically just calls plot_spect
;    valid_file - this keyword will be set to 0 if there was a problem opening the data,
;                 1 otherwise
;
; :Author: wtriplett
;-
pro mas_read_philips_sparsdat_spect


    common scan_data, project
    
    valid_file = 0

    next_index = project.ni
    imnd = project.imndarray[next_index]
    
    if (n_elements(file) eq 0) then file = dialog_pickfile(/read, /must_exist)
    if (file eq '') then return
    
    input_ext = strmid(file, strlen(file)-4)
    
    ;; Try to deduce the name of the other tilfe
    if (strupcase(input_ext) eq 'SPAR') then begin
        spar_file = file
        sdat_file = strmid(file, 0, strlen(file)-4)+'SDAT'
    endif else if (strupcase(input_ext) eq 'SDAT') then begin
        spar_file = strmid(file, 0, strlen(file)-4)+'SPAR'
        sdat_file = file
    endif 

    spar = mas_read_philips_spar(spar_file, error=spar_error)   
    if (spar_error ne 0) then begin
        valid_file = 0
        message, /info, 'Unable to read .SPAR file.'
        return
    endif


    imnd.spect_nucleus        = spar->getValueForKey('nucleus')
    imnd.spect_spectral_width = float(spar->getValueForKey('sample_frequency'))
    imnd.spect_acq_size       = spar->getValueForKey('samples')
    imnd.spect_num_avg        = spar->getValueForKey('averages')
    imnd.spect_bf1            = spar->getValueForKey('synthesizer_frequency') * 1e-6
    at = spar->getValueForKey('dim1_step') * spar->getValueForKey('dim1_pnts')
    imnd.spect_acq_time       = float(at)*100.0

    imnd.fdim = spar->getValueForKey('samples')
    imnd.adim = spar->getValueForKey('rows')
    imnd.n_avg       = imnd.spect_num_avg
    imnd.scan_name   = spar->getValueForKey('scan_id')
    imnd.orientation = [ "----", "----","----" ]
    imnd.spect_rep_time = spar->getValueForKey('repetition_time')
    imnd.dimensions = 1
    imnd.file_Path = file
    imnd.display_name = file
    imnd.image_type = 99
    imnd.acq_matrix = ptr_new(diag_matrix(replicate(1.0D,4)))
    
;    params = create_struct(name='spectroscopy_data')
;    params.manufacturer = 'PHILIPS'
;    params.src_file     = file 
;    params.src_type     = 'SPAR/SDAT' 
;    params.subject_id   = spar->getValueForKey('patient_name') 
;    params.site_name    = 'University of Florida'
;    params.num_points   = spar->getValueForKey('samples')
;    params.img_freq     = spar->getValueForKey('synthesizer_frequency') * 1e-6
;    params.protocol     = spar->getValueForKey('scan_id')
;    params.num_avgs     = spar->getValueForKey('averages')
;    params.pix_bw       = spar->getValueForKey('sample_frequency')
;    params.tr           = spar->getValueForKey('repetition_time')
;    params.te           = spar->getValueForKey('echo_time')
;    params.vx_ofc_AP    = spar->getValueForKey('ap_off_center')
;    params.vx_ofc_FH    = spar->getValueForKey('cc_off_center')
;    params.vx_ofc_RL    = spar->getValueForKey('lr_off_center')
;    params.vx_sz_AP     = spar->getValueForKey('ap_size')
;    params.vx_sz_FH     = spar->getValueForKey('cc_size')
;    params.vx_sz_RL     = spar->getValueForKey('lr_size')
;    params.num_xpts     = spar->getValueForKey('dim2_pnts') ;; For multi voxel
;    params.num_ypts     = spar->getValueForKey('dim3_pnts') ;; For multi voxel
    
;    ;; this format for scan date/time is not final.
    tmp_scandate = strtrim(spar->getValueForKey('scan_date'),2)
    tmp_datetime = strsplit(tmp_scandate, ' ', /extract)
    scan_date = strjoin(strsplit(tmp_datetime[0], '.', /extract), '-')
    scan_time = tmp_datetime[1] ;;strjoin(strsplit(tmp_datetime[1], ':', /extract), ''))
    imnd.scan_date = scan_date+' '+scan_time
    
    obj_destroy, spar
    
    fid = mas_read_philips_sdat(sdat_file, error=sdat_error)
    if (sdat_error ne 0) then begin
        valid_file = 0
        message, /info, 'Unable to read .SPAR file.'
        return
    endif
    
    num_fids = n_elements(fid)/imnd.fdim
    
    fid = reform(fid, imnd.fdim, num_fids)

    project.imndarray[next_index] = imnd

    project.scan_Open_Flag = 1
    
    project.ci = project.ni
    project.ni = project.ni + 1
    
    mas_add_current_scan_to_scan_selection

    project.dataarray[project.ci].state1 = ptr_new(fid, /no_copy)
    project.procpramarray[project.ci].state_1 = 1

end