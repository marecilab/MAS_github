FUNCTION calc_experimental_chi, slice
  
  ;setting up the default values
  B0 = 4.7
  mu0 = 4*!DPI*1e-07;
  

  ; setting up the matrix points
  xdimension = .03
  ydimension = .03

  desired_x_resolution = .03/128.0
  desired_y_resolution = .03/128.0

  ;setting up the angle of main magnetic field vs the model
  zeta = 90.0*(!PI/180)
  cos_zeta = cos( zeta )
  sin_zeta = sin( zeta )

  samples_x = xdimension/desired_x_resolution
  samples_y = ydimension/desired_y_resolution

  sizearray = [samples_x, samples_y]

  xdim_center = float(round(samples_x/2))
  ydim_center = float(round(samples_y/2))
  
  delta_chi_2 = complexarr(samples_x, samples_y)
  Mk = fltarr(samples_x, samples_y)


  ; this is a flag to only calculate the bout once
  ;closed_solution = 0

  ; Assuming susceptibility << 1
  Mz = (slice) * (B0/mu0);
  Mk_temp = (FFT(Mz))

  Mk = shift(Mk_temp, xdim_center,ydim_center)

  ; Creating a progress bar outside of the looping through the core material sizes
  progressbar = Obj_New('progressbar', Color='red', Text='Calculating Field', /NOCANCEL)
  progressbar -> Start

  for j=0.0, samples_y-1 do begin
      for i=0.0, samples_x-1 do begin
          k_x = (float(i)-xdim_center)
        k_y = (float(j)-ydim_center)
        k_z = float(0)

        if k_x eq 0.0 and k_y eq 0.0 and k_z eq 0.0 then begin
          delta_chi_2[i,j] = 0.0
        endif else begin
          cosinepart = ((k_z*cos_zeta - k_y*sin_zeta)^2.0)/(k_z^2.0 + k_x^2.0 + k_y^2.0)
          delta_chi_2[i,j] = -(Mu0/3.0)*Mk[i,j]*(3.0*cosinepart - 1.0)
        endelse

      endfor
    progressbar -> Update, (float(j)/float(samples_y))*100
    endfor

  progressbar -> Destroy

    delta_chi_2 = shift(delta_chi_2, xdim_center,ydim_center)
    results = REAL_PART(FFT(delta_chi_2,/INVERSE))
    

  return, results

END
; 
;helper routines to:
; define structures
; find CI, TI
; change global parameters

;writes a vector to a .dat file that can be read into MATLAB for plotting
PRO write_dat_file, vector, title
  fpath = DIALOG_PICKFILE( Title = title)
  OpenW, lun,fpath + '.dat', /get_lun
  
  sz_vector = size(vector)
  length = sz_vector[1]
  if sz_vector[2] gt sz_vector[1] then begin
    vector = transpose(vector)
    length = sz_vector[2]
  endif
  
  for i=0, length - 1 do begin
    temp = vector[i]
    PrintF, lun, format = '(E0)', temp
  endfor
  
  close, lun
  free_lun, lun
END
; Subroutine name: mas_extract_delays
; Created by: Garrett Astary
; Calling Information: Called by Calculate Delta B Algebra to find free precession delay

; Bugs or Important Comments to Developers:


; Purpose of subroutine: Read ACQP file and return the vector of delays stored therein

; Editing Information:
function create_chi_map, data, threshold

  size_data = size(data)
  max_data = max(data)
  geo_map = dblarr(size_data(1), size_data(2)) - 9.05E-6 
  for i = round(size_data(1)/4.0) - 1, round(size_data(1)/4.0) + round(size_data(1)/2.0) -1 do begin
    for j = round(size_data(2)/4.0) - 1, round(size_data(2)/4.0) + round(size_data(2)/2.0) - 1 do begin
      if data[i,j] lt threshold*max_data then begin
        geo_map[i,j] = -12E-6
      endif
    endfor
   endfor
     
return, geo_map
end

function mas_extract_delays, path, load_procedure, te = te
  
  load_procedure_split = strsplit(load_procedure, '_', /extract)
  varian_test = strmatch(load_procedure_split, 'varian')
  
  if (total(varian_test) EQ 1) then begin
    oProcpar = obj_new('mas_varian_procpar', path)
    
    delay = oProcpar->lookup('fpdelay', is_array=ia, not_found=nf)
     
  obj_destroy, oProcpar
  endif else begin
  
    path = path + 'acqp'
   
    
    ;reading delays from X data set
    OPENR, text_file_pointer, path , /GET_LUN , ERROR = err
    number_lines = FILE_LINES(path)
    TEXT_file_array = strarr(number_lines)
    counter = 0
    
    parameter_temp_string = ''
    
    While (eof(text_file_pointer) NE 1) Do Begin
           Readf, text_file_pointer, parameter_temp_string
           TEXT_file_array[counter] = parameter_temp_string
           counter = counter + 1
    End
    
    array1 = strpos(text_file_array, '##$D=')
    array2 = strpos(text_file_array,'##$P=')
    array3 = strpos(text_file_array,'##$PLW=')
    
    index1 = where(array1 eq 0)
    index2 = where(array2 eq 0)
    index3 = where(array3 eq 0)
       
    delay_temp = text_file_array[index1[0] + 1]
    pulse_temp = text_file_array[index2[0] + 1]
    
    for i=index1[0]+2, index2[0]-1 do delay_temp = [delay_temp, text_file_array[i]]
    for i=index2[0]+2, index3[0]-1 do pulse_temp = [pulse_temp, text_file_array[i]]
    
    delay_temp = strjoin(delay_temp)
    delay_temp = strsplit(delay_temp, /extract)  
    d = float(delay_temp)
    
    pulse_temp = strjoin(pulse_temp)
    pulse_temp = strsplit(pulse_temp, /extract)
    p = float(pulse_temp)

    delay = ((1e-6*(p[0]+p[1])/2.0 + d[1]+d[2]+d[3])*1E3  - te/2)*1E-3
      free_lun,text_file_pointer
  endelse  

  return, delay
end

; Subroutine name: cdb_create_mask
; Created by: Garrett Astary
; Calling Information: Called by mas_calc_delta_B_display and mas_calc_delta_B

; Bugs or Important Comments to Developers:


; Purpose of subroutine: Creates mask (1's and 0's) to filter out noise areas before displaying results of complex division or delta B calculation.
;The masks are generated using a threshold and the magnitude mode of the images.  

; Editing Information:

pro cdb_create_mask, info, output1 = output1, output2 = output2
  COMPILE_OPT IDL2
  COMMON scan_data
  COMMON common_widgets
  CI = project.CI
  ;create a mask for the entire data set
  dim  = project.imndArray[info.scan1].dimensions
  
  if info.display_3D eq 0 then begin
    slice = project.procpramarray[info.scan1].sdim_start
    ;convert data from fft complex mode to magnitude mode
    p_data_x = abs((*project.dataArray[info.scan1].state1)[*,*,slice])
    p_data_y = abs((*project.dataArray[info.scan2].state1)[*,*,slice])
    size_x_y = size(p_data_x)
    mask_x = dblarr(size_x_y[1], size_x_y[2]) + 1.0
    mask_y = dblarr(size_x_y[1], size_x_y[2]) + 1.0
    max_x = max(p_data_x)
    max_y = max(p_data_y)
    threshold_x = where(p_data_x lt info.threshold*max_x)
    threshold_y = where(p_data_y lt info.threshold*max_y)
    mask_x[threshold_x] = 0
    mask_y[threshold_y] = 0
  endif else begin
    ;convert data from fft complex mode to magnitude mode
    p_data_x = abs(*project.dataArray[info.scan1].state1)
    p_data_y = abs(*project.dataArray[info.scan2].state1)
    size_x_y = size(p_data_x)
    mask_x = dblarr(size_x_y[1], size_x_y[2], size_x_y[3]) + 1.0
    mask_y = dblarr(size_x_y[1], size_x_y[2], size_x_y[3]) + 1.0
    max_x = max(p_data_x)
    max_y = max(p_data_y)
    threshold_x = where(p_data_x lt double(info.threshold*max_x))
    threshold_y = where(p_data_y lt double(info.threshold*max_y))
    mask_x[threshold_x] = 0
    mask_y[threshold_y] = 0
  endelse 

output1 = mask_x
output2 = mask_y
info.mask_created=1

end

; Subroutine name: mas_calc_delta_b_display
; Created by: Garrett Astary
; Calling Information: mas_calc_delta_b_gui_event

; Bugs or Important Comments to Developers:


; Purpose of subroutine: After complex division has been completed (and the result stored in the project structure under algebra_result or algebra_result_3D)
; this routine is used to display the magnitude/real/imaginary/phase images of the result.   

; Editing Information:
pro mas_calc_delta_b_display, info
  COMPILE_OPT IDL2
  COMMON scan_data
  COMMON common_widgets
  CI = project.CI
  
  
    image_type = info.image_type
    
    ;2D display (2D data or single slice of 3D data)
    if info.display_3D eq 0 then begin
      if image_type eq 0 then begin ;magnitude
        image = abs(*project.algebra_result)
      endif else if image_type eq 1 then begin ;real
        image= REAL_PART((*project.algebra_result))
        image = abs(image)
      endif else if image_type eq 2 then begin ;imaginary
        image = imaginary((*project.algebra_result))
        image = abs(image)
      endif else if image_type eq 3 then begin ;phase
        image = atan((*project.algebra_result), /phase)
        if info.mask eq 0 then mag_image = make_Array(size(image,/dimensions),value = 1.0, /float) else $
          cdb_create_mask, info, output1 = mag_image
        if info.phase_unwrap eq 1 then image = mas_phase_unwrap(2, mag_image*exp(dcomplex(0,1)*image))
        if info.hipass eq 1 then begin
          image = mas_phase_unwrap(2, mag_image*exp(dcomplex(0,1)*image))
          kernelSize = [3, 3]
          kernel = REPLICATE(-1., kernelSize[0], kernelSize[1])
          kernel[1, 1] = 8
          image = CONVOL(FLOAT(image), kernel, 9, /CENTER, /EDGE_TRUNCATE)
        endif
        endif
     endif else begin
      if image_type eq 0 then begin ;magnitude
        image = abs(*project.algebra_result_3D)
      endif else if image_type eq 1 then begin ;real
        image= REAL_PART((*project.algebra_result_3D))
        image = abs(image)
      endif else if image_type eq 2 then begin ;imaginary
        image = imaginary((*project.algebra_result_3D))
        image = abs(image)
      endif else if image_type eq 3 then begin ;phase
        
        image = atan((*project.algebra_result_3D), /phase)
        if info.mask eq 0 then mag_image = make_Array(size(image,/dimensions),value = 1.0, /float) else $
          cdb_create_mask, info, output1 = mag_image
        if info.phase_unwrap eq 1 then image = mas_phase_unwrap(3, mag_image*exp(dcomplex(0,1)*image))
        if info.hipass eq 1 then begin
          image = mas_phase_unwrap(3, mag_image*exp(dcomplex(0,1)*image))
          kernelSize = [3, 3]
          kernel = REPLICATE(-1., kernelSize[0], kernelSize[1])
          kernel[1, 1] = 8
          image = CONVOL(FLOAT(image), kernel, 9, /CENTER, /EDGE_TRUNCATE)
        endif
      endif
      endelse 
 
  if info.display_3D eq 0 then begin  
    if info.mask eq 0 then begin
      iimage, image
    endif else begin
      cdb_create_mask, info, output1 = mask_x, output2 = mask_y
      image*= mask_x*mask_y
      iimage, image
    endelse
    
  endif else begin
    if info.mask eq 0 then begin
      my_data_ptr = ptr_new(image)
      mas_display_ortho, data_ptr = my_data_ptr
    endif else begin
      cdb_create_mask, info, output1 = mask_x, output2 = mask_y
      image*= mask_x*mask_y
      my_data_ptr = ptr_new(image)
      ;Orthogonal slice viewer created by Bill Triplett
      mas_display_ortho, data_ptr = my_data_ptr
    endelse
  endelse   

end

; Subroutine name: mas_calc_delta_B
; Created by: Garrett Astary
; Calling Information: mas_calc_delta_b_gui_event

; Bugs or Important Comments to Developers:


; Purpose of subroutine: Calcualtes change in B field for 2D or 3D images. The phase image is extracted from the result of the 
; complex division and delta B is calculated by dividing the phase image by  (delta TE * gyromagnetic ratio) 

; Editing Information:
pro mas_calc_delta_b, info
  COMPILE_OPT IDL2
  COMMON scan_data
  COMMON common_widgets
  CI = project.CI
  
  dTE = info.dTE/1000 ;need to convert to seconds
  gyromagnetic_ratio = (info.gyromagnetic_ratio)*1E6 ;convert to units of Hz/Tesla
  B0 = 1E6*project.imndarray[ci].spect_bf1/gyromagnetic_ratio
  if info.display_3D eq 0 then begin
    phase_image = atan((*project.algebra_result), /phase)
    if info.mask eq 0 then mag_image = make_Array(size(phase_image,/dimensions),value = 1.0, /float) else $
      cdb_create_mask, info, output1 = mag_image
    if info.phase_unwrap eq 1 then phase_image = mas_phase_unwrap(2, mag_image*exp(dcomplex(0,1)*phase_image))
    if info.hipass eq 1 then begin
;      phase_image = mas_phase_unwrap(2, mag_image*exp(dcomplex(0,1)*phase_image))
      kernelSize = [3, 3]
      kernel = REPLICATE(-1., kernelSize[0], kernelSize[1])
      kernel[1, 1] = 8
      phase_image = CONVOL(FLOAT(phase_image), kernel, 9, /CENTER, /EDGE_TRUNCATE)  
    endif
    
    ;delta_B is in units of ppm
    
    if info.free_precession_flag eq 0 then begin      
      delta_B = -1e6*phase_image/(dTE*gyromagnetic_ratio*2.0*!DPI*B0)
    endif else begin
      delta_B = 1e6*phase_image/(2.0*dTE*gyromagnetic_ratio*2.0*!DPI*B0)
    endelse
    if info.mask eq 0 then begin
      iimage, delta_B, title = 'delta B (ppm)'
    endif else begin
      cdb_create_mask, info, output1 = mask_x, output2 = mask_y
      delta_B = delta_B*mask_x*mask_y
      iimage, delta_B, title = 'delta B (ppm)'
    endelse
    
  endif else begin
    phase_image = atan(*project.algebra_result_3D, /phase) 
    if info.mask eq 0 then mag_image = make_Array(size(phase_image,/dimensions),value = 1.0, /float) else $
      cdb_create_mask, info, output1 = mag_image
    if info.phase_unwrap eq 1 then phase_image = mas_phase_unwrap(project.imndArray[ci].dimensions,mag_image*exp(dcomplex(0,1)*phase_image))
    if info.hipass eq 1 then begin
;      phase_image = mas_phase_unwrap(project.imndArray[ci].dimensions, mag_image*exp(dcomplex(0,1)*phase_image))
      kernelSize = [3, 3, 3]
      kernel = REPLICATE(-1., kernelSize[0], kernelSize[1], kernelSize[2])
      kernel[1, 1, 1] = 26
      phase_image = CONVOL(FLOAT(phase_image), kernel, 27, /CENTER, /EDGE_TRUNCATE)
    endif
    ;delta_B is in units of ppm
    delta_B = 1e6*phase_image/(dTE*gyromagnetic_ratio*2.0*!DPI*B0)
    if info.mask eq 0 then begin
       my_data_ptr = ptr_new(delta_B)
       ;Orthogonal slice viewer created by Bill Triplett
       mas_display_ortho, data_ptr = my_data_ptr
    endif else begin
      cdb_create_mask, info, output1 = mask_x, output2 = mask_y
      delta_B*= mask_x*mask_y
      my_data_ptr = ptr_new(delta_B)
      ;Orthogonal slice viewer created by Bill Triplett
      mas_display_ortho, data_ptr = my_data_ptr
    endelse
  endelse
end

; Subroutine name: mas_calc_delta_b_gui_event
; Created by: Garrett Astary
; Calling Information: Called when events triggered by main gui (mas_calc_delta_b_gui)

; Bugs or Important Comments to Developers:


; Purpose of subroutine: Event handler

; Editing Information:

pro mas_calc_delta_b_algebra_gui_event, event

 ;bring in global variables
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets
    
    ;grab parameter structure that is stored in the uvalue of the mas_multi_movie_gui
    widget_control, event.top, get_uvalue=info
    CI = project.CI
    sdim_start = project.procpramArray[CI].sdim_start
    adim_start = project.procpramArray[CI].adim_start 
    fdim_start = project.procpramArray[CI].fdim_start
    pdim_start = project.procpramArray[CI].pdim_start
    slice_axis = project.procpramArray[CI].slice_axis
    CASE event.id OF
    
     
      cdb_ge_button: BEGIN
      info.free_precession_flag = 0
      if ptr_valid(project.IMNDArray[info.scan1].echo_Time_ptr) then te1 = *project.IMNDArray[info.scan1].echo_Time_ptr[0] else $
        te1 = project.IMNDArray[info.scan1].echo_Time
        
      if ptr_valid(project.IMNDArray[info.scan2].echo_Time_ptr) then te2 = *project.IMNDArray[info.scan2].echo_Time_ptr[0] else $
        te2 = project.IMNDArray[info.scan2].echo_Time
      
        if info.scan1 gt -1 and info.scan2 gt -1 then begin
          adim1 = project.procPramArray[info.scan1].adim_start
          adim2 = project.procPramArray[info.scan2].adim_start
          if n_elements(te1) ge adim1 then info.TEx = te1[adim1] else info.TEx = te1
          if n_elements(te2) ge adim2 then info.TEy = te2[adim2] else info.TEy = te2          
          info.dTE = info.TEx - info.TEy
          TE_text = strcompress(string(info.dTE),/remove_all)
          widget_control, cdb_cdb_TE_text, set_value = TE_text
          widget_control, cdb_cd_button, sensitive = 1 
         
        endif
        widget_control, cdb_cdb_TE_label, set_value = 'Delta TE (ms):' 
      END
      
      cdb_fp_button: BEGIN
      info.free_precession_flag = 1
        if info.scan1 gt -1 and info.scan2 gt -1 then begin
          delay = mas_extract_delays(project.imndarray[info.scan1].file_path, project.imndarray[info.scan1].state1_load_procedure)
          
          info.TEx = delay
          delay = mas_extract_delays(project.imndarray[info.scan2].file_path, project.imndarray[info.scan2].state1_load_procedure)
         
          info.TEy = delay
          info.dTE = (info.TEx - info.TEy)*1E3
          TE_text = strcompress(string(info.dTE),/remove_all)
          widget_control, cdb_cdb_TE_text, set_value = TE_text
          widget_control, cdb_cd_button, sensitive = 1 
         
        endif
       widget_control, cdb_cdb_TE_label, set_value = 'Delta FP (ms):'
      END
      
      cdb_loadx_button: BEGIN
        mas_load_state_1
        ;want to make sure we grab the data in the orientation the user is expecting
        if slice_axis eq 0 then begin
          project.algebra1 = ptr_new(reform((*project.dataArray[CI].state1)[*,*,sdim_start,adim_start]))
        endif else if slice_axis eq 1 then begin
          project.algebra1 = ptr_new(reform((*project.dataArray[CI].state1)[*,pdim_start,*,adim_start]))
        endif else if slice_axis eq 2 then begin
          project.algebra1 = ptr_new(reform((*project.dataArray[CI].state1)[fdim_start,*,*,adim_start]))
        endif
        project.algebra1_3D = ptr_new(reform((*project.dataArray[CI].state1)[*,*,*,adim_start]))
        info.scan1 = CI
        if ptr_valid(project.IMNDArray[CI].echo_Time_ptr) then te1 = *project.IMNDArray[CI].echo_Time_ptr[0] else $
          te1 = project.IMNDArray[CI].echo_Time
          if info.free_precession_flag eq 0 then begin
            adim = project.procPramArray[CI].adim_start
            if n_elements(te1) ge adim then info.TEx = te1[adim] else info.TEx = te1
          endif else begin
            delay = mas_extract_delays(project.imndarray[CI].file_path,project.imndarray[CI].state1_load_procedure)
            ;delay 16 is the free precession delay
            info.TEx = delay
          endelse
        ;grab scan name to be displayed in GUI
        scan_name_temp = project.scan_list[CI]
        scan_name_split = strsplit(scan_name_temp,PATH_SEP(),/EXTRACT)
        sz_split = size(scan_name_split)
        scan_name = PATH_SEP() + scan_name_split[sz_split[1]-2] + PATH_SEP() + scan_name_split[sz_split[1]-1] + PATH_SEP()
        info.scan1_name = scan_name
        widget_control, cdb_loadx_text, set_value = scan_name  
        ;if both scans are loaded we want to calculate the delta TE, display it, and activate the complex division buttons
        if info.scan2 gt -1 then begin
          if info.free_precession_flag EQ 1 then begin
            delay = mas_extract_delays(project.imndarray[info.scan1].file_path, project.imndarray[info.scan1].state1_load_procedure)
            info.TEx = delay
            delay = mas_extract_delays(project.imndarray[info.scan2].file_path, project.imndarray[info.scan2].state1_load_procedure)
            info.TEy = delay
            info.dTE = (info.TEx - info.TEy)*1E3
          endif else begin
           info.dTE = info.TEx-info.TEy
          endelse 
          TE_text = strcompress(string(info.dTE),/remove_all)
          widget_control, cdb_cdb_TE_text, set_value = TE_text
          widget_control, cdb_cd_button, sensitive = 1 
          dim  = project.imndArray[info.scan1].dimensions
          print, dim
          if dim gt 2 then begin
          widget_control, cdb_cd3D_button, sensitive = 1 
          endif else begin
          widget_control, cdb_cd3D_button, sensitive = 0
          endelse
        endif
        info.complex_division = 0
        widget_control, cdb_cd_display_button, sensitive = 0
        widget_control, cdb_display_button, sensitive = 0
      END
      
       cdb_loady_button: BEGIN
        mas_load_state_1
        ;want to make sure we grab the data in the orientation the user is expecting
        if slice_axis eq 0 then begin
          project.algebra2 = ptr_new(reform((*project.dataArray[CI].state1)[*,*,sdim_start,adim_start]))
        endif else if slice_axis eq 1 then begin
          project.algebra2 = ptr_new(reform((*project.dataArray[CI].state1)[*,pdim_start,*,adim_start]))
        endif else if slice_axis eq 2 then begin
          project.algebra2 = ptr_new(reform((*project.dataArray[CI].state1)[fdim_start,*,*,adim_start]))
        endif
        project.algebra2_3D = ptr_new(reform((*project.dataArray[CI].state1)[*,*,*,adim_start]))
        info.scan2 = CI
        if ptr_valid(project.IMNDArray[CI].echo_Time_ptr) then te2 = *project.IMNDArray[CI].echo_Time_ptr[0] else $
          te2 = project.IMNDArray[CI].echo_Time
         
        if info.free_precession_flag eq 0 then begin
            adim = project.procPramArray[CI].adim_start
            if n_elements(te2) ge adim then info.TEy = te2[adim] else info.TEy = te2
          endif else begin
            delay = mas_extract_delays(project.imndarray[CI].file_path,project.imndarray[CI].state1_load_procedure)
            ;delay 16 is the free precession delay
            info.TEy = delay
          endelse
        ;grab scan name to be displayed in GUI
        scan_name_temp = project.scan_list[CI]
        scan_name_split = strsplit(scan_name_temp,PATH_SEP(),/EXTRACT)
        sz_split = size(scan_name_split)
        scan_name = PATH_SEP() + scan_name_split[sz_split[1]-2] + PATH_SEP() + scan_name_split[sz_split[1]-1] + PATH_SEP()
        info.scan2_name = scan_name
        widget_control, cdb_loady_text, set_value = scan_name  
        ;if both scans are loaded we want to calculate the delta TE, display it, and activate the complex division buttons
        if info.scan1 gt -1 then begin
          if info.free_precession_flag EQ 1 then begin
            delay = mas_extract_delays(project.imndarray[info.scan1].file_path, project.imndarray[info.scan1].state1_load_procedure)
            info.TEx = delay
            delay = mas_extract_delays(project.imndarray[info.scan2].file_path, project.imndarray[info.scan2].state1_load_procedure)
            info.TEy = delay
            info.dTE = (info.TEx - info.TEy)*1E3
          endif else begin
            info.dTE = info.TEx-info.TEy
          endelse 
          TE_text = strcompress(string(info.dTE),/remove_all)
          widget_control, cdb_cdb_TE_text, set_value = TE_text
          widget_control, cdb_cd_button, sensitive = 1 
          widget_control, cdb_cd3D_button, sensitive = 1 
          dim  = project.imndArray[info.scan1].dimensions
          print, dim
          if dim gt 2 then begin
          widget_control, cdb_cd3D_button, sensitive = 1 
          endif else begin
          widget_control, cdb_cd3D_button, sensitive = 0
          endelse 
        endif
        info.complex_division = 0
        widget_control, cdb_cd_display_button, sensitive = 0
        widget_control, cdb_display_button, sensitive = 0
      END
      
      cdb_cd_button: BEGIN
        if ptr_valid(project.algebra1)eq 0 or ptr_valid(project.algebra2) eq 0 then begin
                  update_status_bar,'please load X and Y with data'
                  return
        endif
        
        if project.procPramArray[info.scan1].signal_type NE 9 or project.procPramArray[info.scan1].signal_type NE 9 then begin
                  update_status_bar,'please load X and Y with complex data'
                  return
        endif
        ;check to make sure algebra 1 and algebra 2 are the same size
        size_1 = size(*project.algebra1)
        size_2 = size(*project.algebra2)
        
        if size_1[1] NE size_2[1] or size_1[2] NE size_2[2] then begin
                  update_status_bar, 'X and Y are not the same size'
                  return
        endif
        
        project.algebra_result = ptr_new((*project.algebra1/*project.algebra2))
        info.complex_division = 1
        widget_control, cdb_cd_display_button, sensitive = info.complex_division
        widget_control, cdb_display_button, sensitive = 1
      END
      
      cdb_cd3D_button: BEGIN
        if ptr_valid(project.algebra1)eq 0 or ptr_valid(project.algebra2) eq 0 then begin
                  update_status_bar,'please load X and Y with data'
                  return
        endif
        
        if project.procPramArray[info.scan1].signal_type NE 9 or project.procPramArray[info.scan1].signal_type NE 9 then begin
                  update_status_bar,'please load X and Y with complex data'
                  return
        endif
         ;check to make sure algebra 1 and algebra 2 are the same size
        size_1 = size(*project.dataArray[info.scan1].state1)
        size_2 = size(*project.dataArray[info.scan2].state1)
        
        if size_1[1] NE size_2[1] or size_1[2] NE size_2[2] or size_1[3] NE size_2[3] then begin
                  update_status_bar, 'X and Y are not the same size'
                  return
        endif
        project.algebra_result_3D = ptr_new((*project.algebra1_3D/*project.algebra2_3D))
        info.complex_division_3D = 1
        
        widget_control, cdb_cd3d_display_button, sensitive = info.complex_division_3D
        widget_control, cdb_3d_display_button, sensitive = 1
      END
      
      cdb_filter_toggle:BEGIN
        info.mask = event.select
        widget_control, cdb_threshold_text, sensitive = info.mask
      END
      
      cdb_unwrap_toggle:BEGIN
        info.phase_unwrap = event.select
      END
      
      cdb_threshold_text: BEGIN
        widget_control, cdb_threshold_text, get_value = thresh_temp
        info.threshold = thresh_temp
        info.mask_created = 0
      END
      cdb_hipass_toggle: BEGIN
        info.hipass = event.select        
        END
        
      cdb_cd_display_droplist: BEGIN
        selection = widget_info(cdb_cd_display_droplist, /DROPLIST_SELECT)
        info.image_type = selection
        
       END
      
      cdb_cd_display_button: BEGIN
        info.display_3D = 0
        mas_calc_delta_b_display, info
      END
        
      cdb_cd3D_display_button: BEGIN
        info.display_3D = 1
        mas_calc_delta_b_display, info
       
      END  
      
      cdb_cdb_gr_text:BEGIN
        widget_control, cdb_cdb_gr_text, get_value=gr_temp
        info.gyromagnetic_ratio = gr_temp
      END
      
      cdb_display_button:BEGIN
        widget_control, cdb_cdb_gr_text, get_value=gr_temp
        info.gyromagnetic_ratio = gr_temp
        info.display_3D = 0
        mas_calc_delta_b, info       
      END
      
      cdb_3D_display_button:BEGIN
        widget_control, cdb_cdb_gr_text, get_value=gr_temp
        info.gyromagnetic_ratio = gr_temp
        info.display_3D = 1
        mas_calc_delta_b, info
      END
    ENDCASE
widget_control, cdb_cd_display_button, sensitive = info.complex_division 
widget_control, cdb_display_button, sensitive = info.complex_division
widget_control, event.top, set_uvalue=info  
end
; Subroutine name: mas_calc_delta_b
; Created by: Garrett Astary
; Calling Information: Called from mas.pro when user selects Calculate Delta B (widget button uname = W_MENU_calculate_delta_b)  
;                      from Analyze dropdown

; Bugs or Important Comments to Developers:


; Purpose of subroutine: Creates widget base and widgets that allow user to complex data in order to display
;                        magnitude/real/imaginary/phase images of complex division result and to calculate the change in B field around and object. 
;                        The data being manipulated must be in complex form (Image -> type -> Complex)
;                        See Haacke et al. "Magnetic Resonance Imaging Physical Principles and Sequence Design"
;                        Chapter 25 on Susceptibility Imaging pg. 759

; Editing Information:
pro mas_calc_delta_b_algebra_gui
    HEAP_GC
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets
    CI = project.CI
    
    
    cdb_window_base = widget_base(column=1, TITLE='Calculate Delta B', xoffset=420, /base_align_center)
    
    cdb_note_base = widget_base(cdb_window_base)
    cdb_note_label = widget_label(cdb_note_base, value='Note: Please use complex image type')
    
    cdb_type_base = widget_base(cdb_window_base, /exclusive)
    cdb_ge_button = widget_button(cdb_type_base, value = 'Gradient Echo', uname = 'cdb_ge_button')
    cdb_fp_button = widget_button(cdb_type_base, value = 'Free Precession Delay', uname = 'cdb_fp_button')
    widget_control, cdb_ge_button, set_button = 1
    
    cdb_load_base = widget_base(cdb_window_base)
    
    cdb_loadx_base = widget_base(cdb_load_base, column=3)
    cdb_loadx_button = widget_button(cdb_loadx_base, value = 'Load X')
    cdb_loadx_label = widget_label(cdb_loadx_base, value = 'Scan:')
    cdb_loadx_text = widget_text(cdb_loadx_base, value = 'Empty', xsize = 30, editable = 0)
    
    cdb_load2_base = widget_base(cdb_window_base)
    
    cdb_loady_base = widget_base(cdb_load2_base, column=3)
    cdb_loady_button = widget_button(cdb_loady_base, value = 'Load Y')
    cdb_loady_label = widget_label(cdb_loady_base, value = 'Scan:')
    cdb_loady_text = widget_text(cdb_loady_base, value = 'Empty', xsize = 30, editable = 0)
    
    
    cdb_complex_division_label = widget_label(cdb_window_base, value = 'Complex Division')
    cdb_complex_division_base = widget_base(cdb_window_base, column=2)
    cdb_cd_button = widget_button(cdb_complex_division_base, value = 'Complex Divide X/Y', sensitive=0)
    cdb_cd3D_button = widget_button(cdb_complex_division_base, value = '3D Complex Divide', sensitive=0)
    
    cdb_display_base_label = widget_label(cdb_window_base, value = 'Display Options')
    
    cdb_filter_base = widget_base(cdb_window_base, column = 2)
    cdb_filter_toggle_base = widget_base(cdb_filter_base, /nonexclusive)
    cdb_filter_toggle = widget_button(cdb_filter_toggle_base, value = 'Mask')
    cdb_threshold_base = widget_base(cdb_filter_base, column = 2)
    cdb_threshold_label = widget_label(cdb_threshold_base, value = 'Threshold:')
    cdb_threshold_text = widget_text(cdb_threshold_base, value='0.05', xsize=5, /editable, sensitive = 0)
    cdb_unwrap_toggle_base = widget_base(cdb_filter_base, /nonexclusive)
    cdb_unwrap_toggle = widget_button(cdb_unwrap_toggle_base, value = 'Enable Phase Unwrapping')
    cdb_hipass_toggle_base = widget_base(cdb_filter_base, /nonexclusive)
    cdb_hipass_toggle = widget_button(cdb_hipass_toggle_base, value = 'High pass filter')
    
    cdb_cd_display_base = widget_base(cdb_window_base, column = 4)
    cdb_cd_display_label = widget_label(cdb_cd_display_base, value = 'Image to Display:')
    choices = ['Magnitude', 'Real', 'Imaginary', 'Phase']
    cdb_cd_display_droplist = widget_droplist(cdb_cd_display_base, value = choices, uname='cdb_cd_display_droplist')
    cdb_cd_display_button = widget_button(cdb_cd_display_base, value = 'Display', sensitive=0)
    cdb_cd3D_display_button = widget_button(cdb_cd_display_base, value = 'Ortho Slice Display', sensitive=0)
    
    cdb_calculate_delta_B_label = widget_label(cdb_window_base, value = 'Calculate Delta B')
    cdb_cdb_base = widget_base(cdb_window_base, row =2 )
    cdb_cdb_gr_label = widget_label(cdb_cdb_base, value='Gyromagnetic Ratio (MHz/T):')
    cdb_cdb_gr_text = widget_text(cdb_cdb_base, value='42.58', editable=0, xsize = 5)
    cdb_cdb_TE_label = widget_label(cdb_cdb_base, value = 'Delta TE (ms):', uname = 'cdb_cdb_TE_label')
    cdb_cdb_TE_text = widget_text(cdb_cdb_base, value = '', xsize = 6, editable=0)
  
    
    cdb_action_base = widget_base(cdb_window_base, column = 2)
    cdb_display_button = widget_button(cdb_action_base, value = 'Display Delta B', sensitive=0)
    cdb_3D_display_button = widget_button(cdb_action_base, value = 'Display Delta B 3D', sensitive=0)

    ;Parameter structure
    ;info = {scan1: X data or scan with longer TE
    ;        scan2: Y data or scan with shorter TE
    ;        scan1_name: string specifying X data name to be displayed in GUI
    ;        scan2_name: string specifying Y data name to be displayed in GUI
    ;        TEx: Echo time used in X data acquisition (ms)
    ;        TEy: Echo time used in Y data acquisition (ms)
    ;        dTE: delta TE in ms (TEx - TEy)
    ;        dTE_flag: flag determining if dTE can be calculated (= 1 if X and Y contain data)
    ;        gyromagnetic_ratio: gyromagnetic ratio in MHz/Tesla (user input, defualt setting is for 1H)
    ;        image_type: 0 - magnitude, 1 - Real, 2 - Imaginary, 3 - Phase, used for displaying result of complex division
    ;        complex_division: flag denoting if complex division has been calculated (initially 0)
    ;        complex_division_3D:  flag denoting if 3D complex division has been calculated (initially 0)
    ;        display_3D: Flag, initially set to 0, set to 1 when user selects Ortho slice viewer as display option
    ;        mask: Flag, set to 0, set to 1 by user selecting 'Mask' checkbox
    ;        mask_created: Flag, set to 0 until mask is generated, reset to 0 when threshold is changed (avoids calculating 
    ;                      the mask mulitple times when the same threshold is being used)
    ;        threshold: Initially set to 0.05 this parameter is changed when the user changes the input to the threshold textbox
    ;                   The threshold is multiplied by the max signal in the magnitude mode of the images, anything below this 
    ;                   threshold is masked out of final display routines
    
    info = {scan1:-1, scan2:-1, scan1_name:'', scan2_name:'', TEx:0.0, TEy:0.0, dTE :0.0, dTE_flag :0, gyromagnetic_ratio:0.0 $
    ,image_type:0, complex_division:0,complex_division_3D:0, display_3D:0, mask:0, mask_created:0, threshold: 0.05, free_precession_flag:0, phase_unwrap:0,hipass:0}
    
    widget_control, cdb_window_base, set_uvalue=info, /no_copy
     widget_control, cdb_window_base, /realize
    xmanager, 'mas_calc_delta_b_algebra_gui', cdb_window_base,/NO_BLOCK

end