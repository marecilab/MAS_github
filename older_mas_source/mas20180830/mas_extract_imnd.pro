;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved


function extract_bruker_onepulse, header_string

    common scan_data, project
    forward_function extract_simple_value
    ni = project.ni
    
    my_header_string = header_string
    
    scan_directory = project.imndArray[project.ni].file_Path
    acqp_file = scan_directory+'acqp'
    if (file_test(acqp_file, /read)) then begin
        openr, lun, acqp_file, /get_lun
        tmpstr = ''
        while (not eof(lun)) do begin
            readf, lun, tmpstr
            my_header_string += tmpstr
        endwhile
    endif
    
    temp = extract_simple_value( my_header_string, '##$IMND_nucleus=', '##$', '$$' )
    if temp eq 0 then temp = extract_simple_value( my_header_string, '##$PVM_Nucleus1=( 8 )', '##$', '$$' )
    if temp eq 0 then temp = extract_simple_value( my_header_string, '##$NUC1=', '##$')
    project.imndarray[ni].spect_nucleus=strtrim(temp,2)
    
    temp = extract_simple_value( my_header_string, '##$IMND_sw_h=', '##$', '$$' )
    if temp eq 0 then temp = extract_simple_value( my_header_string, '##$PVM_SpecSWH=( 1 )', '##$', '$$' )
    if temp eq 0 then temp = extract_simple_value( my_header_string, '##$SW_h=', '#')
    project.imndarray[ni].spect_spectral_width=float(temp)
    
    temp = extract_simple_value( my_header_string, '##$IMND_spect_acq_size=', '##$', '$$' )
    if temp eq 0 then temp = extract_simple_value( my_header_string, '##$TD=', '#')
    project.imndarray[ni].spect_acq_size=long(temp)/2
    if temp eq 0 then begin
      temp = extract_simple_value( my_header_string, '##$PVM_SpecMatrix=( 1 )', '##$', '$$' )
      project.imndarray[ni].spect_acq_size=long(temp)
    endif 
    if temp eq 0 then begin
      temp = extract_simple_value( my_header_string, '##$PVM_EncMatrix=( 1 )', '##$', '$$' )
      project.imndarray[ni].spect_acq_size=long(temp)
    endif
    
    temp = extract_simple_value( my_header_string, '##$IMND_rep_time=', '##$', '$$' )
    if temp eq 0 then temp = extract_simple_value( my_header_string, '##$PVM_RepetitionTime=', '##$', '$$' )
    if temp eq 0 then begin
      temp = extract_simple_value( my_header_string, '##$D= (0..63)', '##$',long_value=1)
      a = strsplit(temp,/extract)
      temp = 1e3*a[1]
    endif
    project.imndarray[ni].spect_rep_time=float(temp)
    
    temp = extract_simple_value( my_header_string, '##$IMND_pulse_length=', '##$', '$$' )
    if temp eq 0 then begin
      temp = extract_simple_value( my_header_string, '##$P= (0..63)', '##$',long_value=1)
      if temp ne 0 then begin
        a = strsplit(temp,/extract)
        temp = a[1]
      endif 
    endif
    project.imndarray[ni].spect_pulse_length=float(temp)

    temp = extract_simple_value( my_header_string, '##$IMND_n_averages=', '##$', '$$' )
    if temp eq 0 then temp = extract_simple_value( my_header_string, '##$PVM_NAverages=', '##$', '$$' )
    if temp eq 0 then temp = extract_simple_value( my_header_string, '##$NS=', '#')
    project.imndarray[ni].spect_num_avg = long(temp)
    
    temp = extract_simple_value( my_header_string, '##$IMND_bf1=', '##$', '$$' )
    if temp eq 0 then temp = extract_simple_value( my_header_string, '##$PVM_FrqRef=( 8 )', '##$', '$$' )
    if temp eq 0 then temp = extract_simple_value( my_header_string, '##$BF1=', '#')
    project.imndarray[ni].spect_bf1=float(temp)

    temp = extract_simple_value( my_header_string, '##$IMND_acq_time=', '##$', '$$' )
    if temp eq 0 then temp = extract_simple_value( my_header_string, '##$PVM_SpecAcquisitionTime=', '##$', '$$' )
    if temp eq 0 then begin
      temp = 1e3*project.imndarray[ni].spect_acq_size/project.imndarray[ni].spect_spectral_width
    endif
    project.imndarray[ni].spect_acq_time=float(temp)
    
    temp = extract_simple_value( header_string ,'##$PVM_NRepetitions=','##$', '$$' )
    if temp eq 0 then temp = extract_simple_value( my_header_string, '##$NR=', '##$', '$$' )
    if temp gt 0 then  project.imndarray[ni].adim=long(temp)
  
  
    project.imndarray[ni].image_type = 99
    
    project.imndarray[ni].state1_load_procedure = 'mas_load_state_1_bruker_onepulse'
    
    return, 1

end

function extract_bruker_config

  common scan_data, project
  forward_function extract_simple_value
  ni = project.ni

  my_header_string = ''
  scan_directory = project.imndArray[project.ni].file_Path
  acqp_file = scan_directory+'configscan'
  if (file_test(acqp_file, /read)) then begin
    openr, lun, acqp_file, /get_lun
    tmpstr = ''
    while (not eof(lun)) do begin
      readf, lun, tmpstr
      my_header_string += tmpstr
    endwhile
  endif
  temp = extract_simple_value( my_header_string, '##$IMND_nucleus=', '##$', '$$' )
  project.imndArray[ni].gcoil = temp
  
end

; Subroutine name: extract_simple_value
; Created by:
; Calling Information:
; parameters:
;   string_to_search  [in] the string to search for the smaller string
;   start_string      [in]   the string to start search for
;   stop_string      [in]    the string to stop searching for
;   retruns          [out]  the number between start_string and stop_string

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Use to give a text value out of an imnd file.

; Editing Information:
    ;Edited by HS 2006/09/30.
    ;Fix spelling mistakes and commenting

FUNCTION extract_simple_value ,string_to_search ,start_string, stop_string, stop_string2, long_value=long_value
    COMPILE_OPT IDL2

    index1 = Strpos(string_to_search, start_string)
    ; if i cannot find the string to search, quit and return zero
    if index1 eq -1 then RETURN, '0'

    index1 =  index1 + STRLEN(start_string)
    index2 = Strpos(string_to_search, stop_string , index1 )
    difference =  index2 - index1

    ;print, 'index',+ index1
    ;print, 'difference: ',+ difference

    ; try and search for the alternate string as well too.
    if (difference GT 50 and N_ELEMENTS(stop_string2)) then begin

         index1 = Strpos(string_to_search, start_string)
           index1 =  index1 + STRLEN(start_string)
           index2 = Strpos(string_to_search, stop_string2 , index1 )
           difference =  index2 - index1
    end

    if (not keyword_set(long_value)) then begin
        IF (difference GT 50) then RETURN, '0'
    endif

    RETURN, STRMID(string_to_search, index1, difference)

END



; Subroutine name: extract_slice_scheme
; Created by:
; Calling Information:
;  parameters:
;   string_to_search =  the string to search for the substring
;   imnd_search_pram =  the specific parameter before the array
;   slice = adim
;

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;   This will search an imnd string to look for a specific sub string that contains
;   the string imnd_search_pram then the number array after that it will extract that
;   array and return it as a float array. if there is an error it will return [-1,-1] array
;   that way  you can always check the first element to see if it is equal to -1 and there should
;   always be at least 2 elements also. if your data has a -1 for the first parameter then you
;   can check to see if it is supposed to have -1 -1 for the first 2 entries.
;

; Editing Information:
    ;Edited by HS 2006/09/30.
    ;Fix spelling mistakes and commenting
    ; Edited by Magdoom 2015/04/12
    ; Included PVM_* search keywords
    
FUNCTION extract_slice_scheme, string_to_search, imnd_search_pram, slice
    COMPILE_OPT IDL2

    imnd_search_pram = imnd_search_pram + '=('
    index1 = Strpos(string_to_search, imnd_search_pram )
    if index1 eq -1 then begin
       temp = intArr(2)
       temp = [-1, -1]
       return, temp
    end
    index1 = Strpos(string_to_search, ' )', index1)+2
    ;index1 should be just before the array

    temp1 = Strpos(string_to_search, '##$IMND' , index1 )
    temp2 = Strpos(string_to_search, '##$PVM' , index1 )
    if temp2 ne -1 and temp2 lt temp1 then index2 = temp2 else index2 = temp1
    
    ;index2 should be just after the array
    ;If it is not found then most likely we are dealing with a METHOD file.
    if index2 eq -1 then index2 = Strpos(string_to_search, '$$' , index1 )

    difference =  index2 - index1
;   print, 'in extract slice scheme, difference',difference
    ;print, 'index',+ index1
    ;print, 'difference: ',+ difference
    temp = STRMID(string_to_search, index1, difference)

    slice_scheme = (Strsplit(temp, ' ',/EXTRACT))[0:slice-1]
    slice_scheme = Float(slice_scheme)
    RETURN, slice_scheme

END


; Subroutine name: extract_slice_scheme_methodfiles
; Created by:
; HS on 20061013
;
; Calling Information:
;  parameters:
;   string_to_search =  the string to search for the substring
;   search_pram =  the specific parameter before the array
;   slice = adim
;

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Sames as extract_slice_scheme but specific for METHOD files
;

; Editing Information:
    ;Edited by HS 2006/10/13


FUNCTION extract_slice_scheme_methodfiles, string_to_search, search_pram, slice
    COMPILE_OPT IDL2

    search_pram = search_pram + '=('
    index1 = Strpos(string_to_search, search_pram )
    if index1 eq -1 then begin
       temp = intArr(2)
       temp = [-1, -1]
       return, temp
    end
    index1 = Strpos(string_to_search, ' )', index1)+2
    ;index1 should be just before the array

    index2 = Strpos(string_to_search, '##$PVM' , index1 )
    ;index2 should be just after the array

    difference =  index2 - index1
;   print, 'in extract slice scheme, difference',difference
    ;print, 'index',+ index1
    ;print, 'difference: ',+ difference
    temp = STRMID(string_to_search, index1, difference)

    slice_scheme = (Strsplit(temp, ' ',/EXTRACT))[0:slice-1]
    slice_scheme = Float(slice_scheme)
    RETURN, slice_scheme

END


; Subroutine name: extract_phi_methods
; Created by: HS - 20061214
; Calling Information: Call it by sending the METHOD text in a variable.

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; To extract and convert xyz (rectangular coordinate directions) in METHOD files into
; Polar coordinate phi used by MAS.

; Editing Information:
; 20070504 CD switched phi to physics definitions
; 03/26/15 Fixed bug with method file in PV6 (Magdoom)

FUNCTION extract_phi_methods, string_to_search
    COMPILE_OPT IDL2

    ;;search_pram = '##$PVM_DwDir=('
    ;; This seems to be a better param to look at for 
    ;; diffusion gradient directions. It has an entry for
    ;; the B0.
    search_pram = '##$PVM_DwGradVec=('
    
    index1 = Strpos(string_to_search, search_pram )
    if index1 eq -1 then begin
       temp = intArr(2)
       temp = [-1, -1]
       return, temp
    end

    index1 = Strpos(string_to_search, ' )', index1)+2
    ;index1 should be just before the array
    
    temp1 = Strpos(string_to_search, '##$PVM', index1)
    temp2 = Strpos(string_to_search, '$$', index1)
    if temp2 ne -1 and temp2 lt temp1 then index2 = temp2 else index2 = temp1
     ;index2 should be just after the array

    difference =  index2 - index1
;   print, 'in extract slice scheme, difference',difference
    ;print, 'index',+ index1
    ;print, 'difference: ',+ difference
    temp = STRMID(string_to_search, index1, difference)

	temp = Strsplit(temp, ' ',/EXTRACT)

  temp = reform(temp, 3, n_elements(temp)/3)
  x = reform(float(temp[0,*]))
  y = reform(float(temp[1,*]))
  z = reform(float(temp[2,*]))
  
	rho = sqrt(x^2+y^2+z^2)

	phi = ATAN(y,x)*180/!PI

	; rotate negative angles by 360 degrees
	idx = WHERE((phi LT 0),zcount)
	IF zcount NE 0 THEN phi[idx] = phi[idx] + 360

	; phi   = arctan(y/x)
	; theta = arccos(z/sqrt(x^2+y^2+z^2))

	return, string(phi)

end

; Subroutine name: extract_theta_methods
; Created by: HS - 20061214
; Calling Information: Call it by sending the METHOD text in a variable.

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; To extract and convert xyz (rectangular coordinate directions) in METHOD files into
; Polar coordinate theta used by MAS.

; Editing Information:
; 20070504 CD switched theta to physics definitions
; 03/26/15 Fixed bug with method file in PV6 (Magdoom)

FUNCTION extract_theta_methods, string_to_search
    COMPILE_OPT IDL2

    ;; search_pram = '##$PVM_DwDir=('
    ;; This is a better param to look at for gradient directions.
    search_pram = '##$PVM_DwGradVec=('
    
    index1 = Strpos(string_to_search, search_pram )
    if index1 eq -1 then begin
       temp = intArr(2)
       temp = [-1, -1]
       return, temp
    end
    index1 = Strpos(string_to_search, ' )', index1)+2
    ;index1 should be just before the array

    temp1 = Strpos(string_to_search, '##$PVM', index1)
    temp2 = Strpos(string_to_search, '$$', index1)
    if temp2 ne -1 and temp2 lt temp1 then index2 = temp2 else index2 = temp1
    ;index2 should be just after the array

    difference =  index2 - index1
;   print, 'in extract slice scheme, difference',difference
    ;print, 'index',+ index1
    ;print, 'difference: ',+ difference
    temp = STRMID(string_to_search, index1, difference)

	temp = Strsplit(temp, ' ',/EXTRACT)

        temp = reform(temp, 3, n_elements(temp)/3)

	x = reform(float(temp[0,*]))
	y = reform(float(temp[1,*]))
	z = reform(float(temp[2,*]))
  
	rho = sqrt(x^2+y^2+z^2)
	theta = acos(z/rho)*180/!pi

  ;; get rid of nans caused by B0
  nans = where(finite(theta, /nan), n_nans)
  if (n_nans gt 0) then theta[nans] = 0.0
  
	; phi   = arctan(y/x)
	; theta = arccos(z/sqrt(x^2+y^2+z^2))

	return, string(theta)

end


; Subroutine name: extract_bval_method
; Created by: CD
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/09/30.

FUNCTION extract_bval_method, string_to_search
    COMPILE_OPT IDL2

    ;;search_pram = '##$PVM_DwBvalEach=( 1 )'
    search_pram = '##$PVM_DwEffBval=( 1 )'

    index1 = STRPOS(string_to_search, search_pram )
    IF index1 ne -1 THEN begin

       index1 = index1 + STRLEN(search_pram)
       index2 = Strpos(string_to_search, '##$PVM' , index1)
       difference =  index2 - index1
       
       RETURN, FIX(STRMID(string_to_search, index1, difference))

    ENDIF ELSE BEGIN
       
       search_pram = '##$PVM_DwEffBval=('
       index1 = Strpos(string_to_search, search_pram )
       if index1 eq -1 then begin
          temp = intArr(2)
          temp = [-1, -1]
          return, -1 ;;temp
       end
       index1 = Strpos(string_to_search, ' )', index1)+2
       index2 = Strpos(string_to_search, '##$PVM' , index1 )
       index3 = Strpos(string_to_search, '$$', index1)
       if (index3 gt 0) then begin
          index_min = min([index2, index3])
       endif else begin
          index_min = index2
       endelse
       difference =  index_min - index1
       temp = STRMID(string_to_search, index1, difference)
       
       temp = Strsplit(temp, ' ',/EXTRACT)

       return, temp

    ENDELSE
end


; Subroutine name: extract_big_delta_method
; Created by: CD
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/09/30.

FUNCTION extract_big_delta_method, string_to_search
    COMPILE_OPT IDL2

    search_pram = '##$PVM_DwGradSep=( 1 )'
    index1 = STRPOS(string_to_search, search_pram )
    IF index1 EQ -1 THEN RETURN, -1

    index1 = index1 + STRLEN(search_pram)
    index2 = Strpos(string_to_search, '##$PVM' , index1)
    difference =  index2 - index1

    RETURN, FLOAT(STRMID(string_to_search, index1, difference))
 end


; Subroutine name: extract_small_delta_method
; Created by: CD
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/09/30.

FUNCTION extract_small_delta_method, string_to_search
    COMPILE_OPT IDL2
    
    search_pram = '##$PVM_DwGradDur=( 1 )'
    index1 = STRPOS(string_to_search, search_pram )
    IF index1 EQ -1 THEN RETURN, -1
    
    index1 = index1 + STRLEN(search_pram)
    index2 = Strpos(string_to_search, '##$PVM' , index1)
    difference =  index2 - index1
    
    RETURN, FLOAT(STRMID(string_to_search, index1, difference))
end

; Subroutine name: extract_grad_mat_method
; Created by: BT
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:

FUNCTION extract_grad_mat_method, string_to_search
    COMPILE_OPT IDL2

    ;;search_pram = '##$PVM_DwGradVec=( 1, 3 )'
    search_pram = '##$PVM_DwDir=( 1, 3 )'
    ;;search_pram = '##$PVM_DwSpDir=( 1, 3 )'
    index1 = Strpos(string_to_search, search_pram )
    if index1 eq -1 then begin
       temp = intArr(2)
       temp = [-1, -1]
       return, temp
    end

    index1 = Strpos(string_to_search, ' )', index1)+2
    ;index1 should be just before the array

    index2 = Strpos(string_to_search, '##$PVM' , index1 )
    ;index2 should be just after the array

    difference =  index2 - index1
;   print, 'in extract slice scheme, difference',difference
                                ;print, 'index',+ index1
                                ;print, 'difference: ',+ difference
    temp = STRMID(string_to_search, index1, difference)
    
    temp = Strsplit(temp, ' ',/EXTRACT)
    x = float(temp[0])
    y = float(temp[1])
    z = float(temp[2])

    ;;temp = extract_simple_value( string_to_search , search_pram,'##$', '$$')
    ;;xyz = strsplit(temp, /extract)
    xyz = float([x,y,z])
    print, xyz
    RETURN, xyz
end

; Subroutine name: extract_grad_mat_imad
; Created by: BT
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:

FUNCTION extract_grad_mat_imnd, string_to_search
    COMPILE_OPT IDL2

    search_pram = '##$IMND_diff_grad_mat=( 1, 3 )'
    temp = extract_simple_value( string_to_search , search_pram,'##$', '$$')
    xyz = strsplit(temp, /extract)
    xyz = float(xyz)
    RETURN, xyz
end

; Subroutine name: extract_bval_imnd
; Created by: CD
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/09/30.

FUNCTION extract_bval_imnd, string_to_search
    COMPILE_OPT IDL2

    search_pram = '##$IMND_diff_b_value=( 1 )'
    index1 = STRPOS(string_to_search, search_pram )
    IF index1 EQ -1 THEN RETURN, -1

	index1 = index1 + STRLEN(search_pram)
    index2 = Strpos(string_to_search, '##$IMND' , index1)
    difference =  index2 - index1
    
    RETURN, FIX(STRMID(string_to_search, index1, difference))
end


; Subroutine name: extract_name
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/09/30.
    ;Fix spelling mistakes and commenting
    ;HS 20061016. New way to extract the name of a sequence for METHOD files.

FUNCTION extract_name, string_to_search
    COMPILE_OPT IDL2

; This next block of code
    ;index1 = Strpos(string_to_search, '<')
    ;index2 = Strpos(string_to_search, '>',  index1 + 1)
    ;difference = index2 - index1
    ;temp = STRMID(string_to_search, index1, difference)
    ;RETURN, STRMID(temp,1,difference)

; Start the new way to extract the name of sequence

    index1 = Strpos(string_to_search, '##$IMND_method=(')
    
    if index1 eq -1 then begin ;If it can't find this argument in the file it means we have a METHOD file

    index1 = Strpos(string_to_search, '##$Method=')
; 
; This is how it appears on the METHOD file:
; ##$Method=RARE
; ##$PVM_EchoTime=12.000
; We are looking for the first part of the name.
; Then we have to add the size of ##$Method= minus 1 to grab the correct part of the name.
;
    index1 = index1+Strlen('##$Method=')
    index2 = Strpos(string_to_search, '##$', index1)
    difference = index2-index1
    temp = STRMID(string_to_search, index1, difference)

    endif else begin ;This part happens when we are working with an IMND file.
    ;This is how it appears in the IMND file
    ;##$IMND_method=( 20 )
    ;<DWI_SE>
    ;##$IMND_dummy_method=DWI_SE

    index1 = Strpos(string_to_search, '<') + 1
    index2 = Strpos(string_to_search, '>',  index1+1)
    difference = index2 - index1
    temp = STRMID(string_to_search, index1, difference)

    endelse
    index = Strpos(string_to_search, '##$PULPROG=')
    if index ne -1 then temp = extract_simple_value(string_to_search ,'##$PULPROG=','#')
RETURN, temp


END


; Subroutine name: seconds_since_beginning_of_year
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This function will take in the ascii date and return seconds since beginning of year

; Editing Information:
    ;Edited by HS 2006/09/30.
    ;Fix spelling mistakes and commenting
;
function seconds_since_beginning_of_year, ascii_date

    temp = BIN_DATE(ascii_date)
    ;print, temp

    factor = 0
    if temp[1] eq 1 then factor = 0
    if temp[1] eq 2 then factor = 31
    if temp[1] eq 3 then factor = 59
    if temp[1] eq 4 then factor = 90
    if temp[1] eq 5 then factor = 120
    if temp[1] eq 6 then factor = 151
    if temp[1] eq 7 then factor = 181
    if temp[1] eq 8 then factor = 212
    if temp[1] eq 9 then factor = 243
    if temp[1] eq 10 then factor = 273
    if temp[1] eq 11 then factor = 304
    if temp[1] eq 12 then factor = 334

    month = DOUBLE(factor*86400)
    ;print, month
    day = double((temp[2]-1)*86400)
    ;print, day
    hour = double(temp[3]*3600)
    ;print, hour
    minute = double(temp[4]*60)
    ;print, minute
    second = double(temp[5])
    ;print, second

    return, month+day+hour+minute+second
end


; Subroutine name: extract_Dynamics_time
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:
;   This procedure will take 2 elements from each pre#&pst# imnd file. the 2 elements are
;   time1 =  will be the time the acquisition started
;   time2 =  will be the time the acquisition took to complete
;   after these two elements are collected we will take time2/2 + time1
;   and at that moment is when the image was collected.
;   the last scan will be used to define the last element of the time array
;   the last element will be time1+time2.
;   time = 0 will be time1 from pst1
;   the pre images will be less than 0 in time taking -(time2/2+time1)...i think
;
; Editing Information:
    ;Edited by HS 2006/09/30.
    ;Fix spelling mistakes and commenting

pro extract_Dynamics_time
    COMPILE_OPT IDL2
    COMMON scan_data

    ;get the proper slash depending upon which operation system is in use.
    slash = get_dir_slash()

    ci = project.ni
    multi_scan_file_array = *project.imndArray[ci].multi_scan_file_array
    n_Pre = project.imndArray[ci].n_Pre

    ;store the adim
    adim = (size(multi_scan_file_array))[1]
    project.imndArray[ci].adim = adim

    ;create the time arrays.
    scan_start = dblarr(adim)
    scan_length_str = strarr(adim)
    scanLengthSeconds = dblarr(adim)

    FOR ii=0, adim-1 DO BEGIN

       ;now open the acqp and imnd file and take all the data out of it and store it to a string.
       ;if done compleatly then we are preaty sure that this data is ok so we will increase the
       acqp_file_string = ''
       acqp_temp_string = ''
       OPENR, acqp_file_descriptor, multi_scan_file_array[ii]+slash[1]+'acqp' , /GET_LUN , ERROR = err
       ; If err is nonzero, something happened. Print the error message to
       ; the standard error file (logical unit -2):
       IF (err NE 0) then  begin
         ;PRINTF, -2, !ERROR_STATE.MSG
         update_status_bar,'Please open the directory containing the acqp files'
         return
       end
       ;if there was no error then the file's text can be extracted.
       While (eof(acqp_file_descriptor) NE 1) Do Begin
         Readf, acqp_file_descriptor, acqp_temp_string
         acqp_file_string = acqp_file_string + acqp_temp_string
       EndWhile
       close, acqp_file_descriptor
       FREE_LUN, acqp_file_descriptor


       imnd_file_string = ''
       imnd_temp_string = ''
       OPENR, imnd_file_descriptor, multi_scan_file_array[ii]+slash[1]+'imnd' , /GET_LUN , ERROR = err
       ; If err is nonzero, something happened. Print the error message to
       ; the standard error file (logical unit -2):
       IF (err NE 0) then  begin
         ;PRINTF, -2, !ERROR_STATE.MSG
         update_status_bar,'Please open the directory containing the acqp files'
         return
       end
       ;if there was no error then the file's text can be extracted.
       While (eof(imnd_file_descriptor) NE 1) Do Begin
         Readf, imnd_file_descriptor, imnd_temp_string
         imnd_file_string = imnd_file_string + imnd_temp_string
       EndWhile
       close, imnd_file_descriptor
       FREE_LUN, imnd_file_descriptor

       ;scan_start is in seconds.
       scan_start[ii] = double(extract_simple_value( acqp_file_string, '##$ACQ_abs_time=','##$'))

       ;extract the time string <0h15m12s>
       index1 = Strpos(imnd_file_string, 'total_time')
       index1 = Strpos(imnd_file_string, '<',index1)
        index2 = Strpos(imnd_file_string, '>',  index1 + 1)+1
        difference = index2 - index1
        scan_length_str[ii] = STRMID(imnd_file_string, index1, difference)

    ENDfor

    ;break apart the bruker notation of how long the scan took into it's meaningful parts
    ;<0h15m12s>.;afer the scan length is broken apart we now need to convert it to seconds.

    for ii=0, adim-1 do begin
       ;print,scanLenghtArray
       temp = (strsplit(scan_length_str[ii],'h', ESCAPE='<>',/EXTRACT ))[0]
       scanLengthSeconds[ii] = double(temp)*3600
       ;print, scanHours
       temp = (strsplit(scan_length_str[ii],'m', ESCAPE='<>',/EXTRACT ))[0]
       temp = (strsplit(temp[0],'h', ESCAPE='<>',/EXTRACT ))[1]
       scanLengthSeconds[ii] = scanLengthSeconds[ii] + double(temp)*60
       ;print, scanMinuts
       temp = (strsplit(scan_length_str[ii],'s', ESCAPE='<>',/EXTRACT ))[0]
       temp = (strsplit(temp[0],'m', ESCAPE='<>',/EXTRACT ))[1]
       scanLengthSeconds[ii] = scanLengthSeconds[ii] + double(temp)
       ;print, scanLengthSeconds[ii]
    end


   ;after we have extracted all the scan lengths and scan_start we have to sort the
   ;multi_scan_array, scan_length and scan_start arrays in increasing time.

    ;project.imndArray[project.ni].multi_scan_file_array = ptr_new(parent_file_listing[SORT(echo_times)])

    project.imndArray[ci].multi_scan_file_array = ptr_new(multi_scan_file_array[sort(scan_start)])
    scan_start = scan_start[sort(scan_start)]

    scanLengthSeconds = scanLengthSeconds[sort(scan_start)]


    ;now start all scans at time zero so subtract the earlyist time.
    ;since scans start is marked at the end of the scan we have to subtract the length
     ; of the scan .
     scan_start -= scanLengthSeconds
    scan_start -= scan_start[0]
    ;print, scan_start

    ;i assume that the middle of the scan_length is when the exact second the scan is considered to be accuired.
    ;to find the middle of the scan we take 1/2 of the scan length and then add
    ; that to the scan start time.
    scan_start += scanLengthSeconds/2
    ;print, scan_start


    ;now store all the values we calculated
    project.imndArray[ci].time.array = ptr_new(scan_start)
    project.imndArray[ci].time.length = ptr_new(scanLengthSeconds)
    project.imndArray[ci].time.zero =  scan_start[n_pre] - scanLengthSeconds[n_pre]/2
    project.imndArray[ci].time.min = scan_start[0]- scanLengthSeconds[0]/2
    project.imndArray[ci].time.max = scan_start[adim-1]+scanLengthSeconds[adim-1]/2

;   print, *project.imndArray[ci].time.array
;   print, *project.imndArray[ci].time.length
;   print, project.imndArray[ci].time.zero
;   print, project.imndArray[ci].time.min
;   print, project.imndArray[ci].time.max

    heap_gc
end


; Subroutine name: extract_simple_string
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/09/30.
    ;Fix spelling mistakes and comment

FUNCTION extract_simple_string ,string_to_search ,start_string
    COMPILE_OPT IDL2

    index1 = Strpos(string_to_search, start_string)
    if index1 eq -1 then return, ''
    index1 =  index1 + STRLEN(start_string)
    index2 = Strpos(string_to_search, '##$IMND' , index1 )

    ; since the new systems have $$ in there imnd files now we have
    ; check to see which index is closer index2 or 3 now and use the
    ; shorter of the two.
    index3 = strpos(string_to_search, '$$', index1)
    if index3 gt 0 and index3 lt index2 then begin
       index2 = index3
    end

    difference =  index2 - index1

    ;print, 'index',+ index1
    ;print, 'difference: ',+ difference

    RETURN, STRMID(string_to_search, index1, difference)
END


; Subroutine name: extract_underscore
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/02.
    ;Fix spelling mistakes and comment

;
FUNCTION extract_underscore , string_to_search
    COMPILE_OPT IDL2

    RETURN , STRSPLIT (string_to_search, '_' ,/EXTRACT)
END


; Subroutine name: extract_boolean
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/02.
    ;Fix spelling mistakes and comment

FUNCTION extract_boolean , string_to_search
    COMPILE_OPT IDL2

    if STRMATCH(string_to_search, 'y*', /FOLD_CASE) then $
       return, 1

    if STRMATCH(string_to_search, 'o*', /FOLD_CASE) then $
       return, 0

    return, 0

END


; Subroutine name: extract_date
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/02.
    ;Fix spelling mistakes and comment
    ;Fixed date keyword with PV6 (Magdoom, 09/14/15)

FUNCTION extract_date, string_to_search
    COMPILE_OPT IDL2

    index1 = Strpos(string_to_search, '$$ 2')
    index2 = Strpos(string_to_search, '$$', index1 + 1)
    difference = index2 - index1
    temp = STRMID(string_to_search, index1, difference)
    temp1 = stregex(temp, '^\$\$ ([a-zA-Z]{3}) ([a-zA-Z]{3})[ ]+([0-9]+) (([0-9]{2}):([0-9]{2}):([0-9]{2})) ([0-9]{4}) (.+)$', /extract, /subexp)
    if strlen(temp1[0]) eq 0 then begin
      temp1 = stregex(temp, '^\$\$ (([0-9]{4})-([0-9]{2})-([0-9]{2})) (([0-9]{2}):([0-9]{2}):([0-9]{2}))', /extract, /subexp)
      return, temp1[2]+'-'+temp1[3]+'-'+temp1[4]+' '+temp1[5]
    endif else return, temp1[2]+' '+temp1[3]+' '+temp1[8]+' '+temp1[5]+':'+temp1[6]+':'+temp1[7]
    ;;RETURN, STRMID(temp,7, 20)
END

; Subroutine name: extract_display_name
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/02.
    ;Fix spelling mistakes and comment

 FUNCTION extract_display_name , file_Path
    COMPILE_OPT IDL2

    ;need to change for windows and unix file paths
    ;print,'1 ',file_Path
    length = STRLEN(file_Path)
    temp = STRMID(file_Path, 0 , length-1)
    ;print,'2 ',temp
    index1 = Strpos(temp, '\',/REVERSE_OFFSET, /REVERSE_SEARCH )
    temp = STRMID(temp, index1+1 , length-1-index1)
    ;print,'3 ',temp
    return, temp
 end

; =============================================================================
; ====================== Functions that extract IMND or METHOD data ===========
; =============================================================================

; Function name: mas_extract_method
; Created by: HS - 20061219
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; To read an incoming string of a METHOD file and parse through the parameters.

; Editing Information:
; Fixed bug reading variable TR experiments acquired using PV6 (Magdoom 09/10/15)

function mas_extract_method , header_string , index
	COMPILE_OPT IDL2
    COMMON scan_data

IF N_ELEMENTS(index) EQ 1 THEN $
    ii= index $
ELSE ii= project.ni


; extract the dimensions
temp = extract_simple_value( header_string ,'##$PVM_Matrix=(',')' )
if temp eq 0 then temp = extract_simple_value( header_string ,'##$PVM_SpecMatrix=(',')' )
if project.imndArray[ii].dimensions eq 0 then project.imndArray[ii].dimensions = temp

n = extract_simple_value( header_string ,'##$PVM_SPackArrSliceOrient=(',')' )
temp = extract_simple_value(header_string ,'##$PVM_SPackArrSliceOrient=( '+n+' )', '##', '$')
project.imndArray[ii].orientation[0] = temp

n = extract_simple_value( header_string ,'##$PVM_SPackArrReadOrient=(',')' )
temp = extract_simple_value(header_string ,'##$PVM_SPackArrReadOrient=( '+n+' )', '##', '$')

if (strpos(temp, '_') ne -1) then begin
    splt = extract_underscore(temp)
endif else begin
    splt = [temp, '']
endelse

project.imndArray[ii].orientation[1:2] = splt[0:1]

; Need to add the orientation for DWI files.

;============= EXTRACTING spectral frequency ===============
temp = extract_simple_value( header_string ,'##$PVM_FrqWork=( 8 )','##$', '$$')
if temp eq 0 then temp = extract_simple_value( header_string ,'##$SFO1=','#')
project.imndArray[ii].spect_bf1  = temp

;============= EXTRACTING Recovery Time ===============
temp = extract_simple_value( header_string ,'##$PVM_RepetitionTime=','##$', '$$' )
project.imndArray[ii].recov_time = temp


; ============= EXTRACTING Number of Averages ===============

temp = extract_simple_value( header_string ,'##$PVM_NAverages=','##$', '$$' )
if temp eq 0 then temp = 1
project.imndArray[ii].n_avg = temp
project.imndArray[ii].pn_avg = ptr_new(temp)


; ============= EXTRACTING Number of Slices ===============

temp = extract_simple_value( header_string ,'##$PVM_ObjOrderList=(' ,' )', '##$')
project.imndArray[ii].slices = temp


; ============= EXTRACTING echo time ===============

temp = extract_simple_value( header_string ,'##$PVM_EchoTime=','##$', '$$' )
project.imndArray[ii].echo_time = temp


; ============= EXTRACTING slice thickness ===============

temp = extract_simple_value( header_string ,'##$PVM_SliceThick=','##$', '$$' )
project.imndArray[ii].thick = temp


; ============= EXTRACTING Fat Suppression ===============

temp = extract_simple_value( header_string ,'##$PVM_FatSupOnOff=','##$', '$$' )
project.imndArray[ii].fat_sup = extract_boolean(temp)


; ============= EXTRACTING DATE AND NAME ===============

project.imndArray[ii].scan_date = extract_date (header_string)

project.imndArray[ii].scan_name = extract_name (header_string)


project.imndArray[ii].display_Name = project.imndArray[ii].file_Path


; ============= EXTRACTING Slice Scheme ===============

temp = extract_slice_scheme_methodfiles( header_string,'##$PVM_ObjOrderList' , project.imndArray[ii].slices )
temp = Ptr_New(temp)
project.imndArray[ii].slice_scheme_ptr = temp


; ============= EXTRACTING RARE Factor ===============

temp = extract_simple_value( header_string ,'PVM_RareFactor=','#' )
; If it can't find a rare value it assumes 1.
if ( temp eq 0) then temp = 1
project.imndArray[ii].rare = temp


; ============= EXTRACTING Number of Echo Images ===============

temp = extract_simple_value( header_string ,'##$PVM_NEchoImages=','##', '$$' )
project.imndArray[ii].n_echo = temp

;==================1D DATA EXTRACTION====================
IF (project.imndArray[ii].dimensions EQ 1) THEN BEGIN
  method = extract_name(header_string)
  project.imndArray[ii].scan_name = method
  project.imndArray[ii].scan_date = extract_date (header_string )
  project.imndArray[ii].scan_name = extract_name (header_string )
  project.imndArray[ii].display_Name = project.imndArray[ii].file_Path

  return, extract_bruker_onepulse(header_string)
ENDIF 

;==================2D DATA EXTRACTION====================

IF (project.imndArray[ii].dimensions EQ 2) THEN BEGIN


	; ============= extract Anti-Aliasing Factor (Over-Sampling) ============
	; This will be used later to read in the correct matrix and update the FOV


	text = extract_simple_value( header_string ,'##$PVM_AntiAlias=( 2 )','##$', '$$' )
	anti_alias = float(STRSPLIT(text, /EXTRACT))

	isoversampled = where(anti_alias gt 1)
	size_isoversampled = size(isoversampled)

  text = extract_simple_value( header_string ,'##$PVM_EncMatrix=( 2 )','##$', '$$' )
  enc_matrix_temp = STRSPLIT(text, /EXTRACT)
  temp = enc_matrix_temp
  
  ;; Not sure why you would want *not* want to use the encoding matrix to specify
  ;; dimensions, since MAS does not attempt to zero-fill to the user-requested
  ;; reconstruction dimensions specified on the console.
	if (0 and isoversampled eq -1) then begin

		text = extract_simple_value( header_string ,'##$PVM_Matrix=( 2 )','##$', '$$' )
	  matrix_temp = STRSPLIT(text, /EXTRACT)
	  
  endif

    project.imndArray[ii].fdim = temp[0]
    fdim = project.imndArray[ii].fdim
    
    project.imndArray[ii].pdim = temp[1]
    pdim = project.imndArray[ii].pdim
    
    project.imndArray[ii].sdim = project.imndArray[ii].slices
    sdim = project.imndArray[ii].sdim

    project.imndarray[ii].k_fdim_span_max = (project.imndarray[ii].k_fdim_span = fdim)
    project.imndarray[ii].k_pdim_span_max = (project.imndarray[ii].k_pdim_span = pdim)
    project.imndarray[ii].k_sdim_span_max = (project.imndarray[ii].k_sdim_span = sdim)

    ;=========== Extract Multiple Diffusion Experiments ===========
    ; I'm not quite sure this does anything productive

    temp = extract_simple_value( header_string ,'##$PVM_DwNDiffDir=','##' )
    temp1 = extract_simple_value( header_string ,'##$PVM_DwNDiffExp=','##' )

    ; if there are more than one diffusion experiments then:
    if temp1 gt 1 and temp1 gt temp then project.imndArray[ii].adim = temp1
       
    ;=========== Extract Variable TR Experiments ===========
;       temp = extract_simple_value( header_string ,'##$MultiRepetitionTime=(',')' )
    temp = extract_simple_value( header_string ,'##$MultiRepTime=(',')' )
    if temp gt 0 then project.imndArray[ii].adim = temp
   
    
    ;=========== Extract Variable TI Experiments ===========
    temp = extract_simple_value( header_string ,'##$PVM_FairTIR_Arr=(',')' )
    if temp gt 0 then project.imndArray[ii].adim = temp
   
    ;=========== Extract Repetitions ===========
    ;       temp = extract_simple_value( header_string ,'##$MultiRepetitionTime=(',')' )
    temp = extract_simple_value( header_string ,'##$PVM_NRepetitions=','##$', '$$' )
    if project.imndArray[ii].adim le 1 and temp gt 0 then project.imndArray[ii].adim = temp
       
    ;========== Extract Field of View (FOV) =========
    ; When this happens then we have a METHOD file and we need to determine if the FOV
    ; is in cm or mm. Our default will be to convert mm to cm.

    text = extract_simple_value( header_string ,'##$PVM_FovCm=( 2 )','##$', '$$' )
      ; Cannot find the FovCm, it is probably a DWI image.

        if text eq 0 then begin

	        text = extract_simple_value(header_string, '##$PVM_Fov=( 2 )', '##$', '$$')
	       ; convert to centimeters
	        temp1 = 0.1*long((STRSPLIT(text, /EXTRACT))[0])
	        temp2 = 0.1*long((STRSPLIT(text, /EXTRACT))[1])
	        text = string(temp1) + string(temp2)

		endif

	temp = STRSPLIT(text, /EXTRACT)

    project.imndArray[ii].f_fov = temp[0]*anti_alias[0]
    project.imndArray[ii].p_fov = temp[1]*anti_alias[1]
    pfov = project.imndArray[ii].p_fov

    project.imndArray[ii].s_fov = sdim*project.imndArray[ii].thick/10

    ;; compute the voxel dimensions
    vox_dim = [ project.imndArray[ii].f_fov/fdim, $
                project.imndArray[ii].p_fov/pdim, $
                project.imndArray[ii].s_fov/sdim ]

    project.imndArray[ii].f_voxsz = vox_dim[0]
    project.imndArray[ii].p_voxsz = vox_dim[1]
    project.imndArray[ii].s_voxsz = vox_dim[2]

    temp = extract_simple_value( header_string ,'##$PVM_SPackArrPhase1Offset=( 1 )','##$', '$$')

; HS - 20061220
; I will be removing the negative from the 0.1 for the next instruction.
; This calculates the phase shift and it has always been wrong before except on the 3D case.
; Guess what, the 3D case has a positive 0.1 factor. There is also a curious thing in mas_shift with putting
; a negative in front of this value when it is read in. I am going to leave it in that subroutine
; and change it here. I also changed it on the next function mas_extract_imnd.

    project.imndArray[ii].pdim_shift = round((0.1*temp/project.imndArray[ii].p_fov)* pdim)


ENDIF


;==================3D DATA EXTRACTION===================
IF (project.imndArray[ii].dimensions EQ 3) THEN BEGIN

	; ============= extract Anti-Aliasing Factor (Over-Sampling) ============
	; This will be used later to read in the correct matrix and update the FOV


	text = extract_simple_value( header_string ,'##$PVM_AntiAlias=( 3 )','##$', '$$' )
	anti_alias = float(STRSPLIT(text, /EXTRACT))

	isoversampled = where(anti_alias gt 1)
	size_isoversampled = size(isoversampled)

	if (isoversampled eq -1) then begin

		text = extract_simple_value( header_string ,'##$PVM_Matrix=( 3 )','##$', '$$' )
	    temp = STRSPLIT(text, /EXTRACT)

	endif else begin

		text = extract_simple_value( header_string ,'##$PVM_EncMatrix=( 3 )','##$', '$$' )
	    temp = STRSPLIT(text, /EXTRACT)

	end

    project.imndArray[ii].fdim = temp[0]
    fdim = project.imndArray[ii].fdim
    
    project.imndArray[ii].pdim = temp[1]
    pdim = project.imndArray[ii].pdim

    project.imndArray[ii].sdim = temp[2]
    sdim = project.imndArray[ii].sdim

    project.imndarray[ii].k_fdim_span_max = (project.imndarray[ii].k_fdim_span = fdim)
    project.imndarray[ii].k_pdim_span_max = (project.imndarray[ii].k_pdim_span = pdim)
    project.imndarray[ii].k_sdim_span_max = (project.imndarray[ii].k_sdim_span = sdim)

    project.imndArray[ii].adim = 1

    ;========== Extract FOV array=========


    text = extract_simple_value( header_string ,'##$PVM_FovCm=( 3 )','##$', '$$' )

      ; Cannot find the FovCm, try in mm.
    if text eq 0 then begin

        text = extract_simple_value(header_string, '##$PVM_Fov=( 3 )', '##$', '$$')
	       ; convert to centimeters
        temp1 = 0.1*long((STRSPLIT(text, /EXTRACT))[0])
        temp2 = 0.1*long((STRSPLIT(text, /EXTRACT))[1])
        temp3 = 0.1*long((STRSPLIT(text, /EXTRACT))[2])
        text = string(temp1) + string(temp2) + string(temp3)

    endif

    temp = STRSPLIT(text, /EXTRACT)
    project.imndArray[ii].f_fov = temp[0]*anti_alias[0]
    project.imndArray[ii].p_fov = temp[1]*anti_alias[1]
    project.imndArray[ii].s_fov = temp[2]*anti_alias[2]
    pfov = project.imndArray[ii].p_fov
    sfov = project.imndArray[ii].s_fov

    ;; compute the voxel dimensions
    vox_dim = [ project.imndArray[ii].f_fov/fdim, $
                project.imndArray[ii].p_fov/pdim, $
                project.imndArray[ii].s_fov/sdim]

    project.imndArray[ii].f_voxsz = vox_dim[0]
    project.imndArray[ii].p_voxsz = vox_dim[1]
    project.imndArray[ii].s_voxsz = vox_dim[2]

	; ====== Pshift =====
    temp = extract_simple_value( header_string ,'##$PVM_SPackArrPhase1Offset=( 1 )','##$', '$$')
    pshift = temp
    project.imndArray[ii].pdim_shift = round(  (0.1*pShift/project.imndArray[ii].p_fov) * pdim )


    temp = extract_simple_value( header_string ,'##$PVM_SPackArrSliceOffset=( 1 )','##$', '$$')
    sShift = temp
    project.imndArray[ii].sdim_shift = round( ((0.1*sShift/project.imndArray[ii].s_fov)+0.5) * sdim )


ENDIF

;========== Extract Variable TR Experiments =========
;##$MultiRepetitionTime=( 5 )
;5000.000 2000.000 1000.000 500.000 300.000
varRepTime = extract_slice_scheme( header_string ,'##$MultiRepetitionTime', project.imndArray[ii].adim )
if varRepTime[0] eq -1 then varRepTime = extract_slice_scheme( header_string ,'##$MultiRepTime', project.imndArray[ii].adim )

if varRepTime[0] ne -1 then begin
	project.imndArray[ii].rep_Time_ptr = ptr_new(varRepTime)
	project.imndArray[ii].recov_time = max(varRepTime)
endif

;========== Extract Variable TI Experiments =========
varIRTime = extract_slice_scheme( header_string ,'##$PVM_FairTIR_Arr', project.imndArray[ii].adim )
if varIRTime[0] ne -1 then project.imndArray[ii].inversion_Time_ptr = ptr_new(varIRTime)

; For Method files directions will be input in x,y,z coordinates. We need to convert them to Spherical Coordinates
;##$PVM_DwDir=( 1, 3 )
;0.000000 0.000000 1.000000

temp = extract_theta_methods(header_string)
if temp[0] ne -1 then project.imndArray[ii].angle_theta = ptr_new(float(temp))

temp = extract_phi_methods(header_string)
if temp[0] ne -1 then project.imndArray[ii].angle_phi = ptr_new(float(temp))

; Using the variable varRepTime as a placeholder
varRepTime = extract_b_matrix_method_files(header_string)
;;if varRepTime[0] ne -1 then project.imndArray[ii].bval_Array = ptr_new(varRepTime)
if varRepTime[0] ne -1 then begin 
   varRepTime[3:5, *] *= 2.0
   project.imndArray[ii].b_matrix = ptr_new(float(varRepTime))
endif

varRepTime = extract_bval_method(header_string)
if varRepTime[0] ne -1 then begin
   project.imndArray[ii].bval_array = ptr_new(float(varRepTime))
   project.imndArray[ii].n_bvals = n_elements(*project.imndarray[ii].bval_array)
   
   ;; Note I am going to double check the adim value here since I believe
   ;; that it can be computed wrongly in some cases when the DTIStandard protocol
   ;; is being used.
   project.imndArray[ii].adim = n_elements(*project.imndarray[ii].bval_array)
   
   varRepTime = extract_big_delta_method(header_string)
   if varRepTime[0] ne -1 then project.imndArray[ii].big_delta = ptr_new(float(varRepTime))

   varRepTime = extract_small_delta_method(header_string)
   if varRepTime[0] ne -1 then project.imndArray[ii].small_delta = ptr_new(float(varRepTime))
endif

;if the scan was a dynamics scan the open and extract the time data.
if project.imndArray[ii].image_type eq 1 then begin
    extract_Dynamics_time
end

;set the stop slice for volumetrics measuring.
project.procPramArray[ii].vol_stop_slice = project.imndArray[ii].sdim

; Everything was performed correctly until this point so return a 1.
return, 1

end


;; Function: show_IMND_file
;; Created by: HS 20070612
;; Calling Information: Called by MAS.PRO
;
;Pro show_IMND_file
;	COMPILE_OPT IDL2
;    COMMON scan_data
;	COMMON common_widgets
;    CI = project.CI
;    display_title = ''
;	slash = get_dir_slash()
;	adim = project.imndArray[CI].adim
;
;	if (adim gt 1) then begin
;
;		multiple_file_array = *project.imndArray[CI].multi_scan_file_array
;		display_title = multiple_file_array[0] + slash[1]
;
;	endif else begin
;
;		display_title = project.imndArray[CI].display_Name
;	end
;
;
;
;	; Lets make the window size dependent on the screen that it is being shown.
;	dimensions = GET_SCREEN_SIZE()
;	xsize = dimensions[0]/1.5
;	ysize = dimensions[1]-dimensions[1]*.1
;
;	IMND = project.imndArray[CI].imnd_file
;	IMND_filepath = display_title+'imnd'
;
;
;	IMND_Exists = FILE_SEARCH(IMND_filepath)
;	IMND_size = size(IMND_Exists)
;
;	if IMND_size[0] eq 0 then begin
;
;		IMND_filepath = IMND_filepath + 'method'
;		METHOD_Exists = FILE_SEARCH(IMND_filepath)
;		METHOD_size = size(METHOD_Exists)
;
;			if METHOD_size[0] eq 0 then begin
;				update_status_bar,'The directory does not contain an imnd or method files'
;	       		return
;			endif
;
;	endif
;
;	; If we are here it's because IMND or METHOD file exists
;
;
;    parameter_temp_string = ''
;
;    OPENR, IMND_file_pointer, IMND_filepath , /GET_LUN , ERROR = err
;
;	IMND_WINDOW_BASE = widget_base(TITLE=IMND_filepath , UVALUE = IMND_WINDOW_BASE, XSIZE = xsize, YSIZE = ysize)
;
;
;	; Read in the IMND line by line and store it in an array
;	number_lines = FILE_LINES( 	IMND_filepath )
;	IMND_file_array = strarr(number_lines)
;	counter = 0
;
;	While (eof(IMND_file_pointer) NE 1) Do Begin
;	       Readf, IMND_file_pointer, parameter_temp_string
;	       IMND_file_array[counter] = parameter_temp_string
;	       counter = counter + 1
;	End
;
;	IMND_textbox = widget_text(IMND_WINDOW_BASE, VALUE = IMND_file_array, /SCROLL, scr_xsize = xsize, scr_ysize = ysize)
;	;widget_control, IMND_textbox, SET_VALUE = IMND_file_array
;
;
;	;widget_control, IMND_textbox, SET_VALUE = parameter_file_string
;	widget_control, IMND_WINDOW_BASE, /realize
;
;End


; Function: show_text_file
; Created by: HS 20070802
; Calling Information: Called by MAS.PRO
; Extracts the appropriate (IMND, Method or ACQP) text file to show in a new window.

Pro show_text_file, type

COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets
    CI = project.CI
    display_title = ''
    slash = get_dir_slash()
    adim = project.imndArray[CI].adim

    if (ptr_valid(project.imndArray[CI].multi_scan_file_array) and project.imndArray[CI].DICOMtype EQ 0) then begin

;;    if (adim gt 1) AND (project.imndArray[CI].DICOMtype EQ 0) then begin

        multiple_file_array = *project.imndArray[CI].multi_scan_file_array
        display_title = multiple_file_array[0] + slash[1]

    endif else begin

        display_title = project.imndArray[CI].display_Name
    end

	; Lets make the window size dependent on the screen that it is being shown.
	dimensions = GET_SCREEN_SIZE()
	xsize = dimensions[0]/1.5
	ysize = dimensions[1]-dimensions[1]*.1


if type eq 'imnd' then begin


	IMND = project.imndArray[CI].imnd_file
	textfile_path = display_title+'imnd'


	IMND_Exists = FILE_SEARCH(textfile_path)
	IMND_size = size(IMND_Exists)

	if IMND_size[0] eq 0 then begin

		textfile_path = display_title + 'method'
		METHOD_Exists = FILE_SEARCH(textfile_path)
		METHOD_size = size(METHOD_Exists)

			if METHOD_size[0] eq 0 then begin
				update_status_bar,'The directory does not contain an imnd or method files'
	       		return
			endif

	endif

	; If we are here it's because IMND or METHOD file exists

endif else if type eq 'acqp' then begin

	textfile_path = display_title + 'acqp'

	ACQP_Exists = FILE_SEARCH(textfile_path)
	ACQP_size = size(ACQP_Exists)

	if ACQP_size[0] eq 0 then begin
	  textfile_path = display_title + 'acqu'
	  ACQU_Exists = FILE_SEARCH(textfile_path)
	  ACQU_size = size(ACQU_Exists)
	  if ACQU_size[0] eq 0 then begin 
		  update_status_bar,'The directory does not contain an acqp file'
      return
    endif
	endif

	endif else if type eq 'config' then begin

	  textfile_path = display_title + 'configscan'

	  CONFIG_Exists = FILE_SEARCH(textfile_path)
	  CONFIG_size = size(CONFIG_Exists)
	  if CONFIG_size[0] eq 0 then begin
	    update_status_bar,'The directory does not contain a config file'
	    return
	  endif
	  
endif else if type eq 'PAR' then begin
	textfile_path = project.imndArray[CI].file_Path + '.PAR'

	PAR_Exists = FILE_SEARCH(textfile_path)
	PAR_size = size(PAR_Exists)

	if PAR_size[0] eq 0 then begin
		update_status_bar,'The directory does not contain a PAR file'
		return
	endif

endif else if type eq 'procpar' then begin

    textfile_path = project.imndArray[CI].file_Path
    if (file_test(textfile_path, /read) eq 0) then begin
        update_status_bar, "PROCPAR file not found."
        return
    endif
    
endif else if type eq 'procpar-pretty' then begin
    opp = obj_new('mas_varian_procpar', project.imndarray[CI].file_path)
    opp->DumpParameters
    return

endif

    parameter_temp_string = ''

    OPENR, text_file_pointer, textfile_path , /GET_LUN , ERROR = err

	TEXT_WINDOW_BASE = widget_base(TITLE=textfile_path , UVALUE = TEXT_WINDOW_BASE, XSIZE = xsize, YSIZE = ysize)


	; Read in the IMND line by line and store it in an array
	number_lines = FILE_LINES( 	textfile_path )
	TEXT_file_array = strarr(number_lines)
	counter = 0

	While (eof(text_file_pointer) NE 1) Do Begin
	       Readf, text_file_pointer, parameter_temp_string
	       TEXT_file_array[counter] = parameter_temp_string
	       counter = counter + 1
	End

	TEXT_file_textbox = widget_text(TEXT_WINDOW_BASE, VALUE = TEXT_file_array, /SCROLL, scr_xsize = xsize, scr_ysize = ysize)
	;widget_control, IMND_textbox, SET_VALUE = IMND_file_array


	;widget_control, IMND_textbox, SET_VALUE = parameter_file_string
	widget_control, TEXT_WINDOW_BASE, /realize

END

; Function name: mas_extract_imnd
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; this will take in a string and extract the important data.
; after the data is extracted it will then be saved to the variable.
; index the scan in question to be updated, if index isn't specified it
; will assume that the next scan is the one to update.
; returns 1 if file is good to open 0 if file will be unable to process

; Editing Information:
    ;Edited by HS 2006/10/02.
    ;Fix spelling mistakes and comment

function mas_extract_imnd , header_string , index
    COMPILE_OPT IDL2
    COMMON scan_data

;this is me testing some thing just to see if breaking up the string on ## is a good idea
;temp = STRSPLIT(header_string, '##', /EXTRACT)
;for cc=0 , (size(temp))[1] -1 do $
;   print, temp[cc]


;==================[easy data integers and floats]===========================
IF N_ELEMENTS(index) EQ 1 THEN $
    ii= index $
ELSE ii= project.ni

;; Please see: mas_open.pro, mas_open_onepulse_spect
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; determine if this is onepulse spectroscopy
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
method = extract_name(header_string)
if (method eq 'ONEPULSE') then begin
    project.imndArray[ii].scan_name = method
    project.imndArray[ii].scan_date = extract_date (header_string )
    project.imndArray[ii].scan_name = extract_name (header_string )
    project.imndArray[ii].display_Name = project.imndArray[ii].file_Path
    
    return, extract_bruker_onepulse(header_string)
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; extract the dimensions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
temp = extract_simple_value( header_string ,'##$IMND_dim=', '##$', '$$' )
if temp eq 0 then temp = extract_simple_value( header_string ,'##$IMND_matrix=(',')' )
if temp eq 0 then temp = extract_simple_value( header_string ,'##$PVM_Matrix=(',')' )
project.imndArray[ii].dimensions = temp
;print, 'dimensions',+ project.imndArray[ii].dimensions

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;extract the slice orientation
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; HS 20061014. New scheme to extract the slice orientation for METHOD files.
temp = extract_simple_value(header_string ,'##$IMND_slice_orient=','##$', '$$' )

; HS - if the previous statement finds a string, then it should have a length greater than 1 and
; this next "if" statement does not happen. This returned string will then be split by the
; "endif" statement afterwards. If nothing was found a "0" is returned from the "extract_simple_ value"
; and this statement happens. It searches for two specific lines in the METHODS file and places
; their values in the corresponding imndArray entry.
if (strlen(temp) eq 1) then begin
	temp = extract_simple_value(header_string ,'##$PVM_SPackArrSliceOrient=( 1 )', '##', '$')
	project.imndArray[ii].orientation[0] = temp
	temp = extract_simple_value(header_string ,'##$PVM_SPackArrReadOrient=( 1 )', '##', '$')
	temp = extract_underscore(temp)
	project.imndArray[ii].orientation[1:2] = temp[0:1]

endif else begin
	temp = extract_underscore(temp)
	project.imndArray[ii].orientation = temp
endelse

; Need to add the orientation for DWI files.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;extract the recovery time
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
temp = extract_simple_value( header_string ,'##$IMND_recov_time=','##$', '$$' )
if temp eq 0 then temp = extract_simple_value( header_string ,'##$IMND_diff_rep_time=( 1 )','##$', '$$' )
if temp eq 0 then temp = extract_simple_value( header_string ,'##$PVM_RepetitionTime=','##$', '$$' )
project.imndArray[ii].recov_time = temp
;print, 'recovtime',+ project.imndArray[ii].recov_time

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; ============= EXTRACTING Number of Averages ===============
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
temp = extract_simple_value( header_string ,'##$IMND_n_averages=','##$', '$$' )
if temp eq 0 then temp = extract_simple_value( header_string ,'##$PVM_NAverages=','##$', '$$' )
project.imndArray[ii].n_avg = temp
project.imndArray[ii].pn_avg = ptr_new(temp)
;print, 'navg',+ project.imndArray[ii].n_avg

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; ============= EXTRACTING Number of Slices ===============
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
temp = extract_simple_value( header_string ,'##$IMND_n_slices=' ,'##$', '$$')

; HS 20061013 - Added this extra comparison to find the number of slices in METHOD files.
if temp eq 0 then temp = extract_simple_value( header_string ,'##$PVM_ObjOrderList=(' ,' )', '##$')
project.imndArray[ii].slices = temp
;print, 'slices',+ project.imndArray[ii].slices

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; ============= EXTRACTING echo time ===============
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
temp = extract_simple_value( header_string ,'##$IMND_echo_time=','##$','$$' )
if temp eq 0 then temp = extract_simple_value( header_string ,'##$IMND_diff_echo_time=','##$', '$$' )
if temp eq 0 then temp = extract_simple_value( header_string ,'##$PVM_EchoTime=','##$', '$$' )
project.imndArray[ii].echo_time = temp
;print, 'echotime',+ project.imndArray[ii].echo_time

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; ============= EXTRACTING slice thickness ===============
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
temp = extract_simple_value( header_string ,'##$IMND_slice_thick=','##$' , '$$')
;nr if temp eq 0 then temp = extract_simple_value( header_string ,'##$PVM_SliceThick=','##$', '$$' )
if temp eq 0.0 then temp = extract_simple_value( header_string ,'##$PVM_SliceThick=','##$', '$$' )
project.imndArray[ii].thick = temp
;print, 'thick',+ project.imndArray[ii].thick

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; ============= EXTRACTING Fat Suppression ===============
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;nr temp = extract_simple_string (header_string ,'##$IMND_suppression=')
;nr
;nr if strlen(temp) le 0 then temp = extract_simple_value( header_string ,'##$PVM_FatSupOnOff=','##$', '$$' )
temp = extract_simple_value( header_string ,'##$IMND_fat_mode=','##$', '$$' )
; HS 20061024
; Added for METHOD file extraction of fat suppresion.

if temp eq '0' then begin
temp = extract_simple_value( header_string ,'##$PVM_FatSupOnOff=','##$', '$$' )
endif

project.imndArray[ii].fat_sup = extract_boolean(temp)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; ============= EXTRACTING DATE AND NAME ===============
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
project.imndArray[ii].scan_date = extract_date(header_string )
project.imndArray[ii].scan_name = extract_name(header_string )
project.imndArray[ii].display_Name = project.imndArray[ii].file_Path

; HS 20061013 - Adding a way to read slice scheme for Method files.
; this is the old way it was done (only IMND)

; temp = Ptr_New(extract_slice_scheme( header_string,'##$IMND_slice_list' , project.imndArray[ii].slices ) )
; project.imndArray[ii].slice_scheme_ptr = temp

; The new way calls the extract_slice_scheme and if it finds a -1 (string not found) on the result
; it calls the new version extract_slice_scheme_methodfiles which is specific for METHOD files.
; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; extracting slice list and order
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
temp = extract_slice_scheme( header_string,'##$IMND_slice_list' , project.imndArray[ii].slices )
if temp[0] eq -1 then temp = extract_slice_scheme_methodfiles( header_string,'##$PVM_ObjOrderList' , project.imndArray[ii].slices )

temp = Ptr_New(temp)
project.imndArray[ii].slice_scheme_ptr = temp
;print, 'extract imnd header slice schem', *temp

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; to extract the RARE factor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
temp = FLOAT(extract_simple_value( header_string ,'IMND_rare_factor=','$' ))
if temp eq 0.0 then temp = FLOAT(extract_simple_value( header_string ,'IMND_diff_rare_factor=','$' ))
if temp LT 1 then temp = FLOAT(extract_simple_value( header_string ,'IMND_diff_rare_factor=','$' ))
; To access the rare encoding in Method files.
if temp eq 0 then temp = extract_simple_value( header_string ,'PVM_RareFactor=','#' )

; If it can't find a rare value it assumes 1.
if ( temp eq 0) then temp = 1
project.imndArray[ii].rare = temp

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; number of echo images
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
temp = extract_simple_value( header_string ,'##$IMND_n_echo_images=','$$' )
if temp eq 0 then temp = extract_simple_value( header_string ,'##$PVM_NEchoImages=','##', '$$' )
project.imndArray[ii].n_echo = temp


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================2D DATA EXTRACTION====================
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
IF (project.imndArray[ii].dimensions EQ 2) THEN BEGIN

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;===========extract dim array===========
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    text = extract_simple_string (header_string ,'##$IMND_matrix=( 2 )')
    if strlen(text) le 0 then text = extract_simple_value( header_string ,'##$PVM_Matrix=( 2 )','##$', '$$' )

    project.imndArray[ii].fdim = (STRSPLIT(text, /EXTRACT))[0]
    fdim = project.imndarray[ii].fdim
    
    project.imndArray[ii].pdim = (STRSPLIT(text, /EXTRACT))[1]
    pdim = project.imndArray[ii].pdim
    
    project.imndArray[ii].sdim = project.imndArray[ii].slices
    sdim = project.imndArray[ii].sdim

    project.imndarray[ii].k_fdim_span_max = (project.imndarray[ii].k_fdim_span = fdim)
    project.imndarray[ii].k_pdim_span_max = (project.imndarray[ii].k_pdim_span = pdim)
    project.imndarray[ii].k_sdim_span_max = (project.imndarray[ii].k_sdim_span = sdim)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;===========extract var tr===========
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    temp = extract_simple_value( header_string ,'##$IMND_n_diff_expts=','##$' )
    if temp LT 1 then temp = extract_simple_value( header_string ,'##$IMND_n_diff_expts=','$$' )

    ; HS - 20061024
    ; Adding support for METHOD files
    if temp LT 1 then temp = extract_simple_value( header_string ,'##$PVM_DwNDiffDir=','$$' )

    ; if there are more than one diffusion experiments then:
    if (temp gt 1) then begin
        project.imndArray[ii].adim = temp
    endif else begin
       temp = extract_simple_value( header_string ,'##$IMND_var_rep_time=(',')' )
       if temp eq 0 then temp = extract_simple_value( header_string ,'##$MultiRepetitionTime=(',')' )

       if temp gt 0 then begin
           project.imndArray[ii].adim = temp
       endif else begin
           project.imndArray[ii].adim = 1
       endelse
    endelse

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;==========extract slice/read vector==============
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    read_vector = bytarr(3)
    phase_vector = bytarr(3)
    slice_vector = bytarr(3)

    temp = extract_simple_value(header_string ,'##$IMND_n_slicepacks=','##$')
    n_slicepacks = strcompress(temp, /remove_all)
    nvec_qualifier = '( '+n_slicepacks+', 3 )'
    text = extract_simple_string (header_string ,'##$IMND_read_vector='+nvec_qualifier)
    if strlen(text) gt 0 then begin ;; We have an imnd file
        
        read_vector[0] = (STRSPLIT(text, /EXTRACT))[0]
        read_vector[1] = (STRSPLIT(text, /EXTRACT))[1]
        read_vector[2] = (STRSPLIT(text, /EXTRACT))[2]

        text = extract_simple_string (header_string ,'##$IMND_slicepack_vector='+nvec_qualifier)
        if strlen(text) gt 0 then begin ;; If we found the first, we should find the 2nd
            slice_vector[0] = (STRSPLIT(text, /EXTRACT))[0]
            slice_vector[1] = (STRSPLIT(text, /EXTRACT))[1]
            slice_vector[2] = (STRSPLIT(text, /EXTRACT))[2]
            
            tmp = where(read_vector + slice_vector eq 0)
            phase_vector[tmp] = 1
            
        endif else begin
            print, "Unable to extract the gradient orientation vectors."
        endelse

    endif else begin ;; We have a 'method' file
        
        text = extract_simple_value( header_string ,'##$PVM_SPackArrGradOrient=( 1, 3, 3 )','##$', '$$', /long_value )
        if (strlen(text) gt 0 and text ne '0') then begin
            read_vector[0] = (STRSPLIT(text, /EXTRACT))[0]
            read_vector[1] = (STRSPLIT(text, /EXTRACT))[1]
            read_vector[2] = (STRSPLIT(text, /EXTRACT))[2]
            phase_vector[0] = (STRSPLIT(text, /EXTRACT))[3]
            phase_vector[1] = (STRSPLIT(text, /EXTRACT))[4]
            phase_vector[2] = (STRSPLIT(text, /EXTRACT))[5]
            slice_vector[0] = (STRSPLIT(text, /EXTRACT))[6]
            slice_vector[1] = (STRSPLIT(text, /EXTRACT))[7]
            slice_vector[2] = (STRSPLIT(text, /EXTRACT))[8]
        endif

    endelse

    matrix = double([ [read_vector,0], [phase_vector,0], [slice_vector,0], [0,0,0,1]])

    if (abs(determ(matrix, /check)) eq 0) then begin
        ;; something is wrong with the matrix, it should not be singular
        print, "Acquisition matrix appears to be singular."
        matrix = diag_matrix(replicate(1.0, 4))
    endif

    project.imndarray[ii].acq_matrix = ptr_new(matrix, /no_copy)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;==========extract fov array=========
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    text = extract_simple_string (header_string ,'##$IMND_fov=( 2 )')

    ; HS - 20061013. Added this if statement to be able to read FOV information from Method files.
    ; If the previous statement is not found then an empty array is returned.
    ; When this happens then we have a METHOD file and we need to determine if the FOV
    ; is in cm or mm. Our default will be to convert mm to cm.

    ; temp = strlen(text)

    if strlen(text) le 0 then begin
      text = extract_simple_value( header_string ,'##$PVM_FovCm=( 2 )','##$', '$$' )
      ; Cannot find the FovCm, it is probably a DWI image.
        if text eq 0 then begin
        text = extract_simple_value(header_string, '##$PVM_Fov=( 2 )', '##$', '$$')
       ; convert to centimeters
        temp1 = 0.1*long((STRSPLIT(text, /EXTRACT))[0])
        temp2 = 0.1*long((STRSPLIT(text, /EXTRACT))[1])
        text = string(temp1) + string(temp2)

        endif
    endif

    project.imndArray[ii].f_fov = (STRSPLIT(text, /EXTRACT))[0]
    project.imndArray[ii].p_fov = (STRSPLIT(text, /EXTRACT))[1]
    pfov = project.imndArray[ii].p_fov

    project.imndArray[ii].s_fov = sdim*project.imndArray[ii].thick/10

    ;; compute the voxel dimensions
    vox_dim = [ project.imndArray[ii].f_fov/fdim, project.imndArray[ii].p_fov/pdim, project.imndArray[ii].s_fov/sdim]
    ;;vox_dim /= min(vox_dim)
    project.imndArray[ii].f_voxsz = vox_dim[0]
    project.imndArray[ii].p_voxsz = vox_dim[1]
    project.imndArray[ii].s_voxsz = vox_dim[2]

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; extract phase shift
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    temp = extract_simple_value( header_string ,'##$IMND_phase1_offset=( 1 )','##$' )
    if temp eq 0 then temp = extract_simple_value( header_string ,'##$IMND_phase1_offset=( 1 )','$$' )
    if temp eq 0 then temp = extract_simple_value( header_string ,'##$PVM_SPackArrPhase1Offset=( 1 )','##$', '##$')
    if temp eq 0 then temp = extract_simple_value( header_string ,'##$PVM_SPackArrPhase1Offset=( 1 )','##$', '$$')

    ;; HS - 20061220
    ;; It is interesting that there is a negative in front of the value for the phase shift
    ;; and one more in the mas_shift inside of mas_load_state_1.pro.
    ;; If you notice on the 3D data extraction part on this subroutine the 0.1 is positive and not negative
    ;; like in this part. Therefore, I will remove it here.
	;project.imndArray[ii].pdim_shift = round((-0.1*temp/project.imndArray[ii].p_fov)* pdim)
    project.imndArray[ii].pdim_shift = round((0.1*temp/project.imndArray[ii].p_fov)* pdim)

    ;; tell mas host ot load the state 1 data
    project.imndArray[ii].state1_load_procedure = 'mas_load_state_1_bruker'

ENDIF

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;==================3D DATA EXTRACTION===================
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
IF (project.imndArray[ii].dimensions EQ 3) THEN BEGIN

    ;project.imndArray[ii].necho = extract_simple_value( header_string ,'##$IMND_n_echoes=','##$' )
    ;print, 'necho',+ project.imndArray[ii].necho
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;==========extract dim array=========
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    text = extract_simple_string (header_string ,'##$IMND_matrix=( 3 )')
    if strlen(text) le 0 then text = extract_simple_value( header_string ,'##$PVM_Matrix=( 3 )','##$', '$$' )
    ;print, 'dim matrix ', text

    project.imndArray[ii].fdim = (STRSPLIT(text, /EXTRACT))[0]
    fdim = project.imndarray[ii].fdim
    ;print, 'fdim ', project.imndArray[ii].fdim

    project.imndArray[ii].pdim = (STRSPLIT(text, /EXTRACT))[1]
    pdim = project.imndArray[ii].pdim
    ;print, 'pdim ', pdim

    project.imndArray[ii].sdim = (STRSPLIT(text, /EXTRACT))[2]
    sdim = project.imndArray[ii].sdim
    ;print, 'sdim ', sdim

    project.imndarray[ii].k_fdim_span_max = (project.imndarray[ii].k_fdim_span = fdim)
    project.imndarray[ii].k_pdim_span_max = (project.imndarray[ii].k_pdim_span = pdim)
    project.imndarray[ii].k_sdim_span_max = (project.imndarray[ii].k_sdim_span = sdim)

    project.imndArray[ii].adim = 1
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;==========extract fov array=========
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    text = extract_simple_string (header_string ,'##$IMND_fov=( 3 )')
    if strlen(text) le 0 then text = extract_simple_value( header_string ,'##$PVM_Fov=( 3 )','##$', '$$' )
    ;print, 'FOV matrix ', text

    project.imndArray[ii].f_fov = (STRSPLIT(text, /EXTRACT))[0]
    ;print, 'f_fov ', project.imndArray[ii].f_fov

    project.imndArray[ii].p_fov =(STRSPLIT(text, /EXTRACT))[1]
    pfov = project.imndArray[ii].p_fov
    ;print, 'pfov ', pfov

    project.imndArray[ii].s_fov = (STRSPLIT(text, /EXTRACT))[2]
    sfov = project.imndArray[ii].s_fov
    ;print, 'sfov ', sfov

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; compute the voxel dimensions
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    vox_dim = [ project.imndArray[ii].f_fov/fdim, project.imndArray[ii].p_fov/pdim, project.imndArray[ii].s_fov/sdim]
    ;;vox_dim /= min(vox_dim)
    project.imndArray[ii].f_voxsz = vox_dim[0]
    project.imndArray[ii].p_voxsz = vox_dim[1]
    project.imndArray[ii].s_voxsz = vox_dim[2]

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; extract phase shift
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    temp = extract_simple_value( header_string ,'##$IMND_phase1_offset=( 1 )','##$', '$$')
    if temp eq 0 then temp = extract_simple_value( header_string ,'##$PVM_SPackArrPhase1Offset=( 1 )','##$', '$$')
    pshift = temp
    project.imndArray[ii].pdim_shift = round(  (0.1*pShift/project.imndArray[ii].p_fov) * pdim )
    ;print, 'imnd pdim_shift= ', pShift

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;; extract slice (phase2 shift) shift
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    temp = extract_simple_value( header_string ,'##$IMND_slicepack_position=( 1 )','##$' )
    if temp eq 0 then temp = extract_simple_value( header_string ,'##$PVM_SPackArrSliceOffset=( 1 )','##$', '$$')
    sShift = temp
    project.imndArray[ii].sdim_shift = round( ((0.1*sShift/project.imndArray[ii].s_fov)+0.5) * sdim )
    ;print, 'imnd sdim_shift= ', sShift

ENDIF

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;========== Extract Variable TR Experiments =========
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;Must change for METHODS file

;##$IMND_var_rep_time=( 5 )
;6000.00 2000.00 1000.00 500.00 274.18

;##$MultiRepetitionTime=( 5 )
;5000.000 2000.000 1000.000 500.000 300.000

varRepTime = extract_slice_scheme ( header_string, '##$IMND_var_rep_time' ,project.imndArray[ii].adim )
if varRepTime[0] eq -1 then varRepTime = extract_slice_scheme( header_string ,'##$MultiRepetitionTime', project.imndArray[ii].adim )

if varRepTime[0] ne -1 then project.imndArray[ii].rep_Time_ptr = ptr_new(varRepTime)
;print, *project.imndArray[ii].repTime_ptr

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Extract diffusion angles
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; For IMND files:
;##$IMND_diff_angle_theta=( 1 )
;63.44

; For Method files:
;##$PVM_DwDir=( 1, 3 )
;0.000000 0.000000 1.000000

temp = extract_slice_scheme ( header_string, 'diff_angle_theta' ,project.imndArray[ii].adim )
if temp[0] eq -1 then temp = extract_theta_methods(header_string)
if temp[0] ne -1 then project.imndArray[ii].angle_theta = ptr_new(temp)

;##$IMND_diff_angle_phi=( 1 )
;288.00
temp = extract_slice_scheme ( header_string, 'diff_angle_phi' ,project.imndArray[ii].adim )
if temp[0] eq -1 then temp = extract_phi_methods(header_string)
if temp[0] ne -1 then project.imndArray[ii].angle_phi = ptr_new(temp)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; extract diffusion bvalue and b matrix
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
varRepTime = extract_slice_scheme ( header_string, '##$IMND_diff_b_value' ,project.imndArray[ii].adim )
if varRepTime[0] eq -1 then varRepTime = extract_b_matrix_method_files(header_string)
;;if varRepTime[0] eq -1 then varRepTime = extract_bval_method(header_string)
if varRepTime[0] ne -1 then project.imndArray[ii].bval_Array = ptr_new(varRepTime)
;print, 'bval_Array', *project.imndArray[ii].bval_Array

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; extract dynamics information
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;if the scan was a dynamics scan the open and extract the time data.
if project.imndArray[ii].image_type eq 1 then begin
    extract_Dynamics_time
end

;set the stop slice for volumetrics measuring.
project.procPramArray[ii].vol_stop_slice = project.imndArray[ii].sdim

return, 1

end
