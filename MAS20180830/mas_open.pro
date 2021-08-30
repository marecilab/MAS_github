;; $Id$
;;Copyright 2003 University of Florida. All Rights Reserved

; List of subroutines contained in this file:


; Subroutine name: mas_chk_bruker_version_3
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;   HS 2006/09/29. This function does NOT swaps endians at this point.

; Purpose of subroutine:
;   This function will check what machine the image was collected on.
;   If the machine is not an SGI then swap the endian for now. It may
;   be the other way around sometime soon.

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting.

function mas_chk_bruker_version_3, imnd_file
    compile_opt idl2
    common scan_data


    ver_pos_start = strpos(imnd_file, 'VERSION')

    ; Check to see if the version keyword is present.
    ; If not then return and dont swap the endian.

    if ver_pos_start eq -1 then return, 0 $
    else return, 1

;   ver_pos_end = strpos(imnd_file, '>',ver_pos_start )
;
;   length = ver_pos_end- ver_pos_start
;
;   raw_version = STRSPLIT(strmid(imnd_file, ver_pos_start, length ), /EXTRACT )
;   print, raw_version[1]

;   void = dialog_message('Is this coming from the new Linux box'
;
     ;project.dataArray[CI].state1 =  ptr_new(SWAP_ENDIAN(  *project.dataArray[CI].state1 ) )

end


; Subroutine name: extract_file_pointer
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;

; Purpose of subroutine:
;
;

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting.


function extract_file_pointer , imnd_file_pointer
    COMPILE_OPT IDL2

    difference = STRLEN(imnd_file_pointer) - 4

    return, STRMID(imnd_file_pointer, 0, difference)
end

; Subroutine name: get_dir_slash
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;

; Purpose of subroutine:
; This function will return the proper slashes for directories depending upon
; which operating system MAS is running on.
; slash[0] is the slash that goes in front of a directory
;           this is useful in unix b/c every thing starts at /
;           in windows this is blank.
; slash[1] is the slash to show the difference between directories
;           in unix the slash is / in windows =  \

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting.

function get_dir_slash
    slash = strarr(2)
    if !d.name eq 'X' then begin
        slash[0] = '/'
        slash[1] = '/'
    endif else begin
        slash[0] = ''
        slash[1] = '\'
    end

    return, slash

end


; Subroutine name: open_imnd_file
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;

; Purpose of subroutine:

; This function will open an IMND file in the directory pointer pointed to by
; imnd_file_pointer. If this file does not exist then it will return -1.

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting.

function open_imnd_file, imnd_file_pointer
    COMPILE_OPT IDL2


    ;print,'imnd_file_pointer:',imnd_file_pointer

    imnd_file_pointer = imnd_file_pointer + 'imnd'

    ;print,'imnd_file_pointer:',imnd_file_pointer

    OPENR, imnd_file_descriptor, imnd_file_pointer , /GET_LUN , ERROR = err

    ; If err is nonzero, something happened. Print the error message to
    ; the standard error file (logical unit -2):
    IF (err NE 0) then  begin
       ;PRINTF, -2, !ERROR_STATE.MSG
       update_status_bar,'Please open the directory containing the imnd files'

       return,-1
    end

    ;print,'Able to open the imnd file'

    imnd_file_string = ''
    imnd_temp_string = ''

    While (eof(imnd_file_descriptor) NE 1) Do Begin
       Readf, imnd_file_descriptor, imnd_temp_string
        imnd_file_string = imnd_file_string + imnd_temp_string
    EndWhile

    ;print,'imnd_file_string is: ',imnd_file_string

    close,imnd_file_descriptor
    FREE_LUN, imnd_file_descriptor

    ;print, imnd_file_string
    return, imnd_file_string
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Subroutine name: mas_open_ADC_Multi
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;

; Purpose of subroutine:
; HS 2006/09/29. Opens a series of T2 weighted images.


; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting

pro mas_open_ADC_Multi
    COMPILE_OPT IDL2
    COMMON scan_data

    ;Get the proper slash depending upon which operation system is in use.
    slash = get_dir_slash()

    ;Check to make sure that max_num_scans was not exceeded.
    if project.ni eq project.max_num_scans then begin
       update_status_bar, 'Reached max limit of scans able to open.'
       return
    end

    ;Have the user select a directory that they want to try and open that scan.
    parent_directory = DIALOG_PICKFILE( Title = 'Select parent directory of the multiple ADC scans' ,$
                 path=project.current_path,DIRECTORY=1)

    ;If the directory returned is empty, then they closed the open without a selection then just return.
    if (parent_directory eq '') then return

    ;Get a listing of the files and directories that are in the parent directory.
    parent_file_listing = FILE_SEARCH(parent_directory+'*',  COUN=n_parent_dir)
    file_adim_str = ''

    ;Next we need to extract the echo time so that the we can change the order of the multi_slice_dir
    ;to the proper read in order of fids
    ;First we will open each scan individually and extract the echo time.
    ;Then we will sort the echo time and that will tell us how to read in the T2 fid files
    ;for each scan.

    bvals = fltarr(n_parent_dir)
    n_avg = bytarr(n_parent_dir)
    b_matrix = fltarr(6, n_parent_dir)
    b_matrix_string = strarr(6)
    imnd_check = 1

    for ii=0, n_parent_dir-1 do begin

        imnd_file_string = ''
        imnd_temp_string = ''
        OPENR, imnd_file_descriptor, parent_file_listing[ii]+slash[1]+'imnd' , /GET_LUN , ERROR = err
        
        IF (err NE 0) then  begin
           imnd_check = 0
           open_file  = parent_file_listing[ii]+slash[1] + '\method'
           OPENR, imnd_file_descriptor, open_file, /GET_LUN , ERROR = err
           
           IF (err NE 0) then begin
              ;; If err is nonzero, something happened. Print the error message to
              ;; the standard error file (logical unit -2):
              ;;PRINTF, -2, !ERROR_STATE.MSG
              update_status_bar,'Please open the directory containing imnd or method files'
              return
           end
        end

        ;;if there was no error then the file's text can be extracted.
        While (eof(imnd_file_descriptor) NE 1) Do Begin
           Readf, imnd_file_descriptor, imnd_temp_string
           imnd_file_string = imnd_file_string + imnd_temp_string
        EndWhile
        close, imnd_file_descriptor
        FREE_LUN, imnd_file_descriptor
        
        if (imnd_check eq 1) then begin
           fa_str = extract_simple_value( imnd_file_string ,'##$IMND_n_diff_expts=','$' )
           n_avg_str= extract_simple_value( imnd_file_string ,'##$IMND_n_averages=','$' )
           b_str = extract_bval_imnd(imnd_file_string)
           temp = extract_simple_value( imnd_file_string ,'##$IMND_n_diff_expts=','$')
           b_matrix_temp = extract_b_matrix( imnd_file_string, temp)
        endif else begin
           fa_str = extract_simple_value( imnd_file_string ,'##$PVM_DwNDiffExp=','#' )
           n_avg_str = extract_simple_value( imnd_file_string ,'##$PVM_NAverages=','#' )
           b_str = extract_bval_method(imnd_file_string)
           b_matrix_temp = extract_b_matrix_method_files(imnd_file_string)
        endelse
        
        n_avg[ii] = n_avg_str
        bvals[ii] = b_str
        file_adim_str += string(fa_str)+' '

        for jj=0, 5 do begin
           b_matrix_string[jj] += b_matrix_temp[jj]
        endfor

    endfor

    for jj=0, 5 do begin
       btmp =  rotate(float(strsplit(b_matrix_string[jj], /extract)),1)
       b_matrix[jj,*] = btmp[where(btmp ne -1)]
    endfor
    
    ;;multiply the off diagionals by 2
    b_matrix[3:5,*] = b_matrix[3:5,*] * 2

    CLOSE, /all

    sorted = sort(bvals)
    ;;sort the echo time and use that as the key to the parent_file_listing
    project.imndArray[project.ni].multi_scan_file_array = ptr_new(parent_file_listing[sorted])
    project.imndArray[project.ni].file_adim = ptr_new(fix(strsplit(file_adim_str, /extract)))
    ;store the echo times to the global structure

    ;next we need to save the reciently extracted string to the project structure
    project.imndArray[project.ni].imnd_File = imnd_file_string
    
    ;once the imnd string has been extracted then store where we found it at
    ;then update the path so that it contains the last file open
    ;that way when the next file is open it will be close to the last one.
    ;locality of refrence
    project.imndArray[project.ni].file_Path = parent_directory
    project.current_Path = parent_directory

    ;try and extract the info from the imnd string. if the extraction of the
    ; imnd file returns 0 then the extraction was bad then display a message to the user
    ; and then return from open
    if 0 eq mas_extract_imnd ( imnd_file_string ) then begin
       update_status_bar, "Unable to extract 1D data"
       return
    end

    ;after we extract the imnd file we have to over write the adim
    ;the number of adim is the number of directories in the parent directory
    project.imndArray[project.ni].adim = n_parent_dir
    project.imndArray[project.ni].pn_avg = ptr_new(n_avg[sorted])
    project.imndArray[project.ni].bval_array  = ptr_new(bvals[sorted])
    project.imndArray[project.ni].b_matrix = ptr_new(b_matrix[*,sorted])

    ;store what type of open it was so that later on we can update it more easily.
    ;7 is ADC Multi
    project.imndArray[project.ni].image_type = 7

    ;since the file opened and was read completely change the current index to the
    ;recently opened file. then increment the next index
    project.ci = project.ni
    project.ni = project.ni + 1

    ;now that scan is open we can set the scan_open_flag
    project.scan_open_flag = 1

    cd, parent_directory

    mas_add_current_scan_to_scan_selection

    mas_redraw_GUI

    update_status_bar, ''

    if project.auto_Display eq 1 then mas_display

    HEAP_GC

 end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Subroutine name: mas_open_T2_Multi
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;

; Purpose of subroutine:
; HS 2006/09/29. Opens a series of T2 weighted images.


; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting

pro mas_open_T2_Multi
    COMPILE_OPT IDL2
    COMMON scan_data

    ;Get the proper slash depending upon which operation system is in use.
    slash = get_dir_slash()


    ;Check to make sure that max_num_scans was not exceeded.
    if project.ni eq project.max_num_scans then begin
       update_status_bar, 'Reached max limit of scans able to open.'
       return
    end


    ;Have the user select a directory that they want to try and open that scan.
    parent_directory = DIALOG_PICKFILE( Title = 'Select parent directory of the multiple T2 scans' ,$
                 path=project.current_path,DIRECTORY=1)

    ;If the directory returned is empty, then they closed the open without a selection then just return.
    if (parent_directory eq '') then return

    ;Get a listing of the files and directories that are in the parent directory.
    parent_file_listing = FILE_SEARCH(parent_directory+'*',  COUN=n_parent_dir)



    ;Next we need to extract the echo time so that the we can change the order of the multi_slice_dir
    ;to the proper read in order of fids
    ;First we will open each scan individually and extract the echo time.
    ;Then we will sort the echo time and that will tell us how to read in the T2 fid files
    ;for each scan.
    echo_times = fltarr( n_parent_dir )
    for ii=0, n_parent_dir-1 do begin



        imnd_file_string = ''
        imnd_temp_string = ''
        OPENR, imnd_file_descriptor, parent_file_listing[ii]+slash[1]+'imnd' , /GET_LUN , ERROR = err

        IF (err NE 0) then  begin

           	open_file  = parent_file_listing[ii]+slash[1] + '\method'
			OPENR, imnd_file_descriptor, open_file, /GET_LUN , ERROR = err

			IF (err NE 0) then begin
				; If err is nonzero, something happened. Print the error message to
       			; the standard error file (logical unit -2):
	           ;PRINTF, -2, !ERROR_STATE.MSG
	           update_status_bar,'Please open the directory containing imnd or method files'
	           return
        	end
        end

        ;if there was no error then the file's text can be extracted.
        While (eof(imnd_file_descriptor) NE 1) Do Begin
           Readf, imnd_file_descriptor, imnd_temp_string
            imnd_file_string = imnd_file_string + imnd_temp_string
        EndWhile
        close, imnd_file_descriptor
        FREE_LUN, imnd_file_descriptor

		temp = float(extract_simple_value (imnd_file_string ,'##$IMND_diff_echo_time=', '$'))
		if temp eq 0 then temp = float(extract_simple_value (imnd_file_string ,'##$PVM_EchoTime=', '#'))
		if temp eq 0 then temp = float(extract_simple_value (imnd_file_string ,'##$IMND_echo_time=', '$'))
		
    echo_times[ii] = temp

    end
    CLOSE, /all

    ;sort the echo time and use that as the key to the parent_file_listing
    project.imndArray[project.ni].multi_scan_file_array = ptr_new(parent_file_listing[SORT(echo_times)])

    ;store the echo times to the global structure
    project.imndArray[project.ni].echo_Time_ptr  = ptr_new(echo_times[SORT(echo_times)])

    ;next we need to save the reciently extracted string to the project structure
    project.imndArray[project.ni].imnd_File = imnd_file_string


    ;once the imnd string has been extracted then store where we found it at
    ;then update the path so that it contains the last file open
    ;that way when the next file is open it will be close to the last one.
    ;locality of refrence
    project.imndArray[project.ni].file_Path = parent_directory
    project.current_Path = parent_directory



    ;try and extract the info from the imnd string. if the extraction of the
    ; imnd file returns 0 then the extraction was bad then display a message to the user
    ; and then return from open
    if 0 eq mas_extract_imnd ( imnd_file_string ) then begin
       update_status_bar, "Unable to extract 1D data"
       return
    end

    ;after we extract the imnd file we have to over write the adim
    ;the number of adim is the number of directories in the parent directory
    project.imndArray[project.ni].adim = n_parent_dir

    ;store what type of open it was so that later on we can update it more easily.
    ;2 is T2 weighting across multiple scans
    project.imndArray[project.ni].image_type = 2

    ;tell mas how to load the state 1 data
    project.imndarray[project.ni].state1_load_procedure = 'mas_load_state_1_bruker'
    
    ;since the file opened and was read completely change the current index to the
    ;recently opened file. then increment the next index
    project.ci = project.ni
    project.ni = project.ni + 1

    ;now that scan is open we can set the scan_open_flag
    project.scan_open_flag = 1



    cd, parent_directory

    mas_add_current_scan_to_scan_selection

    mas_redraw_GUI

    update_status_bar, ''

    if project.auto_Display eq 1 then mas_display

    HEAP_GC

end

; Subroutine name: mas_open_onepulse_spect
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;

; Purpose of subroutine:
; BT 2010/09/25 - Read Bruker IMND-based ONEPULSE spect.


; Editing Information:


pro mas_open_onepulse_spect

    common scan_data, project
    
    ;Check to make sure that max_num_scans was not exceeded
    if project.ni eq project.max_num_scans then begin
        void = dialog_message('Maximum number of scans reached:'+strcompress(string(project.max_num_scans)), $
                              /center, /error)
        update_status_bar, ' Reached max limit of scans able to open.'
        return
    end

    ;Have the user select a directory that they want to try and open that scan.
    scan_directory = DIALOG_PICKFILE( FILTER = ['imnd'], $
                                      title = 'Select directory of the scan to open' ,$
                                      path=project.current_path, $
                                      DIRECTORY=1)

    ;If the directory returned is empty, then they closed the open without a selection then just return
    if (scan_directory eq '') then return


    ; Check to see if the directory contains a fid file.
    ; If it does then the directory is valid.
    ; If not then then send the user a message and then return from open.
    OPENR, descriptor, scan_directory + 'fid' , /GET_LUN , ERROR = err
                                ; If err is nonzero, something
                                ; happened. send information to
                                ; user
                                ; or maybe this is a ser file instead of a fid file.
    IF (err NE 0) then  begin
        OPENR, descriptor, scan_directory + 'ser' , /GET_LUN , ERROR = err
        IF (err NE 0) then  begin
                                ;PRINTF, -2, !ERROR_STATE.MSG
            update_status_bar,'Please open a directory containing a FID file or SER'
            return
        end
    end
    close, descriptor
    FREE_LUN, descriptor


    ;Verify is we have an IMND or METHOD file

    parameter_file_pointer = scan_directory + 'imnd'
    IMND_Exists = FILE_SEARCH(parameter_file_pointer)
    IMND_size = size(IMND_Exists)

    if IMND_size[0] eq 0 then begin

        void = dialog_message('Spect processing is only currently supported for IMND data sets.', $
                              /center, /error)
        return
        
    endif
    
    parameter_file_string = ''
    parameter_temp_string = ''
    OPENR, parameter_file_descriptor, parameter_file_pointer , /GET_LUN , ERROR = err
    ; If err is nonzero, something happened. Print the error message to
    ; the standard error file (logical unit -2):

    ;If there was no error then the file's text can be extracted.
    While (eof(parameter_file_descriptor) NE 1) Do Begin
       Readf, parameter_file_descriptor, parameter_temp_string
        parameter_file_string = parameter_file_string + parameter_temp_string
    EndWhile
    close, parameter_file_descriptor
    FREE_LUN, parameter_file_descriptor

    ;Finally we need to save the recently extracted string to the project structure.
    project.imndArray[project.ni].imnd_File = parameter_file_string

    ;Store what type of open it was so that later on we can update it if neccessary.
    project.imndArray[project.ni].image_type = 0
    ;print, 'project.imndArray[project.ni].image_type = ',project.imndArray[project.ni].image_type

    ;; tell mas how to load the state 1 data
    project.imndArray[project.ni].state1_load_procedure = 'mas_load_state_1_bruker'
    
    ;Once the imnd string has been extracted then store where we found it at.
    ;Then update the path so that it contains the last file open.
    ;That way when the next file is open it will be close to the last one.
    ;Locality of reference
    project.imndArray[project.ni].file_Path = scan_directory
    project.current_Path = scan_directory
    
    method = extract_name(parameter_file_string)
    if (method eq 'ONEPULSE') then begin
        project.imndArray[project.ni].scan_name = method
        project.imndArray[project.ni].scan_date = extract_date(parameter_file_string)
        project.imndArray[project.ni].scan_name = extract_name(parameter_file_string)
        project.imndArray[project.ni].display_Name = project.imndArray[project.ni].file_Path
    
        if (0 eq extract_bruker_onepulse(parameter_file_string)) then begin
            update_status_bar, "Unable to extract 1D data."
            return
        endif
        
    endif

    ;since the file opened and was read compleately change the current index to the
    ;recently opened file. Then increment the next index.
    project.ci = project.ni
    project.ni = project.ni + 1

    ;now that scan is open we can set the scan_open_flag
    project.scan_open_flag = 1

    cd, scan_directory

    mas_add_current_scan_to_scan_selection

    turnoff_motion_correction

    mas_redraw_GUI

    update_status_bar, ''

    if project.auto_Display eq 1 then mas_display

    HEAP_GC

end

; Subroutine name: mas_open_project
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; HS 2006/09/29. Opens a *.mas file type.

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting

PRO mas_open_project, Event
    COMPILE_OPT IDL2
    COMMON scan_data, project

    imnd_file_pointer = DIALOG_PICKFILE( FILTER = ['*.mas'], Title = 'Select mas project')

    if (imnd_file_pointer ne '') then begin
       restore, FILENAME=imnd_file_pointer
       ;if there is no open scan then dont do any thing
       mas_redraw_GUI
       update_scan_selection
    END

END


;pro extract_bval_direction_dicom , bvalStringArray
;    COMPILE_OPT IDL2
;
;    update_status_bar, 'extracting b-values and directions'
;    COMMON scan_data
;
;    CI = project.ci
;    bvals = bvalStringArray
;
;    ;project.  bvalArray
;
;    ;print, bvalStringArray
;
;
;    ;this is part of set up we need to find out how many times to loop
;    ; throug to get all the different bval strings out of the array
;    ; basically removing all the repeats.
;
;    sizeOfBval = (size(bvals))[1]
;    result =WHERE(bvals eq bvals [0])
;    nuberOfRepeatsForOneBval = (size(result))[1]
;    numberOfIterations = sizeOfBval/nuberOfRepeatsForOneBval
;    ;print, 'numberOfIterations', numberOfIterations
;
;    ;now we need to make an array to hold the bval strings that is the right size.
;    reducedBvalStringArray = strArr(numberOfIterations)
;
;
;    if numberOfIterations gt 1 then begin
;
;       for ii=0, numberOfIterations-1 do begin
;
;         element = bvals [0]
;         reducedBvalStringArray[ii]= element
;         result=WHERE(bvals ne element)
;         ;dont do the last one b/c nothing from nothing leaves nothing
;         ; and it wont let me condition of result == -1
;         ; so i just wont do the last one.
;         if ii ne numberOfIterations - 1 then $
;          bvals = bvals[result]
;         ;print,bvals
;       endfor
;    endif
;    if numberOfIterations eq 1 then begin
;       ;if they open the b 0 directory then there will only be 1 b-value
;
;       reducedBvalStringArray[0] = bvals [0]
;
;    endif
;
;;   print, 'reducedBvalStringArray ', reducedBvalStringArray
;;
;;   print, 'reducedBvalStringArray ',reducedBvalStringArray[0]
;
;
;
;    project.imndArray[CI].bval_Array = ptr_new(reducedBvalStringArray)
;
;
;    index1 = STRPOS (reducedBvalStringArray[0],'#' )
;;   print, index1
;
;    b_val = strmid(reducedBvalStringArray[0], 4 , index1-4 )
;;   print, b_val
;
;    update_status_bar, ''
;
;end

;this function will return a structure, described bellow, of the
;b_matrix and other geo_info file information at the specfied location
;
;geo_info = {
;        path        :geo_info_string_array[0]
;        f_dim       :f_dim
;        p_dim       :p_dim
;        s_dim       :s_dim
;        r_voxel     :r_voxel
;        p_voxel     :p_voxel
;        s_voxel     :s_voxel
;        euler_psi   :euler_psi
;        euler_theta :euler_theta
;        euler_phi   :euler_phi
;        method      :method
;        number_dwi  :number_dwi
;        number_bval :number_bval
;        b_values    :b_values
;        b_val       :b_val
;    }




; Subroutine name: get_geo_info_file
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; HS 2006/09/29.

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting
    ; Edited by CD 2007/09/11
    ; Overhauled, made into string parser

FUNCTION get_geo_info_file, location

    ;try and open the file
    OPENR, lun, location , /GET_LUN , ERROR = err

    ; If err is nonzero, something happened. Print the error message to the standard error file (logical unit -2):
    IF (err NE 0) THEN BEGIN
        void = DIALOG_MESSAGE(['Error in get_geo_info_file','extracting text from geo_info file!','Aborting open'],/error, /center)
        RETURN,-1
    ENDIF

	diff_data_set = 0
	multi_echo_set = 0
	multi_rep_set = 0
	param_file_text = ''

	line_of_text = ''
	WHILE EOF(lun) NE 1 DO BEGIN
		READF, lun, line_of_text

		; parse string
		str_data = STRSPLIT(line_of_text,'=',/extract)

		CASE str_data[0] OF
			'filename'			: filename = str_data[1]
			'f_dim'   			: f_dim = FLOAT(str_data[1])
			'p_dim'	  			: p_dim = FLOAT(str_data[1])
			's_dim'   			: s_dim = FLOAT(str_data[1])
			'r_voxel' 			: r_voxel = FLOAT(str_data[1])
			'p_voxel' 			: p_voxel = FLOAT(str_data[1])
			's_voxel' 			: s_voxel = FLOAT(str_data[1])
			'image_plane_psi'	: euler_psi = FLOAT(str_data[1])
			'image_plane_theta'	: euler_theta = FLOAT(str_data[1])
			'image_plane_phi'	: euler_phi = FLOAT(str_data[1])
			'method'			: method = str_data[1]
			'big_delta'			: big_delta = FLOAT(str_data[1])
			'small_delta'		: small_delta = FLOAT(str_data[1])
			'a_dim'				: a_dim = FLOAT(str_data[1])

			'Diffusion Weighting Values: ' : BEGIN
				diff_data_set = 1

				n_avgs 		= FLTARR(a_dim)
    			angle_theta = FLTARR(a_dim)
   			 	angle_phi 	= FLTARR(a_dim)
    			b_val 		= FLTARR(6,a_dim)

       			FOR ii=0, a_dim-1 DO BEGIN
					READF, lun, line_of_text

					; parse string
					flt_data = FLOAT(STRSPLIT(line_of_text,',',/extract))

			        n_avgs[ii] 		= flt_data[0]
         			angle_theta[ii] = flt_data[1]
        	 		angle_phi[ii] 	= flt_data[2]
         			b_val[*,ii] 	= flt_data[3:*]
       			ENDFOR
			END

			'Multiple Echo Times: '	: BEGIN
				multi_echo_set = 1

				READF, lun, line_of_text
				echo_times = FLOAT(STRSPLIT(line_of_text,',',/extract))
			END

			'Multiple Repetition Times: ' : BEGIN
				multi_rep_set = 1

				READF, lun, line_of_text
				rep_times = FLOAT(STRSPLIT(line_of_text,',',/extract))
			END

			'Start of PAR file:' : BEGIN
				WHILE EOF(lun) NE 1 DO BEGIN
					param_file_line = ''
					READF, lun, param_file_line
					param_file_text = param_file_text + STRING(13) + STRING(10) + param_file_line
				ENDWHILE
			END

			'Start of IMND or METHOD file of:' : BEGIN
				WHILE EOF(lun) NE 1 DO BEGIN
					param_file_line = ''
					READF, lun, param_file_line
					param_file_text = param_file_text + STRING(13) + STRING(10) + param_file_line
				ENDWHILE
			END

			ELSE: PRINT, 'UNKNOWN CONTENTS IN GEO_INFO FILE'
		ENDCASE
	ENDWHILE

	CLOSE, lun
    FREE_LUN, lun

	; Parse the parameter_file
	position1 = STRPOS(param_file_text,'----------------------')
	position1 = position1 + STRLEN(position1)
	position2 = STRPOS(param_file_text,'----------------------', position1)
	param_text = STRMID(param_file_text,position1,position2-position1)

    ;create and init the data structure that is going to be passed out of the function
    geo_info = { $
        path        :filename	 $
        ,f_dim      :f_dim  	 $
        ,p_dim      :p_dim  	 $
        ,s_dim      :s_dim   	 $
        ,r_voxel    :r_voxel  	 $
        ,p_voxel    :p_voxel   	 $
        ,s_voxel    :s_voxel     $
        ,euler_psi  :euler_psi   $
        ,euler_theta:euler_theta $
        ,euler_phi  :euler_phi   $
        ,method     :method      $
        ,a_dim 		:a_dim 		 $
        ,param_text	:param_text  $
        ,angle_theta:PTR_NEW()	 $
        ,angle_phi	:PTR_NEW()	 $
        ,n_avgs		:PTR_NEW()	 $
        ,b_matrix	:PTR_NEW()	 $
        ,small_delta:PTR_NEW()	 $
        ,big_delta	:PTR_NEW()	 $
        ,echo_times :PTR_NEW()	 $
        ,rep_times	:PTR_NEW()	 $
    }

	IF diff_data_set EQ 1 THEN BEGIN
		geo_info.angle_theta = PTR_NEW(angle_theta)
		geo_info.angle_phi	 = PTR_NEW(angle_phi)
		geo_info.n_avgs		 = PTR_NEW(n_avgs)
		geo_info.b_matrix	 = PTR_NEW(b_val)
		geo_info.small_delta = PTR_NEW(small_delta)
		geo_info.big_delta	 = PTR_NEW(big_delta)
	ENDIF

	IF multi_echo_set EQ 1 THEN BEGIN
		geo_info.echo_times  = PTR_NEW(echo_times)
	ENDIF

	IF multi_rep_set EQ 1 THEN BEGIN
		geo_info.rep_times  = PTR_NEW(rep_times)
	ENDIF

    RETURN, geo_info
END


; Subroutine name: extract_b_matrix
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This procedure will extract the b-matrix out of a set of directories that are in the parent directory.

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting
	;HS - 20061215
    ;Deleted a bunch of commented out code.

function extract_b_matrix , imnd, n_bval
    COMMON scan_data
    ci = project.ci
    slash = get_dir_slash()

              ;find the b-matrix
              search_strings =['IMND_diff_brr' $
                             ,'IMND_diff_bpp' $
                             ,'IMND_diff_bss' $
                             ,'IMND_diff_brp' $
                             ,'IMND_diff_brs' $
                             ,'IMND_diff_bps' ]

              b_matrix = strarr(6)

              for ii=0, 5 do begin

                 b_matrix[ii] = strjoin(extract_slice_scheme (imnd, search_strings[ii], n_bval ),' ')

              endfor

    return, b_matrix
end


; Subroutine name: extract_b_matrix_method_files
; Created by: HS 20061213
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; To extract the b matrix of METHOD files

; Editing Information:

FUNCTION extract_b_matrix_method_files, imnd_file
    COMPILE_OPT IDL2


search_parameter = '##$PVM_DwBMat='
index1 = Strpos(imnd_file, search_parameter)

if index1 eq -1 then begin
	return, -1
endif

index1 = Strpos(imnd_file, ' )', index1)+2
index2 = Strpos(imnd_file, '##$PVM' , index1 )
difference =  index2 - index1

temp = STRMID(imnd_file, index1, difference)
b_matrix_temp = Strsplit(temp, ' ',/EXTRACT)

n_matrices = n_elements(b_matrix_temp)/9

if (n_matrices gt 1) then begin
   bmt = reform(b_matrix_temp, 9, n_matrices)
   b_xx = bmt[0,*]
   b_yy = bmt[4,*]
   b_zz = bmt[8,*]
   b_xy = bmt[1,*]
   b_xz = bmt[2,*]
   b_yz = bmt[5,*]
   
   b_matrix = fltarr(6, n_matrices)
   b_matrix[0,*] = b_xx
   b_matrix[1,*] = b_yy
   b_matrix[2,*] = b_zz
   b_matrix[3,*] = b_xy
   b_matrix[4,*] = b_xz
   b_matrix[5,*] = b_yz

   return, b_matrix

endif else begin
   
   b_matrix = strarr(6)
   
; This b_matrix must be now subdivided into the 6 unique elements.
   
   b_xx = b_matrix_temp[0]
   b_xy = b_matrix_temp[1]
   b_xz = b_matrix_temp[2]
   b_yx = b_matrix_temp[3]
   b_yy = b_matrix_temp[4]
   b_yz = b_matrix_temp[5]
   b_zx = b_matrix_temp[6]
   b_zy = b_matrix_temp[7]
   b_zz = b_matrix_temp[8]
   
   b_matrix[0] = b_xx + ' '
   b_matrix[1] = b_yy + ' '
   b_matrix[2] = b_zz + ' '
   b_matrix[3] = b_xy + ' '
   b_matrix[4] = b_xz + ' '
   b_matrix[5] = b_yz + ' '
   
   
   return, b_matrix
   
endelse

END


; Subroutine name: mas_open_adt
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This procedure will extract the b-matrix out of a set of directories that are in the parent directory.

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting

pro mas_open_adt
    COMPILE_OPT IDL2
    COMMON scan_data

    ;check to make sure that max_num_scans was not exceeded
    if project.ni eq project.max_num_scans then begin
       update_status_bar, ' Reached max limit of scans able to open. 50 '
       return
    end

    slash = get_dir_slash()

    ;;have the user select a parent directory that contains all the ADT scan folders.
    parent_Path = DIALOG_PICKFILE( Title = 'Select parent directory to ADT data' ,$
                                   path=project.current_path, DIRECTORY=1)
    
    if parent_Path eq '' then return

    file_listing = FILE_SEARCH(parent_Path+'*')
    n_scans =  n_elements(file_listing)
    data_dir_check = bytarr(n_scans)
    
    for ii = 0, n_scans-1 do begin
        ok = stregex(file_basename(file_listing[ii]), '^[0-9]+$', /boolean)
        data_dir_check[ii] = ok
    endfor
    
    if (total(data_dir_check) eq 0) then begin
        junk = dialog_message(['This directory does not appear to contain any', $
                               'ADT data. Please make sure that the directory', $
                               'is an ADT scan.'], /error, /center)
        return
    endif

    parent_dir_file_listing = file_listing[where(data_dir_check ne 0)]
    n_scans = n_elements(parent_dir_file_listing)

;;;;;;;;;;;;;;;;;;;;;;;;;;;
    project.current_path = parent_Path

;    parent_dir_file_listing = FILE_SEARCH(parent_Path+'*')
;    n_scans =  (size(parent_dir_file_listing))[1]

    ;Check to see if all sub directories contain a fid and IMND file.
    update_status_bar,'Verifying fid and imnd presence'
    for ii=0, n_scans-1 do begin
        scan_directory = parent_dir_file_listing[ii]

        open_file  = scan_directory + '\fid'
        if !d.name eq 'X' then open_file  = scan_directory + '/fid'
       
        OPENR, descriptor, open_file , /GET_LUN , ERROR = err
        ;; If err is nonzero, something happened. Send information to user
        IF (err NE 0) then  begin
            ;;PRINTF, -2, !ERROR_STATE.MSG
            update_status_bar,'Error in open directory containing ADT scans'
            print,''
            void = dialog_message(['Directory ',scan_directory+ '\fid',$
                                   'Does not contain a fid file.', $
                                   'Aborting open'],/error)
            return
        endif
        close, descriptor
        FREE_LUN, descriptor

        open_file  = scan_directory + '\imnd'
        if !d.name eq 'X' then open_file  = scan_directory + '/imnd'
        
        OPENR, descriptor, open_file, /GET_LUN , ERROR = err
        
        ;; If err is nonzero, IMND does not exist.
        IF (err NE 0) then  begin
            
            open_file  = scan_directory + '\method'
            if !d.name eq 'X' then open_file  = scan_directory + '/method'
            OPENR, descriptor, open_file, /GET_LUN , ERROR = err
            
            ;;IMND or Method file does not exist.
            IF (err NE 0) then  begin
                
                update_status_bar,'Error in open directory containing ADT scans'
                print,''
                print,'Directory ',scan_directory+ ' does not contain an imnd or method file.'
                print,'Aborting open'
                return
            endif
            
        endif
        
        close, descriptor
        FREE_LUN, descriptor

    endfor
    
    update_status_bar,''

    ;Now open the imnd file and take all the data out of it and store it to a string.
    ;If done completely then we are pretty sure that this data is ok so we will increase the
    ;next_index ++
    imnd_file_string = ''
    imnd_temp_string = ''

	; HS
	; no need to destroy the check that we did before this

    ;open_file = scan_directory + '\imnd'
    ;if !d.name eq 'X' then open_file  = scan_directory + '/imnd'

    OPENR, imnd_file_descriptor, open_file , /GET_LUN , ERROR = err
    ; If err is nonzero, something happened. Print the error message to
    ; the standard error file (logical unit -2):

    IF (err NE 0) then  begin
       ;PRINTF, -2, !ERROR_STATE.MSG
       update_status_bar,'Error extracting text from imnd file.'
       return
    end

    ;if there was no error then the file's text can be extracted.
    While (eof(imnd_file_descriptor) NE 1) Do Begin
       Readf, imnd_file_descriptor, imnd_temp_string
        imnd_file_string = imnd_file_string + imnd_temp_string
    EndWhile

    close, imnd_file_descriptor
    FREE_LUN, imnd_file_descriptor



    ; Try and extract the info from the imnd string. If the extraction of the
    ; imnd file returns 0 then the extraction was bad then display a message to the user
    ; and then return from open.
    if 0 eq mas_extract_imnd ( imnd_file_string ) then begin
       update_status_bar, "Unable to extract 1D data"
       return
    end


    ;Once one imnd file is open then we can open all the imnd files and extract just the adim.
    ;the reason we do this is b/c some scans have a b=0 scan which dosn't have a matching adim
    ;to all the other scans in the directory so we have to be a little more careful when we open
    ;the scans. So now we have an adim counter that is incremented by the adim of the current
    ;open scan.


	;since the file opened and was read completely change the current index to the
    ;recently opened file. then increment the next index

    project.ci = project.ni
    project.ni = project.ni + 1

    adim = 0
    file_adim_str = ''
    n_avg = fltarr(n_scans)
    diff_angle_theta_str = ''
    diff_angle_Phi_str = ''
    b_matrix_string = strarr(6)
    bvalues = fltarr(n_scans)
    bdeltas = fltarr(n_scans)
    sdeltas = fltarr(n_scans)

    ;;common diff_grad_test, diff_grad_xyz
    ;;diff_grad_xyz = fltarr(3, n_scans)

    ;; HS 20061214
    ;; Lets check whether the file to be searched contains the IMND term.
    imnd_check = Strpos(imnd_file_string, 'IMND_')

        
    for ii=0, n_scans-1 do begin
        scan_directory = parent_dir_file_listing[ii]

        ;;Now open the imnd file and take all the data out of it and store it to a string.
        ;;If done completly then we are pretty sure that this data is ok so we will increase the
        ;;next_index ++
        imnd_file_string = ''
        imnd_temp_string = ''
        
        open_file = scan_directory + '\imnd'
        
        if !d.name eq 'X' then open_file  = scan_directory + '/imnd'
        OPENR, imnd_file_descriptor, open_file , /GET_LUN , ERROR = err
        ; If err is nonzero, something happened. Print the error message to
        ; the standard error file (logical unit -2):

        IF (err NE 0) then begin

 			open_file  = scan_directory + '\method'
			if !d.name eq 'X' then open_file  = scan_directory + '/method'
			OPENR, descriptor, open_file, /GET_LUN , ERROR = err

         ;IMND or Method file does not exist.
           IF (err NE 0) then begin

           ;PRINTF, -2, !ERROR_STATE.MSG
           update_status_bar,'Error extracting text from imnd or method file.'
           return
           end

        end

        ;if there was no error then the file's text can be extracted.

        While (eof(imnd_file_descriptor) NE 1) Do Begin
           Readf, imnd_file_descriptor, imnd_temp_string
            imnd_file_string = imnd_file_string + imnd_temp_string
        EndWhile

        close, imnd_file_descriptor
        FREE_LUN, imnd_file_descriptor


                                ; Different algorithms to extract the b matrix for METHOD or IMND files.
        if imnd_check eq -1 then begin
           
           tmp_n_avg = extract_simple_value( imnd_file_string ,'##$PVM_NAverages=','#' )
           tmp_n_reps = extract_simple_value( imnd_file_string ,'##$PVM_NRepetitions=','#' )
           temp = extract_simple_value( imnd_file_string ,'##$PVM_DwNDiffExp=','#')
           if (temp eq 0) then begin
                temp = extract_simple_value( imnd_file_string ,'##$PVM_DwNDiffExp=','$$')
           endif
                                ; We need to calculate the correct theta and phi angles based on the new rectangular
                                ; coordinate definition of gradient direction.
           
           phi_tmp = strjoin(extract_phi_methods(imnd_file_string), ' ')
           theta_tmp = strjoin(extract_theta_methods(imnd_file_string), ' ')
           
           if (temp gt 1) then begin
               diff_angle_phi_str += strjoin(replicate(phi_tmp, temp), ' ')
               diff_angle_theta_str += strjoin(replicate(theta_tmp, temp), ' ')
           endif else begin
               diff_angle_phi_str += phi_tmp
               diff_angle_theta_str += theta_tmp
           endelse
                                ; Now lets extract the b matrix from the METHOD file
           b_matrix_temp = extract_b_matrix_method_files(imnd_file_string)
           
                                ; CD extract the b value
           tmp_bvalues = extract_bval_method( imnd_file_string )
           tmp_bdeltas = extract_big_delta_method( imnd_file_string )
           tmp_sdeltas = extract_small_delta_method( imnd_file_string)
           ;;diff_grad_xyz[*,ii] = extract_grad_mat_method(imnd_file_string)

        endif else begin
           
           tmp_n_avg = extract_simple_value( imnd_file_string ,'##$IMND_n_averages=','$' )
           tmp_n_reps = extract_simple_value( imnd_file_string ,'##$IMND_n_repetitions=','$' )
           temp = extract_simple_value( imnd_file_string ,'##$IMND_n_diff_expts=','$')
                                ;Now extract the diff angles
           diff_angle_theta_str += strjoin(extract_slice_scheme (imnd_file_string, 'diff_angle_theta', temp ),' ')
           diff_angle_phi_str += strjoin(extract_slice_scheme (imnd_file_string, 'diff_angle_phi', temp ),' ')
           
                                ;Extract the b_matrix for the given imnd
           b_matrix_temp = extract_b_matrix( imnd_file_string, temp)
           
                                ; CD extract the b value
           tmp_bvalues = extract_bval_imnd( imnd_file_string )
           tmp_bdeltas = extract_simple_value( imnd_file_string ,'##$IMND_big_delta=','$' )
           tmp_sdeltas = extract_simple_value( imnd_file_string ,'##$IMND_diff_grad_dur=','$' )
           ;;diff_grad_xyz[*,ii] = extract_grad_mat_imnd(imnd_file_string)

        endelse
        
        
		; Adding the b matrix extracted to the list of all the b matrices.
		
      	for jj=0, 5 do begin
      	     btt = strjoin(reform(b_matrix_temp[jj,*]), ' ')
        	 b_matrix_string[jj] += btt;;b_matrix_temp[jj]
		endfor

		adim += temp

        if (ii eq 0) then begin
            n_avg = tmp_n_avg
            n_reps = tmp_n_reps
            bvalues = float(tmp_bvalues)
            bdeltas = float(tmp_bdeltas)
            sdeltas = float(tmp_sdeltas)
        endif else begin
            n_avg = [ n_avg, tmp_n_avg ]
            n_reps = [ n_reps, tmp_n_reps ]
            bvalues = [ bvalues, float(tmp_bvalues) ]
            bdeltas = [ bdeltas, float(tmp_bdeltas) ]
            sdeltas = [ sdeltas, float(tmp_sdeltas) ]
        endelse
        
       ;Store the adim for each file so when we open them
       	file_adim_str += string(temp)+' '

    end ;end of the FOR LOOP for every dwi direction



    ;print, 'mas_open_adt', adim

    ;Check the content of the b-matrix, count the number of -1 and if there are any throw an error.
    print, where(b_matrix_string eq -1 , count)
    if count gt 0 then begin
       void = dialog_message (['B-matrix looks funny ', $
                               'Xeve says "Make sure you specify tensor instead of coefficient"', $
                               'Bill says, "Just use the data with a little caution."'],/error)
                               
       ;project.imndArray[project.ci].b_matrix = ptr_new()
       ;return
    end

    ;print, b_matrix_string
    ;now that we have the correct number for the adim and our b_matrix in a string
    ;it is time to convert the b-matrix to a floating point array.
    b_matrix = fltarr(6,adim)
    for jj=0, 5 do begin
        btmp =  rotate(float(strsplit(b_matrix_string[jj], /extract)),1)
        b_matrix[jj,*] = btmp[where(btmp ne -1)]
    endfor

    ;print, b_matrix
    ;multiply the off diagionals by 2
    b_matrix[3:5,*] = b_matrix[3:5,*] * 2

    project.imndArray[project.ci].b_matrix = ptr_new(b_matrix)
    project.imndArray[project.ci].bval_Array = PTR_NEW(bvalues)
    project.imndArray[project.ci].big_delta = PTR_NEW(bdeltas)
	project.imndArray[project.ci].small_delta = PTR_NEW(sdeltas)

    ;store the pointer to the number of averages to each scan.
    project.imndArray[project.ci].pn_avg = ptr_new(n_avg)

    ;store the angles for the diffusion experiment
    ;angle_theta:ptr_new(), angle_phi:

;    angle_theta = rotate(float(strsplit(diff_angle_theta_str, /extract)),1)
;    angle_phi   = rotate(float(strsplit(diff_angle_phi_str, /extract)),1)
 
;    grad = fltarr(3, n_scans)
;    grad[0,*] = angle_phi
;    grad[1,*] = angle_theta
;    grad[2,*] = 1.0
;    grad_rect = cv_coord(from_sphere=grad, /to_rect, /degrees)
;    grad_rect = vert_t3d(grad_rect, matrix=*project.imndarray[project.ci].acq_matrix, /no_copy)
;    grad = cv_coord(from_rect=grad_rect, /to_sphere, /degrees)
    
;    project.imndArray[project.ci].angle_theta = ptr_new(reform(grad[1,*]))
;    project.imndArray[project.ci].angle_phi = ptr_new(reform(grad[0,*]))
    
    project.imndArray[project.ci].angle_theta = ptr_new(rotate(float(strsplit(diff_angle_theta_str, /extract)),1))
    project.imndArray[project.ci].angle_phi = ptr_new(rotate(float(strsplit(diff_angle_phi_str, /extract)),1))


    ;finally we need to save the recently extracted string to the project structure
    project.imndArray[project.ci].imnd_File = imnd_file_string
    project.imndArray[project.ci].multi_scan_file_array = ptr_new(parent_dir_file_listing)
    project.imndArray[project.ci].file_adim = ptr_new(fix(strsplit(file_adim_str, /extract)))
    ;print, (*project.imndArray[project.ci].file_adim)

    ;once the imnd string has been extracted then store where we found it at.
    ;Then update the path so that it contains the last file open. That way when
    ;the next file is open it will be close to the last one.
    ;locality of reference
    project.imndArray[project.ci].file_Path = parent_Path
    project.current_Path = parent_Path

    ;store what type of open it was so that later on we can update it more easily.
    project.imndArray[project.ci].image_type = 3


    ;now that scan is open we can set the scan_open_flag
    project.scan_open_flag = 1


    ;The user is trying to open an ADT image set so the b-matrix is spread across a couple of directories.
    ;after the imnd data has been extracted then we need to go back and build up the
    ;  b_matrix

    project.imndArray[project.ci].adim = adim
    project.imndArray[project.ci].n_bvals = project.imndArray[project.ci].adim
    project.imndArray[project.ci].state1_load_procedure = 'mas_load_state_1_bruker'

    ;openADTMulti, parent_Path

    mas_add_current_scan_to_scan_selection

    mas_redraw_GUI

    update_status_bar, ''

    if project.auto_Display eq 1 then mas_display

end


; Subroutine name:  mas_open_phase_array
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This program will open phase array data folders and check for a valid FID and IMND file.

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting

pro mas_open_phase_array
    COMPILE_OPT IDL2
    COMMON scan_data

    ;check to make sure that max_num_scans was not exceeded
    if project.ni eq project.max_num_scans then begin
       update_status_bar, ' reached max limit of scans able to open. 50  '
       return
    end

    ;have the user select a parent directory that contains all the ADT scan folders.
    parent_Path = DIALOG_PICKFILE( Title = 'Select parent directory to ADT data' ,$
                             path=project.current_path, DIRECTORY=1)
    if parent_Path eq '' then return

    project.current_path = parent_Path

    parent_dir_file_listing = FILE_SEARCH(parent_Path+'*')
    n_scans =  (size(parent_dir_file_listing))[1]

    ;Check to see if all sub directories contain a fid and IMND file.
    update_status_bar,'verifying fid and imnd presence'
    for ii=0, n_scans-1 do begin
       scan_directory = parent_dir_file_listing[ii]

       open_file  = scan_directory + '\fid'
       if !d.name eq 'X' then open_file  = scan_directory + '/fid'

       OPENR, descriptor, open_file , /GET_LUN , ERROR = err
       ; If err is nonzero, something happened. send an information to user
       IF (err NE 0) then  begin
         ;PRINTF, -2, !ERROR_STATE.MSG
         update_status_bar,'Error in open directory containing ADT scans'
         void = dialog_message (['Directory '+scan_directory+'\fid',' does not contain a fid file.','Aborting open.'],/error)
         return
       end
       close, descriptor
       FREE_LUN, descriptor

       open_file  = scan_directory + '\imnd'
       if !d.name eq 'X' then open_file  = scan_directory + '/imnd'

       OPENR, descriptor, open_file, /GET_LUN , ERROR = err
       ; If err is nonzero, something happened. send information to user
       IF (err NE 0) then  begin
         ;PRINTF, -2, !ERROR_STATE.MSG
         update_status_bar,'Error in open directory containing ADT scans'
         void = dialog_message (['Directory ',scan_directory+ '\imnd',' doesnt contain an imnd file.','Aborting open'],/error)
         return
       end
       close, descriptor
       FREE_LUN, descriptor

    endfor

    update_status_bar,''


    ;now open the imnd file and take all the data out of it and store it to a string.
    ;if done completely then we are pretty sure that this data is ok so we will increase the
    ;next_index ++
    imnd_file_string = ''
    imnd_temp_string = ''

    open_file = scan_directory + '\imnd'
    if !d.name eq 'X' then open_file  = scan_directory + '/imnd'
    OPENR, imnd_file_descriptor, open_file , /GET_LUN , ERROR = err
    ; If err is nonzero, something happened. Print the error message to
    ; the standard error file (logical unit -2):
    IF (err NE 0) then  begin
       ;PRINTF, -2, !ERROR_STATE.MSG
       update_status_bar,'Error extracting text from imnd file.'
       return
    end
    ;if there was no error then the file's text can be extracted.
    While (eof(imnd_file_descriptor) NE 1) Do Begin
       Readf, imnd_file_descriptor, imnd_temp_string
        imnd_file_string = imnd_file_string + imnd_temp_string
    EndWhile
    close, imnd_file_descriptor
    FREE_LUN, imnd_file_descriptor

    ; Try and extract the info from the imnd string. If the extraction of the
    ; imnd file returns 0 then the extraction was bad. If this happens, display a message to the user
    ; and then return from open
    if 0 eq mas_extract_imnd ( imnd_file_string ) then begin
       update_status_bar, "Unable to extract 1D data"
       return
    end

    ;since the file opened and was read completely change the current index to the
    ;recently opened file, then increment the next index.
    project.ci = project.ni
    project.ni = project.ni + 1
    ;finally we need to save the recently extracted string to the project structure
    project.imndArray[project.ci].imnd_File = imnd_file_string
    project.imndArray[project.ci].multi_scan_file_array = ptr_new(parent_dir_file_listing)


    ;setting the nuber of adim to the number of scans in the directory
    project.imndArray[project.ci].adim = n_scans

    ;Once the imnd string has been extracted then store where we found it at.
    ;Then update the path so that it contains the last file open.
    ;That way when the next file is open it will be close to the last one.
    ;locality of reference
    project.imndArray[project.ci].file_Path = parent_Path
    project.current_Path = parent_Path

    ;store what type of open it was so that later on we can update it more easily.
    project.imndArray[project.ci].image_type = 4

    ;now that scan is open we can set the scan_open_flag
    project.scan_open_flag = 1

    mas_add_current_scan_to_scan_selection

    mas_redraw_GUI

    update_status_bar, ''

    if project.auto_Display eq 1 then mas_display


end

; Subroutine name:  bruker_image_file_presence
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Checks the existence of an imnd, reco, d3proc, 2dseq file.

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting
    ;Fixed the misspelled "burker_image_file_presence" in mas_project.pro

function bruker_image_file_presence , directory

    slash = get_dir_slash()

    image_directory = directory

    ;now that we have a directory we have to be sure that this directory is a good one
    ;in order to do so we have to check for the presence of 4 different files.
    ;1. imnd to get all the usual parameters for the scan
    ;2. reco to get the word type when opening the file
    ;3. d3proc to get the actual dimensions of the file
    ;4. 2dseq this is where the data is stored.


    ;print, image_directory
    ;checking for an imnd file

    split_directory = STRSPLIT(image_directory, slash[1], /EXTRACT)
    ;print, split_directory
    ;now that the directory is split up we have to go 2 levels and to find the imnd file
    dir_depth = (size(split_directory))[1]
    split_directory = split_directory[0:dir_depth-3]
    ;now that we have the correct directory we can now put the path back
    ;together and tack on the appropriate slash to the end
    scan_dir = slash[0]+(strjoin(split_directory,slash[1]))+slash[1]
    ;print, scan_dir
    ;now we check to see if the imnd file is there if it isn't then quit.
    fileParameters = FILE_INFO(scan_dir+'imnd' )
    if fileParameters.EXISTS eq 0 then begin
       void = dialog_message(['imnd file not found!',scan_dir+'imnd'],/error)
       return,0
    end


    ;now we check for the reco file to see if it exists
    fileParameters = FILE_INFO(image_directory+'reco' )
    if fileParameters.EXISTS eq 0 then begin
       void = dialog_message(['reco file not found!',image_directory+'reco'],/error)
       return,0
    end

    ;now we check for the d3proc file to see if it exists
    fileParameters = FILE_INFO(image_directory+'d3proc' )
    if fileParameters.EXISTS eq 0 then begin
       void = dialog_message(['d3proc file not found!',image_directory+'d3proc'],/error)
       return,0
    end

    ;now we check for the 2dseq file to see if it exists
    fileParameters = FILE_INFO(image_directory+'2dseq' )
    if fileParameters.EXISTS eq 0 then begin
       void = dialog_message(['d3proc file not found!',image_directory+'2dseq'],/error)
       return,0
    end


    ;if we get to this point then the files are ok
    return, 1

end


; Subroutine name: extract_imnd_file
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This function will extract the imnd data out of the file and place it up on the display
; and fill the structure.

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting


function extract_imnd_file, file_pointer
    COMMON scan_data

    imnd_file_pointer = file_pointer
    ;now open the imnd file and take all the data out of it and store it to a string.
    ;if done completely then we are pretty sure that this data is ok so we will increase the
    ;next_index ++
    imnd_file_string = ''
    imnd_temp_string = ''
    OPENR, imnd_file_descriptor, imnd_file_pointer , /GET_LUN , ERROR = err
    ; If err is nonzero, something happened. Print the error message to
    ; the standard error file (logical unit -2):
    IF (err NE 0) then  begin
       ;PRINTF, -2, !ERROR_STATE.MSG
       update_status_bar,'Unable to open the imnd'
       return,0
    end
    ;if there was no error then the file's text can be extracted.
    While (eof(imnd_file_descriptor) NE 1) Do Begin
       Readf, imnd_file_descriptor, imnd_temp_string
        imnd_file_string = imnd_file_string + imnd_temp_string
    EndWhile
    close, imnd_file_descriptor
    FREE_LUN, imnd_file_descriptor

    ;finally we need to save the recently extracted string to the project structure
    project.imndArray[project.ni].imnd_File = imnd_file_string

    ; Try and extract the info from the imnd string. If the extraction of the
    ; imnd file returns 0 then the extraction was bad. If this happens, display a message to the user
    ; and then return from open.
    if 0 eq mas_extract_imnd ( imnd_file_string ) then begin
       update_status_bar, "Unable to extract 1D data"
       return,0
    end

    return, 1
end

; Subroutine name: get_reco_file_word_size
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Returns the word size from the specified file pointer, if there is an error
; return 0

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting

function get_reco_file_word_size , file_pointer

    ;first we are going to have to extract the reco file in its entirety
    reco_file_string = ''
    reco_temp_string = ''
    OPENR, reco_file_descriptor, file_pointer , /GET_LUN , ERROR = err
    ; If err is nonzero, something happened. Print the error message to
    ; the standard error file (logical unit -2):
    IF (err NE 0) then  begin
       ;PRINTF, -2, !ERROR_STATE.MSG
       update_status_bar,'Unable to open the reco file'
       return,0
    end
    ;if there was no error then the file's text can be extracted.
    While (eof(reco_file_descriptor) NE 1) Do Begin
       Readf, reco_file_descriptor, reco_temp_string
        reco_file_string = reco_file_string + reco_temp_string
    EndWhile
    close, reco_file_descriptor
    FREE_LUN, reco_file_descriptor

    ;after we have the reco file we will have try and extract out the word size
    ;of the file.


    ;##$RECO_wordtype=_32BIT_SGN_INT
    ;##$RECO_wordtype=_16BIT_SGN_INT

    index1 = strpos(reco_file_string,'##$RECO_wordtype=')+strlen('##$RECO_wordtype=')
    index2 = strpos(reco_file_string,'##$',index1)
    s_word_size = strmid(reco_file_string, index1, index2-index1)

    ;bruker data types dont match closely to idl's data types so
    ;a conversion is needed to let idl know what idl data type it is
    if s_word_size eq '_32BIT_SGN_INT' then return, 3
    if s_word_size eq '_16BIT_SGN_INT' then return, 2

    return,0
end


; Subroutine name: get_d3proc_dimension
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This function will open the d3proc file and extract the actual freq phase and slice dimensions.

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting

function get_d3proc_dimension, file_pointer

    ;first we are going to have to extract the reco file in it's entirety
    file_string = ''
    temp_string = ''
    OPENR, file_descriptor, file_pointer , /GET_LUN , ERROR = err
    ; If err is nonzero, something happened. Print the error message to
    ; the standard error file (logical unit -2):
    IF (err NE 0) then  begin
       ;PRINTF, -2, !ERROR_STATE.MSG
       update_status_bar,'Unable to open the d3proc file'
       return,[0,0,0]
    end
    ;if there was no error then the file's text can be extracted.
    While (eof(file_descriptor) NE 1) Do Begin
       Readf, file_descriptor, temp_string
        file_string = file_string + temp_string
    EndWhile
    close, file_descriptor
    FREE_LUN, file_descriptor

    ;now that we have the file string we can look for specfic key words and
    ;a extract the values that are at the end of those key words.
    ;first we will extract the freq. direction
    position1 = strpos(file_string,'##$IM_SIX=')+10
    position2 = strpos(file_string,'##$',position1)
    ;print, position1,position2
    freq = fix(strmid(file_string,position1,position2-position1))

    ;now we extract the phase dimension
    position1 = strpos(file_string,'##$IM_SIY=')+10
    position2 = strpos(file_string,'##$',position1)
    ;print, position1,position2
    phase = fix(strmid(file_string,position1,position2-position1))

    ;now we extract the slice dimension
    position1 = strpos(file_string,'##$IM_SIZ=')+10
    position2 = strpos(file_string,'##$',position1)
    ;print, position1,position2
    slice = fix(strmid(file_string,position1,position2-position1))

    return, [freq,phase,slice]
end



; Subroutine name:  mas_open_bruker_image
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This pro will allow the user to open a single bruker image file

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting

pro mas_open_bruker_image
    COMPILE_OPT IDL2
    COMMON scan_data
    ;setting the type of slash to use for the directories
    slash = get_dir_slash()

    update_status_bar, 'Select processed image directory'

    ;check to make sure that max_num_scans was not exceeded
    if project.ni eq project.max_num_scans then begin
       update_status_bar, ' Reached max limit of scans able to open.'
       return
    end

    ;have the user select a directory that they want to try and open that scan
    image_directory = DIALOG_PICKFILE( Title = 'Select directory of the image file to open' ,$
                 path=project.current_path,DIRECTORY=1)

    ;If the directory returned is empty, then they closed the open without a selection then just return
    if (image_directory eq '') then return

    ;print, image_directory

    project.imndArray[project.ni].file_Path = image_directory
    project.current_Path = image_directory

    ;since we have a image file directory lets try and get the scans directory
    split_directory = STRSPLIT(image_directory, slash[1], /EXTRACT)
    ;now that the directory is split up we have to go 2 levels and to find the imnd file
    dir_depth = (size(split_directory))[1]
    split_directory = split_directory[0:dir_depth-3]
    ;now that we have the correct directory we can now put the path back
    ;together and tack on the appropriate slash to the end
    scan_dir = slash[0]+(strjoin(split_directory,slash[1]))+slash[1]

    update_status_bar, 'Verifying file presence'

    ;check the presence of all the files neccessary to open the scan
    if bruker_image_file_presence (image_directory) eq 0 then begin
       void = dialog_message('Error opening. Aborting',/error)
       return
    end


    ;now that we are sure that all the files are there we have to get the parameters out of the files
    ;1. imnd general file parameters
    ;2. reco to get the word size
    ;3. d3proc we need to get the actual size of the processed arrays
    if  extract_imnd_file( scan_dir+'imnd') eq 0 then begin
       void = dialog_message(['Error extracting imnd data from file',scan_dir+'imnd'],/error)
       return
    end

    ;now get the reco word size for the file.
    word_size = get_reco_file_word_size ( image_directory+'reco')
    if word_size eq '0' then begin
       void = dialog_message('Error opening. Aborting',/error)
       return
    end

    dimensions = get_d3proc_dimension ( image_directory+'d3proc')
    if dimensions[0]  eq '0' then begin
       void = dialog_message('Error opening. Aborting',/error)
       return
    end


    update_status_bar,'Opening processed image'


    ;now that we have the dimesion we can start extracting the data
    OPENR, file_descriptor, image_directory+'2dseq' , /GET_LUN , ERROR = err
    ; If err is nonzero, something happened. Print the error message to
    ; the standard error file (logical unit -2):
    IF (err NE 0) then  begin
       ;PRINTF, -2, !ERROR_STATE.MSG
       update_status_bar,'Unable to open the d3proc file'
       return
    end

    sdim = project.imndArray[project.ni].sdim
    ;print, sdim

    ;print, dimensions

    data= fltarr(dimensions[0],dimensions[1],dimensions[2])


    if mas_chk_bruker_version_3( project.imndArray[project.ni].imnd_file) then Begin
       print, 'wordsize', word_size

        data[*,*,*] = READ_BINARY (file_descriptor, $
           DATA_DIMS = dimensions , $
           ;ENDIAN='big' , $
           DATA_TYPE=word_size $
           )

    endif else begin
       print, 'wordsize', word_size

       data[*,*,*] = READ_BINARY (file_descriptor, $
           DATA_DIMS = dimensions , $
           ;ENDIAN='big' , $
           ENDIAN='little' , $
           DATA_TYPE=word_size $
           )
           print,'Commented out the big endian preset.mas_open_bruker_image'

    endelse

    close, file_descriptor
    FREE_LUN, file_descriptor

    adim = dimensions[2]/sdim


    data_out = fltarr(dimensions[0],dimensions[1],sdim, adim)
    for ii=0 , adim -1 do $
       data_out[*,*,*,ii] = data[*,*,ii*sdim:ii*sdim+sdim-1]

    project.dataArray[project.ni].state1 = ptr_new(data_out)

    project.imndArray[project.ni].fdim = dimensions[0]
    project.imndArray[project.ni].pdim = dimensions[1]
    project.imndArray[project.ni].sdim = sdim
    project.imndArray[project.ni].adim = adim

    project.procPramArray[project.ni].state_1 = 1
    project.imndArray[project.ni].image_type = 5

    ;Now that scan is open we can set the scan_open_flag
    project.scan_open_flag = 1

    ;Since the file opened and was read completely change the current index to the
    ;recently opened file. Then increment the next index.
    project.ci = project.ni
    project.ni = project.ni + 1


    cd, image_directory

    mas_add_current_scan_to_scan_selection

    mas_redraw_GUI

    update_status_bar, ''

    if project.auto_Display eq 1 then mas_display

    HEAP_GC

end

; Subroutine name:  get_image_directories
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This function will return the directories of the image files that need to be opened

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting

function get_image_directories , root_dir, image_dir
    ;setting the type of slash to use for the directories
    slash = get_dir_slash()

    ;scan root for directories that may contain scans.
    multi_scan_dir_listing = FILE_SEARCH(root_dir+'*', COUNT=num_scans )
    image_dir_array = strarr(num_scans)

    ;before we begin we have to know what directory number we are dealing with
    split_directory = STRSPLIT(image_dir, slash[1], /EXTRACT)
    dir_depth = (size(split_directory))[1]
    post_dir_wo_slashes = (split_directory[dir_depth-2:dir_depth-1])
    post_dir_w_slashes = slash[1]+(strjoin(post_dir_wo_slashes,slash[1]))+slash[1]
    ;print, 'post_part  ',post_part


    ;now traverse multi_scan_dir_listing piecing together where the image files should be
    for ii=0, num_scans-1 do begin
         multi_scan_dir_listing[ii] = multi_scan_dir_listing[ii]+post_dir_w_slashes
    endfor

    return ,multi_scan_dir_listing
end


; Subroutine name:  mas_open_bruker_multi_image
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This procedure will allow users to open bruker image across multiple scan and put them together.

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting

pro mas_open_bruker_multi_image
    COMPILE_OPT IDL2
    COMMON scan_data
    ;setting the type of slash to use for the directories
    slash = get_dir_slash()

    update_status_bar, 'Select processed image directory'

    ;check to make sure that max_num_scans was not exceeded.
    if project.ni eq project.max_num_scans then begin
       update_status_bar, ' Reached max limit of scans able to open.'
       return
    end

    ;Have the user select a directory that they want to try and open that scan
    image_directory = DIALOG_PICKFILE( Title = 'Select directory of one of the image files to open' ,$
                 path=project.current_path,DIRECTORY=1)

    ;If the directory returned is empty, then they closed the open without a selection then just return
    if (image_directory eq '') then return

    ;print, image_directory

    project.imndArray[project.ni].file_Path = image_directory
    project.current_Path = image_directory

    ;since we have a image file directory lets try and get the scans directory
    split_directory = STRSPLIT(image_directory, slash[1], /EXTRACT)
    ;now that the directory is split up we have to go 2 levels and to find the imnd file
    dir_depth = (size(split_directory))[1]
    split_directory = split_directory[0:dir_depth-4]
    ;now that we have the correct directory we can now put the path back
    ;together and tack on the appropriate slash to the end
    multi_scan_dir = slash[0]+(strjoin(split_directory,slash[1]))+slash[1]
    split_directory = STRSPLIT(image_directory, slash[1], /EXTRACT)
    split_directory = split_directory[0:dir_depth-3]
    scan_dir = slash[0]+(strjoin(split_directory,slash[1]))+slash[1]


    ;now that we have some general directories we need to try and find all
    ;the image file directories that are going to be opened.
    multi_image_directory = get_image_directories(multi_scan_dir, image_directory )
    num_scans = (size(multi_image_directory))[1]


    ;check the presence of all the files neccessary to open the scan
    if bruker_image_file_presence (image_directory) eq 0 then begin
       void = dialog_message('Error opening aborting',/error)
       return
    end


    ;now that we are sure that all the files are there we have to get the parameters out of the files
    ;1. imnd general file parameters
    ;2. reco to get the word size
    ;3. d3proc we need to get the actual size of the processed arrays
    if  extract_imnd_file( scan_dir+'imnd') eq 0 then begin
       void = dialog_message(['Error extracting imnd data from file',scan_dir+'imnd'],/error)
       return
    end

    ;now get the reco word size for the file.
    word_size = get_reco_file_word_size ( image_directory+'reco')
    if word_size eq '0' then begin
       void = dialog_message('Error opening aborting',/error)
       return
    end

    dimensions = get_d3proc_dimension ( image_directory+'d3proc')
    if dimensions[0]  eq '0' then begin
       void = dialog_message('Error opening aborting',/error)
       return
    end

    ;print, 'gathered dims from d3proc',dimensions

    ;now that we have the dimensions we can make an array large enough to hold
    ;all the data.

    sdim = project.imndArray[project.ni].sdim
    adim = (dimensions[2]/sdim)*num_scans


    ;places to hold data extraced from files.
    data= fltarr(dimensions[0],dimensions[1],dimensions[2])


    project.dataArray[project.ni].state1 = ptr_new(fltarr(dimensions[0],dimensions[1],sdim, adim))

    progressbar = Obj_New('progressbar', Color='red', Text='Extracting Bruker image data from files.',/NOCANCEL)
    progressbar -> Start

    ;this is the point where we open the data and save it to
    ;the data structure.
    adim_counter = 0
    for current_scan =0 , num_scans-1 do begin
        ;update_status_bar, 'opening data'+strtrim(current_scan,2)+'/'+strtrim(num_scans,2)

       ;now that we have the dimesion we can start extracting the data
        OPENR, file_descriptor, multi_image_directory[current_scan]+'2dseq' , /GET_LUN , ERROR = err
        ; If err is nonzero, something happened. Print the error message to
        ; the standard error file (logical unit -2):
        IF (err NE 0) then  begin
           ;PRINTF, -2, !ERROR_STATE.MSG
           update_status_bar,'Unable to open the d3proc file'
           void = dialog_message(['Unable to open the d3proc file',multi_image_directory[current_scan]+'2dseq'],/error)
           return
        end

        ;read in the data
        data[*,*,*] = READ_BINARY (file_descriptor, DATA_DIMS = dimensions , ENDIAN='big' , DATA_TYPE=word_size)

       close, file_descriptor
        FREE_LUN, file_descriptor

        ;after the data is opened copy it to the global structure
        for jj=0 , (dimensions[2]/sdim) -1 do begin
            (*project.dataArray[project.ni].state1)[*,*,*,adim_counter] = data[*,*,jj*sdim:jj*sdim+sdim-1]
            adim_counter ++
            progressBar -> Update, (float(adim_counter)/float(adim-1))*100.0
        end

    endfor

    update_status_bar, ''
    progressbar -> Destroy

    project.imndArray[project.ni].multi_scan_file_array = ptr_new( FILE_SEARCH(multi_scan_dir+'*') )
    project.imndArray[project.ni].sdim = sdim
    project.imndArray[project.ni].adim = adim
    project.procPramArray[project.ni].state_1 = 1
    project.imndArray[project.ni].image_type = 6
    project.imndArray[project.ni].n_bvals =  adim / num_scans

    ;now that scan is open we can set the scan_open_flag
    project.scan_open_flag = 1


    ;since the file opened and was read completely change the current index to the
    ;recently opened file. Then increment the next index
    project.ci = project.ni
    project.ni = project.ni + 1


    cd, image_directory


    void = dialog_message('Unable to extract b-matrix using Bruker multi image')
    ;openADTMulti, multi_scan_dir

    mas_add_current_scan_to_scan_selection

    mas_redraw_GUI

    update_status_bar, ''

    print, project.imndArray[project.ci].adim

    if project.auto_Display eq 1 then mas_display

    HEAP_GC
end

; Subroutine name: chk_T2_single_scan_and_set
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This pro will check to see if there are n_echo_images in the imnd file and if so
; will set the image_type causing the scan to be treated differently when opened.
; This will also set the appropriate array to TE for T2 curve fitting.

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting
    ; Edited by Magdoom 2015/04/12
    ; Included PVM_* search keywords
    
pro chk_T2_single_scan_and_set, imnd
COMPILE_OPT IDL2
COMMON scan_data
                                ;i think this is a good idea to calling the index that we are interested in
                                ;pi for project index this way we can change it more easily
pi = project.ni

;print, 'in T2 check and set'

                                ;some times 3d data has some echos and
                                ;needs to be thrown out and not considered as 3d data.
if (fix(extract_simple_value (imnd ,'##$IMND_dim=', '$'))) eq 3 then return

;print,' after dim check'


                                ;this is the nuber of echos in the
                                ;scan. if the scan has any echos then we will treat it as a
                                ;t2 scan else quit out of this pro.
                                ;n_echos = fix(extract_simple_value (imnd ,'##$IMND_n_echoes=', '$'))
n_echos = fix(extract_simple_value (imnd ,'##$IMND_n_echo_images=', '$'))
if n_echos eq 0 then n_echos = fix(extract_simple_value (imnd ,'##$PVM_NEchoImages=', '$'))
;print, 'n_echos',n_echos

                                ;if there are no echos then exit pro
if n_echos LE 1.0 then return

                                ;now that we are sure how many echos
                                ;there are we can set the adim and the echo time pointer
project.imndArray[pi].adim = n_echos
project.imndArray[pi].image_type = 12

TE = float(extract_simple_value (imnd ,'##$IMND_echo_time=', '$'))
if TE eq 0 then TE = float(extract_simple_value (imnd ,'##$PVM_EchoTime=', '#'))
esp = float(extract_simple_value (imnd ,'##$EchoSpacing=', '$'))
if esp eq 0 then esp = TE
;if TE eq 0 then TE = float(strsplit(extract_simple_value (imnd, '##$EffectiveTE=( '+strtrim(n_echos,2)+' )', '$'),/EXTRACT))

p_echo_time = ptr_new(fltarr(n_echos))
for ii=1, n_echos do begin
    (*p_echo_time)[ii-1] = TE + esp*(ii-1)
end

project.imndArray[pi].echo_Time_ptr = p_echo_time

;print, 'out T2 check and set'

end



; Subroutine name: chk_single_movie_and_set
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This pro will check to see if ##$IMND_movie=Yes and if yes is there then
; image_type will be set to 13.

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting


pro chk_single_movie_and_set, imnd
    COMPILE_OPT IDL2
    COMMON scan_data
    ;I think this is a good idea to calling the index that we are interested in
    ;pi for project index this way we can change it more easily
    pi = project.ni

    ;This is the number of frames in the scan. If the scan has any frames then we will treat it as a
    ;movie scan, else quit out of this pro.
    n_frames = fix(extract_simple_value (imnd ,'##$IMND_n_movie_frames=', '$'))


    ;If there are no echos then exit pro
    if n_frames eq 0 then return

    ;Now that we are sure how many echos there are we can set the adim and the echo time pointer
    project.imndArray[pi].adim = n_frames
    project.imndArray[pi].image_type = 13

end


; Subroutine name: chk_ADT_and_set
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This will check an imnd file to see if it has diffusion experiments and try and extract the b-matrix

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting

pro chk_ADT_and_set, imnd
    compile_opt idl2
    common scan_data

    ;i think this is a good idea to calling the index that we are interested in
    ;pi for project index this way we can change it more easily
    pi = project.ni
    adim = project.imndArray[pi].adim


    if project.imndArray[pi].image_type ne 3 then begin
       ;;project.imndArray[pi].b_matrix = ptr_new()
       return
    end

; HS 20061214
; Lets check whether the file to be searched contains the IMND term.
imnd_check = Strpos(imnd, 'IMND_')

; Different algorithms to extract the b matrix for METHOD or IMND files.
if imnd_check eq -1 then begin

	n_avg = extract_simple_value(imnd, '##$PVM_NAverages=', '#')

	; We have a METHOD file

	; HS - 20061213
	; If this matrix is empty then we have a METHOD file and this b-matrix
	; is reported differently.

	; This is how the b matrix is reported:

	; ======= Begin METHOD file example ==============
	;##$PVM_DwDir=( 1, 3 )
	;0.000000 0.000000 1.000000

	;##$PVM_DwBMat=( 1, 3, 3 )
	;30.065724 0.000000 23.073018 0.000000 0.000000 0.000000 23.073018 0.000000
	;102.862187
	;##$PVM_DwEffBval=( 1 )
	;132.927911
	; ======= End of METHOD file example =============

	b_matrix = extract_b_matrix_method_files(imnd)




endif else begin

	n_avg = extract_simple_value(imnd, 'IMND_n_averages=', '$')

 ;find the b-matrix
    search_strings =['IMND_diff_brr' $
                   ,'IMND_diff_bpp' $
                   ,'IMND_diff_bss' $
                   ,'IMND_diff_brp' $
                   ,'IMND_diff_brs' $
                   ,'IMND_diff_bps' ]


    b_matrix = fltarr(6,adim)

    for ii=0, 5 do begin

       b_matrix[ii,*] = rotate(extract_slice_scheme (imnd, search_strings[ii], adim ),1)
       print, b_matrix


					; Old code left by Ty
					;
					;     imnd_index = 0
					;     while imnd_index lt sz_imnd do begin
					;         flag = strpos(imnd[imnd_index], search_strings[jj])
					;         if (flag ne -1) then break
					;         imnd_index ++
					;     end

					;     n_bvals = extract_simple_value( imnd[imnd_index] ,search_strings[jj] ,')' )

					;     imnd_index++
					;     if imnd_index lt sz_imnd then begin
					;
					;      one_b_matrix_string[jj] = imnd[imnd_index]
					;     end

    endfor

endelse ; closes the METHOD vs. IMND comparison IF statement



;   ;        print, 'number of b_vals', n_bvals
;
;   ;now that we have the one_b_matrix we need make it right.
;   ; it needs to be rotated and converted to floats and mabey some other things.
;   one_b_matrix_float = fltarr(6 , n_bvals)
;
;   for jj=0,5 do begin
;
;     one_b_matrix_float[jj,*] = float(rotate(STRSPLIT(one_b_matrix_string[jj], /EXTRACT),1))
;
;   endfor
;
;
;
;   ;print, one_b_matrix_float
;   b_matrix[*,ii*n_bvals:(ii+1)*n_bvals-1]  = one_b_matrix_float
;
;

; ************************************************************************
; HS 20061213
; There is a problem with reporting the b matrix in .RAW files.
;
; ************************************************************************
;
    ;multiply the off diagionals by 2
    b_matrix[3:5,*] *= 2
;
;

	;store the pointer to the number of averages to each scan.
    ;print, 'n_avg',n_avg
    project.imndArray[project.ci].pn_avg = ptr_new(n_avg)

	; OLD CODE
	;       flag = strpos(imnd[imnd_index], '##$IMND_n_averages=')
	;        if (flag ne -1) then break
	;        imnd_index ++
	;       end
	;   ;once we have located where the string is get the value
	;   n_avg[ii]= float((strsplit(imnd[imnd_index],'=',/extract))[1])

	;   imnd_file_string = ''

    project.imndArray[pi].b_matrix = ptr_new(b_matrix)


end


; Subroutine name: open_multi_slice_movie
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Opens and sets up multi slice movies from bruker files.
; Orders the file listing so that slice position is from lowest to greatest.
; so when the scan data is extracted the files are read in properly.
; set image_type to 14

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting

pro open_multi_slice_movie
    COMPILE_OPT IDL2
    COMMON scan_data

    ;get the proper slash depending upon which operation system is in use.
    slash = get_dir_slash()


    ;check to make sure that max_num_scans was not exceeded
    if project.ni eq project.max_num_scans then begin
       update_status_bar, ' Reached max limit of scans able to open.'
       return
    end


    ;Have the user select a directory that they want to try and open that scan
    parent_directory = DIALOG_PICKFILE( Title = 'Select parent directory of the movies to open' ,$
                 path=project.current_path,DIRECTORY=1)

    ;If the directory returned is empty, then they closed the open without a selection then just return
    if (parent_directory eq '') then return

    ;Get a listing of the files and directories that are in the parent directory
    parent_file_listing = FILE_SEARCH(parent_directory+'*',  COUN=n_parent_dir)

    ;Next we need to extract the slice position so that the we can change the multi_slice_dir
    ;to the proper read in order of FIDs.
    ;First we will open each scan individually and extract the slice_offset.
    ;Then we will sort the slice offset and that will tell us how to read in the fid movies files
    ;for each slice.
    slice_offsets = fltarr( n_parent_dir )
    for ii=0, n_parent_dir-1 do begin

        imnd_file_string = ''
        imnd_temp_string = ''
        OPENR, imnd_file_descriptor, parent_file_listing[ii]+slash[1]+'imnd' , /GET_LUN , ERROR = err
        ; If err is nonzero, something happened. Print the error message to
        ; the standard error file (logical unit -2):
        IF (err NE 0) then  begin
           ;PRINTF, -2, !ERROR_STATE.MSG
           update_status_bar,'Please open the directory containing the imnd files'
           return
        end
        ;If there was no error then the file's text can be extracted.
        While (eof(imnd_file_descriptor) NE 1) Do Begin
           Readf, imnd_file_descriptor, imnd_temp_string
            imnd_file_string = imnd_file_string + imnd_temp_string
        EndWhile
        close, imnd_file_descriptor
        FREE_LUN, imnd_file_descriptor

        slice_offsets[ii] = float(extract_simple_value (imnd_file_string ,'##$IMND_slice_offset=', '$'))

    end

    ;Sort the slice offsets and use that as the key to the parent_file_listing.
    project.imndArray[project.ni].multi_scan_file_array = ptr_new(parent_file_listing[SORT(slice_offsets)])

    ;Next we need to save the reciently extracted string to the project structure.
    project.imndArray[project.ni].imnd_File = imnd_file_string


    ;Once the imnd string has been extracted then store where we found it.
    ;Then update the path so that it contains the last file open.
    ;That way when the next file is open it will be close to the last one.
    ;Locality of reference
    project.imndArray[project.ni].file_Path = parent_directory
    project.current_Path = parent_directory



    ; Try and extract the info from the imnd string. If the extraction of the
    ; imnd file returns 0 then the extraction was bad then display a message to the user
    ; and then return from open.
    if 0 eq mas_extract_imnd ( imnd_file_string ) then begin
       update_status_bar, "Unable to extract 1D data"
       return
    end

    ;After we extract the imnd file we have to overwrite an sdim and adim.
    ;The number of slices is the number of directories in the parent directory
    project.imndArray[project.ni].sdim = n_parent_dir
    ;the number of frames in the movie is the size of adim.
    project.imndArray[project.ni].adim = fix(extract_simple_value (imnd_file_string ,'##$IMND_n_movie_frames=', '$'))
    ;Store what type of open it was so that later on we can update it more easily.
    ;14 is a multi slice movie.
    project.imndArray[project.ni].image_type = 14


    ;Since the file opened and was read completely change the current index to the
    ;recently opened file. Then increment the next index.
    project.ci = project.ni
    project.ni = project.ni + 1

    ;now that scan is open we can set the scan_open_flag
    project.scan_open_flag = 1

    cd, parent_directory

    mas_add_current_scan_to_scan_selection

    mas_redraw_GUI

    update_status_bar, ''

    if project.auto_Display eq 1 then mas_display

    HEAP_GC

end


; Subroutine name: open_dynamics_multi
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This will open the dynamics scans.

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting

pro open_dynamics_multi
    COMPILE_OPT IDL2
    COMMON scan_data

    ;Get the proper slash depending upon which operating system is in use.
    slash = get_dir_slash()

    ;Check to make sure that max_num_scans was not exceeded.
    if project.ni eq project.max_num_scans then begin
       update_status_bar, ' Reached max limit of scans able to open.'
       return
    end


    ;Have the user select a directory that they want to try and open that scan.
    parent_directory = DIALOG_PICKFILE( Title = 'Select directory to open' ,$
                 path=project.current_path,DIRECTORY=1)

    ;If the directory returned is empty, then they closed the open without a selection then just return
    if (parent_directory eq '') then return

    ;Get a listing of the files and directories that are in the parent directory
    parent_file_listing = FILE_SEARCH(parent_directory+'*',  COUN=n_parent_dir)



    ;Once the imnd string has been extracted then store where we found it at.
    ;Then update the path so that it contains the last file open.
    ;That way when the next file is open it will be close to the last one.
    ;Locality of reference
    project.imndArray[project.ni].file_Path = parent_directory
    project.current_Path = parent_directory
    cd, parent_directory

    ;Store the parent dir file listing.
    project.imndArray[project.ni].multi_scan_file_array = ptr_new(parent_file_listing)


    ;Now open the imnd file and take all the data out of it and store it to a string.
    ;If done completely then we are pretty sure that this data is ok so we will increase the
    ;next_index ++
    imnd_file_string = ''
    imnd_temp_string = ''
    OPENR, imnd_file_descriptor, parent_file_listing[0]+slash[1]+'imnd' , /GET_LUN , ERROR = err
    ; If err is nonzero, something happened. Print the error message to
    ; the standard error file (logical unit -2):
    IF (err NE 0) then  begin
       ;PRINTF, -2, !ERROR_STATE.MSG
       update_status_bar,'Please open the directory containing the imnd files'
       return
    end
    ;if there was no error then the file's text can be extracted.
    While (eof(imnd_file_descriptor) NE 1) Do Begin
       Readf, imnd_file_descriptor, imnd_temp_string
        imnd_file_string = imnd_file_string + imnd_temp_string
    EndWhile
    close, imnd_file_descriptor
    FREE_LUN, imnd_file_descriptor

    ;Finally we need to save the recently extracted string to the project structure.
    project.imndArray[project.ni].imnd_File = imnd_file_string

    ;Store what type of open it since we are opening a dynamics scan then set the type to 1
    project.imndArray[project.ni].image_type = 1
    ;print, 'project.imndArray[project.ni].image_type = ',project.imndArray[project.ni].image_type

    ; Try and extract the info from the imnd string. If the extraction of the
    ; imnd file returns 0 then the extraction was bad. Then display a message to the user
    ; and then return from open.
    if 0 eq mas_extract_imnd ( imnd_file_string ) then begin
       update_status_bar, "Unable to extract 1D data"
       return
    end

    ;Since the file opened and was read compleately change the current index to the
    ;recently opened file. Then increment the next index.
    project.ci = project.ni
    project.ni = project.ni + 1

    ;Now that scan is open we can set the scan_open_flag.
    project.scan_open_flag = 1


    mas_add_current_scan_to_scan_selection

    mas_redraw_GUI

    update_status_bar, ''

    if project.auto_Display eq 1 then mas_display

    HEAP_GC

end


; Subroutine name: open_flt_file
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Opens flt files. Will open a file the users has selected and the corresponding geoinfo file.

; Editing Information:
    ;Edited by HS 2006/09/29.
    ;Fix spelling mistakes and commenting
    ;Edited by CD 2007/06/28
    ;Added .DOT file extension, DOT files contain a matrix of DOT probability coefficients (plm's)

pro open_flt_file
    compile_opt idl2
    common scan_data
    ni = project.ni

    ;Check to make sure that max_num_scans was not exceeded
    if project.ni eq project.max_num_scans then begin
        void = dialog_message('reached max limit of scans able to open.', /warning)
        return
    end

    filter = ['*.raw','*.dce','*.T2','*.dwi','*.mpa','*.sps', $
       '*.mps','*.csi','*.T1','*.dcm_dwi','*.smov','*.mmov','*.dot']

    ;Have the user select a directory that they want to try and open that scan.
    flt_file = DIALOG_PICKFILE( FILTER = filter, Title = 'Select *.flt file to open' ,$
                 path=project.current_path)

    ;If the directory returned is empty, then they closed the open without a selection then just return.
    if (flt_file eq '') then return

	; if .dot file, process accordingly
	IF STRPOS(flt_file,'.dot') NE -1 THEN BEGIN

		slh_pos = STRPOS(flt_file,'\',/REVERSE_SEARCH)+1

		pth   = STRMID(flt_file,0,slh_pos)
		froot = STRMID(flt_file,slh_pos)
		froot = STRMID(froot,0,STRLEN(froot)-9)

		real_file = pth + froot + '_real.dot'
		imag_file = pth + froot + '_imag.dot'

		project.dataArray[ni].dotcoeffs = PTR_NEW(COMPLEX(RFLT(real_file), $
			RFLT(imag_file)), /NO_COPY)

		sz_flt = SIZE(*project.dataArray[ni].dotcoeffs)

		project.procPramArray[ni].dot_load_flt = 1

    	project.imndArray[ni].fdim = sz_flt[1]
    	project.imndArray[ni].pdim = sz_flt[2]
    	project.imndArray[ni].sdim = sz_flt[3]

		project.imndArray[ni].slices = sz_flt[3]

		project.imndArray[ni].display_Name = flt_file
		project.imndArray[ni].file_Path = STRMID(flt_file,0,STRPOS(flt_file,'\',/REVERSE_SEARCH)+1)

		trash = INTARR(10,10)
		project.imndArray[ni].b_matrix = ptr_new(trash)

		project.imndArray[ni].image_type = 15
		project.imndarray[ni].state1_load_procedure = 'mas_load_state1_flt'
		
	    project.ci = project.ni
    	project.ni = project.ni + 1
    	project.scan_Open_Flag = 1

	    mas_add_current_scan_to_scan_selection
	    mas_redraw_GUI

		RETURN
	ENDIF

    ;Open the flt file
    project.dataArray[ni].state1 = ptr_new( rflt(flt_file), /no_copy )
    sz_flt = size(*project.dataArray[ni].state1,  /STRUCTURE )

    file_name_length = strpos(flt_file, '.', /REVERSE_SEARCH )
    geo_info_location = strmid( flt_file, 0, file_name_length )+'.GEO_INFO'
    geo_info = get_geo_info_file(geo_info_location)


    project.imndArray[ni].fdim = geo_info.f_dim
    project.imndArray[ni].pdim = geo_info.p_dim
    project.imndArray[ni].sdim = geo_info.s_dim
    project.imndArray[ni].adim = geo_info.a_dim
    project.imndArray[ni].scan_name = geo_info.method
    project.imndArray[ni].f_fov = (geo_info.r_voxel*sz_flt.DIMENSIONS[0]/10)
    project.imndArray[ni].p_fov = (geo_info.p_voxel*sz_flt.DIMENSIONS[1]/10)
    project.imndArray[ni].s_fov = (geo_info.s_voxel*sz_flt.DIMENSIONS[2]/10)
    project.imndArray[ni].display_Name = flt_file
    project.imndArray[ni].n_bvals = geo_info.a_dim

    project.imndArray[ni].file_Path = flt_file
    project.current_Path = flt_file

	; if valid diffusion data add it
	IF PTR_VALID(geo_info.b_matrix) THEN BEGIN
	    ;double the off diagonal
    	(*geo_info.b_matrix)[3:*,*] *= 2
    	project.imndArray[ni].b_matrix = geo_info.b_matrix

		;get bvalue from trace of b_matrix
		bvals = REFORM((*geo_info.b_matrix)[0,*] + (*geo_info.b_matrix)[1,*] + (*geo_info.b_matrix)[2,*])
		project.imndArray[ni].bval_Array = PTR_NEW(bvals)

		project.imndArray[ni].big_delta   = geo_info.big_delta
		project.imndArray[ni].small_delta = geo_info.small_delta
		project.imndArray[ni].angle_theta = geo_info.angle_theta
		project.imndArray[ni].angle_phi   = geo_info.angle_phi
		project.imndArray[ni].pn_avg	  = geo_info.n_avgs
	ENDIF

	IF PTR_VALID(geo_info.echo_times) THEN BEGIN
		project.imndArray[ni].echo_time_ptr = geo_info.echo_times
	ENDIF

	IF PTR_VALID(geo_info.rep_times) THEN BEGIN
		project.imndArray[ni].rep_time_ptr = geo_info.rep_times
	ENDIF

    project.scan_Open_Flag = 1

    project.imndArray[ni].image_type = 15

    project.ci = project.ni

    project.ni = project.ni + 1

    mas_add_current_scan_to_scan_selection

    mas_redraw_GUI

    update_status_bar, ''


end


; Subroutine name: open_raw_file
; Created by: HS - 20061026
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Opens raw files. Will open a file the users has selected and the corresponding geoinfo file.
;
; Flow of subroutine:
;
; Reads the common variables and verifies that the maximum number of scans has not been reached.
; Afterwards, it asks the user to select the .raw file and returns that path. If it is empty,
; it returns without any action. If it contains a path then it creates a integer matrix to read in the
; header file of the .raw file. For more information on the structure of .RAW files check into mas_save.pro
; the subroutine called "save_raw". It then reads the complex data in, attempts to read an IMND string in the
; GEO_INFO file and then some parameters. Sets the flag image_type to zero and this will allow
; mas_load_state_1 to FFT, shift and signal select.
;
; Editing Information:
; Edited on 12212007

pro open_raw_file
    compile_opt idl2
    common scan_data

	ci = project.ci
    ni = project.ni

    ;Check to make sure that max_num_scans was not exceeded
    if project.ni eq project.max_num_scans then begin
        void = dialog_message('Reached max limit of scans able to open.', /warning)
        return
    end

    filter = ['*.raw']

    ;Have the user select a directory that they want to try and open that scan.
    raw_file = DIALOG_PICKFILE( FILTER = filter, Title = 'Select *.raw file to open' ,$
                 path=project.current_path)

    ;If the directory returned is empty, then they closed the open without a selection then just return.
    if (raw_file eq '') then return

		file_name_length = strpos(raw_file, '.', /REVERSE_SEARCH )
    	geo_info_location = strmid( raw_file, 0, file_name_length )+'.GEO_INFO'
		file_test_result = FILE_TEST(geo_info_location)

		IF (file_test_result eq 0) then begin
			print, 'GEO_INFO file is not found'
			print, 'File must be in the same directory as .raw file'
			return
		endif

		; Open the file to get the LUN
		OpenR, lun, raw_file, /Get_LUN, /SWAP_IF_LITTLE_ENDIAN

		; The dimensionality of the data is contained in the first byte of the raw file.
		dimensionality = LonArr(1)

		; Read in the byte number defining the dimensionality of the data
		ReadU, lun, dimensionality

		; Set up the array to determine the dimensionality of the data (up to 5)
		if dimensionality eq 3 then begin
			header_raw = LonArr(6)
			header_temp = LonArr(5)
		endif else begin
			header_raw = LonArr(7)
			header_temp = LonArr(6)
		endelse

		header_raw[0] = dimensionality

		; Read in the rest of the header
		ReadU, lun, header_temp

		header_raw[1:*] = header_temp[*]

		; Create the complex data array depending in the dimensions
		if  dimensionality eq 3 then begin

			complex_data_in = ComplexArr(header_raw[1],header_raw[2], header_raw[3])

		endif else begin

			complex_data_in = ComplexArr(header_raw[1],header_raw[2], header_raw[3], header_raw[4])

		endelse

		; Finish reading the data.
		ReadU, lun, complex_data_in
		free_lun, lun

		project.dataArray[ni].state1 = ptr_new( complex_data_in, /no_copy )

		; HS 20071221
		; going to read the IMND file (or method) attached to the end of a GEO_INFO File.
		; dwi series have multiple IMND strings so I must parse the first one only.


		; Read the IMND string to memory.

	parameter_file_pointer = geo_info_location
	parameter_file_string = ''
    parameter_temp_string = ''
    OPENR, parameter_file_descriptor, parameter_file_pointer , /GET_LUN , ERROR = err

	While (eof(parameter_file_descriptor) NE 1) Do Begin
       Readf, parameter_file_descriptor, parameter_temp_string
        parameter_file_string = parameter_file_string + parameter_temp_string
    EndWhile
    close, parameter_file_descriptor
    FREE_LUN, parameter_file_descriptor

		searchstring ='----------------------'
		position_before_imnd = Strpos(parameter_file_string, searchstring)

	; only going to parse the IMND if it's on the GEO info file
	if (position_before_imnd ne -1) then begin

		position_after_imnd = Strpos( parameter_file_string, searchstring, position_before_imnd+strlen(searchstring))
		difference = position_after_imnd-position_before_imnd

		imnd_file_string = strmid(parameter_file_string, position_before_imnd, difference)

		if 0 eq mas_extract_imnd ( imnd_file_string ) then begin
       		update_status_bar, "Unable to extract IMND parameters"
       		PRINT, 'GEO_INFO file might not have the IMND attached'
       		return
    	end

	endif

; This part must come after the reading of the IMND because the the adim is not being read in correctly
; from the imnd part and would overwrite the adim in raw dwi series. therefore, the header is necessary to
; modify project.imndArray[ni].adim (in line 'project.imndArray[ni].adim = header_raw[4]')


	    print, geo_info_location
	    geo_info = get_geo_info_file(geo_info_location)


		; As far as I know, the dimensions of a scan are stored on the first byte of
		; the header. Substract one to get the correct number. e.g. a 2D scan (RARE) will
		; have a 3 in the zero index. A 2D scan like a variable TE will have a 4 in the zero index.

		if (header_raw[0] eq 2) OR (header_raw[0] eq 4) then begin
			project.imndArray[ni].dimensions = 2
		endif else begin
			project.imndArray[ni].dimensions = header_raw[0]
		endelse

		project.imndArray[ni].fdim = geo_info.f_dim
	    project.imndArray[ni].pdim = geo_info.p_dim
	    project.imndArray[ni].sdim = geo_info.s_dim

	    ;if  dimensionality eq 3 then begin

	    ;endif else  begin
	    	project.imndArray[ni].adim = header_raw[4]
		;endelse

	    project.imndArray[ni].scan_name = geo_info.method

	    project.imndArray[ni].f_fov = (geo_info.r_voxel*header_raw[1]/10)
	    project.imndArray[ni].p_fov = (geo_info.p_voxel*header_raw[2]/10)
	    project.imndArray[ni].s_fov = (geo_info.s_voxel*header_raw[3]/10)

	    project.imndArray[ni].display_Name = raw_file

	    project.imndArray[ni].file_Path = raw_file
	    project.current_Path = raw_file

	    project.imndArray[ni].n_bvals = geo_info.a_dim

		if ptr_valid(b_matrix) then begin
			b_val = *geo_info.b_matrix
		    ;double the off diagonal
	    	b_val[3:*,*] *= 2
	    	project.imndArray[ni].b_matrix = ptr_new(b_val)

			; by this point I truly hope this is a dwi image. if not, bad things will happen
			project.imndArray[ni].image_type = 3


		endif

		; This flag is so that the FFT, shift and signal selection can function on the load_mas_state_1
		project.imndArray[ni].image_type = 0
		project.imndArray[ni].imnd_file = geo_info.param_text

		;project.big_endian = 1


		project.ci = project.ni
	    project.ni = project.ni + 1
	    project.scan_Open_Flag = 1

		mas_add_current_scan_to_scan_selection
	    mas_redraw_GUI
	    update_status_bar, ''

end



;*****************************************************************************************************
;+
; NAME:
;   mas_open
;
; PURPOSE:
;
;   Open and load a scan into MAS
;
; ARGUMENTS:
;
;   open_type [in]    integer that represents the different types of scans to open.[0|1|2]
;
; CALLED BY:
; mas.pro
;
; FLOW:
;
;
; MODIFICATION HISTORY:
;
;- HS 2006/09/29. Fixed spelling mistakes and comment structure.
;*****************************************************************************************************


PRO mas_open, open_type, scan_directory=scan_directory

    COMPILE_OPT IDL2

    COMMON scan_data

    print, "mas_open: arrived."

    ;Check to make sure that max_num_scans was not exceeded
    if project.ni eq project.max_num_scans then begin
        void = dialog_message('Maximum number of scans reached:'+strcompress(string(project.max_num_scans)), $
                              /center, /error)
        update_status_bar, ' Reached max limit of scans able to open.'
        return
    end

    if (n_elements(scan_directory) eq 0) then begin
        ;Have the user select a directory that they want to try and open that scan.
        scan_directory = DIALOG_PICKFILE( title = 'Select directory of the scan to open' ,$
                                          path=project.current_path, $
                                          DIRECTORY=1)
    endif
    
    ;If the directory returned is empty, then they closed the open without a selection then just return
    if (scan_directory eq '') then return

    ; Check to see if the directory contains a fid file.
    ; If it does then the directory is valid.
    ; If not then then send the user a message and then return from open.
    OPENR, descriptor, scan_directory + 'fid' , /GET_LUN , ERROR = err
                                ; If err is nonzero, something
                                ; happened. send information to
                                ; user
                                ; or maybe this is a ser file instead of a fid file.
    IF (err NE 0) then  begin
        OPENR, descriptor, scan_directory + 'ser' , /GET_LUN , ERROR = err
        IF (err NE 0) then  begin
                                ;PRINTF, -2, !ERROR_STATE.MSG
            update_status_bar,'Please open a directory containing a FID file or SER'
            return
        end
    end
    close, descriptor
    FREE_LUN, descriptor


	;Verify is we have an IMND or METHOD file

	parameter_file_pointer = scan_directory + 'imnd'

	IMND_Exists = FILE_SEARCH(parameter_file_pointer)
	IMND_size = size(IMND_Exists)

	if IMND_size[0] eq 0 then begin

		parameter_file_pointer = scan_directory + 'method'
		METHOD_Exists = FILE_SEARCH(parameter_file_pointer)
		METHOD_size = size(METHOD_Exists)
    
    if METHOD_size[0] eq 0 then begin
      parameter_file_pointer = scan_directory + 'acqu'
      ACQU_Exists = FILE_SEARCH(parameter_file_pointer)
      ACQU_size = size(ACQU_Exists)
      if ACQU_size[0] eq 0 then begin
        update_status_bar,'Please open the directory containing the imnd or method files'
        return
      endif
    endif
	endif


    ;Now open the parameter (IMND or METHOD) file and take all the data out of it and store it to a string.
    ;If done completly then we are pretty sure that this data is ok so we will increase the
    ;next_index ++
    parameter_file_string = ''
    parameter_temp_string = ''
    OPENR, parameter_file_descriptor, parameter_file_pointer , /GET_LUN , ERROR = err
    ; If err is nonzero, something happened. Print the error message to
    ; the standard error file (logical unit -2):
    ; if there was an error then look for a method file
		;    IF (err NE 0) then  begin
		;
		;       OPENR, parameter_file_descriptor, scan_directory + 'method' , /GET_LUN , ERROR = err2
		;       IF (err2 NE 0) then  begin
		;           update_status_bar,'Please open the directory containing the imnd or method files'
		;           return
		;       end
    	;	 end

    ;If there was no error then the file's text can be extracted.
    While (eof(parameter_file_descriptor) NE 1) Do Begin
       Readf, parameter_file_descriptor, parameter_temp_string
        parameter_file_string = parameter_file_string + parameter_temp_string
    EndWhile
    close, parameter_file_descriptor
    FREE_LUN, parameter_file_descriptor

    ;Finally we need to save the recently extracted string to the project structure.
    project.imndArray[project.ni].imnd_File = parameter_file_string

    ;Store what type of open it was so that later on we can update it if neccessary.
    project.imndArray[project.ni].image_type = open_type
    ;print, 'project.imndArray[project.ni].image_type = ',project.imndArray[project.ni].image_type

    ;; tell mas how to load the state 1 data
    project.imndArray[project.ni].state1_load_procedure = 'mas_load_state_1_bruker'
    
    ;Once the imnd string has been extracted then store where we found it at.
    ;Then update the path so that it contains the last file open.
    ;That way when the next file is open it will be close to the last one.
    ;Locality of reference
    project.imndArray[project.ni].file_Path = scan_directory
    project.current_Path = scan_directory
    
	; Now we will use the appropriate extracting function depending on what file type we have.
	if IMND_size[0] eq 1 then begin
		if 0 eq mas_extract_imnd ( parameter_file_string ) then begin
       		update_status_bar, "Unable to extract 1D data"
       		return
    	end
	endif else begin
	 if isa(ACQU_size) then if ACQU_size[0] eq 1 then project.imndArray[project.ni].dimensions = 1
	 if 0 eq mas_extract_method(parameter_file_string) then begin
       update_status_bar, "Unable to extract 1D data"
       return
   end
	endelse

    ;Setting the T1 image_type pointer
    ;If the pointer to multiple tr's is valid  and type is 0 then the data must be T1
    if PTR_VALID(project.imndArray[project.ni].rep_Time_ptr) and project.imndArray[project.ni].image_type eq 0 then begin
        project.imndArray[project.ni].image_type = 9
    end
    ;print, 'rep time ptr image_type = ',project.imndArray[project.ni].image_type

    ;since we have the imnd string we should see if the image is T2 single scan type
    chk_T2_single_scan_and_set, parameter_file_string
    ;print, 'after T2 set image_type = ',project.imndArray[project.ni].image_type

    chk_single_movie_and_set, parameter_file_string
    ;print, 'after single movie image_type = ',project.imndArray[project.ni].image_type

    ;check to see if the single scan is an ADT scan and then extract the b-matrix and set the image type
    chk_ADT_and_set, parameter_file_string

    ;since the file opened and was read compleately change the current index to the
    ;recently opened file. Then increment the next index.
    project.ci = project.ni
    project.ni = project.ni + 1

    ;now that scan is open we can set the scan_open_flag
    project.scan_open_flag = 1

    cd, scan_directory

    mas_add_current_scan_to_scan_selection

    turnoff_motion_correction

    mas_redraw_GUI

    update_status_bar, ''

    if project.auto_Display eq 1 then mas_display

    HEAP_GC

END
