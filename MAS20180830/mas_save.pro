;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved

; Editing Information:
; HS 2006/09/30
; Renamed the get_file_extenshion function to get_file_extension and
; all associated calls within this .PRO file.

function ok_to_write, fs_item

    if (file_test(fs_item, /REGULAR) or file_test(fs_item, /DIRECTORY)) then begin

        if (file_test(fs_item, /WRITE)) then return, 1

    endif else begin

        parent = file_dirname(fs_item)

        if (file_test(parent, /WRITE)) then return, 1

    endelse

    junk = display_message('Unable to write to '+fs_item)

    return, 0

end

; Subroutine name: save_subset
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/09/30.
    ;Fix spelling mistakes and commenting

pro save_subset
    COMPILE_OPT IDL2

    COMMON scan_data, project
    CI = project.ci

    ;if there is no open scan then don't do any thing
    if project.scan_Open_Flag eq 0 then return

    image_type = project.imndArray[ci].image_type

    file = (Dialog_PickFile( Title = 'Save current scan to flt', path=project.current_path))

    if file eq '' then return
    ;print, file
    mas_load_state_1

    ;get the file extension that goes at the end of the file name
    ext = get_file_extension( image_type )

    ;update_status_bar,'Saving .flt file...'
    sz_data = size(*project.dataArray[CI].state1)

    sig = *project.dataArray[CI].state1

    ;use the roi mask
    IF project.roi.mask gt 0 THEN $
       for jj=0, sz_data[4]-1 do $
       FOR ii=0, sz_data[3]-1 DO BEGIN
         print, ii,jj
         p_image = ptr_new(sig[*,*,ii,jj])
         mas_roi_mask, p_image
         sig[*,*,ii,jj] = *p_image
    END

    OpenW, lun0, file+ext, /SWAP_IF_LITTLE_ENDIAN, /get_lun
    WriteU, lun0, sz_data, sig
    Free_LUN, lun0

    write_geo_info_file, file, sz_data

    update_status_bar,''

end

; ================== OLD SAVE_RAW SUBROUTINE =========================
; Subroutine name: save_raw
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This pro will allow the user to save raw data to a .sav file.

; Editing Information:
    ;Edited by HS 2006/09/30.
    ;Fix spelling mistakes and commenting

;pro save_raw
;    COMPILE_OPT IDL2
;
;    COMMON scan_data, project
;    CI = project.ci
;
;    ;check to see if the data has a k-spce
;    image_type = project.imndArray[CI].image_type
;    if image_type eq 11 or image_type eq 5 or image_type eq 6 then begin
;        void = dialog_message('The type of image selected does not have a native K-space file.',/error)
;        return
;    end
;
;    save_path = (Dialog_PickFile( Title = 'Save Current K-space data to', path=project.current_path))
;    if save_path eq '' then return
;
;    raw = *(mas_extract_data())
;    save, raw ,FILENAME=save_path+'.sav'
;
;end
;size array = size () [dim, fdim, pdim, sdim, adim, type, size]

; ======================= END OF OLD SAVE_RAW SUBROUTINE =============



; Subroutine name: save_raw
; Created by: HS
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This pro will allow the user to save raw data to an .flt file. This new subroutine
; will be used to provide users with the raw data if they want to process it on another
; programming language. The output will be a file that has 4 header bytes with information
; about the data and then the data serialized in floating point format.

; Editing Information:
;  BT, 2008-03-12 - Will now apply motion correction where applicable
;                   before saving .raw

pro save_raw
    COMPILE_OPT IDL2

    COMMON scan_data, project
    CI = project.ci

    ;;check to see if the data has a k-space
    image_type = project.imndArray[CI].image_type

    if image_type eq 11 or image_type eq 5 or image_type eq 6 then begin
        void = dialog_message('The type of image selected does not have '+ $
                              'a native K-space file.',/error)
        return
    end
    
    save_path = (Dialog_PickFile( Title = 'Save current K-space data to:', $
                                  path=project.current_path))
    extension = +'.raw'
    if save_path eq '' then return

    ;; This next call is what extracts the data to this "raw" variable.
    raw = mas_extract_data()

    ;; don't free this pointer until you're sure it doesn't point to
    ;; state1 anymore
    tmp = ptr_new()

    if project.procPramArray[CI].mc_enable eq 2 then begin
 
       ;; this is a hackish for now, but ok.
        
        ;; save state1 pointer
        tmp = project.dataArray[ci].state1
        
        ;; point state1 to the newly loaded raw data
        project.dataArray[ci].state1 = raw
        
        ;; motion correct the raw data
        mas_motion_correct

        ;; point state1 back to its original data
        project.dataArray[ci].state1 = tmp
        

    endif

    ;; We need to get the size of the incoming matrix to write the header of the created
    ;; .RAW file.
    sz_data = size(*raw)

    ;; Now we write the header and the data to the desired path.
    OpenW, lun0, save_path+extension, /SWAP_IF_LITTLE_ENDIAN, /get_lun
    WriteU, lun0, sz_data, *raw
    Free_LUN, lun0
    
    ;; Write the GEO_INFO file which contains the minimum amount of information to
    ;; read in and display the data.
    write_geo_info_file, save_path, sz_data
    
    ptr_free, raw

    HEAP_GC

    update_status_bar,''
    
end


; Subroutine name: write_geo_info_file
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This subroutine writes a text file with necessary information for opening .RAW files.

; Editing Information:
; Edited by HS 2006/09/30.
; Fix spelling mistakes and commenting
; HS - 20061028
; Now this function includes the IMND file when it creates the GEO_INFO file
; per request by Dr. Mareci on 20061027.


PRO write_geo_info_file, file, size_array
    COMPILE_OPT IDL2

    COMMON scan_data, project
    CI = project.ci

	size_of_pn_avg = SIZE((*project.imndArray[ci].pn_avg),/dimension)
	array_dim = project.imndArray[ci].adim

	IF size_of_pn_avg NE 0 THEN n_repeats = array_dim/size_of_pn_avg ELSE n_repeats = array_dim

    OPENW, lun0, file+'.GEO_INFO', /SWAP_IF_LITTLE_ENDIAN, /GET_LUN

    PRINTF, lun0, 'filename='+file
    PRINTF, lun0, 'f_dim='+STRTRIM(size_array[1],2)
    PRINTF, lun0, 'p_dim='+STRTRIM(size_array[2],2)
    PRINTF, lun0, 's_dim='+STRTRIM(size_array[3],2)
    PRINTF, lun0, format='("r_voxel=",E11.4 )',project.imndArray[ci].f_fov*10/size_array[1]
    PRINTF, lun0, format='("p_voxel=",E11.4 )',project.imndArray[ci].p_fov*10/size_array[2]
    PRINTF, lun0, format='("s_voxel=",E11.4 )',project.imndArray[ci].s_fov*10/size_array[3]
    PRINTF, lun0, format='("image_plane_psi=",E11.4 )',0.0
    PRINTF, lun0, format='("image_plane_theta=",E11.4 )',0.0
    PRINTF, lun0, format='("image_plane_phi=",E11.4 )',0.0
    PRINTF, lun0, 'method='+project.imndArray[ci].scan_name

	if ptr_valid(project.imndArray[CI].big_delta) then temp = STRTRIM((*project.imndArray[CI].big_delta)[0],2) else temp = 'N/A'
    ;PRINTF, lun0, 'big_delta='+STRTRIM((*project.imndArray[CI].big_delta)[0],2)
	PRINTF, lun0, 'big_delta='+ temp

	if ptr_valid(project.imndArray[CI].small_delta) then temp = STRTRIM((*project.imndArray[CI].small_delta)[0],2) else temp = 'N/A'
	PRINTF, lun0, 'small_delta=' + temp

    PRINTF, lun0, 'a_dim='+STRTRIM(project.imndArray[ci].adim,2)

	; HS - 20061217
	;When we report the b-matrix we had already scaled it by a factor of 2. Therefore, we need to
	;unscale them so that when they are read back in is correct. Check the subroutine mas_open_adt in mas_open.sav
	;for more information

   ;======================= DIFFUSION SCAN================
   IF PTR_VALID(project.imndArray[ci].b_matrix) THEN BEGIN

		PRINTF, lun0, 'Diffusion Weighting Values: '

        FOR ii=0,  project.imndArray[ci].adim-1 DO BEGIN
            PRINTF, lun0, format='(I,",",E12.4,",",E12.4,",",E12.4,",",E12.4,",",E12.4,",",E12.4,",",E12.4,",",E12.4 )' $
                ,(*project.imndArray[ci].pn_avg)[ii/n_repeats] $
                ,(*project.imndArray[ci].angle_theta)[ii/n_repeats] $
                ,(*project.imndArray[ci].angle_phi)[ii/n_repeats] $
                ,(*project.imndArray[ci].b_matrix)[0,ii] $
                ,(*project.imndArray[ci].b_matrix)[1,ii] $
                ,(*project.imndArray[ci].b_matrix)[2,ii] $
                ,(*project.imndArray[ci].b_matrix)[3,ii]/2 $
                ,(*project.imndArray[ci].b_matrix)[4,ii]/2 $
                ,(*project.imndArray[ci].b_matrix)[5,ii]/2
       	ENDFOR
	ENDIF

	;======================= VAR_TE SCAN========================
	IF PTR_VALID(project.imndArray[ci].echo_Time_ptr) THEN BEGIN
		PRINTF, lun0, 'Multiple Echo Times: '
		PRINTF, lun0, (*project.imndArray[ci].echo_Time_ptr)
	ENDIF

	;======================= VAR_TR SCAN=======================
	IF PTR_VALID(project.imndArray[ci].rep_Time_ptr) THEN BEGIN
		PRINTF, lun0, 'Multiple Repetition Times: '
		PRINTF, lun0, (*project.imndArray[ci].rep_Time_ptr)
	ENDIF

	; if the data is from Philips (DICOM or PAR/REC) this parameter will not be set
	; do not write out IMND data for data acquired on non-Bruker systems
	IF (PTR_VALID(project.imndArray[ci].multi_scan_file_array) NE 1) AND (PTR_VALID(project.imndArray[ci].rep_Time_ptr) EQ 0) THEN BEGIN
		; add PAR file

		textfile_path = project.imndArray[CI].file_Path + '.PAR'

		PAR_Exists = FILE_SEARCH(textfile_path)
		PAR_size = SIZE(PAR_Exists)

		IF PAR_size[0] EQ 0 THEN BEGIN
			update_status_bar,'The directory does not contain a PAR file'
			FREE_LUN, lun0
			RETURN
		ENDIF

		PRINTF, lun0, 'Start of PAR file:'
		PRINTF, lun0, textfile_path
		PRINTF, lun0, '----------------------'

		parameter_temp_string = ''
    	OPENR, text_file_pointer, textfile_path , /GET_LUN

		WHILE (EOF(text_file_pointer) NE 1) DO BEGIN
			READF, text_file_pointer, parameter_temp_string
			PRINTF, lun0, parameter_temp_string
		ENDWHILE

		PRINTF, lun0, '----------------------'

	ENDIF ELSE BEGIN

		; Create a matrix with the array dimensions to store the location of FID files
		location_fids = STRARR(project.imndArray[ci].adim)

		slash =  get_dir_slash()

		; Now generate a matrix with the correct location for all the fids.
		; This will later be used to identify each individual IMND file correctly in the GEO_INFO
		;location_fids = *project.imndArray[ci].multi_scan_file_array

		; only one file to output
		;IF (SIZE(location_fids))[1] EQ 1 THEN BEGIN
		IF (PTR_VALID(project.imndArray[ci].multi_scan_file_array) EQ 0) THEN BEGIN

			PRINTF, lun0, 'Start of IMND or METHOD file of:'
			PRINTF, lun0, project.imndArray[ci].file_Path + 'fid'
			PRINTF, lun0, '----------------------'
			PRINTF, lun0, project.imndArray[ci].imnd_File
			PRINTF, lun0, '----------------------'

		ENDIF ELSE BEGIN
			location_fids = *project.imndArray[ci].multi_scan_file_array

			; Create a matrix with the path to the imnds
			location_imnds = location_fids + slash[1] + 'imnd'

			; Finally attach the correct file to the location of the fids
			location_fids = location_fids +slash[1] + 'fid'
			imnd_temp_string = ''
			imnd_file_string = ''

			; Now it is going to step through each single scan and read in the IMND.
			; This is the same method of reading in the IMND used in the mas_open.pro.
			; Then write the IMND of each DWI scan separated by the FID location.

			FOR ii=0, project.imndArray[ci].adim-1 DO BEGIN
				PRINTF, lun0, 'Start of IMND or METHOD file of:'
				PRINTF, lun0, location_fids[ii]
				PRINTF, lun0, '----------------------'

				; Now we have to verify the existance of the IMND or METHOD file.
				OPENR, imnd_file_descriptor, location_imnds[ii], /GET_LUN , ERROR = err

				IF (err NE 0) THEN BEGIN
				 	scan_directory = *project.imndArray[ci].multi_scan_file_array + slash[1] + 'method'

		    		OPENR, imnd_file_descriptor, scan_directory[ii], /GET_LUN , ERROR = err2

		       		IF (err2 NE 0) THEN BEGIN
		           		update_status_bar,'Please open the directory containing the imnd or method files'
		           		FREE_LUN, lun0
		           		RETURN
		       		ENDIF
			    ENDIF

				; Here we read in all the text of the IMND of METHOD into one variable.

				WHILE (EOF(imnd_file_descriptor) NE 1) DO BEGIN
		       		READF, imnd_file_descriptor, imnd_temp_string
		        	imnd_file_string = imnd_file_string + imnd_temp_string
	    		ENDWHILE

				CLOSE, imnd_file_descriptor
	    		FREE_LUN, imnd_file_descriptor

				PRINTF, lun0, imnd_file_string
				PRINTF, lun0, '----------------------'

				; Clear the IMND string.
				imnd_file_string = ''
			ENDFOR
		ENDELSE
	ENDELSE

	FREE_LUN, lun0
END


; Subroutine name: get_file_extension
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This pro will return the extension depending on a flag.

; Editing Information:
    ;Edited by HS 2006/09/30.
    ;Fix spelling mistakes and commenting
    ;Renamed the function from get_file_extenshion to get_file_extension.

function get_file_extension, type
    if type eq 0 then return,'.raw'
    if type eq 1 then return,'.dce'
    if type eq 2 then return,'.T2'
    if type eq 3 then return,'.dwi'
    if type eq 4 then return,'.mpa'
    if type eq 5 then return,'.sps'
    if type eq 6 then return,'.mps'
    if type eq 7 then return,'.spc'
    if type eq 8 then return,'.csi'
    if type eq 9 then return,'.T1'
    if type eq 10 then return,'.fib'
    if type eq 11 then return,'.dcm_dwi'
    if type eq 12 then return,'.T2'
    if type eq 13 then return,'.smov'
    if type eq 14 then return,'.mmov'
    if type eq 16 then return, '.nii'
    return,-1

end


; Subroutine name: save_flt
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This procedure will save the flt files from the current project
; after they have been FFT'ed

; Editing Information:
    ;Edited by HS 2006/09/30.
    ;Fix spelling mistakes and commenting

pro save_flt
    COMPILE_OPT IDL2

    COMMON scan_data, project
    CI = project.ci

    ;if there is no open scan then don't do any thing
    if project.scan_Open_Flag eq 0 then return

    image_type = project.imndArray[ci].image_type

    file = (Dialog_PickFile( Title = 'Save Current Scan to flt', path=project.current_path))

    if file eq '' then return
    ;print, file
    mas_load_state_1

    pImage = ptr_new(*project.dataArray[CI].state1)
    mas_rotate_flip, pImage

    ;get the file extension that goes at the end of the file name
    ext = get_file_extension( image_type )

    ;update_status_bar,'Saving .flt file...'
    sz_data = size(*pImage)
    OpenW, lun0, file+ext, /SWAP_IF_LITTLE_ENDIAN, /get_lun
    WriteU, lun0, sz_data, *pImage
    Free_LUN, lun0

    write_geo_info_file, file, sz_data

    update_status_bar,''

    ptr_free, pImage

end

;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: save_fid

;; Created 20080207 Bill Triplett

;; Calling Information:



;; Bugs or Important Comments to Developers:

;; Currently only saved single scan files

;; Purpose of subroutine:

;; Saves post-processed time series data to Bruker FID file

;; Editing Information:

PRO save_fid

    COMPILE_OPT IDL2

    COMMON scan_data, project

    ci = project.ci

    image_type = project.imndArray[ci].image_type
    tmp_signal = project.procPramArray[ci].signal_type

    ;; magnitude data is preferred for this operation
    ;;project.procPramArray[ci].signal_type = 0
    ;; reset state_1 data, to guarantee state_1 has magnitude data

    if image_type ne 9 or image_type ne 0 then begin
        void = dialog_message('This image type cannot currently be rewritten as a FID file.', $
                              /ERROR)
        return
    endif

    ;; 1) determine what kind of scan data we're dealing with (Single, DWI ADT
    ;; Dynamics). There seem to be two format types: one is a single scan with
    ;; fid/imnd, and the other is a scan directory with a sequence of
    ;; subdirectories containing imnd/fid pairs.

    ;; 2) if we have a single imnd/fid pair
    ;;         a) then we can just write it out
    ;;    else
    ;;         b) need to make each subdirectory
    ;;            and write and perform 2a).
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    save_path = (dialog_pickfile(Title = 'Save fid file to:', path=project.current_path))
    if save_path eq '' then return

    ;; We want to save post-state1 and pre-state 2
    if project.procPramArray[ci].state_1 eq 0 then begin
        mas_load_state_1
    endif

    export = ptr_new(*project.dataArray[ci].state1)

    ;; data needs to be un-shifted
    mas_shift, export

    fdim = project.imndarray[ci].fdim
    pdim = project.imndarray[ci].pdim
    sdim = project.imndarray[ci].sdim
    adim = project.imndarray[ci].adim
    rare_factor = project.imndarray[ci].rare
    islice = *project.imndarray[ci].slice_scheme_ptr

    ;; will hold the fft^-1'd data
    output = complexarr(fdim, pdim, sdim, adim) ;

    for i = 0, adim-1 do begin
        for j = 0, sdim-1 do begin
            output[*,*,j,i] = fft((*export)[*,*,j,i], /INVERSE)
        endfor
    endfor

    ptr_free, export

    ;; open the output fid (pp = post-processed)
    OPENW, lun, 'fid_pp', ERROR=err, /GET_LUN

    if (err ne 0) then begin
        print, "Error opening file."
        return
    endif

    for oo = 0, adim-1 do begin

        ;; Check data file for size in mod 128
        junk_size = 2 * (fdim MOD 128)
        subdim = pdim/rare_factor
        block = lonarr(2,fdim)

        if junk_size gt 0 then begin
            junk = lonarr(junk_size)
        endif

        For i=0,subdim-1 Do Begin

            For j=0,sdim-1 Do Begin

                For k=0,rare_factor-1 Do Begin

                    ;; split out the re, im parts and write.
                    block[0,*] = real_part(output[*, i+subdim*k, islice[j], oo])
                    block[1,*] = imaginary(output[*, i+subdim*k, islice[j], oo])
                    writeu, lun, block ;; TRANSFER_COUNT=ct

                    ;; pad non-bruker block sized chunks
                    if junk_size gt 0 then begin
                        writeu, lun, junk
                    endif

                EndFor

            EndFor

        Endfor

    endfor

    close, lun

    ;; restore the user settings, state_1 will need to be reloaded
    project.procPramArray[ci].signal_type = tmp_signal
    project.procPramArray[ci].state_1 = 0

    heap_gc

END
;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;


; Subroutine name: save_tiff
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Saves state 2 in a tiff file.

; Editing Information:
    ;Edited by HS 2006/09/30.
    ;Fix spelling mistakes and commenting
	;Edited by HS 20070831
	; Changing the behavior of how the tiff files get generated
	; Only the current slice will be generated unless the user selects "multiple" from
	; the display window.

PRO save_tiff
    COMPILE_OPT IDL2

    COMMON scan_data, project

    mas_load_state_2

    tiffFileLocation = DIALOG_PICKFILE( FILTER = '*.tif', Title = 'Save Current Scan', path=project.current_path)

    if tiffFileLocation eq '' then return

    if (file_test(tiffFileLocation, /directory)) then begin
        project.current_path = tiffFileLocation
    endif else begin
        project.current_path = file_dirname(tiffFileLocation)
    endelse

    ;print, tiffFileLocation
    sizeOfImage = size(*project.dataArray[project.ci].state2 )

    ;setting the correct field of view so that it correctly represents
    ; the rotation of the image
    rotate_Direction = project.procPramArray[project.ci].rotate_Direction
    if rotate_Direction eq 0 or rotate_Direction eq 2 then begin
       x_fov = project.imndArray[project.ci].f_fov
       y_fov = project.imndArray[project.ci].p_fov
    endif else begin
       x_fov = project.imndArray[project.ci].p_fov
       y_fov = project.imndArray[project.ci].f_fov
    end


    XRESOL = sizeOfImage[1]/(x_fov * (1/2.54) )
    YRESOL = sizeOfImage[2]/(y_fov * (1/2.54) )

    ;since we are windowing the data that may not have been what the user wants so we have to reset the
    ;state 2 flag

if (project.procpramarray[project.ci].single_multi_flag eq 0) then begin ; User has 'single' selected
; from the main MAS window

	currentslice = project.procpramarray[project.ci].sdim_start

	if (sizeOfImage[0] eq 2) then begin

		;print, 'inside tiff_write with one slice'
		image = reverse((*project.dataArray[project.ci].state2),2)
	    mas_windowing, image

		WRITE_TIFF,tiffFileLocation + '_' + STRTRIM(currentslice+1,2)+'.tif',image , $
            XRESOL=XRESOL ,YRESOL=YRESOL , orientation = 1

;	endif else if (sizeOfImage[0] eq 3) then begin
;
;		print, 'inside tiff_write with multiple slice'
;		image = reverse((*project.dataArray[project.ci].state2),2)
;	    mas_windowing, image
;
;		currentadim = project.procpramarray[project.ci].adim_start
;
;	    slice = image[*,*,currentslice, currentadim]
;
;		WRITE_TIFF,tiffFileLocation + '_' + STRTRIM(currentslice,2)+'.tif',slice , $
;            XRESOL=XRESOL ,YRESOL=YRESOL , orientation = 1

	endif

endif else if (project.procpramarray[project.ci].single_multi_flag eq 1) then begin ;This would be if the user
; selects 'Multiple' image display from the main MAS window

    if sizeOfImage[0]eq 2 then begin
       ;print,'2dimension'
       ;mas_windowing,*project.dataArray[project.ci].state2
       WRITE_TIFF,tiffFileLocation+'.tif',reverse(*project.dataArray[project.ci].state2,2), $
         XRESOL=XRESOL ,YRESOL=YRESOL , orientation = 1
    end

    if sizeOfImage[0]eq 3 then begin

       progressbar = Obj_New('progressbar', Color='red', Text='Saving tiff files',/NOCANCEL)
       progressbar -> Start

        image = reverse((*project.dataArray[project.ci].state2),2)
        mas_windowing, image

       for ii=0,sizeOfImage[3]-1 do begin
         ;image = reverse((*project.dataArray[project.ci].state2)[*,*,ii],2)

         WRITE_TIFF,tiffFileLocation + '_' + STRTRIM(ii+1,2)+'.tif',image[*,*,ii] , $
            XRESOL=XRESOL ,YRESOL=YRESOL , orientation = 1

         progressBar -> Update, (float(ii)/float(sizeOfImage[3]-1))*100.0
       end

       progressbar -> Destroy

    end

endif 

end

pro save_tiff_series, series, dest=dest

    common scan_data, project
    common common_widgets

    ;; get the destination directory
    if (not keyword_set(dest)) then begin
        dest = dialog_pickfile(/directory, path=project.current_path)
    endif
    
    if (dest eq '') then return
    
    ;; save the sdim to restore later
    prev_sdim = project.procPramArray[project.ci].sdim_start
    
    parent = dest

    for i=0, n_elements(series)-1 do begin

        ;; we want two digit numbers
        if (series[i] lt 10) then begin
            serial = strcompress('0'+string(series[i]), /remove_all)
        endif else begin
            serial = strcompress(string(series[i]), /remove_all)
        endelse

        ;; create the directory
        dirname =  strcompress(parent+'S'+serial+'/')
        if (not file_test(dirname, /directory)) then begin
            print, "making: "+dirname
            file_mkdir, dirname
        endif
        
        ;; increment the sdim slider 
        widget_control, sdim_slider, set_value=series[i]
        project.procPramArray[project.ci].state_2=0

        ;; load up the imagery in state1
        mas_load_state_2

        ;; prepare the image data
        sizeOfImage = size(*project.dataArray[project.ci].state2 )
        
        rotate_Direction = project.procPramArray[project.ci].rotate_Direction
        if rotate_Direction eq 0 or rotate_Direction eq 2 then begin
            x_fov = project.imndArray[project.ci].f_fov
            y_fov = project.imndArray[project.ci].p_fov
        endif else begin
            x_fov = project.imndArray[project.ci].p_fov
            y_fov = project.imndArray[project.ci].f_fov
        endelse

        XRESOL = sizeOfImage[1]/(x_fov * (1/2.54) )
        YRESOL = sizeOfImage[2]/(y_fov * (1/2.54) )

        image = reverse((*project.dataArray[project.ci].state2),2)
        mas_windowing, image

        ;; save the adims as tiffs
        for ii=0,sizeOfImage[3]-1 do begin
            
            tifname = strcompress(dirname+'A'+string(ii+1)+'.tif', /remove_all)
            print, tifname
            WRITE_TIFF, tifname, image[*,*,ii], XRESOL=XRESOL ,YRESOL=YRESOL , orientation = 1
            
        endfor

    end

    widget_control, sdim_slider, set_value=prev_sdim
    project.procPramArray[project.ci].state_2 = 0
    
end


; Subroutine name: save_mpeg
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Saves state 2 in an mpeg file.

; Editing Information:
    ;Edited by HS 2006/09/30.
    ;Fix spelling mistakes and commenting


PRO save_mpeg
    COMPILE_OPT IDL2

    COMMON scan_data, project

    update_status_bar,'MPEG exporter'


    mas_load_state_2

    ;find out how big the image or array of images is
    mpgSize = size(*project.dataArray[project.ci].state2)

    ;if the image is not 3 dimensions then don't make a movie
    if( not ( mpgSize[0] eq 3) ) then begin
       update_status_bar,'not enough images to make a movie'
       return
    end

    ;select a place to save the data.
    mpegFileLocation = DIALOG_PICKFILE( FILTER = '*.mpg', Title = 'Save Current Scan')

    if mpegFileLocation eq '' then return

    mpegFileLocation = mpegFileLocation+'.mpg'


    min_intensity = min(*(project.dataArray[project.ci].state2), max=max_intensity)

    dim = mpgSize[1:2]
    ;print,dim

    ; Open an MPEG sequence:
    mpegID = MPEG_OPEN( dim, FILENAME=mpegFileLocation, QUALITY=100)

    ; Add the frames:
    ii=0
    for ii=0,mpgSize[3]-1 do begin
       image = (*project.dataArray[project.ci].state2)[*,*,ii]
       mas_windowing, image, max_intensity, min_intensity
       MPEG_PUT, mpegID, IMAGE=image, FRAME=ii
       update_status_bar,string('Loading frame ',strTrim(ii,2),'/',strTrim(mpgSize[3]-1,2) )
    end

    update_status_bar, string('Saving:', mpegFileLocation)
    ; Save the MPEG sequence in the file myMovie.mpg:
    MPEG_SAVE, mpegID, FILENAME=mpegFileLocation

    ; Close the MPEG sequence:
    MPEG_CLOSE, mpegID

    update_status_bar,''

end


; Subroutine name: MAS_MPEG_WIDGET_EVENT
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Saves state 2 in an mpeg file.

; Editing Information:
    ;Edited by HS 2006/09/30.
    ;Fix spelling mistakes and commenting

pro MAS_MPEG_WIDGET_EVENT, event
    COMPILE_OPT IDL2
    HEAP_GC
    ;bring in the appropriate global variables
    COMMON scan_data
    CI = project.ci


    wWidget =  Event.top

    case Event.id of


       Widget_Info(wWidget, FIND_BY_UNAME='MAS_MPEG_FILE_BUTTON'):begin
         ;the user clicked the button to select a file
         ;have the users select a file
         file_name = dialog_pickfile(PATH=project.current_Path, /WRITE, FILTER='*.mpg')

         ;check to see if the users selected nothing
         if file_name eq '' then return

         ;check to see if they put a .mpg in the file name.if they didnt then add it
         if strpos(file_name, '.mpg') lt 0  then begin
          file_name +='.mpg'
         end

         ;store the path for later.
         project.current_Path = file_name

         ;display the file name in the text widget
         widget_control, Widget_Info(wWidget, FIND_BY_UNAME='MAS_MPEG_FILE_NAME'), $
          SET_VALUE = file_name
       end

       Widget_Info(wWidget, FIND_BY_UNAME='MAS_MPEG_CREATE'):begin
         ;create the mpeg file

         ;first gather all the info from the mpeg widgets
         ;i should do some error checking to make sure that the values are actually number and not strings
         ;for rep_frame and frame_gap
         WIDGET_control, Widget_Info(wWidget, FIND_BY_UNAME='MAS_MPEG_FILE_NAME'), get_value = file_name
         if file_name eq 'Pick a file' then begin
          void = dialog_message( ['Please select a file name'])
                 return
         end

         WIDGET_control, Widget_Info(wWidget, FIND_BY_UNAME='MAS_MPEG_REP_FRAME'), get_value = sTemp
          rep_frame = (fix(sTemp))[0]

         WIDGET_control, Widget_Info(wWidget, FIND_BY_UNAME='MAS_MPEG_FRAME_GAP'), get_value = sTemp
          frame_gap = (fix(sTemp))[0]

         WIDGET_control, Widget_Info(wWidget, FIND_BY_UNAME='MAS_MPEG_QUALITY'), get_value = sTemp
          quality = (fix(sTemp))[0]

         if quality lt 0 or quality gt 100 then begin
          void = dialog_message( ['Please change quality from 0 to 100'])
                 return
         end

         ;check to see if frame_gap is the within the correct range.
         if frame_gap le 0 then begin
          void = dialog_message( ['The Frame Gap can not be less than 1' $
                              ,'Please change Frame Gap'])
              return

         end

         ;check to see if frame_gap is the within the correct range.
         if rep_frame le 0 then begin
          void = dialog_message( ['The Frame Repetions can not be less than 1' $
                              ,'Please change Frame Gap'])
              return

         end


         mas_load_state_2

         ;find out how big the image or array of images is
         mpgSize = size(*project.dataArray[project.ci].state2)

         ;if the image is not 3 dimensions then don't make a movie
         if( not ( mpgSize[0] eq 3) ) then begin
          void = dialog_message( ['Not enough images to make a movie' $
                           ,'Please select multi slice'])
          return
         end

         min_intensity = min(*(project.dataArray[project.ci].state2), max=max_intensity)

           dim = mpgSize[1:2]
           ;print,dim

         ;make a nice little display bar for the users to see how the loading of the movies is going.
         progressbar = Obj_New('progressbar', Color='red', Text='Loading Frames',/NOCANCEL)
           progressbar -> Start


             ;Open an MPEG sequence:
         mpegID = MPEG_OPEN( dim , FILENAME=file_name, quality=quality , IFRAME_GAP=frame_gap)



         ; Add the frames:
         counter=0
         for ii=0,mpgSize[3]-1 do begin

          image = (*project.dataArray[project.ci].state2)[*,*,ii]

          mas_windowing, image, max_intensity, min_intensity

          progressBar -> Update, (float(ii)/float(mpgSize[3]-1))*100.0

          for jj=0, rep_frame-1 do begin

          MPEG_PUT, mpegID, IMAGE=image , FRAME=counter, /ORDER
          counter++

          endfor


         endfor

         progressbar -> Destroy
         text = strjoin('Saving:'+file_name+'  Does not update')
         progressbar = Obj_New('progressbar', Color='red', Text=text,/NOCANCEL)
         progressbar -> Start

           ; Save the MPEG sequence in the file myMovie.mpg:
           MPEG_SAVE, mpegID, FILENAME=mpegFileLocation

           ; Close the MPEG sequence:
           MPEG_CLOSE, mpegID

         progressbar -> Destroy

       end


    else:
       ;If there is an event that i don't capture then do nothing.
       ;the 3 text boxes i don't care about.
    endcase

end


; Subroutine name: MAS_MPEG_WIDGET
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This will create the widget interface to control the movie creation process.

; Editing Information:
    ;Edited by HS 2006/09/30.
    ;Fix spelling mistakes and commenting


PRO MAS_MPEG_WIDGET

    COMPILE_OPT IDL2
    HEAP_GC
    ;bring in the appropriate global variables
    COMMON scan_data
    COMMON common_widgets
    CI = project.ci


    ;If the ADT ROI window is already on the screen then don't make another one.
    IF N_ELEMENTS(MAS_MPEG_WINDOW_BASE) EQ 1 THEN $
        IF     widget_info( MAS_MPEG_WINDOW_BASE,  /VALID_ID ) eq 1 THEN return

    ;if the current scan does not have a good b-matrix then make every thing insensitive
    SENSITIVE =  ptr_valid (project.imndArray[CI].b_matrix)

    title = string ('MAS MPEG Creation')
    MAS_MPEG_WINDOW_BASE = widget_base(TITLE=title, $
            UVALUE = 'MAS_MPEG_WINDOW_BASE', $
            XOFFSET=420 ,YOFFSET=0  ,ROW=5)

    void = widget_button(MAS_MPEG_WINDOW_BASE, value='Choose file name',UNAME = 'MAS_MPEG_FILE_BUTTON',XSIZE = 100)
    void = widget_text(MAS_MPEG_WINDOW_BASE, value='Pick a file' ,XSIZE = 40, /EDITABLE ,UNAME ='MAS_MPEG_FILE_NAME')

    void = widget_label(MAS_MPEG_WINDOW_BASE, value='Repeat Frames (10)',XSIZE = 100, /ALIGN_RIGHT)
    void = widget_text(MAS_MPEG_WINDOW_BASE, value='10' ,XSIZE = 40, /EDITABLE ,UNAME ='MAS_MPEG_REP_FRAME')

    void = widget_label(MAS_MPEG_WINDOW_BASE, value='Frame Gap (2)',XSIZE = 100, /ALIGN_RIGHT)
    void = widget_text(MAS_MPEG_WINDOW_BASE, value='2' ,XSIZE = 40, /EDITABLE ,UNAME ='MAS_MPEG_FRAME_GAP')

    void = widget_label(MAS_MPEG_WINDOW_BASE, value='Quality (80)',XSIZE = 100, /ALIGN_RIGHT)
    void = widget_text(MAS_MPEG_WINDOW_BASE, value='80' ,XSIZE = 40, /EDITABLE ,UNAME ='MAS_MPEG_QUALITY')

    void = widget_button(MAS_MPEG_WINDOW_BASE, value='Create MPEG',UNAME = 'MAS_MPEG_CREATE',XSIZE = 100)

    widget_control, MAS_MPEG_WINDOW_BASE, /realize

    xmanager, 'MAS_MPEG_WIDGET',MAS_MPEG_WINDOW_BASE,/no_block, GROUP_LEADER=WID_BASE_MAIN

END


; Subroutine name: save_project
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
; Calls the dialog box to ask for filename to save. Appends .mas to
; the end of selection.
;
; Purpose of subroutine:
;
; Editing Information:
    ;Edited by HS 2006/09/30.
    ;Fix spelling mistakes and commenting


PRO save_project, Event
    COMPILE_OPT IDL2

    COMMON scan_data, project
    ;if there is no open scan then don't do any thing
    if project.scan_Open_Flag eq 0 then return

    file = Dialog_PickFile(PATH=project.current_path, FILTER='*.mas', title='Save Current Scan')
    save, project, FILENAME=file+'.mas'

END

pro mas_restore_session

    common scan_data

    if project.scan_Open_Flag eq 0 then return

    ci = project.ci

    file = Dialog_PickFile(PATH=project.current_path, FILTER='*.mas', title='Restore Session')

    progressbar = obj_new('progressbar', title="Restoring Session...", text="Please Wait... (progressbar will not update)")
    progressbar->start

    restore, file
    
    progressbar->update, 0.5

    project.imndarray[ci]      = *data[0]
    project.procpramarray[ci]  = *data[1]
    project.dataarray[ci] = *data[2]

    heap_gc
    mas_redraw

    progressbar->destroy

end

; Subroutine name: MAS_SAVE_SESSION
; Created by: BT
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Saves the pertinent data structures for a the currently selected
; scan for later restore.

pro mas_save_session

    common scan_data

    if project.scan_Open_Flag eq 0 then return
    file = Dialog_PickFile(PATH=project.current_path, FILTER='*.mas', Title='Save Current Session')
    
    ci = project.ci

    data = ptrarr(3)

    data[0] = ptr_new(project.imndarray[ci])
    data[1] = ptr_new(project.procpramarray[ci])
    data[2] = ptr_new(project.dataarray[ci])

    save, data, filename=file, /compress

    print, "Save complete."

end

; Subroutine name: MAS_SAVE
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;

; Editing Information:
    ;Edited by HS 2006/09/30.
    ;Fix spelling mistakes and commenting


PRO MAS_SAVE
    COMPILE_OPT IDL2


END
