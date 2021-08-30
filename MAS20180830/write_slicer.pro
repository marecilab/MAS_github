;; $Id$
;;

; Subroutine name: getSlicerPath
; Created by: CD 7/30/07
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: this routine loads a file selection dialog box that takes
; as input the current image directory and moves up one directory to prevent saving
; in the image directory, the user can then select an output directory and provide
; a filename which will be used as an identifier in constructing filenames for the
; files outputed into the slicer directory
;
; Editing Information:

FUNCTION getSlicerPath, current_path

	dir_path = current_path

	pos = STRPOS(STRMID(dir_path,0,STRLEN(dir_path)-1),'\',/REVERSE_SEARCH)
	dir_path = STRMID(dir_path,0,pos)

	path_pattern = DIALOG_PICKFILE( Title = 'Select directory and provide file pattern to save Slicer data.' ,$
    	path=dir_path)

	RETURN, path_pattern
END


; Subroutine name: write_slicer_diffusion_image
; Created by: CD 7/30/07
; Calling Information:

; Bugs or Important Comments to Developers:
; currently only supports non-interleaved(grayscale) images due to slicer file formats
;
; Purpose of subroutine:
; this routine is called by when pressing the output to slicer button on the ADT regression
; GUI, it checks to make sure that the data for all slices is available and that the
; the user has not selected multiADT image output, if the user has selected a grayscale
; ADT image then images for all slices are exported to slicer files, if the user has selected
; an interleaved(color) image then a warning is displayed.  Color image output is not
; currently supported
;
;
; Editing Information:

PRO write_slicer_diffusion_image

    COMMON scan_data
    CI = project.ci
    
    ffov = project.imndArray[ci].f_fov ;freq field of view (mm)
    pfov = project.imndArray[ci].p_fov ;phase field of view (mm)
    sfov = project.imndArray[ci].s_fov ;slice field of view (mm)
    
    ;; check to make sure the user is not expecting the multiple ADT image display
    IF project.procPramArray[CI].adt_display_multi EQ 0 THEN BEGIN
        
        ;; check to make sure that the tensor has been calculated for all slices
        IF project.procPramArray[CI].adt_proccess_flag EQ 1 THEN BEGIN

            ;; from the ADT GUI, the type of image the user has selected for output
            CASE project.procPramArray[CI].adt_display_type OF
                
                0: BEGIN        ;tensor
                    
                    ADT_Start = project.procPramArray[CI].ADT_Start
                    data_out = *project.dataArray[CI].adt ; the previously calculated tensor
                    
                    CASE ADT_Start OF
                        
                        0: BEGIN ;S0
                            void = DIALOG_MESSAGE('ERROR:Color S0 format not supported in export to slicer.',/ERROR)
                            RETURN
                        END
                        
                        7: BEGIN ;Color Trace
                            void = DIALOG_MESSAGE('ERROR:Color Trace format not supported in export to slicer.',/ERROR)
                            RETURN
                        END
                        
                        8: BEGIN ;Orientation
                            void = DIALOG_MESSAGE('ERROR:Color Orientation format not supported in export to slicer.',/ERROR)
                            RETURN
                        END
                        
                        ELSE: BEGIN ;XX,YY,ZZ,XY,XZ,YZ, ADT_Start 1:6
                            
                                ; retrieve min and max scaling parameters from GUI
                            IF ADT_Start GE 1 AND ADT_Start LE 3 THEN BEGIN
                                min_scale = project.procPramArray[CI].dia_min
                                max_scale = project.procPramArray[CI].dia_max
                            ENDIF ELSE IF ADT_Start GE 4 AND ADT_Start LE 6 THEN BEGIN
                                max_scale = project.procPramArray[CI].off_dia
                                min_scale = -max_scale
                            ENDIF ELSE BEGIN
                                void = DIALOG_MESSAGE('ERROR:attempting to write unsupported image format to slicer.',/ERROR)
                                RETURN
                            ENDELSE
                            
                                ; scale the data
                            data_out = REFORM(data_out[*,*,*,ADT_Start])
                            size_out = SIZE(data_out)
                            IF size_out[0] EQ 2 THEN sdim = 1 ELSE sdim = size_out[3]
                            FOR ii=0,sdim-1 DO BEGIN
                                data_out[*,*,ii] = BYTSCL(data_out[*,*,ii],MIN=min_scale,MAX=max_scale)
                            ENDFOR
                        END
                    ENDCASE
                END
                
                                ; retrieve and scale ADC data, scaling parameters from GUI
                2: data_out = BYTSCL(*project.dataArray[project.CI].Avg_Dif, $
                                     MAX=project.procPramArray[CI].avg_d_max, $
                                     MIN=project.procPramArray[CI].avg_d_min) ;ADC
                
                                ; retrieve and scale FA data, scaling parameters from GUI
                3: data_out = BYTSCL(*project.dataArray[project.CI].frac_Ani, $
                                     MAX=project.procPramArray[CI].fa_max) ;FA
                
                4: BEGIN        ;S0 image (state2 data)
                    void = DIALOG_MESSAGE('ERROR:Please use File->Save/Export to export DWI images to slicer.',/ERROR)
                    RETURN
                END
                
                ELSE: BEGIN
                    void = DIALOG_MESSAGE('ERROR:attempting to write unsupported image format to slicer.',/ERROR)
                    RETURN
                END
            ENDCASE
            
        ENDIF ELSE BEGIN
            void = DIALOG_MESSAGE('Please process the entire dataset before writing to slicer.',/ERROR)
            RETURN
        ENDELSE
        
    ENDIF ELSE BEGIN
        void = DIALOG_MESSAGE('Please deselect multi display option before writing to slicer.',/ERROR)
        RETURN
    ENDELSE
    
    size_out = SIZE(data_out)
    
    fdim = size_out[1]          ; matrix size in freq dimension
    pdim = size_out[2]          ; matrix size in phase dimension
    
    IF size_out[0] EQ 2 THEN sdim = 1 ELSE sdim = size_out[3]
    
                                ; get slicer file output path and pattern and output to slicer
    file_path_pattern = getSlicerPath(project.current_path)
    IF (file_path_pattern EQ '') THEN RETURN
    
    output_slicer_files, data_out, fdim, pdim, sdim, -1, ffov, pfov, sfov, file_path_pattern
    
END


; Subroutine name: write_slicer_image
; Created by: CD and BL 7/30/07
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; this routine is called when File->Save/Export->Output to Slicer is called
; it output image (State1) data to slicer, it the user has selected multiple from
; the single/multiple toggle on the GUI then will output individual datasets for
; each adim, if single is selected then outputs the currently selected adim
;
; Editing Information:

PRO write_slicer_image
    COMMON scan_data
    ci = project.ci

    file_path_pattern = getSlicerPath(project.current_path)
    IF (file_path_pattern EQ '') THEN RETURN

	; make sure multi-slice state2 data has been calculated for currently
	; selected adim
    curr_state_2 = project.procPramArray[project.ci].state_2
    curr_multi = project.procPramArray[CI].single_Multi_flag
    
    project.procPramArray[project.ci].state_2   = 0
    project.procPramArray[CI].single_Multi_flag = 1
    
    mas_load_state_2 ;check to make sure selected data set is loaded, if not this loads it
    
    project.procPramArray[project.ci].state_2   = curr_state_2
    project.procPramArray[CI].single_Multi_flag = curr_multi
    
    
    adim = project.imndArray[ci].adim ;dim4, probably 1 for non DTI
    
    ffov = project.imndArray[ci].f_fov ;freq field of view (mm)
    pfov = project.imndArray[ci].p_fov ;phase field of view (mm)
    sfov = project.imndArray[ci].s_fov ;slice field of view (mm)

    imgdata = *project.dataArray[CI].state2 ; image matrix fdim x pdim x sdim x adim (4D)
    
    size_img = SIZE(imgdata)
    
    fdim = size_img[1]          ; matrix size in freq dimension
    pdim = size_img[2]          ; matrix size in phase dimension
    
    IF size_img[0] EQ 2 THEN sdim = 1 ELSE sdim = size_img[3]
    
    adim_start = project.procpramArray[ci].adim_start
    output_slicer_files, imgdata, fdim, pdim, sdim, adim_start+1, ffov, pfov, sfov, file_path_pattern
    
	; from display tab of GUI, if set to single output only current adim, otherwise
	; output all adims
	;IF project.procPramArray[CI].single_Multi_flag EQ 0 THEN BEGIN
	;	; selects the current adim to output to Slicer
	;	adim_start = project.procpramArray[ci].adim_start
	;	adim_end   = adim_start
	;ENDIF ELSE BEGIN
		;adim_start = 0
		;adim_end   = adim-1
	;ENDELSE

   	;progressbar = OBJ_NEW('progressbar', Color='red', Text='Saving all adims to SLICER format...',/NOCANCEL)
   	;progressbar -> Start

	;FOR ii = adim_start, adim_end DO BEGIN
		;output_slicer_files, imgdata[*,*,*,ii], fdim, pdim, sdim, ii+1, ffov, pfov, sfov, file_path_pattern
		;progressBar -> Update, (float(ii+1)/float(adim_end+1))*100.0
	;ENDFOR

	;progressbar -> Destroy
END


; Subroutine name: output_slicer_files
; Created by: CD and BL on 7/30/07
; Calling Information:

; Bugs or Important Comments to Developers:
; this routine should be called individually for each adim to maintain
; generallity and compatibility with ADT slicer output
;
; Space directions are incorrect and cause Slicer to tag the directions in a wrong way
;
; Purpose of subroutine:
; outputs the 3D dataset in imgdata (freq/phase/slice) to slicer files in
; the directory with name scheme specified in file_path_pattern, adim is
; a numerical label for which adim is contained in imgdata, also outputs
; the nrrd header for slicer files
;
; Editing Information:

PRO output_slicer_files, imgdata, fdim, pdim, sdim, adim, ffov, pfov, sfov, file_path_pattern

	SCALE_FACTOR = 255 ;the maximum pixel intensity

	imgdata = REFORM(FLOAT(imgdata))
	imgdata = imgdata/MAX(imgdata) * SCALE_FACTOR; ; Rescale image matrix from 0-SCALE_FACTOR

	sz_img = SIZE(imgdata)
	IF sz_img[0] EQ 2 THEN sdim = 1 ELSE sdim = sz_img[3]

	FOR jj = 0, sdim - 1 DO BEGIN

		slice_data = REFORM(imgdata[*,*,jj]) ; 2D matrix (xy) for a single image slice

		; allows index to be written to 4 total characters
		IF ((jj+1)/10000.0 LT 1.0) THEN BEGIN      		  ; example: 1024
			leading_zeros = ''
		ENDIF

		IF ((jj+1)/1000.0 LT 1.0) THEN BEGIN   ; example: 0235
			leading_zeros = '0'
		ENDIF

		IF ((jj+1)/100.0 LT 1.0) THEN BEGIN    ; example: 0045
			leading_zeros = '00'
		ENDIF

		IF ((jj+1)/10.0 LT 1.0) THEN BEGIN     ; example: 0002
			leading_zeros = '000'
		ENDIF

		IF adim EQ -1 THEN BEGIN
			;filename = file_pattern + zero-padded slice #
			fname = file_path_pattern + '.' + leading_zeros + STRTRIM(STRING(jj+1),2)
		ENDIF ELSE BEGIN
			;filename = file_pattern + adim + adim# + zero-padded slice #
			fname = file_path_pattern + '_adim_' + STRTRIM(STRING(adim),2) + $
				'.' + leading_zeros + STRTRIM(STRING(jj+1),2)
		ENDELSE

		OPENW, 1, fname  			; Open a new file for writing as IDL file unit number 1
		WRITEU, 1, slice_data  		; Write the data in scaled_slice_data to the file
		CLOSE, 1  					; Close file unit 1
	ENDFOR

	; write nrrd file header for slicer data

	IF adim EQ -1 THEN BEGIN
		;filename = file_pattern + zero-padded slice #
		fname = file_path_pattern + '.' + leading_zeros + STRTRIM(STRING(jj+1),2)
	ENDIF ELSE BEGIN
		;filename = file_pattern + adim + adim# + zero-padded slice #
		fname = file_path_pattern + '_adim_' + STRTRIM(STRING(adim),2) + $
			'.' + leading_zeros + STRTRIM(STRING(jj+1),2)
	ENDELSE

	IF adim EQ -1 THEN BEGIN
		fpattern_nrrd = file_path_pattern
	ENDIF ELSE BEGIN
		fpattern_nrrd = file_path_pattern + '_adim_' + STRTRIM(STRING(adim),2)
	ENDELSE

	slash = get_dir_slash()
    pos = STRPOS(fpattern_nrrd,slash[1],/REVERSE_SEARCH)
    fpattern_nrrd_nodir = STRMID(fpattern_nrrd,pos+1,STRLEN(fpattern_nrrd))

	nrrd_version    = 'NRRD0005'
	nrrd_content    = 'content: 3D Slicer'
	nrrd_data_type  = 'type: float'
	nrrd_dimension  = 'dimension: ' + '3'
	nrrd_sizes      = 'sizes: ' + STRTRIM(STRING(fdim),2) + $
                          ' '       + STRTRIM(STRING(pdim),2) + $
                          ' '       + STRTRIM(STRING(sdim),2)
	nrrd_orientation= 'space: right-anterior-superior'
        ;; This space_dir is incorrectly defined as it only works
        ;; for a specific case. this has to be fixed to a general solution
        nrrd_space_dir  = 'space directions: (' + STRTRIM(STRING(ffov*10/fdim),2) + $
                          ',0,0) (0,'           + STRTRIM(STRING(pfov*10/pdim),2) + $
                          ',0) (0,0,'           + STRTRIM(STRING(sfov*10/sdim),2) + $
                          ')'
	nrrd_centerings = 'centerings: cell cell cell'
	nrrd_kinds      = 'kinds: space space space'
	nrrd_endian     = 'endian: little'
	nrrd_encoding   = 'encoding: raw'
	nrrd_units      = 'space units: "mm" "mm" "mm"'
	nrrd_data_file  = 'data file: ' + fpattern_nrrd_nodir + $
                          '.%04d '      + '1 ' + STRTRIM(STRING(sdim),2) + ' 1'

	nrrd_fname = fpattern_nrrd + '.nrrd' ;filename = file_pattern + adim + adim#

	OPENW, 1, nrrd_fname  			; Open a new file for writing as IDL file unit number 1
	PRINTF, 1, nrrd_version 		; Write the "version" field to the file
	PRINTF, 1, nrrd_content  		; Write the "content" field to the file
	PRINTF, 1, nrrd_data_type  		; Write the "type" field to the file
	PRINTF, 1, nrrd_dimension  		; Write the "dimension" field to the file
	PRINTF, 1, nrrd_sizes  			; Write the "sizes" field to the file
	PRINTF, 1, nrrd_orientation 	        ; Write the "space" field to the file
	PRINTF, 1, nrrd_space_dir  		; Write the "space directions" field to the file
	PRINTF, 1, nrrd_centerings  	        ; Write the "centerings" field to the file
	PRINTF, 1, nrrd_kinds  			; Write the "kinds" field to the file
	PRINTF, 1, nrrd_endian  		; Write the "endian" field to the file
	PRINTF, 1, nrrd_encoding  		; Write the "encoding" field to the file
	PRINTF, 1, nrrd_units  			; Write the "space units" field to the file
	PRINTF, 1, nrrd_data_file  		; Write the "data file" field to the file

	CLOSE, 1  						; Close file unit 1
END

