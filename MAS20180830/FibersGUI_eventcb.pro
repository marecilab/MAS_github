;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved
;
;
; IDL Event Callback Procedures
; FibersGUI_eventcb
;
; Generated on:	01/11/2007 09:30.38
;
;-----------------------------------------------------------------
; Notify Realize Callback Procedure.
; Argument:
;   wWidget - ID number of specific widget.
;
;
;   Retrieve the IDs of other widgets in the widget hierarchy using
;       id=widget_info(Event.top, FIND_BY_UNAME=name)

;-----------------------------------------------------------------
pro WID_BASE_Main_Realized, wWidget
	compile_opt idl2, hidden

	COMMON scan_data
	;COMMON common_widgets
   	COMMON share_fibers_main, pTHR, pFA, pS_0, pAvD, pDIR, loadedthr, loadeddir, resx,$
   	resy, resz, interpol, rfx, ss, thr, ll, sd, fname, tfile, dfile, gfile, path, angle_limit $
   	, selectedslices, ROI_seeds_global, transform_matrix, project_dataset

        CI = project.ci

	cd, current=path
	cd, path

	slash = get_dir_slash()

	id=widget_info(wWidget, FIND_BY_UNAME='WID_TEXT_Output_Filename')
	widget_control, id, SET_VALUE = path+slash[1]+'tracts.trs'

    id=widget_info(wWidget, FIND_BY_UNAME='WID_TEXT_Transform_File')
    widget_control, id, SET_VALUE = ''

    id=widget_info(wWidget, FIND_BY_UNAME='WID_TEXT_Project_Dataset')
    widget_control, id, SET_VALUE = string(CI, format='(I2)')
    
	interpol = '1' ;This is the default interpolation
	tfile = 'Fract_Anisotropy'
	dfile = 'Eigenvectors'
	gfile = 'Fract_Anisotropy'

;	id = widget_info(wWidget, FIND_BY_UNAME='WID_LABEL_ResolutionTitle')
;	widget_control, id, /DYNAMIC_RESIZE
        pixszs = [ project.imndarray[ci].f_voxsz / project.procpramarray[ci].freq_interp, $
                   project.imndarray[ci].p_voxsz / project.procpramarray[ci].phase_interp, $
                   project.imndarray[ci].s_voxsz / project.procpramarray[ci].slice_interp ]
        ;;pixszs /= max(pixszs)

        id=widget_info(wWidget, FIND_BY_UNAME='WID_Text_Res_X_Change')
        
        matrixdimension = project.imndArray[ci].fdim
        directiondimension = project.imndArray[ci].f_fov
        
                                ; the if statement is to avoid any divide by zeros
        if (matrixdimension ne 0) then begin
           temp = directiondimension/matrixdimension
           temp = STRCOMPRESS(string(temp))

           temp = strcompress(string(pixszs[0]), /remove_all)
           widget_control, id, SET_VALUE=(temp)
           
                                ; For the Y Textbox
           matrixdimension = project.imndArray[ci].pdim
           directiondimension = project.imndArray[ci].p_fov
           temp = directiondimension/matrixdimension
           temp = strcompress(string(pixszs[1]), /remove_all)
           id=widget_info(wWidget, FIND_BY_UNAME='WID_Text_Res_Y_Change')
           
           widget_control, id, SET_VALUE=STRCOMPRESS(string(temp))
           
                                ; For the Z Textbox
           matrixdimension = project.imndArray[ci].sdim
           directiondimension = project.imndArray[ci].s_fov
           temp = directiondimension/matrixdimension
           temp = strcompress(string(pixszs[2]), /remove_all)
           id=widget_info(wWidget, FIND_BY_UNAME='WID_Text_Res_Z_Change')
           
           widget_control, id, SET_VALUE=STRCOMPRESS(string(temp))
           
        endif

end
;-----------------------------------------------------------------
; Activate Button Callback Procedure.
; Argument:
;   Event structure:
;
;   {WIDGET_BUTTON, ID:0L, TOP:0L, HANDLER:0L, SELECT:0}
;
;   ID is the widget ID of the component generating the event. TOP is
;       the widget ID of the top level widget containing ID. HANDLER
;       contains the widget ID of the widget associated with the
;       handler routine.

;   SELECT is set to 1 if the button was set, and 0 if released.
;       Normal buttons do not generate events when released, so
;       SELECT will always be 1. However, toggle buttons (created by
;       parenting a button to an exclusive or non-exclusive base)
;       return separate events for the set and release actions.

;   Retrieve the IDs of other widgets in the widget hierarchy using
;       id=widget_info(Event.top, FIND_BY_UNAME=name)

;-----------------------------------------------------------------

pro WID_BUTTON_SpecifyROI_Pressed, Event

    common scan_data
    common share_fibers_main
    common share_fibers_roiobj, roi_objects
    
    ci = project.ci
;;>>>>>>
;    project.dataArray[project.CI].fibers_frac_Ani = ptr_new(*project.dataarray[ci].frac_ani)
;    project.dataArray[project.CI].fibers_Avg_Dif = ptr_new(project.dataArray[project.CI].Avg_Dif)
;;<<<<<<

    pTHR      = project.dataArray[project.CI].fibers_frac_Ani
    pFA       = project.dataArray[project.CI].fibers_frac_Ani
    ;; pFA and pTHR are the SAME pointer!

;    pS_0_temp = (*project.dataArray[project.CI].fibers_adt)[*,*,*,0]
;    pS_0      = ptr_new(pS_0_temp)

    pS_0      = ptr_new((*project.dataArray[project.CI].adt)[*,*,*,0])
    pAvD      = project.dataArray[project.CI].fibers_Avg_Dif
    
    ;; Only the first 5th dimension is read in the eigenvector matrix. This
    ;; is how Evren did it for the files stored in the directory.
    pDIR_TEMP = reform((*project.dataArray[project.CI].fibers_eign_Vec)[*,*,*,*,0])
;;    pDIR_TEMP = reform((*project.dataArray[project.CI].eign_Vec)[*,*,*,*,0])
    pDIR      = ptr_new(pDIR_TEMP)

    
    ;; We only need to rotate pTHR OR pFA because they are the same pointer
    mas_rotate_flip, pTHR

    ;; Rotate the other images
    mas_rotate_flip, pS_0
    mas_rotate_flip, pAvD
    mas_rotate_flip, pDIR
    
    ;; Read the selected values from the GUI.
    ;; resx
    id=widget_info(Event.top, FIND_BY_UNAME='WID_Text_Res_X_Change')
    widget_control, id, GET_VALUE=(temp)
    temp = float(temp)
    resx = temp

    ;; resy
    id=widget_info(Event.top, FIND_BY_UNAME='WID_Text_Res_Y_Change')
    widget_control, id, GET_VALUE=(temp)
    temp = float(temp)
    resy = temp

    ;; resz
    id=widget_info(Event.top, FIND_BY_UNAME='WID_Text_Res_Z_Change')
    widget_control, id, GET_VALUE=(temp)
    temp = float(temp)
    resz = temp

    ;; Visualization Factor (rfx for some reason)
    id=widget_info(Event.top, FIND_BY_UNAME='WID_SLIDER_VisualizationFactor')
    widget_control, id, GET_VALUE=(temp)
    temp = float(temp)
    rfx = temp

    ;; step size
    id=widget_info(Event.top, FIND_BY_UNAME='WID_SLIDER_StepSize')
    widget_control, id, GET_VALUE=(temp)
    temp = float(temp)
    ss = 0.1*temp
    
    ;; threshold
    id=widget_info(Event.top, FIND_BY_UNAME='WID_SLIDER_Threshold')
    widget_control, id, GET_VALUE=(temp)
    temp = float(temp)
    thr = 0.005*temp
    
    
    ;; length limit
    id=widget_info(Event.top, FIND_BY_UNAME='WID_SLIDER_LengthLimit')
    widget_control, id, GET_VALUE=(temp)
    temp = float(temp)
    ll = 0.1*temp
    
    ;; seed density
    id=widget_info(Event.top, FIND_BY_UNAME='WID_SLIDER_SeedDensity')
    widget_control, id, GET_VALUE=(temp)
    temp = float(temp)
    sd = temp
    
    ;; Filename to save ROI
    id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_Output_Filename')
    widget_control, id, GET_VALUE=(temp)
    fname = temp
    lastslash = STRPOS(fname,'.',REVERSE_SEARCH=1)
    fname     = STRMID(fname,0,lastslash)
    roi_fname = fname+'.ROI'
    trs_fname = fname+'.trs'

   ;; Filename to find transformation matrix
    id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_Transform_File')
    widget_control, id, GET_VALUE=(temp)
    if (file_test(temp, /read)) then begin
        transform_matrix = dblarr(4,4)
        openr, lun, temp, /get_lun
        readf, lun, transform_matrix
        close, lun
        free_lun, lun
        print, transform_matrix
    endif else begin
        transform_matrix = diag_matrix(replicate(1D, 4))
    endelse
    
    id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_Project_Dataset')
    widget_control, id, get_value=(temp)
    project_dataset = fix(temp)
    
    id=widget_info(Event.top, FIND_BY_UNAME='WID_SLIDER_AngleLimit')
    widget_control, id, get_value=temp
    angle_limit = cos(float(temp)*!PI/180);;float(temp)/100.0
    
    ;; Need to make sure that we can write tracts.{ROI,trs} files and
    ;; that there isn't already an existing tracts job in here

    ok_to_save = 1
    if (file_test(roi_fname, /regular)) then begin

        if (not file_test(roi_fname, /write, /regular)) then begin
            junk = dialog_message(roi_fname+' exists but is not writable. Please rename '+$
                                  'the existing tracts files.', /center)
            return
        endif else begin
            junk = dialog_message([roi_fname+' exists.', $
                                  'Are you sure you want to overwrite it?'],/question, /center)
            if (junk eq 'No') then return ;;else  ok to overwrite.
        endelse

    endif else if (not file_test(file_dirname(fname), /directory, /write)) then begin
        print, fname
        junk = dialog_message('Unable to create tracts.ROI file. Parent directory is not writable.', /error, /center)
        return
    endif
  
    if (file_test(trs_fname, /regular)) then begin

        if (not file_test(trs_fname, /write, /regular)) then begin
            junk = dialog_message([trs_fname+' exists but is not writable. Please rename ',$
                                  'the existing tracts files.'] , /center, /error)
            return
        endif else begin
            junk = dialog_message([trs_fname+' exists.', $
                                  'Are you sure you want to overwrite it?'] , /question, /center)
            if (junk eq 'No') then return ;;else  ok to overwrite.
        endelse

    endif else if (not file_test(file_dirname(fname), /directory, /write)) then begin
        junk = dialog_message('Unable to create tracts.trs file. Parent directory is not writable.', /center, /error)
        return
    endif
  
    if (n_elements(roi_objects) ne 0) then begin
        ptr_free, roi_objects
    endif 
    
    roi_objects = ptrarr(1)

    ;; Going to call the routine to do the ROI selection
    ptr_seeds=ft_roi(resx, resy, resz, interpol)
    
    num_rois = total(ptr_valid(roi_objects))

    if (fix(num_rois) le 0) then begin

        junk = dialog_message('No ROIs specified.', /information, /center)
        ptr_free, pS_0
        ptr_free, PDIR
        return
        
    endif

    tmp_rois = ptrarr(num_rois)

    ct = 0
    for i=0, n_elements(roi_objects)-1 do begin
        if ptr_valid(roi_objects[i]) then begin
            tmp_rois[ct++] = roi_objects[i]
        endif
    endfor
    roi_objects = tmp_rois

    ;; If the user returns an ROI value then the tracking can proceed
    if total(ptr_valid(ptr_seeds)) then begin
        
        temp = *ptr_seeds[0]
        temp_size = size(temp)
        validpointer = 0
        
        ;; going to check each pointer with ptr_valid. This returns a 0 for a bad pointer, 1 for a valid.
        ;; The total sum of the whole for loop should be the total number of pointers, else one is invalid
        for i=0, temp_size[1]-1 do begin
            badpointer_flag = ptr_valid(temp[i])
            validpointer = validpointer+badpointer_flag
        endfor
        
        if (validpointer eq temp_size[1]) then begin 
            ;;this one checks to see if the user closed the window before selecting an ROI
            ;; Every single pointer in the temp variable has to return valid else it does not execute.

            ;; Now Im saving the ROIs in a global variable that contains a pointer for each
            ;; ROI
            ;lastslash = STRPOS(fname,'.',REVERSE_SEARCH=1)
            ;fname     = STRMID(fname,0,lastslash)
            ROI_matrix = (size(*project.dataarray[project.ci].state1))[1:3]
            if (ok_to_save eq 1) then $
              SAVE, ptr_seeds, roi_objects, ROI_matrix, FILENAME = roi_fname ;+'.ROI'

            temp = *ptr_seeds[0]
            
            num_rois = size(temp)
            temp_roi_array = ptrarr(num_rois[1])
            
            for i=0,num_rois[1]-1 do begin
                roitemp = *temp[i]
                temp_roi_array[i] = ptr_new(roitemp)
            end
            
            ROI_seeds_global = temp_roi_array
            
            print, 'Calculating Tracts'
            ft_track_newversion, ptr_seeds
            
        endif
    endif
    
;    ptr_free, pTHR
;    ptr_free, pS_0
;    ptr_free, pAvD
;    ptr_free, pDIR
    
    ;; Now I have to restore the original values of the ADT data in case the user closes Fibers and
    ;; restarts it later (because we need to reset the rotation of images).
    project.dataArray[project.CI].fibers_adt = ptr_new(*project.dataArray[project.CI].adt)
    project.dataArray[project.CI].fibers_eign_val = ptr_new(*project.dataArray[project.CI].eign_val)
    project.dataArray[project.CI].fibers_eign_Vec = ptr_new(*project.dataArray[project.CI].eign_Vec)
    project.dataArray[project.CI].fibers_frac_Ani = ptr_new(*project.dataArray[project.CI].frac_Ani)
    project.dataArray[project.CI].fibers_Avg_Dif = ptr_new(*project.dataArray[project.CI].Avg_Dif)

end

;; Gets the transformation matrix file from the disk and stores the path in the 
;; appropriate text box

pro WID_BUTTON_Transform_File_Pressed, Event

    common share_fibers_main
    common scan_data
    
    id = widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_Transform_File')
    
    tx_filepath = dialog_pickfile(path=project.current_path, /file, /read)
    
    if (tx_filepath eq '') then return
    
    widget_control, id, set_value=tx_filepath

end

;-----------------------------------------------------------------
; Activate Button Callback Procedure.
; Argument:
;   Event structure:
;
;   {WIDGET_BUTTON, ID:0L, TOP:0L, HANDLER:0L, SELECT:0}
;
;   ID is the widget ID of the component generating the event. TOP is
;       the widget ID of the top level widget containing ID. HANDLER
;       contains the widget ID of the widget associated with the
;       handler routine.

;   SELECT is set to 1 if the button was set, and 0 if released.
;       Normal buttons do not generate events when released, so
;       SELECT will always be 1. However, toggle buttons (created by
;       parenting a button to an exclusive or non-exclusive base)
;       return separate events for the set and release actions.

;   Retrieve the IDs of other widgets in the widget hierarchy using
;       id=widget_info(Event.top, FIND_BY_UNAME=name)

;-----------------------------------------------------------------
pro WID_BUTTON_RepeatTrack_Pressed, Event
	common scan_data
	common share_fibers_main

        ROI_matrix = [0,0,0]
        
	ROI_file = DIALOG_PICKFILE( FILTER = '*.ROI', Title = 'Select ROI file to open')

	; If empty we don't do anything.
        if (ROI_file eq '') then return

	RESTORE, ROI_file

        if (total(ROI_matrix) ne 0) then begin
            size_test = total(ROI_matrix - (size(*project.dataarray[project.ci].state1))[1:3])
            if (size_test ne 0) then begin
                junk = dialog_message('ROI dimensions do not match current project dimensions.')
                return
            endif
        endif
        
	pTHR = project.dataArray[project.CI].fibers_frac_Ani
        pFA  = project.dataArray[project.CI].fibers_frac_Ani
	; pFA and pTHR are the SAME pointer!

	pS_0_temp = (*project.dataArray[project.CI].fibers_adt)[*,*,*,0]
        pS_0      = ptr_new(pS_0_temp)

        pAvD      = project.dataArray[project.CI].fibers_Avg_Dif

        ;; Only the first 5th dimension is read in the eigenvector matrix. This
        ;; is how Evren did it for the files stored in the directory.
        pDIR_TEMP = reform((*project.dataArray[project.CI].fibers_eign_Vec)[*,*,*,*,0])
        pDIR      = ptr_new(pDIR_TEMP)
        
        ;; We only need to rotate pTHR OR pFA because they are the same pointer
	mas_rotate_flip, pTHR
	;mas_rotate_flip, pFA

        ;; Rotate the other images
	mas_rotate_flip, pS_0
	mas_rotate_flip, pAvD
	mas_rotate_flip, pDIR

        ;; Read the selected values from the GUI.
	; resx
	id=widget_info(Event.top, FIND_BY_UNAME='WID_Text_Res_X_Change')
        widget_control, id, GET_VALUE=(temp)
        temp = float(temp)
        resx = temp

	; resy
        id=widget_info(Event.top, FIND_BY_UNAME='WID_Text_Res_Y_Change')
        widget_control, id, GET_VALUE=(temp)
        temp = float(temp)
        resy = temp

	; resz
	id=widget_info(Event.top, FIND_BY_UNAME='WID_Text_Res_Z_Change')
        widget_control, id, GET_VALUE=(temp)
        temp = float(temp)
        resz = temp

	; Visualization Factor (rfx for some reason)
        id=widget_info(Event.top, FIND_BY_UNAME='WID_SLIDER_VisualizationFactor')
        widget_control, id, GET_VALUE=(temp)
        temp = float(temp)
        rfx = temp


	; step size
	id=widget_info(Event.top, FIND_BY_UNAME='WID_SLIDER_StepSize')
		widget_control, id, GET_VALUE=(temp)
		temp = float(temp)
	ss = 0.1*temp


	; threshold
	id=widget_info(Event.top, FIND_BY_UNAME='WID_SLIDER_Threshold')
		widget_control, id, GET_VALUE=(temp)
		temp = float(temp)
	thr = 0.005*temp


	; length limit
	id=widget_info(Event.top, FIND_BY_UNAME='WID_SLIDER_LengthLimit')
		widget_control, id, GET_VALUE=(temp)
		temp = float(temp)
	ll = 0.1*temp

	; seed density
	id=widget_info(Event.top, FIND_BY_UNAME='WID_SLIDER_SeedDensity')
		widget_control, id, GET_VALUE=(temp)
		temp = float(temp)
	sd = temp

	; Filename
	id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_Output_Filename')
		widget_control, id, GET_VALUE=(temp)
	fname = temp

   ;; Filename to find transformation matrix
    id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_Transform_File')
    widget_control, id, GET_VALUE=(temp)
    if (file_test(temp, /read)) then begin
        transform_matrix = dblarr(4,4)
        openr, lun, temp, /get_lun
        readf, lun, transform_matrix
        close, lun
        free_lun, lun
        print, transform_matrix
    endif else begin
        transform_matrix = diag_matrix(replicate(1D, 4))
    endelse
    
    id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_Project_Dataset')
    widget_control, id, get_value=(temp)
    project_dataset = fix(temp)
 
        id=widget_info(Event.top, FIND_BY_UNAME='WID_SLIDER_AngleLimit')
        widget_control, id, get_value=temp
        angle_limit = cos(float(temp)*!PI/180.0); float(temp)/100.0

	lastslash = STRPOS(fname,'.',REVERSE_SEARCH=1)
	fname = STRMID(fname,0,lastslash)


	ft_track_newversion, ptr_seeds

; Now I have to restore the original values of the ADT data in case the user closes Fibers and
; restarts it later (because we need to reset the rotation of images).
project.dataArray[project.CI].fibers_adt = ptr_new(*project.dataArray[project.CI].adt)
project.dataArray[project.CI].fibers_eign_val = ptr_new(*project.dataArray[project.CI].eign_val)
project.dataArray[project.CI].fibers_eign_Vec = ptr_new(*project.dataArray[project.CI].eign_Vec)
project.dataArray[project.CI].fibers_frac_Ani = ptr_new(*project.dataArray[project.CI].frac_Ani)
project.dataArray[project.CI].fibers_Avg_Dif = ptr_new(*project.dataArray[project.CI].Avg_Dif)


end


; Routine that selects the current path for the fiber track mapping

PRO WID_BUTTON_SpecifyDirectory_Pressed, event, recurse=recurse

	common scan_data
	common share_fibers_main

        if (not keyword_set(recurse)) then begin
            recurse = 0
        endif else if recurse gt 4 then begin
            junk = dialog_message('Recurse limit reached. Please try again.', /error, /center)
            return
        endif
        
        sFolder = DIALOG_PICKFILE(PATH=path, /DIRECTORY, TITLE="Choose directory containing save ROI file.")
        if (sFolder eq '') then return
        
        path = sFolder

        duh = file_test(path, /read, /executable, /directory)
        if (duh eq 0) then begin
            junk = dialog_message(['You do not have read permissions in this directory.',$
                                  'Please select a readable directory.'], /center)
            WID_BUTTON_SpecifyDirectory_Pressed, Event, recurse=(recurse+1)
            return
        endif

        cd, path
        
        slash = get_dir_slash()
        
        id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_Output_Filename')
	    widget_control, id, SET_VALUE = path+'tracts.trs'
        id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_Transform_File')
        widget_control, id, SET_VALUE = path+'transform.matrix'

END

pro WID_BUTTON_Revisualize_Pressed, event

    common scan_data
    common share_fibers_main

    id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_Output_Filename')
    widget_control, id, GET_VALUE=tract_file
    
    if not file_test(tract_file, /regular, /read) then begin
        junk = dialog_message('Tract file: '+tract_file+' not found or cannot be read.', /error, /center)
        return
    endif
    
   ;; Filename to find transformation matrix
    id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_Transform_File')
    widget_control, id, GET_VALUE=(temp)
    if (file_test(temp, /read)) then begin
        transform_matrix = dblarr(4,4)
        openr, lun, temp, /get_lun
        readf, lun, transform_matrix
        close, lun
        free_lun, lun
        print, transform_matrix
    endif else begin
        transform_matrix = diag_matrix(replicate(1D, 4))
    endelse
    
    id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_Project_Dataset')
    widget_control, id, get_value=(temp)
    project_dataset = fix(temp)

    ;;fiber_data = make_fibers(tract_file)
    make_fibers, tract_file, fiber_data=fiber_data, transform=invert(transform_matrix)

    ;datasize = size(*project.dataarray[project.ci].state1)
    
    ;result = total([ datasize[1], datasize[2], datasize[3] ] - (*fiber_data).matrix)
    
    ;if (result ne 0) then begin
    ;    junk = dialog_message('Project dimensions do not match '+$
    ;                          'tract file. Aborting', /error, /center)
    ;    ptr_free, (*fiber_data).fibers
    ;    ptr_free, fiber_data
    ;    return
    ;endif

    view_fibers, fiber_data, thick=1, alpha=.5, /show_axis, /in_mas, project_dataset=project_dataset

end

;-----------------------------------------------------------------
; Slider Value Changed Callback Procedure.
; Argument:
;   Event structure:
;
;   {WIDGET_SLIDER, ID:0L, TOP:0L, HANDLER:0L, VALUE:0L, DRAG:0}
;
;   ID is the widget ID of the component generating the event. TOP is
;       the widget ID of the top level widget containing ID. HANDLER
;       contains the widget ID of the widget associated with the
;       handler routine.

;   VALUE returns the new value of the slider. DRAG returns integer 1
;       if the slider event was generated as part of a drag
;       operation, or zero if the event was generated when the user
;       had finished positioning the slider.

;   Retrieve the IDs of other widgets in the widget hierarchy using
;       id=widget_info(Event.top, FIND_BY_UNAME=name)

;-----------------------------------------------------------------
pro WID_SLIDER_StepSize_ChangeValue, Event

	id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_StepSizeValue')
	widget_control, id, SET_VALUE=STRCOMPRESS(string(Event.value*0.1))


end
;-----------------------------------------------------------------
; Slider Value Changed Callback Procedure.
; Argument:
;   Event structure:
;
;   {WIDGET_SLIDER, ID:0L, TOP:0L, HANDLER:0L, VALUE:0L, DRAG:0}
;
;   ID is the widget ID of the component generating the event. TOP is
;       the widget ID of the top level widget containing ID. HANDLER
;       contains the widget ID of the widget associated with the
;       handler routine.

;   VALUE returns the new value of the slider. DRAG returns integer 1
;       if the slider event was generated as part of a drag
;       operation, or zero if the event was generated when the user
;       had finished positioning the slider.

;   Retrieve the IDs of other widgets in the widget hierarchy using
;       id=widget_info(Event.top, FIND_BY_UNAME=name)

;-----------------------------------------------------------------
pro WID_SLIDER_Threshold_ChangeValue, Event

	id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_ThresholdValue')
	widget_control, id, SET_VALUE=STRCOMPRESS(string(Event.value*0.005))

end
;-----------------------------------------------------------------
; Slider Value Changed Callback Procedure.
; Argument:
;   Event structure:
;
;   {WIDGET_SLIDER, ID:0L, TOP:0L, HANDLER:0L, VALUE:0L, DRAG:0}
;
;   ID is the widget ID of the component generating the event. TOP is
;       the widget ID of the top level widget containing ID. HANDLER
;       contains the widget ID of the widget associated with the
;       handler routine.

;   VALUE returns the new value of the slider. DRAG returns integer 1
;       if the slider event was generated as part of a drag
;       operation, or zero if the event was generated when the user
;       had finished positioning the slider.

;   Retrieve the IDs of other widgets in the widget hierarchy using
;       id=widget_info(Event.top, FIND_BY_UNAME=name)

;-----------------------------------------------------------------
pro WID_SLIDER_SeedDensity_ChangeValue, Event

	id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_SeedDensityValue')
	widget_control, id, SET_VALUE=STRCOMPRESS(string(Event.value))

end

pro WID_SLIDER_AngleLimit_ChangeValue, Event

    id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_AngleLimitValue')
    widget_control, id, SET_VALUE=STRCOMPRESS(string(Event.value))

end

;-----------------------------------------------------------------
; Slider Value Changed Callback Procedure.
; Argument:
;   Event structure:
;
;   {WIDGET_SLIDER, ID:0L, TOP:0L, HANDLER:0L, VALUE:0L, DRAG:0}
;
;   ID is the widget ID of the component generating the event. TOP is
;       the widget ID of the top level widget containing ID. HANDLER
;       contains the widget ID of the widget associated with the
;       handler routine.

;   VALUE returns the new value of the slider. DRAG returns integer 1
;       if the slider event was generated as part of a drag
;       operation, or zero if the event was generated when the user
;       had finished positioning the slider.

;   Retrieve the IDs of other widgets in the widget hierarchy using
;       id=widget_info(Event.top, FIND_BY_UNAME=name)

;-----------------------------------------------------------------
pro WID_SLIDER_VisualizationFactor_ChangeValue, Event

	id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_VisualizationFactorValue')
	widget_control, id, SET_VALUE=STRCOMPRESS(string(Event.value))

end
;-----------------------------------------------------------------
; Slider Value Changed Callback Procedure.
; Argument:
;   Event structure:
;
;   {WIDGET_SLIDER, ID:0L, TOP:0L, HANDLER:0L, VALUE:0L, DRAG:0}
;
;   ID is the widget ID of the component generating the event. TOP is
;       the widget ID of the top level widget containing ID. HANDLER
;       contains the widget ID of the widget associated with the
;       handler routine.

;   VALUE returns the new value of the slider. DRAG returns integer 1
;       if the slider event was generated as part of a drag
;       operation, or zero if the event was generated when the user
;       had finished positioning the slider.

;   Retrieve the IDs of other widgets in the widget hierarchy using
;       id=widget_info(Event.top, FIND_BY_UNAME=name)

;-----------------------------------------------------------------
pro WID_SLIDER_LengthLimit_ChangeValue, Event

	id=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_LengthLimitValue')
	widget_control, id, SET_VALUE=STRCOMPRESS(string(Event.value*0.1))

end


; Name: WID_DROPLIST_InterpolationScheme_Change
; Subroutine that gets called when an event happens on the Interpolation Scheme Droplist.

pro WID_DROPLIST_InterpolationScheme_Change, Event

	common share_fibers_main

	id=widget_info(Event.top, FIND_BY_UNAME='WID_DROPLIST_InterpolationScheme')
	selected=widget_info(id, /DROPLIST_SELECT)

	CASE selected of
	0: interpol = '1' ; Linear interpolation
	1: interpol = '0'
	2: interpol = '2'
	ENDCASE

end


; Name: WID_DROPLIST_ROI_Selection_Change
; Subroutine that gets called when an event happens on the ROI selection Droplist.

pro WID_DROPLIST_ROI_Selection_Change, Event

	common share_fibers_main

	id=widget_info(Event.top, FIND_BY_UNAME='WID_DROPLIST_ROI_Selection')
	selected=widget_info(id, /DROPLIST_SELECT)

	CASE selected of
	0: gfile = 'Fract_Anisotropy'
	1: gfile = 'Average Diffusion'
	2: gfile = 'S0'
	3: gfile = 'Orientation'

	ENDCASE

end


; Name: WID_DROPLIST_Thresholding_Change
; Subroutine that gets called when an event happens on the Thresholding Scheme Droplist.

pro WID_DROPLIST_Thresholding_Change, Event

	common share_fibers_main

	id=widget_info(Event.top, FIND_BY_UNAME='WID_DROPLIST_Thresholding')
	selected=widget_info(id, /DROPLIST_SELECT)

	CASE selected of
	0: tfile = 'Fract_Anisotropy'
	1: tfile = 'Generalized_Anisotropy'

	ENDCASE

end


pro WID_TEXT_StepSizeValue_Event, Event

	if (Event.type eq 0) then begin
		if (Event.CH eq 10) then begin

		id_slider=widget_info(Event.top, FIND_BY_UNAME='WID_SLIDER_StepSize')
		id_textbox=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_StepSizeValue')

		widget_control, id_textbox, GET_VALUE=threshold_value
		widget_control, id_slider, SET_VALUE=STRCOMPRESS(string(threshold_value/0.1))
		endif

	endif

end


pro WID_TEXT_ThresholdValue_Event, Event

	;print, Event
	if (Event.type eq 0) then begin
		if (Event.CH eq 10) then begin

		id_slider=widget_info(Event.top, FIND_BY_UNAME='WID_SLIDER_Threshold')
		id_textbox=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_ThresholdValue')

		widget_control, id_textbox, GET_VALUE=threshold_value
		widget_control, id_slider, SET_VALUE=STRCOMPRESS(string(threshold_value/0.005))
		endif

	endif

end


pro WID_TEXT_SeedDensityValue_Event, Event

	if (Event.type eq 0) then begin
		if (Event.CH eq 10) then begin

		id_slider=widget_info(Event.top, FIND_BY_UNAME='WID_SLIDER_SeedDensity')
		id_textbox=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_SeedDensityValue')

		widget_control, id_textbox, GET_VALUE=threshold_value
		widget_control, id_slider, SET_VALUE=STRCOMPRESS(string(threshold_value))
		endif

	endif

end


pro WID_TEXT_LengthLimitValue_Event, Event

	if (Event.type eq 0) then begin
		if (Event.CH eq 10) then begin

		id_slider=widget_info(Event.top, FIND_BY_UNAME='WID_SLIDER_LengthLimit')
		id_textbox=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_LengthLimitValue')

		widget_control, id_textbox, GET_VALUE=threshold_value
		widget_control, id_slider, SET_VALUE=STRCOMPRESS(string(threshold_value/0.1))
		endif

	endif

end

pro WID_TEXT_AngleLimitValue_Event, Event

    if (Event.type eq 0) then begin
        if (Event.CH eq 10) then begin

        id_slider=widget_info(Event.top, FIND_BY_UNAME='WID_SLIDER_AngleLimit')
        id_textbox=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_AngleLimitValue')

        widget_control, id_textbox, GET_VALUE=threshold_value
        widget_control, id_slider, SET_VALUE=STRCOMPRESS(string(threshold_value))
        endif

    endif

end


pro WID_TEXT_VisualizationFactorValue_Event, Event

	if (Event.type eq 0) then begin
		if (Event.CH eq 10) then begin

		id_slider=widget_info(Event.top, FIND_BY_UNAME='WID_SLIDER_VisualizationFactor')
		id_textbox=widget_info(Event.top, FIND_BY_UNAME='WID_TEXT_VisualizationFactorValue')

		widget_control, id_textbox, GET_VALUE=threshold_value
		widget_control, id_slider, SET_VALUE=STRCOMPRESS(string(threshold_value))
		endif

	endif

end


pro WID_BASE_FIBERS_MAIN_EVENT, Event
end

;
; Empty stub procedure used for autoloading.
;
pro FibersGUI_eventcb
end
