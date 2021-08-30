;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved

; Subroutine name: mas_get_multi_slice_selection
; Created by: BT, 2008-06
; Calling Information:
;  event - event structure
; Bugs or Important Comments to Developers:

; Purpose of subroutine:
;
; event handler for mas_get_multi_slice_selected
;
; Editing Information:

pro mas_get_multi_slice_selection_event, event

    id = event.id

    name = widget_info(id, /uname)

    if (name eq 'btn_ok') then begin
        widget_control, event.top, /destroy
        return
    endif else if (name eq 'btn_cancel') then begin
        widget_control, event.top, get_uvalue=toggled
        (*toggled)[*] = 0 ;; zero everything out if canceled
        widget_control, event.top, /destroy
        return
    endif

    widget_control, event.top, get_uvalue=toggled

    if (ptr_valid(toggled)) then begin

        widget_control, event.id, get_uvalue=value

        if (n_elements(value) eq 0) then return
        if (value gt n_elements(*toggled)) then return

        ;; update the pointer's contents
        (*toggled)[value] = event.select

    endif

end

; Subroutine name: mas_get_multi_slice_selection
; Created by: BT, 2008-06
; Calling Information:
;  group_leader  - IN - id of group leader TLB.
;  count         - OUT- set this to a named variable to receive then
;                  number of slices chosen
; Bugs or Important Comments to Developers:

; Purpose of subroutine:
;
; Presents the user with a selection grid containing all of the slices
; in the current slice axis. Returns an array containing the indices
; of the slices within the state1 array. Returns -1 if no slices were selected

; Editing Information:

function mas_get_multi_slice_selection, group_leader, count=count

    common scan_data

    ci = project.ci
    slice_axis = project.procpramarray[ci].slice_axis

    case slice_axis of

        0: sdim = project.imndarray[ci].sdim * project.procpramarray[ci].slice_interp
        1: sdim = project.imndarray[ci].pdim * project.procpramarray[ci].phase_interp
        2: sdim = project.imndarray[ci].fdim * project.procpramarray[ci].freq_interp
        else: begin
            print, 'mas_get_multi_slice_selection: unknown slice axis:'+string(slice_axis)
            count = 0
            return, -1
        end
    endcase

    ;; try to be as square as possible
    x_width = floor(sqrt(sdim))
    y_width = ceil(1.0 * sdim/x_width)

    base = widget_base(title='Batch select slices', /modal, group_leader=group_leader, row=2)

    bt_base = widget_base(base, row=y_width)

    iter = 0

    for y = 0, y_width-1 do begin

        x_base = widget_base(bt_base, /nonexclusive, /row, /grid_layout)

        for x = 0, x_width-1 do begin

            ;; set the uvalue of the button to the scan index
            btn = widget_button(x_base, value=string(iter, format='(I3.3)'), uvalue=iter)

            if (++iter gt sdim) then break ;; no more slices

        endfor

    endfor

    ok_base = widget_base(base, /align_center, column=2)
    btn_ok = widget_button(ok_base, value="Select checked slices", uname='btn_ok', /align_center)
    btn_cancel = widget_button(ok_base, value="Cancel", uname='btn_cancel', /align_center)

    ;; the contents of this pointer will be altered in the event
    ;; handler according to the user's preference
    toggled = ptr_new(bytarr(sdim), /no_copy)
    widget_control, base, set_uvalue=toggled

    ;; this will block since the tlb is modal
    widget_control, base, /realize
    xmanager, 'mas_get_multi_slice_selection', base

    ;; at this point the OK or cancel button has been pressed.
    selected = where(*toggled ne 0, count)

    ptr_free, toggled

    return, count ne 0 ? selected : -1

end

; Subroutine name: mas_get_text_dialog_event
; Created by: BT, 2008-06
; Calling Information:
; event - event structure
; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; event handler for the "rename tab" dialog box
; Editing Information:

pro mas_get_text_dialog_event, event

   name = widget_info(event.id, /uname)

   case name of

       'btn_ok': begin
           text_id = widget_info(event.top, find_by_uname='text')

           widget_control, text_id, get_value=txt
           widget_control, event.top, get_uvalue=new_text

           (*new_text)[0] = txt
       end
       'btn_cancel': begin

       end
       'text': begin
           event_name = tag_names(event, /structure_name)
           if (event_name eq 'WIDGET_TEXT_CH') then begin
               if (event.ch eq 10) then begin
                   widget_control, event.id,  get_value=txt
                   widget_control, event.top, get_uvalue=new_text
                   (*new_text)[0] = txt
               endif
           endif
       end

       else: begin
           print, 'mas_get_text_dialog: unknown event from: '+name
           help, event, /structure
       end
   endcase

   widget_control, event.top, /destroy

end

; Subroutine name: mas_get_text_dialog
; Created by: BT, 2008-06
; Calling Information:
; tlb   -  the "group leader" tlb. Must exist.
; initial_text - IN - text to populate the text field. default: <empty>
; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; Display a dialog to the user asking for text input and returns the
; text string entered

; Editing Information:

function mas_get_text_dialog, $
                              tlb, $
                              title=title, $
                              initial_text=initial_text, $
                              message_text=message_text, $
                              text_size=text_size

    common common_widgets
    compile_opt hidden

    if (not keyword_set(title) or n_elements(title) eq 0) then begin
        title = 'Text Entry'
    endif
    
    if (keyword_set(initial_text) and n_elements(initial_text) ne 0) then begin
        display_text = initial_text
    endif else begin
        display_text = ''
    endelse

    if (keyword_set(message_text) and n_elements(message_text) ne 0) then begin
        display_message = message_text
    endif else begin
        display_message = 'Please enter text:'
    endelse

    if (keyword_set(text_size) and n_elements(text_size) ne 0) then begin
        xsize=text_size < 120
    endif else begin
        xsize=30
    endelse

    if (n_elements(tlb) eq 0) then begin
        tlb = WID_BASE_MAIN
    endif
    
    ;; this travels with the dialog. the event handler will put the
    ;; user's text in this pointer, the contents of which will be returned
    new_value = ptr_new(strarr(1))

    dialog_base = widget_base(title=title, /column, /modal, group_leader=tlb, uvalue=new_value)

    message_base = widget_base(dialog_base, /align_center)
    void = widget_label(message_base, value=message_text)
    entry_base = widget_base(dialog_base, /align_center)
    txt = widget_text(entry_base, value=initial_text, uname='text', xsize=xsize, /editable, /sensitive, /align_center)
    btn_base = widget_base(dialog_base, /row, /grid, /align_center)
    btn_cancel = widget_button(btn_base, value="Cancel", uname='btn_cancel', /align_center)
    btn_ok = widget_button(btn_base, value="Ok", uname='btn_ok', /align_center)

    widget_control, dialog_base, /realize
    xmanager, 'mas_get_text_dialog', dialog_base

    new_str = strcompress((*new_value)[0])
    ptr_free, new_value
    return, new_str

end

; Subroutine name: mas_zoom_tab_callback
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:


; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting


pro mas_zoom_tab_callback, event
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets
    CI = project.ci

    ;if the event was an interpolation change then set state 1 and state 2 to false

    if event.id eq freq_interp_slider then begin
        WIDGET_CONTROL , freq_interp_slider , GET_VALUE=temp
        if (temp ne project.procPramArray[CI].freq_interp) then begin
            project.procPramArray[CI].freq_interp =temp
            project.procPramArray[CI].state_2 = 0
            project.procPramArray[CI].state_1 = 0
            project.procpramArray[CI].adt_proccess_flag = 0
        endif
    end
    if event.id eq phase_interp_slider then begin

        WIDGET_CONTROL , phase_interp_slider , GET_VALUE=temp
        if (temp ne project.procPramArray[CI].Phase_interp) then begin
            project.procPramArray[CI].Phase_interp  =temp
            project.procPramArray[CI].state_2 = 0
            project.procPramArray[CI].state_1 = 0
            project.procpramArray[CI].adt_proccess_flag = 0
        endif
    end

    if event.id eq slice_interp_slider then begin

        WIDGET_CONTROL , slice_interp_slider , GET_VALUE=temp
        if (temp ne project.procPramArray[CI].Slice_interp) then begin
            project.procPramArray[CI].Slice_interp = temp
            project.procPramArray[CI].state_2 = 0
            project.procPramArray[CI].state_1 = 0
            project.procpramArray[CI].adt_proccess_flag = 0
        endif
    end

    if event.id eq x_zoom_slider or event.id eq y_zoom_slider  then begin

        WIDGET_CONTROL , x_zoom_slider , GET_VALUE=temp
        project.procPramArray[CI].x_zoom  =temp


        WIDGET_CONTROL , y_zoom_slider , GET_VALUE=temp
        project.procPramArray[CI].y_zoom =temp

        project.procPramArray[CI].state_2 = 0

    end

    if event.id eq zpad_chk_box then begin
       
       project.procPramArray[CI].zpad_flag  = event.select  
       project.procPramArray[CI].adt_proccess_flag  = 0                    
       project.procPramArray[CI].state_1 = 0
       project.procPramArray[CI].state_2 = 0

    end
    
    if event.id eq lock_interp_chk_box then begin
        project.procPramArray[CI].lock_interp =event.select
    end

	if event.id eq down_sample_chk_box then begin
        project.procPramArray[CI].down_sample =event.select
        project.procPramArray[CI].state_2 = 0
        project.procPramArray[CI].state_1 = 0
        project.procpramArray[CI].adt_proccess_flag = 0
        print, project.procPramArray[CI].down_sample
    end

    ;if the images are locked then set them to the same value
    if project.procPramArray[CI].lock_interp eq 1 then begin
        WIDGET_CONTROL , freq_interp_slider , GET_VALUE=temp
        project.procPramArray[CI].freq_interp  = temp
        project.procPramArray[CI].Phase_interp = temp
        project.procPramArray[CI].Slice_interp = temp

    end


    if project.imndArray[CI].sdim eq 1 and project.procPramArray[CI].Slice_interp ge 2 then begin
        void = dialog_message(['Interpolation in the Slice dimension is not possible',$
                               'when the number of slices is 1.',$
                               'Setting the slice interpolation to 1',$
                               'Removing the Interpolation Lock'],/INFORMATION, /center )
        project.procPramArray[CI].Slice_interp = 1
        project.procPramArray[CI].lock_interp = 0


    end

    mas_redraw_GUI

    if project.auto_Display eq 1 then mas_display

end


; Subroutine name: mas_settings_event
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:


; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

PRO mas_settings_event, Event

    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets
    CI = project.ci

    case Event.id of


       iImage_flag_button: begin
        project.iImage_flag = event.SELECT
        project.procPramArray[CI].state_2 = 0
       end

       p_dim_shift :BEGIN
			;These are here to troubleshoot
	;       	PRINT, 'Event structure type: ', TAG_NAMES(EVENT, /STRUCTURE_NAME)
	;       	PRINT, 'Event Value: ', event.value
	;
	         project.imndArray[CI].pdim_shift = event.value
	         project.procPramArray[CI].state_1 = 0
	         project.procPramArray[CI].state_2 = 0
	         ;mas_redraw_GUI
	; We should only redraw the actual widget and not the whole display

	         widget_control    ,p_dim_shift,SENSITIVE=1   $
	            ,SET_SLIDER_Min=-project.imndArray[ci].pdim $
	            ,SET_SLIDER_MAX=project.imndArray[ci].pdim $
	            ,SET_VALUE=project.imndArray[ci].pdim_shift

	         if project.auto_Display eq 1 then mas_display

       end

       s_dim_shift :BEGIN
	         project.imndArray[CI].sdim_shift = event.value
	         project.procPramArray[CI].state_1 = 0
	         project.procPramArray[CI].state_2 = 0

			;HS 20061221
			;Don't need to redraw every time the slider is moved.
	         ;mas_redraw_GUI

			widget_control    ,s_dim_shift,SENSITIVE=1   $
	            ,SET_SLIDER_Min=-project.imndArray[ci].sdim $
	            ,SET_SLIDER_MAX=project.imndArray[ci].sdim $
	            ,SET_VALUE=project.imndArray[ci].sdim_shift

	         if project.auto_Display eq 1 then mas_display

       end

	       f_dim_shift :BEGIN
	         project.imndArray[CI].fdim_shift = event.value
	         project.procPramArray[CI].state_1 = 0
	         project.procPramArray[CI].state_2 = 0

			;HS 20061221
			;Don't need to redraw every time the slider is moved.
	         ;mas_redraw_GUI

			widget_control    ,f_dim_shift,SENSITIVE=1   $
		            ,SET_SLIDER_Min=-project.imndArray[ci].fdim $
		            ,SET_SLIDER_MAX=project.imndArray[ci].fdim $
		            ,SET_VALUE=project.imndArray[ci].fdim_shift


         if project.auto_Display eq 1 then mas_display

       end

       TIME_ZERO_SLIDER : begin
         project.imndArray[ci].time.zero = event.value
       end
       SLIDER_PREIMAGES : begin
         project.imndArray[ci].n_pre = event.value
       end
       n_avg_flag_button : begin
        project.procPramArray[CI].n_avg_flag = event.SELECT
         project.procPramArray[CI].state_1 = 0
         project.procPramArray[CI].state_2 = 0
         mas_redraw_GUI
       end
      phi_unwrap_flag_button : begin
        project.procPramArray[CI].phi_unwrap_flag = event.SELECT
        project.procPramArray[CI].state_1 = 0
        project.procPramArray[CI].state_2 = 0
        mas_redraw_GUI
        end
        
       kspace_fdim_shift :BEGIN
           project.imndArray[CI].k_fdim_shift = event.value
           project.procPramArray[CI].state_1 = 0
           project.procPramArray[CI].state_2 = 0

      ;HS 20061221
      ;Don't need to redraw every time the slider is moved.
           ;mas_redraw_GUI

      widget_control    ,kspace_fdim_shift,SENSITIVE=1   $
                ,SET_SLIDER_Min=-project.imndArray[ci].fdim $
                ,SET_SLIDER_MAX=project.imndArray[ci].fdim $
                ,SET_VALUE=project.imndArray[ci].k_fdim_shift


         if project.auto_Display eq 1 then mas_display

       end
       
        kspace_pdim_shift :BEGIN
           project.imndArray[CI].k_pdim_shift = event.value
           project.procPramArray[CI].state_1 = 0
           project.procPramArray[CI].state_2 = 0

      ;HS 20061221
      ;Don't need to redraw every time the slider is moved.
           ;mas_redraw_GUI

      widget_control    ,kspace_pdim_shift,SENSITIVE=1   $
                ,SET_SLIDER_Min=-project.imndArray[ci].pdim $
                ,SET_SLIDER_MAX=project.imndArray[ci].pdim $
                ,SET_VALUE=project.imndArray[ci].k_pdim_shift


         if project.auto_Display eq 1 then mas_display

       end
       
       kspace_sdim_shift :BEGIN
           project.imndArray[CI].k_sdim_shift = event.value
           project.procPramArray[CI].state_1 = 0
           project.procPramArray[CI].state_2 = 0

      ;HS 20061221
      ;Don't need to redraw every time the slider is moved.
           ;mas_redraw_GUI

      widget_control    ,kspace_sdim_shift,SENSITIVE=1   $
                ,SET_SLIDER_Min=-project.imndArray[ci].sdim $
                ,SET_SLIDER_MAX=project.imndArray[ci].sdim $
                ,SET_VALUE=project.imndArray[ci].k_sdim_shift


         if project.auto_Display eq 1 then mas_display

       end
       
       kspace_fdim_subsamp: begin
           project.procPramArray[CI].state_1 = 0
           project.procPramArray[CI].state_2 = 0
           project.imndarray[CI].k_fdim_span = (event.value > 1) < (project.imndarray[CI].k_fdim_span_max)
           project.imndarray[ci].f_voxsz = project.imndarray[ci].f_fov/project.imndarray[CI].k_fdim_span
       end
       
       kspace_pdim_subsamp: begin
           project.procPramArray[CI].state_1 = 0
           project.procPramArray[CI].state_2 = 0
           project.imndarray[CI].k_pdim_span = (event.value > 1) < (project.imndarray[CI].k_pdim_span_max)
           project.imndarray[ci].p_voxsz = project.imndarray[ci].p_fov/project.imndarray[CI].k_pdim_span
       end
       
       kspace_sdim_subsamp: begin
           project.procPramArray[CI].state_1 = 0
           project.procPramArray[CI].state_2 = 0
           project.imndarray[CI].k_sdim_span = (event.value > 1) < (project.imndarray[CI].k_sdim_span_max)
           project.imndarray[ci].s_voxsz = project.imndarray[ci].s_fov/project.imndarray[CI].k_sdim_span
       end
       
    endcase
end


; Subroutine name: auto_display_toggle
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; This changes the flag for auto_Display and changes the message on the button

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

pro auto_display_toggle, event
    COMMON scan_data
    project.auto_Display = event.SELECT


end


; Subroutine name: set_flip_Direction
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; This updates the status bar on the bottom of the main display.

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting


PRO set_flip_Direction , direction
    COMPILE_OPT IDL2

    COMMON scan_data
    ;print, direction
    CI = project.ci
    project.procPramArray[CI].flip_Direction = direction
    project.procPramArray[CI].state_2 = 0

    ;mas_redraw_GUI
    if project.auto_Display eq 1 then mas_display

end


;***************************************
; Subroutine name: set_smooth_Direction
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; This updates the status bar on the bottom of the main display for the smooth part.

; Editing Information:
    ;Created by AMV 2006/12/1.
    ;Created the subroutine


PRO set_smooth_Direction , direction
    COMPILE_OPT IDL2

    COMMON scan_data
    ;print, direction
    CI = project.ci
    project.procPramArray[CI].smooth_Direction = direction
    project.procPramArray[CI].state_1 = 0
    project.procPramArray[CI].state_2 = 0

    ;mas_redraw_GUI
    if project.auto_Display eq 1 then mas_display

end





; Subroutine name: set_signal_type
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:


; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting
    ; Added optional keyword to specify the scan index (Magdoom 05/14/2015)
    
PRO set_signal_type , type, ci = ci
    COMPILE_OPT IDL2

    COMMON scan_data
    COMMON common_widgets

    if ~keyword_set(ci) then CI = project.ci
    project.procPramArray[CI].signal_type = type
    project.procPramArray[CI].state_1 = 0
    project.procPramArray[CI].state_2 = 0

    mas_redraw_GUI
    if project.auto_Display eq 1 then mas_display

end


; Subroutine name: toggle_motion_correction
; Created by: Chad Durgin
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:


; Editing Information:


;; project.use_motion_correction constants:
;;    0   -  All Off
;;    1   -  Edge On
;;    2   -  MI On

PRO toggle_motion_correction, TYPE=type
    COMPILE_OPT IDL2

    COMMON scan_data
    COMMON common_widgets

    ci = project.ci

    project.procPramArray[ci].mc_enable = type
    project.procPramArray[ci].state_1 = 0
    project.procPramArray[ci].state_2 = 0

    mas_redraw_GUI

end


; Subroutine name: turnoff_motion_correction
; Created by: Hector Sepulveda
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
	; Turning off the motion correction toggle whenever a new scan is loaded
; Editing Information:


PRO turnoff_motion_correction
    COMPILE_OPT IDL2

    COMMON scan_data
    COMMON common_widgets

    project.procpramArray[project.ci].mc_enable = 0

    mas_redraw_GUI
end

; Subroutine name: set_slice_axis
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:


; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

PRO set_slice_axis , direction
    COMPILE_OPT IDL2

    COMMON scan_data
    COMMON common_widgets

    CI = project.ci
    project.procPramArray[CI].slice_axis = direction
    project.procPramArray[CI].state_2 = 0


    mas_redraw_GUI
    if project.auto_Display eq 1 then mas_display

end


; Subroutine name: set_rotate_Direction
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:


; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

PRO set_rotate_Direction , direction
    COMPILE_OPT IDL2

    COMMON scan_data
    ;print, direction
    CI = project.ci
    project.procPramArray[CI].rotate_Direction = direction
    project.procPramArray[CI].state_2 = 0

    ;mas_redraw_GUI
    if project.auto_Display eq 1 then mas_display
end


; Subroutine name: update_status_bar
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; This updates the status bar on the bottom of the main display

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

PRO update_status_bar , text
    COMPILE_OPT IDL2

    COMMON common_widgets
    widget_control, LABEL_STATUS_BAR , SET_VALUE=text
end


; Subroutine name: single_Multi_flag_toggle
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; This changes the flag for sliderControlDisplay and changes the message on the button.

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

pro single_Multi_flag_toggle

    COMPILE_OPT IDL2

    COMMON scan_data
    COMMON common_widgets

    ;if there is no open scan then dont do any thing
    if project.scan_Open_Flag eq 0 then return
    CI = project.ci

    ;1 is for multiple and 0 is for single
    if project.procPramArray[CI].single_Multi_flag eq 1 then begin
       project.procPramArray[CI].single_Multi_flag = 0
       widget_control, SLIDER_DISPLAY_CONTROL , SET_VALUE='Single' ,SENSITIVE=1
    endif else begin
       project.procPramArray[CI].single_Multi_flag  = 1
       widget_control, SLIDER_DISPLAY_CONTROL , SET_VALUE='Multiple' ,SENSITIVE=1
    endelse

    project.procPramArray[CI].state_2 = 0
    mas_redraw_GUI
    if project.auto_Display eq 1 then mas_display
end

; Subroutine name: sort_dir_toggle
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; This changes the flag for sort_dir and changes the message on the button

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

pro sort_dir_toggle
    COMPILE_OPT IDL2

    COMMON scan_data
    COMMON common_widgets
    ;if there is no open scan then dont do any thing
    if project.scan_Open_Flag eq 0 then return

    ;1 is for adim and 0 is for sdim

    CI = project.ci
    if project.procPramArray[CI].sort_dir eq 0 then begin
       project.procPramArray[CI].sort_dir = 1
       widget_control, SLIDER_DIMENSION_CONTROL , SET_VALUE='  Adim  ',SENSITIVE=1
    endif else begin
       project.procPramArray[CI].sort_dir = 0
       widget_control, SLIDER_DIMENSION_CONTROL , SET_VALUE='  Sdim  ',SENSITIVE=1
    endelse

    project.procPramArray[CI].state_2 = 0
    mas_redraw_GUI
    if project.auto_Display eq 1 then mas_display

end

; Subroutine name: display_threshold_callback
; Created by: wtriplett
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; Handles changes to the overall display scale min/max thresholds

; Editing Information:

pro display_threshold_callback, event

    common scan_data, project
    common common_widgets
    
    uname = widget_info(event.id, /uname)

    case uname of
    
        'DISPLAY_SCALE_RESET_BUTTON': begin
            min_txt_box = widget_info(event.top, find_by_uname='DISPLAY_DATA_MIN')
            max_txt_box = widget_info(event.top, find_by_uname='DISPLAY_DATA_MAX')
            
            min_data_value = min(*project.dataarray[project.ci].state1, max=max_data_value) 
            
            ;; convert to float from compound values (like complex  if applicable)
            min_data_value = float(min_data_value)
            max_data_value = float(max_data_value)
            
            project.procpramarray[project.ci].min_display_thresh = min_data_value
            project.procpramarray[project.ci].max_display_thresh = max_data_value
            
            widget_control, min_txt_box, set_value=strtrim(string(min_data_value),2)
            widget_control, max_txt_box, set_value=strtrim(string(max_data_value),2)        
        end
        'DISPLAY_DATA_MIN': begin
            widget_control, event.id, get_value=value
            project.procpramarray[project.ci].min_display_thresh = float(value)
        end
        
        'DISPLAY_DATA_MAX': begin
            widget_control, event.id, get_value=value
            project.procpramarray[project.ci].max_display_thresh = float(value)
        end
        
        else: message, /info, "Rec'd event from untracked widget: "+uname

    endcase
    
end

; Subroutine name: intensity_slider_callback
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; This takes the values out of the sliders and updates the intensityMax, intensityCen, intensityMin

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

pro intensity_slider_callback
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets
    CI = project.ci

    ;since we change a parameter belonging to state 2 then we need to clear state 2 so that it is reprocessed
    project.procPramArray[CI].state_2 = 0

    widget_control, INTENSITY_SLIDER_MAX , GET_VALUE=temp
    project.procPramArray[CI].intensity_Max = temp

    widget_control, INTENSITY_SLIDER , GET_VALUE=temp
    project.procPramArray[CI].intensity_Cen = temp

    widget_control, INTENSITY_SLIDER_MIN , GET_VALUE=temp
    project.procPramArray[CI].intensity_Min = temp

    if project.auto_Display eq 1 then mas_display

end


; Subroutine name: adim_slider_callback
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; This takes the values out of the sliders and updates the adimStartSlice

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

pro adim_slider_callback
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets

    WIDGET_CONTROL , ADIM_SLIDER , GET_VALUE=temp
    CI = project.CI
    ;if there is no open scan then dont do any thing
    if project.scan_Open_Flag eq 0 then return
    ;since we chage a parameter belonging to state 2 then we need to clear state 2 so that it is reprocessed

    project.procPramArray[CI].state_2 = 0

    project.procPramArray[CI].adim_Start = (temp) > 0

    widget_control, ADIM_SLIDER, SET_VALUE = project.procPramArray[CI].adim_start $
      , SET_SLIDER_MAX=project.imndArray[CI].adim-1

    if project.auto_Display eq 1 then mas_display
end


; Subroutine name: sdim_slider_callback
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; This takes the values out of the sliders and updates the sdim_Start

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

pro sdim_slider_callback
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets

    WIDGET_CONTROL , SDIM_SLIDER , GET_VALUE=slider_value
    CI = project.ci

    ;if there is no open scan then dont do any thing
    IF project.scan_Open_Flag EQ 0 THEN RETURN

    ;since we chage a parameter belonging to state 2 then we need to clear state 2 so that it is reprocessed
    project.procPramArray[CI].state_2 = 0

    ;depending on what dimension we are dealing with we need to store the appropriate value in
    ; the structure. [0|1|2]  0=freq phase, 1=freq slice 2=phase slice

	slice_axis = project.procPramArray[CI].slice_axis
  is_3d = project.imndarray[CI].dimensions eq 3
  
	CASE slice_axis OF
		0: BEGIN
			 project.procPramArray[CI].sdim_start = slider_value ;traverse the sdim
			 if isa(project.imndArray[CI].dicomfiles) then slidermax = project.imndArray[CI].sdim * project.procPramArray[CI].Slice_interp*(1+is_3d*project.procpramArray[CI].zpad_flag) - 1 else $
			   slidermax = project.imndArray[CI].k_sdim_span * project.procPramArray[CI].Slice_interp*(1+is_3d*project.procpramArray[CI].zpad_flag) - 1
   	   ;if the sdim slider is greater then the slider max then cut back on the sdim
   	   if project.procPramArray[CI].sdim_start gt slidermax then project.procPramArray[CI].sdim_start = slidermax
			 sliderval = project.procPramArray[CI].sdim_start
			 sliderlab = 'Slice'
       END

		1: BEGIN
			 project.procPramArray[CI].pdim_start = slider_value ;traverse the pdim
       if isa(project.imndArray[CI].dicomfiles) then slidermax = project.imndArray[CI].pdim * project.procPramArray[CI].Phase_interp*(1+is_3d*project.procpramArray[CI].zpad_flag) - 1 else $
        slidermax = project.imndArray[CI].k_pdim_span * project.procPramArray[CI].Phase_interp*(1+is_3d*project.procpramArray[CI].zpad_flag) - 1
       ;if the sdim slider is greater then the slider max then cut back on the sdim
       if project.procPramArray[CI].pdim_start gt slidermax then project.procPramArray[CI].pdim_start = slidermax
       sliderval = project.procPramArray[CI].pdim_start
       sliderlab = 'Phase'
       END

		2: BEGIN
			 project.procPramArray[CI].fdim_start = slider_value ;traverse the fdim
       if isa(project.imndArray[CI].dicomfiles) then slidermax = project.imndArray[CI].fdim * project.procPramArray[CI].freq_interp*(1+is_3d*project.procpramArray[CI].zpad_flag) - 1 else $
        slidermax = project.imndArray[CI].k_fdim_span * project.procPramArray[CI].freq_interp*(1+is_3d*project.procpramArray[CI].zpad_flag) - 1
       ;if the sdim slider is greater then the slider max then cut back on the sdim
   		 if project.procPramArray[CI].fdim_start gt slidermax then project.procPramArray[CI].fdim_start = slidermax
			 sliderval = project.procPramArray[CI].fdim_start
			 sliderlab = 'Freq'
       END

	ENDCASE

	if project.imndArray[CI].sdim le 1 then sensitivity = 0 else sensitivity = 1

   	widget_control, SDIM_SLIDER_MAX, SET_VALUE= STRTRIM( slidermax,1), SENSITIVE=sensitivity
   	widget_control, SDIM_SLIDER, SET_VALUE = sliderval, SET_SLIDER_MAX=slidermax, SENSITIVE=sensitivity
   	widget_control, SDIM_SLIDER_TITLE, SET_VALUE=slicerlab, SENSITIVE=sensitivity

    if project.auto_Display eq 1 then mas_display
    project.procPramArray[CI].image_fit_flag=0
end


; Subroutine name: update_scan_selection_event
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting


pro update_scan_selection_event, Event
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets

    ;whatever you click "assiant" (HS 20061003: I have no idea what he means) that to ci
    ;update label looks at the ci  and then displays that entirely.

    project.ci = Event.INDEX
    mas_redraw_GUI
    ;Widget_control, WID_SCAN_LIST, SET_LIST_SELECT=project.ci
    if project.auto_display eq 1 then mas_display
end


; Subroutine name: update_scan_selection
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting

PRO update_scan_selection
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets
    newListStringArray = strarr(project.ni)
    newListStringArray = project.scan_list[0:project.ni-1]
    Widget_control, WID_SCAN_LIST , SET_VALUE= newListStringArray ,SENSITIVE=1
    mas_redraw_GUI

END


; Subroutine name: mas_add_current_scan_to_scan_selection
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:

;=========================
    ;summary
    ; we have just added a new scan to the project and we need to manage
    ;  the scan selection list widget by adding a new scan. since the widget
    ;  wont allow us to read the data from the value directly we need to
    ;  maintain a global array to hold the values to update to the array.
    ;==========================
    ;algorithm
    ; 1) make a new string array that is 1 bigger than than the current list
    ;      string array. call this new array newListStringArray
    ; 2) copy all the old values to the new list string array. the old array
    ;      is maintained by the mas_extract_imnd. the only value that needs
    ;         to be updated it the value displayed on the widget. so we will assign
    ;   this later
    ; 3) after we have this new array compiled from all across the imnd array
    ;   and located in one place we just store the new array to old
    ;      array location. project.listStringArray
    ; 4) after we hace backed up the data in step 3 we can now assign the value
    ;   directly to the WID_SCAN_LIST's value
    ;==========================
    ;note
    ;  since the project.listStringArray has a cap of 50 items in the array we
    ;  need pass only what is needed to the widget so it doesn't have 50 empty
    ;  entries displayed and that are selectable.
    ;==========================


; Editing Information:
    ;Edited by HS 2006/10/03.
    ;Fix spelling mistakes and commenting
    ;Edited by CD 2007/06/28
    ;Use displayname is set, otherwise display path

pro mas_add_current_scan_to_scan_selection
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets

    ci = project.ci
    display_name = project.imndArray[ci].display_Name
    if (display_name eq '') then begin
        display_name = project.imndArray[ci].file_Path
    endif

;    new_name = mas_get_text_dialog(WID_BASE_MAIN, $
;                                   initial_text=display_name, $
;                                   message_text='Please enter a name for this scan:', $
;                                   text_size=100)
;    if (new_name ne '') then begin
;        if (new_name ne display_name) then begin
;            project.imndArray[ci].display_Name = new_name + ' ['+display_name+']'
;        endif
;    endif else begin
;        project.imndArray[ci].display_Name = display_name
;    endelse

    ;1
    newListStringArray = strarr(project.ni)

    ;2
    FOR ii=0, project.ni-1 DO BEGIN
    	IF STRCMP(project.imndArray[ii].display_Name,'') EQ 1 THEN BEGIN

            newListStringArray[ii] = project.imndArray[ii].file_Path

        ENDIF ELSE BEGIN

            newListStringArray[ii] = project.imndArray[ii].display_Name

        ENDELSE
    ENDFOR

    ;3
    project.scan_list = newListStringArray

    ;4
    Widget_control, WID_SCAN_LIST , SET_VALUE= newListStringArray ,SENSITIVE=1
    mas_redraw_GUI
END


;*****************************************************************************************************
;
; NAME:
;   mas_callback
;
; PURPOSE:
;   most events that are generated by the main mas display can be handled here if they are small.
;
; ARGUMENTS:
;
;
; MODIFICATION HISTORY:
;
;
;*****************************************************************************************************
pro mas_callback
    COMPILE_OPT IDL2
end
