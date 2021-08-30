;+
; :Description:
;    Objective function for autophase.
;
; :Params:
;    param_array - this array contains the parameters (0-order, x1, y1, etc.) to 
;                  compute the phase correction.
;                  
; :Author: wtriplett
;-
function mas_image_phase_autophase_objective, param_array

    common mas_image_phase_powell, slice
    
    params = { phase_0:  param_array[0], $
               pivot_x:  param_array[1]*(size(slice,/dim))[0], $
               phase_x1: param_array[2], $
               pivot_y:  param_array[3]*(size(slice,/dim))[1], $
               phase_y1: param_array[4] }

    temp = slice
    mas_image_phase_apply_correction, temp, params

    return, total(abs(imaginary(temp)));; - 0.001*total(real_part(temp))

end

;+
; :Description:
;    Attempts to autophase the image set based on the selected array element.
;    The phase parameters found for that array element will be applied to each
;    array element in turn for a selected slice.
;
; :Params:
;    state - widget state pointer
;
; :Keywords:
;    single_slice - set this keyword to a slice to perform the auto phase for a single subject.
;
; :Author: wtriplett
;-
pro mas_image_phase_autophase, state, single_slice=single_slice

    common scan_data, project
    common mas_image_phase_powell, slice
    
    ci = project.ci
    
    if (keyword_set(single_slice)) then begin
        slice_start = single_slice
        slice_end   = single_slice
    endif else begin
        slice_start = 0
        slice_end   = project.imndarray[ci].sdim-1
    endelse
    
    widget_control, (*state).sl_aid, get_value = a_num
    
    progressbar = obj_new('progressbar', color='red', text='Autophasing...')
    progressbar->Start
    
    ;; number of times to iterate Powell for improved result.
    max_iter = 5
        
    for s = slice_start, slice_end do begin 
        
        wait, 0.0001 ;; needed for some reason so that the progress bar will update
        progressbar->Update, float(s)/float((slice_end-slice_start)) * 100.0
        
        for n = 0, max_iter-1 do begin

            slice = reform((*project.dataarray[ci].state1)[*,*,s,a_num])

            Xi    = diag_matrix([1.0,1.0,1.0,1.0,1.0])
            start = [ 0.0, 0.5, 0.0, 0.5, 0.0 ]
            powell, start, Xi, 1e-7, fmin, 'mas_image_phase_autophase_objective'
    
            params = { phase_0: start[0], $
                       pivot_x: start[1]*(size(slice,/dim))[0], $
                       phase_x1: start[2], $
                       pivot_y: start[3]*(size(slice,/dim))[1], $
                       phase_y1: start[4] }
        
            for a = 0, project.imndarray[ci].adim-1 do begin
            
                temp = reform((*project.dataarray[ci].state1)[*,*,s,a])

                mas_image_phase_apply_correction, temp, params

                ;; correct for negative maxima in real part
                if (n+1 eq max_iter and a eq 0 and total(real_part(temp)) lt 0) then begin
                    temp = reform((*project.dataarray[ci].state1)[*,*,s,a])
                    params.phase_0 = params.phase_0 + 0.5
                    mas_image_phase_apply_correction, temp, params
                endif
                
                (*project.dataarray[ci].state1)[*,*,s,a] = temp
        
            endfor

            ;; Reset these to defaults, since we want the user to be
            ;; able to fine tune the powell results 
            (*state).phase_corr_params[s,*].phase_0  = 0
            (*state).phase_corr_params[s,*].phase_x1 = 0
            (*state).phase_corr_params[s,*].pivot_x  = project.imndarray[ci].fdim/2
            (*state).phase_corr_params[s,*].phase_y1 = 0
            (*state).phase_corr_params[s,*].pivot_y  = project.imndarray[ci].pdim/2

        endfor

    endfor
    
    progressbar->Destroy
    
    project.procpramarray[0].state_2 = 0
    
end

;+
; :Description:
;    Refreshes the surface plots and image windows.
;
; :Params:
;    state  - widget state pointer
;    slice  - 2D slice information to draw
;
;
; :Author: wtriplett
;-
pro mas_image_phase_refresh_display, state, slice


    if (n_elements(slice) eq 0) then begin
        widget_control, (*state).sl_aid, get_value   = a_num
        widget_control, (*state).sl_slice, get_value = s_num
        
        pdata = (*state).pdata
        slice = (*pdata)[*,*,s_num,a_num]
        mas_image_phase_apply_correction, slice, (*state).phase_corr_params[s_num,a_num]
    endif

    ;; show results
    az = (*state).rotation_az
    ax = (*state).rotation_ax
    zrange = (*state).data_range * [0.95,1.05]
    wset, (*state).surf_window
    !P.multi = [0, 2, 1]
    xstyle = 1+4+8
    ystyle = 1+4+8
    zstyle = 1+8
    !Z.charsize = 1.35
    !Z.minor = 1
    if (1) then begin
        shade_surf, real_part(slice), zrange=zrange, az=az, ax=ax, xstyle=xstyle, ystyle=ystyle, zstyle=zstyle
        shade_surf, imaginary(slice), zrange=zrange, az=az, ax=ax, xstyle=xstyle, ystyle=ystyle, zstyle=zstyle
    endif else begin
        surface, real_part(slice), zrange=zrange, az=az, ax=ax, xstyle=xstyle, ystyle=ystyle, zstyle=zstyle
        surface, imaginary(slice), zrange=zrange, az=az, ax=ax, xstyle=xstyle, ystyle=ystyle, zstyle=zstyle
    endelse
    !Z.minor = 0
    !Z.gridstyle = 0
    total_imag = total(abs(imaginary(slice)))
    xyouts, 0.75, 0.95, string(total_imag, format='(%"Total Imag: %0.3f")'), /normal, align=0.5, charsize=1.5
    xyouts, 0.25, 0.95, string(total(abs(real_part(slice))), format='(%"Total Real: %0.3f")'), /normal, align=0.5, charsize=1.5
    !P.multi = 0
    
    wset, (*state).img_window
    sz_slice = size(slice, /dimensions)
    xsize = sz_slice[0] * 175.0/sz_slice[0]
    ysize = sz_slice[1] * 175.0/sz_slice[0]
    tvscl, congrid([real_part(slice), imaginary(slice)], 2*xsize, ysize, /interp)

end
;+
; :Description:
;    Applies the correction parameters to a single slice. Does not replace 
;    any of the main data. Alters the slice data in place
;
; :Params:
;    slice  - 2D slice data
;    params - structure containing the phase parameters
;
;
;
; :Author: wtriplett
;-
pro mas_image_phase_apply_correction, slice, params

    phase_0  = params.phase_0
    pivot_x  = params.pivot_x
    phase_x1 = params.phase_x1
    pivot_y  = params.pivot_y
    phase_y1 = params.phase_y1

    dims = size(slice, /dimensions)
    xdim = dims[0]
    ydim = dims[1]
    
    two_pi_i = 2*!PI*complex(0,1)
    
    ;; scale phase ramp 0..1 in image space
    index = (findgen(xdim) - pivot_x)/xdim
    
    ;; phase = X1 * ramp + X0
    arg = ( float(phase_x1)*index + float(phase_0) )
    phase_x = exp(two_pi_i * arg) 
    
    ;; create a phase correction map
    phase = complex(rebin(real_part(phase_x), xdim,ydim, /sample), $
                    rebin(imaginary(phase_x), xdim,ydim, /sample))
                    
    ;; apply correction
    slice *= phase

    ;; scale phase ramp 0..1 in image space
    index = transpose((findgen(ydim) - pivot_y)/ydim)
    
    ;; phase = X1 * ramp + X0
    arg = ( float(phase_y1)*index )
    phase_y = exp(two_pi_i * arg) 

    ;; create a phase correction map
    phase = complex(rebin(real_part(phase_y), xdim,ydim, /sample), $
                    rebin(imaginary(phase_y), xdim,ydim, /sample))

    ;; apply correction
    slice *= phase
    
end

;+
; :Description:
;    Handles the "Apply" button click. Applies correction values to the data set.
;
; :Params:
;    event
;
;
;
; :Author: wtriplett
;-
pro mas_image_phase_apply_event, event

    common scan_data, project
    
    widget_control, event.top, get_uvalue=state
    
    widget_control, (*state).sl_aid, get_value   = a_num
    widget_control, (*state).sl_slice, get_value = s_num

    pdata = (*state).pdata

    xdim = (*state).data_dim[0]
    ydim = (*state).data_dim[1]
    sdim = (*state).data_dim[2]
    adim = (*state).data_dim[3]
    
    for s_num = 0, sdim-1 do begin
    
        for a_num = 0, adim-1 do begin
        
            slice = (*pdata)[*,*,s_num,a_num]
            mas_image_phase_apply_correction, slice, (*state).phase_corr_params[s_num,a_num]
            (*pdata)[*,*,s_num,a_num] = slice
            
        endfor
    endfor
    
    ptr_free, state
    
    widget_control, event.top, /destroy

    project.procpramarray[project.ci].state_2 = 0
    
end

;+
; :Description:
;    Handles the Autophase Button click.
;
; :Params:
;    event
;
;
;
; :Author: wtriplett
;-
pro mas_image_phase_auto_event, event

    widget_control, event.top, get_uvalue=state

    YesNo = dialog_message(['Warning: Autophase is experimental.', $
                            'Autophase works on the data in memory, not a copy.', $
                            'If you are unhappy with the results, please reload the scan from disk.', $
                            'Do you want to continue?' ], $
                            /center, /question)
    if (YesNo eq 'No') then return
    
    widget_control, (*state).sl_aid,   get_value = a_num
    widget_control, (*state).sl_slice, get_value = s_num
    
    mas_image_phase_autophase, state
    
    mas_image_phase_slice_change_event, { top: event.top, value: 0 }
    mas_image_phase_refresh_display, state

end

;+
; :Description:
;    Handles the Cancel Button click.
;
; :Params:
;    event
;
;
;
; :Author: wtriplett
;-
pro mas_image_phase_cancel_event, event

    widget_control, event.top, get_uvalue=state
    ptr_free, state
    widget_control, event.top, /destroy
    
end


;+
; :Description:
;    Handes events generated by slice-change slider.
;
; :Params:
;    event
;
;
;
; :Author: wtriplett
;-
pro mas_image_phase_slice_change_event, event

    widget_control, event.top, get_uvalue=state
    
    widget_control, (*state).sl_aid, get_value   = a_num
    widget_control, (*state).sl_slice, get_value = s_num

    widget_control, (*state).sl_x0, set_value=(*state).phase_corr_params[s_num,a_num].phase_0*360.0
    widget_control, (*state).sl_x1, set_value=(*state).phase_corr_params[s_num,a_num].phase_x1*1e3
    widget_control, (*state).sl_xp, set_value=(*state).phase_corr_params[s_num,a_num].pivot_x
    
    widget_control, (*state).sl_y1, set_value=(*state).phase_corr_params[s_num,a_num].phase_y1*1e3
    widget_control, (*state).sl_yp, set_value=(*state).phase_corr_params[s_num,a_num].pivot_y

    mas_image_phase_event, event
    
end

;+
; :Description:
;    Event handler to handle click/drag in the main drawing window. Events
;    are used to rotate the display in 3D
;
; :Params:
;    event
;
;
;
; :Author: wtriplett
;-
pro mas_image_phase_draw_event, event

    ;help, event, /str
    widget_control, event.id, get_uvalue=uv
    widget_control, event.top, get_uvalue=state
    
    if (event.release) then begin
        uv.down = 0
    endif else if (event.press) then begin
        uv.down = 1
        uv.down_x = event.x
        uv.down_y = event.y
    endif
    
    if (uv.down) then begin
        dx = (event.x - uv.down_x)/2
        dy = (event.y - uv.down_y)/2
        
        vy = (*state).rotation_ax
        vx = (*state).rotation_az
        
        nz = ((vx+dx)/1.0) mod 360
        if (nz le 0) then nz = 360
        
        nx = ((vy-dy)/1.0) mod 360
        if (nx le 0) then nx = 360
        
        (*state).rotation_az = nz
        (*state).rotation_ax = nx
                
        mas_image_phase_refresh_display, state
        
        uv.down_x = event.x
        uv.down_y = event.y

    endif
        
    widget_control, event.id, set_uvalue=uv
    
end

;+
; :Description:
;    Handes event generated by main GUI.
;
; :Params:
;    event
;
;
;
; :Author: wtriplett
;-
pro mas_image_phase_event, event

    widget_control, event.top, get_uvalue=state
    
    widget_control, (*state).sl_aid, get_value   = a_num
    widget_control, (*state).sl_slice, get_value = s_num
    
    widget_control, (*state).sl_x0, get_value=sl_x0
    widget_control, (*state).sl_x1, get_value=sl_x1
    widget_control, (*state).sl_xp, get_value=sl_xp
    
    widget_control, (*state).sl_y1, get_value=sl_y1
    widget_control, (*state).sl_yp, get_value=sl_yp

    xdim = (*state).data_dim[0]
    ydim = (*state).data_dim[1]
    sdim = (*state).data_dim[2]
    adim = (*state).data_dim[3]

    ;; divide the 0, X1, Y1 by max slider value to scale between [0,1]
    if (widget_info((*state).btn_apply_across_adim, /button_set)) then begin
        (*state).phase_corr_params[s_num,*].phase_0  = sl_x0/360.0
        (*state).phase_corr_params[s_num,*].pivot_x  = sl_xp
        (*state).phase_corr_params[s_num,*].phase_x1 = sl_x1/1e3
        (*state).phase_corr_params[s_num,*].pivot_y  = sl_yp
        (*state).phase_corr_params[s_num,*].phase_y1 = sl_y1/1e3
    endif else begin
        (*state).phase_corr_params[s_num,a_num].phase_0  = sl_x0/360.0
        (*state).phase_corr_params[s_num,a_num].pivot_x  = sl_xp
        (*state).phase_corr_params[s_num,a_num].phase_x1 = sl_x1/1e3
        (*state).phase_corr_params[s_num,a_num].pivot_y  = sl_yp
        (*state).phase_corr_params[s_num,a_num].phase_y1 = sl_y1/1e3
    endelse
    
    pdata = (*state).pdata
    slice = (*pdata)[*,*,s_num,a_num]
    ;; apply
    mas_image_phase_apply_correction, slice, (*state).phase_corr_params[s_num,a_num]

    mas_image_phase_refresh_display, state, slice
    
end

;+
; :Description:
;    Procedure to correct phase of a complex 2D image.
;
; :Params:
;    pdata - pointer to complex image data in 2D, 3D or 4D format
;
;
; :Author: wtriplett
;-
pro mas_image_phase, pdata
    
    common scan_data, project
    
    set_plot, !VERSION.OS_FAMILY eq 'unix' ? 'X' : 'WIN'
    device, get_screen_size=screen_size
    
    if (not ptr_valid(pdata)) then begin
        mas_load_state_1
        pdata = project.dataarray[project.ci].state1
    endif
    
    sz_data = size(*pdata, /dimensions)

    xdim = sz_data[0]
    ydim = sz_data[1]
    sdim = (n_elements(sz_data) gt 2) ? sz_data[2] : 1
    adim = (n_elements(sz_data) gt 3) ? sz_data[3] : 1
    
    order1_bound = 5000
    
    base = widget_base(title='Image Phase Adjustment', /column)
    
    imbase   = widget_base(base)
    drawmain = widget_draw(imbase, scr_xsize=screen_size[0]*0.8, $
                                   scr_ysize=screen_size[1]*0.6, $
                                   graphics_level=1, /button_events, /motion_events, $
                                   event_pro='mas_image_phase_draw_event')
    widget_control, drawmain, set_uvalue={ down: 0, down_x : 0, down_y : 0 }
    controlbase = widget_base(base, /row)
    
    ;; Holds two slider bases
    slider_parent = widget_base(controlbase, /column)
    
    ;; Holds the sliders that control phase correction
    slbase = widget_base(slider_parent, /row, /frame)
    sl_x0 = cw_fslider(slbase, min=-360, max=360, value=0, title='0-phase', /drag, /edit, scroll=1.0)
    
    sl_xp = widget_slider(slbase, min=0, max=xdim, value=xdim/2, title='X-pivot', /drag, scroll=1)
    sl_x1 = widget_slider(slbase, min=-order1_bound, max=order1_bound, value=0, title='X_1-phase', /drag, scroll=1)
    
    sl_yp = widget_slider(slbase, min=0, max=ydim, value=ydim/2, title='Y-pivot', /drag, scroll=1)    
    sl_y1 = widget_slider(slbase, min=-order1_bound, max=order1_bound, value=0, title='Y_1-phase', /drag, scroll=1)
    
    ;; Holds slice/array select, and display axis rotation sliders
    slbase = widget_base(slider_parent, /row, /frame, /base_align_center)
    
    btn_base = widget_base(slbase, /column, /nonexclusive)
    btn_apply_across_adim = widget_button(btn_base, value='Same correction per slice for all array elements')
    widget_control, btn_apply_across_adim, /set_button
    
    if (adim gt 1) then begin
        sl_aid = widget_slider(slbase, min=0, max=adim-1, value=0, title='Array', /drag, event_pro='mas_image_phase_slice_change_event')
    endif else begin
        sl_aid = widget_slider(slbase, min=0, max=100, value=0, title='Array', sensitive=0)
    endelse
    if (sdim gt 1) then begin
        sl_slice = widget_slider(slbase, min=0, max=sdim-1, value=0, title='Slice', /drag, event_pro='mas_image_phase_slice_change_event')
    endif else begin
        sl_slice = widget_slider(slbase, min=0, max=100, value=0, title='Slice', sensitive=0)
    endelse
    
    ;; apply/cancel
    btn_base   = widget_base(controlbase, /align_top, /base_align_top, row=3, /frame, /grid)
    btn_apply  = widget_button(btn_base, value='Apply', event_pro='mas_image_phase_apply_event')
    btn_auto   = widget_button(btn_base, value='Auto Phase', event_pro='mas_image_phase_auto_event')
    btn_cancel = widget_button(btn_base, value='Cancel', event_pro='mas_image_phase_cancel_event')
    
    ;; display windows for graphics
    ;max_xdim = 
    dbase   = widget_base(controlbase, /align_right, /base_align_right, /row)
    drawimg = widget_draw(dbase, scr_xsize=175*2, scr_ysize=175, graphics_level=1, /align_right)
    ;;plotimg = widget_draw(dbase, scr_xsize=300, scr_ysize=ydim, graphics_level=1)
    
    ;; open window and grab the display window IDs
    widget_control, base, /realize
    widget_control, drawmain, get_value=wid
    widget_control, drawimg, get_value=imgwid
    ;;widget_control, plotimg, get_value=plotwid
    
    ;; Set up some phase correction parameters for each slice
    ;; Right now, one phase correction is appled per slice, across all array elements.
    params = replicate({ phase_0: float(0), $
                         pivot_x: xdim/2, $
                         phase_x1: float(0), $
                         pivot_y: ydim/2, $ 
                         phase_y1: float(0) }, sdim*adim)
    params = reform(params, sdim, adim)
    
    ;; Save some widget stuff
    state = { data_dim: [xdim, ydim, sdim, adim], $
              data_range: max(abs(*pdata))*[-1.01, 1.01], $
              sl_x0: sl_x0, $
              sl_x1: sl_x1, $
              sl_xp: sl_xp, $
              sl_y1: sl_y1, $
              sl_yp: sl_yp, $
              sl_aid: sl_aid, $
              sl_slice: sl_slice, $
              btn_apply_across_adim: btn_apply_across_adim, $
              surf_window: wid, $
              img_window: imgwid, $
              $;;plot_window: plotwid, $
              rotation_az: 45, $
              rotation_ax: 45, $
              btn_apply: btn_apply, $
              pdata: pdata, $
              imag_history: fltarr(500), $
              imag_hist_pt: 0, $
              phase_corr_params: params }
     
    widget_control, base, set_uvalue=ptr_new(state)
    
    ;; This just causes the screen to be refreshed one time 
    mas_image_phase_event, { top: base }
    
    xmanager, 'mas_image_phase', base
    
end
