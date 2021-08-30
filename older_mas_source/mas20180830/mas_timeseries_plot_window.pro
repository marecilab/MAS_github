
function mas_timeseries_plot_init, proj_id, group_leader=group_leader

    common scan_data, project
    
    if (n_elements(proj_id) eq 0) then return, 0

    ci = proj_id
    
    ;; get parameters from MAS
    TR = project.imndarray[ci].recov_time*1e-3
    if (TR eq 0) then begin
        void = dialog_message(['Warning: TR parameter for this scan is 0!', $
                               'Setting it to 1, but please set to proper value in GUI.'], $
                              /error, /center)
        TR = 1.0
    endif
    
    EV_time = 60
    start_index = project.imndarray[ci].epi_num_ref_scans
    NA = project.imndarray[ci].n_avg
    
    if (n_elements(group_leader) eq 0) then begin
        base = widget_base(title='mas_timeseries_plot', row=2)
    endif else begin
        base = widget_base(title='mas_timeseries_plot', group_leader=group_leader, row=2)
    endelse
    
    ;; the rest of this just sets up the GUI for the plot window.
    base_row = widget_base(base, /row)
    
    b = widget_base(base_row, /column, /frame)
    lbl = widget_label(b, value='Delta t (s):')
    txt_delta_t = widget_text(b, xsize=4, value=strtrim(TR*NA,2), /editable)

    b = widget_base(base_row, /column, /frame)
    lbl = widget_label(b, value='Event time (s):')
    txt_key_interval = widget_text(b, xsize=4, value=strtrim(EV_time,2), /editable)

    b = widget_base(base_row, /column, /frame)
    lbl = widget_label(b, value='Start Index:')
    txt_start_index = widget_text(b, xsize=4, value=strtrim(start_index,2), /editable)
    
    b = widget_base(base_row, /column, /frame)
    lbl = widget_label(b, value='Time Series Correction')
    dl_ts_correction = widget_droplist(b, value=['Uncorrected', 'De-mean', 'Linear Drift Correction'])

;    b = widget_base(base_row, /column, /frame)
;    lbl = widget_label(b, value='TS Smoothing (avg.)')
;    txt_smoothing = widget_text(b, xsize=4, value='0', /editable)

    b = widget_base(base_row, /column, /frame)
    lbl = widget_label(b, value='Spectrum Type')
    dl_fq_type = widget_droplist(b, value=['Raw', 'Coherence'])
    
    b = widget_base(base_row, /column, /frame)
    lbl = widget_label(b, value='Create Coherence Map')
    b1 = widget_base(b, /row)
    btn_coh_map_all = widget_button(b1, value='All Freqs.')
    btn_coh_map_act = widget_button(b1, value='Activation Freq.')

    b = widget_base(base_row, /column, /frame, /align_center)
    lbl = widget_label(b, value='Colors')
    b1 = widget_base(b, /row)
    btn_invert_colors = widget_button(b1, value='Invert')
    
    device, get_screen_size=screen_size
    maxdim = (screen_size[0]*0.75) - 50
    b = widget_draw(base, xsize=875, ysize=375)
    
    widget_control, base, /realize
    widget_control, b, get_value=window_id
    tool_state = { tlb: base, $
                   proj_id: 0L, $
                   current_data: ptr_new(), $
                   current_coord: lonarr(3), $
                   color_scheme: 0B, $
                   is_active: 1B, $
                   window_id: window_id, $
                   txt_delta_t: txt_delta_t, $
                   txt_key_interval: txt_key_interval, $
                   txt_start_index: txt_start_index, $
                   btn_coh_map_all: btn_coh_map_all, $
                   btn_coh_map_act: btn_coh_map_act, $
                   btn_invert_colors: btn_invert_colors, $
                   $;txt_smoothing: txt_smoothing, $
                   dl_ts_correction: dl_ts_correction, $
                   dl_fq_type: dl_fq_type }
    widget_control, base, set_uvalue=ptr_new(tool_state)
    
    xmanager, 'mas_timeseries_plot', base, /no_block, cleanup='mas_timeseries_plot_cleanup'
    return, base

end

pro mas_timeseries_plot_cleanup, tlb
    
    widget_control, tlb, get_uvalue=tool_state
    
    ptr_free, (*tool_state).current_data
    ptr_free, tool_state


end

pro mas_timeseries_plot_event, event
    
    widget_control, event.top, get_uvalue=tool_state
    
    ;; Specific event handlers exist for the two button events. The rest of the events
    ;; are handled by mas_timeseries_plot_update. The values from the GUI are read 
    ;; inside that procedure rather than here (for now).
    case event.id of
    
        (*tool_state).btn_coh_map_all: begin

            widget_control, (*tool_state).txt_delta_t, get_value=delta_t
            delta_t = float(delta_t[0])
    
            widget_control, (*tool_state).txt_key_interval, get_value=key_interval
            key_interval = float(key_interval[0])

            widget_control, (*tool_state).txt_start_index, get_value=start_index
            start_index = long(start_index[0])
        
            mas_timeseries_coherence_map, delta_t, start_index, key_interval, proj_id=(*tool_state).proj_id, /all_freq
            
        end
        
        (*tool_state).btn_coh_map_act: begin
        
            widget_control, (*tool_state).txt_delta_t, get_value=delta_t
            delta_t = float(delta_t[0])
    
            widget_control, (*tool_state).txt_key_interval, get_value=key_interval
            key_interval = float(key_interval[0])

            widget_control, (*tool_state).txt_start_index, get_value=start_index
            start_index = long(start_index[0])
        
            mas_timeseries_coherence_map, delta_t, start_index, key_interval, proj_id=(*tool_state).proj_id
        
        end
        
        (*tool_state).btn_invert_colors: begin
            (*tool_state).color_scheme = 1 - (*tool_state).color_scheme
            mas_timeseries_plot_update, tool_state
        end
        
        else: mas_timeseries_plot_update, tool_state

    endcase

end

pro mas_timeseries_plot, coord, state

    common scan_data
    
    if (project.imndarray[(*state).proj_index].adim le 1) then return
    
    ;; if no plot window exists, create one
    if (widget_info((*state).time_series_gui_win, /valid_id) eq 0) then begin
        (*state).time_series_gui_win = mas_timeseries_plot_init((*state).proj_index, group_leader=(*state).base_main)
    endif
    
    ;; get the scan id from the ortho viewer and set it in the plot window.
    widget_control, (*state).time_series_gui_win, get_uvalue=tool_state
    (*tool_state).proj_id = (*state).proj_index
    if ((*tool_state).is_active eq 0) then return

    ;; get the data values from MAS's data array for the correct scan
    values = reform((*project.dataarray[(*state).proj_index].state1)[coord[0], $
                                                                     coord[1], $
                                                                     coord[2], *])
    (*tool_state).current_coord = coord
    if (ptr_valid((*tool_state).current_data)) then begin
        *((*tool_state).current_data) = values
    endif else begin
        (*tool_state).current_data = ptr_new(values)
    endelse
    
    ;; call update to refresh the window.
    mas_timeseries_plot_update, tool_state

end

pro mas_timeseries_plot_update, tool_state
    
    common scan_data, project
    
    ;; Get plot parameters from widgets
    widget_control, (*tool_state).txt_delta_t, get_value=delta_t
    delta_t = float(delta_t[0])
    
    widget_control, (*tool_state).txt_key_interval, get_value=key_interval
    key_interval = float(key_interval[0])

    widget_control, (*tool_state).txt_start_index, get_value=start_index
    start_index = long(start_index[0])

    ;; smoothing is not enabled at the moment.
    ;widget_control, (*tool_state).txt_smoothing, get_value=smooth_width
    smooth_width = 0;long(smooth_width[0])
    
    ts_corr_type    = widget_info((*tool_state).dl_ts_correction, /droplist_select)
    fq_display_type =  widget_info((*tool_state).dl_fq_type, /droplist_select)
    
    subtract_mean = 0
    deline        = 0
    case ts_corr_type of 
        0:
        1: subtract_mean  = 1
        2: deline = 1
        else:
    endcase 
    
    TR_sec = delta_t
    
    values = (*((*tool_state).current_data))[start_index:*]

    ;; setup for split screen plot
    !P.multi = [0,1,2]
    nvalues = n_elements(values)
    
    ;; perform demean or linear drift correction
    if (subtract_mean) then begin
        values -= mean(values)
        yfit = 0.0
    endif else if (deline) then begin
        values -= mean(values)
        lin = linfit(findgen(nvalues)*TR_sec, values, yfit=yfit)
    endif else begin
        yfit = 0.0
    endelse
    
    color_scheme = (*tool_state).color_scheme
    
    ;; create axes for the plots.
    ts_axis = findgen(nvalues)*TR_sec
    fq_axis = (findgen(nvalues/2-1)/(nvalues/2-1)) * 1.0/(2*TR_sec)
    set_plot, (!VERSION.OS_FAMILY eq 'unix') ? 'X' : 'WIN'
    wset, (*tool_state).window_id
    device, decomposed=1
    
    ;; compute spectrum and find the activation freq.
    freq = fft(values-yfit)
    freq = freq[0:nvalues/2-1]
    freq = abs(freq)
    interval = key_interval
    act_freq = 1.0/interval
    xthk = 1
    ythk = 1
    
    ;; Set up color scheme
    if (color_scheme eq 0) then begin
        ;; white on black
        plot_bgcolor = '000000'x
        plot_fgcolor = 'FFFFFF'x
        marker_color = '0000FF'x
    endif else begin
        plot_bgcolor = 'FFFFFF'x
        plot_fgcolor = '000000'x
        marker_color = 'FF0000'x    
    endelse
    
    min_values = min(values, max=max_values)
    
    ;; plot the time series
    plot,  ts_axis, values-yfit, $
           background=plot_bgcolor, color=plot_fgcolor, $
           xstyle=1, ystyle=1, $
           yticklen=-0.01, xthick=xthk, ythick=ythk, $
           yrange=[min_values,max_values], $
           title=string((*tool_state).current_coord, format='(%"Voxel Time Series: [%d,%d,%d]")'), $
           xtitle='Time (seconds)', ytitle='Signal Intensity'
    ;; overplot the data point symbols
    oplot, ts_axis, values-yfit, psym=1, color=plot_fgcolor
    
    ;; overplot the activation / stimulation invervals
    for i = 1, ceil(max(ts_axis)/interval) do begin
        oplot, interval*[i,i], [min_values,max_values], linestyle=1, color=marker_color, thick=2
    endfor
    
    ;; prepare to compute the spectrum
    min_fq_axis = min(fq_axis, max=max_fq_axis)
    void = min(abs(fq_axis-act_freq), min_loc)
    
    ;; create spectrum plot title and range
    if (fq_display_type eq 1) then begin
        freq /= sqrt(total(freq^2))
        fq_title = string((*tool_state).current_coord, format='(%"Coherence value: [%d,%d,%d]")')
        fq_yrange = [0, 1]
    endif else begin
        fq_title = string((*tool_state).current_coord, format='(%"Raw spectrum: [%d,%d,%d]")')
        fq_yrange = [0, 0.5]
    endelse
    fq_title = string(act_freq, format='(%"Freq. @ marker: %0.3f")')+'           '+$
                fq_title+'       '+$
                string(freq[min_loc], min_loc, format='(%"Value @ marker: %0.3f (index=%d)")')
    
    ;; plot the spectrum
    plot, fq_axis, freq, $
          background=plot_bgcolor, color=plot_fgcolor, $
          xstyle=1, ystyle=1, yticklen=-0.01, xticklen=-0.025, xthick=xthk, ythick=ythk, $
          yrange=fq_yrange, title=fq_title, xtitle='Frequency (Hz)', xticks=15
    ;; overplot the activation marker
    oplot, act_freq*[1,1], [0,1000], linestyle=1, color=marker_color, thick=2
    
    ;;compute the position of the marker text
    ;; disabled in favor of putting marker info in plot title
;    if (min_loc lt (n_elements(fq_axis)/2)) then begin
;        marker_xpos = act_freq+(max_fq_axis-min_fq_axis)*0.01
;        marker_align = 0.0
;    endif else begin
;        marker_xpos = act_freq-(max_fq_axis-min_fq_axis)*0.01
;        marker_align = 1.0
;    endelse
;    
;    ;; show the marker text
;    xyouts, marker_xpos, fq_yrange[1]*0.9, $
;            string(freq[min_loc], min_loc, format='(%"Value @ marker: %0.3f (index=%d)")'), $
;            alignment=marker_align
            
end

pro mas_timeseries_coherence_map, eff_tr, start_idx, period, proj_id=proj_id, all_freq=all_freq, coherence=coherence

    common scan_data, project

    if (n_elements(start_idx) eq 0) then start_idx = 0
    
    if (n_elements(proj_id) eq 0) then proj_id = project.ci
    
    subset   = (*project.dataarray[proj_id].state1)[*,*,*,start_idx:*]
    dims     = size(subset, /dimensions)
    TR_sec   = float(eff_tr)
    nvalues  = float(dims[3])
    interval = float(period)
    
    fq_axis = (findgen(nvalues/2-1)/(nvalues/2-1)) * 1.0/(2*TR_sec)
    void = min(abs(fq_axis-1.0/interval), min_loc)

    pbar = obj_new('progressbar', title='Computing Coherence...', text='Computing Coherence...')
    pbar->Start

    for x = 0, dims[0]-1 do begin
        for y = 0, dims[1]-1 do begin
            for z = 0, dims[2]-1 do begin
                mn = mean(reform(subset[x,y,z,*]))
                subset[x,y,z,*] -= mn
            endfor
        endfor
        
    endfor
    
    pbar->Update, 33.0
    
    fsubset = fft(subset, dimension=4)
    fsubset = abs(fsubset[*,*,*,0:dims[3]/2-1])        

    block = (sqrt(total(fsubset^2, 4)))

    pbar->Update, 66.0
    
    for x = 0, dims[0]-1 do begin
        for y = 0, dims[1]-1 do begin
            for z = 0, dims[2]-1 do begin
                temp = reform(fsubset[x,y,z,*])
                temp = temp / replicate(reform(block[x,y,z]), dims[3]/2)
                fsubset[x,y,z,*] = temp
            endfor
        endfor
    endfor

    pbar->update, 100.0
    pbar->Destroy
    
    if (keyword_set(all_freq)) then begin
        coherence = temporary(fsubset)
    endif else begin
        coherence = reform(fsubset[*,*,*,min_loc])
        print, "Extracting index: "+strtrim(string(min_loc),2)
    endelse
   
    
    nifti_file = dialog_pickfile(title='Enter a filename to save the coherence:', filter='*.nii', $
                                default_extension='nii')
    if (nifti_file eq '') then return
   
    p = ptr_new(coherence)
    ci = project.ci
    project.ci = proj_id
    mas_export_nifti, data_ptr=p, file_name=nifti_file
    project.ci = ci
    ;;mas_display_ortho, data_ptr=p
    ptr_free, p
   
end

