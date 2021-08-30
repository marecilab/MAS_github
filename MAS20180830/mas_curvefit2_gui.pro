;+
; :Description:
;    Rotates and image using MAS's mas_rotate_flip procedure.
;    This applies the user's settings for rotate and flip found
;    in the MAS main window to the displayed image
;
; :Params:
;    image_data - the image data to rotate and flip
;
; :Author: wtriplett
;-
function mas_cfit2_rotate_flip, image_data

    ptmp = ptr_new(image_data)
    mas_rotate_flip, ptmp
    image_out = *ptmp
    ptr_free, ptmp
    
    return, image_out

end

;
;+
; :Description:
;    Creates the default state pointer, which is saved
;    in the mas common block (scan_data). 
;
;
;
;
;
; :Author: wtriplett
;-
function mas_cfit2_make_default_state

    common scan_data, project
    ci = project.ci
    
    ;; models of fitting, all 3 arrays need to be updated when a model is added or removed!
    model_names = [ 'Please Select', $
                    'T1: Inversion Recovery SE', $
                    'T1: Inversion Recovery SE (ABS w/ Polarity Correction)', $
                    'T1: Inversion Recovery SE (ABS w/ Powell Optimizer)', $
                    'T1: Saturation Recovery SE', $
                    'T2: Variable TE w/o Baseline', $
                    'T2: Variable TE w/ Baseline', $
                    'ADC: Multi B-value' ]
    
    model_pros_roi = [ '', $
                       'mas_cfit2_ROI_FIT_T1_VarTI', $
                       'mas_cfit2_ROI_FIT_T1_VarTI_ABS', $
                       'mas_cfit2_ROI_FIT_T1_VarTI_POWELL', $
                       'mas_cfit2_ROI_FIT_T1_VarTR', $
                       'mas_cfit2_ROI_FIT_T2_VarTE', $
                       'mas_cfit2_ROI_FIT_T2_wBaseline_VarTE', $
                       'mas_cfit2_ROI_FIT_ADC_VarB' ]
                       
    model_pros_img = [ '', $
                       'mas_cfit2_IMG_FIT_T1_VarTI', $
                       'mas_cfit2_IMG_FIT_T1_VarTI_ABS', $
                       'mas_cfit2_IMG_FIT_T1_VarTI_POWELL', $
                       'mas_cfit2_IMG_FIT_T1_VarTR', $
                       'mas_cfit2_IMG_FIT_T2_VarTE', $
                       'mas_cfit2_IMG_FIT_T2_wBaseline_VarTE', $
                       'mas_cfit2_IMG_FIT_ADC_VarB' ]

    state = { gui: ptr_new(), $                  ;; pointer to GUI elements
              project_ci: ci, $                  ;; index of the current scan
              project_id: randomu(systime(1)), $ ;; unique identifier for this scan
              model_names: model_names, $        ;; display names of the models that can be fit
              model_pros_roi: model_pros_roi, $  ;; model procedures for ROI fitting
              model_pros_img: model_pros_img, $  ;; model procedures for image fitting
              plot_layout: 'MULTIPLOT', $        ;; controls how ROI plotting is handled
              show_histogram: 0, $               ;; 0,1 whether or not to display the histogram
              combined_yscale: 0, $              ;; 0,1 whether yscale is kept fixed between region plots
              roi_show_stats_window: 1, $        ;; 0,1 whether or not to show stats window
              roi_save_ps: 0, $                  ;; 0,1 whether to save plots to ps or display on screen
              roi_ps_output_file: '', $          ;; filename and path to .ps destination file
              roi_stats_base: 0L, $              ;; widget id of ROI stats window
              roi_append_stats: 1B, $            ;; 0,1 whether or not to append to ROI stats or replace
              fit_model: 0, $                    ;; index into the above three arrays of which model we are fitting
              meas_type: 0, $                    ;; 0,1,2... control how the signal is created from ROI
              img_histogram_bins: 50, $          ;; # of bins in histogram when displaying 
              img_fit_threshold: 0.10, $         ;; the percentage of max signal cutoff for voxel fitting
              img_fit_result: ptr_new(), $       ;; pointer to the fit result from an image fit
              img_stats_base: 0L, $              ;; widget id of img stats window
              img_save_ps: 0, $                  ;; 0,1 whether or not to save PS of image fit plots
              img_ps_output_file: '', $          ;; filename for above
              vol_fit_result: ptr_new(), $       ;; pointer to hold volume fit results
              vol_fit_threshold: 0.10, $         ;; the percentage of max signal cutoff for volume fitting
              current_path: project.current_path $ ;; current project path 
             }
 
    return, ptr_new(state, /no_copy)
end

pro mas_cfit2_vb_display_vol_event, event

    widget_control, event.id, get_uvalue=param_state
    if (n_elements(param_state) eq 0) then begin
        return
    endif

    widget_control, event.top, get_uvalue=state
    fr = (*state).vol_fit_result

    param_name = widget_info(event.id, /uname)
    parinfo = (*fr).parinfo
    parnum = 0L
    
    widget_control, param_state.min_wid, get_value=min_scl
    widget_control, param_state.max_wid, get_value=max_scl
    
    for i = 0, n_elements(parinfo)-1 do begin
        if (parinfo[i].fixed) then continue
        
        if (param_name eq parinfo[i].parname) then begin
        
            temp = ptr_new(reform((*fr).param_maps[*,*,*,parnum]))
            
            lo = where(*temp lt min_scl, n_lo)
            hi = where(*temp gt max_scl, n_hi)
            if (n_lo gt 0) then (*temp)[lo] = min_scl
            if (n_hi gt 0) then (*temp)[hi] = max_scl
            
            mas_display_ortho, data_ptr=temp;, /block
            ;ptr_free, temp
            return
            
        endif
        
        parnum++
        
    endfor
    
    if (param_name eq 'CHI2') then begin
        temp = ptr_new(reform((*fr).chisq_map))
        lo = where(*temp lt min_scl, n_lo)
        hi = where(*temp gt max_scl, n_hi)
        if (n_lo gt 0) then (*temp)[lo] = min_scl
        if (n_hi gt 0) then (*temp)[hi] = max_scl
        
        mas_display_ortho, data_ptr=temp, /block
        ptr_free, temp
    endif
end

pro mas_cfit2_vb_export_vol_event, event

    
    widget_control, event.id, get_uvalue=param_state
    if (n_elements(param_state) eq 0) then begin
        return
    endif

    widget_control, event.top, get_uvalue=state
    fr = (*state).vol_fit_result

    param_name = widget_info(event.id, /uname)
    parinfo = (*fr).parinfo
    parnum = 0L

    file_name = dialog_pickfile(title='Please select a file name for the '+param_name+'map.', $
                                default_extension='nii', filter='*.nii;*.nii.gz', /overwrite_prompt)
    if (file_name eq '') then return
    
    parnum = 0L
    for i = 0, n_elements(parinfo)-1 do begin
        if (parinfo[i].fixed) then continue
        
        if (param_name eq parinfo[i].parname) then begin
        
            temp = ptr_new(reform((*fr).param_maps[*,*,*,parnum]))
                        
            mas_export_nifti, file_name=file_name, data_ptr=temp
            
            ptr_free, temp
            
            return
        endif
        
        parnum++
        
    endfor
    
    if (param_name eq 'CHI2') then begin
            temp = ptr_new(reform((*fr).chisq_map))
            mas_export_nifti, file_name=file_name, data_ptr=temp
            ptr_free, temp
    endif
    
end

;+
; :Description:
;    Updates the variable base which holds the widgets that control
;    display and stats related to image fitting.
;
; :Params:
;    event
;
; :Author: wtriplett
;-
pro mas_cfit2_vb_stats_img_event, event

    widget_control, event.top, get_uvalue=state
    fr = (*state).img_fit_result
    gui = *((*state).gui)
    
    widget_control, event.id, get_uvalue=param_state
    if (n_elements(param_state) eq 0) then begin
        return
    endif

    widget_control, param_state.min_wid, get_value=min_scl
    widget_control, param_state.max_wid, get_value=max_scl

    hist_bins = (*state).img_histogram_bins
    
    param_name = widget_info(event.id, /uname)

    parinfo = (*fr).parinfo
    
    if (param_name eq 'CHI2') then begin
        temp = reform((*fr).chisq_map)
        parunits = ''
        parname = 'Chisquared'
    endif else begin    
        parnum = 0L
        for i = 0, n_elements(parinfo)-1 do begin
    
            if (parinfo[i].fixed) then continue
            
            if (param_name eq parinfo[i].parname) then begin                
                temp = reform((*fr).param_maps[*,*,parnum])
                parname = param_name
                parunits = parinfo[i].parunits
                break
            endif
    
            parnum++
            
        endfor
    endelse
    
    if (n_elements(temp) eq 0) then return

    region_names = mas_roi_get_current_name()
    num_regions = n_elements(region_names)
    if (num_regions eq 1 and region_names[0] eq '') then begin
        void = dialog_message(['Selected ROI set has no regions.', $
                               'Please create an ROI to use this feature.'], $
                               /error, /center)
        return
    endif

    save_ps = (*state).img_save_ps
    psfile  = (*state).img_ps_output_file
    psfile  = file_dirname(psfile, /mark_directory)+file_basename(psfile, '.ps')+'_'+param_name+'.ps'
    !P.multi = [0, 1, num_regions]
    if (save_ps ne 0 and n_elements(psfile) ne '') then begin
        PS=1
        set_plot, 'PS'
        filename=psfile
        device, filename=filename, bits_per_pixel=8, font_size=11, language_level=2
        device,  /landscape, /color, /inches, /TIMES, /BOLD, $
                xsize=11, ysize=8.5, xoffset=0.25, decomposed=1, $
                yoffset=11-(0.25), scale_factor=0.95
    endif else begin
        PS=0
        set_plot, (!VERSION.OS_FAMILY eq 'unix' ) ? 'X' : 'WIN'
        device, decomposed=0
        device, get_screen_size=scrn_size
        win_xsize = scrn_size[0]*0.75
        win_ysize = scrn_size[1]*0.75
        window, 8, xsize=win_xsize, ysize=win_ysize
    endelse
    color_bg = 255L
    color_fg = 0L
    
    !P.multi = [0, 1, num_regions]
    
    indep_scale = widget_info(widget_info(event.top, find_by_uname=param_name+'_scale'), /button_set)
    
    stats = replicate({ region_name: '', mean: 0.0, stdev: 0.0, min: 0.0, max: 0.0, count: 0L }, num_regions)
    
    for r = 0, num_regions-1 do begin
    
        region_mask = mas_roi_get_current_mask(region_num=r,/no_transform)
        active_ind  = where(*region_mask ne 0, n_active)
        
        if (n_active eq 0) then continue
        
        data = temp[active_ind]
        if (indep_scale) then begin
            hist_min = min(data, max=hist_max)
        endif else begin
            hist_min = min_scl
            hist_max = max_scl
        endelse
        
        stats[r].region_name = region_names[r]
        stats[r].mean = mean(data)
        stats[r].stdev = stdev(data)
        stats[r].min = min(data)
        stats[r].max = max(data)
        stats[r].count = n_elements(data)

        hist = histogram(data, nbins=hist_bins, min=hist_min, max=hist_max, locations=loc)
        min_hist = min(hist, max=max_hist)
        min_loc  = min(loc, max=max_loc)
        !X.charsize = 1.5
        !Y.charsize = 1.5
        !P.charsize = 1.25
        plot, loc, replicate(0, n_elements(loc)), /nodata, color=color_fg, background=color_bg, $
              title='Histogram for '+parname+' (Region: '+region_names[r]+')', $
              xtitle=parname+' '+parunits, $
              yrange=[0, max_hist], ticklen=-!P.ticklen, $
              ytitle='Count', xstyle=1+8, ystyle=1+8, $
              xmargin=[17,35] ,ymargin=[5,5]
              
        oplot, loc, hist, color=color_fg
        oplot, loc, hist, psym=1, color=color_fg

        xyouts, max_loc, max_hist*0.95, string(stats[r].mean, format='(%"Mean: %0.3e")'), /data, /noclip, color=color_fg
        xyouts, max_loc, max_hist*0.85, string(stats[r].stdev, format='(%"Stdev: %0.3e")'), /data, /noclip, color=color_fg
        xyouts, max_loc, max_hist*0.75, string(stats[r].min, format='(%"Min: %0.3e")'), /data, /noclip, color=color_fg
        xyouts, max_loc, max_hist*0.65, string(stats[r].max, format='(%"Max: %0.3e")'), /data, /noclip, color=color_fg
        xyouts, max_loc, max_hist*0.55, string(stats[r].count, format='(%"Count: %d")'), /data, /noclip, color=color_fg
        !X.charsize = 0
        !Y.charsize = 0    
        !P.charsize = 0    
        ptr_free, region_mask 

    endfor

    if (PS) then begin
        device, /close_file
        set_plot, (!VERSION.OS_FAMILY eq 'unix' ) ? 'X' : 'WIN'
    endif
    !P.multi = 0
    
    line = 'Paramter: '+parname+string(09B)+'Units: '+parunits
    line = [line, strjoin(tag_names(stats), string(09B))]
    for r = 0, num_regions-1 do begin
        temp = string(region_names[r], stats[r].mean, stats[r].stdev, stats[r].min, stats[r].max, stats[r].count, $
                      format='(%"%s\t%0.3e\t%0.3e\t%0.3e\t%0.3e\t%d")')
        line = [ line, temp ]
    endfor
    line = [ line, '' ]
    
    if (not widget_info((*state).img_stats_base, /valid_id)) then begin
        stat_base = widget_base(title='Stats', group_leader=event.top)
        stat_txt = widget_text(stat_base, value=line, xsize=max(strlen(line))+45, ysize=n_elements(line)+7, /scroll)
        widget_control, stat_base, /realize, set_uvalue=stat_txt
        (*state).img_stats_base = stat_base
    endif else begin
        widget_control, (*state).img_stats_base, get_uvalue=stat_txt
        widget_control, stat_txt, get_value=existing_stats
        line = [ line, existing_stats ]
        widget_control, stat_txt, set_value=line
    endelse
    
end

;+
; :Description:
;    Handles the displaying of images and maps from an image fit..
;
; :Params:
;    event
;
;
;
; :Author: wtriplett
;-
pro mas_cfit2_vb_display_img_event, event

    widget_control, event.id, get_uvalue=param_state
    if (n_elements(param_state) eq 0) then begin
        return
    endif

    widget_control, event.top, get_uvalue=state
    fr = (*state).img_fit_result
    
    ;; check if the user wants to update the existing map display
    want_update = 0;widget_info(param_state.btn_update_map, /button_set)
    if (want_update) then begin
        if (not ptr_valid(param_state.display_state)) then begin
            want_update = 0B
        endif
    endif
    param_name = widget_info(event.id, /uname)

    parinfo = (*fr).parinfo
    
    parnum = 0L
    for i = 0, n_elements(parinfo)-1 do begin
        if (parinfo[i].fixed) then continue
        
        if (param_name eq parinfo[i].parname) then begin
        
            temp = reform((*fr).param_maps[*,*,parnum])
            widget_control, param_state.min_wid, get_value=min_scl
            widget_control, param_state.max_wid, get_value=max_scl
            image = bytscl(temp, min=min_scl, max=max_scl)
            
            image = mas_cfit2_rotate_flip(image)
            temp  = mas_cfit2_rotate_flip(temp)
            
            if (0 && want_update) then begin
                mmd_replace_image_data, param_state.display_state, image, new_fpvals=temp, new_fprange=[min_scl,max_scl]
            endif else begin
                mas_display_multi, image, fp_vals=temp, standalone=1, $
                                   tab_title=parinfo[i].parname+' map', $
                                   drawstate_handle=ds, fp_range=[min_scl,max_scl]
                param_state.display_state = ds
                widget_control, event.id, set_uvalue=param_state
            endelse
            return
        
        endif
        
        parnum++
        
    endfor
    
    if (param_name eq 'CHI2') then begin
        temp = reform((*fr).chisq_map)
        widget_control, param_state.min_wid, get_value=min_scl
        widget_control, param_state.max_wid, get_value=max_scl
        image = bytscl(temp, min=min_scl, max=max_scl)

        image = mas_cfit2_rotate_flip(image)
        temp  = mas_cfit2_rotate_flip(temp)

        if (0 && want_update) then begin
            mmd_replace_image_data, param_state.display_state, image, new_fpvals=temp, new_fprange=[min_scl,max_scl]
        endif else begin
            mas_display_multi, image, fp_vals=temp, standalone=1, $
                               tab_title='Chisquare map', drawstate_handle=ds, $
                               fp_range=[min_scl,max_scl]
            param_state.display_state = ds
            widget_control, event.id, set_uvalue=param_state
        endelse
    endif
    
end

pro mas_cfit2_gui_update_variable_base_vol, state

    gui = *((*state).gui)
    
    ;; if this happens the maybe the main curvefit window went away?
    if (not widget_info(gui.variable_base_vol, /valid_id)) then return
    
    ;; destroy existing widgets and prepare to recreate
    widget_control, gui.variable_base_vol, get_uvalue=vb_child
    if (widget_info(vb_child, /valid_id)) then begin
        widget_control, vb_child, /destroy
        if (ptr_valid((*state).vol_fit_result)) then begin
            ;; commented out due to annoyance factor.
;            YesNo = dialog_message(['Volume fitting results exist for this scan,', $
;                                    ' do you want to reload them?'], /center, /question)
;            if (YesNo eq 'No') then begin
;                ptr_free, (*state).vol_fit_result
;            endif
        endif
    endif
    
    vb_child = widget_base(gui.variable_base_vol, /column)
    widget_control, gui.variable_base_vol, set_uvalue=vb_child
    
    fr = (*state).vol_fit_result
    if (not ptr_valid(fr)) then return
    
    parinfo = (*fr).parinfo
    parnum  = 0
    for i = 0, n_elements(parinfo)-1 do begin
    
        if (parinfo[i].fixed) then continue
        
        b = widget_base(vb_child, /row)
        parname = parinfo[i].parname
        parunits = parinfo[i].parunits
        
         ;; It was requested to hard limit T1 and T2 to max of 3 S.
        if (parname eq 'T1' or parname eq 'T2') then begin
            min_p = 0.0
            max_p = 6.0
        endif else begin
            min_p = min((*fr).param_maps[*,*,parnum], max=max_p)
        endelse
                
        ;; create min/max display scaling widgets
        min_title = parname+' scale min'
        if (parunits ne '') then min_title += ' ('+parunits+')'
        cfs_display_min = cw_fslider(b, title=min_title, /edit, $
                                     min=min_p*0.8, max=max_p*1.6, value=min_p)
                                     
        max_title = parname+' scale max'
        if (parunits ne '') then max_title += ' ('+parunits+')'
        cfs_display_max = cw_fslider(b, title=max_title, /edit, $
                                     min=min_p*0.8, max=max_p*1.6, value=max_p)
        
        ;; create buttons to show map or hist/stats
        bb = widget_base(b, /column)                             
        display_btn = widget_button(bb, value='Display '+parname+' map', uname=parname, $
                                    event_pro='mas_cfit2_vb_display_vol_event')
                                    
        export_btn   = widget_button(bb, value='Export '+parname+' to .nii', uname=parname, $
                                    event_pro='mas_cfit2_vb_export_vol_event')
        
        param_state = { min_wid: cfs_display_min, $
                        max_wid: cfs_display_max, $
                        display_state: ptr_new() }
        
        widget_control, display_btn, set_uvalue=param_state                
        widget_control, export_btn,  set_uvalue=param_state
        parnum++
        
    endfor 
    
    ;; do the above for the chi2 map
    b = widget_base(vb_child, /row)
    min_p = min((*fr).chisq_map, max=max_p)
    cfs_display_min = cw_fslider(b, title='Chisquare scale min', /edit, $
                                 min=min_p*0.8, max=max_p*1.7, value=min_p)
    cfs_display_max = cw_fslider(b, title='Chisquare scale max', /edit, $
                                 min=min_p*0.8, max=max_p*1.7, value=max_p)
    bb = widget_base(b, /column) 
    display_btn = widget_button(bb, value='Display CHI2 map', uname='CHI2', $
                                event_pro='mas_cfit2_vb_display_vol_event')
                                
    export_btn  = widget_button(bb, value='Export CHI2 map to .nii', uname='CHI2', $
                                event_pro='mas_cfit2_vb_export_vol_event')
    
    param_state = { min_wid: cfs_display_min, $
                    max_wid: cfs_display_max, $
                    display_state: ptr_new() }
    
    widget_control, display_btn, set_uvalue=param_state                
    widget_control, export_btn,  set_uvalue=param_state
    
    widget_control, vb_child, set_uvalue=fr, /realize
    
    xmanager, 'mas_cfit2_vb_display_vol', vb_child, /no_block



end


;+
; :Description:
;    Creates the widgets necessary to control an image fit result.
;
; :Params:
;    state
;
;
;
; :Author: wtriplett
;-
pro mas_cfit2_gui_update_variable_base_img, state

    gui = *((*state).gui)
    
    ;; if this happens the maybe the main curvefit window went away?
    if (not widget_info(gui.variable_base, /valid_id)) then return
    
    ;; destroy existing widgets and prepare to recreate
    widget_control, gui.variable_base, get_uvalue=vb_child
    if (widget_info(vb_child, /valid_id)) then begin
        widget_control, vb_child, /destroy
        if (ptr_valid((*state).img_fit_result)) then begin
            ;; commented out due to annoyance factor.
;            YesNo = dialog_message(['Image fitting results exist for this scan,', $
;                                    ' do you want to reload them?'], /center, /question)
;            if (YesNo eq 'No') then begin
;                ptr_free, (*state).img_fit_result
;            endif
        endif
    endif
    
    vb_child = widget_base(gui.variable_base, /column)
    widget_control, gui.variable_base, set_uvalue=vb_child
    
    fr = (*state).img_fit_result
    if (not ptr_valid(fr)) then return
    
    parinfo = (*fr).parinfo
    parnum = 0
    for i = 0, n_elements(parinfo)-1 do begin
    
        if (parinfo[i].fixed) then continue
        
        b = widget_base(vb_child, /row)
        parname = parinfo[i].parname
        parunits = parinfo[i].parunits
        
        ;; It was requested to hard limit T1 and T2 to max of 3 S.
        if (parname eq 'T1' or parname eq 'T2') then begin
            min_p = 0.0
            max_p = 6.0
        endif else begin
            min_p = min((*fr).param_maps[*,*,parnum], max=max_p)
        endelse
        
        ;; create min/max display scaling widgets
        min_title = parname+' scale min'
        if (parunits ne '') then min_title += ' ('+parunits+')'
        cfs_display_min = cw_fslider(b, title=min_title, /edit, $
                                     min=min_p*0.8, max=max_p*1.6, value=min_p)
                                     
        max_title = parname+' scale max'
        if (parunits ne '') then max_title += ' ('+parunits+')'
        cfs_display_max = cw_fslider(b, title=max_title, /edit, $
                                     min=min_p*0.8, max=max_p*1.6, value=max_p)
        
        ;; create buttons to show map or hist/stats
        bb = widget_base(b, /column)                             
        display_btn = widget_button(bb, value='Display '+parname+' map', uname=parname, $
                                    event_pro='mas_cfit2_vb_display_img_event')
                                    
        stats_btn   = widget_button(bb, value='Display '+parname+' stats', uname=parname, $
                                    event_pro='mas_cfit2_vb_stats_img_event')
        
        ;; a couple of other options
        bbb = widget_base(b, /column, /nonexclusive)
        ;btn_update_map = widget_button(bbb, value='Update existing map', uname=parname+'_update')
        ;widget_control, btn_update_map, /set_button
        btn_indep_scale = widget_button(bbb, value='Independent plot scaling', uname=parname+'_scale')
        
        
        param_state = { min_wid: cfs_display_min, $
                        max_wid: cfs_display_max, $
                        display_state: ptr_new(), $
                        $;btn_update_map: btn_update_map, $
                        btn_indep_scale: btn_indep_scale }
        
        widget_control, display_btn, set_uvalue=param_state                
        widget_control, stats_btn,   set_uvalue=param_state
        parnum++
        
    endfor 
    
    ;; do the above for the chi2 map
    b = widget_base(vb_child, /row)
    min_p = min((*fr).chisq_map, max=max_p)
    cfs_display_min = cw_fslider(b, title='Chisquare scale min', /edit, $
                                 min=min_p*0.8, max=max_p*1.7, value=min_p)
    cfs_display_max = cw_fslider(b, title='Chisquare scale max', /edit, $
                                 min=min_p*0.8, max=max_p*1.7, value=max_p)
    bb = widget_base(b, /column) 
    display_btn = widget_button(bb, value='Display CHI2 map', uname='CHI2', $
                                event_pro='mas_cfit2_vb_display_img_event', $
                                uvalue=[cfs_display_min, cfs_display_max])
                                
    stats_btn   = widget_button(bb, value='Display CHI2 stats', uname='CHI2', $
                                event_pro='mas_cfit2_vb_stats_img_event', $
                                uvalue=[cfs_display_min, cfs_display_max])
    
    bbb = widget_base(b, /column, /nonexclusive)
;    btn_update_map = widget_button(bbb, value='Update existing map', uname=parname+'_update')
;    widget_control, btn_update_map, /set_button
    btn_indep_scale = widget_button(bbb, value='Independent plot scaling', uname='CHI2_scale')

    param_state = { min_wid: cfs_display_min, $
                    max_wid: cfs_display_max, $
                    display_state: ptr_new(), $
                    $;btn_update_map: btn_update_map, $
                    btn_indep_scale: btn_indep_scale }
    
    widget_control, display_btn, set_uvalue=param_state                
    widget_control, stats_btn,   set_uvalue=param_state
    
    widget_control, vb_child, set_uvalue=fr, /realize
    
    xmanager, 'mas_cfit2_vb_display_img', vb_child, /no_block

end

;+
; :Description:
;    Displays the threshold mask given the current threshold percentage.
;    Displays it as a binary image.
;
; :Params:
;    event
;    
; :Author: wtriplett
;-
pro mas_cfit2_gui_show_thrmask_event, event

    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        void = dialog_message(['Unable to compute threshold mask:', $
                               !ERROR_STATE.msg], /error, /center)
        return
    endif
    
    widget_control, event.top, get_uvalue=state
    gui = *((*state).gui)
    
    if (event.id eq gui.btn_img_show_thrmask) then begin
        widget_control, gui.sl_img_fit_threshold, get_value=thr_pct
    endif else begin
        widget_control, gui.sl_vol_fit_threshold, get_value=thr_pct
    endelse
    
    widget_control, event.id, get_uvalue=drawstate
    
    thr_pct /= 100.0
    
    slice = mas_cfit2_get_slice_data()
    dims = size(slice, /dimensions)
    
    thr_mask   = bytarr(dims[0], dims[1]) + 1B

    ;; compute threshold
    thr_pct = n_elements(thr_pct) eq 0 ? 0.0001 : float(thr_pct)
    thr_max = max(abs(slice), dimension=3)
    thr_val = max(thr_max) * thr_pct
  
    ;; update mask
    excl = where(thr_max lt thr_val, n_excl)
    if (n_excl gt 0) then begin
        thr_mask[excl] = 0
    endif
    
    ;; display it
    thr_img = bytscl(thr_mask, min=0, max=1)
    if (not ptr_valid(drawstate)) then begin
        mas_display_multi, thr_img, fp_vals=thr_mask, /standalone, tab_title='Image Fit Threshold Mask', $
                           drawstate_handle=drawstate
        widget_control, event.id, set_uvalue=drawstate
    endif else begin
        mmd_replace_image_data, drawstate, thr_img, new_fpvals=thr_mask
    endelse

end


;+
; :Description:
;    Event handler for misc. image fit widgets.
;
; :Params:
;    event
;
;
;
; :Author: wtriplett
;-
pro mas_cfit2_gui_img_event, event

;    widget_control, event.top, get_uvalue=state
;    
;    gui = *((*state).gui)
;    
;    case event.id of 
;
;            gui.sl_img_fit_threshold: begin
;                widget_control, gui.sl_img_fit_threshold, get_value=val
;                (*state).img_fit_threshold = val
;            end
;            
;            gui.sl_histogram_bins: begin
;                widget_control, gui.sl_histogram_bins, get_value=val
;                (*state).img_histogram_bins = val
;            end
;            
;            gui.btn_img_save_ps: (*state).img_save_ps = event.select
;        else: 
;
;    endcase

end


;+
; :Description:
;    Event handler for misc. ROI related widgets.
;
; :Params:
;    event
;
;
;
; :Author: wtriplett
;-
pro mas_cfit2_gui_roi_event, event

    widget_control, event.top, get_uvalue=state
    
    gui = *((*state).gui)
    
    case event.id of 
    
        gui.dl_model: begin
            (*state).fit_model = event.index
        end
        
        gui.dl_meas_type: begin
            (*state).meas_type = event.index
        end
        
        gui.btn_roi_save_ps: begin
            widget_control, gui.txt_roi_ps_file, sensitive=event.select
            widget_control, gui.btn_roi_browse_ps, sensitive=event.select
            (*state).roi_save_ps = event.select
        end

        gui.btn_img_save_ps: begin
            widget_control, gui.txt_img_ps_file, sensitive=event.select
            widget_control, gui.btn_img_browse_ps, sensitive=event.select
            (*state).img_save_ps = event.select
        end
        
        gui.btn_show_stats: begin
            (*state).roi_show_stats_window = event.select
        end
        
        gui.btn_append_stats: begin
            (*state).roi_append_stats = event.select
        end
        
        gui.btn_no_plot: begin
            (*state).plot_layout = 'NO_PLOT'
        end
        
        gui.btn_multiplot: begin
            (*state).plot_layout = 'MULTIPLOT'
        end

        gui.btn_all_in_one: begin
            (*state).plot_layout = 'ALL_IN_ONE'
            widget_control, gui.btn_combined_y, sensitive=1-event.select
        end

        gui.btn_singleplot: begin
            (*state).plot_layout = 'ONE_PER_PAGE'
        end

        gui.btn_combined_y: begin
            (*state).combined_yscale = event.select
        end

        gui.btn_show_hist: begin
            (*state).show_histogram = event.select
        end
        
        else: begin
        
        end
        
    endcase
    
end

;+
; :Description:
;    Event handler that runs the image fit and produces the results..
;
; :Params:
;    event
;
;
;
; :Author: wtriplett
;-
pro mas_cfit2_gui_run_vol_fit_event, event

    common scan_data, project
    ci = project.ci
    
    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        void = dialog_message(['Unable to perform volume fit:', $
                               !ERROR_STATE.msg], /error, /center)
        return
    endif
    catch, /cancel
                      
    widget_control, event.top, get_uvalue=state
    
    gui = *((*state).gui)

    model_index = widget_info(gui.dl_model, /droplist_select)
    
    model_name = (*state).model_names[model_index]
    model_pro  = (*state).model_pros_img[model_index]
    
    if (model_pro eq '') then return
        
    thr_pct = (*state).vol_fit_threshold
    
    n_slices   = project.imndarray[ci].sdim
    fit_result = ptrarr(n_slices)
    
    old_slice_axis = project.procpramarray[ci].slice_axis
    project.procpramarray[ci].slice_axis = 0
    
    pbar = obj_new('progressbar', text='Volume Fit Running...', color='Red')
    pbar->Start
    pbar_cancelled = 0
    for s = 0, n_slices-1 do begin
        ;; this runs the fit!
        call_procedure, model_pro, fit_result=fr, $
                        for_slice=s, $
                        threshold_pct=thr_pct/100.0, $
                        progressbar=0
        
        fit_result[s] = fr
        
        if (pbar->CheckCancel()) then begin
            pbar_cancelled = 1
            break
        endif else begin
            message, /info, 'Finished slice '+string(s, n_slices-1, format='(%"%03d of %03d")')
        endelse
        
        pbar->update, float(s+1)/n_slices * 100.0
        
    endfor
    
    project.procpramarray[ci].slice_axis = old_slice_axis
    pbar->Destroy
    
    if (pbar_cancelled eq 1) then begin
        for s = 0, n_slices-1 do begin
            ptr_free, fit_result[s]
        endfor
        return
    endif else begin
        parinfo = (*fit_result[0]).parinfo
        param_maps_size = size((*fit_result[0]).param_maps, /dimensions)
        param_maps_vol  = fltarr([param_maps_size[0:1], n_slices, param_maps_size[2]])
        thresh_mask_vol = bytarr([param_maps_size[0:1], n_slices])
        chisq_map_vol   = fltarr([param_maps_size[0:1], n_slices])
        
        for s = 0, n_slices-1 do begin
            param_maps_vol[*,*,s,*] = (*fit_result[s]).param_maps
            thresh_mask_vol[*,*,s]  = (*fit_result[s]).thresh_mask
            chisq_map_vol[*,*,s]    = (*fit_result[s]).chisq_map
            ptr_free, fit_result[s]
        endfor
    endelse
    
    (*state).vol_fit_result = ptr_new({ parinfo: parinfo, $
                                        param_maps: temporary(param_maps_vol), $
                                        thresh_mask: temporary(thresh_mask_vol), $
                                        chisq_map: temporary(chisq_map_vol) }, /no_copy)

    mas_cfit2_gui_update_variable_base_vol, state
   
end


;+
; :Description:
;    Event handler that runs the image fit and produces the results..
;
; :Params:
;    event
;
;
;
; :Author: wtriplett
;-
pro mas_cfit2_gui_run_img_fit_event, event

    common scan_data, project
    ci = project.ci
    
    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        void = dialog_message(['Unable to perform image fit:', $
                               !ERROR_STATE.msg], /error, /center)
        return
    endif
                      
    widget_control, event.top, get_uvalue=state
    
    gui = *((*state).gui)

    widget_control, gui.variable_base, get_uvalue=vb_child
    if (widget_info(vb_child, /valid_id)) then begin
        widget_control, vb_child, /destroy
    endif

    model_index = widget_info(gui.dl_model, /droplist_select)
    
    model_name = (*state).model_names[model_index]
    model_pro  = (*state).model_pros_img[model_index]
    
    if (model_pro eq '') then return
    
    if (ptr_valid((*state).img_fit_result)) then ptr_free, (*state).img_fit_result
    
    thr_pct = (*state).img_fit_threshold
    
    ;; this runs the fit!
    call_procedure, model_pro, fit_result=fit_result, $
                    threshold_pct=thr_pct/100.0, $
                    progressbar=1
    
    (*state).img_fit_result = fit_result

    mas_cfit2_gui_update_variable_base_img, state
    
end

;+
; :Description:
;    Event handler that runs the ROI fit.
;
; :Params:
;    event
;
;
;
; :Author: wtriplett
;-
pro mas_cfit2_gui_run_roi_fit_event, event

    common scan_data, project
    ci = project.ci
    
    widget_control, event.top, get_uvalue=state

    gui = *((*state).gui)
    
    no_plot    = (*state).plot_layout eq 'NO_PLOT' ? 1 : 0
    multiplot  = (*state).plot_layout eq 'MULTIPLOT' ? 1 : 0
    all_in_one = (*state).plot_layout eq 'ALL_IN_ONE' ? 1 : 0
    if ( (*state).plot_layout eq 'ONE_PER_PAGE' ) then begin
        multiplot  = 0
        all_in_one = 0
    endif
    
    combined_yscale = (*state).combined_yscale
    no_stats        = 1-(*state).roi_show_stats_window            

    if ((*state).roi_save_ps eq 1 and (*state).roi_ps_output_file ne '') then begin
        psfile = (*state).roi_ps_output_file
    endif

    model_index = (*state).fit_model; widget_info(gui.dl_model, /droplist_select)
    
    model_name = (*state).model_names[model_index]
    model_pro  = (*state).model_pros_roi[model_index]

    meas_type = (*state).meas_type
    
    mas_cfit2_roi_fit, model_pro, psfile=psfile, $
                       multiplot=multiplot, $
                       all_in_one=all_in_one, $
                       combined_yscale=combined_yscale, $
                       no_legend=no_legend, $
                       no_plot=no_plot, $
                       roi_stats=roi_stats, $
                       measurement_type=meas_type
    
    if (n_elements(roi_stats) eq 0) then return
    
    ;; produces the stats
    if ((*state).roi_show_stats_window) then begin
        if (not widget_info((*state).roi_stats_base, /valid_id)) then begin
            stat_base = widget_base(title='ROI Stats', group_leader=event.top)
            stat_txt = widget_text(stat_base, value=roi_stats, $
                                   xsize=max(strlen(roi_stats))+20, $
                                   ysize=n_elements(roi_stats)+7, /scroll)
            widget_control, stat_base, /realize, set_uvalue=stat_txt
            (*state).roi_stats_base = stat_base
        endif else begin
            widget_control, (*state).roi_stats_base, get_uvalue=stat_txt
            if ((*state).roi_append_stats eq 1) then begin
                widget_control, stat_txt, get_value=existing_stats
                roi_stats = [ roi_stats, '', existing_stats ]
            endif
            widget_control, stat_txt, set_value=roi_stats
        endelse
    endif
end

;+
; :Description:
;    Procedure to browser and obtain an file path to a destination file
;    for postscript plotting.
;
; :Params:
;    event
;
;
;
; :Author: wtriplett
;-
pro mas_cfit2_gui_browse_roi_psfile_event, event

    widget_control, event.top, get_uvalue=state
    
    gui = *((*state).gui)
    
    file = dialog_pickfile(path=(*state).current_path, $
                           get_path=path, $
                           title='Please select an output file.', $
                           filter='*.ps', default_extension='ps', /overwrite_prompt)
    if (file eq '') then return
    
    if (file_basename(file) eq '.ps') then begin
        file = file_dirname(file, /mark_directory) + 'cf2_plot.ps'
    endif
    
    (*state).roi_ps_output_file = file
    (*state).current_path = path                           
    widget_control, gui.txt_roi_ps_file, set_value=file
    
end

pro mas_cfit2_gui_browse_img_psfile_event, event

    widget_control, event.top, get_uvalue=state
    
    gui = *((*state).gui)
    
    file = dialog_pickfile(path=(*state).current_path, $
                           get_path=path, $
                           title='Please select an output file.', $
                           filter='*.ps', default_extension='ps', /overwrite_prompt)
    if (file eq '') then return
    
    if (file_basename(file) eq '.ps') then begin
        file = file_dirname(file, /mark_directory) + 'cf2_plot.ps'
    endif
    
    (*state).img_ps_output_file = file
    (*state).current_path = path                           
    widget_control, gui.txt_img_ps_file, set_value=file
    
end

pro mas_cfit2_gui_cleanup, tlb

    common scan_data, project
    
    widget_control, tlb, get_uvalue=state
    
    ptr_free, project.curvefit2_gui
    ptr_free, (*state).img_fit_result
    ptr_free, state

end

pro mas_cfit2_gui_event, event

    widget_control, event.top, get_uvalue=state
    
    gui = *((*state).gui)
    
    case event.id of 
            gui.s_array: select_data_gui, index = (*state).project_ci
            
            gui.sl_img_fit_threshold: begin
                widget_control, gui.sl_img_fit_threshold, get_value=val
                (*state).img_fit_threshold = val
            end
            
            gui.sl_vol_fit_threshold: begin
                widget_control, gui.sl_vol_fit_threshold, get_value=val
                (*state).vol_fit_threshold = val
            end
            
            gui.sl_histogram_bins: begin
                widget_control, gui.sl_histogram_bins, get_value=val
                (*state).img_histogram_bins = val
            end
            
            gui.btn_img_save_ps: (*state).img_save_ps = event.select
        else: 

    endcase

end

pro mas_cfit2_gui_redraw, void

    common scan_data, project
    ci = project.ci

    pgui = project.curvefit2_gui
    if (not ptr_valid(pgui)) then return
    
    new_state = project.procpramarray[ci].curvefit2_state
    
    if (not ptr_valid(new_state)) then begin
        new_state = mas_cfit2_make_default_state()
        project.procpramarray[ci].curvefit2_state = new_state
    endif else begin
        widget_control, (*pgui).tlb, get_uvalue=old_state
        if (ptr_valid(old_state)) then begin
            if ((*old_state).project_id eq (*new_state).project_id) then return
        endif
    endelse

    (*new_state).gui = pgui
    
    gui = *pgui
    
    widget_control, gui.tlb, set_uvalue=new_state
    
    widget_control, gui.dl_model, set_droplist_select=(*new_state).fit_model
    widget_control, gui.dl_meas_type, set_droplist_select=(*new_state).meas_type
    widget_control, gui.btn_show_stats, set_button=(*new_state).roi_show_stats_window
    widget_control, gui.btn_append_stats, set_button=(*new_state).roi_append_stats
    case (*new_state).plot_layout of
        'MULTIPLOT': widget_control, gui.btn_multiplot, /set_button
        'ALL_IN_ONE': widget_control, gui.btn_all_in_one, /set_button
        'ONE_PER_PAGE': widget_control, gui.btn_singleplot, /set_button
        'NO_PLOT': widget_control, gui.btn_no_plot, /set_button
    endcase
    widget_control, gui.btn_combined_y, set_button=(*new_state).combined_yscale
    widget_control, gui.btn_roi_save_ps, set_button=(*new_state).roi_save_ps
    widget_control, gui.txt_roi_ps_file, set_value=(*new_state).roi_ps_output_file
    
    widget_control, gui.sl_img_fit_threshold, set_value=(*new_state).img_fit_threshold
    widget_control, gui.sl_histogram_bins, set_value=(*new_state).img_histogram_bins
    widget_control, gui.btn_img_save_ps, set_button=(*new_state).img_save_ps
    widget_control, gui.txt_img_ps_file, set_value=(*new_state).img_ps_output_file
    
    mas_cfit2_gui_update_variable_base_img, new_state
    mas_cfit2_gui_update_variable_base_vol, new_state
end

pro mas_cfit2_gui

    common scan_data, project
    common common_widgets

    ci = project.ci
    if (project.procpramarray[ci].state_1 eq 0) then begin
        mas_load_state_1
    endif
    
    ;; models of fitting, all 3 arrays need to be updated when a model is added or removed!
    ;; models of fitting, all 3 arrays need to be updated when a model is added or removed!
    model_names = [ 'Please Select', $
                    'T1: Inversion Recovery SE', $
                    'T1: Inversion Recovery SE (ABS w/ Polarity Correction)', $
                    'T1: Inversion Recovery SE (ABS w/ Powell Optimizer)', $
                    'T1: Saturation Recovery SE', $
                    'T2: Variable TE w/o Baseline', $
                    'T2: Variable TE w/ Baseline', $
                    'ADC: Multi B-value' ]
    
    model_pros_roi = [ '', $
                       'mas_cfit2_ROI_FIT_T1_VarTI', $
                       'mas_cfit2_ROI_FIT_T1_VarTI_ABS', $
                       'mas_cfit2_ROI_FIT_T1_VarTI_POWELL', $
                       'mas_cfit2_ROI_FIT_T1_VarTR', $
                       'mas_cfit2_ROI_FIT_T2_VarTE', $
                       'mas_cfit2_ROI_FIT_T2_wBaseline_VarTE', $
                       'mas_cfit2_ROI_FIT_ADC_VarB' ]
                       
    model_pros_img = [ '', $
                       'mas_cfit2_IMG_FIT_T1_VarTI', $
                       'mas_cfit2_IMG_FIT_T1_VarTI_ABS', $
                       'mas_cfit2_IMG_FIT_T1_VarTI_POWELL', $
                       'mas_cfit2_IMG_FIT_T1_VarTR', $
                       'mas_cfit2_IMG_FIT_T2_VarTE', $
                       'mas_cfit2_IMG_FIT_T2_wBaseline_VarTE', $
                       'mas_cfit2_IMG_FIT_ADC_VarB' ]

    case 1 of
    
        ptr_valid(project.imndarray[ci].inversion_time_ptr): begin
            model_selection = 2
        end
        ptr_valid(project.imndarray[ci].echo_time_ptr): begin
            model_selection = 5
        end
        ptr_valid(project.imndarray[ci].rep_time_ptr): begin
            model_selection = 4
        end
        ptr_valid(project.imndarray[ci].bval_array): begin
            model_selection = 7
        end
        else: model_selection = 0
        
    endcase
    
                       
    base = widget_base(group_leader=WID_BASE_MAIN, title='Curve Fit 2 (BETA)', /column)
    
    eh_roi = 'mas_cfit2_gui_roi_event'
    eh_img = 'mas_cfit2_gui_img_event'
    
    type_base = widget_base(base, /row)
    lbl_type  = widget_label(type_base, value='Curve fit model:')
    dl_model  = widget_droplist(type_base, value=model_names, event_pro=eh_roi)    
    widget_control, dl_model, set_droplist_select=model_selection
    s_array   = widget_button(type_base,value = 'Select Data to Fit')
    
    tab_base = widget_tab(base)
    
    tab_roi = widget_base(tab_base, title='ROI Curve Fit', /column)

    type_base      = widget_base(tab_roi, /row)
    lbl_roi_method = widget_label(type_base, value='Fit Method:')
    dl_meas_type   = widget_droplist(type_base, event_pro=eh_roi, $
                                     value=['mean(measured data)', $
                                            'mean(abs(measured_data))', $
                                            'mean(signal(measured_data))'])

    out_base       = widget_base(tab_roi, /column, /nonexclusive)
    btn_show_stats = widget_button(out_base, value='Show statistics window', event_pro=eh_roi)
    widget_control, btn_show_stats, /set_button
    btn_append_stats = widget_button(out_base, value='Append results to statistics window', event_pro=eh_roi)
    widget_control, btn_append_stats, /set_button
    
    btn_show_hist  = 1;widget_button(out_base, value='Display parameter histograms')
    btn_combined_y = widget_button(out_base, value='Unified y-scale for all ROIs', event_pro=eh_roi)
    out_base       = widget_base(tab_roi, /column, /exclusive)
    btn_no_plot    = widget_button(out_base, value="Don't show plot window", event_pro=eh_roi)
    btn_singleplot = widget_button(out_base, value='One plot per page (many windows)', event_pro=eh_roi)
    btn_multiplot  = widget_button(out_base, value='Multiple single plots on page', event_pro=eh_roi)
    btn_all_in_one = widget_button(out_base, value='All plots on same page', event_pro=eh_roi)
    widget_control, btn_multiplot, /set_button
    
    colbase = widget_base(tab_roi, /row)
    b1  = widget_base(colbase, /row, /nonexclusive, xpad=0, ypad=0)
    btn_roi_save_ps = widget_button(b1, value='Save plots to PS file:', event_pro=eh_roi)
    b2  = widget_base(colbase, /row, xpad=0, ypad=0)
    txt_roi_ps_file = widget_text(b2, xsize=60, /editable)
    btn_roi_browse_ps = widget_button(b2, value='Browse...', event_pro='mas_cfit2_gui_browse_roi_psfile_event')
    
    btn_base  = widget_base(tab_roi, /row)
    btn_run_roi_fit  = widget_button(btn_base, value='Run Fit', event_pro='mas_cfit2_gui_run_roi_fit_event')
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    tab_img = widget_base(tab_base, title='IMAGE Curve Fit', /column)

    colbase = widget_base(tab_img, /row)
    sl_img_fit_threshold = cw_fslider(colbase, title='Signal Threshold (%)', min=0.01, max=100, value=0.10, /edit)
    
    btn_base  = widget_base(colbase, /column)
    btn_img_show_thrmask = widget_button(btn_base, value='Show Threshold Mask', event_pro='mas_cfit2_gui_show_thrmask_event')
    btn_run_img_fit = widget_button(btn_base, value='Run Image Fit', event_pro='mas_cfit2_gui_run_img_fit_event')
    
    sl_base = widget_base(colbase, /row)
    sl_histogram_bins = widget_slider(sl_base, title='Histogram Bins', min=10, max=200, value=50)

    variable_base = widget_base(tab_img, /column)
    variable_base_child = widget_base(variable_base, /column)
    widget_control, variable_base, set_uvalue=variable_base_child

    ps_base           = widget_base(tab_img, /row)
    btn_base          = widget_base(ps_base, xpad=0, ypad=0, /row, /nonexclusive)
    btn_img_save_ps   = widget_button(btn_base, value='Save plots to PS file:')
    txt_img_ps_file   = widget_text(ps_base, value='', xsize=55, /edit)
    btn_img_browse_ps = widget_button(ps_base, value='Browse...', uvalue=txt_img_psfile, $
                                      event_pro='mas_cfit2_gui_browse_img_psfile_event')

    sdim = project.imndarray[ci].sdim
    
    if (sdim gt 1) then begin
        tab_volume = widget_base(tab_base, title='Volume Curve Fit', /column)
        colbase = widget_base(tab_volume, /row)
        sl_vol_fit_threshold = cw_fslider(colbase, title='Signal Threshold (%)', min=0.01, max=100, value=0.10, /edit)
        
        btn_base  = widget_base(colbase, /column)
        btn_vol_show_thrmask = widget_button(btn_base, value='Show Threshold Mask', event_pro='mas_cfit2_gui_show_thrmask_event')
        btn_run_vol_fit = widget_button(btn_base, value='Run Volume Fit', event_pro='mas_cfit2_gui_run_vol_fit_event')
        btn_base   = widget_base(tab_volume, /row)
                
        variable_base_vol = widget_base(tab_volume, /column)
        variable_base_vol_child = widget_base(variable_base_vol, /column)
        widget_control, variable_base_vol, set_uvalue=variable_base_vol_child
    endif else begin
        btn_run_vol_fit       = 0
        variable_base_vol     = 0
        ;;txt_vol_load_nii_mask = 0
        ;;btn_vol_load_nii_mask = 0
        sl_vol_fit_threshold  = 0
        btn_vol_show_thrmask  = 0
        ;;btn_vol_apply_nii_mask = 0
    endelse
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    widget_control, base, /realize
    
    gui =  ptr_new({ tlb: base, $
                     stats_base: 0L, $
                     dl_model: dl_model, $
                     s_array:s_array, $
                     dl_meas_type: dl_meas_type, $
                     btn_show_stats: btn_show_stats, $
                     btn_append_stats: btn_append_stats, $
                     btn_show_hist: btn_show_hist, $
                     btn_multiplot: btn_multiplot, $
                     btn_no_plot: btn_no_plot, $
                     btn_singleplot: btn_singleplot, $
                     btn_all_in_one: btn_all_in_one, $
                     btn_combined_y: btn_combined_y, $
                     btn_roi_save_ps: btn_roi_save_ps, $
                     btn_roi_browse_ps: btn_roi_browse_ps, $
                     txt_roi_ps_file: txt_roi_ps_file, $
                     btn_run_roi_fit: btn_run_roi_fit, $
                     $
                     variable_base: variable_base, $
                     dl_noise_roiset: 0, $;dl_noise_roiset, $            ;; unused
                     dl_noise_region: 0, $;dl_noise_region, $            ;; unused
                     btn_refresh_roi_list: 0, $; btn_refresh_roi_list, $ ;; unused
                     sl_img_fit_threshold: sl_img_fit_threshold, $
                     sl_histogram_bins: sl_histogram_bins, $
                     btn_img_save_ps: btn_img_save_ps, $
                     btn_img_browse_ps: btn_img_browse_ps, $
                     txt_img_ps_file: txt_img_ps_file, $
                     btn_img_show_thrmask: btn_img_show_thrmask, $
                     $
                     $;;txt_vol_load_nii_mask: txt_vol_load_nii_mask, $
                     $;;btn_vol_browse_nii_mask: btn_vol_browse_nii_mask, $
                     $;;btn_vol_apply_nii_mask: btn_vol_apply_nii_mask, $
                     sl_vol_fit_threshold: sl_vol_fit_threshold, $
                     btn_vol_show_thrmask: btn_vol_show_thrmask, $
                     btn_run_vol_fit: btn_run_vol_fit, $
                     variable_base_vol: variable_base_vol })

    state = mas_cfit2_make_default_state()
    (*state).gui = gui
    (*state).project_ci =  ci                  ;; the project scan id
    (*state).project_id =  randomn(systime(1)) ;; this is a unique id to identify this scan
    (*state).model_names =  model_names        ;; the names of the available models
    (*state).model_pros_roi =  model_pros_roi  ;; the procedures for each model's ROI fit
    (*state).model_pros_img =  model_pros_img  ;; the procedures for each models's image fit
    (*state).show_histogram =  0               ;; 1,0 - whether or not to display the histogram
    (*state).plot_layout =  'MULTIPLOT'        ;; a string indicating the type of plot to show
    (*state).combined_yscale =  0              ;; 1,0 - whether or not to unify the y-scale between regions
    (*state).roi_show_stats_window =  1        ;; 1,0 - show the ROI stats window when fit is run
    (*state).roi_save_ps =  0                  ;; 1,0 - save the plots to postscript
    (*state).roi_ps_output_file =  ''          ;; postscript output file name
    (*state).roi_stats_base =  0L              ;; the widget ID of the ROI stats window
    (*state).roi_append_stats =  1B            ;; 1,0 - append stats or replace window contents with new stats
    (*state).fit_model =  model_selection      ;; the current fit model (index into models array)
    (*state).meas_type =  0                    ;; the numeric measurement type (see: mas_cfit2_roi_get_data)
    (*state).img_histogram_bins =  50          ;; number of bins in image histogram fit
    (*state).img_fit_threshold =  0.10         ;; the pct threshold in image fitting
    (*state).img_fit_result =  ptr_new()       ;; poiner to the fit result (last run)
    (*state).img_stats_base =  0L              ;; widget id of the window displaying the image fit stats
    (*state).img_save_ps =  0                  ;; 1,0 - save the image histogram plot to postscript
    (*state).img_ps_output_file =  ''          ;; the path to the image fit postscript output file
    (*state).vol_fit_result =  ptr_new()       ;; pointer to hold volume fit results
    (*state).vol_fit_threshold =  0.10         ;; the pct threshold in volume fitting
    (*state).current_path = project.current_path
    
;    state = { , $
;              project_ci: ci, $                  ;; the project scan id
;              project_id: randomn(systime(1)), $ ;; this is a unique id to identify this scan
;              model_names: model_names, $        ;; the names of the available models
;              model_pros_roi: model_pros_roi, $  ;; the procedures for each model's ROI fit
;              model_pros_img: model_pros_img, $  ;; the procedures for each models's image fit
;              show_histogram: 0, $               ;; 1,0 - whether or not to display the histogram
;              plot_layout: 'MULTIPLOT', $        ;; a string indicating the type of plot to show
;              combined_yscale: 0, $              ;; 1,0 - whether or not to unify the y-scale between regions
;              roi_show_stats_window: 1, $        ;; 1,0 - show the ROI stats window when fit is run
;              roi_save_ps: 0, $                  ;; 1,0 - save the plots to postscript
;              roi_ps_output_file: '', $          ;; postscript output file name
;              roi_stats_base: 0L, $              ;; the widget ID of the ROI stats window
;              roi_append_stats: 1B, $            ;; 1,0 - append stats or replace window contents with new stats
;              fit_model: model_selection, $      ;; the current fit model (index into models array)
;              meas_type: 0, $                    ;; the numeric measurement type (see: mas_cfit2_roi_get_data)
;              img_histogram_bins: 50, $          ;; number of bins in image histogram fit
;              img_fit_threshold: 0.10, $         ;; the pct threshold in image fitting
;              img_fit_result: ptr_new(), $       ;; poiner to the fit result (last run)
;              img_stats_base: 0L, $              ;; widget id of the window displaying the image fit stats
;              img_save_ps: 0, $                  ;; 1,0 - save the image histogram plot to postscript
;              img_ps_output_file: '', $          ;; the path to the image fit postscript output file
;              vol_fit_result: ptr_new(), $       ;; pointer to hold volume fit results
;              current_path: project.current_path } 
    
;    pstate = ptr_new(state, /no_copy)
    widget_control, base, set_uvalue=state
    
    project.procpramarray[ci].curvefit2_state = state
    project.curvefit2_gui = (*state).gui

    xmanager, 'mas_cfit2_gui', base, /no_block, cleanup='mas_cfit2_gui_cleanup'
    
end
