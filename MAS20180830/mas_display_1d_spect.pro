;+
; :Description:
;    Inversion recovery fit function.
;
; :Params:
;    TI - inversion times
;    A  - model parameters
;
; :Author: wtriplett
;-
function mas_1ds_inv_rec_function, TI, A

    ;; A[0] - steady state signal
    ;; A[1] - T1
    
    return, A[0]*( 1.0 - 2.0*exp(-TI/A[1]) ) 

end

;+
; :Description:
;    T2 fit parameters.
;
; :Params:
;    TE - Echo times
;    A  - model parameters
;
; :Author: wtriplett
;-
function mas_1ds_t2_function, TE, A

    ;; A[0] - steady state signal
    ;; A[1] - T2
    
    return, A[0]*( exp(-TE/A[1]) ) 

end

function mas_1ds_fft, data

    dims = size(data)
    if (dims[0] gt 1) then begin
        adim = dims[2]
    endif else begin
        adim = 1
    endelse
    
    acq_size = dims[1]
    
    if (adim gt 1) then begin
        fdata = shift(fft(data, dimension=1), acq_size/2,0)
    endif else begin
        fdata = shift(fft(data, dimension=1), acq_size/2)
    endelse
    
    fdata = reverse(fdata)

    return, fdata

end

function mas_1ds_compute_fwhm, state, which, index = index, which_fid = which_fid

    if n_elements(which_fid) ne 0 then curr_fid = which_fid else curr_fid = (*state).curr_fid
    (*state).ofq_plot->getCursorInfo, cursor01_info=c1, cursor02_info=c2
    
    if (n_elements(which) ne 0 && which eq 2) then begin
        data_index = c2.data_index
    endif else begin
        data_index = c1.data_index
    endelse
    if keyword_set(index) then data_index = index
    peak_val = real_part((*state).fq_data[data_index, curr_fid])
    half_peak = 0.5 * peak_val
    
    lt_fwhm_pt = data_index
    while (lt_fwhm_pt gt 0) do begin
        tempval = real_part((*state).fq_data[lt_fwhm_pt, curr_fid])
        if (tempval le half_peak) then break
        lt_fwhm_pt--
    endwhile
    m_lt = -( tempval-real_part((*state).fq_data[lt_fwhm_pt+1, curr_fid]) )
    k_lt = half_peak/m_lt
    x_lt = (((half_peak - tempval) + m_lt*lt_fwhm_pt)/m_lt)
                    
    rt_fwhm_pt = data_index
    while (rt_fwhm_pt lt (*state).acq_size) do begin
        tempval = real_part((*state).fq_data[rt_fwhm_pt,curr_fid])
        if (tempval le half_peak) then break
        rt_fwhm_pt++
    endwhile
    m_rt = ( tempval-real_part((*state).fq_data[rt_fwhm_pt-1, curr_fid]) )
    k_rt = half_peak/m_rt
    x_rt = (((half_peak - tempval) + m_rt*rt_fwhm_pt)/m_rt)
    
    fwhm = abs(x_rt-x_lt)/(*state).acq_size * (*state).spect_width

    return, fwhm
        
end


function mas_1ds_truncating_shift, data, shift_amt

    data_size = size(data)
    
    if (shift_amt ne 0) then begin
        if (data_size[0] gt 1) then begin
            data = shift(data, shift_amt, 0)
            if (shift_amt lt 0) then begin
                data[data_size[1]-1+shift_amt:*,*] = complex(0,0)
            endif else begin
                data[0:shift_amt-1,*] = complex(0,0)
            endelse
        endif else begin
            data = shift(data, shift_amt)
            if (shift_amt lt 0) then begin
                data[data_size[1]-1+shift_amt:*] = complex(0,0)
            endif else begin
                data[0:shift_amt-1] = complex(0,0)
            endelse
        endelse
    endif

    return, data
    
end

function mas_1ds_apply_baseline_corr, data, npts

    if (npts gt 0) then begin
        rmean = mean(data[n_elements(data)-1-npts:*])
        data = data - replicate(rmean, n_elements(data))    
    endif
    
    return, data
    
end


function mas_1ds_apply_ts_weight, data, type, acq_time, val

    data_size = size(data)
    if (data_size[0] gt 1) then begin
        num_fids = data_size[2]
    endif else begin
        num_fids = 1
    endelse
    
    acq_size = float(data_size[1])
    
    KdeltaT = findgen(acq_size)/acq_size * acq_time
    W = float(val[0])
    
    if (type eq 'EM') then begin
        corr = exp(-!PI*W*KdeltaT)
    endif else if (type eq 'GM') then begin
        corr = exp(-W*KdeltaT^2)
    endif else begin
        return, data
    endelse
    
    if (num_fids gt 1) then begin
        corr = rebin(corr, n_elements(corr), num_fids, /sample)
    endif
    
    data *= complex(corr)
 
    return, data

end

function mas_1ds_compute_phase_corr, state, slope_intercept=slope_intercept

    ;;(*state).ofq_plot->getCursorInfo, cursor01_info=c1
    ;;pivot = float(c1.data_index)/(*state).acq_size
    
    pivot = float((*state).fq_phase_pivot)/(*state).acq_size
    
    y1 = float((*state).fq_0th_corr_widval)*!DTOR
    x1 = pivot
    m  = float((*state).fq_1st_corr_widval)*!DTOR
    b = y1 - m*x1

    phase_corr = [ m, b ]

    return, phase_corr

end

function mas_1ds_convert_data_to_window, state, x

    (*state).ofq_window->getProperty, dimensions=window_dim
    return, (x + (*state).zero_ref)/(*state).spect_width * window_dim[0]

end

function mas_1ds_convert_window_to_data, state, x

    (*state).ofq_window->getProperty, dimensions=window_dim
    return, x/float(window_dim[0]) * (*state).spect_width - (*state).zero_ref

end

pro mas_1ds_compute_fq_baseline, state, region_size=region_size, A=A, B=B, $
                                 checkcancel=checkcancel

    if (keyword_set(region_size)) then begin
        region_size = long(region_size)
    endif else begin
        region_size=32
    endelse
    
    checkcancel = 0
    
    prog = obj_new('progressbar', title='Estimating Baseline', text='Estimating baseline...')
    prog->start
    progmax = (*state).num_fids * 30.0
    
    for fid = 0, (*state).num_fids-1 do begin
    
        prog->SetProperty, text=string(fid, (*state).num_fids,$
                                       format='(%"Estimating Baseline, FID %03d of %03d")')
                                       
        fq_data = reform((*state).fq_data[*,fid])
        
        regions = dcomplexarr(region_size, n_elements(fq_data)/region_size)
        nregions = (size(regions, /dimensions))[1]
        for i = 0, nregions-1 do begin 
            regions[*,i] = (fq_data)[i*region_size:(i+1)*region_size-1]
        endfor
        r_regions = real_part(regions)
        mean_regions = total(r_regions,1)/region_size
        expanded_rmean = transpose(rebin(mean_regions, nregions, region_size, /sample))
        var_regions = total((r_regions - expanded_rmean)^2, 1)/(region_size-1)
        min_mean = where(abs(mean_regions) eq min(abs(mean_regions)))
        psuedo_nvar = sqrt(var_regions[min_mean])
        
        s = psuedo_nvar[0]
        spec = real_part(fq_data)
        L = n_elements(spec)
        
        A = (5e-8^(0.25)*L)^4
        B = 1.25D
        
        A = ((5e-8)^(0.25)*L)^4
        B = 1.25D
            
        Bs= -B/s
        As= -A/s
        
        bd=(dblarr(L)+1)*median(spec)
        bd0=spec
        nm=norm(bd-bd0)
        nm0=1e20
        
        M0 = -(dblarr(L)+1)/As
        e  = (dblarr(L)+1)
        D0 = dblarr(L,L)
        
        for i = 0, L-1 do begin
            D0[i,i] = 12.0
            if (i lt (L-1)) then begin
                D0[i,i+1] = -8 & D0[i+1,i] = -8
            endif
            
            if (i lt (L-2)) then begin
                D0[i,i+2] =  2 & D0[i+2,i] =  2
            endif
            
        endfor
        
        D0[0,0] = 2      &  D0[L-1,L-1] = 2
        D0[1,0] = -4     &  D0[0,1] = -4
        D0[L-1,L-2] = -4 &  D0[L-2,L-1] = -4
        D0[1,1] = 10     &  D0[L-2,L-2] = 10
        It=0
            
        while (nm gt 1e-3 and It lt 30) do begin
            prog->update, ((fid*30.0) + It)/progmax * 100.0
            if (prog->CheckCancel()) then begin
                checkcancel=1
                (*state).fq_baseline[*] = 0.0
                goto, FINISHED
            endif
            It++
            M = M0
            D = D0
            bd0 = bd
            nm0 = nm
            
            for i = 0, L-1 do begin
                if (bd[i] gt spec[i]) then begin
                    M[i]   = M[i]+2*Bs*spec[i]/As
                    D[i,i] = D[i,i]+2*Bs/As
                endif
            endfor
            
            Dsp = sprsin(D, /double)
            bd  = linbcg(Dsp, M, bd0, /double)
            nm  = norm(bd0-bd);
            print, "ITER: ", it, "Norm: ", nm
        end
    
        (*state).fq_baseline[*,fid] = bd
        
    endfor
    
    FINISHED:
    prog->Destroy
    print, "Ok Done"


end


pro mas_1ds_run_processing_chain, state, data

    my_data = data
    
    my_data = mas_1ds_apply_baseline_corr(my_data, (*state).baseline_corr)
    
    my_data = mas_1ds_truncating_shift(my_data, (*state).ts_shift)
    
    if ((*state).weight_active eq 1) then begin
        my_data = mas_1ds_apply_ts_weight(my_data, $
                                          (*state).weight_type, $
                                          (*state).acq_time, $
                                          (*state).weight_param)
    endif
    
    (*state).ts_data = my_data
    (*state).ots_plot->setData, (*state).ts_data
    
    (*state).fq_data = mas_1ds_fft(my_data)
    ;; fq will be set during phase correction
    ;;(*state).ofq_plot->setData, (*state).fq_data

    mas_1ds_apply_phase_corr, state, /global
 
    if ((*state).fq_baseline_correction_on) then begin
        (*state).fq_data -= (*state).fq_baseline
        (*state).ofq_plot->setData, (*state).fq_data
    endif
    
end

pro mas_1ds_window_init, data, title=title, tlb=tlb

    common scan_data, project
    common common_widgets
    
    ci = project.ci
    color_orange = [255,128,0]
    color_blue   = [0,0,255]
    
    if (n_elements(title) eq 0) then begin
        title = 'mas_1d_spect'
    endif
    
    sw = project.imndarray[ci].spect_spectral_width
    sf = project.imndarray[ci].spect_bf1
    at = project.imndarray[ci].spect_acq_time
    adim_start = project.procpramarray[ci].adim_start
    adim = project.imndarray[ci].adim
    
    acq_size = project.imndarray[ci].spect_acq_size
    
    data = *project.dataarray[ci].state1
    
    max_data = max(abs(data), max_ind)
    print, "Max FID of ("+string(max_data, format='(G0)')+') occurs at index: '+string(max_ind, format='(G0)')
    tlb = widget_base(title=title, column=2, ypad=0, mbar=menu_bar)
    file = widget_button(menu_bar, value='File', /menu)
    view = widget_button(menu_bar, value='View', /menu)
    analyze = widget_button(menu_bar, value='Analyze', /menu)
    
    view_frq = widget_button(view, value='Spectrum', /menu)
    view_ts = widget_button(view, value='Free Induction Decay', /menu)
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    left_panel = widget_base(tlb, /column)
    main_panel = widget_base(tlb, /column)
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    file_save_tif = widget_button(file, value="Save main plots as TIF...", $
                                  event_pro='mas_1ds_save_main_plots', $
                                  uname='file_save_tif')
    file_save_ps = widget_button(file, value="Save main plots as Postscript...", $
                                  event_pro='mas_1ds_save_main_plots', $
                                  uname='file_save_ps')                              
    file_exit = widget_button(file, value="Exit", event_pro='mas_1ds_exit_event')
    
    apply_freq_baseline = widget_button(analyze, value="Turn ON Frequency Baseline Correction", $
                                        event_pro='mas_1ds_apply_baseline_event')
    show_multiplot = widget_button(analyze, value="Plot Multiple Spectra", $
                                   event_pro='mas_1ds_plot_spect_multi_window')
    curve_fit = widget_button(analyze, value='Curve fitting', /menu)
    curve_fit_IR = widget_button(curve_fit, accelerator='CTRL+SHIFT+1', value='T1 (Inversion Recovery)', event_pro='mas_1ds_cfit_T1_TI')
    curve_fit_T2 = widget_button(curve_fit, accelerator='CTRL+SHIFT+2', value='T2 (Multi-Echo)', event_pro='mas_1ds_cfit_T2_TE')
    curve_fit_T2F = widget_button(curve_fit, accelerator='CTRL+SHIFT+3', value='T2 (FID)', event_pro='mas_1ds_cfit_t2_fid')

    time_Series = widget_button(analyze, value='Stability', /menu)
    time_Series_freq = widget_button(time_Series, value='Frequency', event_pro='mas_stability_plot')
    time_Series_phase = widget_button(time_Series, value='Phase', event_pro='mas_stability_plot')
    time_Series_amp = widget_button(time_Series, value='Amplitude', event_pro='mas_stability_plot')
    time_Series_amp = widget_button(time_Series, value='Linewidth', event_pro='mas_stability_plot')
    time_Series_amp = widget_button(time_Series, value='Export', event_pro='mas_stability_plot')
    
    calctemp = widget_button(analyze, value='Temperature', /menu)
    calctemp_MeOH = widget_button(calctemp, value='Pure Methanol', event_pro='mas_calctemp')
    calctemp_Ethyl_Glycol = widget_button(calctemp, value='Pure Ethylene Glycol', event_pro='mas_calctemp')
                                      
    freq_base = widget_base(left_panel, /frame, /column)
    lbl_freq_base = widget_label(freq_base, value="Spectrum")  
    freq_choose = widget_base(freq_base, /column, /nonexclusive)
    btn_show_freq_real = widget_button(view_frq, value="Hide Real", uvalue=1)
    btn_show_freq_imag = widget_button(view_frq, value="Show Imaginary", uvalue=0)
    btn_show_freq_abs  = widget_button(view_frq, value="Show Magnitude", uvalue=0)
    btn_show_data_pts = widget_button(freq_choose, value="Show Data Points")
    btn_show_zero_line = widget_button(freq_choose, value="Show Zero Line")
    widget_control, btn_show_zero_line, /set_button
    btn_show_trend_line = widget_button(freq_choose, value="Show Trend Line")
    
    ;freq_phase_corr = widget_base(freq_base, /column, ypad=0)
    ;void = widget_base(freq_choose, /row, /nonexclusive, ypad=0)
    btn_enable_phase_corr = widget_button(freq_choose, value="Phase Correction", $
                                          event_pro='mas_1ds_act_phase_corr_event')
    
    widget_control, btn_show_freq_real, /set_button
    
    freq_scale = widget_base(freq_base, /column)
    xzoom_base = widget_base(freq_scale, /row)                     
    void = widget_button(xzoom_base, value="X Zoom In", event_pro='mas_1ds_freq_xscale', uname='btn_xzoom_in')
    void = widget_button(xzoom_base, value="X Zoom Out", event_pro='mas_1ds_freq_xscale', uname='btn_xzoom_reset')
    units_base = widget_base(freq_scale, /row)
    lbl = widget_label(units_base, value='Axis Units:')
    aaa = widget_base(units_base, /row, /exclusive)
    btn_ppm = widget_button(aaa, value='ppm', uname='btn_ppm', event_pro='mas_1ds_axis_units_event')
    btn_hz  = widget_button(aaa, value='Hz', uname='btn_hz', event_pro='mas_1ds_axis_units_event')
    widget_control, btn_ppm, /set_button
    
    void = widget_slider(freq_scale, title="Plot Y Scaling (2^n)", min=-10, max=10, /drag,$
                         value=0, event_pro='mas_1ds_freq_yscale', uname='fq_yscale')
    ;void = widget_slider(freq_scale, title="Plot Y Offset", min=-10000, max=10000, /drag,$
    ;                     value=0, event_pro='mas_1ds_freq_yoffset', uname='fq_yoffset')
                            
    if (project.imndarray[ci].adim gt 1) then begin
        void = widget_slider(freq_scale, title="Array Select", min=0, max=adim-1, /drag,$
                             value=0, event_pro='mas_1ds_change_fid_event', uname='fid_select')
    endif else begin
        ;void = widget_slider(freq_scale, title="Array Select", min=0, max=1, sensitive=0, /drag,$
        ;                     value=0, event_pro='mas_1ds_change_fid_event', uname='fid_select')
    endelse
    
    weight_base = widget_base(freq_base, /frame, /column)
    lbl = widget_label(weight_base, value="Weight Function")
    btn_base = widget_base(weight_base, /nonexclusive)
    btn = widget_button(btn_base, value="Apply Weighting", uname='weight_activate', event_pro='mas_1ds_apply_weight_function')
    weight_type_base = widget_base(weight_base, /row, /exclusive)
    btn = widget_button(weight_type_base, value='Exponential', uname='weight_exp', event_pro='mas_1ds_apply_weight_function', sensitive=0)
    widget_control, btn, /set_button
    btn = widget_button(weight_type_base, value='Gaussian', uname='weight_gauss', event_pro='mas_1ds_apply_weight_function', sensitive=0)
    lbl_base = widget_base(weight_base, /row)
    sl = cw_fslider(lbl_base, title='Linewidth (Hz)', /drag, /edit, $
                    uname='weight_slider', minimum=0, maximum=200, value=0)
    
   ;; lbl = widget_label(lbl_base, value="Linewidth:")
    ;;txt = widget_text(lbl_base, xsize=9, uname='weight_param', value='-0.003', /editable)
    ;;lbl = widget_label(lbl_base, value="Hz")
   ;; btn = widget_button(weight_base, value="Apply Weight", event_pro='mas_1ds_apply_weight_function')
    
    ts_base = widget_base(left_panel, /frame, /column)
    lbl_ts_title = widget_label(ts_base, value="Free Induction Decay")
    ts_choose = widget_base(ts_base, /column, /nonexclusive)
    btn_show_ts_real = widget_button(view_ts, value="Hide Real", uvalue=1)
    btn_show_ts_imag = widget_button(view_ts, value="Show Imaginary", uvalue=0)
    btn_show_ts_abs  = widget_button(view_ts, value="Show Magnitude", uvalue=0)
    widget_control, btn_show_ts_real, /set_button

    ts_shift_base = widget_base(ts_base, /row)
    sl_ts_shift = widget_slider(ts_shift_base, title='Time Shift', $
                                minimum=-acq_size/4, maximum=acq_size/4, value=0, scroll=1, /drag, $
                                event_pro='mas_1ds_display_time_shift_event')
    lbl = widget_label(ts_shift_base, value='(Max @ '+string(max_ind, format='(G0)')+')')

    ts_shift_base = widget_base(ts_base, /row)
    sl_baseline = widget_slider(ts_shift_base, title='Baseline Correction', $
                                minimum=0, maximum=acq_size/4, value=0, scroll=1, /drag, $
                                event_pro='mas_1ds_display_baseline_event')
                                 
;    void = widget_button(xzoom_base, value="X Zoom In", event_pro='mas_1ds_ts_xscale', uname='btn_xzoom_in')
;    void = widget_button(xzoom_base, value="X Zoom Out", event_pro='mas_1ds_ts_xscale', uname='btn_xzoom_reset')

    freq_plot_base = widget_base(main_panel, /column, /frame)
    frq_combined = widget_draw(freq_plot_base, scr_xsize=850, scr_ysize=375, $
                               graphics_level=2, retain=2, $
                               /button_events, /motion_events, $
                               event_pro='mas_1ds_display_freq_button_event')
    vdim = [850,375]
    status_base = widget_base(main_panel, /row, /align_center)
    b1 = widget_base(status_base, /row, /align_left, /frame)
    void = widget_label(b1, value="Red Cursor [X,Real,Imag]")
    txt = widget_text(b1, xsize=8, uname='MARKER_01_XPOS', /all_events, event_pro='mas_1ds_keyboard_adjust_event')
    txt = widget_text(b1, xsize=8, uname='MARKER_01_REAL')
    txt = widget_text(b1, xsize=8, uname='MARKER_01_IMAG')
    b1 = widget_base(status_base, /row, /align_right, /frame)
    void = widget_label(b1, value="Green Cursor [X,Real,Imag]")
    txt = widget_text(b1, xsize=8, uname='MARKER_02_XPOS', /all_events, event_pro='mas_1ds_keyboard_adjust_event')
    txt = widget_text(b1, xsize=8, uname='MARKER_02_REAL')    
    txt = widget_text(b1, xsize=8, uname='MARKER_02_IMAG')    
    btn = widget_button(status_base, value="More Info...", event_pro='mas_1ds_stats_button_event')
    ;;btn = widget_button(status_base, value="Jump to Peak", event_pro='mas_1ds_stats_button_event')
    
    ts_plot_base = widget_base(main_panel, /column, /frame)
    ts_combined = widget_draw(ts_plot_base, scr_xsize=850, scr_ysize=250, $
                               graphics_level=2, retain=2, $
                               /button_events, /motion_events, $
                               event_pro='mas_1ds_display_ts_button_event')

    
    widget_control, tlb, /realize
    
    widget_control, frq_combined, get_value=freq_combined_draw
    widget_control, ts_combined, get_value=ts_combined_draw
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ts_plot = obj_new('mas_complex_plot', data, fq_window=ts_combined_draw, $
                      xlim=at)
    ts_plot->setHideFreqPlot, real=1, imag=1
    ts_plot->setCursorStyle, 1, /hide
    ts_plot->setCursorStyle, 2, /hide
    ts_plot->setHideZeroLine, 0
    ts_plot->setHideIntervalTrendLine, 1
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    fdata = mas_1ds_fft(data)
    
    ;;xlim = sw
    xlim = sw/(sf * 1e0)
;    max_abs = max(abs(fdata), max_abs_index)
    xzero =  4.75/xlim*acq_size + acq_size/2
    
    fq_plot = obj_new('mas_complex_plot', fdata, fq_window=freq_combined_draw, $
                      xlim=xlim, /reverse_xaxis)
    fq_plot->setHideFreqPlot, abs=1, imag=1
    fq_plot->setCursorStyle, 1, hide=0, thick=2, linestyle=2
    fq_plot->setCursorStyle, 2, hide=0, thick=2, linestyle=5
    fq_plot->setHideZeroLine, 0
    fq_plot->setHideIntervalTrendLine, 1
    fq_plot->setZeroReference, xzero
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    
    state = { project_id: ci, $
              num_fids: adim, $
              curr_fid: adim_start, $
              smooth_factor: 0, $
              acq_size: acq_size, $
              ts_data: data, $
              fq_data: fdata, $
              fq_baseline: dblarr(size(fdata, /dimensions)), $
              fq_baseline_correction_on: 0, $
              fq_button_pressed: 0, $
              ts_button_pressed: 0, $
              tlb: tlb, $
              phase_corr_tlb: 0, $
              base_stats_window: -1L, $
              btn_show_freq_real: btn_show_freq_real, $
              btn_show_freq_imag: btn_show_freq_imag, $
              btn_show_freq_abs: btn_show_freq_abs, $
              btn_show_zero_line: btn_show_zero_line, $
              btn_show_trend_line: btn_show_trend_line, $
              btn_show_ts_real: btn_show_ts_real, $
              btn_show_ts_imag: btn_show_ts_imag, $
              btn_show_ts_abs: btn_show_ts_abs, $
              btn_show_data_pts: btn_show_data_pts, $
              btn_enable_phase_corr: btn_enable_phase_corr, $
              fq_xaxis_units: 'ppm', $
              fq_marker_01_xfreqpos: sw/4, $
              fq_marker_01_xpixelpos: (sw/4)/sw*vdim[0], $ 
              fq_marker_02_xfreqpos: sw*0.75, $
              fq_marker_02_xpixelpos: (sw*0.75)/sw*vdim[0], $ 
              fq_marker_moving: 0L, $
              fq_phase_corr_active: 0, $
              fq_phase_corr_slope: 0.0, $
              fq_phase_corr_intercept: 0.0, $
              fq_phase_pivot: project.imndarray[ci].spect_phase_corr[2], $
              fq_0th_corr: float(project.imndarray[ci].spect_phase_corr[0]), $
              fq_1st_corr: float(project.imndarray[ci].spect_phase_corr[1]), $
              fq_0th_corr_widval: float(project.imndarray[ci].spect_phase_corr[0]), $
              fq_1st_corr_widval: float(project.imndarray[ci].spect_phase_corr[1]), $
              fq_phase_corr: 0, $
              fq_scale: [1, 1], $
              fq_xzoom_limits: [0, acq_size-1], $
              ts_scale: [1.15, 1.15], $
              ts_shift: 0L, $
              baseline_corr: 0L, $
              weight_type: '', $
              weight_param: 0.0, $
              weight_active: 0, $
              zero_ref: xzero, $
              ofq_plot: fq_plot, $
              ots_plot: ts_plot, $
              spect_width: sw, $
              base_freq: sf, $
              acq_time: at }
    state = ptr_new(state, /no_copy)
        
    widget_control, tlb, set_uvalue=state
    
    xmanager, "mas_display_1d_spect", tlb, /no_block, cleanup='mas_display_1d_spect_cleanup'

    phase_corr = mas_1ds_compute_phase_corr(state)
    (*state).fq_phase_corr_slope = phase_corr[0]
    (*state).fq_phase_corr_intercept = phase_corr[1]
    
    mas_1ds_apply_phase_corr, state, phase_corr, /global
    
;    mas_1ds_update_marker_values, state
    
end

pro mas_display_1d_spect_cleanup, tlb

    widget_control, tlb, get_uvalue=state
    
    obj_destroy, (*state).ofq_plot
    obj_destroy, (*state).ots_plot
    
    ptr_free, state
    
end

pro mas_1ds_freq_xscale, event

    widget_control, event.top, get_uvalue=state
    uname = widget_info(event.id, /uname)
    
    if (uname eq 'btn_xzoom_reset') then begin
        (*state).ofq_plot->applyXZoom, /reset
        (*state).fq_xzoom_limits = [0, (*state).acq_size-1]
        return
    endif

    (*state).ofq_plot->getCursorInfo, cursor01_info=c1, cursor02_info=c2
    
    left = min([c1.freq_xval, c2.freq_xval], max=right)
    (*state).ofq_plot->applyXZoom, left-100, right+100
    (*state).fq_xzoom_limits = [min([c1.data_index, c2.data_index]), $
                                max([c1.data_index, c2.data_index])]
    
end

pro mas_1ds_ts_xscale, event

    widget_control, event.top, get_uvalue=state
    uname = widget_info(event.id, /uname)
    
    if (uname eq 'btn_xzoom_reset') then begin
        (*state).ots_plot->applyXZoom, /reset
        return
    endif

    (*state).ots_plot->getCursorInfo, cursor01_info=c1, cursor02_info=c2
    
    left = min([c1.freq_xval, c2.freq_xval], max=right)
    (*state).ots_plot->applyXZoom, left-100, right+100
    
end

pro mas_1ds_freq_yoffset, event

    val = event.value
    
    widget_control, event.top, get_uvalue=state
    
    (*state).ofq_plot->applyYoffset, val

end

pro mas_1ds_freq_yscale, event

    val = event.value
    
    widget_control, event.top, get_uvalue=state
    
    (*state).ofq_plot->applyYZoom, val

end

;; See: Chen, et. al., An efficient algorithm for automatic
;;      phase correction of NMR spectra based on entropy
;;      minimization, Journ. Magnetic Resonance 158 (2002) 164-168

function mas_1ds_autophase_objective, params, feedback=feedback

    common autophase, data

    II = dcomplex(0,1)

    ;;;mas_1ds_apply_phase_corr, data.state, reverse(params)
    if (n_elements(data) eq 0) then return, -1
    
    phi = params[0] + params[1]*findgen(n_elements(data.spectra))/(n_elements(data.spectra)-1)
    
    S = exp(II*phi)*data.spectra
    R = real_part(S)
    
    R_prime = R;[1:*] - R[0:*]
    ;;R_prime = R_prime[1:*] - R_prime[0:*]
    
    h = abs(R_prime)/total(abs(R_prime))
    P = total( float(R lt 0)*R^2 ) * data.gam

    ent = -(total(h * alog(h))) + P
    
    if (keyword_set(feedback)) then begin
        print, "Penalty  : "+string(P, format='(G0)')
        print, "Entropy  : "+string(ent, format='(G0)')
        print, "Int. Imag: "+string(total(imaginary(S)), format="(G0)")
        print, "Int. Real: "+string(total(R), format="(G0)")
        print, "-----------------------------"
    endif
    
;    if (data.objective_niter eq n_elements(data.objective_history)) then begin
;        data.objective_niter = 0
;    endif
;    data.objective_history[data.objective_niter++] = ent
;    plot, data.objective_history[0:data.objective_niter-1]

    return, ent

end

pro mas_1ds_autophase01, event

    common autophase, data
    
    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        help, /traceback
        print, !error_state.MSG
        print, !error_state.CODE
        if (!error_state.CODE eq -637) then begin
            void = dialog_message("Phase correction failed to converge.", /error, /center)
        endif else begin
            void = dialog_message(["An error occurred during automatic phasing:", $
                                   !error_state.MSG], /error, /center)
        endelse
        return
    endif
    
    widget_control, event.top, get_uvalue=state
    
    active_fid = (*state).curr_fid
    data = { spectra: reform((*state).fq_data[*,active_fid]), $
             gam: 5e-5, $
             m: 1.0D, $
             state:state, $
             objective_history: dblarr(1), $
             objective_niter: 0L }
            
    start = [ 0.0D, 0.0D ]
    best = start
    bestmin = 1e7
    for j = 1, 6 do begin
    
        Xi = diag_matrix([1.0D,1.0D])
        POWELL, start, Xi, 1.0e-5, fmin, 'mas_1ds_autophase_objective', itmax=300, /double
        
        if (fmin lt bestmin) then begin
            best = start
            bestmin = fmin
        endif
        
        start[0] = 360/float(j) * !DTOR
        start[1] = 0
        
    endfor

    start = best

    max = max(abs((*state).fq_data), max_ind)
    pivot = float(max_ind)/(*state).acq_size
    
    m = start[1]
    b = start[0]
    
    y1 = m*pivot + b

    (*state).fq_0th_corr_widval = (y1*1.0/!DTOR) mod 360
    (*state).fq_1st_corr_widval =  m*1.0/!DTOR
    (*state).fq_phase_corr_slope = start[1]
    (*state).fq_phase_corr_intercept = start[0]
    (*state).fq_phase_pivot = max_ind
    
    mas_1ds_apply_phase_corr, state, reverse(start), /global
    mas_1ds_update_master_phase_corr, state, reverse(start)
    
    widget_control, widget_info(event.top, find_by_uname='fq_0th_order'), set_value=(*state).fq_0th_corr_widval
    widget_control, widget_info(event.top, find_by_uname='fq_1st_order'), set_value=(*state).fq_1st_corr_widval
    
end

pro mas_1ds_autophase0, event

    widget_control, event.top, get_uvalue=state

    fdata = (*state).fq_data
    
    props = [90-45.0, 90, 90+45.]
    areas = [0.0,0,0]
    dx = 5.0
    tol = 0.01
    max_iter = 40
    sl_event = { top: (*state).tlb, $
                 value: 1.0 }
    global_max = 0.0
    for i = 0, 360-1,5 do begin
        sl_event.value = i
        mas_1ds_0phase_corr, sl_event
        tmpmax = max(real_part((*state).fq_data))
        if (tmpmax gt global_max) then begin
            global_max = tmpmax
            props[*] = float(i)
        endif
    endfor
    
    props = [props[1]-dx, props[1], props[1]+dx]
    iter = 0
    while (iter++ lt max_iter) do begin
    
        for i = 0, n_elements(props)-1 do begin
            
            sl_event.value = props[i]
            mas_1ds_0phase_corr, sl_event
            areas[i] = abs(total(imaginary((*state).fq_data)))
           
        endfor
        
        min_area = min(areas, min_index)
        if (min_index eq 1) then begin
            dx /= 2.0
            if (dx lt tol) then break
            props[0] += dx
            props[2] -= dx
            ;; trough
        endif else if (min_index eq 0) then begin
            dx /= 2.0
            props[2] = props[0]
            props[1] = props[2]-dx
            props[0] = props[1]-dx
        endif else if (min_index eq 2) then begin
            dx /= 2.0
            props[0] = props[2]
            props[1] = props[0]+dx
            props[2] = props[1]+dx
        endif
    
    end
    
    if (total(real_part((*state).fq_data)) lt 0) then begin
        sl_event.value = (180.0 + sl_event.value) mod 360.0
        mas_1ds_0phase_corr, event
    endif
    
    phase_0_widget = widget_info(event.top, find_by_uname='fq_0th_order')
    if (widget_info(phase_0_widget, /valid_id)) then begin
        widget_control, phase_0_widget, set_value=sl_event.value
    endif
    
    if (iter ge max_iter) then begin
        print, "Convergence not reached in "+string(max_iter, format='(G0)')+" iterations."
    endif

    print, "Min imaginary area:", min_area

end

pro mas_1ds_0phase_corr, event

    common scan_data
    
    widget_control, event.top, get_uvalue=state
    
    val = float(event.value)
    (*state).fq_0th_corr_widval = val
    ;theta = val*!DTOR
    
    ;(*state).fq_0th_corr = theta
    
    phase_corr = mas_1ds_compute_phase_corr(state)
    (*state).fq_phase_corr_slope = phase_corr[0]
    (*state).fq_phase_corr_intercept = phase_corr[1]
    
    ;;void = mas_1ds_autophase_objective(reverse(phase_corr), /feedback)
    
    mas_1ds_apply_phase_corr, state, phase_corr
    mas_1ds_update_master_phase_corr, state, phase_corr
    
end

pro mas_1ds_1phase_corr, event

    widget_control, event.top, get_uvalue=state
    
    val = float(event.value)
    (*state).fq_1st_corr_widval = val
    ;theta = val*!DTOR

    ;(*state).fq_1st_corr = theta
    
    phase_corr = mas_1ds_compute_phase_corr(state)
    (*state).fq_phase_corr_slope = phase_corr[0]
    (*state).fq_phase_corr_intercept = phase_corr[1]

    ;;void = mas_1ds_autophase_objective(reverse(phase_corr), /feedback)
    
    mas_1ds_apply_phase_corr, state, phase_corr
    mas_1ds_update_master_phase_corr, state, phase_corr
    
end

pro mas_1ds_update_master_phase_corr, state, phase_corr

    common scan_data, project
    
    if (n_elements(project) ne 0) then begin
        
        phase_corr_n = [ (*state).fq_0th_corr_widval, $
                        (*state).fq_1st_corr_widval, $
                        (*state).fq_phase_pivot ]
                        
        project.imndarray[(*state).project_id].spect_phase_corr = phase_corr_n
        
    endif

end

pro mas_1ds_apply_weight_function, event;, function_name, params

    common scan_data, project
    
    widget_control, event.top, get_uvalue=state
    uname = widget_info(event.id, /uname)
    if (uname eq 'weight_activate') then begin
        widget_control, widget_info(event.top, find_by_uname='weight_exp'), sensitive=event.select
        widget_control, widget_info(event.top, find_by_uname='weight_gauss'), sensitive=event.select
        (*state).weight_active = event.select
    endif
        
    type_exp   = widget_info(widget_info(event.top, find_by_uname='weight_exp'), /button_set)
    type_gauss = widget_info(widget_info(event.top, find_by_uname='weight_gauss'), /button_set)
    
    if (type_exp) then begin
        (*state).weight_type = 'EM'
    endif else if (type_gauss) then begin
        (*state).weight_type = 'GM'
    endif else begin
        (*state).weight_type = ''
        return
    endelse
    
    widget_control, widget_info(event.top, find_by_uname='weight_slider'), get_value=val
    (*state).weight_param = float(val[0])/1000.0
    
    (*state).ts_data = *project.dataarray[(*state).project_id].state1
    mas_1ds_run_processing_chain, state, (*state).ts_data
    
end

pro mas_1ds_apply_baseline_event, event

    widget_control, event.top, get_uvalue=state
    
    old_baseline = (*state).fq_baseline
    (*state).fq_data += old_baseline
    (*state).fq_baseline[*] = 0
    
    checkcancel=0
    
    if ((*state).fq_baseline_correction_on eq 0) then begin
    
        mas_1ds_compute_fq_baseline, state, checkcancel=checkcancel
    
        if (checkcancel eq 0) then begin
            (*state).fq_data -= (*state).fq_baseline
            (*state).fq_baseline_correction_on = 1
            widget_control, event.id, set_value='Turn OFF Frequency Baseline Correction'
        endif else begin
            (*state).fq_baseline_correction_on = 0
            widget_control, event.id, set_value='Turn ON Frequency Baseline Correction'
            ;;widget_control, event.id, set_button=0
        endelse
        ;;window, /free, xsize=400, ysize=300
        ;;plot, (*state).fq_baseline
        
    endif else begin
    
        (*state).fq_baseline_correction_on = 0
        widget_control, event.id, set_value='Turn ON Frequency Baseline Correction'
        
    endelse

    (*state).ofq_plot->setData, (*state).fq_data
    
end


pro mas_1ds_apply_phase_corr, state, phase_corr, global=global

    common scan_data, project
    
    nfids = (*state).num_fids
    acq_size = (*state).acq_size
    
    if (keyword_set(global)) then begin
        m = (*state).fq_phase_corr_slope
        b = (*state).fq_phase_corr_intercept
    endif else begin
        m = phase_corr[0]
        b = phase_corr[1]
    endelse

    if (m ne 0) then begin
        freqs = findgen(acq_size)/(acq_size)
    endif else begin
        freqs = replicate(1.0, acq_size)
    endelse
    
    if (keyword_set(global)) then begin

        data_size = size((*state).ts_data)
        
        (*state).fq_data = mas_1ds_fft((*state).ts_data)

        freqs = rebin(freqs, n_elements((*state).fq_data)/nfids, nfids, /sample)
        correction = (complex(0,1.0))*(m*freqs + b)
        (*state).fq_data *= exp(correction)
        (*state).ofq_plot->setData, (*state).fq_data
        
    endif else begin
    
        correction = (complex(0,1.0))*(m*freqs + b)
        fdata = mas_1ds_fft((*state).ts_data[*,(*state).curr_fid])        
        fdata *= exp(correction)
        (*state).ofq_plot->setDisplayedData, fdata
        
    endelse
        
    mas_1ds_update_marker_values, state
    if (widget_info((*state).base_stats_window, /valid_id)) then begin
        mas_1ds_stats_button_event, { top: (*state).tlb }
    endif
    
    if (widget_info((*state).btn_show_trend_line, /button_set)) then begin
        (*state).ofq_plot->setCursorPosition, 0, -1
    endif
    
end

pro mas_1ds_swap_data, state, new_fid=new_fid

    common scan_data, project
    
    if (n_elements(project) eq 0) then return
    
    (*state).ots_plot->selectDataArrayElement, new_fid
    (*state).ofq_plot->selectDataArrayElement, new_fid
    (*state).curr_fid = new_fid
    
    mas_1ds_update_marker_values, state
    if (widget_info((*state).base_stats_window, /valid_id)) then begin
        mas_1ds_stats_button_event, { top: (*state).tlb }
    endif
    
end

pro mas_1ds_set_pivot, event

    widget_control, event.top, get_uvalue=state

    (*state).ofq_plot->getCursorInfo, cursor01_info=c1

    (*state).fq_phase_pivot = c1.data_index
    
    (*state).ofq_plot->removeMarker, 'PHASE_PIVOT'
    (*state).ofq_plot->addMarker, c1.data_index, 'PHASE_PIVOT'
    
    mas_1ds_update_master_phase_corr, state
    mas_1ds_apply_phase_corr, state, mas_1ds_compute_phase_corr(state)
    
    sl_1st_order = widget_info(event.top, find_by_uname='fq_1st_order')
    if (widget_info(sl_1st_order, /valid_id)) then begin
        widget_control, sl_1st_order, sensitive=1
    endif
    (*state).ofq_plot->setHeadsupMessageText, ['Phase Correction Active.', $
                                               'Pivot @ ' + string(c1.freq_xval, format='(G0)') + ' '+(*state).fq_xaxis_units]
    
end

pro mas_1ds_set_ref_value_event, event

    widget_control, event.top, get_uvalue=state
    widget_control, event.id, get_uvalue=txt_ref_value
    
    widget_control, txt_ref_value, get_value=ref_value
    ref_value = float(ref_value)
       
    (*state).ofq_plot->getCursorInfo, cursor01_info=c1
    
    xlim = (*state).spect_width/((*state).base_freq * 1e0)    
    xzero = c1.data_index + ref_value/xlim*(*state).acq_size
    if ((*state).fq_xaxis_units ne 'pts') then begin
        (*state).ofq_plot->setZeroReference, xzero;;c1.data_index
    endif
    (*state).zero_ref = xzero;;c1.data_index
    
    mas_1ds_update_marker_values, state
    
end

pro mas_1ds_set_zero, event

    widget_control, event.top, get_uvalue=state
    
    (*state).ofq_plot->getCursorInfo, cursor01_info=c1
    (*state).ofq_plot->setZeroReference, c1.data_index
    (*state).zero_ref = c1.data_index
    mas_1ds_update_marker_values, state
    if (widget_info((*state).base_stats_window, /valid_id)) then begin
        mas_1ds_stats_button_event, { top: (*state).tlb }
    endif
    
end

pro mas_1ds_reset_spectrum, event

    common scan_data
    
    widget_control, event.top, get_uvalue=state
    
    (*state).ts_data = *project.dataarray[(*state).project_id].state1
    data_size = size((*state).ts_data)
    
    if ((*state).ts_shift ne 0) then begin
        if (data_size[0] gt 1) then begin
            (*state).ts_data = shift((*state).ts_data, (*state).ts_shift, 0)
            if ((*state).ts_shift lt 0) then begin
                (*state).ts_data[data_size[1]-1+(*state).ts_shift:*,*] = complex(0,0)
            endif else begin
                (*state).ts_data[0:(*state).ts_shift-1,*] = complex(0,0)
            endelse
        endif else begin
            (*state).ts_data = shift((*state).ts_data, (*state).ts_shift)
            if ((*state).ts_shift lt 0) then begin
                (*state).ts_data[data_size[1]-1+(*state).ts_shift:*] = complex(0,0)
            endif else begin
                (*state).ts_data[0:(*state).ts_shift-1] = complex(0,0)
            endelse
        endelse
    endif
    
    fdata = mas_1ds_fft((*state).ts_data)    
    (*state).fq_data = fdata
 
    wid = widget_info((*state).tlb, find_by_uname='weight_activate')
    widget_control, wid, set_button=0
    
    wid = widget_info(event.top, find_by_uname='fq_0th_order')
    if (widget_info(wid, /valid_id)) then begin
        widget_control, wid, set_value=0
    endif
    (*state).fq_0th_corr = 0
    (*state).fq_phase_corr_slope = 0
    (*state).fq_0th_corr_widval = 0
    
    wid = widget_info(event.top, find_by_uname='fq_1st_order')
    if (widget_info(wid, /valid_id)) then begin
        widget_control, wid, set_value=0
    endif
    (*state).fq_1st_corr = 0
    (*state).fq_phase_corr_intercept = 0
    (*state).fq_1st_corr_widval = 0
    
    (*state).ofq_plot->setData, (*state).fq_data
    (*state).ots_plot->setData, (*state).ts_data
    
    mas_1ds_swap_data, state, new_fid=(*state).curr_fid
    
    mas_1ds_update_marker_values, state
    
end

pro mas_display_1d_spect_event, event

    widget_control, event.top, get_uvalue=state
    uname = widget_info(event.id, /uname)
    
    if (uname eq 'weight_slider') then begin
        mas_1ds_apply_weight_function, event
        return
    endif
    
    if (uname eq 'fq_0th_order') then begin
        mas_1ds_0phase_corr, event
        return
    endif
    
    print, mas_1ds_compute_fwhm(state)
    case event.id of 
        (*state).btn_show_freq_abs: begin
            widget_control, event.id, get_uvalue=ison
            (*state).ofq_plot->setHideFreqPlot, abs=ison
            ison = 1-ison
            new_value = (ison eq 1) ? 'Hide' : 'Show'
            new_value += ' Magnitude'
            widget_control, event.id, set_uvalue=ison, set_value=new_value
        end
        (*state).btn_show_freq_imag: begin
            widget_control, event.id, get_uvalue=ison
            (*state).ofq_plot->setHideFreqPlot, imag=ison
            ison = 1-ison
            new_value = (ison eq 1) ? 'Hide' : 'Show'
            new_value += ' Imaginary'
            widget_control, event.id, set_uvalue=ison, set_value=new_value
        end
        (*state).btn_show_freq_real: begin
            widget_control, event.id, get_uvalue=ison
            (*state).ofq_plot->setHideFreqPlot, real=ison
            ison = 1-ison
            new_value = (ison eq 1) ? 'Hide' : 'Show'
            new_value += ' Real'
            widget_control, event.id, set_uvalue=ison, set_value=new_value
        end
        (*state).btn_show_ts_abs: begin
            widget_control, event.id, get_uvalue=ison
            (*state).ots_plot->setHideFreqPlot, abs=ison
            ison = 1-ison
            new_value = (ison eq 1) ? 'Hide' : 'Show'
            new_value += ' Magnitude'
            widget_control, event.id, set_uvalue=ison, set_value=new_value
        end
        (*state).btn_show_ts_imag: begin
            widget_control, event.id, get_uvalue=ison
            (*state).ots_plot->setHideFreqPlot, imag=ison
            ison = 1-ison
            new_value = (ison eq 1) ? 'Hide' : 'Show'
            new_value += ' Imaginary'
            widget_control, event.id, set_uvalue=ison, set_value=new_value
        end
        (*state).btn_show_ts_real: begin
            widget_control, event.id, get_uvalue=ison
            (*state).ots_plot->setHideFreqPlot, real=ison
            ison = 1-ison
            new_value = (ison eq 1) ? 'Hide' : 'Show'
            new_value += ' Real'
            widget_control, event.id, set_uvalue=ison, set_value=new_value
        end
        else:
    endcase
    
;    (*state).ofq_plot->setHideFreqPlot, $
;                abs=1-widget_info((*state).btn_show_freq_abs , /button_set), $
;                real=1-widget_info((*state).btn_show_freq_real , /button_set), $
;                imag=1-widget_info((*state).btn_show_freq_imag , /button_set)
;     
;    (*state).ots_plot->setHideFreqPlot, $
;                abs=1-widget_info((*state).btn_show_ts_abs , /button_set), $
;                real=1-widget_info((*state).btn_show_ts_real , /button_set), $
;                imag=1-widget_info((*state).btn_show_ts_imag , /button_set)

    (*state).ofq_plot->setHidePlotSymbols, $
                1-widget_info((*state).btn_show_data_pts , /button_set)
                
    (*state).ofq_plot->setHideZeroLine, $
                1-widget_info((*state).btn_show_zero_line , /button_set)

    (*state).ofq_plot->setHideIntervalTrendLine, $
                1-widget_info((*state).btn_show_trend_line , /button_set)
    
end

pro mas_1ds_display_freq_button_event, event

    widget_control, event.top, get_uvalue=state
    
    (*state).ofq_plot->windowClickEvent, event, updated=updated
    
    if (updated eq 1) then begin
        mas_1ds_update_marker_values, state
        if (widget_info((*state).base_stats_window, /valid_id)) then begin
            mas_1ds_stats_button_event, { top: (*state).tlb }
        endif
    endif
    
    return
     
end

pro mas_1ds_display_ts_button_event, event

;    widget_control, event.top, get_uvalue=state
    
;    (*state).ots_plot->windowClickEvent, event, updated=updated
    
;    if (updated eq 1) then begin
;        mas_1ds_update_marker_values, state
;        if (widget_info((*state).base_stats_window, /valid_id)) then begin
;            mas_1ds_stats_button_event, { top: (*state).tlb }
;        endif
;    endif
    
    return
     
end

pro mas_1ds_change_fid_event, event

    widget_control, event.top, get_uvalue=state
    
    new_Fid = event.value
    
    mas_1ds_swap_data, state, new_fid=new_fid

end

pro mas_1ds_axis_units_event, event

    widget_control, event.top, get_uvalue=state
    
    uname = widget_info(event.id, /uname)
    
    case uname of 
        'btn_ppm': begin
            if ((*state).fq_xaxis_units eq 'ppm') then return
            ppm = (*state).spect_width/(*state).base_freq
            (*state).ofq_plot->setXAxislimit, ppm
            (*state).ofq_plot->refreshFrequencyPlot
            (*state).fq_xaxis_units = 'ppm'
        end
        
        'btn_hz': begin
            if ((*state).fq_xaxis_units eq 'Hz') then return
            hz = (*state).spect_width
            (*state).ofq_plot->setXAxislimit, hz
            (*state).ofq_plot->refreshFrequencyPlot
            (*state).fq_xaxis_units = 'Hz'
        
        end
    endcase
    
end

;+
; :Description:
;    Performs basic T2 fitting to multi echo experiment.
;
; :Params:
;    event - event from xmanager
;
; :Author: wtriplett
;-
pro mas_1ds_cfit_t2_te, event

    common scan_data, project
    ci = project.ci
    
    widget_control, event.top, get_uvalue=state
    
    pTE = project.imndarray[ci].echo_time_ptr
   
    if (not ptr_valid(pTE)) then begin
       void = dialog_message(['No echo times found for this scan.'], /center, /error)
       return
    endif else begin
        X = *pTE * 1e-3      
    endelse
    
    (*state).ofq_plot->getCursorInfo, cursor01_info=c1, cursor02_info=c2
    min_x_ind = min([c1.data_index, c2.data_index], max=max_x_ind)
    Y = total(real_part((*state).fq_data[min_x_ind:max_x_ind, *]),1)
    
    params = mpfitfun('mas_1ds_t2_function', X, Y, err, [max(Y), mean(X)], weights=1d, yfit=yfit, /quiet)
    resid = Y - yfit
    chi2 = total(resid^2)
    
    device, decomposed=0
    loadct, 0, /silent
    old_p_color = !P.color
    old_p_background = !P.background
    old_p_thick = !P.thick
    !P.color = 0
    !P.background = 255
    
    plot, X, Y, title='T2 Fit (Echoes)', xtitle='Echo Time (s)', ytitle='Signal', psym=1
    
    ;oplot, X, Y, psym=1
    
    X_interp = findgen(1000)/1000 * max(X)
    oplot, X_interp, mas_1ds_t2_function(X_interp, params), linestyle=0
    
    xyouts, 0.75, 0.800, string(params[0], format='(%"M0 = %0.3f")'), /normal
    xyouts, 0.75, 0.830, string(params[1], format='(%"T2 = %0.3f s")'), /normal
    xyouts, 0.75, 0.860, string(chi2, format='(%"CHI2 = %0.3e")'), /normal

    !P.color = old_p_color
    !P.background = old_p_background
    !P.thick = old_p_thick
    
end

;+
; :Description:
;    Performs basic T2 fitting to multi echo experiment.
;
; :Params:
;    event - event from xmanager
;
; :Author: wtriplett
;-
pro mas_1ds_cfit_t2_fid, event

  common scan_data, project
  ci = project.ci

  widget_control, event.top, get_uvalue=state

  pTE = project.imndarray[ci].echo_time_ptr
  sw = project.imndarray[ci].spect_spectral_width
  np = project.imndarray[ci].spect_acq_size
  
  X = (1+findgen(np))/sw
  data = (*state).ts_data
  p = (*state).curr_fid
  Y = abs(data[*,p])
  params = mpfitfun('mas_1ds_t2_function', X, Y, err, [max(Y), mean(X)], weights=1d, yfit=yfit, /quiet)
  resid = Y - yfit
  chi2 = total(resid^2)

  device, decomposed=0
  loadct, 0, /silent
  old_p_color = !P.color
  old_p_background = !P.background
  old_p_thick = !P.thick
  !P.color = 0
  !P.background = 255
  plot, X, Y, title='T2 Fit (FID)', xtitle='Time (s)', ytitle='Signal', psym=1 

  ;oplot, X, Y, psym=1

  X_interp = findgen(1000)/1000 * max(X)
  oplot, X_interp, mas_1ds_t2_function(X_interp, params), linestyle=0

  xyouts, 0.75, 0.800, string(params[0], format='(%"M0 = %0.3f")'), /normal
  xyouts, 0.75, 0.830, string(params[1], format='(%"T2 = %0.3f s")'), /normal
  xyouts, 0.75, 0.860, string(chi2, format='(%"CHI2 = %0.3e")'), /normal

  !P.color = old_p_color
  !P.background = old_p_background
  !P.thick = old_p_thick

end

;+
; :Description:
;    Performs basic T1 fitting to inversion recovery experiment.
;
; :Params:
;    event - event from xmanager
;
; :Author: wtriplett
;-
pro mas_1ds_cfit_T1_TI, event

    common scan_data, project
    ci = project.ci
    
    widget_control, event.top, get_uvalue=state
    
    pTI = project.imndarray[ci].inversion_time_ptr
    
    if (not ptr_valid(pTI)) then begin
        void = dialog_message(['No inversion times found for this scan.'], /center, /error)
        return
    endif else begin
        X = *pTI * 1e-3
    endelse
    
    (*state).ofq_plot->getCursorInfo, cursor01_info=c1, cursor02_info=c2    
    min_x_ind = min([c1.data_index, c2.data_index], max=max_x_ind)
    
    Y = total(real_part((*state).fq_data[min_x_ind:max_x_ind, *]),1)

    params = mpfitfun('mas_1ds_inv_rec_function', X, Y, err, [max(Y), mean(X)], weights=1d, yfit=yfit, /quiet)
    resid = Y - yfit
    chi2  = total(resid^2)

    device, decomposed=0
    loadct, 0, /silent
    old_p_color = !P.color
    old_p_background = !P.background
    old_p_thick = !P.thick
    !P.color = 0
    !P.background = 255
    
    plot, X, Y, title='T1 fit by Inversion Recovery', xtitle='Inversion Time (s)', ytitle='Signal', psym=1 ;, linestyle=1
    ;oplot, X, Y, psym=1
    
    X_interp = findgen(1000)/1000 * max(X)
    oplot, X_interp, mas_1ds_inv_rec_function(X_interp, params), linestyle=0
    xyouts, 0.75, 0.20, string(params[0], format='(%"M0 = %0.3f")'), /normal
    xyouts, 0.75, 0.23, string(params[1], format='(%"T1 = %0.3f s")'), /normal
    xyouts, 0.75, 0.26, string(chi2, format='(%"CHI2 = %0.3e")'), /normal

    !P.color = old_p_color
    !P.background = old_p_background
    !P.thick = old_p_thick
    
end

pro mas_1ds_plot_spect_multi_window_event, event

    widget_control, event.top, get_uvalue=stash
    
    state = (*stash).state

    widget_control, (*stash).txt_start, get_value=temp
    start_spect = (long(temp))[0]
    widget_control, (*stash).txt_end, get_value=temp
    end_spect = (long(temp))[0] 
    widget_control, (*stash).txt_step, get_value=temp
    ind_step = (long(temp))[0] 
    
    widget_control, (*stash).sl_yposition, get_value=yposition
    widget_control, (*stash).sl_yscale, get_value=yscale
    widget_control, (*stash).sl_yoffset, get_value=yoffset
    
    widget_control, (*stash).sl_xposition, get_value=xposition
    widget_control, (*stash).sl_xscale, get_value=xscale
    widget_control, (*stash).sl_xoffset, get_value=xoffset
    color_scheme_id = widget_info((*stash).dl_color_scheme, /droplist_select)
    
    case event.id of 
    
    
        (*stash).btn_save_tif: begin
            tiffile = dialog_pickfile(/write, filter='*.tif;*.tiff;*.TIF;*.TIFF', $
                                      default_extension='tif', /overwrite_prompt)
            if (tiffile eq '') then return
            mas_1ds_plot_spect_multi, state, start_spect, end_spect, ind_step, $
                                      yscale=yscale, xscale=xscale, $
                                      xoffset=xoffset, yoffset=yoffset, $
                                      xposition=xposition, yposition=yposition, $
                                      tiffile=tiffile, color_scheme=color_scheme_id
            return
        end

        (*stash).btn_save_ps: begin
            psfile = dialog_pickfile(/write, filter='*.ps', $
                                      default_extension='ps', /overwrite_prompt)
            if (psfile eq '') then return
            mas_1ds_plot_spect_multi, state, start_spect, end_spect, ind_step, $
                                      yscale=yscale, xscale=xscale, $
                                      xoffset=xoffset, yoffset=yoffset, $
                                      xposition=xposition, yposition=yposition, $
                                      psfile=psfile, color_scheme=color_scheme_id
            return
        end
        
        (*stash).btn_cancel: begin
            ptr_free, stash
            widget_control, event.top, /destroy
            return
        end
        
        else: begin
        
        end
        
    endcase
    
    mas_1ds_plot_spect_multi, state, start_spect, end_spect, ind_step, $
                              yscale=yscale, xscale=xscale, $
                              xoffset=xoffset, yoffset=yoffset, $
                              xposition=xposition, yposition=yposition, $
                              color_scheme=color_scheme_id
    
end


pro mas_1ds_plot_spect_multi_window, event

    widget_control, event.top, get_uvalue=state

    wind_base = widget_base(group_leader=(*state).tlb, title='Plot Multiple Spectra', /column);, /modal)
    
    base = widget_base(wind_base, /column, /base_align_right, /align_center)
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    base_top = widget_base(base, column=2)
    base_top_l = widget_base(base_top, /column, /base_align_right, /align_center, /frame)
    base_top_r = widget_base(base_top, /column, /base_align_right, /align_center, /frame)
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    b1 = widget_base(base_top_l, /row)
    lbl = widget_label(b1, value='Starting Spectrum:')
    txt_start = widget_text(b1, xsize=5, value='0', /edit)
    
    b2 = widget_base(base_top_l, /row)
    lbl = widget_label(b2, value='Ending Spectrum:')
    txt_end = widget_text(b2, xsize=5, value=strtrim(string((*state).num_fids-1),2), /edit)

    b2 = widget_base(base_top_l, /row)
    lbl = widget_label(b2, value='Step Increment:')
    txt_step = widget_text(b2, xsize=5, value='1', /edit)
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    b2 = widget_base(wind_base, /row, xpad=0, /base_align_center, /align_center)
    lbl = widget_label(b2, value='Color Scheme:')
    loadct, get_names=ct_names
    ct_names = [ 'Solid Black', ct_names ]
    dl_color_scheme = widget_droplist(b2, value=ct_names)
    
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    bot_base = widget_base(wind_base, column=2)
    
    base = widget_base(bot_base, /column, /base_align_center, /align_center, /frame)
    
    b3p5 = widget_base(base, /row)
    sl_yposition = cw_fslider(b3p5, title='YPosition', minimum=0, maximum=1, value=0, /edit)
    
    b3p5 = widget_base(base, /row)
    sl_yscale = cw_fslider(b3p5, title='YScale',  minimum=0, maximum=1, value=0.5, /edit)

    b3p5 = widget_base(base, /row)
    sl_yoffset = cw_fslider(b3p5, title='YOffset',  minimum=-1, maximum=1, value=0, /edit)

    base = widget_base(bot_base, /column, /base_align_center, /align_center, /frame)
    
    b3p5 = widget_base(base, /row)
    sl_xposition = cw_fslider(b3p5, title='XPosition', minimum=-1, maximum=1, value=0, /edit)

    b3p5 = widget_base(base, /row)
    sl_xscale = cw_fslider(b3p5, title='XScale', minimum=0, maximum=5, value=1, /edit)

    b3p5 = widget_base(base, /row)
    sl_xoffset = cw_fslider(b3p5, title='XOffset', minimum=-1, maximum=1, value=0, /edit)
    
    b4 = widget_base(wind_base, /row, /base_align_center, /align_center)
    
    btn_cancel   = widget_button(b4, value='Close Window')
    btn_save_tif = widget_button(b4, value='Save as TIF')
    btn_save_ps  = widget_button(b4, value='Save as Postscript')
    
    stash = ptr_new({ state: state, $
                      start_spectrum: 0, $
                      end_spectrum:(*state).num_fids-1, $
                      sl_yscale: sl_yscale, $
                      sl_yposition: sl_yposition, $
                      sl_yoffset: sl_yoffset, $
                      sl_xscale: sl_xscale, $
                      sl_xposition: sl_xposition, $
                      sl_xoffset: sl_xoffset, $
                      dl_color_scheme: dl_color_scheme,$
                      txt_start: txt_start, $
                      txt_end: txt_end, $
                      txt_step: txt_step, $
                      btn_save_tif: btn_save_tif, $
                      btn_save_ps: btn_save_ps, $
                      btn_cancel: btn_cancel })
                      
    widget_control, wind_base, /realize, set_uvalue=stash
    
    xmanager, 'mas_1ds_plot_spect_multi_window', wind_base, /no_block
    
    
end


pro mas_1ds_plot_spect_multi, state, start_ind, end_ind, ind_step, $
                              yscale=yscale, yoffset=yoffset, yposition=yposition, $
                              xscale=xscale, xoffset=xoffset, xposition=xposition, $
                              color_scheme=color_scheme, $
                              xsize=xsize, ysize=ysize, tiffile=tiffile, psfile=psfile

    if (n_elements(color_scheme) eq 0) then begin
        color_scheme = 0 
    endif
    
    if (n_elements(start_ind) eq 0) then start_ind = 0
    if (n_elements(end_ind) eq 0) then end_ind = (*state).num_fids-1
    if (n_elements(ind_step) eq 0) then ind_step = 1
    if (n_elements(yscale) eq 0) then begin
        yscale = .9
    endif else begin
        yscale = ROUND((yscale*100))/100.0 ;; (ysqueeze > 0.0) < 1.0
    endelse
    if (n_elements(yoffset) eq 0) then yoffset = 0.0
    if (n_elements(xoffset) eq 0) then xoffset = 0.0
    if (n_elements(xscale) eq 0) then xscale = 0.5
    if (n_elements(xposition) eq 0) then xposition = 0.0
    yscale = (1.0-yscale)^3 > 0.001
    
    spect_start = ( min([start_ind,end_ind]) > 0 ) < ((*state).num_fids-1)
    spect_end   = ( max([start_ind,end_ind]) > 0 ) < ((*state).num_fids-1)
    
    spect_start = float(spect_start)
    spect_end   = float(spect_end)
    
    ywindow_size = 600
    num_fids = (spect_end-spect_start)+1
    acq_size = (*state).acq_size
    
    spect = (*state).fq_data
    
    xlim = (*state).spect_width/(*state).base_freq
    
    (*state).ofq_plot->getCursorInfo, $
                       cursor01_info=cur01, $
                       cursor02_info=cur02

    normalized_ts_indices = findgen(acq_size)/acq_size
    if ((*state).zero_ref ne 0) then begin
        normalized_fq_indices = (findgen(acq_size)-acq_size+(*state).zero_ref)/acq_size
    endif else begin
        normalized_fq_indices = (findgen(acq_size))/acq_size
    endelse
    
    xaxis_ppm = reverse(normalized_fq_indices) * xlim
    
    old_device  = !D.name
    
    if (keyword_set(tiffile)) then begin
        set_plot, 'Z'
        device, set_pixel_depth=24;;set_pixel_depth=8
        device, set_resolution=[11*150,8.5*150], decomposed=1
    endif else if (keyword_set(psfile)) then begin
        set_plot, 'PS', /interpolate, /copy
        device, /landscape, /color, BITS_PER_PIXEL = 8, FONT_SIZE=12, decomposed=1
        device, FILENAME = psfile, /inches;, xsize=11, ysize=8.5, /TIMES, /BOLD
    endif else begin
        device, window_state=ws, decomposed=1
        if (ws[23] eq 0) then begin
            window, 23, xsize=800, ysize=ywindow_size
        endif else wset, 23
        
    endelse
    
    spect = real_part(spect)
    spect /= max(spect)
    
    ymax = max(spect, min=ymin)
    yrange = yscale*abs(ymax - ymin)
    yrange_total = yrange*num_fids

    left  = (*state).fq_xzoom_limits[0]
    right = (*state).fq_xzoom_limits[1]
    
    t3d, /reset
    t3d, scale=[xscale,1,1]
    t3d, translate=[xposition,yposition,0]
    
    case color_scheme of
        0: begin
            loadct, 0
            rgb_table = transpose([ 0,0,0 ])
        end
        else: begin
            loadct, color_scheme-1, rgb_table=rgb_table
        end
    endcase 

    num_colors = n_elements(rgb_table)/3
    
    for p = spect_start, spect_end, ind_step do begin
    
        data = reform(spect[left:right,p])
        axis = xaxis_ppm[left:right]
        
        fidnum = p-spect_start

        colnum = ( fidnum/((spect_end-spect_start+1)) * num_colors ) mod num_colors
        color = rgb_table[colnum,*]
        color = color[0] + color[1]*256L + color[2]*(256L)^2
        
        if (fidnum eq 0) then begin
            title = 'Multiple Spectra'
            plot, /t3d, axis, data+fidnum*yrange-ymin, xticklayout=1, $
                  xstyle=1+8, ystyle=1+4, xtitle='ppm',$; title=title, $
                  yrange=[0, yrange*num_fids], $
                  xrange=[max(axis),min(axis)], $
                  background='FFFFFF'x, color='000000'x
            
        endif else begin
            t3d, translate=[xoffset/( (spect_end-spect_start)/float(ind_step) ), $
                            yoffset/( (spect_end-spect_start)/float(ind_step) ), 0]
            oplot, /t3d, axis-(fidnum*0), data+fidnum*yrange-ymin, color=color; '000000'x;color
        endelse
        
        empty
;        id_string = string(p, format='(%"#%03d")')
;        xyouts, /t3d, charsize=0.75, max(axis)+0.01*(abs(max(axis)-min(axis))), (fidnum*yrange-ymin)+yrange*0.25, id_string, /data, color='000000'x
         
    endfor
    t3d, /reset
    loadct, 0
    if (keyword_set(tiffile)) then begin
        data = reverse(tvrd(TRUE=3),2)
        write_tiff, tiffile, data, units=2, planarconfig=2
        set_plot, old_device
    endif else if (keyword_set(psfile)) then begin
        device, /close_file
        set_plot, old_device
    endif

end

pro mas_1ds_save_main_plots, event

    common scan_data, project
    
    old_device  = !D.name
    old_p_multi = !P.multi
    
    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        old_device = !D.name
        old_p_multi = !P.multi
        print, !error_state.MSG
        help, /traceback
        return
    endif
    
    widget_control, event.top, get_uvalue=state
    uname = widget_info(event.id, /uname)
    
    widget_control, widget_info(event.top, find_by_uname='fq_yscale'), get_value=val
    zoom_val = 1.0/float((2.0^val))
    
    (*state).ofq_plot->getCursorInfo, $
                       cursor01_info=cur01, $
                       cursor02_info=cur02
                       
    print, zoom_val
    if (uname eq 'file_save_tif') then begin
        filter='*.tif;*.tiff'
        def_ext = 'tif
    endif else begin
        filter='*.ps'
        def_ext = 'ps'
    endelse
    
    filename = dialog_pickfile(path=project.current_path, filter=filter, /write, $
                               /overwrite_prompt, default_extension=def_ext, /fix_filter, $
                               title='Where would you like to save the plot?')
                               
    if (filename eq '') then return

    if (uname eq 'file_save_tif') then begin
        set_plot, 'Z'
        device, set_pixel_depth=8
        device, set_resolution=[11*150,8.5*150], decomposed=1
    endif else begin
        set_plot, 'PS', /interpolate, /copy
        device, /landscape, /color, BITS_PER_PIXEL = 8, FONT_SIZE=12
        device, FILENAME = filename, /inches, xsiz=9, ysiz=7, /TIMES, /BOLD
    endelse
    
    !P.multi = [0,1,2]
    !p.multi = [0,0,0,0,0]
    
    acq_size    = (*state).acq_size
    acq_time    = (*state).acq_time
    spect_width = (*state).spect_width
    
    fq_yrange = [min(real_part((*state).fq_data)), $
                 max(real_part((*state).fq_data))]
    fq_yrange *= zoom_val
   
    fq_xrange = [min([cur01.freq_xval, cur02.freq_xval]), $
                 max([cur01.freq_xval, cur02.freq_xval])]
    
    if ((*state).fq_xaxis_units eq 'ppm') then begin
        xlim = (*state).spect_width/(*state).base_freq
    endif else begin
        xlim = (*state).spect_width
    endelse
           
    normalized_ts_indices = findgen(acq_size)/acq_size
    if ((*state).zero_ref ne 0) then begin
        normalized_fq_indices = (findgen(acq_size)-acq_size+(*state).zero_ref)/acq_size
    endif else begin
        normalized_fq_indices = (findgen(acq_size))/acq_size
    endelse
    
    normalized_fq_indices = reverse(normalized_fq_indices)
    
    left  = (*state).fq_xzoom_limits[0]
    right = (*state).fq_xzoom_limits[1]
    
    ;;fq_yrange = [min(real_part((*state).fq_data)),max(real_part((*state).fq_data))]*10
    
    for p = (*state).curr_fid, (*state).curr_fid do begin ;; (*state).num_fids-1 do begin
    ;;for p = 0, (*state).num_fids-1 do begin
        scan_name = project.imndarray[project.ci].display_name
        data = reverse(real_part(reform((*state).fq_data[left:right,p])))
        axis = normalized_fq_indices[left:right] * xlim
        plot,  axis, $
               reverse(data), xrange=[max(axis),min(axis)],$
               xstyle=1, xtitle=(*state).fq_xaxis_units, $
               title=string(scan_name,p, format='(%"%s!CFrequency Spectrum, FID %03d")'), $
               background='FFFFFF'x, color='000000'x, yrange=fq_yrange, thick=4.0, Xthick = 4.0, Ythick = 4.0, CHARTHICK = 3.0, xcharsize = 2, ycharsize = 2
        oplot, normalized_fq_indices * spect_width, replicate(0,acq_Size), linestyle=1, color='000000'x, thick=4.0

;        plot, normalized_ts_indices * acq_time, $
;              (*state).ts_data[*,p], xrange=[0,acq_time], xstyle=1, $
;              title=string(p, format='(%"Time Series, FID %03d")'), $
;              xtitle='Time (ms)', $
;              background='FFFFFF'x, color='000000'x, thick=2.0
;        oplot,[0, acq_time], [0,0], linestyle=1, color='000000'x, thick=2.0

        empty
        
    endfor
    
    if (uname eq 'file_save_tif') then begin
       data = reverse(tvrd(),2)
       write_tiff, filename, data, units=2;;, $
                   ;;xresol=(*state).save_resolution, $
                  ;; yresol=(*state).save_resolution    
    endif else begin
        device, /close_file       
    endelse
    
    set_plot, old_device
    !P.multi = old_p_multi

end


pro mas_1ds_display_baseline_event, event

    common scan_Data, project

    widget_control, event.top, get_uvalue=state

    (*state).ts_data = *project.dataarray[(*state).project_id].state1
    
    (*state).baseline_corr = event.value
    
    mas_1ds_run_processing_chain, state, (*state).ts_data
    
end


pro mas_1ds_display_time_shift_event, event

    common scan_data, project
    
    widget_control, event.top, get_uvalue=state
    
    if ((*state).ts_shift ne event.value) then begin
        ;(*state).fq_phase_corr_slope = 0
        ;(*state).fq_1st_corr = 0
        ;(*state).fq_1st_corr_widval = 0
        ;mas_1ds_update_master_phase_corr, state
    endif
        
    (*state).ts_shift = event.value

    (*state).ts_data = *project.dataarray[(*state).project_id].state1
    
    mas_1ds_run_processing_chain, state, (*state).ts_data
    
    ;;mas_1ds_reset_spectrum, event
    
end

pro mas_1ds_update_marker_values, state

    (*state).ofq_plot->getCursorInfo, cursor01_info=c1, cursor02_info=c2
    widget_control, widget_info((*state).tlb, find_by_uname='MARKER_01_XPOS'), set_value=string(c1.freq_xval, format='(G0.6)')
    widget_control, widget_info((*state).tlb, find_by_uname='MARKER_01_REAL'), set_value=string(real_part(c1.freq_yval), format='(G0.6)')
    widget_control, widget_info((*state).tlb, find_by_uname='MARKER_01_IMAG'), set_value=string(imaginary(c1.freq_yval), format='(G0.6)')
   
    widget_control, widget_info((*state).tlb, find_by_uname='MARKER_02_XPOS'), set_value=string(c2.freq_xval, format='(G0.6)')
    widget_control, widget_info((*state).tlb, find_by_uname='MARKER_02_REAL'), set_value=string(real_part(c2.freq_yval), format='(G0.6)')
    widget_control, widget_info((*state).tlb, find_by_uname='MARKER_02_IMAG'), set_value=string(imaginary(c2.freq_yval), format='(G0.6)')
    
end

pro mas_1ds_act_weight_func_event, event

    widget_control, event.top, get_uvalue=state
    
    weight_tlb = widget_base(title="Weight Function", group_leader=(*state).tlb, /column)
    
    b = widget_base(weight_tlb, /row)
    void = widget_label(b, value='Function: ')
    
    void = widget_droplist(b, value=['Exponential'])
   
    widget_control, weight_tlb, /realize
    

end


pro mas_1ds_act_phase_corr_event, event

    widget_control, event.top, get_uvalue=state
    
    sen = event.select
        
    if ((*state).num_fids gt 1) then begin
        widget_control, widget_info(event.top, find_by_uname='fid_select'), sensitive=1-sen
    endif
    
    if (sen eq 1) then begin
    
        (*state).ofq_plot->setCursorStyle, 1, hide=0, thick=2, linestyle=0, alpha=1.0
        (*state).ofq_plot->setCursorStyle, 2, hide=0, thick=2, alpha=0.3
        (*state).ofq_plot->setHideMarker, 'PHASE_PIVOT', 0
        
        phase_corr_tlb = widget_base(title="Phase Correction", group_leader=(*state).tlb, /column, $
                                     kill_notify='mas_1ds_phase_corr_destroy_event')
        void = widget_label(phase_corr_tlb, value="Phase Correction Controls")
;        sl_0th_order = widget_slider(phase_corr_tlb, title='0th Order', min=-360, max=360, /drag,$
 ;                                    value=(*state).fq_0th_corr_widval, event_pro='mas_1ds_0phase_corr',$
 ;                                    uname='fq_0th_order', sensitive=1)
        tb = widget_base(phase_corr_tlb, /base_align_center, /align_center, /column)
        sl_0th_order = cw_fslider(tb, title='0th Order', min=-360, max=360, scroll=1, /drag,$
                                  value=(*state).fq_0th_corr_widval,uname='fq_0th_order') 
        sl_1st_order = widget_slider(phase_corr_tlb, title='1st Order', min=-7200, max=7200, /drag,$
                                     value=(*state).fq_1st_corr_widval, event_pro='mas_1ds_1phase_corr', $
                                     uname='fq_1st_order', sensitive=0)
        btn_0_auto = widget_button(phase_corr_tlb, value="Automatic Phasing", event_pro='mas_1ds_autophase01' );;event_pro='mas_1ds_autophase0')

        if ((*state).fq_phase_pivot ne 0) then begin
            widget_control, sl_1st_order, sensitive=1
            (*state).ofq_plot->removeMarker, 'PHASE_PIVOT'
            (*state).ofq_plot->addMarker, (*state).fq_phase_pivot, 'PHASE_PIVOT'
        endif
        
        btn = widget_button(phase_corr_tlb, value="Set Pivot Point", uname='btn_set_pivot', $
                            event_pro='mas_1ds_set_pivot', sensitive=1)
        base_ref_value = widget_base(phase_corr_tlb, /row)
        lbl_ref_value = widget_label(base_ref_value, value="Ref. Value (ppm):")
        txt_ref_value = widget_text(base_ref_value, xsize=4, value='4.75', uname='txt_ref_value', /edit)
        btn_set_zero = widget_button(base_ref_value, value="Set", uname='btn_set_zero', $
                                     event_pro='mas_1ds_set_ref_value_event', sensitive=1, uvalue=txt_ref_value)
        
;        btn = widget_button(phase_corr_tlb, value="Set Zero Reference", uname='btn_set_zero', $
;                            event_pro='mas_1ds_set_zero', sensitive=1)
        btn = widget_button(phase_corr_tlb, value="Reset Spectrum", uname='btn_reset_spectrum', $
                            event_pro='mas_1ds_reset_spectrum')
        (*state).phase_corr_tlb = phase_corr_tlb
        widget_control, phase_corr_tlb, set_uvalue=state
        widget_control, phase_corr_tlb, /realize
        (*state).fq_phase_corr_active=sen        
        
        if ((*state).fq_phase_pivot ne 0) then begin
            (*state).ofq_plot->getCursorInfo, cursor01_info=c1
            (*state).ofq_plot->setHeadsupMessageText, ['Phase Correction Active.', $
                                                       'Pivot @ ' + string(c1.freq_xval, format='(G0)')+' '+(*state).fq_xaxis_units]
        endif else begin
            (*state).ofq_plot->setHeadsupMessageText, ['Phase Correction Active.',$
                                                       'Note: pivot must be set to enable 1st order correction.', $
                                                       'Position the solid red cursor at the desired pivot', $
                                                       'location and click "Set Pivot Point"']
        endelse
        
        xmanager, 'mas_display_1d_spect', phase_corr_tlb, /no_block
        
    endif else begin
        
        if (widget_info((*state).phase_corr_tlb, /valid_id)) then begin
            widget_control, (*state).phase_corr_tlb, /destroy
        endif
        
        (*state).fq_phase_corr_active=sen
        (*state).phase_corr_tlb = 0
        (*state).ofq_plot->setHeadsupMessageText, ''
        (*state).ofq_plot->setHideMarker, 'PHASE_PIVOT', 1
        
    endelse

    
end

pro mas_1ds_phase_corr_destroy_event, event

    widget_control, event, get_uvalue=state

    if (not ptr_valid(state)) then begin
        return
    endif
    
    (*state).ofq_plot->setCursorStyle, 1, hide=0, thick=2, linestyle=2, alpha=1.0
    (*state).ofq_plot->setCursorStyle, 2, hide=0, thick=2, alpha=1.0    
    if ((*state).fq_phase_corr_active) then begin
        mas_1ds_apply_phase_corr, state, mas_1ds_compute_phase_corr(state), /global
    endif    
    
    if ((*state).num_fids gt 1) then begin
        widget_control, widget_info((*state).tlb, find_by_uname='fid_select'), sensitive=1
    endif
    widget_control, (*state).btn_enable_phase_corr, set_button=0
    (*state).fq_phase_corr_active = 0
    (*state).phase_corr_tlb = 0
    (*state).ofq_plot->setHideMarker, 'PHASE_PIVOT', 1
    (*state).ofq_plot->setHeadsupMessageText, ''
    
end

pro mas_1ds_keyboard_adjust_event, event

    if (tag_names(event, /structure_name) ne 'WIDGET_TEXT_CH') then return
    
    if (event.ch eq 43B or event.ch eq 61B) then begin
        amt = 1; '+' or '='
    endif else if (event.ch eq 45) then begin
        amt = -1; '-'
    endif else return
    
    widget_control, event.top, get_uvalue=state
    uname = widget_info(event.id, /uname)
    
    (*state).ofq_plot->getCursorInfo, cursor01_info=c1, cursor02_info=c2
    
    if (uname eq 'MARKER_01_XPOS') then begin
        (*state).ofq_plot->setCursorPosition, c1.data_index + amt, 1
    endif else if (uname eq 'MARKER_02_XPOS') then begin
        (*state).ofq_plot->setCursorPosition, c2.data_index + amt, 2
    endif
    
    mas_1ds_update_marker_values, state

    return
    
    if (widget_info((*state).base_stats_window, /valid_id)) then begin
        mas_1ds_stats_button_event, { top: (*state).tlb }
    endif

end
pro mas_stability_plot, event

  common scan_data

  widget_control, event.top, get_uvalue=state
  widget_control, event.id, get_value = plot_type
  
  (*state).ofq_plot->getCursorInfo, cursor01_info=c1, cursor02_info=c2
  (*state).ofq_plot->getDataValues, c1.data_index, real=all_red_real, /all_fids
  (*state).ofq_plot->getDataValues, c2.data_index, real=all_grn_real, /all_fids
  
  if c1.data_index lt c2.data_index then ft_trunc = (*state).fq_data[c1.data_index:c2.data_index,*] else $
    ft_trunc = (*state).fq_data[c2.data_index:c1.data_index,*]
  N =  (*state).num_fids
  freq = make_Array(N,/float)
  amp = make_Array(N,/float)
  phi = make_Array(N,/float)
  fwhm = make_Array(N,/float)
  for i = 0,N-1 do begin
    if c1.data_index lt c2.data_index then ind = c1.data_index + where(real_part(ft_trunc[*,i]) eq max(real_part(ft_trunc[*,i]))) else $
      ind = c2.data_index + where(real_part(ft_trunc[*,i]) eq max(real_part(ft_trunc[*,i])))
    (*state).ofq_plot->getDataValues, ind, freqval = freqval, which_fid = i
    (*state).ofq_plot->getDataValues, ind, real = real, which_fid = i
    (*state).ofq_plot->getDataValues, ind, imag = imag, which_fid = i

    fwhm[i] = mas_1ds_compute_fwhm(state, index = ind, which_fid = i)
    freq[i] = freqval
    amp[i] = sqrt(real^2 + imag^2)
    phi[i] = 180*atan(complex(real,imag),/phase)/!Pi 
  endfor
  tex = project.imndArray[project.ci].spect_rep_time*project.imndArray[project.ci].spect_num_avg/(60.0*60.0*1000.0)
  time = tex*findgen(N)
  case plot_type of
    'Frequency' : begin
                  ymax = mean(freq) + 10*stddev(freq)
                  ymin = mean(freq) - 10*stddev(freq)
                  h1 = plot(time,freq, 'ro--', ytitle = 'Frequency (' +strtrim((*state).fq_xaxis_units,2) +' )',xtitle = 'Time (hrs)',title = 'Frequency vs Time',yrange = [ymin,ymax])
;                  h2 = plot(time,deriv(time,freq), 'bo--',xtitle = 'Time (hrs)', ytitle = '$df/dt$ (' +strtrim((*state).fq_xaxis_units,2) +'/hr)',title = 'Frequency Change over Time',layout = [1,2,2],/current)
                  end
    'Phase'     : begin
                  ymax = mean(phi) + 10*stddev(phi)
                  ymin = mean(phi) - 10*stddev(phi)
                  h1 = plot(time,phi,'ro--', ytitle = 'Phase (degrees)',xtitle = 'Time (hrs)',title = 'Phase vs Time',yrange = [ymin,ymax])
;                  h2 = plot(time,deriv(time,phi),'bo--',xtitle = 'Time (hrs)', ytitle = '$d\phi/dt$ (degrees/hr)',title = 'Phase Change over Time',layout = [1,2,2],/current)
                  end
    'Amplitude' : begin
                  ymax = mean(amp) + 10*stddev(amp)
                  ymin = mean(amp) - 10*stddev(amp)
                  h1 = plot(time,amp, 'ro--', ytitle = 'Amplitude',xtitle = 'Time (hrs)',title = 'Amplitude vs Time',yrange = [ymin,ymax])
;                  h2 = plot(time,deriv(time,amp), 'bo--',xtitle = 'Time (hrs)', ytitle = '$dA/dt (hr^{-1})$',title = 'Amplitude Change over Time',layout = [1,2,2],/current)
                  end
    'Linewidth' : begin
                  ymax = mean(fwhm) + 10*stddev(fwhm)
                  ymin = mean(fwhm) - 10*stddev(fwhm)
                  h1 = plot(time,fwhm,'ro--', ytitle = 'Linewidth (' +strtrim((*state).fq_xaxis_units,2) +' )',xtitle = 'Time (hrs)',title = 'Linewidth vs Time',yrange = [ymin,ymax])
;                  h2 = plot(time,deriv(time,fwhm),'bo--',xtitle = 'Time (hrs)', ytitle = 'Linewidth (' +strtrim((*state).fq_xaxis_units,2) +'/hr)',title = 'Linewidth Change over Time',layout = [1,2,2],/current)
                  end
    'Export'    : begin
                  file_name = DIALOG_PICKFILE(/WRITE, FILTER = '*.csv')
                  header_string = ['Time (s)', 'Frequency ('+ strtrim((*state).fq_xaxis_units,2) + ')', 'Phase (degrees)', 'Amplitude', 'Linewidth ('+strtrim((*state).fq_xaxis_units,2) +' )']
                  WRITE_CSV,file_name, 60*60*time, freq, phi, amp, fwhm, header = header_string 
                  end
  endcase   
  
end
pro mas_calctemp, event

  common scan_data

  widget_control, event.top, get_uvalue=state
  widget_control, event.id, get_value = solvent
  
  ppm = (*state).spect_width/(*state).base_freq
  hz = (*state).spect_width
  if (*state).fq_xaxis_units eq 'Hz' then (*state).ofq_plot->setXAxislimit, ppm
  freq1 = [2.8, 3.8]        ; 1-H methyl chemical shift : 3.31 ppm
  freq2 = [4.3,5.3]         ; 1-H hydroxyl chemical shift : 4.78 ppm
  (*state).ofq_plot->getFreqPos, freq1 , pos1  ; Local search index for freq 1
  (*state).ofq_plot->getFreqPos, freq2 , pos2  ; Local search index for freq 2
  
  ft = real_part((*state).fq_data) 
  N =  (*state).num_fids
  T = make_Array(N,/float)
  flag = 1
  for i = 0,N-1 do begin
    ind1 =  pos1[0]+where(ft[pos1[0]:pos1[1],i] eq max(ft[pos1[0]:pos1[1],i]))
    ind2 =  pos2[0]+where(ft[pos2[0]:pos2[1],i] eq max(ft[pos2[0]:pos2[1],i]))
    if ft[ind1,i] lt 0.1*max(ft[*,i]) or ft[ind2,i] lt 0.1*max(ft[*,i])  then begin
      h = dialog_message('No peak found',/error)
      flag = 0
      break
    endif
    
    (*state).ofq_plot->getDataValues, ind1, freqval = freq11, which_fid = i
    (*state).ofq_plot->getDataValues, ind2, freqval = freq21, which_fid = i
    delta = freq21 - freq11
    case solvent of
    'Pure Methanol'         : T[i] = 409.0 - 36.54*delta - 21.85*delta^2 -273.15 
    'Pure Ethylene Glycol'  : T[i] = 466.5 - 102.00*delta - 273.15
    endcase    
  endfor
  if flag eq 1 then begin
    if N gt 1 then begin
      tex = project.imndArray[project.ci].spect_rep_time*project.imndArray[project.ci].spect_num_avg/(60.0*60.0*1000.0)
      time = tex*findgen(N)
      ymax = mean(T) + 10*stddev(T)
      ymin = mean(T) - 10*stddev(T)
      h = plot(time,T, 'ro--', ytitle = 'Temperature (celsius)',xtitle = 'Time (hrs)', title = 'Temperature vs Time',yrange = [ymin,ymax])
;      h2 = plot(time,deriv(T), 'bo--',xtitle = 'Time (hrs)', ytitle = '$dT/dt$ (celsius/hr)', title = 'Temperature Change over Time',layout = [1,2,2],/current)  
    endif else h = dialog_message(['Chemical Shift Difference :'+ strtrim(delta,2) + ' ppm','Temperature '+strtrim(T[0],2)+'Celsius'],/information,/center)
      
  endif
  if (*state).fq_xaxis_units eq 'Hz' then (*state).ofq_plot->setXAxislimit, hz else (*state).ofq_plot->setXAxislimit, ppm

end

pro mas_1ds_make_stats_plot, state, plotwindow=plotwindow, psfile=psfile, tiffile=tiffile

    old_device  = !D.name
    old_p_multi = !P.multi
    
    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        old_device = !D.name
        old_p_multi = !P.multi
        print, !error_state.MSG
        help, /traceback
        return
    endif
    catch, /cancel

    if ((*state).fq_xaxis_units eq 'Hz') then begin
        native_units = 'Hz'
        alt_units = 'ppm'
        native2alt_factor = 1.0/(*state).base_freq
    endif else begin
        native_units = 'ppm'
        alt_units = 'Hz'
        native2alt_factor = (*state).base_freq
    endelse

    (*state).ofq_plot->getCursorInfo, cursor01_info=c1, cursor02_info=c2    
    min_x = min([c1.data_index, c2.data_index], max=max_x)
    
    if (abs(max_x-min_x) eq 0) then return
    if ((*state).num_fids eq 1) then return
    
    if (keyword_set(tiffile)) then begin
        set_plot, 'Z'
        device, set_pixel_depth=8
        device, set_resolution=[11*150,8.5*150], decomposed=1
        thick = 2.0
    endif else if (keyword_set(psfile)) then begin
        set_plot, 'PS', /interpolate, /copy
        device, /landscape, /color, BITS_PER_PIXEL=8, FONT_SIZE=12
        device, FILENAME=psfile, /inches, xsiz=9, ysiz=7, /TIMES, /BOLD
        thick = 2.0
    endif else if (keyword_set(plotwindow)) then begin
        wset, plotwindow
        thick = 1.0
    endif else begin
        message, /info, 'No plot output device selected.'
    endelse
    
    (*state).ofq_plot->getDataValues, c1.data_index, real=all_red_real, /all_fids
    (*state).ofq_plot->getDataValues, c2.data_index, real=all_grn_real, /all_fids
                
    old_p_multi = !p.multi
    !p.multi = [0,0,2,0,0]
    plottitle = ' at '+string(c1.freq_xval, format='(G0.5)')+' '+native_units +$
                ' ('+string(c1.freq_xval*native2alt_factor, format='(G0.5)')+' '+alt_units+')'
    plot, indgen(n_elements(all_red_real)), all_red_real, background='FFFFFF'x, color='000000'x, $
          title="Red Cursor Array"+plottitle, xrange=[0,n_elements(all_red_real)-1], $
          xtitle='Array Element', ytitle='Value', thick=thick
    plottitle = ' at '+string(c2.freq_xval, format='(G0.5)')+' '+native_units +$
                ' ('+string(c2.freq_xval*native2alt_factor, format='(G0.5)')+' '+alt_units+')'
    plot, indgen(n_elements(all_grn_real)), all_grn_real, background='FFFFFF'x, color='000000'x, $
          title="Green Cursor Array"+plottitle, xrange=[0,n_elements(all_red_real)-1], $
          xtitle='Array Element', ytitle='Value', thick=thick
          
    if (keyword_set(tiffile)) then begin
       data = reverse(tvrd(),2)
       write_tiff, tiffile, data, units=2;;, $
                   ;;xresol=(*state).save_resolution, $
                  ;; yresol=(*state).save_resolution    
    endif else if (keyword_set(psfile)) then begin
        device, /close_file       
    endif
    
    set_plot, old_device
    !P.multi = old_p_multi
    
end

pro mas_1ds_save_stats_event, event

    uname = widget_info(event.id, /uname)
    widget_control, event.top, get_uvalue=temp
    
    stats_text_wid = temp[0]
    
    if (widget_info(temp[2], /valid_id) eq 0) then begin
        message, /info, 'Bad widget ID for top-level window.'
        return
    endif
    
    widget_control, temp[2], get_uvalue=state
    
    case uname of
    
        'save_stats_text': begin
            widget_control, stats_text_wid, get_value=stats_text
            output_file = dialog_pickfile(title="Where would you like to save the text file?", $
                                          /write, filter='*.txt', default_extension='txt', $
                                          /overwrite_prompt)
            if (output_file eq '') then return
            
            openw, lun, output_file, /get_lun
            printf, lun, transpose(stats_text)
            close, lun & free_lun, lun
        end
        
        'save_stat_plot_tif': begin
            output_file = dialog_pickfile(title="Where would you like to save the TIFF file?", $
                                          /write, filter='*.tiff;*.tif', default_extension='tif', $
                                          /overwrite_prompt)
            if (output_file eq '') then return
            mas_1ds_make_stats_plot, state, tiffile=output_file
        
        end
        
        'save_stat_plot_ps': begin
            output_file = dialog_pickfile(title="Where would you like to save the PS file?", $
                                          /write, filter='*.ps', default_extension='ps', $
                                          /overwrite_prompt)
            if (output_file eq '') then return
            mas_1ds_make_stats_plot, state, psfile=output_file
        
        end            
        
    endcase
    
end

pro mas_1ds_stats_button_event, event

    common scan_data, project
    
    widget_control, event.top, get_uvalue=state
    
    nfids = project.imndarray[(*state).project_id].adim
    max_line = ''
    datastr = ''
    
    if (not widget_info((*state).base_stats_window, /valid_id)) then begin
        base_stats_window = widget_base(title="Spectroscopy Statistics", $
                                        group_leader=(*state).tlb, /column, mbar=mbar)
                                        
        file_menu = widget_button(mbar, value='File', /menu)
        
        save_stats_text = widget_button(file_menu, value='Save Stats Text...', $
                                        uname='save_stats_text', $
                                        event_pro='mas_1ds_save_stats_event')
        save_stat_plot_tif = widget_button(file_menu, value='Save Plots as TIFF...', $
                                           uname='save_stat_plot_tif', $
                                           event_pro='mas_1ds_save_stats_event')
        save_stat_plot_ps = widget_button(file_menu, value='Save Plots as PS...', $
                                          uname='save_stat_plot_ps', $
                                          event_pro='mas_1ds_save_stats_event')
        
        stats_text = widget_text(base_stats_window, xsize=60, ysize=13, /scroll, /frame)
        widget_control, base_stats_window, /realize
        xmanager, '', base_stats_window, /no_block
        
        if (nfids gt 1) then begin
            geom = widget_info(stats_text, /geometry)
            draw = widget_draw(base_stats_window, graphics_level=1, scr_xsize=geom.SCR_XSIZE, scr_ysize=2*250)
            widget_control, draw, get_value=plot_window
            widget_control, base_stats_window, set_uvalue=[stats_text, plot_window, (*state).tlb]
        endif else begin
            widget_control, base_stats_window, set_uvalue=[stats_text, 0, (*state).tlb]
        endelse
        
        (*state).base_stats_window = base_stats_window
        
    endif else begin
        base_stats_window = (*state).base_stats_window
        widget_control, base_stats_window, get_uvalue=temp
        stats_text = temp[0]
        if (nfids gt 1) then begin
            plot_window = temp[1]
        endif
    endelse
    
    (*state).ofq_plot->getCursorInfo, cursor01_info=c1, cursor02_info=c2    
    min_x_ind = min([c1.data_index, c2.data_index], max=max_x_ind)
    
    (*state).ofq_plot->getExtremeValues, real=real_extrema
    global_max_real=real_extrema.max_yval
    global_max_freqval=real_extrema.max_xval
    
    if ((*state).fq_xaxis_units eq 'Hz') then begin
        native_units = 'Hz'
        alt_units = 'ppm'
        native2alt_factor = 1.0/(*state).base_freq
    endif else begin
        native_units = 'ppm'
        alt_units = 'Hz'
        native2alt_factor = (*state).base_freq
    endelse
    
    if (abs(min_x_ind-max_x_ind) ne 0) then begin
        (*state).ofq_plot->getCursorIntervalStats, $
                  min_x=min_x, max_x=max_x, npts=npts, $
                  realmean=realmean, realstdev=realstdev, realmin=realmin, realmax=realmax, realint=realint, $
                  imagmean=imagmean, imagstdev=imagstdev, imagmin=imagmin, imagmax=imagmax, imagint=imagint, $
                  trendline=trendline, realsum=realsum, imagsum=imagsum
                
        max_line = 'Peak (real): '+string(global_max_real, format='(G0.6)')+$
                            ' at '+string(global_max_freqval, format='(G0.6)')+' '+native_units+$
                              ' ('+string(global_max_freqval*native2alt_factor, format='(G0.6)')+' '+alt_units+')'
        red_line   = 'Red cursor loc. (real)  : '+string(real_part(c1.freq_yval), format='(G0.5)')+ $
                                      ' at '+string(c1.freq_xval, format='(G0.5)')+' '+native_units +$
                                      ' ('+string(c1.freq_xval*native2alt_factor, format='(G0.5)')+' '+alt_units+')'+string(10B)
        red_line  += 'Red cursor FWHM (real)  : '+string(mas_1ds_compute_fwhm(state, 1), format='(G0.5)')+' Hz'
        green_line = 'Green cursor loc. (real): '+string(real_part(c2.freq_yval), format='(G0.5)')+ $
                                      ' at '+string(c2.freq_xval, format='(G0.5)')+' '+native_units+$
                                      ' ('+string(c2.freq_xval*native2alt_factor, format='(G0.5)')+' '+alt_units+')'+string(10B)
        green_line+= 'Green cursor FWHM (real): '+string(mas_1ds_compute_fwhm(state, 2), format='(G0.5)')+' Hz'
         mas_1ds_make_stats_plot, state, plotwindow=plot_window
         
        if (nfids gt 1) then begin

            (*state).ofq_plot->getDataValues, c1.data_index, real=all_red_real, /all_fids
            (*state).ofq_plot->getDataValues, c2.data_index, real=all_grn_real, /all_fids
            
            datastr = strarr(n_elements(all_red_real))
            
            if (min_x_ind eq max_x_ind) then begin
                cursor_integral = replicate(0.0, n_elements(all_red_real))
            endif else begin
                cursor_integral = total(real_part((*state).fq_data[min_x_ind:max_x_ind, *]),1)
            endelse
            
            for j = 0, n_elements(all_red_real)-1 do begin
                datastr[j] = strjoin([string(j, format='(I0)'), $
                                      string(all_red_real[j], format='(G0.5)'), $
                                      string(all_grn_real[j], format='(G0.5)'), $
                                      string(cursor_integral[j], format='(G0.5)')], $
                                  string(09B))
            endfor
            
        endif
       
        stats = strjoin(['Type','Min', 'Max', 'N', 'Mean', 'Std Dev', 'Sum', 'Integral'], string(09B))
        
        valstr_real = [ 'Real', string(realmin, format='(G0.5)'), $
                        string(realmax, format='(G0.5)'), $
                        string(npts, format='(I0)'), $
                        string(realmean, format='(G0.5)'), $
                        string(realstdev, format='(G0.5)'), $
                        string(realsum, format='(G0.8)'), $
                        string(realint, format='(G0.8)') ]

        valstr_imag = [ 'Imag', string(imagmin, format='(G0.5)'), $
                        string(imagmax, format='(G0.5)'), $
                        string(npts, format='(I0)'), $
                        string(imagmean, format='(G0.5)'), $
                        string(imagstdev, format='(G0.5)'), $
                        string(imagsum, format='(G0.8)'), $
                        string(imagint, format='(G0.8)') ]
        
        adim = project.imndarray[(*state).project_id].adim
        adim_start = project.procpramarray[(*state).project_id].adim_start

        if (nfids gt 1) then begin
            widget_control, stats_text, set_value=[max_line, '', 'Cursor Information', '-------------------------',$
                                                   red_line, green_line, '', 'Interval Statistics', '-------------------------', $
                                                   'Trendline slope: '+string(trendline[1], format='(G0.5)')+', y-intercept: '+string(trendline[0], format='(G0.5)'), $
                                                   stats, $
                                                   strjoin(valstr_real, string(09B)), $
                                                   strjoin(valstr_imag, string(09B)), $
                                                   '', 'Multi-fid Statistics', '-------------------------',$
                                                   strjoin(['Array', 'Red Cursor Real', 'Green Cursor Real', 'Interval Sum Real'], string(09B)), $
                                                   datastr ]
         
        
        endif else begin
            widget_control, stats_text, set_value=[max_line, '', 'Cursor Information', '-------------------------',$
                                                   red_line, green_line, '', 'Interval Statistics', '-------------------------',$
                                                   'Trendline slope: '+string(trendline[1], format='(G0.5)')+', y-intercept: '+string(trendline[0], format='(G0.5)'), $
                                                   stats, $
                                                   strjoin(valstr_real, string(09B)), $
                                                   strjoin(valstr_imag, string(09B)) ]
        endelse
        widget_control, base_stats_window, /realize
    endif

end

pro mas_1ds_exit_event, event

    widget_control, event.top, /destroy

end

pro mas_display_1d_spect, ts_data, title=title

    mas_1ds_window_init, ts_data, title=title, tlb=tlb

end
