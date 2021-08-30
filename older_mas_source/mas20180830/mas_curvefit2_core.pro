;+
; :Description:
;    Uses image histogram to discard outlier intensity values for better image
;    scaling for display purposes.
;
; :Params:
;    image   - the image data, unscaled
;    low_pct - the low intensity percent. intensity values less than this
;              may be set to zero
;    hi_pct  - the high intensity percent. values greater than this may be
;              set to 255
;
; :Author: wtriplett
;-
function mas_cfit2_robust_intensity, image, low_pct, hi_pct, mask=mask
    
    if (n_elements(low_pct) eq 0) then low_pct = 0.03
    if (n_elements(hi_pct) eq 0) then hi_pct = 0.97
    
    if (n_elements(mask) ne 0) then begin
        keep = where(mask ne 0, n_keep)
        if (n_keep gt 0) then begin
        endif
    endif else begin
        keep = lindgen(n_elements(image))
    endelse
    
    hist = histogram((float(image))[keep], nbins=900, locations=loc)
    hist /= total(hist)
    pdf = fltarr(900)
    for h = 0, 900-1 do pdf[h] = total(hist[0:h])
    min_cut = min(loc[where(pdf ge low_pct)])
    max_cut = min(loc[where(pdf ge hi_pct)])

    return, [min_cut, max_cut]

end

;+
; :Description:
;    Creates an MPFIT-compatible parinfo struct from fitting input parameters. Each
;    argument should have the same number of values.
;
; :Params:
;    values - the initial values for each parameter 
;    fixed  - a 1,0 array indicating shich parameters are fixed and which are not.
;    names  - a string array containing the common names of the parameters (T2, T1, ADC, etc.)
;    units  - the units attached to the parameters (ms, etc).
;
; :Author: wtriplett
;-
function mas_cfit2_make_basic_parinfo, values, fixed, names, units

    nvalues = n_elements(values)
    nfixed  = n_elements(fixed)
    nnames  = n_elements(names)
    nunits  = n_elements(units)
    
    if (nvalues eq 0) then message, 'VALUE parameter is required'
    if (nfixed eq 0) then message,  'FIXED parameter is required'
    if (nnames  eq 0) then message, 'NAME parameter is required'
    if (nunits  eq 0) then message, 'UNITS parameter is required'

    if (nvalues ne nfixed or nfixed ne nnames or nnames ne nunits) then begin
        message, 'values, fixed, names, and units must have the same number of elements.'
    endif
    
    parinfo = replicate({ VALUE: 0.0, $
                          FIXED: 0, $
                          PARNAME: '', $
                          PARUNITS: '', $
                          LIMITED: [0, 0], $
                          LIMITS: [0.0, 0.0] }, nvalues)
    
    for v = 0, nvalues-1 do begin
        parinfo[v].value    = values[v]
        parinfo[v].fixed    = fixed[v]
        parinfo[v].parname  = names[v]
        parinfo[v].parunits = units[v]
    endfor

    return, parinfo
    
end

;+
; :Description:
;    Returns a single slice, based on the selection in MAS. Alternatively,
;    It is possible to specify which slice. It returns all of the array
;    data for the slice.
;
; :Params:
;    for_slice - a number to indicate which slice to return, overriding the
;                slice selected in MAS gui.
;
; :Author: wtriplett
; 
; Modifications (mkulam)
; Added a feature which allows to select which arrays to include in the fit

function mas_cfit2_get_slice_data, for_slice, ci = ci

    common scan_data, project
    
    if ~n_elements(ci) then ci = project.ci
    
    pdata = project.dataarray[ci].state1
    
    if (not ptr_valid(pdata)) then begin
        message, 'Scan data has not been loaded.'
    endif
    
    r = project.procpramarray[ci].fdim_start
    p = project.procpramarray[ci].pdim_start
    s = project.procpramarray[ci].sdim_start
    
    case project.procpramarray[ci].slice_axis of
        0: sl_num = s
        1: sl_num = p
        2: sl_num = r
    endcase

    if (n_elements(for_slice) ne 0) then begin
        sl_num = for_slice
    endif
    
    a = project.procpramarray[ci].array_select
    
    case project.procpramarray[ci].slice_axis of
        0: if ptr_valid(a) then slice = reform((*pdata)[*,*,sl_num,*a]) else slice = reform((*pdata)[*,*,sl_num,*])        
        1: if ptr_valid(a) then slice = reform((*pdata)[*,sl_num,*,*a]) else slice = reform((*pdata)[*,sl_num,*,*])
        2: if ptr_valid(a) then slice = reform((*pdata)[sl_num,*,*,*a]) else slice = reform((*pdata)[sl_num,*,*,*])
    endcase

    return, slice

end

;+
; :Description:
;    Executes a curvefitting operation and returns the result.
;
; :Params:
;    func - a string name of the function to fit to.
;    X    - the X values 
;    Y    - the Y values
;    err  - the measurement errors for each data point. (optional)
;    parinfo - The parinfo struct containing the information about
;              whether the parameter is fixed, the name of the 
;              parameter, the starting value and the units.
;
; :Keywords:
;    fit_error   - this will be set to 1 if something happened in the
;                  fitting process that caused an error.
;    use_powell  - Set this to use the POWELL optimizer instead of MPFIT.
;                  Note that some features may not be available, like
;                  sigmas on parameters.
;    powell_objective - if use_powell is enabled, this must be set to the
;                       string name of an objective function to minimize.
;                       This is generally different from the model function, since
;                       it will return a scalar value.
;
; :Author: wtriplett
;-
function mas_cfit2_run_fit, func, X, Y, err, parinfo, fit_error=fit_error, $
                            use_powell=use_powell, powell_objective=powell_objective

    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        fit_error = 1
        return, ptr_new()
    endif
    
    fit_error = 0

    if (keyword_set(use_powell)) then begin
        if (n_elements(powell_objective) eq 0) then begin
            message, 'Use of POWELL optimizer requires a specified objective function.'
        endif
    endif
    
    if (n_elements(err) eq 0) then err = replicate(1.0, n_elements(X))
    
    if (keyword_set(use_powell)) then begin
    
        ;; we will set up a common block for powell
        common __curvefit2_powell, fit_params
        fit_params = { measurements: Y, params: X, err: err }
        A   = float(parinfo.value)
        Xi  = diag_matrix(1-float(parinfo.fixed))
        dof = n_elements(X)-n_elements(a)

        !EXCEPT=0
        powell, A, Xi, 1e-7, fmin, powell_objective
        !EXCEPT=1

        expected  = call_function(func, X, A)
        if (total(1-finite(expected)) ne 0) then begin
            message, /info, 'Fit parameters from POWELL caused infinte y-value.'
        endif
        resid = Y - expected
        
        sigma1 = abs(replicate(!VALUES.F_NAN, n_elements(parinfo))) ;; make these meaningless
        params = A
        
    endif else begin
        
        params = mpfitfun(func, X, Y, err, void, /quiet, $
                          parinfo=parinfo, perror=sigma1)
        
        resid = Y - call_function(func, X, params)
    
    endelse
    
    ;; apply errors to the reported chi^2
    if (n_elements(err) ne 0) then begin
        chisq = total( (resid/err)^2 )
    endif else begin
        chisq = total(resid^2)
    endelse
    
    ;; set up a fit results structure. some input
    ;; parameters are repeated for convienuence
    fit_result = { func: func, $             ;; model function 
                   X : X, $                  ;; X values
                   Y : Y, $                  ;; Y values
                   params : params, $        ;; result of the fit
                   sigma1 : sigma1, $        ;; 1-sigma errors
                   resid  : resid, $         ;; residuals
                   chisq  : chisq, $         ;; chi^2 value
                   parinfo: parinfo, $       ;; parinfo (copied from input)
                   xtitle: 'X', $            ;; Some titles, you can change after.
                   ytitle: 'Y', $
                   title: 'Fit to '+func, $
                   alt_title: '' }

    return, ptr_new(fit_result, /no_copy)

end

;+
; :Description:
;    Executes an image fit operation.
;
; :Params:
;    func - a string name of the function to fit to.
;    X    - the X values 
;    slice - the image slice data to fit
;    err  - the measurement errors for each data point. (optional)
;    parinfo - The parinfo struct containing the information about
;              the starting values for each voxel
;              whether the parameter is fixed, the name of the 
;              parameter, the starting value and the units.
;    start_map - this is a 2D array with the same dimensions as slice containing
;
; :Keywords:
;    thr_pct  - a percent value (0 .. 1.0) indicated the threshold value to 
;               fit in the image. don't fit if value is < max(slice)*thr_pct 
;    progressbar - set this keyword to have a progress bar pop up
;    use_powell  - set this to use the POWELL optimizer. See 'mas_cfit2_run_fit'
;                  for information about this parameter
;    powell_objective - must be set to the string name of the powell objective function
;
; :Author: wtriplett
; Modification (Magdoom)
; - Added functionality to perform the fit in selected voxels given by mask = 1

function mas_cfit2_run_image_fit, func, X, slice, start_map, err, parinfo, $
                                  thr_pct=thr_pct, progressbar=progressbar, $
                                  use_powell=use_powell, powell_objective=powell_objective, mask = mask

    dims = size(slice, /dimensions)
    
    n_unfixed = total(parinfo.fixed eq 0)
    
    param_maps = fltarr(dims[0], dims[1], n_unfixed)
    chisq_map  = fltarr(dims[0], dims[1])
    thr_mask   = bytarr(dims[0], dims[1]) + 1B
    
    thr_pct = n_elements(thr_pct) eq 0 ? 0.0001 : float(thr_pct)
    thr_max = max(abs(slice))
    thr_val = thr_max * thr_pct
   
    sorted = sort(X)
    X_srt = X[sorted]
    
    if (keyword_set(progressbar)) then begin
        pbar = obj_new('progressbar', text='Image Fit Running...', color='Red', /fast_loop)
        pbar->Start
    endif
    n = 0L
    
    n_parinfo = n_elements(parinfo)
    t = systime(1)
    for slx = 0, dims[0]-1 do begin
        for sly = 0, dims[1]-1 do begin
            if isa(mask) then if mask[slx,sly] eq 0 then continue
            n++
            Y = (reform(slice[slx, sly, *]))[sorted]
            if (max(abs(y)) lt thr_val) then begin
                thr_mask[slx, sly] = 0 ;; record the thr mask 
                continue
            endif
            
            for p = 0, n_parinfo-1 do begin
                parinfo[p].value = start_map[slx,sly,p]
            endfor

            temp = mas_cfit2_run_fit(func, X_srt, Y, err, parinfo, fit_error=fe, $
                                     use_powell=use_powell, $
                                     powell_objective=powell_objective)
            if (fe eq 1) then continue
            
            chisq_map[slx, sly] = (*temp).chisq
            param_num = 0
            for f = 0, n_parinfo-1 do begin
                if (parinfo[f].fixed eq 0) then begin
                    param_maps[slx, sly, param_num++] = (*temp).params[f]
                endif
            endfor
            
            ptr_free, temp
            if (obj_valid(pbar)) then begin
                pbar->update, 100.0 * float(n)/(dims[0]*dims[1])
                if (n mod 100 eq 0) then begin
                    if (pbar->CheckCancel()) then begin
                        pbar->Destroy
                        return, ptr_new()
                    endif
                endif
            endif
            
        endfor
    endfor
    message, /info, 'Elapsed time: '+string(systime(1)-t, format='(%"%0.4f")')

    if (obj_valid(pbar)) then pbar->Destroy
    
    return, ptr_new({ parinfo: parinfo, $
                      param_maps: param_maps, $
                      thresh_mask: thr_mask, $
                      chisq_map: chisq_map }, /no_copy)
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;+
; :Description:
;    Frees the ROI data structure.
;
; :Params:
;    roi_data - data structure to free
;
; :Author: wtriplett
;-
pro mas_cfit2_roi_free_data, roi_data

    for i = 0, n_elements(roi_data)-1 do begin
        ptr_free, roi_data[i].pmask
        ptr_free, roi_data[i].pdata
    endfor

end

;+
; :Description:
;    Retrieves data from MAS which is contained in the currently
;    selected ROI set.
;
; :Params:
;    measurement_type - this is a keyword describing the type of function
;                       to apply to the pixels in the ROI. For example,
;                       0 - the mean of the pixels for each array element 
;                           will be computed
;                       1 - the mean of the absolute value will be computed
;                       2 - the rician noise correction function will be
;                           applied.
;
; :Keywords:
;    for_slice - set this to the slice number to apply the ROIs to. Otherwise
;                the ROIs are appliede to the currently selected slice in MAS.
;
; :Author: wtriplett
;-
function mas_cfit2_roi_get_data, measurement_type, for_slice=for_slice

    common scan_data, project
    ci = project.ci
    
    num_rois = mas_roi_get_number_of_rois()
    if (num_rois eq 0) then begin
        message, 'The current ROI set contains no regions.'
    endif
    
    mtype     = n_elements(measurement_type) ? measurement_type : 0
    if ptr_valid(project.procpramArray[ci].array_select) then adim = (size(*project.procpramArray[ci].array_select,/Dimension))[0] else adim = project.imndarray[ci].adim
    p_masks   = ptrarr(num_rois)
    roi_names = strarr(num_rois)
    
    case project.procpramarray[ci].slice_axis of
        0: active_slice = project.procpramarray[ci].sdim_start
        1: active_slice = project.procpramarray[ci].pdim_start
        2: active_slice = project.procpramarray[ci].fdim_start
    endcase

    if (n_elements(for_slice) ne 0) then begin
        if (for_slice lt 0 or for_slice ge sdim_start) then begin
            message, 'Requested slice number if out of range: ['+ $
                      string(0, sdim_start-1, format='(%"%d,$d")')
        endif
    endif else begin
        for_slice = active_slice
    endelse
    
    roi_data_str = { name: '', $
                     pmask: ptr_new(), $
                     pdata: ptr_new(), $
                     slice_num: for_slice, $
                     measurement_func: '', $
                     measurement: fltarr(adim), $
                     measure_err: fltarr(adim), $
                     data_mean: fltarr(adim), $
                     data_stdev: fltarr(adim), $
                     data_signal : fltarr(adim), $
                     data_absmean : fltarr(adim), $
                     indep_var: ptr_new(), data_type: '' }
                     
    output       = replicate(roi_data_str, num_rois)
    
    slice    = mas_cfit2_get_slice_data(for_slice)
    size_str = size(slice, /structure)
    if (size_str.type_name eq 'COMPLEX') then begin
        message, /info, 'WARNING: data loaded is of COMPLEX type -- REAL_PART will be used.'
        slice = real_part(slice)
        size_str = size(slice, /structure)
    endif
    
    for r = 0, num_rois-1 do begin
    
        output[r].name  = mas_roi_get_current_name(region_num=r)
        
        pmask = mas_roi_get_current_mask(region_num=r, /no_transform)
        output[r].pmask = pmask
        
        mask_active = where(*pmask ne 0, active_count)
        if (active_count eq 0) then begin
            message, /info, string(r, output[r].name, format='(%"ROI %s (#%d) has empty mask")')
            continue
        endif
        
        pdata = ptr_new(fltarr(active_count, adim))

        for a = 0,adim-1 do begin
        
            (*pdata)[*,a] = reform( (slice[*,*,a])[mask_active] )
            
            output[r].data_mean[a]  = mean((*pdata)[*,a])
            output[r].data_stdev[a] = stdev((*pdata)[*,a])
            
            temp = abs( (*pdata)[*,a] )
            var = variance(temp)
            n   = float(n_elements(temp)) 
            ssq = total( temp^2 )
            output[r].data_signal[a] = sqrt(ssq/n - 2.0*var) 
            output[r].data_absmean[a] = mean(temporary(temp))
            
        endfor
        
        case mtype of
            0: begin
                message, /info, 'Using mean(data)'
                output[r].measurement = output[r].data_mean
                output[r].measurement_func = 'mean(data)'
            end
            1: begin
                message, /info, 'Using mean(abs(data))'
                output[r].measurement = output[r].data_absmean
                output[r].measurement_func = 'mean(abs(data))'
               end
            2: begin
                message, /info, 'Using signal(abs(data))'
                output[r].measurement = output[r].data_signal
                output[r].measurement_func = 'signal(abs(data))'
               end
            
            else: message, 'Unknown measurement type: '+string(mtype, format='(%"%d")')
        endcase
        
        output[r].pdata = pdata
        output[r].data_type = size_str.type_name

    endfor

    return, output
    
end


;+
; :Description:
;    Using the results of the fit, and the input ROI data, this procedure
;    computes some summary statistics and returns them as an array of lines.
;
; :Params:
;    fit_result - a fit result structure as returned from mas_cfit2_run_fit
;    roi_data - an roi data structure as returned from mas_cfit2_roi_get_data
;
; :Keywords:
;    stat_lines - the stat lines are returned in this keyword parameter.
;
; :Author: wtriplett
;-
pro mas_cfit2_roi_make_stats, fit_result, roi_data, stat_lines=stat_lines

    common scan_data, project
    ci = project.ci
    
    nlines = 0
    tabch = string(09B)
    slice_axes = [ 'Freq-Phase', 'Freq-Slice', 'Phase-Slice' ]
    slice_axis_str = slice_axes[project.procpramarray[ci].slice_axis]
    
    interp_str = string(project.procpramarray[ci].freq_interp, $
                        project.procpramarray[ci].phase_interp,$
                        project.procpramarray[ci].slice_interp,$
                        format='(%"Interpolation (RPS): %d,%d,%d")')
    
    fit_func = (*fit_result[0]).func
    
    all_lines = [ 'Curve fitting results for '+project.imndarray[ci].display_name ]
    all_lines = [ all_lines, 'Slice #'+string(roi_data[0].slice_num, format='(%"%03d")')+tabch+'Model function: '+fit_func ]
    all_lines = [ all_lines, 'Slice Axis: '+slice_axis_str+tabch+interp_str ]
    
    parinfo = (*fit_result[0]).parinfo    
    unfixed = where(parinfo.fixed eq 0)
    title   = 'Region'+tabch
    for t = 0, n_elements(parinfo)-1 do begin
        if (parinfo[t].fixed eq 0) then begin
            name = parinfo[t].parname
            unit = parinfo[t].parunits
            if (unit ne '') then unit = ' ('+unit+')'
            title += strjoin([name+unit, 'Sigma ('+parinfo[t].parname+')'], tabch) + tabch
        endif
    endfor
    title += 'Chisquared'
    all_lines = [all_lines, title]
    nlines++
    
    nlines = 3
    for f = 0, n_elements(fit_result)-1 do begin

        line = [ roi_data[f].name+tabch]
        
        parinfo = (*fit_result[f]).parinfo
                
        for p = 0, n_elements(parinfo)-1 do begin
        
            if (parinfo[p].fixed) then continue
            
            line += string((*fit_result[f]).params[p], (*fit_result[f]).sigma1[p], $
                           format='(%"%0.3e\t%0.3e\t")')
        
        endfor
    
        line_out = line + string((*fit_result[f]).chisq, format='(%"%0.3e")')
        
        all_lines = (n_elements(all_lines) ne 0) ? [ all_lines, line_out ] : line_out
        
    endfor

    stat_lines = all_lines    

end

;+
; :Description:
;    Given a fit result, this procedure will create a plot.
;
; :Params:
;    fit_result - the fit result as returned from 'mas_cfit2_run_fit'
;
; :Keywords:
;    yrange      - the y range of the plot [min,max]
;    region_num  - The region number index into fit_result to show
;    region_name - The name associate with this region number
;    no_legend   - Suppress the legend on the plot
;    overplot    - set this to keep whatever is on the current plot
;                  without erasing
;    all_in_one  - this keyword will plot all regions on a single plot
;                  the layout is slightly different in this configuration.
;    _EXTRA      - Any extra keyword are passed to plot 
;
; :Author: wtriplett
;-
pro mas_cfit2_roi_make_plot, fit_result, yrange=c_yrange, $
                             region_num=region_num, region_name=region_name, $
                             no_legend=no_legend, overplot=overplot, $
                             all_in_one=all_in_one, _EXTRA=ex

    colors   = [ 'FF0000'x, '00FF00'x, '0000FF'x, $
                 'FF00FF'x, '0066FF'x, '00FFFF'x, $
                 '67EB15'x ]
    pscolors = [ 0, 91, 22, 79, 132, 50, 82 ]
    
    if (1 or !D.name eq 'PS') then begin ;; force white background color scheme for PS and window plots.
        decomp = 0
    endif else begin
        device, get_decomposed=decomp
    endelse
    
    if (decomp eq 0) then begin
        loadct, 5, /silent
        use_colors = pscolors
        color_white = 255
        color_black = 0
        color_fixed = 91
    endif else begin
        use_colors = colors
        color_white = 'ffffff'x
        color_black = '000000'x
        color_fixed = '0000ff'x
    endelse
    
    psyms  = [ 1, 2, 4, 5, 6, 7, 5 ]
    n_colors = n_elements(use_colors)

    parinfo = (*fit_result).parinfo
    
    X_hi_res = findgen(200)/200 * max((*fit_result).X)*1.25
    Y_hi_res = call_function((*fit_result).func, X_hi_res, (*fit_result).params)
    
    max_x_hi_res = max(X_hi_res, min=min_x_hi_res)
    max_y_hi_res = max(Y_hi_res, min=min_y_hi_res)
    
    if (n_elements(c_yrange) ne 2) then begin
        yrange = [min_y_hi_res,max_y_hi_res*1.1]
    endif else begin
        yrange = c_yrange
    endelse
    
    if (not keyword_set(no_legend)) then xmargin = [9,23]

    if (1 or !D.name eq 'PS') then begin ;; forcing white BG 
        plot_fg = color_black
        plot_bg = color_white
    endif else begin
        plot_fg = color_white
        plot_bg = color_black
    endelse

    plot_color = (keyword_set(all_in_one)) ? use_colors[region_num mod n_colors] : plot_fg
    plot_thick = 2.0
    
    if (not keyword_set(overplot)) then begin
        plot, X_hi_res, Y_hi_res, /nodata, color=plot_fg, background=plot_bg, $
              xrange=[min_x_hi_res,min_x_hi_res], $
              yrange=yrange, $
              xticklen=-0.015, yticklen=-0.015, $
              xmargin=xmargin, thick=plot_thick, _EXTRA=ex
    endif
    
    oplot, X_hi_res, Y_hi_res, color=plot_color, thick=plot_thick, _EXTRA=ex
    
    oplot, (*fit_result).X, (*fit_result).Y, psym=psyms[region_num mod n_colors], symsize=1.5, thick=plot_thick, color=plot_color

    if (keyword_set(no_legend)) then return

    if (not keyword_set(overplot)) then begin
        oplot, 1.01*replicate(max_x_hi_res,2), $
               1.05*[min_y_hi_res,max_y_hi_res], linestyle=1, /noclip
    endif
        
    X_loc  = !P.clip[2]+!D.X_CH_SIZE*3
    Y_loc  = !P.clip[3]-!D.Y_CH_SIZE*2
    y_mult = !D.Y_CH_SIZE * 1.5

    if (keyword_set(all_in_one)) then begin
        xyouts, X_loc, Y_loc - (1+region_num)*y_mult, region_name, color=plot_color, /device
        return
    endif
    
    for p = 0, n_elements(parinfo)-1 do begin

        name      = strmid(parinfo[p].parname, 0, 14)
        fit_value = (*fit_result).params[p]
        sigma1    = (*fit_result).sigma1[p]
        
        format = '(%"%s%0.2e %s")'
        color  = (parinfo[p].fixed) ? color_fixed : plot_fg

        if (p gt 0) then begin
            xyouts, X_loc, Y_loc, strjoin(replicate('_', 20), ''), /device, color=color, charsize=charsize, /noclip
            Y_loc -= 1.25*y_mult
        endif
        
        xyouts, X_loc, Y_loc, 'Parameter: '+name, color=color, /device, charsize=charsize, /noclip
        Y_loc -= y_mult
        
        xyouts, X_loc, Y_loc, $
                string('Value: ', fit_value, parinfo[p].parunits, format=format),$
                color=color, /device, charsize=charsize, /noclip
        Y_loc -= y_mult
        
        xyouts, X_loc, Y_loc, $
                string('Sigma: ', sigma1, parinfo[p].parunits, format=format),$
                color=color, /device, charsize=charsize, /noclip
        Y_loc -= 1.25*y_mult

    endfor
    
    xyouts, X_loc, Y_loc, strjoin(replicate('_', 20), ''), /device, charsize=charsize, color=color, /noclip
    Y_loc -= 1.25*y_mult
    
    xyouts, X_loc, Y_loc, /device, string((*fit_result).chisq, format='(%"ChiSQ: %0.2e")'), color=color, charsize=charsize, /noclip

end

;pro mas_cfit2_img_make_plot, fit_result, yrange=c_yrange, $
;                             region_num=region_num, region_name=region_name, $
;                             no_legend=no_legend, overplot=overplot, $
;                             all_in_one=all_in_one, _EXTRA=ex
;
;    colors = [ 'FF0000'x, '00FF00'x, '0000FF'x, $
;               'FF00FF'x, '0066FF'x, '00FFFF'x, $
;               '67EB15'x ]
;    psyms  = [ 1, 2, 4, 5, 6, 7, 5 ]
;    n_colors = n_elements(colors)
;
;    parinfo = (*fit_result).parinfo
;    
;    if (not keyword_set(no_legend)) then xmargin = [9,23]
;
;    if (!D.name eq 'PS') then begin
;        plot_fg = '000000'X
;    endif else begin
;        plot_fg = 'FFFFFF'X
;    endelse
;
;    plot_color = (keyword_set(all_in_one)) ? colors[region_num mod n_colors] : plot_fg
;    
;    if (not keyword_set(overplot)) then begin
;        plot, X_hi_res, Y_hi_res, /nodata, color=plot_fg, $
;              xrange=[min_x_hi_res,min_x_hi_res], $
;              yrange=yrange, $
;              xticklen=-0.015, yticklen=-0.015, $
;              xmargin=xmargin, _EXTRA=ex
;    endif
;    
;    oplot, X_hi_res, Y_hi_res, color=plot_color, _EXTRA=ex
;    
;    oplot, (*fit_result).X, (*fit_result).Y, psym=psyms[region_num mod n_colors], symsize=1.0, color=plot_color
;
;    if (keyword_set(no_legend)) then return
;
;    if (not keyword_set(overplot)) then begin
;        oplot, 1.01*replicate(max_x_hi_res,2), $
;               1.05*[min_y_hi_res,max_y_hi_res], linestyle=1, /noclip
;    endif
;        
;    X_loc  = !P.clip[2]+!D.X_CH_SIZE*3
;    Y_loc  = !P.clip[3]-!D.Y_CH_SIZE*2
;    y_mult = !D.Y_CH_SIZE * 1.5
;
;    if (keyword_set(all_in_one)) then begin
;        xyouts, X_loc, Y_loc - (1+region_num)*y_mult, region_name, color=plot_color, /device
;        return
;    endif
;    
;    for p = 0, n_elements(parinfo)-1 do begin
;
;        name      = strmid(parinfo[p].parname, 0, 14)
;        fit_value = (*fit_result).params[p]
;        sigma1    = (*fit_result).sigma1[p]
;        
;        format = '(%"%s%0.2e %s")'
;        color  = (parinfo[p].fixed) ? '0000FF'x : (!D.name eq 'PS') ? '000000'x : 'FFFFFF'x
;
;        if (p gt 0) then begin
;            xyouts, X_loc, Y_loc, strjoin(replicate('_', 20), ''), /device, charsize=charsize, /noclip
;            Y_loc -= 1.25*y_mult
;        endif
;        
;        xyouts, X_loc, Y_loc, 'Parameter: '+name, color=color, /device, charsize=charsize, /noclip
;        Y_loc -= y_mult
;        
;        xyouts, X_loc, Y_loc, $
;                string('Value: ', fit_value, parinfo[p].parunits, format=format),$
;                color=color, /device, charsize=charsize, /noclip
;        Y_loc -= y_mult
;        
;        xyouts, X_loc, Y_loc, $
;                string('Sigma: ', sigma1, parinfo[p].parunits, format=format),$
;                color=color, /device, charsize=charsize, /noclip
;        Y_loc -= 1.25*y_mult
;
;    endfor
;    
;    xyouts, X_loc, Y_loc, strjoin(replicate('_', 20), ''), /device, charsize=charsize, /noclip
;    Y_loc -= 1.25*y_mult
;    
;    xyouts, X_loc, Y_loc, /device, string((*fit_result).chisq, format='(%"ChiSQ: %0.2e")'), charsize=charsize, /noclip
;
;end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;+
; :Description:
;    This will exectute and display the results of an ROI fit, and send
;    the results to the appropriate plot device.
;
; :Params:
;    model_pro - a string containing the procedure model to fit.
;
; :Keywords:
;    psfile           - some options controlling the output behavior,
;    multiplot          mostly they should be self explanatory
;    all_in_one
;    combined_yscale
;    no_legend
;    no_plot
;    roi_stats
;    measurement_type
;
; :Author: wtriplett
;-
pro mas_cfit2_roi_fit, model_pro, psfile=psfile, $
                       multiplot=multiplot, $
                       all_in_one=all_in_one, $
                       combined_yscale=combined_yscale, $
                       no_legend=no_legend, $
                       no_plot=no_plot, $
                       roi_stats=roi_stats, $
                       measurement_type=measurement_type

    common scan_data, project
    ci = project.ci

    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        loadct, 0, /silent
        void = dialog_message(['Unable to perform ROI curve fit:', $
                               !ERROR_STATE.msg], /error, /center)
        return
    endif
    
    if (n_elements(model_pro) eq 0) then begin
        message, 'No model specified.'
    endif
    
    set_plot, (!VERSION.OS_FAMILY eq 'unix' ) ? 'X' : 'WIN'

    m_type     = n_elements(measurement_type) ? measurement_type : 0
    roi_data   = mas_cfit2_roi_get_data(measurement_type)
    roi_wanted = (n_elements(roi_number) gt 0) ? roi_number > 0 : 0
    num_rois   = n_elements(roi_data)
    
    if (num_rois-1 lt roi_wanted) then begin
        message, string(roi_wanted, num_rois-1, $
                        format='(%"Chosen ROI out of range. (wanted = %d, max = %d)")')
    endif
    
    call_procedure, model_pro, roi_data, fit_result=fit_result
    
    if (not arg_present(roi_stats) and keyword_set(no_plot)) then begin
        goto, DONE
    endif else if (keyword_set(no_plot)) then begin
        goto, STATS
    endif
    
    if (n_elements(psfile) ne 0) then begin
        PS=1
        set_plot, 'PS'
        filename=psfile
        device, filename=filename, bits_per_pixel = 8, font_size=11, language_level=2
        device,  /landscape, /color, /inches, /TIMES, /BOLD, $
                xsize=11, ysize=8.5, xoffset=0.25, $;;decomposed=0, $
                yoffset=11-(0.25), scale_factor=0.95
    endif else begin
        PS=0
        set_plot, (!VERSION.OS_FAMILY eq 'unix' ) ? 'X' : 'WIN'
        device, decomposed=0
        wid = 5L
        device, get_screen_size=scrn_size
        win_xsize = scrn_size[0]*0.6
        win_ysize = scrn_size[1]*0.7
        window, wid++, xsize=win_xsize, ysize=win_ysize
    endelse
    
    if (keyword_set(multiplot) and keyword_set(all_in_one)) then begin
        message, /info, 'all_in_one and multiplot are mutually exclusive options, disregarding multiplot.'
        a = temporary(multiplot)
    endif
        
    if (keyword_set(multiplot) and num_rois gt 1) then begin
        case 1 of
            num_rois ge 3: !P.multi = [0, 2, 2]
            num_rois gt 1: !P.multi = [0, 1, 2]
            else :         !P.multi = 0
        endcase
    endif else begin
        !P.multi = 0
    endelse
    
    if (keyword_set(all_in_one)) then begin
        combined_yscale = 1
        noerase         = 1
        !P.multi        = 0
    endif else begin
        noerase = 0
    endelse
    
    if (keyword_set(combined_yscale)) then begin
        yrange_check = fltarr(2, num_rois)
        for r = 0, num_rois-1 do begin
            mn = min((*fit_result[r]).y, max=mx)
            yrange_check[*,r] = [ mn, mx ]
        endfor
        yrange = [ 0.9*min(yrange_check[0,*]), max(yrange_check[1,*])*1.1 ]
    endif
    
    for r = 0, num_rois-1 do begin

        xtitle   = (*fit_result[r]).xtitle
        ytitle   = (*fit_result[r]).ytitle + ' ['+roi_data[r].measurement_func+']'

        if (keyword_set(all_in_one)) then begin
            title = (*fit_result[r]).title + ' (All ROIs)'
        endif else begin
            title = (*fit_result[r]).title + ' ('+roi_data[r].name+')'
        endelse

        if (not keyword_set(multiplot) and PS eq 0 and r gt 0) then begin
            if (not keyword_set(all_in_one)) then begin
                window, wid++, xsize=win_xsize, ysize=win_ysize
            endif
        endif 
        
        if (keyword_set(all_in_one) and r gt 0) then begin
            title    = ''
            xtitle   = ''
            ytitle   = ''
            overplot = 1
            noerase  = 1
        endif else begin
            overplot = 0
            noerase  = 0
        endelse

        mas_cfit2_roi_make_plot, fit_result[r], yrange=yrange, noerase=noerase, $
                                 region_num=r, region_name=roi_data[r].name, $
                                 title=title, xtitle=xtitle, ytitle=ytitle, color=color, $
                                 xstyle=1+8, ystyle=1+8, linestyle=2, overplot=overplot, $
                                 no_legend=keyword_set(no_legend), $
                                 all_in_one=keyword_set(all_in_one)
                                   
        if (keyword_set(multiplot) and !P.multi[0] eq 0 and r ne num_rois-1) then begin
            if (PS eq 0 and not keyword_set(all_in_one)) then begin
                window, wid++, xsize=win_xsize, ysize=win_ysize
            endif
        endif
        
    endfor

    if (PS) then begin
        device, /close_file
        set_plot, (!VERSION.OS_FAMILY eq 'unix' ) ? 'X' : 'WIN'
    endif else begin
        loadct, 0, /silent
    endelse

    STATS:
    
    if (arg_present(roi_stats)) then begin
        mas_cfit2_roi_make_stats, fit_result, roi_data, $
                                  stat_lines=roi_stats
    endif
    
    DONE:
    
    mas_cfit2_roi_free_data, roi_data
    ptr_free, fit_result
    !P.multi = 0
    
end
