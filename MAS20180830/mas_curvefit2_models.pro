function CF2_T1_InvRecovFunction_SE_ABS_POWELL, A

    common __curvefit2_powell, fit_params
    
    TI  = fit_params.params
    Y   = fit_params.measurements
    err = fit_params.err
    
    ;; A[0] - steady state signal
    ;; A[1] - T1
    ;; A[2] - TR (fixed)
    
    return, total( ( (Y - abs(A[0]*( 1.0 - 2.0*exp(-TI/A[1]) + exp(-A[2]/A[1]) )))/err )^2 )

end

function CF2_T1_InvRecovFunction_SE_ABS, TI, A

    ;; A[0] - steady state signal
    ;; A[1] - T1
    ;; A[2] - TR (fixed)
    
    return, abs(A[0]*( 1.0 - 2.0*exp(-TI/A[1]) + exp(-A[2]/A[1]) )) 

end

function CF2_T1_InvRecovFunction_SE, TI, A

    ;; A[0] - steady state signal
    ;; A[1] - T1
    ;; A[2] - TR (fixed)
    
    return, A[0]*( 1.0 - 2.0*exp(-TI/A[1]) + exp(-A[2]/A[1]) ) 

end

function CF2_T2_wBaseline_MultiEchoFunction_SE, TE, A

    ;; A[0] - M0
    ;; A[1] - T2
    ;; A[2] - C (baseline)
    
    return, (A[0] * exp(-TE/A[1])) + A[2]

end

function CF2_T2_MultiEchoFunction_SE, TE, A

    ;; A[0] - M0
    ;; A[1] - T2
    
    return, A[0] * exp(-TE/A[1]) 

end

function CF2_T1_SatRecovFunction_SE, TR, A

    ;; A[0] - Steady state signal
    ;; A[1] - T1
    ;; A[2] - TE (fixed)
    
    return, A[0]*(1-2*exp(-(TR-A[2]/2)/A[1])+exp(-TR/A[1]))
    
end

function CF2_ADC_MultiBVAL_SE, B, A

    ;; A[0] - S0
    ;; A[1] - Diffusion Coeff.
    
    return, A[0] * exp(-B*A[1])

end

;+
; :Description:
;    Model to perform T1 Inversion Recovery fitting for imagery using
;    POWELLs method of optimization.
;
;
; :Keywords:
;    fit_result    - returned, contains the results of the fit
;    threshold_pct - the percent image value threshold 
;    for_slice     - reconstruct a slice that is not currently selected
;    progressbar   - set this to display a progress bar during the fit
;
; :Author: wtriplett
; 
; Modifications (mkulam)
; Added a feature which allows to select which arrays to include in the fit
pro mas_cfit2_IMG_FIT_T1_VarTI_POWELL, fit_result=fit_result, $
                                       threshold_pct=threshold_pct, $
                                       for_slice=for_slice, $
                                       progressbar=progressbar, mask = mask, ci = ci

    common scan_data, project
    if ~n_elements(ci) then ci = project.ci

    ptr_free, project.dataArray[ci].imgFit1
    project.procPramArray[CI].image_fit_flag = 0
    
    TI = project.imndarray[ci].inversion_time_ptr
    if (not ptr_valid(TI)) then begin
        message, 'This model requires an inversion time array.'
    endif

    names = [ 'M0', 'T1', 'TR' ]
    units = [ '', 'sec', 'sec' ]
    fixed = [ 0, 0, 1 ]
    func  = 'CF2_T1_InvRecovFunction_SE_ABS'
    X     = *TI * 1e-3
    if ptr_valid(project.procpramArray[ci].array_select) then X = X[*project.procpramArray[ci].array_select]    
    TR    = project.imndarray[ci].recov_time * 1e-3
    
    parinfo = mas_cfit2_make_basic_parinfo([0, 0, TR], fixed, names, units) 
    parinfo[1].limited = [1, 0]
    parinfo[1].limits = [0, 1e4]
    
    if (n_elements(for_slice) eq 0) then begin
        slice   = abs(mas_cfit2_get_slice_data())
    endif else begin
        slice   = abs(mas_cfit2_get_slice_data(for_slice,ci = ci))
    endelse
    
    dims = size(slice, /dimensions)
    start_map = fltarr(dims[0], dims[1], 3)
    start_map[*,*,0] = max(slice,dimension=3)
    start_map[*,*,1] = mean(X)
    start_map[*,*,2] = TR

    threshold_pct = n_elements(threshold_pct) ? threshold_pct : 0.001
    if isa(mask) then fit_result = mas_cfit2_run_image_fit(func, X, slice, start_map, err, parinfo, $
                                         thr_pct=threshold_pct, progressbar=progressbar, /use_powell, $
                                         powell_objective='CF2_T1_InvRecovFunction_SE_ABS_POWELL', mask = mask) else $
                      fit_result = mas_cfit2_run_image_fit(func, X, slice, start_map, err, parinfo, $
                                         thr_pct=threshold_pct, progressbar=progressbar, /use_powell, $
                                         powell_objective='CF2_T1_InvRecovFunction_SE_ABS_POWELL')                                         
    if (not ptr_valid(fit_result)) then return
    
    p_image = ptr_new(fltarr([2, dims]))
    (*p_image)[0,*,*] = (*fit_result).param_maps[*,*,0]
    (*p_image)[1,*,*] = (*fit_result).param_maps[*,*,1]
    project.dataArray[ci].imgFit1 = p_image
    project.procPramArray[CI].image_fit_flag = 1

end

;+
; :Description:
;    Model to perform T1 fitting by inversion recovery to the signal
;    from an ROI. This uses the POWELL optimization method, and is generally
;    suitable for fitting the IR curve to magnitude mode data, which has a
;    discontinuous derivative at the zero crossing.
;
; :Params:
;    roi_data  - a data structure obtained from the mas_cfit2_roi_get_data() function
;
; :Keywords:
;    fit_result    - returned fit result
;
; :Author: wtriplett
; 
; Modifications (mkulam)
; Added a feature which allows to select which arrays to include in the fit
pro mas_cfit2_ROI_FIT_T1_VarTI_POWELL, roi_data, fit_result=fit_result

    common scan_data, project
    ci = project.ci
    
    TI = project.imndarray[ci].inversion_time_ptr
    if (not ptr_valid(TI)) then begin
        message, 'This model requires an inversion time array.'
    endif
    
    num_rois = n_elements(roi_data)

    names = [ 'M0', 'T1', 'TR' ]
    units = [ '', 'sec', 'sec' ]
    fixed = [ 0, 0, 1 ]
    func  = 'CF2_T1_InvRecovFunction_SE_ABS'
    X     = *TI * 1e-3   
    if ptr_valid(project.procpramArray[ci].array_select) then X = X[*project.procpramArray[ci].array_select]    
    X_sorted = sort(x)
    X = X[X_sorted]
    
    if (n_elements(fit_result) ne 0) then temp = temporary(fit_result)
    
    for roi = 0, num_rois-1 do begin
    
        Y = abs(roi_data[roi].measurement)
        Y = Y[X_sorted]
                
        ;; Error estimates applied to POWELL method not completely understood.
        err   = replicate(1.0, n_elements(X)) ;  roi_data[roi].data_stdev
        start = [ max(Y), mean(X), project.imndarray[project.ci].recov_time * 1e-3 ]
        
        parinfo = mas_cfit2_make_basic_parinfo(start, fixed, names, units)
    
        fit_res_temp = mas_cfit2_run_fit(func, X, Y, err, parinfo, /use_powell, $
                                         powell_objective='CF2_T1_InvRecovFunction_SE_ABS_POWELL')
                                        
        (*fit_res_temp).title  = 'T1 (POWELL) fit to variable TI'
        (*fit_res_temp).xtitle = 'TI (sec.)'
        (*fit_res_temp).ytitle  = 'Measured data'
        
        if (n_elements(fit_result) eq 0) then begin
            fit_result = fit_res_temp
        endif else begin
            fit_result = [ fit_result, fit_res_temp ]
        endelse
        
    endfor
    
end

;+
; :Description:
;    Model to perform T1 Inversion Recovery fitting for imagery using
;    a method to estimate the polarity of each data point before
;    the fitting.
;    
;    1. Points that can be, are estimated by determining if the data values are
;       increasing or decreasing.
;    2. Ambiguous points (should only be one) are thrown out and a preliminary
;       fit is run to determine just the sign of the thrown-out points.
;    3. Sign is corrected and full data set are fit
;
;
; :Keywords:
;    fit_result    - returned, contains the results of the fit
;    threshold_pct - the percent image value threshold 
;    for_slice     - reconstruct a slice that is not currently selected
;    progressbar   - set this to display a progress bar during the fit
;
; :Author: wtriplett
; 
; Modifications (mkulam)
; Added a feature which allows to select which arrays to include in the fit
pro mas_cfit2_IMG_FIT_T1_VarTI_ABS, fit_result=fit_result, $
                                    threshold_pct=threshold_pct, $
                                    for_slice=for_slice, $
                                    progressbar=progressbar

    common scan_data, project
    ci = project.ci

    ptr_free, project.dataArray[ci].imgFit1
    project.procPramArray[CI].image_fit_flag = 0
    
    TI = project.imndarray[ci].inversion_time_ptr
    if (not ptr_valid(TI)) then begin
        message, 'This model requires an inversion time array.'
    endif

    names = [ 'M0', 'T1', 'TR' ]
    units = [ '', 'sec', 'sec' ]
    fixed = [ 0, 0, 1 ]
    func  = 'CF2_T1_InvRecovFunction_SE'
    X     = *TI * 1e-3
    if ptr_valid(project.procpramArray[ci].array_select) then X = X[*project.procpramArray[ci].array_select]       
    TR    = project.imndarray[project.ci].recov_time * 1e-3
    
    parinfo = mas_cfit2_make_basic_parinfo([0, 0, TR], fixed, names, units) 
    parinfo[1].limited = [1, 0]
    parinfo[1].limits = [0, 1e4]
    
    if (n_elements(for_slice) eq 0) then begin
        slice   = abs(mas_cfit2_get_slice_data())
    endif else begin
        slice   = abs(mas_cfit2_get_slice_data(for_slice))
    endelse
    
    dims = size(slice, /dimensions)
    start_map = fltarr(dims[0], dims[1], 3)
    start_map[*,*,0] = max(slice,dimension=3)
    start_map[*,*,1] = mean(X)
    start_map[*,*,2] = TR

    threshold_pct = n_elements(threshold_pct) ? threshold_pct : 0.001
    thr_max = max(abs(slice))
    thr_val = thr_max * threshold_pct
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    if (1 or keyword_set(progressbar)) then begin
        pbar = obj_new('progressbar', /nocancel, $
                       title='Polarity Correction...', $
                       text='Polarity correcting slice data')
        pbar->Start
    endif
    message, /info, 'Polarity correcting magnitude slice data...'
    for xi = 0, dims[0]-1 do begin
        for yi = 0, dims[1]-2 do begin

            y = reform(slice[xi, yi, *])
            
            if (max(abs(y)) lt thr_val) then begin
                continue
            endif
            
            y_abs_min  = min(abs(y), y_abs_min_ind)
            y_diff_map = (y[1:*] - y[0:*])
            
            ;; find the pivot where the signal crosses x axis. this is the ambiguous point
            y_decrease_test = where(y_diff_map lt 0, n_dec, complement=y_increase_test, ncomplement=n_inc)
            if (n_dec ne 0) then begin
                pivot = y_increase_test[0]
            endif else begin
                ;; if for some there is no change in direction, use the abs(min()) point.
                pivot = y_abs_min_ind
            endelse
            
            ;; throw out the pivot 
            lmap = where(lindgen(n_elements(y)) ne pivot)
            dummy_x = x[lmap]
            dummy_y = y[lmap]
            
            ;; flip the sign of the decreasing points
            if (n_dec gt 0) then  dummy_y[y_decrease_test] *= -1
            err   = replicate(1.0, n_elements(dummy_x))
            start = [ max(dummy_y), mean(dummy_x), project.imndarray[project.ci].recov_time * 1e-3 ]
            parinfo = mas_cfit2_make_basic_parinfo(start, fixed, names, units)
            
            ;; fit the curve without the pivot
            fit_res_temp = mas_cfit2_run_fit(func, dummy_x, dummy_y, err, parinfo)
            temp = CF2_T1_InvRecovFunction_SE(X, (*fit_res_temp).params)
            temp = call_function(func,X, (*fit_res_temp).params)
            
            ;; obtain the sign of the pivot from the preliminary fit.
            sign_map = float(fix(temp/abs(temp)))
            y *= sign_map
            slice[xi, yi, *] = y
        endfor
        if (obj_valid(pbar)) then pbar->Update, float(xi)/dims[0] * 100
    endfor
    if (obj_valid(pbar)) then pbar->Destroy
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    err   = replicate(1.0, n_elements(X))
    fit_result = mas_cfit2_run_image_fit(func, X, slice, start_map, err, parinfo, $
                                         thr_pct=threshold_pct, progressbar=progressbar)
                                         
    if (not ptr_valid(fit_result)) then return
    
    p_image = ptr_new(fltarr([2, dims]))
    (*p_image)[0,*,*] = (*fit_result).param_maps[*,*,0]
    (*p_image)[1,*,*] = (*fit_result).param_maps[*,*,1]
    project.dataArray[ci].imgFit1 = p_image
    project.procPramArray[CI].image_fit_flag = 1

end

;+
; :Description:
;    Model to perform T1 Inversion Recovery fitting for ROI using
;    a method to estimate the polarity of each data point before
;    the fitting.
;    
;    1. Points that can be, are estimated by determining if the data values are
;       increasing or decreasing.
;    2. Ambiguous points (should only be one) are thrown out and a preliminary
;       fit is run to determine just the sign of the thrown-out points.
;    3. Sign is corrected and full data set are fit
;
;
; :Keywords:
;    fit_result    - returned, contains the results of the fit
;    threshold_pct - the percent image value threshold 
;    for_slice     - reconstruct a slice that is not currently selected
;    progressbar   - set this to display a progress bar during the fit
;
; :Author: wtriplett
; 
; Modifications (mkulam)
; Added a feature which allows to select which arrays to include in the fit
pro mas_cfit2_ROI_FIT_T1_VarTI_ABS, roi_data, fit_result=fit_result

    common scan_data, project
    ci = project.ci
    
    TI = project.imndarray[ci].inversion_time_ptr
    if (not ptr_valid(TI)) then begin
        message, 'This model requires an inversion time array.'
    endif
    
    num_rois = n_elements(roi_data)

    names = [ 'M0', 'T1', 'TR' ]
    units = [ '', 'sec', 'sec' ]
    fixed = [ 0, 0, 1 ]
    func = 'CF2_T1_InvRecovFunction_SE'
    X     = *TI * 1e-3   
    if ptr_valid(project.procpramArray[ci].array_select) then X = X[*project.procpramArray[ci].array_select]        
    X_sorted = sort(x)
    X = X[X_sorted]
    
    if (n_elements(fit_result) ne 0) then temp = temporary(fit_result)
    
    for roi = 0, num_rois-1 do begin
    
        Y = abs(roi_data[roi].measurement)
        Y = Y[X_sorted]
        
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        y_abs_min  = min(abs(y), y_abs_min_ind)
        y_diff_map = (y[1:*] - y[0:*])
        
        ;; find the pivot where the signal crosses x axis. this is the ambiguous point
        y_decrease_test = where(y_diff_map lt 0, n_dec, complement=y_increase_test, ncomplement=n_inc)
        if (n_dec ne 0) then begin
            pivot = y_increase_test[0]
        endif else begin
            ;; if for some there is no change in direction, use the abs(min()) point.
            pivot = y_abs_min_ind
        endelse
        
        ;; throw out the pivot 
        lmap = where(lindgen(n_elements(y)) ne pivot)
        dummy_x = x[lmap]
        dummy_y = y[lmap]
        
        ;; flip the sign of the decreasing points
        if (n_dec gt 0) then  dummy_y[y_decrease_test] *= -1
        err   = replicate(1.0, n_elements(dummy_x))
        start = [ max(dummy_y), mean(dummy_x), project.imndarray[project.ci].recov_time * 1e-3 ]
        parinfo = mas_cfit2_make_basic_parinfo(start, fixed, names, units)
        
        ;; fit the curve without the pivot
        fit_res_temp = mas_cfit2_run_fit(func, dummy_x, dummy_y, err, parinfo)
        temp = CF2_T1_InvRecovFunction_SE(X, (*fit_res_temp).params)
        temp = call_function(func,X, (*fit_res_temp).params)
        
        ;; obtain the sign of the pivot from the preliminary fit.
        sign_map = float(fix(temp/abs(temp)))
        y *= sign_map
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        
        ;; Error estimates applied to POWELL method not completely understood.
        err   = roi_data[roi].data_stdev
        start = [ max(Y), mean(X), project.imndarray[project.ci].recov_time * 1e-3 ]
        
        parinfo = mas_cfit2_make_basic_parinfo(start, fixed, names, units)
        
        fit_res_temp = mas_cfit2_run_fit(func, X, Y, err, parinfo)                                         
        (*fit_res_temp).title  = 'T1 (ABS) fit to variable TI'
        (*fit_res_temp).xtitle = 'TI (sec.)'
        (*fit_res_temp).ytitle  = 'Measured data'
        
        if (n_elements(fit_result) eq 0) then begin
            fit_result = fit_res_temp
        endif else begin
            fit_result = [ fit_result, fit_res_temp ]
        endelse
        
    endfor
    
end


;+
; :Description:
;    Model to perform T1 Inversion Recovery fitting for imagery using
;    MPFIT optimization.
;
;
; :Keywords:
;    fit_result    - returned, contains the results of the fit
;    threshold_pct - the percent image value threshold 
;    for_slice     - reconstruct a slice that is not currently selected
;    progressbar   - set this to display a progress bar during the fit
;
; :Author: wtriplett
; 
; Modifications (mkulam)
; Added a feature which allows to select which arrays to include in the fit
pro mas_cfit2_IMG_FIT_T1_VarTI, fit_result=fit_result, $
                                threshold_pct=threshold_pct, $
                                for_slice=for_slice, $
                                progressbar=progressbar

    common scan_data, project
    ci = project.ci

    ptr_free, project.dataArray[ci].imgFit1
    project.procPramArray[CI].image_fit_flag = 0
    
    TI = project.imndarray[ci].inversion_time_ptr
    if (not ptr_valid(TI)) then begin
        message, 'This model requires an inversion time array.'
    endif

    names = [ 'M0', 'T1', 'TR' ]
    units = [ '', 'sec', 'sec' ]
    fixed = [ 0, 0, 1 ]
    func  = 'CF2_T1_InvRecovFunction_SE'
    X     = *TI * 1e-3
    if ptr_valid(project.procpramArray[ci].array_select) then X = X[*project.procpramArray[ci].array_select]        
    TR    = project.imndarray[project.ci].recov_time * 1e-3
    
    parinfo = mas_cfit2_make_basic_parinfo([0, 0, TR], fixed, names, units) 
;    parinfo[1].limited = [1, 0]
;    parinfo[1].limits = [0, 1e4]
    
    if (n_elements(for_slice) eq 0) then begin
        slice   = real_part(mas_cfit2_get_slice_data())
    endif else begin
        slice   = real_part(mas_cfit2_get_slice_data(for_slice))
    endelse

    dims = size(slice, /dimensions)
    start_map = fltarr(dims[0], dims[1], 3)
    start_map[*,*,0] = max(slice,dimension=3)
    start_map[*,*,1] = mean(X)
    start_map[*,*,2] = TR

    threshold_pct = n_elements(threshold_pct) ? threshold_pct : 0.001
    fit_result = mas_cfit2_run_image_fit(func, X, slice, start_map, err, parinfo, $
                                         thr_pct=threshold_pct, progressbar=progressbar)
    if (not ptr_valid(fit_result)) then return
    
    p_image = ptr_new(fltarr([2, dims]))
    (*p_image)[0,*,*] = (*fit_result).param_maps[*,*,0]
    (*p_image)[1,*,*] = (*fit_result).param_maps[*,*,1]
    project.dataArray[ci].imgFit1 = p_image
    project.procPramArray[CI].image_fit_flag = 1

end


;+
; :Description:
;    Inversion recovery T1 fit to ROI data..
;
; :Params:
;    roi_data
;
; :Keywords:
;    fit_result
;
; :Author: wtriplett
; 
; Modifications (mkulam)
; Added a feature which allows to select which arrays to include in the fit
pro mas_cfit2_ROI_FIT_T1_VarTI, roi_data, fit_result=fit_result

    common scan_data, project
    ci = project.ci
    
    TI = project.imndarray[ci].inversion_time_ptr
    if (not ptr_valid(TI)) then begin
        message, 'This model requires an inversion time array.'
    endif
    
    num_rois = n_elements(roi_data)

    names = [ 'M0', 'T1', 'TR' ]
    units = [ '', 'sec', 'sec' ]
    fixed = [ 0, 0, 1 ]
    func  = 'CF2_T1_InvRecovFunction_SE'
    X     = *TI * 1e-3   
    if ptr_valid(project.procpramArray[ci].array_select) then X = X[*project.procpramArray[ci].array_select]    
            
    if (n_elements(fit_result) ne 0) then temp = temporary(fit_result)
    
    for roi = 0, num_rois-1 do begin
    
        Y     = roi_data[roi].measurement
        err   = roi_data[roi].data_stdev
        start = [ max(Y), mean(X), project.imndarray[project.ci].recov_time * 1e-3 ]
        
        parinfo = mas_cfit2_make_basic_parinfo(start, fixed, names, units)
    
        fit_res_temp = mas_cfit2_run_fit(func, X, Y, err, parinfo)
        (*fit_res_temp).title  = 'T1 fit to variable TI'
        (*fit_res_temp).xtitle = 'TI (sec.)'
        (*fit_res_temp).ytitle  = 'Measured data'
        
        if (n_elements(fit_result) eq 0) then begin
            fit_result = fit_res_temp
        endif else begin
            fit_result = [ fit_result, fit_res_temp ]
        endelse
        
    endfor
    
end

;+
; :Description:
;    Model to perform T1 fitting for imagery using arrayed TR.
;
; :Keywords:
;    fit_result    - returned, contains the results of the fit
;    threshold_pct - the percent image value threshold 
;    for_slice     - reconstruct a slice that is not currently selected
;    progressbar   - set this to display a progress bar during the fit
;
; :Author: wtriplett
; 
; Modifications (Magdoom)
; 1) Fixed bug (Included a piece of code to save the fits, useful for proton density maps)
; 2) Added a feature which allows to select which arrays to include in the fit
; 3) Added functionality to perform the fit in selected voxels given by mask = 1

pro mas_cfit2_IMG_FIT_T1_VarTR, fit_result=fit_result, $
                                threshold_pct=threshold_pct, $
                                for_slice=for_slice, $
                                progressbar=progressbar, mask = mask, ci = ci

    common scan_data, project
    
    if ~n_elements(ci) then ci = project.ci
    
    TR = project.imndarray[ci].rep_time_ptr
    if (not ptr_valid(TR)) then begin
        message, 'This model requires an repetition time array.'
    endif

    names = [ 'M0', 'T1', 'TE' ]
    units = [ '', 'sec', 'sec' ]
    fixed = [ 0, 0, 1 ]
    func  = 'CF2_T1_SatRecovFunction_SE'
    X     = (*TR) * 1e-3
    if ptr_valid(project.procpramArray[ci].array_select) then X = X[*project.procpramArray[ci].array_select]
        
    TE    = project.imndarray[CI].echo_time * 1e-3
    
    parinfo = mas_cfit2_make_basic_parinfo([0, 0, TE], fixed, names, units) 

    if (n_elements(for_slice) eq 0) then begin
        slice   = real_part(mas_cfit2_get_slice_data(ci = ci))
    endif else begin
        slice   = real_part(mas_cfit2_get_slice_data(for_slice,ci = ci))
    endelse

    dims = size(slice, /dimensions)
    start_map = fltarr(dims[0], dims[1], 3)
    start_map[*,*,0] = max(slice,dimension=3)
    start_map[*,*,1] = mean(X)
    start_map[*,*,2] = TE

    threshold_pct = n_elements(threshold_pct) ? threshold_pct : 0.001
    if isa(mask) then fit_result = mas_cfit2_run_image_fit(func, X, slice, start_map, err, parinfo, $
                                         thr_pct=threshold_pct, progressbar=progressbar,mask=mask) else $
                      fit_result = mas_cfit2_run_image_fit(func, X, slice, start_map, err, parinfo, $
                                         thr_pct=threshold_pct, progressbar=progressbar)                   
    if (not ptr_valid(fit_result)) then return
    
    p_image = ptr_new(fltarr([2, dims]))
    (*p_image)[0,*,*] = (*fit_result).param_maps[*,*,0]
    (*p_image)[1,*,*] = (*fit_result).param_maps[*,*,1]
    project.dataArray[ci].imgFit1 = p_image
    project.procPramArray[CI].image_fit_flag = 1
    
end


;+
; :Description:
;    Arrayed TR,  T1 fit to ROI data..
;
; :Params:
;    roi_data
;
; :Keywords:
;    fit_result
;
; :Author: wtriplett
; 
; Modifications (mkulam)
; Added a feature which allows to select which arrays to include in the fit
pro mas_cfit2_ROI_FIT_T1_VarTR, roi_data, fit_result=fit_result, $
                                signal_method=signal_method

    common scan_data, project
    ci = project.ci
    
    TR = project.imndarray[ci].rep_time_ptr
    if (not ptr_valid(TR)) then begin
        message, 'This model requires an repetition time array.'
    endif

    if (not ptr_valid(TR)) then begin
        message, 'Data set has invalid TR values.'
    endif
    
    num_rois = n_elements(roi_data)
    
    TE = project.imndarray[project.ci].echo_time * 1e-3 

    names = [ 'M0', 'T1', 'TE' ]
    units = [ '', 'sec', 'sec' ]
    fixed = [ 0, 0, 1 ]
    func  = 'CF2_T1_SatRecovFunction_SE'
    X = *TR * 1e-3
    if ptr_valid(project.procpramArray[ci].array_select) then X = X[*project.procpramArray[ci].array_select]
               
    if (n_elements(fit_result) ne 0) then temp = temporary(fit_result)
    
    for roi = 0, num_rois-1 do begin
    
        Y     = roi_data[roi].measurement
        err   = roi_data[roi].data_stdev
        start = [ max(Y), mean(X), TE ]
        
        parinfo = mas_cfit2_make_basic_parinfo(start, fixed, names, units)
        
        fit_res_temp = mas_cfit2_run_fit(func, X, Y, err, parinfo)
        (*fit_res_temp).title  = 'T1 fit to variable TR'
        (*fit_res_temp).xtitle = 'TR (sec.)'
        (*fit_res_temp).ytitle  = 'Measured data'
        
        if (n_elements(fit_result) eq 0) then begin
            fit_result = fit_res_temp
        endif else begin
            fit_result = [ fit_result, fit_res_temp ]
        endelse
        
    endfor
    
end

;+
; :Description:
;    Model to perform T2 fitting for imagery using from arrayed TEs.
;
;
; :Keywords:
;    fit_result    - returned, contains the results of the fit
;    threshold_pct - the percent image value threshold 
;    for_slice     - reconstruct a slice that is not currently selected
;    progressbar   - set this to display a progress bar during the fit
;
; :Author: wtriplett
; 
; Modifications (mkulam)
; Added a feature which allows to select which arrays to include in the fit
pro mas_cfit2_IMG_FIT_T2_VarTE, fit_result=fit_result, $
                                threshold_pct=threshold_pct, $
                                for_slice=for_slice, $
                                progressbar=progressbar

    common scan_data, project
    ci = project.ci
    
    TE = project.imndarray[ci].echo_time_ptr
    if (not ptr_valid(TE)) then begin
        message, 'This model requires an echo time array.'
    endif

    names = [ 'M0', 'T2' ]
    units = [ '', 'sec' ]
    fixed = [ 0, 0 ]
    func  = 'CF2_T2_MultiEchoFunction_SE'
    X     = (*TE) * 1e-3
    if ptr_valid(project.procpramArray[ci].array_select) then X = X[*project.procpramArray[ci].array_select]    
    
    parinfo = mas_cfit2_make_basic_parinfo([0, 0], fixed, names, units)
    if (n_elements(for_slice) eq 0) then begin
        slice   = real_part(mas_cfit2_get_slice_data())
    endif else begin
        slice   = real_part(mas_cfit2_get_slice_data(for_slice))
    endelse

    dims = size(slice, /dimensions)
    start_map = fltarr(dims[0], dims[1], 3)
    start_map[*,*,0] = max(slice,dimension=3)
    start_map[*,*,1] = mean(X)

    threshold_pct = n_elements(threshold_pct) ? threshold_pct : 0.001
    fit_result = mas_cfit2_run_image_fit(func, X, slice, start_map, err, parinfo, $
                                         thr_pct=threshold_pct, progressbar=progressbar)
    if (not ptr_valid(fit_result)) then return
    
    p_image = ptr_new(fltarr([2, dims]))
    (*p_image)[0,*,*] = (*fit_result).param_maps[*,*,0]
    (*p_image)[1,*,*] = (*fit_result).param_maps[*,*,1]
    project.dataArray[ci].imgFit1 = p_image
    project.procPramArray[CI].image_fit_flag = 1
    
end

;+
; :Description:
;    T2 ROI fit to multiple echo times
;
; :Params:
;    roi_data
;
; :Keywords:
;    fit_result
;
; :Author: wtriplett
; 
; Modifications (mkulam)
; Added a feature which allows to select which arrays to include in the fit
pro mas_cfit2_ROI_FIT_T2_VarTE, roi_data, fit_result=fit_result

    common scan_data, project
    ci = project.ci
        
    TE = project.imndarray[ci].echo_time_ptr
    if (not ptr_valid(TE)) then begin
        message, 'This model requires an echo time array.'
    endif

    num_rois = n_elements(roi_data)

    names = [ 'M0', 'T2' ]
    units = [ '', 'sec' ]
    fixed = [ 0, 0 ]
    func  = 'CF2_T2_MultiEchoFunction_SE'
    X     = *TE * 1e-3
    if ptr_valid(project.procpramArray[ci].array_select) then X = X[*project.procpramArray[ci].array_select]    
            
    if (n_elements(fit_result) ne 0) then temp = temporary(fit_result)
    
    for roi = 0, num_rois-1 do begin
    
        Y     = roi_data[roi].measurement
        err   = roi_data[roi].data_stdev
        start = [ max(Y), mean(X) ]
        
        parinfo    = mas_cfit2_make_basic_parinfo(start, fixed, names, units)

        fit_res_temp = mas_cfit2_run_fit(func, X, Y, err, parinfo)
        (*fit_res_temp).title  = 'T2 fit to variable TE'
        (*fit_res_temp).xtitle = 'TE (sec.)'
        (*fit_res_temp).ytitle  = 'Measured data'
        
        if (n_elements(fit_result) eq 0) then begin
            fit_result = fit_res_temp
        endif else begin
            fit_result = [ fit_result, fit_res_temp ]
        endelse

    endfor
        
end

;+
; :Description:
;    Model to perform T2 fitting for imagery using arrayed TEs and
;    also fit a baseline or noise floor
;
; :Keywords:
;    fit_result    - returned, contains the results of the fit
;    threshold_pct - the percent image value threshold 
;    for_slice     - reconstruct a slice that is not currently selected
;    progressbar   - set this to display a progress bar during the fit
;
; :Author: wtriplett
; 
; Modifications (mkulam)
; Added a feature which allows to select which arrays to include in the fit
pro mas_cfit2_IMG_FIT_T2_wBaseline_VarTE, fit_result=fit_result, $
                                          threshold_pct=threshold_pct, $
                                          for_slice=for_slice, $
                                          progressbar=progressbar

    common scan_data, project
    ci = project.ci
    
    TE = project.imndarray[ci].echo_time_ptr
    if (not ptr_valid(TE)) then begin
        message, 'This model requires an echo time array.'
    endif

    names = [ 'M0', 'T2', 'C' ]
    units = [ '', 'sec', '' ]
    fixed = [ 0, 0, 0 ]
    func  = 'CF2_T2_wBaseline_MultiEchoFunction_SE'
    X     = *TE * 1e-3
    if ptr_valid(project.procpramArray[ci].array_select) then X = X[*project.procpramArray[ci].array_select]    
    
    parinfo = mas_cfit2_make_basic_parinfo([0, 0, 0], fixed, names, units) 

    if (n_elements(for_slice) eq 0) then begin
        slice   = real_part(mas_cfit2_get_slice_data())
    endif else begin
        slice   = real_part(mas_cfit2_get_slice_data(for_slice))
    endelse

    dims = size(slice, /dimensions)
    start_map = fltarr(dims[0], dims[1], 3)
    start_map[*,*,0] = max(slice,dimension=3)
    start_map[*,*,1] = mean(X)
    start_map[*,*,0] = 0.0
    
    threshold_pct = n_elements(threshold_pct) ? threshold_pct : 0.001
    fit_result = mas_cfit2_run_image_fit(func, X, slice, start_map, err, parinfo, $
                                         thr_pct=threshold_pct, progressbar=progressbar)
    if (not ptr_valid(fit_result)) then return

    p_image = ptr_new(fltarr([2, dims]))
    (*p_image)[0,*,*] = (*fit_result).param_maps[*,*,0]
    (*p_image)[1,*,*] = (*fit_result).param_maps[*,*,1]
    project.dataArray[ci].imgFit1 = p_image
    project.procPramArray[CI].image_fit_flag = 1
    
end

;+
; :Description:
;    T2 fit to multiple TEs with a baseline or noise floor..
;
; :Params:
;    roi_data
;
; :Keywords:
;    fit_result
;
; :Author: wtriplett
; 
; Modifications (mkulam)
; Added a feature which allows to select which arrays to include in the fit
pro mas_cfit2_ROI_FIT_T2_wBaseline_VarTE, roi_data, fit_result=fit_result

    common scan_data, project
    ci = project.ci
        
    TE = project.imndarray[ci].echo_time_ptr
    if (not ptr_valid(TE)) then begin
        message, 'This model requires an echo time array.'
    endif

    num_rois = n_elements(roi_data)

    names = [ 'M0', 'T2', 'C' ]
    units = [ '', 'sec', '' ]
    fixed = [ 0, 0, 0 ]
    func  = 'CF2_T2_wBaseline_MultiEchoFunction_SE'
    X     = *TE * 1e-3
    if ptr_valid(project.procpramArray[ci].array_select) then X = X[*project.procpramArray[ci].array_select]    
            
    if (n_elements(fit_result) ne 0) then temp = temporary(fit_result)
    
    for roi = 0, num_rois-1 do begin
    
        Y     = roi_data[roi].measurement
        err   = roi_data[roi].data_stdev
        start = [ max(Y), mean(X), 0.0 ]
        
        parinfo    = mas_cfit2_make_basic_parinfo(start, fixed, names, units)

        fit_res_temp = mas_cfit2_run_fit(func, X, Y, err, parinfo)
        (*fit_res_temp).title  = 'T2 fit to variable TE (w/ Baseline)'
        (*fit_res_temp).xtitle = 'TE (sec.)'
        (*fit_res_temp).ytitle  = 'Measured data'
        
        if (n_elements(fit_result) eq 0) then begin
            fit_result = fit_res_temp
        endif else begin
            fit_result = [ fit_result, fit_res_temp ]
        endelse

    endfor
        
end

;+
; :Description:
;    Model to perform ADC fitting for imagery using arrayed b-values
;
; :Keywords:
;    fit_result    - returned, contains the results of the fit
;    threshold_pct - the percent image value threshold 
;    for_slice     - reconstruct a slice that is not currently selected
;    progressbar   - set this to display a progress bar during the fit
;
; :Author: wtriplett
; 
; Modifications (mkulam)
; 1) Fixed ADC unit display
; 2) Added a feature which allows to select which arrays to include in the fit
pro mas_cfit2_IMG_FIT_ADC_VarB, fit_result=fit_result, $
                                threshold_pct=threshold_pct, $
                                for_slice=for_slice, $
                                progressbar=progressbar

    common scan_data, project
    ci = project.ci
    
    B = project.imndarray[ci].bval_array
    if (not ptr_valid(B)) then begin
        message, 'This model requires an b-value array.'
    endif

    names = [ 'S0', 'D' ]
    units = [ '', 'um^2/ms' ]
    fixed = [ 0, 0 ]
    func  = 'CF2_ADC_MultiBVAL_SE'
    X     = *B * 1e-3
    if ptr_valid(project.procpramArray[ci].array_select) then X = X[*project.procpramArray[ci].array_select]    
    
    parinfo = mas_cfit2_make_basic_parinfo([0, 0], fixed, names, units) 
    if (n_elements(for_slice) eq 0) then begin
        slice   = real_part(mas_cfit2_get_slice_data())
    endif else begin
        slice   = real_part(mas_cfit2_get_slice_data(for_slice))
    endelse

    dims = size(slice, /dimensions)
    start_map = fltarr(dims[0], dims[1], 3)
    start_map[*,*,0] = max(slice,dimension=3)
    start_map[*,*,1] = mean(X)

    threshold_pct = n_elements(threshold_pct) ? threshold_pct : 0.001
    fit_result = mas_cfit2_run_image_fit(func, X, slice, start_map, err, parinfo, $
                                         thr_pct=threshold_pct, progressbar=progressbar)
    if (not ptr_valid(fit_result)) then return
    
    p_image = ptr_new(fltarr([2, dims]))
    (*p_image)[0,*,*] = (*fit_result).param_maps[*,*,0]
    (*p_image)[1,*,*] = (*fit_result).param_maps[*,*,1]
    project.dataArray[ci].imgFit1 = p_image
    project.procPramArray[CI].image_fit_flag = 1
    
end

;+
; :Description:
;    ADC fit to multiple b-values
;
; :Params:
;    roi_data
;
; :Keywords:
;    fit_result
;
; :Author: wtriplett
; 
; Modifications (mkulam)
; 1) Fixed ADC unit display
; 2) Added a feature which allows to select which arrays to include in the fit
pro mas_cfit2_ROI_FIT_ADC_VarB, roi_data, fit_result=fit_result

    common scan_data, project
    ci = project.ci
        
    B = project.imndarray[ci].bval_array
    if (not ptr_valid(B)) then begin
        message, 'This model requires an b-value array.'
    endif

    num_rois = n_elements(roi_data)
    
    names = [ 'S0', 'D' ]
    units = [ '', 'um^2/ms' ]
    fixed = [ 0, 0 ]
    func  = 'CF2_ADC_MultiBVAL_SE'
    X     = *B * 1e-3
    if ptr_valid(project.procpramArray[ci].array_select) then X = X[*project.procpramArray[ci].array_select]    
    
    if (n_elements(fit_result) ne 0) then temp = temporary(fit_result)
    
    for roi = 0, num_rois-1 do begin
    
        Y     = roi_data[roi].measurement
        err   = roi_data[roi].data_stdev
        start = [ max(Y), mean(X) ]
        
        parinfo = mas_cfit2_make_basic_parinfo(start, fixed, names, units)

        fit_res_temp = mas_cfit2_run_fit(func, X, Y, err, parinfo)
        (*fit_res_temp).title  = 'ADC fit to variable B'
        (*fit_res_temp).xtitle = 'B (ms/um^2)'
        (*fit_res_temp).ytitle  = 'Measured data'
        
        if (n_elements(fit_result) eq 0) then begin
            fit_result = fit_res_temp
        endif else begin
            fit_result = [ fit_result, fit_res_temp ]
        endelse

    endfor
        
end
