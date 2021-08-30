function mas_hardi_tractography_loadPrepDirs, event

    txt = widget_info(event.top, find_by_uname='txt_prep_dir')
    widget_control, txt, get_value=dir, get_uvalue=files
    if (dir eq '') then begin
        return, -1
    endif

    nfiles = n_elements(files)    
    if (nfiles eq 0) then begin
        void = dialog_message(['Count not find any direction files in',$
                              dir, $
                              'Appropriate files should be named: "V_XX.nii",', $
                              'where XX is a number (WITH leading zero.)'], /error, /center)
        return, -1
    endif else begin
        progressBar = obj_new('progressbar', text="Loading Directions...", title="Please Wait...")
        progressBar->Start
        dir_ptrs = ptrarr(nfiles)
        for f = 0, nfiles-1 do begin
            nif = mas_read_nifti(nifti_filename=files[f], read_status=rs)
            if (rs ne 0 and ptr_valid(nif.voxel_data)) then begin
                dir_ptrs[f] = nif.voxel_data
            endif else begin
                print, "Bad nifti file: "+files[f]+' (read_status: '+string(rs, format='(I0)')+')'
            endelse
            progressBar->Update, float(f)/float(nfiles) * 100.0
        endfor
        progressBar->Destroy
        return, dir_ptrs
    endelse
end

pro mas_hardi_tractography_gui_togglevf, event

    ;; desensitize some of the user interface elements depending on whether
    ;; or not we are tracking from preprocessed directions or from live data
    txt_prep_dir    = widget_info(event.top, find_by_uname='txt_prep_dir')
    dl_prob_method  = widget_info(event.top, find_by_uname='dl_prob_method')
    txt_odf_val     = widget_info(event.top, find_by_uname='txt_odf_val')
    btn_odf_prep    = widget_info(event.top, find_by_uname='btn_odf_prep')
    btn_config_prob = widget_info(event.top, find_by_uname='btn_config_prob')
    dl_tess_level   = widget_info(event.top, find_by_uname='dl_tess_level')
    dl_interp_method= widget_info(event.top, find_by_uname='dl_interp_method')
    
    txt_prep_dir = widget_info(event.top, find_by_uname='txt_prep_dir')
    btn_prep_dir = widget_info(event.top, find_by_uname='btn_prep_dir')

    widget_control, dl_tess_level,    sensitive=1-event.select
    widget_control, txt_prep_dir,     sensitive=1-event.select
    widget_control, dl_prob_method,   sensitive=1-event.select
    widget_control, txt_odf_val,      sensitive=1-event.select
    widget_control, btn_odf_prep,     sensitive=1-event.select
    widget_control, btn_config_prob,  sensitive=1-event.select
    if (event.select eq 1) then begin
        old_select = widget_info(dl_interp_method, /droplist_select)
        widget_control, dl_interp_method, set_uvalue=old_select, set_droplist_select=0
    endif else begin
        widget_control, dl_interp_method, get_uvalue=old_select
        widget_control, dl_interp_method, set_droplist_select=old_select
    endelse
    widget_control, dl_interp_method, sensitive=1-event.select
    
    widget_control, txt_prep_dir, sensitive=event.select
    widget_control, btn_prep_dir, sensitive=event.select
    
end

pro mas_hardi_tractography_gui_preprocess_event, event

    forward_function mas_hardi_tractography_gui_getODFObj
    forward_function mas_hardi_tractography_gui_getTessLevel
    
    base = widget_base(title='Preprocessing Parameters', group_leader=event.top, /modal, /column)
    
    b = widget_base(base, /row)
    lbl = widget_label(b, value='Number of directions to resolve:')
    ndirs = widget_droplist(b, value=['1','2','3','4','5'], uname='num_dirs')
    
    b = widget_base(base, /row, /align_center)
    btn_cancel = widget_button(b, value='Cancel', event_pro='mas_hardi_tractography_gui_prep_btn_cancel')
    btn_ok = widget_button(b, value='Start', event_pro='mas_hardi_tractography_gui_prep_btn_ok')
    
    result = ptr_new({num_dirs: 0, dest_dir: ''})
    
    widget_control, base, set_uvalue=result, /realize
    xmanager, 'mas_hardi_tractography_gui', base
    ;; At this point the result structure will contain the user-entered parameters
    
    if ((*result).dest_dir eq '') then begin
        ptr_free, result
        return
    endif else begin
        tess_level = mas_hardi_tractography_gui_getTessLevel(event)
        odf_obj   = mas_hardi_tractography_gui_getODFObj(event)    
        mas_hardi_tractography_preprocess, odf_obj, threshold=1,$
                                           destdir=(*result).dest_dir, $
                                           ndirections=(*result).num_dirs, $
                                           tess_level=tess_level
        ptr_free, result
    endelse 
end

pro mas_hardi_tractography_gui_prep_btn_ok, event

    common scan_data
    widget_control, event.top, get_uvalue=result
    dl_ndirs = widget_info(event.top, find_by_uname='num_dirs')
    (*result).num_dirs = fix(widget_info(dl_ndirs, /droplist_select))+1
    
    dest_dir = dialog_pickfile(path=project.current_path, title='Select a destination directory', $
                               /directory, /write)
    (*result).dest_dir = dest_dir
    
    widget_control, event.top, /destroy
    
end

pro mas_hardi_tractography_gui_prep_btn_cancel, event

    widget_control, event.top, get_uvalue=result
    (*result).dest_dir = ''
    widget_control, event.top, /destroy
    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mas_hardi_tractography_gui_getAnisThr, event

    txt_fa_halt = widget_info(event.top, find_by_uname='txt_fa_halt')
    widget_control, txt_fa_halt, get_value=fa_halt
    return, float(fa_halt[0])
    
end

function mas_hardi_tractography_gui_getProbThr, event

    txt_odf_val = widget_info(event.top, find_by_uname='txt_odf_val')
    widget_control, txt_odf_val, get_value=odf_val
    return, float(odf_val[0])
    
end

function mas_hardi_tractography_gui_getAngleThr, event

    txt_turn_angle = widget_info(event.top, find_by_uname='txt_turn_angle')
    widget_control, txt_turn_angle, get_value=angle_thr
    return, float(angle_thr[0])
    
end

function mas_hardi_tractography_gui_getStepSize, event

    txt_stepsize = widget_info(event.top, find_by_uname='txt_stepsize')
    widget_control, txt_stepsize, get_value=stepsize
    return, float(stepsize[0])
    
end

function mas_hardi_tractography_gui_getTessLevel, event

    dl_tess_level = widget_info(event.top, find_by_uname='dl_tess_level')
    tess_level = widget_info(dl_tess_level, /droplist_select)
    
    return, fix(tess_level)+3
    
end

function mas_hardi_tractography_gui_getInterp, event

    dl_interp = widget_info(event.top, find_by_uname='dl_interp_method')
    interp_meth = widget_info(dl_interp, /droplist_select)
    return, interp_meth[0]
    
end

function mas_hardi_tractography_gui_getSeedAlgorithm, event

    dl_seed_algo  = widget_info(event.top, find_by_uname='dl_seed_algorithm')
    algo_type = widget_info(dl_seed_algo, /DROPLIST_SELECT)
    return, fix(algo_type[0])
    
end

function mas_hardi_tractography_gui_getSeedDensity, event

    txt_seed_dens = widget_info(event.top, find_by_uname='txt_seed_dens')
    widget_control, txt_seed_dens, get_value=seed_density
    return, fix(seed_density[0])

end

function mas_hardi_tractography_gui_isSplittingEnabled, event

    btn_splitting_enabled = widget_info(event.top, find_by_uname='btn_splitting_enabled')
    return, widget_info(btn_splitting_enabled, /button_set)

end

function mas_hardi_tractography_gui_isVectorFieldEnabled, event

    btn_use_vector_field = widget_info(event.top, find_by_uname='btn_use_vector_field')
    return, widget_info(btn_use_vector_field, /button_set)

end

function mas_hardi_tractography_gui_getMaximumSplit, event

    txt_max_split = widget_info(event.top, find_by_uname='txt_max_split')
    widget_control, txt_max_split, get_value=max_split
    return, fix(max_split[0])

end

function mas_hardi_tractography_gui_getFiberFlip, event

    btn_flip_x = widget_info(event.top, find_by_uname='btn_flip_x')
    btn_flip_y = widget_info(event.top, find_by_uname='btn_flip_y')
    btn_flip_z = widget_info(event.top, find_by_uname='btn_flip_z')
    
    flip_x = widget_info(btn_flip_x, /button_set) ? -1 : 1
    flip_y = widget_info(btn_flip_y, /button_set) ? -1 : 1
    flip_z = widget_info(btn_flip_z, /button_set) ? -1 : 1
    
    return,  float([flip_x, flip_y, flip_z, 1.0])

end

function mas_hardi_tractography_gui_getWorkDir, event

    txt_working_dir = widget_info(event.top, find_by_uname='txt_working_dir')
    widget_control, txt_working_dir, get_value=working_dir
    return, strtrim(working_dir, 2)

end

;function mas_hardi_tractography_gui_getOutputFile, event
;
;    txt_trk_output = widget_info(event.top, find_by_uname='txt_trk_output')
;    widget_control, txt_trk_output, get_value=trk_output
;    return, strtrim(trk_output, 2)
;
;end

;function mas_hardi_tractography_gui_getMatrixFile, event
;
;    txt_tx_file = widget_info(event.top, find_by_uname='txt_tx_file')
;    widget_control, txt_tx_file, get_value=matrix_file
;    return, strtrim(matrix_file, 2)
;
;end

function mas_hardi_tractography_gui_getODFObj, event

    common scan_data
    
    dl_odf_select = widget_info(event.top, find_by_uname='dl_prob_method')
    odf_type = widget_info(dl_odf_select, /droplist_select)
    
    case odf_type of
    
        0: begin
            method_state = project.procpramarray[project.ci].odf_state
            
            if (not ptr_valid(method_state)) then begin
                method_state = mas_odf_gui_make_default_state()
                project.procpramarray[project.ci].odf_state = method_state
            endif
            
            mas_odf_init_obj, $
                lmax=(*method_state).odf_lmax, $
                odf_obj=method, /bulk
        end
        
        1: begin
            method_state = project.procpramarray[project.ci].dot_state
            
            if (not ptr_valid(method_state)) then begin
                method_state = mas_dot_gui_make_default_state()
                project.procpramarray[project.ci].dot_state = method_state
            endif
            
            mas_dot_init_obj, $
                lmax=(*method_state).dot_lmax, $
                radius=(*method_state).dot_rzero, $
                dot_obj=method, /bulk
        end
        
        2: begin
            method_state = project.procpramarray[project.ci].mow_state
            
            if (not ptr_valid(method_state)) then begin
                method_state = mas_mow_gui_make_default_state()
                project.procpramarray[project.ci].mow_state = method_state
            endif
            
            mas_mow_init_obj, $
                deco_level=(*method_state).mow_deco_order, $
                radius=(*method_state).mow_radius, $
                comp_meth=(*method_state).mow_comp_method, $
                nnls_so_location=(*method_state).mow_use_nnls_so eq 1 ? (*method_state).mow_nnls_so_location : '', $
                mow_obj=method, /bulk
        end
        
    endcase
    
    return, method
    
end

pro mas_hardi_tractography_gui_event, event

    name = widget_info(event.id, /uname)
    print, name
    ;;help, event, /structure
    
end

pro mas_hardi_tractography_gui_seeddens_event, event
    
    if (event.type eq 0 or event.type eq 2) then begin
    
        widget_control, event.id, get_value=val
        val = long(val)
        
        seeding_algo = mas_hardi_tractography_gui_getSeedAlgorithm(event)
        
        if (seeding_algo eq 0) then begin
            seeds_per_vxl = val^3
        endif else begin
            seeds_per_vxl = val
        endelse
        
        txt_seeds_vxl = widget_info(event.top, find_by_uname='txt_seeds_vxl')
        widget_control, txt_seeds_vxl, set_value=string(seeds_per_vxl, format='(I0)')
    endif

end

pro mas_hardi_tractography_gui_btn_event, event

    widget_control, event.top, /destroy

end

pro mas_hardi_tractography_gui_btn, event

    catch, error_state
    if (error_state ne 0) then begin
        ;; If one of the reconstruction objects fail to materialize,
        ;; we need to know why.
        catch, /cancel
        void = dialog_message(['An error occurred when trying to create the reconstruction object.', $
                              'Please make sure that a scan is loaded, and that it is a diffusion-weighed scan.', $
                              'If there is still problem, please let a developer know.'], $
                              /center, /error)
        print, !ERROR_STATE.MSG
        return
    endif
    
    base = widget_base(title="Configure Probability", group_leader=event.top, /modal, /column)
    
    dl_odf_select = widget_info(event.top, find_by_uname='dl_prob_method')
    odf_type = widget_info(dl_odf_select, /droplist_select)
    
    case odf_type of
        0: begin
            mas_odf_gui_make, base, 0
            void = widget_button(base, value="Apply Changes", event_pro='mas_hardi_tractography_gui_btn_event')
            widget_control, base, /realize
            xmanager, 'mas_odf_gui', base, cleanup='mas_odf_gui_cleanup'
        end
        1: begin
            mas_dot_gui_make, base, 0
            void = widget_button(base, value="Apply Changes", event_pro='mas_hardi_tractography_gui_btn_event')
            widget_control, base, /realize
            xmanager, 'mas_dot_gui', base, cleanup='mas_dot_gui_cleanup'
        end
        2: begin
            mas_mow_gui_make, base, 0
            void = widget_button(base, value="Apply Changes", event_pro='mas_hardi_tractography_gui_btn_event')
            widget_control, base, /realize
            xmanager, 'mas_mow_gui', base, cleanup='mas_mow_gui_cleanup'
        end
        else: print, "Reconstruction method not available."
        
    endcase
    catch, /cancel
         
end

pro mas_hardi_tractography_gui_anroi_manage, event

    common scan_data
    
    name = widget_info(event.id, /uname)
    
    an_roi_list = widget_info(event.top, find_by_uname='list_anal_rois')
    
    widget_control, an_roi_list, get_uvalue=plist_values

    work_dir = mas_hardi_tractography_gui_getWorkDir(event)
    if (work_dir eq '') then begin
        work_dir = project.current_path
    endif
    
    case name of
    
        'btn_add_anal_roi': begin
        
            file = dialog_pickfile(path=work_dir, /read, /multiple, filter='*.nii;*.nii.gz')
            if (n_elements(file) eq 0) then return
            if (file[0] eq '') then return
            
            for f = 0, n_elements(file)-1 do begin
                if (file_test(file[f], /directory)) then begin
                    void = dialog_message([file[f], 'is a directory.'], /error, /center)
                    endif else if (not file_test(file[f], /read, /regular)) then begin
                    void = dialog_message(['Cannot read:', file[f], 'Please check permissions'], /error, /center)
                    endif
            
                if (ptr_valid(plist_values) && (*plist_values)[0] ne '[None]') then begin
                    (*plist_values) = [ *plist_values, file[f] ]
                endif else begin
                    plist_values = ptr_new([ file[f] ])
                    endelse
            
                widget_control, an_roi_list, set_value=*plist_values, set_uvalue=plist_values
            endfor
        end
        
        'btn_rem_anal_roi': begin
            selected = widget_info(an_roi_list, /list_select)
            for n = 0, n_elements(*plist_values)-1 do begin
                if (selected eq n) then continue
                new_values = (n_elements(new_values) eq 0) ? [ (*plist_values)[n] ] : [ new_values, (*plist_values)[n] ]
            endfor
            if (n_elements(new_values) eq 0) then new_values = ['[None]']
            ptr_free, plist_values
            widget_control, an_roi_list, set_value=new_values, set_uvalue=ptr_new(new_values)
        end
        
    endcase
    
end

pro mas_hardi_tractography_gui_setfibfile, event

    common scan_data
    
    name = widget_info(event.id, /uname)
    
    txt_fibfile_wid = widget_info(event.top, find_by_uname='txt_fibfile')
    
    if (not widget_info(txt_fibfile_wid, /valid_id)) then return

    work_dir = mas_hardi_tractography_gui_getWorkDir(event)
    if (work_dir eq '') then begin
        work_dir = project.current_path
    endif
    
    case name of 
    
        'btn_fibfile_choose' : begin
            file = dialog_pickfile(path=work_dir, /read, filter='*.sav;*.trk')
            if (file eq '') then return
            if (file_test(file, /directory)) then begin
                junk = dialog_message([file,'is a directoy. Please choose a file.'], /error, /center)
            endif else if (not file_test(file, /read, /regular)) then begin
                junk = dialog_message(['Unable to read:', file], /error, /center)
                return
            endif
            widget_control, txt_fibfile_wid, set_value=file
        end
        
    endcase

end

pro mas_hardi_tractography_gui_setmatfile, event

    name = widget_info(event.id, /uname)
    
    txt_matfile_wid = widget_info(event.top, find_by_uname='txt_matfile')
    
    if (not widget_info(txt_matfile_wid, /valid_id)) then return
    
    case name of 
    
        'btn_matfile_choose' : begin
            file = dialog_pickfile(/read)
            if (file eq '') then return
            if (not file_test(file, /read)) then begin
                junk = dialog_message(['Unable to read:', file], /error, /center)
                return
            endif
            widget_control, txt_matfile_wid, set_value=file
        end
        
        'btn_matfile_clear' : begin
            widget_control, txt_matfile_wid, set_value='<Identity>'
        end
        
    endcase

end

pro mas_hardi_tractography_gui_setwd, event

    txt = widget_info(event.top, find_by_uname='txt_working_dir')
    if (not widget_info(txt, /valid_id)) then return
    
    oldwd = mas_hardi_tractography_gui_getWorkDir(event)
    dir = dialog_pickfile(path=oldwd, /directory, /read)
    if (dir eq '') then return
    
    if (not file_test(dir, /write, /directory)) then begin
        junk = dialog_message(['Unable to write to:', dir], /error, /center)
        return
    endif 
    
    widget_control, txt, set_value=dir

end

pro mas_hardi_tractography_gui_prepd, event

    txt = widget_info(event.top, find_by_uname='txt_prep_dir')
    if (not widget_info(txt, /valid_id)) then return
    
    wd = mas_hardi_tractography_gui_getWorkDir(event)
    dir = dialog_pickfile(/directory, /read, title='Choose the directory containing the preprocessed directions.', $
                            path=wd)
    if (dir eq '') then return
    
    files = file_search(dir+'V_[0-9][0-9]*.{nii,nii.gz}', /test_read, /test_regular, count=nfiles)
    if (nfiles eq 0) then begin
        void = dialog_message(['Count not find any direction files in',$
                              dir, $
                              'Appropriate files should be named: "V_XX.nii",', $
                              'where XX is a number (WITH leading zero.)'], /error, /center)
        widget_control, txt, set_value='', set_uvalue=['']
    endif else begin
        base_files = replicate('    ', n_elements(files)) + file_basename(files)
        void = dialog_message(['Found the following direction files:', $
                              base_files], /information, /center)
        widget_control, txt, set_value=dir, set_uvalue=files
    endelse
    
    

end

pro mas_hardi_tractography_gui_roi_choose, event

    common scan_data
    
    name = widget_info(event.id, /uname)
    identifiers = stregex(name, '^btn_roi_choose_(.+$)', /extract, /subexpr)
    txt = 'txt_roi_path_'+identifiers[1]
    wid_txt = widget_info(event.top, find_by_uname=txt)
    if (not widget_info(wid_txt, /valid_id)) then return
    work_dir = mas_hardi_tractography_gui_getWorkDir(event)
    if (work_dir eq '') then begin
        work_dir = project.current_path
    endif
    txt_val = dialog_pickfile(path=work_dir, /read, filter='*.nii;*.nii.gz')   
    if (txt_val eq '') then begin
        return
    endif
    widget_control, wid_txt, set_value=txt_val, sensitive=1

end


pro mas_hardi_tractography_gui_roi_enable, event

    common scan_data
    name = widget_info(event.id, /uname)
    identifiers = stregex(name, '^btn_roi_enable_(.+$)', /extract, /subexpr)
    txt = 'txt_roi_path_'+identifiers[1]
    wid_txt = widget_info(event.top, find_by_uname=txt)
    if (not widget_info(wid_txt, /valid_id)) then return
    
    btn = 'btn_roi_choose_'+identifiers[1]
    wid_btn = widget_info(event.top, find_by_uname=btn)
    if (not widget_info(wid_btn, /valid_id)) then return

    txt = 'dl_roi_type_'+identifiers[1]
    wid_dl = widget_info(event.top, find_by_uname=txt)
    if (not widget_info(wid_dl, /valid_id)) then wid_dl = wid_btn
    
    work_dir = mas_hardi_tractography_gui_getWorkDir(event)
    if (work_dir eq '') then begin
        work_dir = project.current_path
    endif
    
    widget_control, wid_txt, get_value=txt_val
    if (event.select) then begin
        if (txt_val eq '') then begin
            txt_val = dialog_pickfile(path=work_dir, /read, filter='*.nii;*.nii.gz')
            if (txt_val eq '') then begin
                widget_control, event.id, set_button=0
                return
            endif
        endif
        widget_control, wid_txt, set_value=txt_val, sensitive=1
        widget_control, wid_btn, sensitive=1
        widget_control, wid_dl,  sensitive=1
    endif else begin
        widget_control, wid_txt, sensitive=0
        widget_control, wid_btn, sensitive=0
        widget_control, wid_dl,  sensitive=0
    endelse
    
end

;pro mas_hardi_tractography_gui_bootstrap, event
;
;    void = dialog_message('Are you sure?', /question, /center)
;    if (void eq 'No') then return
;    
;    mas_odf_init_obj, lmax=4, /bulk, odf_obj=odf
;    
;    mas_probtrack_bootstrap_data, odf, nbootstrap=100
;
;    obj_destroy, odf
;    
;end

pro mas_hardi_tractography_gui_estimate, event

    common scan_data
    name = widget_info(event.id, /uname)
    
    dl_seed_algo  = widget_info(event.top, find_by_uname='dl_seed_algorithm')
    txt_seed_dens = widget_info(event.top, find_by_uname='txt_seed_dens')
    
    algo_type = widget_info(dl_seed_algo, /DROPLIST_SELECT)
    widget_control, txt_seed_dens, get_value=dens_str
    seed_density = fix(dens_str[0])
    
    ;; 0 = "uniform n^3 method", 1 = "Random method"
    multiplier = algo_type eq 0 ? seed_density^3 : seed_density
    
    case name of
    
        'btn_sd_roi_est': begin
            junk = dialog_message(['Not implemented yet...'], $
                /center, title='Seed point estimate', /information)
                
        end
        
        'btn_sd_fa_est': begin
            if (not ptr_valid(project.dataarray[project.ci].frac_ani)) then begin
                void = dialog_message(['This dataset does not appear to be a diffusion', $
                                        'weighed dataset.'], /center, /error)
                return
            endif
            txt_fa = widget_info(event.top, find_by_uname='txt_sd_fa')
            if (not widget_info(txt_fa, /valid_id)) then return
            widget_control, txt_fa, get_value=fa_thr
            fa_thr = float(fa_thr[0])
            ind = where(*project.dataarray[project.ci].frac_ani ge fa_thr, count)
            
            nseeds = count * multiplier
            
            junk = dialog_message(['Number of voxels: '+string(count, format='(I0)'), $
                'Number of seed points / voxel: '+string(multiplier, format='(I0)'), $
                'Total number of seed points: '+string(nseeds, format='(I0)')], $
                /center, title='Seed point estimate', /information)
        end
        
    endcase
    
end

pro mas_hardi_tractography_gui_showvoxel, event

    common scan_data, project
    ci = project.ci
    
    txt_sd_vxl_x = widget_info(event.top, find_by_uname='txt_sd_vxl_x')
    txt_sd_vxl_y = widget_info(event.top, find_by_uname='txt_sd_vxl_y')
    txt_sd_vxl_z = widget_info(event.top, find_by_uname='txt_sd_vxl_z')

    widget_control, txt_sd_vxl_x, get_value=seed_x
    widget_control, txt_sd_vxl_y, get_value=seed_y
    widget_control, txt_sd_vxl_z, get_value=seed_z
    
    if (seed_x eq '' and seed_y eq '' and seed_z eq '') then begin
        seed_x = project.procpramarray[ci].fdim_start
        seed_y = project.procpramarray[ci].pdim_start
        seed_z = project.procpramarray[ci].sdim_start
    endif
    
    vox = [ fix(seed_x[0]), fix(seed_y[0]), fix(seed_z[0]) ]
    proj_dims = (size(*project.dataarray[project.ci].state1, /dimensions))[0:2]
    
    bad_voxel = 0
    if (vox[0] lt 0 or vox[0] ge proj_dims[0]) then bad_voxel = 1
    if (vox[1] lt 0 or vox[1] ge proj_dims[1]) then bad_voxel = 1
    if (vox[2] lt 0 or vox[2] ge proj_dims[2]) then bad_voxel = 1
    
    if (bad_voxel) then begin
        junk = dialog_message('Voxel out of range: ['+strjoin(string(proj_dims, format='(I0)'), ', ')+']', /error, /center)
        return
    endif

    
    prob_thr = mas_hardi_tractography_gui_getProbThr(event)
    tess_level = mas_hardi_tractography_gui_getTessLevel(event)
    
    btn_compare_adt = widget_info(event.top, find_by_uname='btn_compare_adt')
    compare_adt = widget_info(btn_compare_adt, /button_set)
    
    ;if (mas_hardi_tractography_gui_isVectorFieldEnabled(event)) then begin


    ;endif else begin
        odf_obj   = mas_hardi_tractography_gui_getODFObj(event)    
        mas_hardi_tractography_show_voxel, odf_obj, vox, $
                                           prob_thr=prob_thr, $
                                           recon_tess_level=tess_level, $
                                           restart_tess_level=0, $
                                           nbootstrap=0, compare_adt=compare_adt
    ;endelse
    
    obj_destroy, odf_obj
    
end


pro mas_hardi_tractography_gui_trackvoxel, event
    
    common scan_data

    if (not ptr_valid(project.imndarray[project.ci].bval_array)) then begin
        if (not mas_hardi_tractography_gui_isVectorFieldEnabled(event)) then begin
            void = dialog_message(['Data set does not appear to be a diffusion', $
                                'weighted data set.'], /error, /center)
            return
        endif
    endif
    
    txt_sd_vxl_x = widget_info(event.top, find_by_uname='txt_sd_vxl_x')
    txt_sd_vxl_y = widget_info(event.top, find_by_uname='txt_sd_vxl_y')
    txt_sd_vxl_z = widget_info(event.top, find_by_uname='txt_sd_vxl_z')

    widget_control, txt_sd_vxl_x, get_value=seed_x
    widget_control, txt_sd_vxl_y, get_value=seed_y
    widget_control, txt_sd_vxl_z, get_value=seed_z
    
    vox = [ fix(seed_x[0]), fix(seed_y[0]), fix(seed_z[0]) ]
    proj_dims = (size(*project.dataarray[project.ci].state1, /dimensions))[0:2]
    
    bad_voxel = 0
    if (vox[0] lt 0 or vox[0] ge proj_dims[0]) then bad_voxel = 1
    if (vox[1] lt 0 or vox[1] ge proj_dims[1]) then bad_voxel = 1
    if (vox[2] lt 0 or vox[2] ge proj_dims[2]) then bad_voxel = 1
    
    if (bad_voxel) then begin
        junk = dialog_message('Voxel out of range: ['+strjoin(string(proj_dims, format='(I0)'), ', ')+']', /error, /center)
        return
    endif

    mask = ptr_new(bytarr(proj_dims))
    (*mask)[vox[0], vox[1], vox[2]] = 1
    
    interp   = mas_hardi_tractography_gui_getInterp(event)
    density  = mas_hardi_tractography_gui_getSeedDensity(event)
    seed_alg = mas_hardi_tractography_gui_getSeedAlgorithm(event)
    step_size= mas_hardi_tractography_gui_getStepSize(event)
    anis_thr = mas_hardi_tractography_gui_getAnisThr(event)
    prob_thr = mas_hardi_tractography_gui_getProbThr(event)
    angle_thr = mas_hardi_tractography_gui_getAngleThr(event)
    tess_level = mas_hardi_tractography_gui_getTessLevel(event)
    work_dir = mas_hardi_tractography_gui_getWorkDir(event)
    if (mas_hardi_tractography_gui_isSplittingEnabled(event)) then begin
        max_split = mas_hardi_tractography_gui_getMaximumSplit(event)
    endif else begin
        max_split = 0
    endelse
    
    if (mas_hardi_tractography_gui_isVectorFieldEnabled(event)) then begin
        prepdirs = mas_hardi_tractography_loadPrepDirs(event)
        if (total(ptr_valid(prepdirs)) ne n_elements(prepdirs)) then begin
            void = dialog_message(['Error reading preprocessed directions.', $
                                   'Please check to make sure that the directory', $
                                   'is set and contains preprocessed direction files.'], /error, /center)
            return
        endif
        prepdirs = ptr_new(prepdirs, /no_copy)
        mas_hardi_tractography_track, mask, handle=state, density=density, odf_obj=odf_obj, working_dir=work_dir, $
                                    randomize=(seed_alg eq 0) ? 0 : 1, interp_method=interp, prob_thr=prob_thr, $
                                    step_size=step_size, anis_thr=anis_thr, angle_thr=angle_thr, $
                                    maximum_split=max_split, preprocessed_directions=prepdirs, /track_vector_field
    endif else begin
        odf_obj   = mas_hardi_tractography_gui_getODFObj(event)    
        mas_hardi_tractography_track, mask, handle=state, density=density, odf_obj=odf_obj, working_dir=work_dir, tess_level=tess_level, $
                                    randomize=(seed_alg eq 0) ? 0 : 1, interp_method=interp, prob_thr=prob_thr, $
                                    step_size=step_size, anis_thr=anis_thr, angle_thr=angle_thr, $
                                    maximum_split=max_split
        mas_hardi_tractography_addROI, state, mask
        obj_destroy, odf_obj
    endelse
    
    ptr_free, mask
    obj_destroy, odf_obj
    
end

pro mas_hardi_tractography_gui_trackfa, event

    common scan_data

    if (not ptr_valid(project.imndarray[project.ci].bval_array)) then begin
        if (not mas_hardi_tractography_gui_isVectorFieldEnabled(event)) then begin
            void = dialog_message(['Data set does not appear to be a diffusion', $
                                'weighted data set.'], /error, /center)
            return
        endif
    endif
    
    txt_fa = widget_info(event.top, find_by_uname='txt_sd_fa')
    if (not widget_info(txt_fa, /valid_id)) then return
    widget_control, txt_fa, get_value=fa_thr
    fa_thr = float(fa_thr[0])
    ind = where(*project.dataarray[project.ci].frac_ani ge fa_thr, count)
    
    if (count gt 100) then begin
        void = dialog_message(['FA threshold tractography generally results in', $
                               'a large number of seed points and may take a while.', $
                                 'Do you want to continue?'], /center, /question)
        if (void eq 'No') then return
    endif
    
    proj_dims = (size(*project.dataarray[project.ci].state1, /dimensions))[0:2]
    mask = ptr_new(bytarr(proj_dims))
    (*mask)[ind] = 1
    
    interp   = mas_hardi_tractography_gui_getInterp(event)
    density  = mas_hardi_tractography_gui_getSeedDensity(event)
    seed_alg = mas_hardi_tractography_gui_getSeedAlgorithm(event)
    step_size= mas_hardi_tractography_gui_getStepSize(event)
    anis_thr = mas_hardi_tractography_gui_getAnisThr(event)
    prob_thr = mas_hardi_tractography_gui_getProbThr(event)
    angle_thr = mas_hardi_tractography_gui_getAngleThr(event)
    tess_level = mas_hardi_tractography_gui_getTessLevel(event)
    work_dir = mas_hardi_tractography_gui_getWorkDir(event)    
    if (mas_hardi_tractography_gui_isSplittingEnabled(event)) then begin
        max_split = mas_hardi_tractography_gui_getMaximumSplit(event)
    endif else begin
        max_split = 0
    endelse

    if (mas_hardi_tractography_gui_isVectorFieldEnabled(event)) then begin
        prepdirs = mas_hardi_tractography_loadPrepDirs(event)
        if (total(ptr_valid(prepdirs)) ne n_elements(prepdirs)) then begin
            void = dialog_message(['Error reading preprocessed directions.', $
                                   'Please check to make sure that the directory', $
                                   'is set and contains preprocessed direction files.'], /error, /center)
            return
        endif
        prepdirs = ptr_new(prepdirs, /no_copy)    
        mas_hardi_tractography_track, mask, handle=state, density=density, odf_obj=odf_obj, working_dir=work_dir, $
                                    randomize=(seed_alg eq 0) ? 0 : 1, interp_method=interp, prob_thr=prob_thr, $
                                    step_size=step_size, anis_thr=anis_thr, angle_thr=angle_thr, $
                                    maximum_split=max_split, preprocessed_directions=prepdirs, /track_vector_field
    endif else begin
        odf_obj   = mas_hardi_tractography_gui_getODFObj(event)    
        mas_hardi_tractography_track, mask, handle=state, density=density, odf_obj=odf_obj, working_dir=work_dir, tess_level=tess_level,$
                                    randomize=(seed_alg eq 0) ? 0 : 1, interp_method=interp, prob_thr=prob_thr, $
                                    step_size=step_size, anis_thr=anis_thr, angle_thr=angle_thr, $
                                    maximum_split=max_split
        obj_destroy, odf_obj
    endelse
            
    ptr_free, mask

end

pro mas_hardi_tractography_gui_trackroi, event

    common scan_data
    
    if (not ptr_valid(project.imndarray[project.ci].bval_array)) then begin
        if (not mas_hardi_tractography_gui_isVectorFieldEnabled(event)) then begin
            void = dialog_message(['Data set does not appear to be a diffusion', $
                                'weighted data set.'], /error, /center)
            return
        endif
    endif

    bound_roi = ptr_new()
    exclusion_roi = ptr_new()

    progressBar = obj_new('progressbar', title="Please Wait...", text="Loading ROIs")
    progressBar->Start

    for i = 1, 7 do begin
    
        if (i eq 6) then begin
            ident = 'bound'
        endif else if (i eq 7) then begin
            ident = 'exclusion'
        endif else begin
            ident = string(i, format='(I0)')
        endelse
        
        btn_roi_enable = widget_info(event.top, find_by_uname='btn_roi_enable_'+ident)
        if (widget_info(btn_roi_enable, /button_set) eq 0) then continue
        
        txt_file_name = widget_info(event.top, find_by_uname='txt_roi_path_'+ident)
        widget_control, txt_file_name, get_value=file_name
        
        if (not file_test(file_name, /read)) then begin
            void = dialog_message(["Skipping: "+file_name, 'Cannot read file.'], /error, /center)
            continue
        endif
        
        nif = mas_read_nifti(nifti_filename=file_name, read_status=status)
        
        if (status eq 0) then begin
            void = dialog_message(["Skipping: "+file_name, 'Cannot read file.'], /error, /center)
        endif else begin
            if (i eq 6) then begin ;; bound ROI
                bound_roi = nif.voxel_data
            endif else if (i eq 7) then begin ;; exclusion ROI
                exclusion_roi = nif.voxel_data
            endif else begin
                dl_roi_type = widget_info(event.top, find_by_uname='dl_roi_type_'+ident)
                roi_type = widget_info(dl_roi_type, /droplist_select)
                if (not widget_info(dl_roi_type, /valid_id)) then continue          
                case roi_type of 
                    0: begin
                        seed_rois = (n_elements(seed_rois) eq 0) ? [ nif.voxel_data ] : [ seed_rois, nif.voxel_data ]
                        anal_rois = (n_elements(anal_rois) eq 0) ? [ nif.voxel_data ] : [ anal_rois, nif.voxel_data ]
                    end 
;; These are no longer relavent since we visualize in TrackVis, or use the "Visualization" tab.
;;                    1: seed_rois = (n_elements(seed_rois) eq 0) ? [ nif.voxel_data ] : [ seed_rois, nif.voxel_data ]
;;                    2: anal_rois = (n_elements(anal_rois) eq 0) ? [ nif.voxel_data ] : [ anal_rois, nif.voxel_data ]
                    else: print, "Unknown ROI type: "+string(roi_type, format='(I0)')
                 endcase
            endelse
        endelse
    
        progressBar->Update, float(i)/7.0 * 100.0
        
    endfor

    progressBar->Destroy
    
    if (n_elements(seed_rois) eq 0) then begin
        void = dialog_message("No seeding rois selected.", /error, /center)
        return
    endif
    
    interp    = mas_hardi_tractography_gui_getInterp(event)
    density   = mas_hardi_tractography_gui_getSeedDensity(event)
    seed_alg  = mas_hardi_tractography_gui_getSeedAlgorithm(event)
    step_size = mas_hardi_tractography_gui_getStepSize(event)
    anis_thr  = mas_hardi_tractography_gui_getAnisThr(event)
    prob_thr  = mas_hardi_tractography_gui_getProbThr(event)
    angle_thr = mas_hardi_tractography_gui_getAngleThr(event)
    tess_level = mas_hardi_tractography_gui_getTessLevel(event)
    work_dir = mas_hardi_tractography_gui_getWorkDir(event)    
    if (mas_hardi_tractography_gui_isSplittingEnabled(event)) then begin
        max_split = mas_hardi_tractography_gui_getMaximumSplit(event)
    endif else begin
        max_split = 0
    endelse

    if (mas_hardi_tractography_gui_isVectorFieldEnabled(event)) then begin
        prepdirs = mas_hardi_tractography_loadPrepDirs(event)
        if (total(ptr_valid(prepdirs)) ne n_elements(prepdirs)) then begin
            void = dialog_message(['Error reading preprocessed directions.', $
                                   'Please check to make sure that the directory', $
                                   'is set and contains preprocessed direction files.'], /error, /center)
            return
        endif
        prepdirs = ptr_new(prepdirs, /no_copy)    
        mas_hardi_tractography_track, seed_rois, handle=state, density=density, odf_obj=odf_obj, working_dir=work_dir, $
                                    randomize=(seed_alg eq 0) ? 0 : 1, interp_method=interp, prob_thr=prob_thr, $
                                    step_size=step_size, anis_thr=anis_thr, bound_roi=bound_roi, exclusion_roi=exclusion_roi, angle_thr=angle_thr, $
                                    maximum_split=max_split, preprocessed_directions=prepdirs, /track_vector_field
    endif else begin
        odf_obj   = mas_hardi_tractography_gui_getODFObj(event)    
        mas_hardi_tractography_track, seed_rois, handle=state, density=density, odf_obj=odf_obj, working_dir=work_dir, tess_level=tess_level,$
                                    randomize=(seed_alg eq 0) ? 0 : 1, interp_method=interp, prob_thr=prob_thr, $
                                    step_size=step_size, anis_thr=anis_thr, bound_roi=bound_roi, exclusion_roi=exclusion_roi, angle_thr=angle_thr, $
                                    maximum_split=max_split
        obj_destroy, odf_obj
        mas_hardi_tractography_addROI, state, anal_rois    
    endelse
    
    if (ptr_valid(bound_roi)) then ptr_free, bound_roi
    if (ptr_valid(exclusion_roi)) then ptr_free, exclusion_roi
    if (n_elements(anal_rois) gt 0) then ptr_free, anal_rois
    ptr_free, seed_rois

end

pro mas_hardi_tractography_gui_launchvis, event

    common scan_data
    
    txt_fib_file = widget_info(event.top, find_by_uname='txt_fibfile')
    widget_control, txt_fib_file, get_value=fib_filename
    
    if (fib_filename eq '') then begin
        void = dialog_message('Please select a fiber data file.', /error, /center)
        return
    endif else if (not file_test(fib_filename, /read, /regular)) then begin
        void = dialog_message(['Cannot read:', fib_filename], /error, /center)
        return
    endif
    
    list_anal_roi = widget_info(event.top, find_by_uname='list_anal_rois')
    widget_control, list_anal_roi, get_uvalue=plist_values
    
    if (not ptr_valid(plist_values)) then return
    
    nrois = 0
    
    for r = 0, n_elements(*plist_values)-1 do begin
    
        file = (*plist_values)[r]
        if (not file_test(file, /read, /regular)) then begin
            print, "Skipping: "+file
            continue
        endif
        
        nif = mas_read_nifti(nifti_filename=file, read_status=rs)
        if (rs eq 0) then begin
            print, "mas_read_nifti: returned read status=0"
            continue
        endif
        
        oroi = obj_new('mas_hardi_ROI', nif.voxel_data, name="ROI_N")
    
        subs =  oroi->getSubRegions()
    
        for i = 0, n_elements(subs)-1 do begin
            tmp_seeds = oroi->makeSeeds(subs[i], density=1)
            if (tmp_seeds[0] eq -1) then continue
            roi_arr = (nrois eq 0) ? ptr_new(tmp_seeds) : [ roi_arr, ptr_new(tmp_seeds) ]
            roi_names = (nrois eq 0) ? file : [ roi_names, file ]
            nrois++
            endfor
    
        obj_destroy, oroi
        ptr_free, nif.voxel_data
    
    endfor

    working_dir = mas_hardi_tractography_gui_getWorkDir(event)
    if (working_dir eq '') then begin
        working_dir = project.current_path
    endif
    
    fiber_flip =  mas_hardi_tractography_gui_getFiberFlip(event)
    
    if (strmid(fib_filename, 2, /reverse_offset) eq 'sav') then begin
        restore, fib_filename
    endif else begin
        
        openr, lun, fib_filename, /get_lun
        tv_hdr = create_struct(name='MAS_HARDI_TRACTOGRAPHY_TRACKVIS_HEADER')
        readu, lun, tv_hdr

        if (string(tv_hdr.id_string[0:4]) ne 'TRACK') then begin
            void = dialog_message([fib_filename, 'does not seem to be a valid trackvis file.'], /center,/error)
            return
        endif
        
        tv2mas_mat = diag_matrix(replicate(1.0, 4)) ;; transform from trackvis to mas orientation
        print, "Track Voxel Order: "+string(tv_hdr.voxel_order)
        if (1 or project.procpramarray[project.ci].have_orientation eq 1) then begin
            orient_mat = float(mas_orient_get_matrix())
            orient_code = strjoin(mas_orient_from_matrix(orient_mat, /lettercode))
            print, "Data Voxel Order: "+orient_code
            if (string(tv_hdr.voxel_order) ne orient_code) then begin
                tv_orient_mat = mas_orient_get_matrix(lettercode=string(tv_hdr.voxel_order)) ;; get the matrix for the trackvis orientation
                our_orient_mat = mas_orient_get_matrix() ;; the matrix for our orientation
                tv2mas_mat[0:2,0:2] = invert(float(tv_orient_mat)) # float(our_orient_mat)
                tv_hdr.dim = tv_hdr.dim
                print, "Transform: " & print, tv2mas_mat
            endif
        endif
        progressbar = obj_new('progressbar', title='Reading Track File...', text='Reading Track File...', /fast_loop)
        progressbar->Start

        ct = tv_hdr.n_count & print, "Track Count: ", ct
        nf = 0L
        seeds = lonarr(ct)
        fibers = ptrarr(ct)
        ;; to perform the transform, we need to move the tracks to the
        ;; data set center, apply the transform, then shift back to the
        ;; new orientation's center
        to_origin   = diag_matrix(replicate(1.0,4)) & to_origin[3,0:2]   = (-tv_hdr.dim/2.0) 
        from_origin = diag_matrix(replicate(1.0,4)) & from_origin[3,0:2] = abs((tv_hdr.dim/2.0) # tv2mas_mat[0:2,0:2])
        final_mat = diag_matrix([[1.0/tv_hdr.voxel_size], 1.0]) # to_origin # diag_matrix([fiber_flip]) # tv2mas_mat # from_origin 
        for f = 0L, ct-1 do begin

            if (eof(lun)) then break
       
            readu, lun, nf
            tmp_fib = fltarr(3,nf)
            readu, lun, tmp_fib
            
            if (nf lt 3) then continue
            
            tmp_fib = vert_t3d(tmp_fib, matrix=final_mat, /no_divide, /no_copy)

            seeds[f]  = n_elements(tmp_fib)/3/2.
            fibers[f] = ptr_new(tmp_fib, /no_copy)
            
            progressbar->Update, float(f)/float(ct) * 100.0
            
            if (f mod 20 eq 0) then begin
                if (progressbar->CheckCancel()) then begin
                    ptr_free, fibers
                    close, lun
                    free_lun, lun
                    progressbar->Destroy
                    return
                endif
            endif
        endfor
        
        progressbar->Destroy
        close, lun
        free_lun, lun
        
        fiber_data = ptr_new({ fibers: fibers, $
            seed_points: seeds, $
            transform: diag_matrix(replicate(1.0,4)), $
            voxel_dims: tv_hdr.voxel_size, $
            branch_levels: [0], $
            dimensions: tv_hdr.dim, $
            nfibers: tv_hdr.n_count, $
            filedir: '', $
            tractfile: '', $
            roifile: '' })
    
    endelse
    
    if (n_elements(roi_arr) ne 0) then begin
        roi_data = ptr_new({ vertices: roi_arr, names: roi_names })
        view_fibers_hardi, fiber_data, roi_data=roi_data, $
                           /show_axis, /in_mas, handle=state, $
                           working_dir=working_dir;, &
                           ;fstep=10;;, project_dataset=project_dataset

    endif else begin
        view_fibers_hardi, fiber_data, /show_axis, /in_mas, handle=state, $
                           working_dir=working_dir;, &
                           ;fstep=10;;, project_dataset=project_dataset
    endelse
end

pro mas_hardi_tractography_gui

    common common_widgets
    common scan_data
    
    base = widget_base(title="HARDI Tractography", column=2, $
                       group_leader=WID_BASE_MAIN)
    
    base_left = widget_base(base, row=1)
    base_right = widget_base(base, /column)
    
    tabbase = widget_tab(base)
    
    ;;tab_prep = widget_base(tabbase, title="Preprocessing", /column)
    
    tab_param = widget_base(tabbase, title="Tracking Parameters", column=2)
    
    tab_seeds = widget_base(tabbase, title="Regions of Interest")
    
    tab_visual = widget_base(tabbase, title="Visualization", row=2)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Preprocessing Tab
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
;    prep_param = widget_base(tab_prep, /column)
;    b = widget_base(prep_param, /column, /frame)
;    base_odf_type = widget_base(b, /row)
;    
;    lbl_prob_method = Widget_Label(base_odf_type,  $
;        UNAME='lbl_prob_method', VALUE='ODF Computation:')
;        
;    dl_prob_method = Widget_Droplist(base_odf_type,  $
;        UNAME='dl_prob_method',VALUE=[ 'QBI', 'DOT', 'MOW' ])
;        
;    btn_config_prob = widget_button(base_odf_type, $
;        UNAME='btn_config_prob', VALUE='Configure', $
;        event_pro='mas_hardi_tractography_gui_btn')
;        
;    base_odf_val = widget_base(b, /row)
;    lbl_odf_val = widget_label(base_odf_val, value='ODF Value Thr:', uname='lbl_odf_val')
;    txt_odf_val = widget_text(base_odf_val, value='0.5', xsize=8, /editable, uname='txt_odf_val')
;    
;    b = widget_base(tab_prep, /row, /frame)
;    lbl_num_directions = widget_label(b, value='Number of directions to resolve per voxel:', /align_left)
;    txt_num_directions = widget_text(b, xsize=4, value='4', uname='txt_num_directions', /editable)
;   
;    b = widget_base(tab_prep, /column, /align_left, /frame)
;    lbl_prep_output_dir = widget_label(b, value='Save direction output files to:', /align_left, uname='lbl_prep_output_dir')
;    b1 = widget_base(b, /row)
;    txt_prep_output_dir = widget_Text(b1, value=project.current_path, xsize=50, uname='txt_prep_output_dir', /editable)
;    btn_prep_output_dir = widget_button(b1, value="Choose", uname='btn_prep_output_dir', $
;                                    event_pro='mas_hardi_tractography_gui_setwd')

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Tracking Parameters Tab
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    param_left = widget_base(tab_param, /column);, /frame)
    
    b = widget_base(param_left, /frame, /column)
    base_fa_halt = widget_base(b, /row)
    lbl_fa_halt = widget_label(base_fa_halt, value="FA Halting Thr:", uname='lbl_fa_halt')
    txt_fa_halt = widget_text(base_fa_halt, value='0.05', xsize=8, /editable, uname='txt_fa_halt')
     
    base_turn_angle = widget_base(b, /row)
    lbl_turn_angle = widget_label(base_turn_angle, value='Turn Angle Thr:', uname='lbl_turn_angle')
    txt_turn_angle = widget_text(base_turn_angle, value='50.', xsize=8, /editable, uname='txt_turn_angle')
    
    b = widget_base(param_left, /frame, /column)
    base_seed_dens = widget_base(b, /row)
    lbl_seed_type = widget_label(base_seed_dens, value='Seed Algorithm:', uname='lbl_seed_type')
    dl_seed_algorithm = widget_droplist(base_seed_dens, value=['Uniform (density^3)', 'Random (Uniform Dist.)'], uname='dl_seed_algorithm')
    base_seed_dens = widget_base(b, /row)
    lbl_seed_dens = widget_label(base_seed_dens, value='Seed Density:', uname='lbl_seed_dens')
    txt_seed_dens = widget_text(base_seed_dens, value='1', xsize=4, /editable, uname='txt_seed_dens', /all_events, event_pro='mas_hardi_tractography_gui_seeddens_event')
    lbl_nseeds = widget_label(base_seed_dens, value="Seeds/Vxl:")
    txt_seeds_vxl = widget_text(base_seed_dens, xsize=10, uname='txt_seeds_vxl', value='1')
    
    b = widget_base(param_left, /column, /frame)
    base_odf_type = widget_base(b, /row)
    
    lbl_prob_method = Widget_Label(base_odf_type,  $
        uname='lbl_prob_method', VALUE='ODF Computation:')
        
    dl_prob_method = Widget_Droplist(base_odf_type,  $
        uname='dl_prob_method',VALUE=[ 'QBI', 'DOT', 'MOW' ])
        
    btn_config_prob = widget_button(base_odf_type, $
        uname='btn_config_prob', VALUE='Configure', $
        event_pro='mas_hardi_tractography_gui_btn')
        
    base_odf_val = widget_base(b, /row)
    lbl_odf_val = widget_label(base_odf_val, value='ODF Value Thr:', uname='lbl_odf_val')
    txt_odf_val = widget_text(base_odf_val, value='0.5', xsize=8, /editable, uname='txt_odf_val')
    dl_tess_level = widget_droplist(base_odf_val, title="Tessellation Level:", value=['3','4','5'], uname='dl_tess_level')
    widget_control, dl_tess_level, set_droplist_select=0
    btn_prep_base = widget_base(b, /base_align_center)
    btn_odf_prep = widget_button(btn_prep_base, value='Preprocess Dataset', uname='btn_odf_prep', $
                                 event_pro='mas_hardi_tractography_gui_preprocess_event')
    
    b = widget_base(param_left, /column, /frame)
    base_interpolation = widget_base(b, /row)
    
    lbl_interp = Widget_Label(base_interpolation,  $
        uname='lbl_interp', VALUE='Interpolation:')
        
    dl_interp_method = Widget_Droplist(base_interpolation,  $
        uname='dl_interp_method', VALUE=[ 'Nearest Neighbor', 'Tricubic', 'Trilinear' ])
        
    base_stepsize = widget_base(b, /row)
    lbl_stepsize = widget_label(base_stepsize, value="Step Size:", uname='lbl_stepsize')
    txt_stepsize = widget_text(base_stepsize, value='0.5', xsize=8, /editable, uname='txt_stepsize')
    
    param_right = widget_base(tab_param, /column);, /frame)
    
    b = widget_base(param_right, /column, /align_left, /frame)
    b1 = widget_base(b, /nonexclusive)
    btn_splitting_enabled = widget_button(b1, value="Copy and branch on multidirection", uname="btn_splitting_enabled")
    b2 = widget_base(b, /row)
    lbl = widget_label(b2, value="Maximum branches:")
    txt_max_split = widget_text(b2, value='5', xsize=8, /editable, uname="txt_max_split")

    b = widget_base(param_right, /column, /align_left, /frame)
    b1 = widget_base(b, /nonexclusive)
    btn_use_vector_field = widget_button(b1, value="Track using preprocessed directions", uname="btn_use_vector_field",$
                                         event_pro='mas_hardi_tractography_gui_togglevf')
    lbl_prep_dir = widget_label(b, value='Directory containing preprocessed directions:', /align_left, uname='lbl_prep_dir', sensitive=1)
    b2 = widget_base(b, /row)
    txt_prep_dir = widget_text(b2, xsize=50, uname='txt_prep_dir', sensitive=0)
    btn_prep_dir = widget_button(b2, value="Choose", uname='btn_prep_dir', sensitive=0, $
                                    event_pro='mas_hardi_tractography_gui_prepd')
    
    b = widget_base(param_right, /column, /align_left, /frame)
    lbl_working_dir = widget_label(b, value='Working Directory:', /align_left, uname='lbl_working_dir')
    b1 = widget_base(b, /row)
    txt_working_dir = widget_Text(b1, value=project.current_path, xsize=50, uname='txt_working_dir', /editable)
    btn_working_dir = widget_button(b1, value="Choose", uname='btn_working_dir', $
                                    event_pro='mas_hardi_tractography_gui_setwd')
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; ROI  Parameteres Tab
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    roi_base = widget_base(tab_seeds, /column)
    
    voxel_base = widget_base(roi_base, /row, /frame)
    lbl_sd_vxl_x = widget_label(voxel_base, value="Track Single Voxel X:", uname='lbl_sd_vxl_x')
    txt_sd_vxl_x = widget_text(voxel_base, value='', xsize=5, uname='txt_sd_vxl_x', /editable)
    lbl_sd_vxl_y = widget_label(voxel_base, value=" Y:", uname='lbl_sd_vxl_y')
    txt_sd_vxl_y = widget_text(voxel_base, value='', xsize=5, uname='txt_sd_vxl_y', /editable)
    lbl_sd_vxl_z = widget_label(voxel_base, value=" Z:", uname='lbl_sd_vxl_z')
    txt_sd_vxl_z = widget_text(voxel_base, value='', xsize=5, uname='txt_sd_vxl_z', /editable)
    btn_sd_track = widget_button(voxel_base, value="Track", uname='btn_sd_track', $
                                event_pro='mas_hardi_tractography_gui_trackvoxel')
    btn_sd_show = widget_button(voxel_base, value="Show Voxel", uname='btn_sd_show', $
                                event_pro='mas_hardi_tractography_gui_showvoxel')
    bb = widget_base(voxel_base, /nonexclusive)
    btn_comp_adt = widget_button(bb, value='Compare to ADT', uname='btn_compare_adt')
                                
    wm_base = widget_base(roi_base, /row, /frame)
    lbl_sd_fa = widget_label(wm_base, value="Whole Brain (seed all voxels with FA greater than):", uname='lbl_sd_fa')
    txt_sd_fa = widget_text(wm_base, value='0.25', xsize=7, uname='txt_sd_fa', /editable)
    btn_sd_fa_est = widget_button(wm_base, value="Estimate # of Seeds", uname='btn_sd_fa_est', $
                                  event_pro='mas_hardi_tractography_gui_estimate')
    btn_sd_fa_trk = widget_button(wm_base, value="Track", uname='btn_sd_fa_trk', $
                                event_pro='mas_hardi_tractography_gui_trackfa')
    b = widget_base(roi_base, /column, /frame)
    l = widget_label(b, value="Tractography using 3D volume ROIs in NIFTI format:")
    for roi = 1, 5 do begin
    
        base_roi_1 = widget_base(b, /row)
        i = string(roi, format='(I0)')
        
        non_ex_base = widget_base(base_roi_1, /nonexclusive)
        btn = widget_button(non_ex_base, value='3D ROI #'+i+': ', $
            uname='btn_roi_enable_'+i, sensitive=1, $
            event_pro='mas_hardi_tractography_gui_roi_enable')
            
        txt = widget_text(base_roi_1, xsize=50, uname='txt_roi_path_'+i)
        
;        dl  = widget_droplist(base_roi_1, value=['Seeding+Analysis', 'Seeding Only', 'Analysis Only'], $
;                              uname='dl_roi_type_'+i, sensitive=0)

        dl  = widget_droplist(base_roi_1, value=['Seeding'], $
                              uname='dl_roi_type_'+i, sensitive=0)
                             
        btn = widget_button(base_roi_1, sensitive=0, $
            uname='btn_roi_choose_'+i, value='Choose...', $
            event_pro='mas_hardi_tractography_gui_roi_choose')
                        
    endfor
    
    base_roi_1 = widget_base(b, /row)
    
    non_ex_base = widget_base(base_roi_1, /nonexclusive)
    btn = widget_button(non_ex_base, value='Termination ROI: ', $
        uname='btn_roi_enable_bound', sensitive=1, $
            event_pro='mas_hardi_tractography_gui_roi_enable')
        
    txt = widget_text(base_roi_1, xsize=60, uname='txt_roi_path_bound')
    
    btn = widget_button(base_roi_1, sensitive=0, $
        uname='btn_roi_choose_bound', value='Choose...', $
        event_pro='mas_hardi_tractography_gui_roi_choose')

    ;;;;;;;;;;;
    
    base_roi_1 = widget_base(b, /row)
    
    non_ex_base = widget_base(base_roi_1, /nonexclusive)
    btn = widget_button(non_ex_base, value='Exclusion ROI: ', $
        uname='btn_roi_enable_exclusion', sensitive=1, $
            event_pro='mas_hardi_tractography_gui_roi_enable')
        
    txt = widget_text(base_roi_1, xsize=60, uname='txt_roi_path_exclusion')
    
    btn = widget_button(base_roi_1, sensitive=0, $
        uname='btn_roi_choose_exclusion', value='Choose...', $
        event_pro='mas_hardi_tractography_gui_roi_choose')
    
    ;;;;;;;;;;;
    
    btn_base = widget_base(b, /row)
    btn_sd_roi_est = widget_button(btn_base, value="Estimate # of Seeds", uname='btn_sd_roi_est', $
                                  event_pro='mas_hardi_tractography_gui_estimate')
    btn_sd_roi_trk = widget_button(btn_base, value="Track", uname='btn_sd_roi_trk', $
                                event_pro='mas_hardi_tractography_gui_trackroi')
        
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Visualization Parameteres Tab
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    bv = widget_base(tab_visual, column=2, /align_left)
    b_main = widget_base(bv, /column, /align_left)
    b = widget_base(b_main, /column, /align_left, /frame)
    lbl_fib_file = widget_label(b, value='Fiber Data File:', /align_left, uname='')
    b1 = widget_base(b, /row)
    txt_fib_file = widget_Text(b1, value='', xsize=75, uname='txt_fibfile')
    btn_fib_choose = widget_button(b1, value="Choose", uname='btn_fibfile_choose', $
                                  event_pro='mas_hardi_tractography_gui_setfibfile')
    
    b = widget_base(b_main, /column, /align_left, /frame)
    lbl_roi_list = widget_label(b, value='Analysis ROIs to load:')
    list_anal_rois = widget_list(b, xsize=80, ysize=13, uname='list_anal_rois', value=['[None]'], uvalue=ptr_new(['[None]']))
    btn_base = widget_base(b, /row)
    btn_add_anal_roi = widget_button(btn_base, value='Add...', uname='btn_add_anal_roi', $
                        event_pro='mas_hardi_tractography_gui_anroi_manage')
    btn_rem_anal_roi = widget_button(btn_base, value='Remove', uname='btn_rem_anal_roi', $
                        event_pro='mas_hardi_tractography_gui_anroi_manage')
    btn_launch_vis = widget_button(btn_base, value='Launch Visualization', event_pro='mas_hardi_tractography_gui_launchvis')
 
    br = widget_base(bv, /column, /align_left)
    lbl = widget_label(br, value="TrackVis Fiber Flipping:")
    b = widget_base(br, /nonexclusive, /column, /frame)
    btn_flip_x = widget_button(b, value="Flip fibers' X axis", uname="btn_flip_x")
    btn_flip_y = widget_button(b, value="Flip fibers' Y axis", uname="btn_flip_y")
    btn_flip_z = widget_button(b, value="Flip fibers' Z axis", uname="btn_flip_z")
    ;boot_btn = widget_button(btn_base, value='Do not press', event_pro='mas_hardi_tractography_gui_bootstrap')
    
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    widget_control, base, /realize
    xmanager, 'mas_hardi_tractography_GUI', base, /no_block
    
    
end
