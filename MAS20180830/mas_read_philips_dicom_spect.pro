;+
; :Description:
;    Takes a parameter structure (or array) and puts the appropriate params 
;    into the MAS imndarray format. NOTE this also sets up state1 with the
;    data.
;
; :Params:
;    params - measured parameters from the scan
;    fid    - the fid data
;
;
;
; :Author: btt
;-
pro mas_philips_spect_to_imnd, params, fid

    common scan_data, project
    
    next_index = project.ni
    imnd = project.imndarray[next_index]
    
    ;;;;;;;;;;;;;;;;;;;; Compute dimensions of data array ;;;;;;;;;;;;;;;;;;;;;

    imnd.spect_nucleus        = '1H'
    imnd.spect_spectral_width = params[0].pix_bw
    imnd.spect_acq_size       = params[0].num_points
    imnd.spect_num_avg        = long(params[0].num_avgs)
    imnd.spect_bf1            = params[0].img_freq
    imnd.spect_acq_time       = 400.0

    imnd.fdim = params[0].num_points
    imnd.adim = n_elements(params)
    imnd.n_avg       = long(params[0].num_avgs)
    imnd.scan_name   = params[0].protocol
    imnd.orientation = [ "----", "----","----" ]
    scantime         = params[0].scan_date + ' ' + params[0].scan_time
    imnd.scan_date = scantime  

    imnd.spect_rep_time=params[0].tr
    imnd.dimensions = 1
        
    imnd.file_Path    = params[0].src_file
    imnd.display_Name = params[0].src_file
    imnd.image_type   = 99
    
    imnd.echo_time    = params[0].te
    
    if (n_elements(uniq(params.ti)) ne 1) then begin
        imnd.inversion_time_ptr = ptr_new(params.ti)
    endif else if (n_elements(uniq(params.te)) ne 1) then begin
        imnd.echo_time_ptr = ptr_new(params.te)
    endif
    
    acq_matrix = ptr_new(diag_matrix(replicate(1.0D,4)))
    imnd.acq_matrix = acq_matrix

    ;; matrix to transform a unit cube centered at the origin to SVS voxel location
    ;; in patient coordinate frame.
    t3d, /reset, $
         rotate=[-params[0].VX_ANG_RL, params[0].VX_ANG_AP, params[0].VX_ANG_FH], $
         scale=[params[0].VX_SZ_RL, params[0].VX_SZ_AP, params[0].VX_SZ_FH ], $
         matrix=mat
    mat[3,*] = [-params[0].VX_OFC_RL, params[0].VX_OFC_AP, params[0].VX_OFC_FH, 1.0 ]
    ptr_free, imnd.acq_matrix
    imnd.acq_matrix = ptr_new(mat)
    ;;abs(invert(*project.imndarray[0].sform_matrix) ## *project.imndarray[1].acq_matrix)
    
    imnd.state1_load_procedure = 'mas_philips_dicom_spect_load_state1'
    
    ni = project.ni
    project.imndarray[ni] = imnd
    project.scan_Open_Flag = 1
    
    project.ci = project.ni
    project.ni = project.ni + 1
    
    mas_add_current_scan_to_scan_selection
    
    project.dataarray[project.ci].state1 = ptr_new(fid, /no_copy)
    project.procpramarray[project.ci].state_1 = 1
    
end

;+
; :Description:
;    Dummy procedure to load state 1, this is all done when the scan
;    is opened.
;
; :Author: wtriplett
;-
pro mas_philips_dicom_spect_load_state1

    message, /info, 'Called, but nothing to do.'
    
end

;
;
;+
; :Description:
;    Read a Philips Spectroscopy DICOM file.
;
; :Params:
;    file - input file name. If empty, a file select
;           box will open.
;
; :Keywords:
;    params - on output, this keyword will have the spectroscopy information.
;             if there are most than one FID in the file, then params will be
;             an array of structres of type SPECTROSCOPY_DATA (read_spect.pro)
;    fid    - on output this keyword will contain the FID data (num_points x num_fids)
;    valid_file - this keyword will be set to 0 if there was a problem opening the data,
;                 1 otherwise
;
;-
pro read_philips_dicom_spect, file, fid=fid, params=params, valid_file=valid_file

    valid_file = 0
    
    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        message, /info, !error_state.MSG
        help, /traceback
        if (obj_valid(odicom)) then obj_destroy, odicom
        valid_file = 0
        return
    endif
    ;catch, /cancel
    if (n_elements(file) eq 0) then begin
        file = dialog_pickfile(/must_exist, /read)
        if (file eq '') then return
    endif
    
    mas_dicom_get_dicom_obj, dicom_obj=odicom
    
    result = odicom->Read(file)
    
    if (result eq 0) then begin
        message, /info, 'Could not read DICOM file ('+file+').', /noprint
        obj_destroy, odicom
        return
    endif
    
    ;; If this tag is present, then the DICOM file doesn't have
    ;; acquired data. This seems to indicate that the file
    ;; contains only parameter data.
    p_non_data = odicom->GetValue('2005'x, '1132'x)
    if (ptr_valid(p_non_data)) then begin
        message, /info, 'This DICOM (parameter dicom) file contains no pertinent data.', /noprint
        obj_destroy, odicom
        return
    endif
    
    p_seq_name = odicom->GetValue('0018'x, '9005'x)
    if (ptr_valid(p_seq_name)) then begin
      if (*p_seq_name[0] ne 'SPECTROSCOPY') then begin
        message, /info, 'DICOM file does not seem to be spectroscopy', /noprint
        obj_destroy, odicom
        return
      endif
    endif else begin
      message, /info, 'Could not find pulse sequence name.', /noprint
      obj_destroy, odicom
      return      
    endelse
    
    ;; read relavent values
    p_num_pts  = odicom->GetValue('0018'x, '9127'x) ;; # of complex points
    p_img_freq = odicom->GetValue('2001'x, '1083'x) ;; Img Freq.
    p_pix_bw   = odicom->GetValue('2005'x, '1357'x) ;; not sure about this
    p_tr       = odicom->GetValue('0018'x, '0080'x) ;; TR
    p_te       = odicom->GetValue('2005'x, '1310'x) ;; TE
    p_ti       = odicom->GetValue('2005'x, '1312'x) ;; TI ;2005,1312
    p_proto_nm = odicom->GetValue('0018'x, '1030'x) ;; Protocol Name
    p_subj_id  = odicom->GetValue('0010'x, '0020'x) ;; Patient ID
    p_site_name = odicom->GetValue('0008'x, '0080'x) ;; Institution Name
    p_vx_an_AP = odicom->GetValue('2005'x, '1054'x) ;; Voxel Angulation AP (deg)
    p_vx_an_FH = odicom->GetValue('2005'x, '1055'x) ;; Voxel Angulation FH (deg)
    p_vx_an_RL = odicom->GetValue('2005'x, '1056'x) ;; Voxel Angulation RL (deg)
    p_vx_sz_AP = odicom->GetValue('2005'x, '1057'x) ;; Voxel Size AP (mm)
    p_vx_sz_FH = odicom->GetValue('2005'x, '1058'x) ;; Voxel Size FH (mm)
    p_vx_sz_RL = odicom->GetValue('2005'x, '1059'x) ;; Voxel Size RL (mm)
    p_vx_of_AP = odicom->GetValue('2005'x, '105a'x) ;; Voxel Offset AP (P+mm)
    p_vx_of_FH = odicom->GetValue('2005'x, '105b'x) ;; Voxel Offset FH (H+mm)
    p_vx_of_RL = odicom->GetValue('2005'x, '105c'x) ;; Voxel Offset RL (L+mm)
    p_num_avgs = odicom->GetValue('2001'x, '1088'x) ;; Number of Averages
    p_num_fids = odicom->GetValue('2001'x, '1081'x) ;; Total # fids
    p_fid_num  = odicom->GetValue('2005'x, '1313'x) ;; This fid's number
    p_series_num = odicom->GetValue('0020'x, '0011'x) ;; series number
    p_instance_num = p_fid_num ;; odicom->GetValue('0020'x, '0013'x) ;; instance number
    p_acqu_num = odicom->GetValue('2005'x, '1316'x) ;; Acquisition Number
    p_phase_rows = odicom->GetValue('0018'x,'9095'x) ;; SpectroscopyAcquisitionPhaseRows
    p_phase_cols = odicom->GetValue('0018'x,'9234'x) ;; SpectroscopyAcquisitionPhaseColumns
    
    ;; casting - note not all of the above are set here
    num_fid_in_file = n_elements(p_te)
    params = replicate(create_struct(name='SPECTROSCOPY_DATA'), num_fid_in_file)
    
    for f = 0, num_fid_in_file-1 do begin

        params[f].manufacturer = 'PHILIPS'
        params[f].src_file     = file
        params[f].series_num   = (*p_series_num[0])[0]
        params[f].instance_num = (*p_instance_num[0])[0]
        params[f].acqu_num     = (*p_acqu_num[0])[0]
        params[f].src_type     = 'DICOM'
        params[f].subject_id   = strtrim((*p_subj_id[0])[0],2)
        params[f].site_name    = strtrim((*p_site_name[0])[0],2)
        params[f].num_points   = (long(*p_num_pts[0]))[0]
        params[f].img_freq     = (double(*p_img_freq[0]))[0]
        params[f].pix_bw       = (float(*p_pix_bw[0]))[0]
        params[f].tr           = (float(*p_tr[0]))[0]
        params[f].te           = (float(*p_te[f]))[0]
        params[f].ti           = (float(*p_ti[f]))[0]
        params[f].protocol     = strtrim((*p_proto_nm[0])[0], 2)
        params[f].vx_ang_AP    = (float(*p_vx_an_AP[0]))[0]
        params[f].vx_ang_FH    = (float(*p_vx_an_FH[0]))[0]
        params[f].vx_ang_RL    = (float(*p_vx_an_RL[0]))[0]
        params[f].vx_sz_AP     = (float(*p_vx_sz_AP[0]))[0]
        params[f].vx_sz_FH     = (float(*p_vx_sz_FH[0]))[0]
        params[f].vx_sz_RL     = (float(*p_vx_sz_RL[0]))[0]
        params[f].vx_ofc_AP    = (float(*p_vx_of_AP[0]))[0]
        params[f].vx_ofc_FH    = (float(*p_vx_of_FH[0]))[0]
        params[f].vx_ofc_RL    = (float(*p_vx_of_RL[0]))[0]
        
        params[f].num_avgs     = (long(*p_num_avgs[0]))[0]
        params[f].num_xpts     = (long(*p_phase_cols[0]))[0]
        params[f].num_ypts     = (long(*p_phase_rows[0]))[0]
        
        p_scandate = odicom->GetValue('0008'x,'0021'x) ;; Series Instance Date
        tmp_scandate = strtrim(string(*(p_scandate)[0],2))
        tmp_scandate = strjoin([strmid(tmp_scandate, 0, 4), $
                                strmid(tmp_scandate, 4, 2), $
                                strmid(tmp_scandate, 6, 2) ], '-')
        params[f].scan_date = tmp_scandate
    
        p_scantime = odicom->GetValue('0008'x,'0031'x) ;; Series instance Time
        tmp_scantime = strtrim(string(*(p_scantime)[0],2))
        tmp_scantime = strjoin([strmid(tmp_scantime, 0, 2), $
                                strmid(tmp_scantime, 2, 2), $
                                strmid(tmp_scantime, 4, 2) ], ':')
        params[f].scan_time = tmp_scantime
    
    endfor

    ptr_free, [ p_num_pts, p_img_freq, p_pix_bw, p_tr, p_te, p_proto_nm, p_subj_id, $
                p_vx_an_AP, p_vx_an_FH, p_vx_an_RL, p_vx_sz_AP, p_vx_sz_FH, p_vx_sz_RL, $
                p_vx_of_AP, p_vx_of_FH, p_vx_of_RL, p_num_avgs, p_num_fids, p_fid_num, $
                p_series_num, p_instance_num, p_acqu_num, p_phase_rows, p_phase_cols ]

    bad_dicom = 1B
    p_fid_temp = odicom->GetValue('2005'x, '1270'x) ;; fid as short (2 byte) int.
    if (ptr_valid(p_fid_temp)) then begin
        ;; newer single-fid DICOMS
        ;; Cast short data to float
        num_points = params[0].num_points
        fid_temp   = *(p_fid_temp[0])
        fid_data   = fltarr(num_points*2)
        for p = 0, num_points*2-1 do begin
            fid_data[p] = float(fid_temp,p*4)
        endfor
        ;; returned as output variable
        fid   = complex(fid_data[0:*:2], fid_data[1:*:2])
        bad_dicom = 0
    endif else begin
        ;; Older multi-fid DICOMs.
        ;; we will read the data outselves since gdlffdicom has trouble
        ;; with this dicom type.
        offset = odicom->GetOffset('5600'x,'0020'x)
        length = odicom->GetLength('5600'x,'0020'x)
        if (offset ne -1 and length ne -1) then begin
            fid_data = fltarr(length/4) ;; gdlffdicom divides by 2, we need 4
            openr, lun, file, /get_lun
            point_lun, lun, offset
            readu,lun, fid_data
            free_lun, lun
            fid = complex(fid_data[0:*:2], fid_data[1:*:2])
            fid = reform(reform(fid, params[0].num_points, num_fid_in_file))
            ;; second reform eliminates a [npoints, 1] situation
            bad_dicom = 0
        endif
    endelse
    
    ;; check to see if we're good, otherwise bail.
    if (bad_dicom eq 1) then begin
        message, /info, 'DICOM file does not seem to be spectroscopy', /noprint
        obj_destroy, odicom
        return
    endif

    valid_file = 1
    
    obj_destroy, odicom
    
end


pro read_philips_dicom_spect_dir_event, event

    widget_control, event.top, get_uvalue=uvalue
    
    case event.id of
    
        (*uvalue).btn_ok: begin

            select = widget_info((*uvalue).table, /table_select)
            select = reform(select[1,*])
            select = select[uniq(select[sort(select)])]
            
            (*uvalue).selected_rows[*] = 0
            (*uvalue).selected_rows[select] = 1
            (*uvalue).do_import = 1
            
            widget_control, event.top, /destroy
            
        end
        
        (*uvalue).btn_cancel: begin
        
            (*uvalue).selected_rows[*] = 0
            widget_control, event.top, /destroy
        
        end
        
        (*uvalue).table: begin
            ;;help, event, /str
            if ((event.type eq 9 or event.type eq 4) && event.sel_top ge 0) then begin

                ;; extract just the rows that are selected, regardless
                ;; of whatever columns are selected. selection will be
                ;; expanded to cover all columns in each row.
                select = widget_info((*uvalue).table, /table_select)
                n_cols = n_elements(widget_info((*uvalue).table, /column_widths))
                select = reform(select[1,*])
                select = select[uniq(select[sort(select)])]
                ;; the trick is to adjust the selection and then to keep the view the same.
                table_view = widget_info((*uvalue).table, /table_view)
                
                if (event.type eq 9) then begin
                    exclude = lindgen((event.sel_bottom - event.sel_top + 1)) + event.sel_bottom
                    ;;print, exclude 
                endif
                
                for c = 0, n_elements(select)-1 do begin
                
                    if (event.type eq 9 && select[c] eq exclude) then continue
                    
                    ;; fill the selection out so that entire rows appear to be selected.
                    if ( n_elements(new_select) eq 0 ) then begin
                        new_select = transpose([[lindgen(n_cols)], [replicate(select[c],n_cols)]])
                    endif else begin
                        new_select = [ [new_select], [transpose([[lindgen(n_cols)], [replicate(select[c],n_cols)]])] ]
                    endelse
                endfor

                widget_control, (*uvalue).table, set_table_select=new_select
                widget_control, (*uvalue).table, set_table_view=table_view
            
            endif
        
        end
        
        else: 
        
     endcase
     
end

;+
; :Description:
;    Opens and parsed a DICOM directory thought to contain philips spectroscopy.
;
; :Params:
;    dir - the directory to open, will ask if not provided
;
; :Keywords:
;    params - a structure (or array of structures) containing parameters from the experimet
;    fid    - the [npoints,nfids] array of measured data
;
; :Author: wtriplett
;-
pro read_philips_dicom_spect_dir, dir, params=params, fid=fid
    
    common scan_data, project
    
    if (n_elements(dir) eq 0) then begin
        dir = dialog_pickfile(/directory, /must_exist, /read, path=project.current_path)
        if (dir eq '') then return
    endif
    
    files = file_search(dir, 'XX_*', count=n_files)
    if (n_files eq 0) then begin
        void = dialog_message(['No DICOM files found in selected directory', $
                               dir], /center)
        return
    endif
    
    params_allfiles = replicate(create_struct(name='SPECTROSCOPY_DATA'), n_files)
    num_valid = 0L
    valid = 0L
    
    progressbar = obj_new('progressbar', text='Reading directory...', color='Red')
    progressbar->Start
    
    for i = 0L, n_files-1 do begin
    
        ;; the fid is discarded each time, since only the parameters are of interest
        ;; at the moment. Later, the fids will be read when the selection is made.
        read_philips_dicom_spect, files[i], fid=junk, params=p_tmp, valid_file=vf
        
        if (i mod 20 eq 0) then begin
            progressbar->Update, 100.0*float(i)/n_files
        endif
        
        if (vf eq 0) then continue

        if (p_tmp[0].src_type eq 'SPAR/SDAT') then begin
            message, /info, 'Excluding '+p_tmp[0].src_file+' (Use single file method.)'
            continue
        endif

        if (n_elements(p_tmp) gt 1) then begin
            message, /info, 'Excluding '+p_tmp[0].src_file+' (Multiple FID file not supported. Use single file method.)'
            continue
        endif
         
        params_allfiles[num_valid] = p_tmp[0]
        num_valid++
        
    endfor

    progressbar->Destroy

    if (num_valid eq 0) then begin
        message, /info, 'No eligible spectroscopy files found in '+dir
        return
    endif
    
    params_allfiles = params_allfiles[0:num_valid-1]
    
    if (not keyword_set(title)) then title='Spectroscopy Directory Import'
    
    device, get_screen_size=scr_size
    base       = widget_base(title=title, row=3)
    lbl        = widget_label(base, value='Select the rows to import, and click OK to proceed', /align_center)
    table      = widget_table(base, value=params_allfiles, column_labels=tag_names(params_allfiles), $
                              /disjoint, /all_events, /resizeable_columns, scr_ysize=scr_size[1]*0.80)
    btn_base   = widget_base(base, /row, /align_center, /base_align_center)
    btn_ok     = widget_button(btn_base, value='Import Selected')
    btn_cancel = widget_button(btn_base, value='Cancel')
    
    uvalue = ptr_new({ selected_rows: bytarr(n_elements(params_allfiles)), $
                       do_import: 0B, $
                       table: table, $
                       base: base, $
                       btn_ok: btn_ok, $
                       btn_cancel: btn_cancel })
                       
    widget_control, base, /realize, set_uvalue=uvalue
    xmanager, 'read_philips_dicom_spect_dir', base
    
    if ((*uvalue).do_import eq 0) then begin
        ;; the cancel button was pressed, or closebox was clicked.
        ptr_free, uvalue
        return
    endif

    selected_rows = where((*uvalue).selected_rows ne 0, n_selected)
    ptr_free, uvalue
    
    if (n_selected eq 0) then return
    
    params = params_allfiles[selected_rows]
    fid = complexarr(params[0].num_points, n_selected)
    
    for i = 0, n_selected-1 do begin
    
        read_philips_dicom_spect, params[i].src_file, fid=fid_tmp, params=param_tmp
        
        if (n_elements(param_tmp) eq 1) then begin
            fid[*,i] = fid_tmp
        endif

    endfor
    
    ;; Sort by instance number, which appears to be consistent with the ordering
    ;; in the spar file. Only sort if they have different instance numbers.
    ;; This is because some Siemens data may be read in the correct order, but
    ;; since all the instance numbers are the same, the sort() function will 
    ;; order them strangely.
    temp = params.instance_num
    if (n_elements(temp) gt 1) then begin
        if (total(temp[1:*] - temp[0:*]) ne 0) then begin
            sorted = sort(params.instance_num)
            params = params[sorted]
            fid    = fid[*,sorted]
        endif
    endif
    
    project.current_path = dir
    mas_philips_spect_to_imnd, params, fid
    
    
end

;; main data structre for spectroscopy parameters.
pro spectroscopy_data__define

        struct = { SPECTROSCOPY_DATA, $
                   src_file     : '', $
                   series_num   : 0L, $
                   instance_num : 0L, $
                   acqu_num     : 0L, $
                   src_type     : '', $
                   manufacturer : '', $
                   subject_id   : '', $
                   site_name    : '', $
                   scan_date    : '', $
                   scan_time    : '', $
                   num_points   : long(0), $  
                   num_xpts     : long(0), $
                   num_ypts     : long(0), $ 
                   img_freq     : double(0), $ ;; MHz
                   pix_bw       : float(0), $  ;; Hz
                   tr           : float(0), $  ;; ms
                   te           : float(0), $  ;; ms
                   ti           : float(0), $  ;; ms
                   num_avgs     : long(0), $
                   protocol     : '', $
                   vx_ang_AP    : float(0), $ ;; deg
                   vx_ang_FH    : float(0), $ ;; deg
                   vx_ang_RL    : float(0), $ ;; deg
                   vx_ofc_AP    : float(0), $ ;; mm
                   vx_ofc_FH    : float(0), $ ;; mm
                   vx_ofc_RL    : float(0), $ ;; mm
                   vx_sz_AP     : float(0), $  ;; mm
                   vx_sz_FH     : float(0), $  ;; mm
                   vx_sz_RL     : float(0), $  ;; mm
                   phase0_corr  : float(0), $
                   phase1_corr  : float(0)  }
end
