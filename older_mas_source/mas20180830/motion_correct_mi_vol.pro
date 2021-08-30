function motion_correct_mi_vol_size_reduce, vol, new_size

    old_size = size(vol, /dimensions)
    
    tmp = fltarr(new_size[0], new_size[1], old_size[2]) 
    for i = 0, old_size[2]-1 do begin
        tmp[*,*,i] = congrid(reform(vol[*,*,i]), new_size[0], new_size[1], cubic=-0.5)
    endfor
    outvol = fltarr(new_size)
    for i = 0, new_size[0]-1 do begin
        outvol[i,*,*] = congrid(reform(tmp[i,*,*]), new_size[1], new_size[2], cubic=-0.5)
    endfor

    return, outvol

end

; Subroutine name: motion_correct_mi_vol_objective
; Created by: Bill Triplett
; Calling Information:
;
;   psrc_data: Pointer to a 3d array that will be used as the source data
;              for interpolation
;
;   sz_src: the size of the source data
;
;   interp_ind: [x,y,z] coordinates representing the interpolated data points
;
;   orig_ind: [x,y,z] coordinates representing the untransformed data points
;
;   pvol_out: pointer to hold the interpolated data
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Filled in the interpolated data array with the data values interpolated using 
; the original indices and the transformed indices
;
; Editing Information:

function motion_correct_mi_vol_interp, psrc_data, sz_src, interp_ind, orig_ind, pvol_out, tricubic=tricubic

    k=0L
    st = systime(1)
        
    for x = 0, sz_src[0]-1 do begin
        for y = 0, sz_src[1]-1 do begin
            for z = 0, sz_src[2]-1 do begin
                inn = interp_ind[*,k++]
                if (inn[0] lt 0 or inn[0] gt sz_src[0]-1) then continue
                if (inn[1] lt 0 or inn[1] gt sz_src[1]-1) then continue
                if (inn[2] lt 0 or inn[2] gt sz_src[2]-1) then continue
                (*pvol_out)[x,y,z] = mas_tricubic_3d(inn, data_ptr=psrc_data)
            endfor
        endfor
    endfor
    ;;print, "Elapsed Time:", systime(1)-st

    return, pvol_out

end

; Subroutine name: motion_correct_mi_vol_objective
; Created by: Bill Triplett
; Calling Information:
;
;    tx_vec: vector of transformation data provided by the optimizer
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;   This function is called by an optimizing procedure and should not be
;   called manually
;
; Editing Information:

function motion_correct_mi_vol_objective, tx_vec

    common mi_vol, work

    center_voxel_sm = work.center_voxel_sm 
    sz_ref_sm = work.sz_ref_sm
    
    t3d, /reset
    ;;t3d, scale=replicate(sm_factor, 3)
    t3d, translate=-center_voxel_sm
    t3d, translate=[tx_vec[0],tx_vec[1],tx_vec[2]]
    !P.t = !P.t # [ [ 1.0, tx_vec[3], tx_vec[4], 0.0 ], $
                    [ 0.0,       1.0, tx_vec[5], 0.0 ], $
                    [ 0.0,       0.0,       1.0, 0.0 ], $
                    [ 0.0,       0.0,       0.0, 1.0 ] ]
    t3d, scale=[1.0, 1.0, 1.0] + [tx_vec[6], tx_vec[7], tx_vec[8]]
    t3d, rotate=[tx_vec[9],tx_vec[10],tx_vec[11]]
    t3d, translate=center_voxel_sm
    ;;t3d, scale=replicate(1.0/sm_factor, 3)
    
    ind_sm_tx = vert_t3d(*work.indices_in_sm, matrix=invert(!P.t), /no_divide)
    
    ;vol_out_sm = ptr_new(interpolate(*work.inp_vol_sm, ind_sm_tx[0,*], ind_sm_tx[1,*], ind_sm_tx[2,*]))
    ;*vol_out_sm = reform(*vol_out_sm, sz_ref_sm[0], sz_ref_sm[1], sz_ref_sm[2], /overwrite)

    *work.vol_out_sm = reform(interpolate(*work.inp_vol_sm, ind_sm_tx[0,*], ind_sm_tx[1,*], ind_sm_tx[2,*]), $
                        sz_ref_sm[0], sz_ref_sm[1], sz_ref_sm[2], /overwrite)
    
    ;for i = 0L, n_elements(ind_sm_tx)/3-1 do begin
    ;    pt_or = (*work.indices_in_sm)[*,i]
    ;    pt_tx = (reform(ind_sm_tx[*,i]) > [0.,0.,0.]) < (sz_ref_sm-1)
    ;    (*work.vol_out_sm)[pt_or[0], pt_or[1], pt_or[2]] = (*work.inp_vol_sm)[pt_tx[0], pt_tx[1], pt_tx[2]]
    ;endfor
    
    mi = mc_compute_mi(*work.ref_vol_sm, *work.vol_out_sm, work.mi_binsize)
    ;;;ptr_free, vol_out_sm
    
    return, -mi
        
end

; Subroutine name: motion_correct_mi_makeindices
; Created by: Bill Triplett
; Calling Information:
;
;  dimensions: the dimensions of the array
;
;  zyx: set this keyword for the indices to be looped in z,y,x order rather than
;       x,y,z
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
;  Makes [x,y,z] indices for an array of given dimensions. 
;
; Editing Information:

function motion_correct_mi_makeindices, dimensions, zyx=zyx

    ind = lonarr(3, long(product(dimensions)))
    k = 0L

    if (keyword_set(zyx)) then begin

        for z = 0, dimensions[2]-1 do begin
            for y = 0, dimensions[1]-1 do begin
                for x = 0, dimensions[0]-1 do begin
                    ind[*,k++] = [x,y,z]
                endfor
            endfor
        endfor

    endif else begin

        for x = 0, dimensions[0]-1 do begin
            for y = 0, dimensions[1]-1 do begin
                for z = 0, dimensions[2]-1 do begin
                    ind[*,k++] = [x,y,z]
                endfor
            endfor
        endfor
        
    endelse

    return, ind
 
end

; Subroutine name: motion_correct_mi_alignvol
; Created by: Bill Triplett
; Calling Information:
;
;  pref_vol:  The reference volume; that is the volume that should not change
;
;  pinp_vol:  The input volume: the volume that should be changed to be in alignment
;             with pref_vol
;
;  mi_binsize: The hisogram bin size for the mutial information computation.
;
;  center_voxel: The [x,y,z] location of the center voxel. Rotations and shears will be done
;                using this voxel as the center
;
;  use_work_struct: not implemented
;
; Purpose of subroutine:
; applies a homogenous coordinate transformation matrix to data_ptr (or state1 of project.ci).
;
; Editing Information:

function motion_correct_mi_alignvol, pref_vol=pref_vol, $
                                           pinp_vol=pinp_vol, $
                                           mi_binsize=mi_binsize, $
                                           center_voxel=center_voxel, $
                                           size_reduce=size_reduce, $
                                           use_work_struct=use_work_struct

    common mi_vol, work
    
    work = create_struct(name='MOTION_CORRECT_MI_VOL')
    
    sz_inp = size(*pinp_vol, /dimensions)
    sz_ref = size(*pref_vol, /dimensions)
    
    if (not keyword_set(center_voxel)) then begin
        center_voxel = sz_inp/2.0
    endif
    
    if (keyword_set(size_reduce)) then begin
        pref_vol_sm = ptr_new(congrid(*pref_vol, sz_ref[0]*size_reduce, sz_ref[1]*size_reduce, sz_ref[2]*size_reduce))
        pinp_vol_sm = ptr_new(congrid(*pinp_vol, sz_inp[0]*size_reduce, sz_inp[1]*size_reduce, sz_inp[2]*size_reduce))
    endif else begin
        pref_vol_sm = pref_vol
        pinp_vol_sm = pinp_vol
        size_reduce = 1.0
    endelse

    work.vol_out_sm      = ptr_new(fltarr(sz_ref), /no_copy)    
    work.use_transformations = replicate(1.0, 12)
    work.reduce_factor = size_reduce
    work.mi_binsize = keyword_set(mi_binsize) ? mi_binsize : 1
    work.ref_vol_sm = pref_vol_sm
    work.inp_vol_sm = pinp_vol_sm
    work.sz_ref_sm = sz_inp * size_reduce
    ;;work.sz_inp_sm = sz_inp
    work.center_voxel_sm = center_voxel * size_reduce
    work.indices_in_sm = ptr_new(motion_correct_mi_makeindices(work.sz_ref_sm, /zyx))

    dir_vec    = float([1.0, 1.0, 1.0, $ ;; trans
                         1.0, 1.0, 1.0, $ ;; shear
                         1.0, 1.0, 1.0, $ ;; scale perterbation
                         1.0, 1.0, 1.0]) * work.use_transformations
    dir_vec_wk  = diag_matrix(dir_vec)
    
    start_tx_sm = float([0,0,0,0,0,0,0,0,0,0,0,0])
    
    powell, start_tx_sm, dir_vec_wk, 1e-4, max_mi,$
             'motion_correct_mi_vol_objective', ITER=n_iter
    
    end_trans  = start_tx_sm[0:2]
    end_scale  = start_tx_sm[6:8] + [1.,1.,1.] 
    end_shear  = start_tx_sm[3:5]
    end_rotate = start_tx_sm[9:11]
    
    t3d, /reset
    t3d, scale=replicate(work.reduce_factor,3)
    t3d, translate=-center_voxel
    t3d, translate=end_trans
    !P.t = !P.t # [ [ 1.0, end_shear[0], end_shear[1], 0.0 ], $
                    [ 0.0,          1.0, end_shear[2], 0.0 ], $
                    [ 0.0,          0.0,          1.0, 0.0 ], $
                    [ 0.0,          0.0,          0.0, 1.0 ] ]
    t3d, scale=end_scale
    t3d, rotate=end_rotate
    t3d, translate=center_voxel
    t3d, scale=replicate(1.0/work.reduce_factor,3)

    return, !P.t

end

; Subroutine name: motion_correct_mi_vol
; Created by: Bill Triplett
; Calling Information:
;
;  ref_index: The index into the data that should be used as a reference (standard)
;
;  use_translations: 12 element array indicating which transformations are admittable:
;                    [ 0:2 -> translations, 3:5 -> shear, 6:8 -> scale, 9:11 -> rotate ]
;
;  size_reduce: attempt to save time by downsampling the arrays before the computation
;
;  center_voxel: the center voxel with respect to rotations and shears
;
;  mi_binsize: the histogram bin size for the mutual information calculation
;
;  save_intermediate: set this keyword to save each registered volume when complete as a nifti file
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Performs volume-wise registration if a multi-volume dataet using mutual information..
;
; Editing Information:

pro motion_correct_mi_vol, ref_index=ref_index, $
                               use_transformations=use_transformations, $
                               size_reduce=size_reduce, $
                               center_voxel=center_voxel, $
                               mi_binsize=mi_binsize, $
                               save_intermediate=save_intermediate

    common scan_data
    common mi_vol, work
    
    ci = project.ci
    
    ref_idx = -1
    inp_idx =  2
    work = create_struct(name='MOTION_CORRECT_MI_VOL')
    
    mi_binsize = keyword_set(mi_binsize) ? mi_binsize : 1
    print, "MI_BINSIZE: ", mi_binsize
    
    save_intermediate = (keyword_set(save_intermediate)) ? 1 : 0
     
    if n_elements(ref_index) gt 0 then ref_idx = ref_index
    if keyword_set(use_transformations) then begin
        use_tx = use_transformations
    endif else begin
        use_tx = replicate(1.0, 12)
    endelse
     
    nvolumes = project.imndarray[ci].adim

    ;; get the original datas
    if (ref_idx eq -1) then begin
       print, "REF_IDX: "+string(ref_idx)+" (Using average of all volumes as reference)"
       ref_vol = total(*project.dataarray[ci].state1, 4)
    endif else begin
       print, "REF_IDX: "+string(ref_idx)
       ref_vol = reform((*project.dataarray[ci].state1)[*,*,*,ref_idx])
    endelse
    
    if (keyword_set(size_reduce)) then begin
       sm_factor = float(size_reduce)
    endif else begin
       sm_factor = 1.0
    endelse
    print, "SIZE_REDUCE: "+string(sm_factor)
    
    ;; reduce size for mi search
    sz_ref = size(ref_vol, /dimensions)
    if (not keyword_set(center_voxel)) then begin
        center_voxel = sz_ref/2.
    endif
    print, "CENTER_VOXEL: ", center_voxel
    
    if (sm_factor ne 1.0) then begin
        ref_vol_sm = motion_correct_mi_vol_size_reduce(ref_vol, sz_ref*sm_factor)
    endif else begin
        ref_vol_sm = ref_vol
    endelse
    ;ref_vol_sm = congrid(ref_vol, $
    ;                     sz_ref[0]*sm_factor, $
    ;                     sz_ref[1]*sm_factor, $
    ;                     sz_ref[2]*sm_factor)

    sz_ref_sm = size(ref_vol_sm, /dimensions)

    ind    = ptr_new(motion_correct_mi_makeindices(sz_ref), /no_copy)
    ind_sm = ptr_new(motion_correct_mi_makeindices(sz_ref_sm, /zyx), /no_copy)
    
    work.ref_vol         = ptr_new(ref_vol, /no_copy)
    work.sz_ref          = sz_ref
    work.vol_out_sm      = ptr_new(fltarr(sz_ref_sm), /no_copy)
    work.ref_vol_sm      = ptr_new(ref_vol_sm, /no_copy)
    work.sz_ref_sm       = sz_ref_sm
    work.reduce_factor   = sm_factor
    work.center_voxel    = center_voxel
    work.center_voxel_sm = center_voxel * sm_factor
    work.mi_binsize      = mi_binsize
    work.indices_in      = ind
    work.indices_in_sm   = ind_sm
    
    pbar = obj_new('progressbar', color='red')
    pbar->Start
    
    st_total = systime(1)
    
    ;;nvolumes = 4
    for inp_idx = 0, nvolumes-1 do begin
    
        if (inp_idx eq ref_idx) then begin
            (*project.dataarray[ci].state1)[*,*,*,inp_idx] = $
               (*project.dataarray[ci].state1)[*,*,*,inp_idx]
            continue
        endif
        
        inp_vol    = reform((*project.dataarray[ci].state1)[*,*,*,inp_idx])
        
        if (sm_factor ne 1.0) then begin
            inp_vol_sm = motion_correct_mi_vol_size_reduce(inp_vol, sz_ref_sm)
        endif else begin
            inp_vol_sm = inp_vol
        endelse
        
        work.inp_vol    = ptr_new(inp_vol, /no_copy)
        work.inp_vol_sm = ptr_new(inp_vol_sm, /no_copy)
        
        dir_vec    = float([1.0, 1.0, 1.0, $ ;; trans
                             2.0, 2.0, 2.0, $ ;; shear
                             1.5, 1.5, 1.5, $ ;; scale perterbation
                             2.0, 2.0, 2.0]) * use_tx * sm_factor
        dir_vec_wk  = diag_matrix(dir_vec)
        
        start_tx_sm = float([0,0,0,0,0,0,0,0,0,0,0,0])      
        
        st = systime(1)
        
        powell, start_tx_sm, dir_vec_wk, 1e-5, max_mi,$
                 'motion_correct_mi_vol_objective', ITER=n_iter
        
        end_trans  = start_tx_sm[0:2]
        end_scale  = start_tx_sm[6:8] + [1.,1.,1.] 
        end_shear  = start_tx_sm[3:5]
        end_rotate = start_tx_sm[9:11]
        
        t3d, /reset
        t3d, scale=replicate(sm_factor,3)
        t3d, translate=-center_voxel * sm_factor
        t3d, translate=end_trans
        !P.t = !P.t # [ [ 1.0, end_shear[0], end_shear[1], 0.0 ], $
                        [ 0.0,          1.0, end_shear[2], 0.0 ], $
                        [ 0.0,          0.0,          1.0, 0.0 ], $
                        [ 0.0,          0.0,          0.0, 1.0 ] ]
        t3d, scale=end_scale
        t3d, rotate=end_rotate
        t3d, translate=center_voxel * sm_factor
        t3d, scale=replicate(1.0/sm_factor,3)
        
        final_tx = invert(!P.t) ;;tx_matrix)
        tiny = where(abs(final_tx) lt 1e-8, n_tiny)
        if (n_tiny gt 0) then begin
            final_tx[tiny] = 0.0
        endif
        
        print, "--------------------FINAL-----------------------"
        print, "VOL     : ", inp_idx
        print, "ST TRANS: ", start_tx_sm[0:2]
        print, "EN TRANS: ", end_trans
        print, "ST SCALE: ", start_tx_sm[6:8] + [1.,1.,1]
        print, "EN SCALE: ", end_scale
        print, "ST SHEAR: ", start_tx_sm[3:5]
        print, "EN SHEAR: ", end_shear
        print, "ST ROT  : ", start_tx_sm[9:11]
        print, "EN ROT  : ", end_rotate
        print, "MI      : ", -max_mi
        print, "NITER   : ", n_iter
        print, "Matrix  : "
        print, reform(string(final_tx, format='(%"%15.11f")'), 4, 4)
        print, "VOL_TIME:", (systime(1) - st), (systime(1) - st)/60.
        print, "------------------------------------------------"
        
        ind_tx = vert_t3d(*ind, matrix=final_tx, /no_divide)
        
        final_vol = motion_correct_mi_vol_interp(work.inp_vol, sz_ref, ind_tx, *ind, ptr_new(fltarr(sz_ref)))
        (*project.dataarray[ci].state1)[*,*,*,inp_idx] = *final_vol
        
        if (save_intermediate) then begin
            fname = 'mc_mi_vol_interm_'+string(inp_idx, format='(I04)')+'.nii'
            interm_fname = project.current_path + (get_dir_slash())[0] + fname
            mas_export_nifti, data_ptr=final_vol, file_name=interm_fname
            fname = 'mc_mi_vol_interm_'+string(inp_idx, format='(I04)')+'.matrix'
            interm_fname = project.current_path + (get_dir_slash())[0] + fname
            openw, lun, interm_fname, /get_lun
            ;;printf, lun, reform(string(!P.t, format='(%"%14.11f")'), 4, 4)
            printf, lun, reform(string(!P.t, format='(F0.7)'), 4, 4)
            close, lun & free_lun, lun
        endif
        
        ptr_free, final_vol
        ptr_free, [work.inp_vol, work.inp_vol_sm]
        
        pbar->update, float(inp_idx+1)/float(nvolumes) * 100.
        
        if (pbar->CheckCancel()) then begin
            break
        endif
        
    endfor

    pbar->Destroy

    ptr_free, [work.ref_vol, work.ref_vol_sm, work.indices_in, work.indices_in_sm]
    
    print, "TOTAL TIME: ", systime(1)-st_total
    
end

pro motion_correct_mi_vol_GUI_event, event

    name = widget_info(event.id, /uname)
    widget_control, event.top, get_uvalue=state
    
    case name of 
        
        'btn_cancel': begin
            widget_control, event.top, /destroy
        end
        
        'btn_start': begin
        
            txt_binsz = widget_info(event.top, find_by_uname='txt_binsz')
            widget_control, txt_binsz, get_value=val
            (*state).mi_binsize = fix(val)
            
            motion_correct_mi_vol, $
                ref_index=(*state).ref_index, $
                use_transformations=(*state).use_transformations, $
                mi_binsize=(*state).mi_binsize, $
                save_intermediate=(*state).save_intermediate
            
            ptr_free, state
            widget_control, event.top, /destroy
            
        end
        
        'btn_use_translations': begin
            (*state).use_transformations[0:2] = event.select
        end
                                     
        'btn_use_rotations': begin
            (*state).use_transformations[9:11] = event.select
        end

        'btn_use_shears': begin
            (*state).use_transformations[3:5] = event.select
        end

        'btn_use_scale': begin
            (*state).use_transformations[6:8] = event.select
        end

        'btn_save_interm_vol': begin
            (*state).save_intermediate = event.select
        end
        
        'sldr_ref_index': begin
            (*state).ref_index = event.value
        end
        
        'txt_binsz': begin
            widget_control, event.id, get_value=val
            (*state).mi_binsize = fix(val)
         end
         
     else: print, 'Recvd event from unknown widget: '+name
     
   endcase
   
end


pro motion_correct_mi_vol_GUI

    common scan_data
    common common_widgets
    
    base = widget_base(title="Volume-based motion correction", /column, group_leader=WID_BASE_MAIN, /modal)
    
    lbl = widget_label(base, value="Include transformations:")
    base_tx = widget_base(base, /nonexclusive, row=4)
    btn_use_translations = widget_button(base_tx, value="Translations", uname='btn_use_translations')
    btn_use_rotations = widget_button(base_tx, value="Rotations", uname='btn_use_rotations')
    btn_use_shears = widget_button(base_tx, value="Shears", uname='btn_use_shears')
    btn_use_scale = widget_button(base_tx, value="Scaling", uname='btn_use_scale')
    widget_control, btn_use_translations, /set_button
    widget_control, btn_use_rotations, /set_button
    widget_control, btn_use_shears, /set_button
    widget_control, btn_use_scale, /set_button
    
    base_mi_opts = widget_base(base, /column)
    base_binsz = widget_base(base_mi_opts, /row)
    lbl = widget_label(base_binsz, value="MI Binsize:")
    txt_binsz = widget_text(base_binsz, xsize=5, value='2', uname='txt_binsz', /editable)
        
    base_extra = widget_base(base, /column)
    max_vol = project.imndarray[project.ci].adim
    sel_vol = project.procpramarray[project.ci].adim_start
    sldr_ref_index = widget_slider(base_extra, title="Reference volume", minimum=0, maximum=max_vol-1, value=sel_vol, uname='sldr_ref_index')
    
    base_interp = widget_base(base, row=2)
    lbl = widget_label(base_interp, value="Interpolation: ")
    dl_interp = widget_droplist(base_interp, value=['Tricubic', 'Nearest Neightbor', 'Trilinear'], uname='dl_interp', sensitive=0)

    base_save_opts = widget_base(base, /column, /nonexclusive)
    btn_save_interm_vol = widget_button(base_save_opts, value="Save intermediate files?", uname='btn_save_interm_vol')
    
    base_buttons = widget_base(base, /row, /align_center)
    btn_cancel = widget_button(base_buttons, value="Cancel", uname='btn_cancel')
    btn_start = widget_button(base_buttons, value="Start", uname='btn_start')
    
    state = ptr_new({ use_transformations: replicate(1.0, 12), $
                       mi_binsize: 2, $
                       ref_index: sel_vol, $
                       save_intermediate: 0, $
                       interp_method: 0 })
    
    widget_control, base, /realize, set_uvalue=state
    
    xmanager, 'motion_correct_mi_vol_GUI', base
    
end


; Subroutine name: motion_correct_mi_vol__define
; Created by: Bill Triplett
; Calling Information:
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; 
;  Work structure definition for volume-based MI registration
;
; Editing Information:

pro motion_correct_mi_vol__define

    struct = { MOTION_CORRECT_MI_VOL, $
               use_transformations: fltarr(12), $
               ref_vol: ptr_new(), $
               inp_vol: ptr_new(), $
               ref_vol_sm: ptr_new(), $
               inp_vol_sm: ptr_new(), $
               vol_out_sm: ptr_new(), $
               sz_ref: fltarr(3), $
               sz_ref_sm: fltarr(3), $
               center_voxel: fltarr(3), $
               center_voxel_sm: fltarr(3), $
               reduce_factor: float(0), $
               indices_in: ptr_new(), $
               indices_in_sm: ptr_new(), $
               mi_binsize: 0L }
               
end

