;; $Id$
;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Subroutine name: mnf_bmatrix2tx
; Created by: BT, 2008-10
; Calling Information:
;
;   NUM:  The index in the b_matrix array of the individual b_matrix
;         that should be returned
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
;    Returns a properly formatted b-matrix from the MAS project. Note
;    that this entails dividing the off-diagonal elements by 2 since
;    it is a mas convention to double them when they are read.
;
; Editing Information:

;function mnf_bmatrix2tx, num
;
;    common scan_data
;    forward_function mnf_bmatrix2grad
;
;    b = *project.imndarray[project.ci].b_matrix
;    b[3:5,*] /= 2.0
;
;    ii = num
;
;    bbb = [ [ b[0,ii], b[3,ii], b[4,ii] ], $
;        [ b[3,ii], b[1,ii], b[5,ii] ], $
;        [ b[4,ii], b[5,ii], b[2,ii] ] ]
;    bbb /= (*project.imndarray[project.ci].bval_array)[num]
;
;    eval = eigenql(bbb, eigenvectors=evecs)
;
;    tx = fltarr(4,4)
;    tx[0:2, 0:2] = evecs
;    tx[3,3] = 1.0
;
;    return, tx
;
;end

; Subroutine name: mnf_show
; Created by: BT, 2008-10
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:

;pro mnf_show, grad
;
;    grad_sph = cv_coord(from_rect=grad, /to_sphere, /degrees)
;    adt_show_gradient_profile, gradients=[grad_sph[0,*], 90-grad_sph[1,*]]
;
;end

; Subroutine name: mnf_bmatrix2grad
; Created by: BT, 2008-10
; Calling Information:
;
;   SPHERE: set this keyword to make the function return the gradients
;           in spherical coordinates.
;
;   degrees: set this keyword to have the function return the
;            gradients in degrees rather than in radians
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
;   Recovers the gradient directions used to create the b-matrix in
;   MAS.
;
; Editing Information:

;function mnf_bmatrix2grad, sphere=sphere, degrees=degrees
;
;    common scan_data
;    ci = project.ci
;
;    b_matrix = *project.imndarray[ci].b_matrix
;    b_matrix[3:5,*] /= 2.0
;
;    bvals = *project.imndarray[ci].bval_array
;
;    nbvals = n_elements(bvals)
;
;    grad = fltarr(3, nbvals)
;
;    for ii = 0, nbvals-1 do begin
;
;        bi_mat = [ [ b_matrix[0,ii], b_matrix[3,ii], b_matrix[4, ii] ], $
;            [ b_matrix[3,ii], b_matrix[1,ii], b_matrix[5, ii] ], $
;            [ b_matrix[4,ii], b_matrix[5,ii], b_matrix[2, ii] ] ]
;
;        evals = eigenql(bi_mat, eigenvectors=evecs)
;
;
;        grad[*,ii] = reform(evecs[0,*])
;
;    endfor
;
;    if (keyword_set(sphere)) then begin
;        degrees = keyword_set(degrees) ? 1 : 0
;        grad = cv_coord(from_rect=grad, /to_sphere, degrees=degrees)
;    endif
;
;    return, grad
;
;end

; Subroutine name: mnf_qtcompose
; Created by: BT, 2008-10
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:

function mnf_qtcompose, axis, phi

    if n_params() EQ 0 then begin
        info = 1
        USAGE:
        message, 'USAGE:', /info
        message, 'Q = QTCOMPOSE(AXIS, PHI)', info=1
        return, 0
    endif
    
    nph = n_elements(phi)
    nv  = n_elements(axis)/3
    if nph LT 1 OR nv LT 1 then goto, USAGE
    
    nq  = nv > nph
    q = make_array(value=axis[0]*phi[0]*0., 4,nq)
    
    sph = sin(phi/2) & cph = cos(phi/2)
    if nph EQ 1 AND nv EQ 1 then return, [ axis[0:2] * sph[0], cph[0] ]
    if nph GT 1 AND nv EQ 1 then begin
        ;; Single axis, multiple rotation angles
        q[0,*] = axis[0]*sph
        q[1,*] = axis[1]*sph
        q[2,*] = axis[2]*sph
        q[3,*] = cph
    endif else if nph EQ 1 AND nv GT 1 then begin
        ;; Multiple axis, single rotation
        q[0:2,*] = axis*sph[0]
        q[3,*]   = cph[0]
    endif else if nph EQ nv then begin
        ;; Multiple axes, multiple rotations
        q[0:2,*] = axis*rebin(reform(temporary(sph),1,nq),3,nq)
        q[3,*]   = temporary(cph)
    endif else begin
        message, 'ERROR: number of axes and angles do not match'
    endelse
    
    return, q
end


; Subroutine name: mnf_mat2qt
; Created by: BT, 2008-10
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: convert a transformation matrix to quaternion representation.
;
; Editing Information:

function mnf_mat2qt, mat
    amat = transpose(mat)
    r11 = amat[0,0] & r12 = amat[0,1] & r13 = amat[0,2]
    r21 = amat[1,0] & r22 = amat[1,1] & r23 = amat[1,2]
    r31 = amat[2,0] & r32 = amat[2,1] & r33 = amat[2,2]
    a = 0.5 * sqrt(1 + r11 + r22 + r33)
    if (a ne 0) then begin
        b = 0.25 * (r32 - r23)/a
        c = 0.25 * (r13 - r31)/a
        d = 0.25 * (r21 - r12)/a
    endif else begin
        xd = 1.0 + r11 - (r22+r33) ;  /* 4*b*b */
        yd = 1.0 + r22 - (r11+r33) ;  /* 4*c*c */
        zd = 1.0 + r33 - (r11+r22) ;  /* 4*d*d */
        if ( xd gt 1.0 ) then begin
            b = 0.5 * sqrt(xd) ;
            c = 0.25 * (r12+r21) / b ;
            d = 0.25 * (r13+r31) / b ;
            a = 0.25 * (r32-r23) / b ;
        endif else if( yd gt 1.0 ) then begin
            c = 0.5  * sqrt(yd) ;
            b = 0.25 * (r12+r21) / c ;
            d = 0.25 * (r23+r32) / c ;
            a = 0.25 * (r13-r31) / c ;
        endif else begin
            d = 0.5  * sqrt(zd) ;
            b = 0.25 * (r13+r31) / d ;
            c = 0.25 * (r23+r32) / d ;
            a = 0.25 * (r21-r12) / d ;
        endelse
        
        if ( a lt 0.0 ) then begin
            b=-b
            c=-c
            d=-d
            a=-a
        endif
    endelse
    
    return, [a, b, c, d]

end

; Subroutine name: mnf_qt2mat
; Created by: BT, 2008-10
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: convert a quaterion representation into a transformation matrix
;     note that only the bcd part of the [a,b,c,d] is needed.
;
; Editing Information:

function mnf_qt2mat, bcd, invert=invert

    ; THIS IS ADAPTED FROM CHAPTER 12 BY F.L.MARKLEY[]
    if n_params() EQ 0 then begin
        info = 1
        USAGE:
        message, 'USAGE:', /info
        message, 'MATRIX = QTMAT(Q)', info=info
        return, 0
    endif
    
    
    nbcd = n_elements(bcd)/3
    if nbcd LT 1 then goto, USAGE
    
    q = [ sqrt(1.0 - (total(bcd^2))), bcd[0], bcd[1], bcd[2] ]
    
    if NOT keyword_set(invert) then begin
        q1 = q[0,*] & q2 = q[1,*] & q3 = q[2,*] & q4 = q[3,*]
    endif else begin
        q1 = -q[0,*] & q2 = -q[1,*] & q3 = -q[2,*] & q4 = q[3,*]
    endelse
    
    a = q[0] & b = q[1] & c = q[2] & d = q[3]
    
    mat = dblarr(4,4) 
    mat[0,0] = a*a+b*b-c*c-d*d
    mat[0,1] = 2*b*c+2*a*d
    mat[0,2] = 2*b*d-2*a*c
    mat[1,0] = 2*b*c-2*a*d
    mat[1,1] = a*a+c*c-b*b-d*d
    mat[1,2] = 2*c*d+2*a*b
    mat[2,0] = 2*b*d+2*a*c
    mat[2,1] = 2*c*d-2*a*b
    mat[2,2] = a*a+d*d-c*c-b*b
    mat[3,3] = 1.0
    return, mat[0:2,0:2]


end


;+
; :Description:
;    Performs a NIFTI export using the details indicated in the
;    MAS_EXPORT_NIFTI_EXPORT_SPEC structure. This structure
;    contains information about what should be exported, where to
;    save the files, and whether or not to use compression.
;
; :Params:
;    export_spec - see MAS_EXPORT_NIFTI_EXPORT_SPEC structure definition
;
; :Author: wtriplett
;-
pro mas_export_nifti_run, export_spec

    common scan_data
    
    on_error, 2
    
    ci = project.ci

    SL = path_sep()
    
    file_prefix = export_spec.file_prefix
    use_compression = export_spec.use_compression
    
    if (strpos(file_prefix, SL) ne -1) then begin
        message, 'Filename prefix cannot contain any "'+SL+'"s.'
        return
    endif

    dest_dir = export_spec.dest_dir
    if (dest_dir eq '') then return
    if (not file_test(dest_dir)) then begin
        if (file_test(file_dirname(dest_dir), /write)) then begin
            file_mkdir, dest_dir
        endif else begin
            message, "Unable to create destination directory."
            return
        endelse
    endif
    
    file_path_prefix = dest_dir+SL+file_prefix

    ;; only FLOAT is supported at the moment.
    data_type = 'FLOAT'
  
    if (export_spec.export_scan_data) then begin
        data_ptr = project.dataarray[ci].state1
        if (ptr_valid(data_ptr)) then begin
            message, /info, "Exporting NIFTI: Scan Data"
            mas_export_nifti, data_ptr=data_ptr, data_type=data_type, $
                file_name=file_path_prefix+'data.nii', compress=use_compression
        endif
    endif
    
    if (export_spec.export_bvals) then begin
        message, /info, "Exporting NIFTI: b-balues"
        bvals = *project.imndarray[ci].bval_array
        openw, lun, file_path_prefix+'bvals.txt', /get_lun
        printf, lun, strjoin(string(bvals, format='(%"%0.3f")'), '  ')
        close, lun
        free_lun, lun
    endif
    
    if (export_spec.export_bvecs) then begin
        message, /info, "Exporting NIFTI: gradient directions"
        data_dims = size(*project.dataarray[ci].state1)
        gr = fltarr(3, n_elements(*project.imndarray[ci].angle_theta))
        gr[0,*] = *project.imndarray[ci].angle_phi         ;; longitude
        gr[1,*] = 90. - *project.imndarray[ci].angle_theta ;; latitude
        gr[2,*] = 1.0
        
        gr_rect = cv_coord(from_sphere=gr, /to_rect, /degrees)
        acq_matrix = *project.imndarray[ci].acq_matrix

        if (n_elements(gr_rect)/3 lt n_elements(bvals)) then begin
            gr_rect = transpose([ [replicate(gr_rect[0], n_elements(bvals))], $
                                  [replicate(gr_rect[1], n_elements(bvals))], $
                                  [replicate(gr_rect[2], n_elements(bvals))] ])
        endif             
        gr_rect = vert_t3d(gr_rect, matrix=acq_matrix, /no_copy)
        
        openw, lun, dest_dir+SL+file_prefix+'bvecs.txt', /get_lun
        printf, lun, strjoin(string(reform(gr_rect[0,*]), format='(%"%0.6f")'), '  ')
        printf, lun, strjoin(string(reform(gr_rect[1,*]), format='(%"%0.6f")'), '  ')
        printf, lun, strjoin(string(reform(gr_rect[2,*]), format='(%"%0.6f")'), '  ')
        close, lun
        free_lun, lun
    endif
    
    if (export_spec.export_te) then begin
        message, /info, "Exporting NIFTI: arrayed echo times"
        if (ptr_valid(project.imndarray[ci].echo_time_ptr)) then begin
            tes = *project.imndarray[ci].echo_time_ptr
        endif else begin
            tes = [ project.imndarray[ci].echo_time ]
        endelse
        openw, lun, file_path_prefix+'te.txt', /get_lun
        printf, lun, strjoin(string(tes, format='(%"%0.3f")'), '  ')
        close, lun
        free_lun, lun
    endif

    if (export_spec.export_ti) then begin
        message, /info, "Exporting NIFTI: arrayed inversion times"
        if (ptr_valid(project.imndarray[ci].inversion_time_ptr)) then begin
            tis = *project.imndarray[ci].inversion_time_ptr
        endif else begin
            tis = [ project.imndarray[ci].inversion_time ]
        endelse
        openw, lun, file_path_prefix+'ti.txt', /get_lun
        printf, lun, strjoin(string(tis, format='(%"%0.3f")'), '  ')
        close, lun
        free_lun, lun
    endif

    if (export_spec.export_tr) then begin
        message, /info, "Exporting NIFTI: arrayed repetition times"
        if (ptr_valid(project.imndarray[ci].rep_time_ptr)) then begin
            trs = *project.imndarray[ci].rep_time_ptr
        endif else begin
            trs = [ project.imndarray[ci].recov_time ]
        endelse
        openw, lun, file_path_prefix+'tr.txt', /get_lun
        printf, lun, strjoin(string(trs, format='(%"%0.3f")'), '  ')
        close, lun
        free_lun, lun
    endif

    if (export_spec.export_ad) then begin
        data_ptr = project.dataarray[ci].avg_dif
        if (ptr_valid(data_ptr)) then begin
            message, /info, "Exporting NIFTI: AD"
            mas_export_nifti, data_ptr=data_ptr, data_type=data_type, $
                              file_name=file_path_prefix+'AD.nii', compress=use_compression
        endif
    endif

    if (export_spec.export_fa) then begin
        data_ptr = project.dataarray[ci].frac_ani
        if (ptr_valid(data_ptr)) then begin
            message, /info, "Exporting NIFTI: FA"
            mas_export_nifti, data_ptr=data_ptr, data_type=data_type, $
                              file_name=file_path_prefix+'FA.nii', compress=use_compression
        endif
    endif

    if (export_spec.export_s0) then begin
        data_ptr = ptr_new(reform((*project.dataarray[ci].adt)[*,*,*,0]))
        if (ptr_valid(data_ptr)) then begin
            message, /info, "Exporting NIFTI: S0"
            mas_export_nifti, data_ptr=data_ptr, data_type=data_type, $
                              file_name=file_path_prefix+'S0.nii', compress=use_compression
            ptr_free, data_ptr
        endif
    endif

    if (export_spec.export_v1) then begin
        data_ptr = project.dataarray[ci].eign_vec
        if (ptr_valid(data_ptr)) then begin
            message, /info, "Exporting NIFTI: V1"
            vec = reform((*data_ptr)[*,*,*,*,0])
            mas_export_nifti, data_ptr=ptr_new(vec), file_name=file_path_prefix+'V1.nii', $
                              compress=use_compression
        endif
    endif
    
    if (export_spec.export_v2) then begin
        data_ptr = project.dataarray[ci].eign_vec
        if (ptr_valid(data_ptr)) then begin
            message, /info, "Exporting NIFTI: V2"
            vec = reform((*data_ptr)[*,*,*,*,1])
            mas_export_nifti, data_ptr=ptr_new(vec), file_name=file_path_prefix+'V2.nii', $
                              compress=use_compression
        endif
    endif

    if (export_spec.export_v3) then begin
        data_ptr = project.dataarray[ci].eign_vec
        if (ptr_valid(data_ptr)) then begin
            message, /info, "Exporting NIFTI: V3"
            vec = reform((*data_ptr)[*,*,*,*,2])
            mas_export_nifti, data_ptr=ptr_new(vec), file_name=file_path_prefix+'V3.nii', $
                              compress=use_compression
        endif
    endif

    if (export_spec.export_e1) then begin
        data_ptr = project.dataarray[ci].eign_val
        if (ptr_valid(data_ptr)) then begin
            message, /info, "Exporting NIFTI: E1"
            vec = reform((*data_ptr)[*,*,*,0])
            mas_export_nifti, data_ptr=ptr_new(vec), data_type=data_type, $
                              file_name=file_path_prefix+'E1.nii', compress=use_compression
        endif
    endif

    if (export_spec.export_e2) then begin
        data_ptr = project.dataarray[ci].eign_val
        if (ptr_valid(data_ptr)) then begin
            message, /info, "Exporting NIFTI: E2"
            vec = reform((*data_ptr)[*,*,*,1])
            mas_export_nifti, data_ptr=ptr_new(vec), data_type=data_type, $
                              file_name=file_path_prefix+'E2.nii', compress=use_compression
        endif
    endif

    if (export_spec.export_e3) then begin
        data_ptr = project.dataarray[ci].eign_val
        if (ptr_valid(data_ptr)) then begin
            message, /info, "Exporting NIFTI: E3"
            vec = reform((*data_ptr)[*,*,*,2])
            mas_export_nifti, data_ptr=ptr_new(vec), data_type=data_type, $
                              file_name=file_path_prefix+'E3.nii', compress=use_compression
        endif
    endif

    if (export_spec.export_curvefit) then begin
        temp = *project.dataarray[ci].imgfit1
        sz = size(temp)
        if (sz[0] eq 3) then begin
            ;; Curvefit is [n, xdim, ydim]
            data = fltarr(sz[2], sz[3], 1, sz[1])
            for v = 0, sz[1]-1 do begin
                data[*,*,0,v] = temp[v,*,*]
            endfor
        endif else begin
            ;; curvefit is 2D [ xdim, ydim ]
            data = temporary(temp)
        endelse
        data_ptr = ptr_new(data, /no_copy)
        if (ptr_valid(data_ptr)) then begin
            message, /info, "Exporting NIFTI: curvefit.nii"
            mas_export_nifti, data_ptr=data_ptr, data_type=data_type, $
                              file_name=file_path_prefix+'curvefit.nii', compress=use_compression
        endif
        ptr_free, data_ptr
        
    endif

    if (export_spec.export_sig_enhance) then begin
        data_ptr = project.dataarray[ci].signal_enhancement
        if (ptr_valid(data_ptr)) then begin
            message, /info, "Exporting NIFTI: sig_enhance.nii"
            mas_export_nifti, data_ptr=data_ptr, data_type=data_type, $
                              file_name=file_path_prefix+'sig_enhance.nii', compress=use_compression
        endif
    endif
    
end

; Subroutine name: mas_export_nifti_event
; Created by: BT, 2008-10
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
;    Event handler for the nifti export choice window.
;
; Editing Information:

pro mas_export_nifti_event, event

    common scan_data
    
    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        void = dialog_message(['An error has ocurred:', $
                               !ERROR_STATE.MSG], /error, /center)
        return
    endif
    
    ;; user wants to go or cancel...
    name = widget_info(event.id, /uname)
    if (name eq 'nif_cancel') then begin
        widget_control, event.top, /destroy
        return
    endif else if (name ne 'nif_export') then begin
        return
    endif
    
    ci = project.ci
        
    SL = (get_dir_slash())[1]

    ;; fill in export_specification with user-supplied info from GUI
    
    export_spec = create_struct(name='MAS_EXPORT_NIFTI_EXPORT_SPEC')
    
    txt_file_prefix = widget_info(event.top, find_by_uname='txt_file_prefix')
    widget_control, txt_file_prefix, get_value=file_prefix
    export_spec.file_prefix = strtrim(file_prefix, 2)

    btn_use_compression = widget_info(event.top, find_by_uname='btn_use_compression')
    export_spec.use_compression = widget_info(btn_use_compression, /button_set)
    
    if (strpos(file_prefix, SL) ne -1) then begin
        void = dialog_message('Filename prefix cannot contain any "'+SL+'"s.', /center, /error)
        return
    endif

    dest_dir = dialog_pickfile(path=project.current_path, /directory, /write, title="Select the destination directory:")
    if (dest_dir eq '') then return
    if (not file_test(dest_dir)) then begin
        if (file_test(file_dirname(dest_dir), /write)) then begin
            file_mkdir, dest_dir
        endif else begin
            void = dialog_message("Unable to create destination directory.", /error, /center)
            return
        endelse
    endif
    export_spec.dest_dir = dest_dir

    export_spec.data_type = 'FLOAT'
    ;; NOTE NOTE This widget is deactivated, only FLOAT is supported.
    dl_datatype = widget_info(event.top, find_by_uname='nif_datatype')
    if (widget_info(dl_datatype, /valid_id)) then begin
        datatype_val = widget_info(dl_datatype, /droplist_select)
        case datatype_val of
            0: data_type = 'FLOAT'
            1: data_type = 'SHORT'
            else: data_type = 'FLOAT'
        endcase
    endif
    
    id = widget_info(event.top, find_by_uname='nif_data')
    export_spec.export_scan_data = widget_info(id, /button_set)

    id = widget_info(event.top, find_by_uname='nif_bval')
    export_spec.export_bvals = widget_info(id, /button_set)
    
    id = widget_info(event.top, find_by_uname='nif_bvec')
    export_spec.export_bvecs = widget_info(id, /button_set)
    
    id = widget_info(event.top, find_by_uname='nif_te')
    export_spec.export_te = widget_info(id, /button_set)

    id = widget_info(event.top, find_by_uname='nif_ti')
    export_spec.export_ti = widget_info(id, /button_set)

    id = widget_info(event.top, find_by_uname='nif_tr')
    export_spec.export_tr = widget_info(id, /button_set)

    id = widget_info(event.top, find_by_uname='nif_ad')
    export_spec.export_ad = widget_info(id, /button_set)
    
    id = widget_info(event.top, find_by_uname='nif_fa')
    export_spec.export_fa = widget_info(id, /button_set)

    id = widget_info(event.top, find_by_uname='nif_s0')
    export_spec.export_s0 = widget_info(id, /button_set)

    id = widget_info(event.top, find_by_uname='nif_v1')
    export_spec.export_v1 = widget_info(id, /button_set)
    
    id = widget_info(event.top, find_by_uname='nif_v2')
    export_spec.export_v2 = widget_info(id, /button_set)
    
    id = widget_info(event.top, find_by_uname='nif_v3')
    export_spec.export_v3 = widget_info(id, /button_set)
    
    id = widget_info(event.top, find_by_uname='nif_e1')
    export_spec.export_e1 = widget_info(id, /button_set)

    id = widget_info(event.top, find_by_uname='nif_e2')
    export_spec.export_e2 = widget_info(id, /button_set)

    id = widget_info(event.top, find_by_uname='nif_e3')
    export_spec.export_e3 = widget_info(id, /button_set)
        
    id = widget_info(event.top, find_by_uname='nif_curvefit')
    export_spec.export_curvefit = widget_info(id, /button_set)
    
    id = widget_info(event.top, find_by_uname='nif_sigenhance')
    export_spec.export_sig_enhance = widget_info(id, /button_set)

    mas_export_nifti_run, export_spec
    
    widget_control, event.top, /destroy
    
end

;+
; :Description:
;    Structure definition for the NIFTI export specification.
;    
;    This is used by mas_export_nifti_run to determine what to export
;    and where to place the files.
;
; :Author: wtriplett
;-
pro mas_export_nifti_export_spec__define

    struct = { MAS_EXPORT_NIFTI_EXPORT_SPEC, $
               file_prefix: '', $       ;; string that is prepended to each file that is exported
               dest_dir: '', $          ;; The destination direction where the files will be placed (created)
               data_type: '', $         ;; the data type to export the data (not currently used; always FLOAT)
               use_compression: 0B, $   ;; set = 1 to use gz compression, 0 uncompressed
               export_scan_data: 0B, $  ;; set = 1 to export state1 scan data
               export_bvals: 0B, $      ;; set = 1 to export diffusion bvals
               export_bvecs: 0B, $      ;; set = 1 to export diffusion gradient vecs.
               export_te: 0B, $         ;; set = 1 to export (possibly arrayed) echo times
               export_ti: 0B, $         ;; set = 1 to export (possibly arrayed) repetition times
               export_tr: 0B, $         ;; set = 1 to export (possibly arrayed) inversion times
               export_ad: 0B, $         ;; set = 1 to export DTI "FA" data
               export_fa: 0B, $         ;; set = 1 to export DTI "AD" data
               export_s0: 0B, $         ;; set = 1 to export DTI "S0" data
               export_v1: 0B, $         ;; set = 1 to export DTI 1st eigenvector data
               export_v2: 0B, $         ;; set = 1 to export DTI 2nd eigenvector data
               export_v3: 0B, $         ;; set = 1 to export DTI 3rd eigenvector data
               export_e1: 0B, $         ;; set = 1 to export DTI 1st eigenvalue data
               export_e2: 0B, $         ;; set = 1 to export DTI 2nd eigenvalue data
               export_e3: 0B, $         ;; set = 1 to export DTI 3rd eigenvalue data
               export_curvefit: 0B, $   ;; set = 1 to export curvefit data
               export_sig_enhance: 0B } ;; set = 1 to export DCE signal enhancement data
               
end


; Subroutine name: mas_export_nifti
; Created by: BT, 2008-10
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:
pro mas_export_nifti, $
    diffusion=diffusion, $
    data_ptr=data_ptr, $
    file_name=file_name, $
    data_type=data_type, $
    compress=compress, $
    is_mask=is_mask

    common scan_data, project
    forward_function get_dir_slash
    
    ci = project.ci
    
    if (keyword_set(data_ptr)) then begin
        if (not ptr_valid(data_ptr)) then begin
            junk = dialog_message('Data pointer is not a valid pointer.', /error, /center)
            return
        endif
    endif else begin
        mas_load_state_1
        data_ptr = project.dataarray[ci].state1
    endelse

    nif = create_struct(name='MAS_NIFTI_HEADER')
        
    if (keyword_set(is_mask)) then begin
        data_type = 'BYTE'
    endif else begin
        is_mask = 0
    endelse
    
    if (not keyword_set(data_type)) then begin
        data_type = 'FLOAT'
    endif else begin
        data_type = strupcase(data_type)
        if (data_type ne 'FLOAT' and data_type ne 'SHORT' and data_type ne 'BYTE') then begin
            print, "Unsupported data type: "+data_type+". Using FLOAT."
            data_type = 'FLOAT'
        endif
    endelse
    
    if (not keyword_set(file_name)) then begin
        dest = dialog_pickfile(path=project.current_path, /directory, /write)
        if (dest eq '') then return
    endif else begin
        dest = file_name
    endelse
    
    if (keyword_set(compress)) then begin
        if (stregex(dest, "\.gz$", /boolean) eq 0) then begin
            dest += ".gz"
        endif
    endif
    
    dest_dir = file_dirname(dest)
    
    if (file_test(dest)) then begin
        if (file_test(dest, /write)) then begin
            void = dialog_message(['File exists:', $
                dest, $
                'Do you want to overwrite?'], /question, /center)
            if (void eq 'No') then return
        endif else begin
            void = dialog_message('Unable to save file, permission denied.', /error, /center)
            return
        endelse
    endif else if (not file_test(dest_dir, /write)) then begin
        void = dialog_message('Unable to save file, permission denied.', /error, /center)
        return
    endif
        
    openw, lun, dest, /get_lun, error=err, compress=compress
    
    if (err ne 0) then begin
        void = dialog_message('Unable to save file, permission denied.', /error, /center)
        return
    endif
  
    ;; Construct NIFTI header
    ;; hdr size
    nif.sizeof_hdr = 348
    
    ;; dim info
    nif.dim_info = 0
    
    ;; data dim
    data_dims = size(*data_ptr)
    tmp = intarr(8) + 1
    tmp[0:data_dims[0]] = data_dims[0:data_dims[0]]
    nif.dim = tmp
    
    ;; data type (16=float), (4=short)
    tmp = intarr(1)
    case data_type of
        'FLOAT': tmp[*] = 16
        'SHORT': tmp[*] = 4
        'BYTE':  tmp[*] = 2
        else :  tmp[*] = 16
    endcase
    nif.datatype = tmp

    ;; bits / pix (32=float), (16=short)
    case data_type of
        'FLOAT': tmp[*] = 32
        'SHORT': tmp[*] = 16
        'BYTE':  tmp[*] = 8
        else :  tmp[*] = 32
    endcase
    nif.bitpix = tmp
    
    ;; slice start
    tmp[*] = 0
    nif.slice_start = tmp
        
    ;; vox data offset
    tmp = float(352)
    nif.vox_offset = tmp
    
    ;; scale slope
    if (data_type eq 'SHORT' and not is_mask) then begin
        scale_fac = float(32766.0)/float(max(*data_ptr))
    endif else if (data_type eq 'BYTE' and not is_mask) then begin
        scale_fac = float(255.0)/float(max(*data_ptr))
    endif else begin
        scale_fac = float(1)
    endelse
    tmp = 1.0/scale_fac
    nif.scl_slope = tmp
    
    ;; scale intercept
    tmp = float(0)
    nif.scl_inter = tmp
    
    ;; xyzt units
    tmp = byte(10)
    nif.xyzt_units = tmp
    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    voxel_dim = [ project.imndarray[ci].f_voxsz/project.procpramarray[ci].freq_interp , $
                  project.imndarray[ci].p_voxsz/project.procpramarray[ci].phase_interp, $
                  project.imndarray[ci].s_voxsz/project.procpramarray[ci].slice_interp ]
    
    if (project.procpramarray[ci].have_orientation eq 1) then begin
        orient_mat = mas_orient_get_matrix([ project.procpramarray[ci].orientation_Imin, $
                                             project.procpramarray[ci].orientation_Jmin, $
                                             project.procpramarray[ci].orientation_Kmin ])
    endif else begin
        orient_mat = diag_matrix([-1.0,1.0,1.0])
    endelse
    
    determ_mat = determ(orient_mat)
    qmat = orient_mat
    if (determ_mat lt 0) then begin
        qmat[2,*] = -qmat[2,*]
    endif
    quatern = mnf_mat2qt(qmat)
    quatern = quatern[1:3]

    ;; qform code
    if (ptr_valid(project.imndarray[ci].qform_bcd)) then begin
        tmp = 1B
    endif else begin
        tmp = byte(1)
    endelse
    nif.qform_code = tmp
    
    ;; quatern_bcd
    if (ptr_valid(project.imndarray[ci].qform_bcd)) then begin
        quatern = *project.imndarray[ci].qform_bcd
    endif
    nif.quatern_b = quatern[0]
    nif.quatern_c = quatern[1]
    nif.quatern_d = quatern[2]
    
    ;; qoffset_xyz
    if (ptr_valid(project.imndarray[ci].qoffset_xyz)) then begin
        qoffset = *project.imndarray[ci].qoffset_xyz
    endif else begin
        fovs = [project.imndarray[ci].f_fov, $
                project.imndarray[ci].p_fov, $
                project.imndarray[ci].s_fov] # orient_mat * 10.
        qoffset = [-fovs[0]/2., -fovs[1]/2., -fovs[2]/2.]
    endelse
    nif.qoffset_x = qoffset[0]
    nif.qoffset_y = qoffset[1]
    nif.qoffset_z = qoffset[2]
    
    ;; sform code
    tmp = byte(1)
    nif.sform_code = tmp
    
    if (ptr_valid(project.imndarray[ci].sform_matrix)) then begin
        sform_scl = diag_matrix([ 1.0/project.procpramarray[ci].freq_interp , $
                                  1.0/project.procpramarray[ci].phase_interp, $
                                  1.0/project.procpramarray[ci].slice_interp ])
        sform_mat = *project.imndarray[ci].sform_matrix
        sform_mat[0:2,0:2] = sform_mat[0:2,0:2] # sform_scl
        nif.srow_x = sform_mat[*,0]
        nif.srow_y = sform_mat[*,1]
        nif.srow_z = sform_mat[*,2]
        
    endif else begin
        ;; Sform matrix
        sform_mat = diag_matrix(voxel_dim*10.) # orient_mat
        fovs = [project.imndarray[ci].f_fov, $
                project.imndarray[ci].p_fov, $
                project.imndarray[ci].s_fov] # orient_mat * 10.
        nif.srow_x = [sform_mat[*,0], -fovs[0]/2.]
        nif.srow_y = [sform_mat[*,1], -fovs[1]/2.]
        nif.srow_z = [sform_mat[*,2], -fovs[2]/2.]
    endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        
    tmp = fltarr(8)+1.0
    tmp[0] = (determ_mat)*1.0 ;; "qfac"
    tmp[1] = voxel_dim[0]*10. ;; -> mm
    tmp[2] = voxel_dim[1]*10.
    tmp[3] = voxel_dim[2]*10.
    nif.pixdim = tmp

    ;; magic
    tmp = 'n+1'
    nif.magic = byte(tmp)
    
    ;; voxel data
    writeu, lun, nif
    point_lun, lun, 352
    
    case data_type of
        'BYTE': begin
            if (is_mask) then begin
                tmp = fix(*data_ptr, type=1)
                high = where(tmp gt 0, n_high, complement=low, ncomplement=n_low)
                if (n_high ne 0) then begin
                    tmp[high] = 1
                endif
            endif else begin
                tmp = fix((*data_ptr * scale_fac), type=1)
            endelse
            writeu, lun, tmp
        end
        'SHORT': begin
            print, "DATA SCALE FACTOR for "+file_name+": ", scale_fac
            tmp = fix((*data_ptr * scale_fac), type=2)
            writeu, lun, tmp
        end
        'FLOAT': begin
        
            progressbar = obj_new('progressbar', $
                                  title="Writing NIFTI Data...",$
                                  text="Saving: "+(file_basename(dest))[0], $
                                  /fast_loop)
            progressbar->start
            
            tmpdim = long(nif.dim)
            blksize = tmpdim[1]*tmpdim[2]
            nvoxels = (tmpdim[1]*tmpdim[2]*tmpdim[3]*( (tmpdim[4] eq 0L) ? 1L : tmpdim[4]) )
            for i = 0L, nvoxels-1, blksize do begin
                writeu, lun, fix((*data_ptr)[i:i+blksize-1], type=4)
                progressbar->update, float(i)/float(nvoxels) * 100.0
            endfor
            
            progressbar->destroy
            ;writeu, lun, fix((*data_ptr), type=4)
            
        end
        
        else: writeu, lun, fix((*data_ptr), type=4)
    endcase
    
    close, lun
    free_lun, lun
    
    print, "End of export..."

end

; Subroutine name: mas_export_nifti_gui
; Created by: BT, 2008-10
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
;    Creates the nifti export GUI
;
; Editing Information:

pro mas_export_nifti_gui, DICOM=DICOM

    common scan_data
    common common_widgets
    
    ci = project.ci
    
    if (keyword_set(DICOM)) then begin
        title = 'DICOM export'
    endif else begin
        title = 'NIFTI export'
        if (project.procpramarray[ci].have_orientation eq 0) then begin
            yesno = dialog_message(['This data set is currently unoriented. Would you like to orient it before exporting?',$
                                    '(Orienting the data set will improve interoperability with other software.)'], $
                                    /question, /center)
            if (yesno eq 'Yes') then begin
                mas_orient_data
            endif
        endif
    endelse
    
    base = widget_base(group_leader=WID_BASE_MAIN, title=title, /column, /modal)
    
    ; base_lbl = widget_base(base, /align_center)
    ; lbl = widget_label(base_lbl, value="Select the data to export:", /align_center)
    
    base_2COL = widget_base(base, /row)
    base_L = widget_base(base_2COL, /column)
    base_R = widget_base(base_2COL, /column)
    
    ;; Scan Data
    base_typegroup = widget_base(base_L, /frame, /column)
    lbl = widget_label(base_typegroup, value='Scan Data:')
    
    base_nonex = widget_base(base_typegroup, /column, /align_right, /nonexclusive)
    
    data = widget_button(base_nonex, uname='nif_data', value='Scan Data', $
        sensitive=ptr_valid(project.dataarray[ci].state1))
        
    if (not keyword_set(DICOM)) then begin
        bvec = widget_button(base_nonex, uname='nif_bvec', value="Diffusion Directions", $
            sensitive=ptr_valid(project.imndarray[ci].angle_theta))
            
        bval = widget_button(base_nonex, uname='nif_bval', value="Diffusion Weights", $
            sensitive=ptr_valid(project.imndarray[ci].bval_array))
            
        te = widget_button(base_nonex, uname='nif_te', value="Echo Times", $
            sensitive=1)

        ti = widget_button(base_nonex, uname='nif_ti', value="Inversion Times", $
            sensitive=ptr_valid(project.imndarray[ci].inversion_time_ptr))

        tr = widget_button(base_nonex, uname='nif_tr', value="Repetition Times", $
            sensitive=1)
            
        ;; Misc. Fitted Data
        base_typegroup = widget_base(base_L, /frame, /column)
        lbl = widget_label(base_typegroup, value='Misc. Data:')
        
        base_nonex = widget_base(base_typegroup, /column, /align_right, /nonexclusive)
        
        curvefit = widget_button(base_nonex, uname='nif_curvefit', value='Curve Fit', $
            sensitive=ptr_valid(project.dataarray[ci].imgFit1))
            
        sigenhance = widget_button(base_nonex, uname='nif_sigenhance', value="Signal Enhancement", $
            sensitive=ptr_valid(project.dataarray[ci].signal_Enhancement))
            
        base_datatype = widget_base(base_L, /frame, /row)
        lbl = widget_label(base_datatype, value='Datatype: ')
        dl = widget_droplist(base_datatype, sensitive=0, value=['32 Bit Float', '16 Bit Short'], uname='nif_datatype')
        
    endif
    
    ;; Diffusion Data
    base_typegroup = widget_base(base_R, /frame, /column)
    lbl = widget_label(base_typegroup, value="Processed Diffusion Data:")
    
    base_nonex = widget_base(base_typegroup, /column, /align_right, /nonexclusive)
    
    fa = widget_button(base_nonex, uname='nif_fa', value='Fractional Anisotropy', $
        sensitive=ptr_valid(project.dataarray[ci].frac_ani))
        
    ad = widget_button(base_nonex, uname='nif_ad', value='Average Diffusivity', $
        sensitive=ptr_valid(project.dataarray[ci].avg_dif))

    s0 = widget_button(base_nonex, uname='nif_s0', value='S0 Image', $
        sensitive=ptr_valid(project.dataarray[ci].adt))
        
    if (not keyword_set(DICOM)) then begin
    
        v1 = widget_button(base_nonex, uname='nif_v1', value='1st Eigenvector', $
            sensitive=ptr_valid(project.dataarray[ci].eign_vec))
            
        v2 = widget_button(base_nonex, uname='nif_v2', value='2nd Eigenvector', $
            sensitive=ptr_valid(project.dataarray[ci].eign_vec))
            
        v3 = widget_button(base_nonex, uname='nif_v3', value='3rd Eigenvector', $
            sensitive=ptr_valid(project.dataarray[ci].eign_vec))
            
        e1 = widget_button(base_nonex, uname='nif_e1', value='1st Eigenvalue', $
            sensitive=ptr_valid(project.dataarray[ci].eign_val))
            
        e2 = widget_button(base_nonex, uname='nif_e2', value='2nd Eigenvalue', $
            sensitive=ptr_valid(project.dataarray[ci].eign_val))
            
        e3 = widget_button(base_nonex, uname='nif_e3', value='3rd Eigenvalue', $
            sensitive=ptr_valid(project.dataarray[ci].eign_val))
            
    endif
    
    if (not keyword_set(dicom)) then begin
        base_prefix = widget_base(base, /row, /align_center, /frame)
        lbl = widget_label(base_prefix, value="Prefix for output file names:")
        txt_file_prefix = widget_text(base_prefix, xsize=25, value="", /editable, uname='txt_file_prefix')
    
        base_comp = widget_base(base, /row, /align_center, /nonexclusive)
        btn_use_compression = widget_button(base_comp, value="Use gzip compression", uname="btn_use_compression")
    endif else begin
        base_comp = widget_base(base, /row, /align_center, /nonexclusive)
        btn_include_rescale = widget_button(base_comp, value="Include real-world rescale info", uname="btn_include_rescale")
        widget_control, btn_include_rescale, /set_button
    endelse
    
    
    base_btn = widget_base(base, /row, /align_center)
    btn_export = widget_button(base_btn, uname='nif_export', value="Export")
    btn_cancel = widget_button(base_btn, uname='nif_cancel', value="Cancel")
    
    widget_control, base, /realize
    
    if (keyword_set(DICOM)) then begin
        event_prefix = 'mas_export_dicom'
    endif else begin
        event_prefix = 'mas_export_nifti'
    endelse
    
    xmanager, event_prefix, base
    
end

;+
; :Description:
;    Reads all of the data in a text file and returns it in an array.
;    As the name suggests, it assume the data are in floating point numbers,
;    so if you are expecting ints or something, you can convert it after the fact.
;
; :Params:
;    textfile - input:  a string path to a text file
;    array - output: the resulting data array of floats read from textfile
;
;
; :Author: wtriplett
;-
pro mnf_read_floats_from_text_file, textfile, array

    openr, lun, textfile, /get_lun
    
    line = ''
    
    if (n_elements(array) ne 0) then begin
        junk = temporary(array)
    endif
    
    while (not eof(lun)) do begin
        
        readf, lun, line
        temp = float(strsplit(line, ' ', /extract))
        
        if (n_elements(array) eq 0) then begin
            array = temp
        endif else begin
            array = [ array, temp ]
        endelse
        
    endwhile
    
    free_lun, lun

end

; Subroutine name: mnf_make_data_array
; Created by: BT, 2008-10
; Calling Information:
;
;   NIFTI: a nifti structure
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
;   Given a nifti structure, this function will create the data array
;   of the proper type and size for the scan data.
;
; Editing Information:

function mnf_make_data_array, nifti

    dim = nifti.dim
    type = nifti.datatype
        
    dims = dim[1:dim[0]]
    
    case type of
    
        0:  data_array = intarr(dims)       ;; unknown datatype
        1:  data_array = bytarr(dims)       ;; binary (1 bit/voxel)  (BT: unsure of IDLs format for this)
        2:  data_array = bytarr(dims)       ;; unsigned char (8 bits/voxel)
        4:  data_array = intarr(dims)       ;; signed short (16 bits/voxel)
        8:  data_array = lonarr(dims)       ;; signed int (32 bits/voxel)
        16: data_array = fltarr(dims)       ;; float (32 bits/voxel)
        32: data_array = complex_arr(dims)  ;; complex (64 bits/voxel)
        64: data_array = dblarr(dims)       ;; double (64 bits/voxel)
        256: data_array = bytarr(dims)      ;; signed char (8 bits)
        512: data_array = intarr(dims)      ;; unsigned short (16 bits)
        else: begin
            print, "Data type not supported."
            return, ptr_new()
        end
        
    endcase
    
    return, ptr_new(data_array, /no_copy)
    
end

; Subroutine name: mas_read_nifti
; Created by: BT, 2008-10
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:

function mas_read_nifti, $
        nifti_filename  = nifti_filename, $
        default_dir = default_dir, $
        bval_filename  = bval_filename, $
        bvec_filename  = bvec_filename, $
        tr_filename = tr_filename, $
        ti_filename = ti_filename, $
        te_filename = te_filename, $
        progress  = progress, $
        hdr_only  = hdr_only, $
        read_status=read_status
        
    common scan_data, project
    
    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        void = dialog_message(['An error occurred while reading the NIFTI file:', $
                               !ERROR_STATE.msg], /center, /error)
        read_status = 0
        if (obj_valid(progressbar)) then progressbar->Destroy
        return, -1
    endif
    
    read_status = 0
    
    if (not keyword_set(nifti_filename)) then begin
    
        if (not keyword_set(default_dir)) then begin
            default_dir = n_elements(project) ne 0 ? project.current_path : '.'
        endif
        
        nifti_filename = dialog_pickfile(title='Select an .nii or .nii.gz file', $
            path=default_dir, filter='*.nii;*.nii.gz;*.hdr')
        if (nifti_filename eq '') then begin
            read_status = 0
            return, ptr_new()
        endif
        
    endif else begin
        if (not file_test(nifti_filename, /read)) then begin
            junk = dialog_message('Cannot read: '+nifti_filename, /error, /center)
            read_status = 0
            return, 0
        endif
    endelse
    
    if keyword_set(progress) and not keyword_set(hdr_only) then begin
        progressbar = obj_new('progressbar', title='MAS NIFTI READER', $
            text='Reading NIFTI file...', color='Red')
        progressbar->start
    endif
    
    if (not keyword_set(bval_filename)) then bval_filename = ''
    if (not keyword_set(bvec_filename)) then bvec_filename = ''
    if (not keyword_set(tr_filename)) then tr_filename = ''
    if (not keyword_set(ti_filename)) then ti_filename = ''
    if (not keyword_set(te_filename)) then te_filename = ''
    
    openr, lun, nifti_filename, /get_lun, err=err, /compress
    
    if err ne 0 then begin
        read_status = 0
        return, 0
    endif
    
    read_status = 1
    
    nif_hdr = create_struct(name='MAS_NIFTI_HEADER')
    nifti   = create_struct(name='MAS_NIFTI_FILE')
    readu, lun, nif_hdr
    struct_assign, nif_hdr, nifti, /verbose
    nifti.nifti_filename = nifti_filename
    nifti.bval_filename = bval_filename
    nifti.bvec_filename = bvec_filename
    nifti.tr_filename = tr_filename
    nifti.ti_filename = ti_filename
    nifti.te_filename = te_filename            
    
    if (keyword_set(hdr_only)) then begin
        close, lun
        free_lun, lun
        return, nifti
    endif
    
    pdata = mnf_make_data_array(nifti)
    
    if (strmid(nifti_filename, strlen(nifti_filename)-3) eq 'hdr') then begin
        img_filename = strmid(nifti_filename, 0, strlen(nifti_filename)-3)+'img'
        close, lun
        free_lun, lun
        if (not file_test(img_filename, /read)) then begin
            message, 'Corresponding .img file not found for .hdr file.', /info
            read_status = 0
            return, 0
        endif
        openr, lun, img_filename, /get_lun, err=err
        if (err ne 0) then begin
            read_status = 0
            return, 0
        endif
        
    endif else begin
        point_lun, lun, nifti.vox_offset
    endelse
    
    blksize = long(nifti.dim[1])*long(nifti.dim[2])
    blk = (*pdata)[0:blksize-1] 
    nvoxels = n_elements(*pdata)
    for i = 0L, nvoxels-1, blksize do begin
    
        readu, lun, blk
        
        if (nifti.scl_slope ne 0) then begin
            blk_scaled = (blk * nifti.scl_slope) + nifti.scl_inter
        endif else begin
            blk_scaled = blk
        endelse

        (*pdata)[i:i+blksize-1] = blk_scaled
        if (obj_valid(progressbar)) then begin
            progressbar->update, float(i)/float(nvoxels) * 100.0
        endif
    endfor
        
    nifti.voxel_data = pdata
    
    close, lun
    free_lun, lun
    
    if (obj_valid(progressbar)) then begin
        progressbar->update, 100.0
        progressbar->destroy
    endif
    
    return, nifti
    
end


; Subroutine name: mnf_nifti_hdr_to_mas
; Created by: BT, 2008-10
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: converts information contained in NIFTI header to MAS's imnd array
;
; Editing Information:

pro mnf_nifti_hdr_to_mas, nifti

    forward_function get_dir_slash
    
    common scan_data
    
    error_state = 0
    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        void = dialog_message(['A error ocurred while reading the NIFTI data.', $
            'Please make sure that the NIFTI data is not corrupted, or if you', $
            'Are also importing bval/bvec/tr/ti/te text file, make sure that', $
            'the number of entries in those files match the expected data dimension.'], $
            /error, /center)
        print, !ERROR_STATE.msg
        help, /traceback
        ptr_free, nifti.voxel_data
        return
    endif
    
    ni = project.ni
    
    project.imndArray[ni].fdim = nifti.dim[1]
    project.imndArray[ni].pdim = nifti.dim[2]
    project.imndArray[ni].sdim = nifti.dim[3]
    project.imndArray[ni].adim = nifti.dim[4] > 1
    
    project.imndArray[ni].k_fdim_span = nifti.dim[1]
    project.imndArray[ni].k_pdim_span = nifti.dim[2]
    project.imndArray[ni].k_sdim_span = nifti.dim[3]
    
    project.imndArray[ni].scan_date = 'Unknown'
    project.imndArray[ni].scan_name = 'Unknown Name'
    project.imndArray[ni].recov_time = 0.0
    project.imndArray[ni].echo_time = 0.0
    project.imndArray[ni].n_echo = 0.0
    project.imndArray[ni].n_avg = 1
    
    project.imndArray[ni].pn_avg = PTR_NEW(1)
    
    project.imndArray[ni].f_fov = nifti.dim[1] * nifti.pixdim[1] / 10.
    project.imndArray[ni].p_fov = nifti.dim[2] * nifti.pixdim[2] / 10.
    project.imndArray[ni].s_fov = nifti.dim[3] * nifti.pixdim[3] / 10.
    
    if (nifti.sform_code ne 0) then begin
    
        ;; find the axes permutation matrix for this transformation
        sform_mat = [ [nifti.srow_x], [nifti.srow_y], [nifti.srow_z] ]
        A = sform_mat[0:2,0:2]
        sgn_A = fltarr(size(A, /dimensions))
        
        zro = where(abs(A) lt 1e-3, n_zro, complement=non_zro, ncomplement=n_non_zro)
        if (n_zro gt 0) then begin 
            A[zro] = 0.0
        endif 
        
        if (n_non_zro eq 0) then begin
            void = dialog_message(['Note: SFORM code is nonzero, but SFORM matrix', $
                                   'contains all zeros. This should not happen.', $
                                   'Will set SFORM to identity.'], /error, /center)
            A = diag_matrix([1.0,1.0,1.0,1.0])
        endif
          
        sgn_A[non_zro] = A[non_zro]/(abs(A[non_zro]))
        
        svdc, A, sigma, U, V
        ;; A = U ## diag_matrix(sigma) ## transpose(V)
        print, 'SFORM:'
        print, A
        axes_tx = fltarr(3,3)
        for column = 0, 2 do begin
            tmp = reform(abs(A[column,*]))
            max_ind = where(tmp eq max(tmp))
            axes_tx[column, max_ind[0]] = 1.0
        endfor
        ;;axes_tx = round(U##V)
        ;; polar decomposition to get cardinal axes to determine orientation.
        axes_tx = abs(axes_tx) * sgn_A
        print, 'POLAR_DECOMP:'
        print, axes_tx
        ijk = mas_orient_from_matrix(axes_tx)
        print, 'AXES: ', ijk
        
        project.procPramArray[ni].orientation_Imin = ijk[0]
        project.procPramArray[ni].orientation_Jmin = ijk[1]
        project.procPramArray[ni].orientation_Kmin = ijk[2]
        project.procpramarray[ni].have_orientation = 1
    endif else begin
        axes_tx = diag_matrix(float([1,1,1]))
        sform_mat = fltarr(4,3)
        project.procpramarray[ni].have_orientation = 0
    endelse
    
    project.imndArray[ni].slices = nifti.dim[3]
    project.imndArray[ni].thick = nifti.pixdim[3]
    
    project.imndarray[ni].f_voxsz = nifti.pixdim[1]/10.
    project.imndarray[ni].p_voxsz = nifti.pixdim[2]/10.
    project.imndarray[ni].s_voxsz = nifti.pixdim[3]/10.
    
    project.imndArray[ni].dimensions = nifti.dim[0]
    
    project.imndArray[ni].fat_sup = 0
    
    project.imndArray[ni].file_Path = nifti.nifti_filename
    
    project.imndArray[ni].display_Name = nifti.nifti_filename
    
    project.imndArray[ni].orientation = ['Unknown', '']
    project.imndArray[ni].sform_matrix = ptr_new(sform_mat)
    if (nifti.qform_code ne 0) then begin
        project.imndArray[ni].qform_bcd = ptr_new([nifti.quatern_b, nifti.quatern_c, nifti.quatern_d])
        project.imndArray[ni].qoffset_xyz = ptr_new([nifti.qoffset_x, nifti.qoffset_y, nifti.qoffset_z])
        project.imndArray[ni].qfac = nifti.pixdim[0]
    endif
    
    ;; Note that we do not infer any kind of acq matrix from NIFTI.
    ;; In other words, the diffusion gradients incoming from a text file
    ;; are assumes to be correctly specified.
    acq_matrix = ptr_new(diag_matrix(replicate(1.0D,4)))
    project.imndArray[ni].acq_matrix  = acq_matrix
    
    project.imndArray[ni].image_type = 18
    
    project.current_path = file_dirname(nifti.nifti_filename)
    
    bval_file = nifti.bval_filename
    bvec_file = nifti.bvec_filename
    if (file_test(bval_file, /read)) then begin
    
        ;; read bvals & bvecs
        bval_array = fltarr(nifti.dim[4])
        n_bvals = n_elements(bval_array)
        
        openr, lun, bval_file, /get_lun
        readf, lun, bval_array
        if (not eof(lun)) then begin
            junk = dialog_message(['bval file appears to contain more information ', $
                'than the expected number of b-values.'], $
                /error, /center)
            message, 'bval file appears to contain more information '+ $
                'than the expected number of b-values.'
        endif
        project.imndarray[ni].bval_array = ptr_new(bval_array)
        project.imndarray[ni].n_bvals = n_bvals
        close, lun
        free_lun, lun
    endif
    
    if (file_test(bval_file, /read) and file_test(bvec_file, /read)) then begin
        bvec_array = fltarr(nifti.dim[4],3)
        bvec_temp  = fltarr(nifti.dim[4] * 3)
        openr, lun, bvec_file, /get_lun
        readf, lun, bvec_temp
        bvec_array = reform(bvec_temp, nifti.dim[4], 3)
        close, lun
        free_lun, lun
        
        bvec_array = transpose(bvec_array)
        
        ;; normalize if not already
        for ii = 0, n_bvals-1 do begin
            bvec_norm = norm(bvec_array[*,ii])
            if (abs(bvec_norm) lt 1e-5) then continue
            
            bvec_array[*, ii] /= bvec_norm
            
        endfor
        
        ;; Save the original vectors. 
        bvec_array_save = bvec_array
        bvec_array_sph = cv_coord(from_rect=bvec_array, /to_sphere, /degrees)
        
        ;; adjust for negative phi => [0, 359]
        phi = reform(bvec_array_sph[0,*])
        neg_phi = where(phi lt 0, num_neg_phi)
        if (num_neg_phi gt 0) then begin
            phi[neg_phi] += 360.0
            bvec_array_sph[0,*] = phi
        endif
        
        ;; adjust for IDL's spherical coord convention
        project.imndarray[ni].angle_theta = ptr_new(90. - reform(bvec_array_sph[1,*]))
        project.imndarray[ni].angle_phi   = ptr_new(phi)
        
        ;; These come from the philips 3T assumed values. We need to verify that
        ;; these are still correct.
        project.imndArray[ni].big_delta   = ptr_new(replicate(29.73, n_bvals))
        project.imndArray[ni].small_delta = ptr_new(replicate(16.91, n_bvals))
        
        bvec_normalized = bvec_array_save
        
        bmatrix = fltarr(6, n_bvals)
        for ii = 0, n_bvals-1 do begin
        
            bval = bval_array[ii]
            x_part = bvec_normalized[0,ii]
            y_part = bvec_normalized[1,ii]
            z_part = bvec_normalized[2,ii]
            
            xx = bval * x_part * x_part
            yy = bval * y_part * y_part
            zz = bval * z_part * z_part
            xy = bval * x_part * y_part * 2 ;; The factor of two is a MAS convention
            xz = bval * x_part * z_part * 2
            yz = bval * z_part * y_part * 2
            
            bmatrix[*,ii] = [xx, yy, zz, xy, xz, yz]
            
        endfor
        project.imndarray[ni].image_type = 16
        project.imndarray[ni].b_matrix = ptr_new(bmatrix)
        
    endif
    
    if (file_test(nifti.tr_filename)) then begin
        mnf_read_floats_from_text_file, nifti.tr_filename, tr_array
        n_trs = n_elements(tr_array)
        
        if (n_trs eq 1) then begin
            project.imndarray[ni].recov_time = tr_array[0]
        endif else if (n_trs eq nifti.dim[4]) then begin
            project.imndarray[ni].rep_time_ptr = ptr_new(tr_array)
            project.imndarray[ni].recov_time = tr_array[0]
        endif else begin
            message, 'Number of TRs in file does not match with data dimensions'
        endelse
    endif

    if (file_test(nifti.ti_filename)) then begin
        mnf_read_floats_from_text_file, nifti.ti_filename, ti_array
        n_tis = n_elements(ti_array)
        
        if (n_tis eq 1) then begin
            ;project.imndarray[ni].recov_time = tr_array[0]
        endif else if (n_tis eq nifti.dim[4]) then begin
            project.imndarray[ni].inversion_time_ptr = ptr_new(ti_array)
        endif else begin
            message, 'Number of TIs in file does not match with data dimensions'
        endelse
    endif

    if (file_test(nifti.te_filename)) then begin
        mnf_read_floats_from_text_file, nifti.te_filename, te_array
        n_tes = n_elements(te_array)
        
        if (n_tes eq 1) then begin
            project.imndarray[ni].echo_time = te_array[0]
        endif else if (n_tes eq nifti.dim[4]) then begin
            project.imndarray[ni].echo_time_ptr = ptr_new(te_array)
            project.imndarray[ni].echo_time = te_array[0]
        endif else begin
            message, 'Number of TEs in file does not match with data dimensions'
        endelse
    endif
    
    project.procPramArray[ni].signal_type = 8
    
    project.imndarray[ni].state1_load_procedure = 'mnf_load_state_1'
    
    catch, /cancel
    
    project.scan_Open_Flag = 1
    
    project.ci = project.ni
    
    project.ni = project.ni + 1
    
    mas_add_current_scan_to_scan_selection
    
    update_status_bar, ''
    
end

; Subroutine name: mas_open_nifti_GUI_browse
; Created by: BT, 2008-10
; Calling Information:
;
;   STATE:  the mas nifti GUI state pointer
;
;   NIFTI:  set this keyword to browse for a nifti data file
;
;   BVALS:  set this keyword to browse for a bval text file
;
;   BVECS:  set this keyword to browse for a gradient direction text
;           file
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
;   Locates and performs checks on the various data file associated
;   with a nifti import.
;
; Editing Information:

pro mas_open_nifti_GUI_browse, state, nifti=nifti, bvals=bvals, bvecs=bvecs, other=other, $
                               selected_file=selected_file

    common scan_data
    
    nkeyword = keyword_set(nifti) + keyword_set(bvals) + keyword_set(bvecs)
;    if (nkeyword gt 1) then begin
;        void = dialog_message('Only one NIFTI file may be opened at one time.', /error, /center)
;        return
;    endif else if (nkeyword eq 0) then begin
;        void = dialog_message('A type of file must be specified.', /error, /center)
;        return
;    endif
    
    if ((*state).file_path eq '') then begin
        file_path = project.current_path
    endif else begin
        file_path = (*state).file_path
    endelse
    
    if keyword_set(nifti) then begin
    
        filename = dialog_pickfile(title='Select an .nii or .nii.gz file', $
            path=file_path, filter='*.nii;*.nii.gz;*.hdr')
        if (filename eq '') then begin
            return
        endif
        
        ;; compressed nifti is partially supported, so this is commented out.
        ;if (stregex('\.gz$', filename, /boolean)) then begin
        ;    void = dialog_message(['This file looks to be compressed (.gz).', $
        ;        'MAS can only read uncompressed NIFTI files', $
        ;        'Please uncompress the file and try again.'], $
        ;        /error, /center)
        ;    return
        ;endif
        
        if (not file_test(filename, /read)) then begin
            void = dialog_message(['The file could not be read', $
                'Check the file permissions and try again.'], $
                /error, /center)
            return
        endif
        
        (*state).nifti_filename = filename
        if ((*state).file_path eq '') then begin
            (*state).file_path = file_dirname(filename)
        endif
        selected_file = filename
        return
    endif
    
    if (keyword_set(bvals)) then begin
    
        filename = dialog_pickfile(title='Select the text file containing the b-values', $
            path=file_path)
        if (filename eq '') then begin
            return
        endif
        
        if (not file_test(filename, /read)) then begin
            void = dialog_message(['The file could not be read', $
                'Check the file permissions and try again.'], $
                /error, /center)
            return
        endif
        
        (*state).bval_filename = filename
        if ((*state).file_path eq '') then begin
            (*state).file_path = file_dirname(filename)
        endif
        selected_file = filename
        return
    endif
    
    if (keyword_set(bvecs)) then begin
    
        filename = dialog_pickfile(title='Select the text file containing the b-vectors', $
            path=file_path)
        if (filename eq '') then begin
            return
        endif
        
        if (not file_test(filename, /read)) then begin
            void = dialog_message(['The file could not be read', $
                'Check the file permissions and try again.'], $
                /error, /center)
            return
        endif
        
        (*state).bvec_filename = filename
        if ((*state).file_path eq '') then begin
            (*state).file_path = file_dirname(filename)
        endif
        selected_file = filename
        return
    endif
    
    if (keyword_set(other)) then begin
        
        filename = dialog_pickfile(title='Select the text file containing the relavent data', $
            path=file_path)
        if (filename eq '') then begin
            return
        endif
        
        if (not file_test(filename, /read)) then begin
            void = dialog_message(['The file could not be read', $
                'Check the file permissions and try again.'], $
                /error, /center)
            return
        endif
        
        if ((*state).file_path eq '') then begin
            (*state).file_path = file_dirname(filename)
    endif
        selected_file = filename
        return
    
    endif
end

; Subroutine name: mas_open_nifti_GUI_event
; Created by: BT, 2008-10
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
;    Event handler for nifti import GUI
;
; Editing Information:

pro mas_open_nifti_GUI_event, ev

    widget_control, ev.top, get_uvalue=state
    
    if (not ptr_valid(state)) then begin
        widget_control, ev.top, /destroy
        return
    endif
    
    uname = widget_info(ev.id, /uname)
    
    case uname of
    
        'choose_nifti': begin
            mas_open_nifti_GUI_browse, state, /nifti
            widget_control, (*state).nifti_txt, set_value=(*state).nifti_filename
        end
        
        'choose_bval': begin
            mas_open_nifti_GUI_browse, state, /bvals
            widget_control, (*state).bval_txt, set_value=(*state).bval_filename
        end
        
        'choose_bvec': begin
            mas_open_nifti_GUI_browse, state, /bvecs
            widget_control, (*state).bvec_txt, set_value=(*state).bvec_filename
        end
        
        'choose_tr': begin
            mas_open_nifti_GUI_browse, state, /other, selected_file=filename
            (*state).tr_filename = filename
            widget_control, (*state).tr_txt, set_value=filename
        end
        
        'choose_ti': begin
            mas_open_nifti_GUI_browse, state, /other, selected_file=filename
            (*state).ti_filename = filename
            widget_control, (*state).ti_txt, set_value=filename     
        end
        
        'choose_te': begin
            mas_open_nifti_GUI_browse, state, /other, selected_file=filename
            (*state).te_filename = filename
            widget_control, (*state).te_txt, set_value=filename
        end
        
        'start_import': begin
            read_status = 0
            nifti = mas_read_nifti(nifti_filename=(*state).nifti_filename, $
                bval_filename=(*state).bval_filename, $
                bvec_filename=(*state).bvec_filename, $
                tr_filename=(*state).tr_filename, $
                ti_filename=(*state).ti_filename, $
                te_filename=(*state).te_filename, $
                /hdr_only, $
                read_status=read_status)
                
            if (not read_status) then return
            
            mnf_nifti_hdr_to_mas, nifti
            ptr_free, state
            widget_control, ev.top, /destroy
        end
        
        'cancel': begin
            ptr_free, state
            widget_control, ev.top, /destroy
        end
        
        else:
        
    endcase
    
end

; Subroutine name: mas_open_nifti_GUI_cleanup
; Created by: BT, 2008-10
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
;    Cleans up state pointer when window is closed
;
; Editing Information:

pro mas_open_nifti_GUI_cleanup, tlb

    widget_control, tlb, get_uvalue=state
    if (ptr_valid(state)) then ptr_free, state
    return
    
end

; Subroutine name: mas_open_nifti_GUI
; Created by: BT, 2008-10
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
;   Creates the GUI window for nifti import
;
; Editing Information:

pro mas_open_nifti_GUI

    common common_widgets
    
    base = widget_base(title="NIFTI Import Options", group_leader=WID_BASE_MAIN, /column, /modal)
    
    b1 = widget_base(base, /row)
    lbl = widget_label(b1, value="NIFTI data file:            ")
    nifti_txt = widget_text(b1, value="", xsize=40)
    btn = widget_button(b1, value="Browse...", uname='choose_nifti')
    
    b2 = widget_base(base, /row)
    lbl = widget_label(b2, value="Diffusion b-values:         ")
    bval_txt = widget_text(b2, value="", xsize=40)
    btn = widget_button(b2, value="Browse...", uname='choose_bval')
    
    b3 = widget_base(base, /row)
    lbl = widget_label(b3, value="Diffusion b-vectors:        ")
    bvec_txt = widget_text(b3, value="", xsize=40)
    btn = widget_button(b3, value="Browse...", uname='choose_bvec')

    b4 = widget_base(base, /row)
    lbl = widget_label(b4, value="Repetition time values (ms):")
    tr_txt = widget_text(b4, value="", xsize=40)
    btn = widget_button(b4, value="Browse...", uname='choose_tr')

    b5 = widget_base(base, /row)
    lbl = widget_label(b5, value="Inversion time values (ms): ")
    ti_txt = widget_text(b5, value="", xsize=40)
    btn = widget_button(b5, value="Browse...", uname='choose_ti')

    b6 = widget_base(base, /row)
    lbl = widget_label(b6, value="Echo time values (ms):      ")
    te_txt = widget_text(b6, value="", xsize=40)
    btn = widget_button(b6, value="Browse...", uname='choose_te')
    
    b99 = widget_base(base, /row, /align_center)
    btn = widget_button(b99, value="Cancel", uname='cancel')
    btn = widget_button(b99, value="Start Import", uname='start_import')
    
    state = ptr_new({ nifti_filename: '', $
        bval_filename: '', $
        bvec_filename: '', $
        tr_filename: '', $
        ti_filename: '', $
        te_filename: '', $
        file_path: '', $
        tr_txt: tr_txt, $
        ti_txt: ti_txt, $
        te_txt: te_txt, $
        nifti_txt: nifti_txt, $
        bval_txt: bval_txt, $
        bvec_txt: bvec_txt })
        
    widget_control, base, set_uvalue=state
    
    widget_control, base, /realize
    
    xmanager, 'mas_open_nifti_GUI', base, cleanup='mas_open_nifti_GUI_cleanup'
    
end

;+
; :Description:
;    The "state1" loader for nifti files..
;
;
; :Author: wtriplett
;-Modifications
; - Zeropadding option added (Magdoom, 07/13)

pro mnf_load_state_1

    common scan_data, project
    
    ci = project.ci
    
    nifti_file = project.imndarray[ci].file_Path

    nif = mas_read_nifti(nifti_filename=nifti_file, read_status=rs, /progress)

    if (rs ne 0) then begin
        if (nif.datatype ne 16) then begin
            project.dataarray[ci].state1 = ptr_new(float(*nif.voxel_data))
            ptr_free, nif.voxel_data
        endif else begin
            project.dataarray[ci].state1 = nif.voxel_data
        endelse
    endif else begin
        ;; error message is generated in mas_read_nifti
        return
    endelse
    if project.procPramArray[CI].zpad_flag eq 1 then  begin
      dim = project.imndarray[CI].dimensions
      if dim eq 3 then *project.dataarray[CI].state1 = fft(*project.dataarray[CI].state1,/center,/overwrite) else begin
        *project.dataarray[CI].state1 = fft(*project.dataarray[CI].state1,dimension = 1,/overwrite)
        *project.dataarray[CI].state1 = fft(*project.dataarray[CI].state1,dimension = 2,/overwrite)
        *project.dataarray[CI].state1 = shift(*project.dataarray[CI].state1, project.imndArray[CI].fdim/2,project.imndArray[CI].pdim/2)
      endelse
      mas_zeropad, project.dataarray[CI].state1
      if dim eq 3 then *project.dataarray[CI].state1 = abs(fft(*project.dataarray[CI].state1,/center,/inverse,/overwrite)) else begin
        *project.dataarray[CI].state1 = abs(fft(*project.dataarray[CI].state1,dimension = 1,/center,/inverse,/overwrite))
        *project.dataarray[CI].state1 = abs(fft(*project.dataarray[CI].state1,dimension = 2,/center,/inverse,/overwrite))
      endelse
    endif
    
    IF project.procPramArray[CI].mc_enable NE 0 THEN BEGIN
        mas_motion_correct
    ENDIF
    
    project.dataArray[CI].state1 = temporary(mas_interpolate( project.dataArray[CI].state1 ))
    
    ;;mas_shift,            project.dataArray[CI].state1, ISDICOM=1
    ;;mas_signal_selection, project.dataArray[CI].state1
    mas_smoothing,        project.dataArray[CI].state1
    
    IF (project.procPramArray[ci].image_domain_filterflag eq 'Image Domain') THEN BEGIN
        image_domain_filter
    ENDIF

    project.procpramarray[ci].state_1 = 1
    
end

; Subroutine name: mas_open_nifti
; Created by: BT, 2008-10
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:

pro mas_open_nifti

    mas_open_nifti_GUI
    return
    
end

; Subroutine name: mas_nifti_header__define
; Created by: BT, 2008-10
; Calling Information:
;
; not generally called explicitly
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Serves as a definition for the MAS_NIFTI_HEADER struct.
; 
; See: http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields
; 
; Editing Information:

pro mas_nifti_header__define

    struct = { MAS_NIFTI_HEADER, $
                sizeof_hdr: long(0), $
                data_type:  bytarr(10), $
                db_name: bytarr(18), $
                extents:  long(0), $
                session_error: 0, $
                regular: 0B, $
                dim_info: 0B, $
                dim: intarr(8), $
                intent_p1: 0.0, $
                intent_p2: 0.0, $
                intent_p3: 0.0, $
                intent_code: 0, $
                datatype: 0, $
                bitpix: 0, $
                slice_start: 0, $
                pixdim: fltarr(8), $
                vox_offset: 0.0, $
                scl_slope: 0.0, $
                scl_inter: 0.0, $
                slice_end: 0, $
                slice_code: 0B, $
                xyzt_units: 0B, $
                cal_max: 0.0, $
                cal_min: 0.0, $
                slice_duration: 0.0, $
                toffset: 0.0, $
                glmax: 0L, $
                glmin: 0L, $
                descrip: bytarr(80), $
                aux_file: bytarr(24), $
                qform_code: 0, $
                sform_code: 0, $
                quatern_b: 0.0, $
                quatern_c: 0.0, $
                quatern_d: 0.0, $
                qoffset_x: 0.0, $
                qoffset_y: 0.0, $
                qoffset_z: 0.0, $
                srow_x: fltarr(4), $
                srow_y: fltarr(4), $
                srow_z: fltarr(4), $
                intent_name: bytarr(16), $
                magic: bytarr(4) $
             }
        
end

; Subroutine name: mas_nifti_file__define
; Created by: BT, 2008-10
; Calling Information:
;
; not generally called explicitly
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Serves as a definition for the MAS_NIFTI_FILE struct.
; 
; 
; Editing Information:

pro mas_nifti_file__define

    struct = { MAS_NIFTI_FILE, inherits MAS_NIFTI_HEADER, $
        nifti_filename: '', $
        bval_filename: '', $
        bvec_filename: '', $
        tr_filename: '', $
        te_filename: '', $
        ti_filename: '', $
        voxel_data: ptr_new() }
        
end
