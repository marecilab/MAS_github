;; $Id$
;;

function mas_mow_test_nnls_so, location
    
    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        print, !ERROR_STATE.MSG
        void = dialog_message("C based nnls solver not working. Falling back to IDL solver.", /center, /error)
        return, 0
    endif
    
    test_vector = double([0,1,0])
    A = double([[1,0,0], [0,0,1], [0,1,0]])
    m = (size(A, /dimensions))[0]; = mda in nnls.c
    n = (size(A, /dimensions))[1]
    w = dblarr(n+1)
    x = w ;; same size
    mode = 0L
    rnorm = 0D
    index = lonarr(n+2)
    zz = dblarr(m+1)
    ok = call_external(location, 'nnls_c',$
                       A, m, m, n, test_vector, x, rnorm, w, zz, index, mode,$
                       /cdecl, /auto_glue)
    w = x[0:n-1]
    
    if (total(abs(w - double([0,0,1]))) ne 0) then begin
        print, "mas_mow_test_nnls_so: expected: 0,0,1, got "+strjoin(string(w, format='(I0)'), ',')+"."
        void = dialog_message("C based nnls solver not working. Falling back to IDL solver.", /center, /error)
        return, 0
    endif
    
    return, 1

end

;; Subroutine name: mas_mow_gui_make_default_state
;; Created by: BT, 2008-07
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; creates a pointer to a state structure containing the 
;; default MOW parametrs
;;
;; Editing Information:

function mas_mow_gui_make_default_state

    common scan_data

    return, ptr_new({ MOW_STATE, $
                      mow_signal_thr:5,  $
                      mow_diff_time: (*project.imndarray[project.ci].big_delta)[0], $
                      mow_diff_sep:  (*project.imndarray[project.ci].small_delta)[0], $
                      mow_anis_thr:0.001, $
                      mow_reco_order:2, $
                      mow_deco_order:2, $
                      mow_show_tess:0,  $
                      mow_bval:400,     $
                      mow_comp_method: 1, $
                      mow_radius:10.0,  $
                      mow_use_nnls_so: 0B, $
                      mow_nnls_so_location: '' $
                    }, /no_copy)
end

;; Subroutine name: mas_mow_gui_event
;; Created by: BT, 2008-07
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Event handler for MOW gui window or diff tools tab
;;
;; Editing Information:

pro mas_mow_gui_event, event

    common scan_data

    ci = project.ci
    
    mow_state = project.procpramarray[ci].mow_state
    difftools_state = project.procpramarray[ci].difftools_state
    
    if (not ptr_valid(mow_state)) then return

    uname = widget_info(event.id, /uname)
    widget_control, event.id, get_value=val

    case uname of

        'sl_glyph_size'  : (*mow_state).mow_glyph_size  = long(val)
        'sl_radius'      : (*mow_state).mow_radius      = double(val)
        'sl_b_thr'       : (*mow_state).mow_bval        = long(val)
        'sl_thr'         : (*mow_state).mow_signal_thr  = long(val)
        'sl_fa_thr'      : (*mow_state).mow_anis_thr    = float(val)
        'btn_show_tess'  : (*mow_state).mow_show_tess   = event.select
        'dl_mow_method'  : (*mow_state).mow_comp_method = event.index
        'txt_diff_sep'   : begin
           temp = float(val); & print, temp
           (*mow_state).mow_diff_sep = temp
        end
        'txt_diff_time'   : begin
           temp = float(val); & print, temp
           (*mow_state).mow_diff_time = temp
        end
        'tf_desc_order_deco': begin
            val = event.value
            if (val gt 3 or val lt 0) then begin
                val = (val < 3) > 0
                widget_control, event.id, set_value=val
            endif
            (*mow_state).mow_deco_order = val
        end
        'tf_desc_order_reco': begin
            val = event.value
            if (val gt 3 or val lt 0) then begin
                val = (val < 3) > 0
                widget_control, event.id, set_value=val
            endif
            (*mow_state).mow_reco_order = val
        end
        'btn_use_nnls_so': begin
            if (event.select eq 1) then begin
                if ((*mow_state).mow_nnls_so_location eq '' or not file_test((*mow_state).mow_nnls_so_location)) then begin

                    so_loc = dialog_pickfile(path=project.current_path, filter='*.so;*.dylib;*.dll', $
                                             title="Please locate the shared library nnls.so", /file)

                    if (so_loc eq '') then begin
                        (*mow_state).mow_use_nnls_so = 0B
                        widget_control, event.id, set_button=0
                        return
                    endif else if (mas_mow_test_nnls_so(so_loc) eq 0) then begin
                        (*mow_state).mow_use_nnls_so = 0B
                        widget_control, event.id, set_button=0
                        return
                    endif
                    
                    (*mow_state).mow_nnls_so_location = so_loc
                    (*mow_state).mow_use_nnls_so = 1B
                    return
                    
                endif 
                
                (*mow_state).mow_use_nnls_so = 1B
                return
                
            endif 
            
            (*mow_state).mow_use_nnls_so = 0B
            
        end
        
        'btn_create': begin
           
           mas_mow_init_obj, mow_obj=reco_obj, $
                             diff_time=(*mow_state).mow_diff_time, $
                             diff_sep=(*mow_state).mow_diff_sep, $
                             bval_thr=(*mow_state).mow_bval, $
                             radius=(*mow_state).mow_radius, $
                             deco_level=(*mow_state).mow_deco_order, $
                             nnls_so_location=(*mow_state).mow_use_nnls_so eq 1 ? (*mow_state).mow_nnls_so_location : '', $
                             comp_meth=(*mow_state).mow_comp_method, $
                             /bulk
           if ((*difftools_state).disp_3d) then begin
              ;; 3D glyph view
              case project.procpramarray[ci].slice_axis of
                 0: slice = project.procpramarray[ci].sdim_start
                 1: slice = project.procpramarray[ci].pdim_start
                 2: slice = project.procpramarray[ci].fdim_start
              endcase
              
                mdt_make_glyph3d, reco_obj, [ slice ], threshold=(*mow_state).MOW_SIGNAL_THR
                
           endif else begin
           
               mdt_make_glyphimage, reco_obj, $
                                anis_threshold=(*mow_state).mow_anis_thr, $
                                threshold=(*mow_state).mow_signal_thr, $
                                /display
           endelse
           
               obj_destroy, reco_obj

        end

        'btn_batch': begin
            slices = mas_get_multi_slice_selection(event.top, count=count)
            if (count eq 0) then return

            mas_mow_init_obj, mow_obj=reco_obj, $
                              diff_time=(*mow_state).mow_diff_time, $
                              diff_sep=(*mow_state).mow_diff_sep, $
                              bval_thr=(*mow_state).mow_bval, $
                              radius=(*mow_state).mow_radius, $
                              deco_level=(*mow_state).mow_deco_order, $
                              comp_meth=(*mow_state).mow_comp_method, $
                              /bulk
            if ((*difftools_state).disp_3d) then begin

               mdt_make_glyph3d, reco_obj, [ slices ], threshold=(*mow_state).MOW_SIGNAL_THR

            endif else begin
            
               for sl = 0, n_elements(slices)-1 do begin
                  slice=slices[sl]
                  mdt_make_glyphimage, reco_obj, $
                                    vis_slice=slice, $
                                    anis_threshold=(*mow_state).mow_anis_thr, $
                                    threshold=(*mow_state).mow_signal_thr, $
                                    /display
               endfor
               
            endelse
            
            obj_destroy, reco_obj

        end
            
        else: begin
            print, 'mas_mow_gui_event: recv''d event from unknown widget: '+uname
            help, event, /structure
        end

    endcase


end

;; Subroutine name: mas_mow_redraw
;; Created by: BT, 2008-06
;; Calling Information:
;
;; Bugs or Important Comments to Developers:

;; Purpose of subroutine:

;; Redraws the GUI, triggered by mas control

;; Editing Information:

pro mas_mow_redraw, mow_base

    common scan_data, project

    ci = project.ci

    mow_state = project.procpramarray[ci].mow_state
    
    if (not ptr_valid(mow_state)) then begin
        
        if (project.ni eq 0) then begin

            sen = 0
            widget_control, mow_base, get_uvalue=widget_state
            
            widget_control, (*widget_state).sl_radius       , sensitive=sen
            widget_control, (*widget_state).sl_b_thr        , sensitive=sen
            widget_control, (*widget_state).sl_mow_thr      , sensitive=sen
            widget_control, (*widget_state).sl_fa_thr       , sensitive=sen
            widget_control, (*widget_state).tf_desc_order_deco, sensitive=sen
            widget_control, (*widget_state).tf_desc_order_reco, sensitive=sen
            ;widget_control, (*widget_state).btn_show_tess   , sensitive=sen
            
            return

        endif else begin

            mow_state = mas_mow_gui_make_default_state()
            project.procpramarray[ci].mow_state = mow_state

        endelse

    endif

    sen = 1
    widget_control, mow_base, get_uvalue=widget_state
    have_adt = project.procpramarray[ci].adt_proccess_flag

    widget_control, (*widget_state).sl_radius         , sensitive=sen, set_value=(*mow_state).mow_radius
    widget_control, (*widget_state).sl_b_thr          , sensitive=sen, set_value=(*mow_state).mow_bval
    widget_control, (*widget_state).sl_mow_thr        , sensitive=sen, set_value=(*mow_state).mow_signal_thr
    widget_control, (*widget_state).sl_fa_thr        , sensitive=sen*have_adt, set_value=(*mow_state).mow_anis_thr
    widget_control, (*widget_state).tf_desc_order_deco, sensitive=sen, set_value=(*mow_state).mow_deco_order
    widget_control, (*widget_state).tf_desc_order_reco, sensitive=sen, set_value=(*mow_state).mow_reco_order
    widget_control, (*widget_state).dl_mow_method   , sensitive=sen, set_combobox_select=(*mow_state).mow_comp_method
    ;widget_control, (*widget_state).btn_show_tess   , sensitive=sen, set_button=(*mow_state).mow_show_tess

end

;; Subroutine name: mas_mow_gui_cleanup
;; Created by: BT, 2008-07
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; cleans up pointers. called by xmanager when window
;; is destroyed
;;
;; Editing Information:

pro mas_mow_gui_cleanup, tlb

    common scan_data

    widget_control, tlb, get_uvalue=widget_state 
    if (ptr_valid(widget_state)) then begin
        ptr_free, widget_state
    endif

    if (project.mow_tlb ne tlb) then begin
        ;print, "MOW tlb not consistent!"
    endif

    project.mow_tlb = -1

end

;; Subroutine name: mas_mow_gui_make
;; Created by: BT, 2008-07
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Note that this may be called to create model
;; windows, like it is within probtrack GUI or 
;; not as in the difftools tabbed window
;;
;; Purpose of subroutine:
;;
;; creates the MOW gui window
;;
;; Editing Information:

pro mas_mow_gui_make, base, state

    common scan_data
    common common_widgets

    ci = project.ci

    if (ptr_valid(project.procpramarray[ci].mow_state)) then begin
        mow_state = project.procpramarray[ci].mow_state
    endif else begin
        mow_state = mas_mow_gui_make_default_state()
        project.procpramarray[ci].mow_state = mow_state
    endelse
 
    param_base = widget_base(base, column=2, xpad=0, /grid);;widget_base(mow_base_c1, /column, /frame)

    param_base_c1 = widget_base(param_base, /column, /frame)

    junk = widget_label(param_base_c1, value='Computation Method:', /align_left)

    dl_mow_method = widget_combobox(param_base_c1, $
                                    /dynamic_resize, $
                                    event_pro='mas_mow_gui_event', $
                                    uname='dl_mow_method', $
                                    value=['Damped LS', $
                                           'Non-negative LS', $
                                           'SVD Pseudoinverse'])

    widget_control, dl_mow_method, set_combobox_select=(*mow_state).mow_comp_method

    junk = widget_base(param_base_c1, /nonexclusive)
    btn_use_nnls_so = widget_button(junk, value="Use C-based NNLS solver", uname='btn_use_nnls_so')
    widget_control, btn_use_nnls_so, set_button=(*mow_state).mow_use_nnls_so
    
    junk = widget_label(param_base_c1, value='Discretization order:', /align_left)
    
    tf_desc_order_deco = cw_field(param_base_c1, title='Deconvolution:', $
                                  uname='tf_desc_order_deco', $
                                  ;event_pro='mas_mow_gui_event', $
                                  /all_events, $
                                  /integer, $
                                  xsize=2, $
                                  value=(*mow_state).mow_deco_order)

    tf_desc_order_reco = cw_field(param_base_c1, title='Reconstruction:',$
                                  uname='tf_desc_order_reco', $
                                  ;event_pro='mas_mow_gui_event', $
                                  /all_events, $
                                  /integer,$
                                  xsize=2, $
                                  value=(*mow_state).mow_reco_order)

;;    non_ex_base = widget_base(param_base_c1, /align_left, /nonexclusive)
;;    btn_show_tess = widget_button(non_ex_base, value='Display vertices', uname='btn_show_tess', event_pro='mas_mow_gui_event')
;;    widget_control, btn_show_tess, set_button=(*mow_state).mow_show_tess

    sl_radius = cw_fslider(param_base_c1, $
                           uname='sl_radius', $
                           title='Radius (um)', $
                           minimum=.00001, $
                           maximum=25.0, $
                           /edit, $
                           value=(*mow_state).mow_radius, $
                           scroll=1)

    param_base_c2 = widget_base(param_base, /column, /frame)

    sl_b_thr = widget_slider(param_base_c2, title='Low B-value Threshold', $
                             minimum=100,  $
                             maximum=1000, $
                             event_pro='mas_mow_gui_event', $
                             value=(*mow_state).mow_bval,    $
                             scroll=100,   $
                             uname='sl_b_thr')


    sl_mow_thr  = widget_slider(param_base_c2, $
                                uname='sl_thr', $
                                event_pro='mas_mow_gui_event', $
                             title='% Signal Threshold', $
                             minimum=0, $
                             maximum=20, $
                             value=(*mow_state).mow_signal_thr, $
                             scroll=1)

    sl_fa_thr = cw_fslider(param_base_c2, $
                         uname='sl_fa_thr', $
                         title='FA Threshold', $
                         minimum=.0001, $
                         maximum=1.0, $
                         /edit, $
                         value=(*mow_state).mow_anis_thr, $
                         scroll=.1)
    widget_control, sl_fa_thr, sensitive=project.procpramarray[ci].adt_proccess_flag
    
    dt_base = widget_base(param_base_c2, /row)
    lbl_dt = widget_label(dt_base, value="Diffusion Time:")
    txt_diff_time = widget_text(dt_base, xsize=7, $
                                uname='txt_diff_time', $
                                event_pro='mas_mow_gui_event', $
                                /editable, /all_events, $
                                value=strcompress((*mow_state).mow_diff_time, /remove_all))

    ds_base = widget_base(param_base_c2, /row)
    lbl_ds = widget_label(ds_base, value="Diff Pulse Length:")
    txt_diff_sep = widget_text(ds_base, xsize=7, $
                               uname='txt_diff_sep', $
                               event_pro='mas_mow_gui_event', $
                               /editable, /all_events, $
                               value=strcompress((*mow_state).mow_diff_sep, /remove_all))

    widget_state = ptr_new({ MOW_WIDGET_STATE, $
                             dl_mow_method: dl_mow_method, $
                             tf_desc_order_deco: tf_desc_order_deco, $
                             tf_desc_order_reco: tf_desc_order_reco, $
                             txt_diff_sep: txt_diff_sep, $
                             txt_diff_time: txt_diff_time, $
$;;                             btn_show_tess: btn_show_tess, $
                             sl_b_thr: sl_b_thr, $
                             sl_fa_thr:sl_fa_thr, $
                             sl_radius: sl_radius, $
                             sl_mow_thr: sl_mow_thr } $
                             , /no_copy) 

    widget_control, base, set_uvalue=widget_state

end

;; Subroutine name: mas_mow_init_obj
;; Created by: BT, 2009-01
;; Calling Information:
;
;; Bugs or Important Comments to Developers:

;; Purpose of subroutine:

;; creates a MOW object from parameters in MAS

;; Editing Information:

pro mas_mow_init_obj, diff_time=diff_time, $
                      diff_sep=diff_sep, $
                      comp_meth=comp_meth, $
                      bval_thr=bval_thr, $
                      deco_level=deco_level, $
                      radius=radius, $
                      mow_obj=mow_obj, $
                      bulk=bulk, $
                      nnls_so_location=nnls_so_location

    common scan_data

    ci = project.ci

    ;;; IMAGING PARAMETERS ;;;;
    ;; convert to s/m^2, appears to be what Evren uses
    bvals = reform(*project.imndarray[ci].bval_array) * 1e-3
    
    if (not keyword_set(comp_meth)) then begin
       comp_meth = 0
    endif

    if (not keyword_set(bval_thr)) then begin
        bval_thr = 400.
    endif
    bval_thr *= 1e-3
    
    if (not keyword_set(deco_level)) then begin
       deco_level = 0
    endif

    ;; diffusion time parameters
    if (keyword_set(diff_time)) then begin
       big_delta = float(diff_time)
    endif else begin
       big_delta = ((*project.imndArray[CI].big_delta)[0]) 
    endelse

    if (keyword_set(diff_sep)) then begin
       small_delta = float(diff_sep)
    endif else begin
       small_delta = ((*project.imndArray[CI].small_delta)[0])
    endelse

    if (keyword_set(nnls_so_location)) then begin
        nnls_so = nnls_so_location
    endif else begin
        nnls_so = ''
    endelse
    
    b_matrix = double(*project.imndarray[ci].b_matrix)
    b_matrix *= 1e-3

    td = big_delta - small_delta/3.
   
    if (keyword_set(bulk)) then begin

       mow_obj = obj_new('mas_mow_bulk_reconstructor', $
                         diff_time=td, $
                         computation_method=comp_meth, $
                         bval_thr=bval_thr, $
                         deconvolution_level=deco_level, $
                         radius=radius, $
                         nnls_so_location=nnls_so, $
                         b_matrix=b_matrix)

    endif else begin

       mow_obj = obj_new('mas_mow', $
                         diff_time=td, $
                         computation_method=comp_meth, $
                         bval_thr=bval_thr, $
                         deconvolution_level=deco_level, $
                         radius=radius, $
                         nnls_so_location=nnls_so, $
                         b_matrix=b_matrix)

    endelse
    
 end

pro mas_mow_tools

end
