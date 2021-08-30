;; $Id$
;;
;; Subroutine name: mas_odf_gui_event
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; event handler for the ODF gui window
;;
;; Editing Information:

pro mas_odf_gui_event, event

    common scan_data

    ci = project.ci
    
    odf_state = project.procpramarray[ci].odf_state
    difftools_state = project.procpramarray[ci].difftools_state

    if (not ptr_valid(odf_state)) then return

    uname = widget_info(event.id, /uname)
    widget_control, event.id, get_value=val
    
    ;; in case of being called from probtrack gui
    if (not ptr_valid(difftools_state)) then begin
        widget_control, event.top, get_uvalue=widget_state
    endif else begin
        ;;widget_control, (*difftools_state).odf_base, get_uvalue=widget_state
    endelse

    case uname of

        'sl_lmax'        : (*odf_state).odf_lmax        = long(val)
        'sl_b_thr'       : (*odf_state).odf_bval        = long(val)
        'sl_thr'         : (*odf_state).odf_signal_thr  = long(val)
        'sl_fa_thr'      : (*odf_state).odf_anis_thr    = float(val)
        'sl_anis_scl'    : (*odf_state).odf_anis_scl    = float(val)
        'sl_sharpen'     : (*odf_state).odf_sharpen     = long(val)
        'btn_anis'       : begin
            (*odf_state).odf_use_derived_anis_image = event.select
        end

        'btn_create': begin
            
           mas_odf_init_obj, odf_obj=reco_obj, $
                             lmax=(*odf_state).odf_lmax, $
                             b_thr=(*odf_state).odf_bval, $
                             lambda=0.01, $
                             /bulk

           if ((*difftools_state).disp_3d) then begin
              ;; 3D glyph view
              case project.procpramarray[ci].slice_axis of
                 0: slice = project.procpramarray[ci].sdim_start
                 1: slice = project.procpramarray[ci].pdim_start
                 2: slice = project.procpramarray[ci].fdim_start
              endcase
              
                mdt_make_glyph3d, reco_obj, [ slice ], threshold=(*odf_state).odf_signal_thr
                
           endif else begin
              
;               if ((*odf_state).odf_use_derived_anis_image ne 0) then begin
;                  ga = mas_odf_compute_ga(reco_obj, scaling=(*odf_state).odf_anis_scl)
;               endif
              
              mdt_make_glyphimage, reco_obj, $
                                   $;image_in=ga, $
                                   anis_threshold=(*odf_state).odf_anis_thr, $
                                   threshold=(*odf_state).odf_signal_thr, $
                                   /display
              
              obj_destroy, reco_obj
              
           endelse
        end
        
        'btn_batch': begin

            slices = mas_get_multi_slice_selection(event.top, count=count)
            if (count eq 0) then return

            mas_odf_init_obj, odf_obj=reco_obj, $
                              lmax=(*odf_state).odf_lmax, $
                              lambda=0.01, $
                              /bulk

            if ((*difftools_state).disp_3d) then begin

               mdt_make_glyph3d, reco_obj, [ slices ], threshold=(*odf_state).odf_signal_thr                

            endif else begin

               mas_odf_init_obj, odf_obj=reco_obj, $
                                 lmax=(*odf_state).odf_lmax, $
                                 lambda=0.01, $
                                 /bulk
               
               for sl = 0, n_elements(slices)-1 do begin
                  slice=slices[sl]
                  mdt_make_glyphimage, reco_obj, $
                                       anis_threshold=(*odf_state).odf_anis_thr, $
                                       vis_slice=slice, $
                                       threshold=(*odf_state).odf_signal_thr, $
                                       /display

               endfor
               
               obj_destroy, reco_obj

            endelse

        end
            
        else: begin
            print, 'mas_odf_gui_event: recv''d event from unknown widget: '+uname
        end

    endcase

end

;; Subroutine name: mas_odf_gui_cleanup
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Cleans up pointers associated with the ODF gui window
;;
;; Editing Information:

pro mas_odf_gui_cleanup, tlb

    common scan_data

    widget_control, tlb, get_uvalue=widget_state 
    if (ptr_valid(widget_state)) then begin
        ptr_free, widget_state
    endif

    if (project.odf_tlb ne tlb) then begin
        ;;print, "ODF tlb not consistent!"
    endif

    project.odf_tlb = -1

end

;; Subroutine name: mas_odf_gui_make_default_state
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; creates a pointer to a state structure containing the default ODF parametrs
;;
;; Editing Information:

function mas_odf_gui_make_default_state

    return, ptr_new({ ODF_STATE, $
                      odf_signal_thr:5,  $
                      odf_anis_thr:0.001, $
                      odf_anis_scl:0.175, $
                      odf_lmax:4,       $
                      odf_sharpen:1,    $
                      odf_use_derived_anis_image: 0, $
                      odf_bval:400      $
                    }, /no_copy)

end

;; Subroutine name: mas_odf_redraw
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Redraws the GUI, triggered by mas control
;;
;; Editing Information:

pro mas_odf_redraw, base_id

    common scan_data, project

    ci = project.ci

    odf_state = project.procpramarray[ci].odf_state
    widget_control, base_id, get_uvalue=widget_state

    if (not ptr_valid(odf_state)) then begin
        
        if (project.ni eq 0) then begin

            sen = 0

            widget_control, (*widget_state).sl_lmax   , sensitive=sen
            widget_control, (*widget_state).sl_b_thr   , sensitive=sen
            widget_control, (*widget_state).sl_odf_thr, sensitive=sen
            widget_control, (*widget_state).sl_anis_scl, sensitive=sen
            widget_control, (*widget_state).sl_sharpen, sensitive=sen
            widget_control, (*widget_state).sl_fa_thr , sensitive=sen
            widget_control, (*widget_state).btn_anis  , sensitive=sen
            
            return

        endif else begin

            odf_state = mas_odf_gui_make_default_state()
            project.procpramarray[ci].odf_state = odf_state
            return

        endelse

    endif

    sen = 1

    have_adt = project.procpramarray[ci].adt_proccess_flag

    widget_control, (*widget_state).sl_lmax   ,  sensitive=sen, set_value=(*odf_state).odf_lmax
    widget_control, (*widget_state).sl_b_thr  ,  sensitive=sen, set_value=(*odf_state).odf_bval
    widget_control, (*widget_state).sl_sharpen,  sensitive=sen, set_value=(*odf_state).odf_sharpen
    widget_control, (*widget_state).sl_odf_thr,  sensitive=sen, set_value=(*odf_state).odf_signal_thr
    widget_control, (*widget_state).sl_fa_thr ,  sensitive=sen*have_adt, set_value=(*odf_state).odf_anis_thr
;    widget_control, (*widget_state).sl_anis_scl, sensitive=sen*(*odf_state).odf_use_derived_anis_image, $
;                                                 set_value=(*odf_state).odf_anis_scl
;    widget_control, (*widget_state).btn_anis   , sensitive=sen, set_button=(*odf_state).odf_use_derived_anis_image

end

;;
;; Subroutine name: mas_odf_gui_make
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Creates the GUI window and all of its widgets, prepares the
;; scan's default parameters
;;
;; Editing Information:

pro mas_odf_gui_make, base, state

    common scan_data
    common common_widgets

    ci = project.ci

    if (ptr_valid(project.procpramarray[ci].odf_state)) then begin
        odf_state = project.procpramarray[ci].odf_state
    endif else begin
        odf_state = mas_odf_gui_make_default_state()
        project.procpramarray[ci].odf_state = odf_state
    endelse
 
    param_base = widget_base(base, /row, xpad=0)
    pbase_c1 = widget_base(param_base, /column, /frame)
    sl_lmax = widget_slider(pbase_c1, $
                            uname='sl_lmax', $
                            title='Harmonic Order (L)', $
                            minimum=0, $
                            maximum=10, $
                            value=(*odf_state).odf_lmax, $
                            event_pro='mas_odf_gui_event', $
                            scroll=2)

    sl_sharpen = widget_slider(pbase_c1, $
                            uname='sl_sharpen', $
                            title='Sharpening Factor', $
                            minimum=1, $
                            maximum=10, $
                            value=(*odf_state).odf_sharpen, $
                            event_pro='mas_odf_gui_event', $
                            scroll=1)

    sl_b_thr = widget_slider(pbase_c1, $
                             title='Low B-value Threshold', $
                             minimum=100,  $
                             maximum=1000, $
                             value=(*odf_state).odf_bval,    $
                             scroll=100,   $
                             event_pro='mas_odf_gui_event', $
                             uname='sl_b_thr')

    sl_odf_thr  = widget_slider(pbase_c1, $
                             uname='sl_thr', $
                             title='% Signal Threshold', $
                             minimum=0, $
                             maximum=20, $
                             value=(*odf_state).odf_signal_thr, $
                             event_pro='mas_odf_gui_event', $
                             scroll=1)

    pbase_c2 = widget_base(param_base, /column, /frame)

    sl_fa_thr = cw_fslider(pbase_c2, $
                           title='FA Threshold', $
                           value=(*odf_state).odf_anis_thr,  $
                           minimum=0.001, $
                           maximum=1.0,$
                           scroll=0.5,  $
                           /edit,     $
                           uname='sl_fa_thr')
    widget_control, sl_fa_thr, sensitive=project.procpramarray[ci].adt_proccess_flag

;    da_base = widget_base(pbase_c2, /column, ypad=0, /frame)
;
;    non_ex_base = widget_base(da_base, /nonexclusive)
;    btn_anis = widget_button(non_ex_base, $
;                             value="Use DA Background", $
;                             uname='btn_anis', $
;                             event_pro='mas_odf_gui_event')
;    widget_control, btn_anis, set_button=(*odf_state).odf_use_derived_anis_image
;
;    sl_anis_scl = cw_fslider(da_base, title='DA Scaling (Max)', $
;                       value=(*odf_state).odf_anis_scl,  $
;                       minimum=0.001, $
;                       maximum=1.0,$
;                       scroll=0.1,  $
;                       /edit,     $
;                       uname='sl_anis_scl')
;    widget_control, sl_anis_scl, sensitive=(*odf_state).odf_use_derived_anis_image

    widget_state = ptr_new({ ODF_WIDGET_STATE, $
                             sl_lmax: sl_lmax, $
                             sl_b_thr: sl_b_thr, $
                             sl_anis_scl: 0, $;;sl_anis_scl, $
                             sl_sharpen: sl_sharpen, $
                             btn_anis: 0, $;;btn_anis, $
                             sl_fa_thr:sl_fa_thr, $
                             sl_odf_thr: sl_odf_thr}, /no_copy)
    widget_control, base, set_uvalue=widget_state

end

;; Subroutine name: mas_dot_compute_ga
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;     DOT_OBJ: DOT reconstructor object
;;   
;;     VIS_SLICE: slice to use. If not set, current MAS selection
;;                according to slice axis is used
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;; 
;;     Computes Generalized Anisotropy from DOT weights
;;
;; Editing Information:

function mas_odf_compute_ga, odf_obj, scaling=scaling, vis_slice=vis_slice
    
    common scan_data
    
    ci = project.ci
    n_slices = n_elements(vis_slice)

    if (not keyword_set(scaling)) then scaling = 1.0
    
    if (n_slices eq 0) then begin
       case project.procpramarray[ci].slice_axis of
          0: slice = project.procpramarray[ci].sdim_start
          1: slice = project.procpramarray[ci].pdim_start
          2: slice = project.procpramarray[ci].fdim_start
       endcase
    endif else begin
       slice = vis_slice[0]
    endelse
    data_sz = size(*project.dataarray[ci].state1, /dimensions)

    case project.procpramarray[ci].slice_axis of
       0: image = reform((*project.dataarray[ci].state1)[*,*,slice,*])
       1: image = reform((*project.dataarray[ci].state1)[slice,*,*,*])
       2: image = reform((*project.dataarray[ci].state1)[*,slice,*,*])
    endcase

    xdim = (size(image, /dimensions))[0]
    ydim = (size(image, /dimensions))[1]

    ga  = FLTARR(xdim,ydim)

    for x = 0, xdim-1 do begin
       for y = 0, ydim-1 do begin

          tmp = odf_obj->getDisplacementProbability(data=image[x,y,*])
          if total(finite(tmp, /nan)) ne 0 then continue
          tot = sqrt(total(tmp^2))
          if (tot eq 0) then continue
          mean_tmp = mean(tmp)
          ga[x,y] = sqrt(n_elements(tmp)/(n_elements(tmp-1))) * (sqrt( total( (tmp - mean_tmp)^2 ))/tot)
          
       endfor
    endfor
    
    return, ga/scaling

end

;; Subroutine name: mas_odf_init_obj
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;    LMAX:   The highest order of the harmonic basis.
;;
;;    LAMBDA: Tikhonov constant, not used. Here for backward compat.
;;
;;    BULK:   Set this keyword to get a "BULK_RECONSTRUCTOR" version of
;;            the object
;;
;;   ODF_OBJ: set this keyword to a named variable that will contain
;;            the prepared object.
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Creates an instance of the ODF reconstructor object using imaging
;; parameters from MAS's data structures.
;;
;; Editing Information:

pro mas_odf_init_obj, lmax=lmax, $
                      lambda=lambda, $
                      b_thr=b_thr, $
                      odf_obj=odf_obj, $
                      bulk=bulk

   common scan_data
   
   ci = project.ci

   if (n_elements(lmax) eq 0) then begin
      lmax = 4
   endif

   ;; not used.
   ;;if (n_elements(lambda) eq 0) then begin
   ;;   lambda = 0.1
   ;;endif

   b_values = reform(*project.imndarray[ci].bval_array)

   the = reform(*project.imndarray[ci].angle_theta)
   phi = reform(*project.imndarray[ci].angle_phi)
   
   gradients = [ [the], [phi] ]
   sz_grad = size(gradients, /dimension)
   if (sz_grad[0] eq 2) then begin
      gradients = transpose(gradients)
   endif

   if (keyword_set(bulk)) then begin

      odf_obj = obj_new('mas_odf_bulk_reconstructor', $
                        l_max=lmax, $
                        bval_thr=b_thr, $
                        b_values=b_values, $
                        gradients=gradients, /spherical, /degrees, $
                        /debug, /use_tikhonov)

   endif else begin
      
      odf_obj = obj_new('mas_odf', $
                        l_max=lmax, $
                        bval_thr=b_thr, $
                        b_values=b_values, $
                        gradients=gradients, /spherical, /degrees, $
                        /use_tikhonov)
   endelse
   
end

pro mas_odf_tools

end
