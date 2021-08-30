;; $Id$
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_refresh_data
;; Created by:
;; Calling Information:
;;
;;    state: The state pointer containing preferences
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;    Usually called when the user selects a different image array to
;;    view. This procedure sets the proper pointer and performs and
;;    data operations necessary to display it in the window. 
;;
;; Editing Information:
;; Magdoom (01/14/15)
;; Included flow speed imagery (PC-MRI)

pro ortho_refresh_data, state

    common scan_data
    forward_function get_orient_rgb_axis
    
    adim = (*state).adim
    
    case (*state).current_imagery of        
        0: begin
           data_ptr = project.dataarray[(*state).proj_index].state1
           (*state).data_max = max((*project.dataarray[(*state).proj_index].state1)[*,*,*,adim])
        end
        1: begin
           data_ptr = project.dataarray[(*state).proj_index].frac_ani
           (*state).data_max = 1.
           adim = 0
        end
        2: begin
           data_ptr = project.dataarray[(*state).proj_index].avg_dif
           (*state).data_max = max(*project.dataarray[(*state).proj_index].avg_dif)
           adim = 0
         end
        3: begin ;; will be combined with evectors later in procedure
           data_ptr = project.dataarray[(*state).proj_index].frac_ani
           (*state).data_max = 1.
           adim = 0
        end
        4 : begin
            data_ptr = project.dataarray[(*state).proj_index].flow
            (*state).data_max = sqrt(3)
            adim = 3
          end   
        5: begin
           adim = 0
           data_ptr = (*state).data_ptr
           end
        else: begin
           data_ptr = project.dataarray[(*state).proj_index].state1
           (*state).data_max = max((*project.dataarray[(*state).proj_index].state1)[*,*,*,adim])
        end
    endcase

    if (not ptr_valid(data_ptr)) then begin
        (*state).data_available = 0
        return
    endif
    
    zoom =(*state).zoom
        
    if ((*state).current_imagery eq 3) then begin
       ;; hangle the eigenvector image a bit differently
       ;; grab the data only for the selected volume
       eigv = abs(reform((*project.dataarray[(*state).proj_index].eign_vec)[*,*,*,*,0], $
                         (*state).matrix[0], (*state).matrix[1], (*state).matrix[2],3))
       
       datar = reform( (*data_ptr)[*,*,*,adim], $
                       (*state).matrix[0], (*state).matrix[1], (*state).matrix[2] )
       datar_size = (size(datar))[1:3]
       
       axisid =  get_orient_rgb_axis()
       temp_eig = eigv
       temp_eig[*,*,*,axisid[0]] = eigv[*,*,*,0] * datar
       temp_eig[*,*,*,axisid[1]] = eigv[*,*,*,1] * datar
       temp_eig[*,*,*,axisid[2]] = eigv[*,*,*,2] * datar
       eigv = temporary(temp_eig)

       ;; bytescale + intensity adjust
       datar = bytscl(bytscl(temporary(eigv), (*state).min_data_thr, (*state).max_data_thr), $
                      top = (*state).center_intensity,$
                      max = (*state).max_intensity, $
                      min = (*state).min_intensity)
                      
       sz_pdata = size(*data_ptr, /dimensions)
       
       datar = reform(datar, (*state).matrix[0], (*state).matrix[1], (*state).matrix[2], 3)
            
    end else begin
       ;; grab the data only for the selected volume
       datar = reform( (*data_ptr)[*,*,*,adim], (*state).matrix[0], (*state).matrix[1], (*state).matrix[2] )
       datar_size = (size(datar))[1:3]
       
       ;; bytescale + intensity adjust
       datar = bytscl(bytscl(datar, (*state).min_data_thr, (*state).max_data_thr), $
                      top = (*state).center_intensity,$
                      max = (*state).max_intensity, $
                      min = (*state).min_intensity)

       sz_pdata = size(*data_ptr, /dimensions)

       datar = reform(datar, (*state).matrix[0], (*state).matrix[1], (*state).matrix[2])

    endelse

    ;; replace state data ptr
    ptr_free, (*state).datar
    (*state).datar = ptr_new(datar, /no_copy)
    (*state).data_available = 1
    (*state).data_ptr = data_ptr

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_update_imagery
;; Created by:
;; Calling Information:
;;
;;     state:       The state pointer containing the preferences
;;
;;     want_adim:   The index of the volume requested for 4D scans
;;
;;     ignore:     This keyword does nothing
;;
;; Bugs or Important Comments to Developers:
;;
;;     This procedure updates the image windows as the user moves
;;     around the data. It slices the data along the three axes to
;;     create a 2D image for each plane, then resizes the image
;;     according to the voxel dimensions and zoom factor.
;;
;; Purpose of subroutine:
;;
;; Editing Information:
;; Magdoom (01/15/15)
;; Added colormap to flow speed imagery

pro ortho_update_imagery, state, want_adim=want_adim, ignore=ignore, $
                          angulation=angulation, $
                          break_chain=break_chain, link_history=link_history

    common scan_data

    if ((*state).linked_scanid ne -1 and not keyword_set(break_chain)) then begin
    
        if (widget_info(project.procpramarray[(*state).linked_scanid].ortho_tlb, /valid_id)) then begin

            link_history = n_elements(link_history) eq 0 ? (*state).linked_scanid : [link_history, (*state).linked_scanid]    
            circulant = where(link_history eq (*state).proj_index, count)
            if (count ne 0) then begin
                break_chain = 1B
            endif else begin
                break_chain = 0B
            endelse
            
            widget_control, project.procpramarray[(*state).linked_scanid].ortho_tlb, get_uvalue=linked_state
            vec = [(*state).freq, (*state).phase, (*state).slice, 1.0] # (diag_matrix([(*state).voxel_dims_actual,0]) ## (*state).link_matrix ## diag_matrix([1.0/(*linked_state).voxel_dims_actual,0]) )

            (*linked_state).freq  = (vec[0] > 0) < ((*linked_state).matrix[0]-1)
            (*linked_state).phase = (vec[1] > 0) < ((*linked_state).matrix[1]-1)
            (*linked_state).slice = (vec[2] > 0) < ((*linked_state).matrix[2]-1)
            
            ortho_update_imagery, linked_state, break_chain=break_chain, link_history=link_history
        endif
    endif
    
    if ((*state).navigation_callback ne '') then begin
        void = where(strmatch(routine_info(), strupcase((*state).navigation_callback)), nmatch)
        if (nmatch ne 0) then begin
            call_procedure, (*state).navigation_callback, [(*state).freq, (*state).phase, (*state).slice], state
        endif
    endif

    if (not arg_present(want_adim)) then begin
       want_adim = max([(*state).adim,0])
    endif

    if (not keyword_set(ignore)) then ignore=''

    ;; make sure want_adim is in range
    if (want_adim ne (*state).adim and $
        want_adim ge 0 and $
        want_adim le project.imndarray[(*state).proj_index].adim) then begin

       (*state).adim = want_adim
       ortho_refresh_data, state
    endif

    if ((*state).current_imagery ne 0) then want_adim = 0
    matrix = (*state).matrix
    pdatar = (*state).datar
    
    ;; set the zoom scaling for non-isotropic voxel dimensions
    zoom = float((*state).zoom) * (*state).voxel_dims
    sz_datar = size(*pdatar, /dimensions)
    
    if ((*state).current_imagery eq 3) then begin
       ;; eigenvector imagery is handled differently
       RP = fltarr(sz_datar[0], sz_datar[1], 3)
       RS = fltarr(sz_datar[0], sz_datar[2], 3)
       PS = fltarr(sz_datar[1], sz_datar[2], 3)

;       for chan = 0, 2 do begin
;              RP[*,*,chan] = reform((*pdatar)[*, *, (*state).slice, chan])
;              RS[*,*,chan] = reform((*pdatar)[*, (*state).phase, *, chan])
;              PS[*,*,chan] = reform((*pdatar)[(*state).freq, *,  *, chan])
;       endfor
       
       RP[*,*,0] = reform((*pdatar)[*, *, (*state).slice, 0])
       RP[*,*,1] = reform((*pdatar)[*, *, (*state).slice, 1])
       RP[*,*,2] = reform((*pdatar)[*, *, (*state).slice, 2])
       
       RS[*,*,0] = reform((*pdatar)[*, (*state).phase, *, 0])
       RS[*,*,1] = reform((*pdatar)[*, (*state).phase, *, 1])
       RS[*,*,2] = reform((*pdatar)[*, (*state).phase, *, 2])

       PS[*,*,0] = reform((*pdatar)[(*state).freq, *,  *, 0])
       PS[*,*,1] = reform((*pdatar)[(*state).freq, *,  *, 1])
       PS[*,*,2] = reform((*pdatar)[(*state).freq, *,  *, 2])

       interleave = 2

    endif else begin

;       RP = extract_slice(*pdatar, matrix[0], matrix[1], $
;                          matrix[0]/2.0,matrix[1]/2.0,(*state).slice, $
;                          ang[0], ang[1], ang[2], $
;                          sample=0, out_val=0);, anisotropy=(*state).voxel_dims)
;                                
;       RS = extract_slice(*pdatar, matrix[0], matrix[2], $
;                          matrix[0]/2.0,(*state).phase,matrix[2]/2.0, $
;                          90, 0, 0, sample=1, out_val=0);, anisotropy=(*state).voxel_dims)
;                          
;       PS = extract_slice(*pdatar, matrix[1], matrix[2], $
;                          (*state).freq,matrix[1]/2.0,matrix[2]/2.0, $
;                          90, 90, 0, sample=1, out_val=0);, anisotropy=(*state).voxel_dims)

       ;; grab imagery to display                          
       RP = reform((*pdatar)[*, *, (*state).slice], sz_datar[0], sz_datar[1])
       RS = reform((*pdatar)[*, (*state).phase, *], sz_datar[0], sz_datar[2])
       PS = reform((*pdatar)[(*state).freq, *,  *], sz_datar[1], sz_datar[2])

       interleave = 1

    endelse

    if (*state).show_line_profile then begin
       ;; update the line profile plot
       (*state).line_prof_selected_window = ignore

       case ignore of
          'RP': begin
             mas_update_line_data_window, window_id=(*state).line_data_gui_win, fov=[(*state).fov_fps[[0,1]]], $
                                          horiz_data=(*((*state).data_ptr))[*, (*state).phase, (*state).slice, want_adim], $
                                          vert_data=(*((*state).data_ptr))[(*state).freq, *, (*state).slice, want_adim]
          end
          'PS': begin
             mas_update_line_data_window, window_id=(*state).line_data_gui_win, fov=[(*state).fov_fps[[1,2]]], $
                                          horiz_data=(*((*state).data_ptr))[(*state).freq, *, (*state).slice, want_adim], $
                                          vert_data=(*((*state).data_ptr))[(*state).freq, (*state).phase, *, want_adim]
          end
          'RS': begin
             mas_update_line_data_window, window_id=(*state).line_data_gui_win, fov=[(*state).fov_fps[[0,2]]], $
                                          horiz_data=(*((*state).data_ptr))[*, (*state).phase , (*state).slice, want_adim], $
                                          vert_data=(*((*state).data_ptr))[(*state).freq, (*state).phase, *, want_adim]
          end
          else:
       endcase

    endif

    ;; draw crosshairs
    if (*state).draw_crosshairs then begin

       freq_co  = ((*state).freq  + 0.5)
       phase_co = ((*state).phase + 0.5)
       slice_co = ((*state).slice + 0.5)
       
       ;; change the color based on whether the line profile is enabled
       ((*state).ov_RP->getByName('main_model/CH_model/crosshair_H'))->setProperty, $
          hide=0, data=[ [freq_co, -1], [freq_co, sz_datar[1]*zoom[1]-1] ], $
          color=(ignore eq 'RP' and (*state).show_line_profile) ? [255,0,0] : [255,255,255]
       ((*state).ov_RP->getByName('main_model/CH_model/crosshair_V'))->setProperty, $
          hide=0, data=[ [-1, phase_co], [sz_datar[0]*zoom[0]-1, phase_co] ], $
          color=(ignore eq 'RP' and (*state).show_line_profile) ? [0,0,255] : [255,255,255]

       ((*state).ov_RS->getByName('main_model/CH_model/crosshair_H'))->setProperty, $
          hide=0, data=[ [freq_co, -1], [freq_co, sz_datar[2]*zoom[2]-1] ], $
          color=(ignore eq 'RS' and (*state).show_line_profile) ? [255,0,0] : [255,255,255]
       ((*state).ov_RS->getByName('main_model/CH_model/crosshair_V'))->setProperty, $
          hide=0, data=[ [-1, slice_co], [sz_datar[0]*zoom[0]-1, slice_co] ], $
          color=(ignore eq 'RS' and (*state).show_line_profile) ? [0,0,255] : [255,255,255]

       ((*state).ov_PS->getByName('main_model/CH_model/crosshair_H'))->setProperty, $
          hide=0, data=[ [phase_co, -1], [phase_co, sz_datar[2]*zoom[2]-1] ], $
          color=(ignore eq 'PS' and (*state).show_line_profile) ? [255,0,0] : [255,255,255]
       ((*state).ov_PS->getByName('main_model/CH_model/crosshair_V'))->setProperty, $
          hide=0, data=[ [-1, slice_co], [sz_datar[1]*zoom[1]-1, slice_co] ], $
          color=(ignore eq 'PS' and (*state).show_line_profile) ? [0,0,255] : [255,255,255]

    endif else begin

       ((*state).ov_RP->getByName('main_model/CH_model/crosshair_H'))->setProperty, hide=1
       ((*state).ov_RP->getByName('main_model/CH_model/crosshair_V'))->setProperty, hide=1
       ((*state).ov_RS->getByName('main_model/CH_model/crosshair_H'))->setProperty, hide=1
       ((*state).ov_RS->getByName('main_model/CH_model/crosshair_V'))->setProperty, hide=1
       ((*state).ov_PS->getByName('main_model/CH_model/crosshair_H'))->setProperty, hide=1
       ((*state).ov_PS->getByName('main_model/CH_model/crosshair_V'))->setProperty, hide=1

    endelse

    ;; replace imagery and update view
    im_PS = (*state).ov_PS->getByName('main_model/IM_model/image')
    im_RS = (*state).ov_RS->getByName('main_model/IM_model/image')
    im_RP = (*state).ov_RP->getByName('main_model/IM_model/image')
   
    im_RP->setProperty, DATA=RP, DIMENSIONS=(size(RP, /dimensions))[0:1], INTERLEAVE=interleave
    im_RS->setProperty, DATA=RS, DIMENSIONS=(size(RS, /dimensions))[0:1], INTERLEAVE=interleave
    im_PS->setProperty, DATA=PS, DIMENSIONS=(size(PS, /dimensions))[0:1], INTERLEAVE=interleave
    
    if (ptr_valid((*state).overlay_data_ptr)) then begin
        ortho_update_overlay, state
    endif
    
    ortho_update_coords_and_int, state
    if ((*state).current_imagery eq 4) then begin
        oPalette = OBJ_NEW('IDLgrPalette')
        oPalette->LOADCT, 13
        im_RP->setProperty, PALETTE = oPalette
        im_RS->setProperty, PALETTE = oPalette
        im_PS->setProperty, PALETTE = oPalette
    endif
        
    (*state).owin_RP->draw, (*state).ov_RP
    ;(*state).owin_RP->draw, (*state).sv_RP
    (*state).owin_RS->draw, (*state).ov_RS
    ;(*state).owin_RS->draw, (*state).sv_RS
    (*state).owin_PS->draw, (*state).ov_PS
    ;(*state).owin_PS->draw, (*state).sv_PS
 end

pro ortho_update_overlay, state

    if (not ptr_valid((*state).overlay_data_ptr)) then return
    
    ;; replace imagery and update view
    im_overlay_PS = (*state).ov_PS->getByName('main_model/IM_OVL_model/image')
    im_overlay_RS = (*state).ov_RS->getByName('main_model/IM_OVL_model/image')
    im_overlay_RP = (*state).ov_RP->getByName('main_model/IM_OVL_model/image')

    ol_dim = size(*((*state).overlay_data_ptr))
    dt_dim = size(*(*state).data_ptr)
    adim = 0
    if (ol_dim[0] eq 4 and ol_dim[0] eq 4) then begin
        if (ol_dim[4] eq dt_dim[4]) then begin
            adim = (*state).adim
        endif
    endif
    poverlay = (*state).overlay_data_ptr
    RP = reform((*poverlay)[*, *, (*state).slice,adim], (*state).matrix[0], (*state).matrix[1])
    RS = reform((*poverlay)[*, (*state).phase, *,adim], (*state).matrix[0], (*state).matrix[2])
    PS = reform((*poverlay)[(*state).freq, *,  *,adim], (*state).matrix[1], (*state).matrix[2])

    interleave = 1
    
    im_overlay_RP->setProperty, DATA=bytscl(RP, min=(*state).overlay_threshold[0], max=(*state).overlay_threshold[1]), $
                                DIMENSIONS=(size(RP, /dimensions))[0:1], INTERLEAVE=interleave, hide=0
    im_overlay_RS->setProperty, DATA=bytscl(RS, min=(*state).overlay_threshold[0], max=(*state).overlay_threshold[1]), $
                                DIMENSIONS=(size(RS, /dimensions))[0:1], INTERLEAVE=interleave, hide=0
    im_overlay_PS->setProperty, DATA=bytscl(PS, min=(*state).overlay_threshold[0], max=(*state).overlay_threshold[1]), $
                                DIMENSIONS=(size(PS, /dimensions))[0:1], INTERLEAVE=interleave, hide=0
    
    (*state).owin_RP->draw, (*state).ov_RP
    (*state).owin_RS->draw, (*state).ov_RS
    (*state).owin_PS->draw, (*state).ov_PS

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_null_event
;; Created by:
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;; 
;;    Null event that does nothing.
;;
;; Editing Information:

pro ortho_null_event, ev
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_view_cleanup
;; Created by:
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;     Cleans up the window and frees any associated pointers and
;;     graphics objects.
;;
;; Editing Information:

pro ortho_view_cleanup, tlb

    common scan_data

    ;; print, "ortho_refresh_cleanup: arrived."

    widget_control, tlb, get_uvalue=state

    ptr_free, (*state).datar

    obj_destroy, [ (*state).owin_RS, (*state).owin_PS, (*state).owin_RP, $
                   (*state).ov_RP, (*state).ov_PS, (*state).ov_RS ]    

    if (ptr_valid((*state).loaded_ROI)) then begin
        ptr_free, (*state).loaded_ROI
       ;;; not implemented!
    endif

    ptr_free, (*state).marker_list
    ptr_free, state

    if (ptr_valid(project.dataarray[project.ci].state1)) then begin
       project.procpramarray[project.ci].ortho_tlb = -1
    endif

    heap_gc

end

;+
; :Description:
;    Describe the procedure.
;
; :Params:
;    state
;
;
;
; :Author: wtriplett
;-
pro ortho_update_coords_and_int, state

     ;; update intensity display
     if ((*state).data_available) then begin
       int = (*((*state).data_ptr))[(*state).freq, (*state).phase, (*state).slice, $
                                    ((*state).current_imagery eq 0) ? (*state).adim : 0]
       widget_control, widget_info((*state).base_main, find_by_uname='txt_int'),$
         set_value=strcompress(string(int),/remove_all)
                         
       ;; update the freq/phase/slice textboxes with new values
       widget_control, widget_info((*state).base_main, find_by_uname='txt_sl'),$
         set_value=strcompress(string(fix((*state).slice)), /remove_all)
       widget_control, widget_info((*state).base_main, find_by_uname='txt_fq'),$
         set_value=strcompress(string(fix((*state).freq)), /remove_all)
       widget_control, widget_info((*state).base_main, find_by_uname='txt_ph'),$
         set_value=strcompress(string(fix((*state).phase)), /remove_all)
     endif

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_view_update_marker_list
;; Created by:
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;     Updates the marker list and synchronizes with project's
;;     marker list.
;;
;; Editing Information:

pro ortho_update_marker_list, state

    common scan_data
    
    marker_list = (*state).marker_list

    if (not ptr_valid(marker_list)) then begin
       widget_control, (*state).marker_list_wid, set_value=''
       ptr_free, project.procpramarray[(*state).proj_index].marker_list
       return
    endif

    project.procpramarray[(*state).proj_index].marker_list = ptr_new(*marker_list) 

    n_markers = n_elements(*marker_list)

    marker_str = strarr(n_markers)

    for i = 0, n_markers-1 do begin

       marker_name = string(i, format='(I3)') + ': ' + (*marker_list)[i].title
       marker_str[i] = marker_name

    endfor
   
    widget_control, (*state).marker_list_wid, set_value=marker_str

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_marker_event
;; Created by:
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;     Event handler for add/remove/jump markers
;;
;; Editing Information:

pro ortho_marker_event, event

    widget_control, event.top, get_uvalue=state
    
    name = widget_info(event.id, /uname)
    
    case name of
   
        'btn_add_marker': begin
            ;; create a new marker structure
            new_marker = { point: [(*state).freq, (*state).phase, (*state).slice] , $
                           title: '' }
            new_marker.title = '[' + strjoin(string(new_marker.point, format='(I3)'), ',') + ']'
            
            if (ptr_valid((*state).marker_list)) then begin
                marker_list = *((*state).marker_list)
                marker_list = [ marker_list, new_marker ]
                ptr_free, (*state).marker_list
                (*state).marker_list = ptr_new(marker_list)
            endif else begin
                marker_list = [ new_marker ]
                (*state).marker_list = ptr_new(marker_list)
            endelse
            
            ortho_update_marker_list, state
            
            widget_control, (*state).marker_list_wid, set_list_select=n_elements(marker_list)-1
            
        end
        
        'btn_remove_marker': begin
        
            if (not ptr_valid((*state).marker_list)) then return
            
            nmarkers = n_elements(*((*state).marker_list))
            if (nmarkers eq 0) then return
            
            index = widget_info((*state).marker_list_wid, /list_select)
            
            if (nmarkers eq 1) then begin
                ptr_free, (*state).marker_list
                ortho_update_marker_list, state
            endif else begin
                marker = 0
                curr_markers = (*((*state).marker_list))
                for m = 0, nmarkers-1 do begin
                    if (m eq index) then continue
                    new_markers = (marker eq 0) ? curr_markers[m] : [ new_markers, curr_markers[m] ]
                    marker++
                endfor
                ptr_free, (*state).marker_list
                (*state).marker_list = ptr_new(new_markers)
                ortho_update_marker_list, state
                widget_control, (*state).marker_list_wid, set_list_select=(index+1 > 0) < marker-1
           endelse 
                    
        end
        
        'btn_jump_marker': begin
            
            index = widget_info((*state).marker_list_wid, /list_select)
            if (not ptr_valid((*state).marker_list)) then return
            
            point = (*((*state).marker_list))[index].point
            (*state).freq  = point[0]
            (*state).phase = point[1]
            (*state).slice = point[2]
            
            widget_control, widget_info(event.top, find_by_uname='txt_sl'),$
              set_value=strcompress(string(fix((*state).slice)), /remove_all)
            widget_control, widget_info(event.top, find_by_uname='txt_fq'),$
              set_value=strcompress(string(fix((*state).freq)), /remove_all)
            widget_control, widget_info(event.top, find_by_uname='txt_ph'),$
              set_value=strcompress(string(fix((*state).phase)), /remove_all)
            ortho_update_imagery, state
            
         end
        
        else: print, "ortho_marker_event: recv'd event from unknown widget: "+name
        
    endcase
    
end

pro ortho_overlay_event, event

    common scan_data, project
    
    widget_control, event.top, get_uvalue=state
    uname = widget_info(event.id, /uname)
    
    im_overlay_PS = (*state).ov_PS->getByName('main_model/IM_OVL_model/image')
    im_overlay_RS = (*state).ov_RS->getByName('main_model/IM_OVL_model/image')
    im_overlay_RP = (*state).ov_RP->getByName('main_model/IM_OVL_model/image')
    
    case uname of
    
        'slider_overlay_opacity': begin
            
            val = float(event.value)/100.0
        
            im_overlay_RP->SetProperty, alpha_channel=val
            im_overlay_RS->SetProperty, alpha_channel=val
            im_overlay_PS->SetProperty, alpha_channel=val
        
            (*state).owin_RP->draw, (*state).ov_RP
            (*state).owin_RS->draw, (*state).ov_RS
            (*state).owin_PS->draw, (*state).ov_PS
        end
        
        'dl_overlay_setpalette': begin
            loadct, event.index, rgb_table=rgb_table
            rgb_table[0,*] = 0
            opalette = obj_new('idlgrpalette', reform(rgb_table[*,0]), reform(rgb_table[*,1]), reform(rgb_table[*,2]))
            im_overlay_RP->SetProperty, palette=opalette
            im_overlay_RS->SetProperty, palette=opalette
            im_overlay_PS->SetProperty, palette=opalette
            obj_destroy, (*state).overlay_opalette
            (*state).overlay_opalette = opalette
            ortho_update_imagery, state
        end
        
        'dl_set_scan_as_overlay': begin
            id = event.index-1
            if (id lt 0) then begin
                (*state).overlay_data_ptr = ptr_new()
                im_overlay_RP->SetProperty, hide=1
                im_overlay_RS->SetProperty, hide=1
                im_overlay_PS->SetProperty, hide=1
                txt_min = widget_info(event.top, find_by_uname='txt_overlay_min')
                txt_max = widget_info(event.top, find_by_uname='txt_overlay_max')
                widget_control, txt_min, set_value=''
                widget_control, txt_max, set_value=''
                
                ortho_update_imagery, state
                widget_control, event.id, set_uvalue=event.index
                return
            endif
            
            if (not ptr_valid(project.dataarray[id].state1)) then begin
                void = dialog_message('There is no open scan having that ID', /error, /center)
                widget_control, event.id, get_uvalue=old_index
                widget_control, event.id, set_droplist_select=old_index
                return
            endif
            
            overlay_ptr = project.dataarray[id].state1
        
            overlay_dims = size(*overlay_ptr, /dimensions)
            
            if (overlay_dims[0] ne (*state).matrix[0] or $
                overlay_dims[1] ne (*state).matrix[1] or $
                overlay_dims[2] ne (*state).matrix[2]) then begin
                
                void = dialog_message(['The requested scan must have the same x,y,z dimensions', $
                                       'as the currently displayed scan.'], /error, /center)
                widget_control, event.id, get_uvalue=old_index
                widget_control, event.id, set_droplist_select=old_index
                
                return
            endif
            
            (*state).overlay_data_ptr = overlay_ptr
            ol_min = min(*overlay_ptr, max=ol_max)
            (*state).overlay_threshold = [ol_min, ol_max]
            txt_min = widget_info(event.top, find_by_uname='txt_overlay_min')
            txt_max = widget_info(event.top, find_by_uname='txt_overlay_max')
            widget_control, txt_min, set_value=string(ol_min, format='(%"%0.3f")')
            widget_control, txt_max, set_value=string(ol_max, format='(%"%0.3f")')
            ortho_update_imagery, state
            widget_control, event.id, set_uvalue=event.index
            
        end
        
        'txt_overlay_min': begin
            widget_control, event.id, get_value=value
            (*state).overlay_threshold[0] = float(value[0])
            ortho_update_imagery, state
         end

        'txt_overlay_max': begin
            widget_control, event.id, get_value=value
            (*state).overlay_threshold[1] = float(value[0])
            ortho_update_imagery, state
         end
         
         'txt_overlay_reset': begin
            if (not ptr_valid((*state).overlay_data_ptr)) then return
            
            ol_min = min(*(*state).overlay_data_ptr, max=ol_max)
            (*state).overlay_threshold = real_part([ol_min, ol_max])
            txt_min = widget_info(event.top, find_by_uname='txt_overlay_min')
            txt_max = widget_info(event.top, find_by_uname='txt_overlay_max')
            widget_control, txt_min, set_value=string(real_part(ol_min), format='(%"%0.3f")')
            widget_control, txt_max, set_value=string(real_part(ol_max), format='(%"%0.3f")')
            ortho_update_imagery, state
         
         end
         
        else:
        
     endcase
     
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_imagery_event
;; Created by:
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;     Handles events generated by the "Imagery Selection" drop list.
;;
;; Editing Information:

pro ortho_imagery_event, event

    common scan_data

    ;; print, "ortho_imagery_event: arrived."

    widget_control, event.top, get_uvalue=state

    if (event.index gt 0 and event.index lt 4) then begin
       if (project.procpramarray[(*state).proj_index].adt_proccess_flag eq 0) then begin
          junk = dialog_message('Please process ADT before selecting this imagery.', /center, /error)
          widget_control, event.id, SET_DROPLIST_SELECT=0
          return
       endif
    endif
    if event.index eq 4  then begin
      if (project.procpramarray[(*state).proj_index].flow_proccess_flag eq 0) then begin
        junk = dialog_message('Please process flow data before selecting this imagery.', /center, /error)
        widget_control, event.id, SET_DROPLIST_SELECT=0
        return
      endif
    endif
    
    (*state).current_imagery = event.index
    ortho_refresh_data, state
    ortho_update_imagery, state

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_axis_event
;; Created by:
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Editing Information:

pro ortho_axis_event, event

    return

    widget_control, event.top, get_uvalue=state

    case event.index of 
        0: begin
            if ((*state).axis_units eq 'CM') then return
            
            (*state).axis_units = 'CM'            
            tick_multiplier = 10.0

        end

        1: begin
            if ((*state).axis_units eq 'MM') then return
            
            (*state).axis_units = 'MM'
            tick_multiplier = 1.0/10.0

        end

        else: return

    endcase

    ox_axis = (*state).ov_RP->getByName('main_model/RP_xaxis')
    ox_axis->getProperty, tickinterval=tick
    ox_axis->setProperty, tickinterval=tick*tick_multiplier
    
    oy_axis = (*state).ov_RP->getByName('main_model/RP_yaxis')
    oy_axis->getProperty, tickinterval=tick
    oy_axis->setProperty, tickinterval=tick*tick_multiplier
    
    ox_axis = (*state).ov_RS->getByName('main_model/RS_xaxis')
    ox_axis->getProperty, tickinterval=tick
    ox_axis->setProperty, tickinterval=tick*tick_multiplier
    
    oy_axis = (*state).ov_RS->getByName('main_model/RS_yaxis')
    oy_axis->getProperty, tickinterval=tick
    oy_axis->setProperty, tickinterval=tick*tick_multiplier
    
    ox_axis = (*state).ov_PS->getByName('main_model/PS_xaxis')
    ox_axis->getProperty, tickinterval=tick
    ox_axis->setProperty, tickinterval=tick*tick_multiplier
    
    oy_axis = (*state).ov_PS->getByName('main_model/PS_yaxis')
    oy_axis->getProperty, tickinterval=tick
    oy_axis->setProperty, tickinterval=tick*tick_multiplier

    ortho_update_imagery, state

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_zoom_event
;; Created by:
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;     Handles events generated by the zoom slider.
;;
;; Editing Information:

pro ortho_zoom_event, event

    ;print, "ortho_zoom_event: arrived."

    widget_control, event.id, get_value=value
    widget_control, event.top, get_uvalue=state

    name = widget_info(event.id, /uname)
    tmp = stregex(name, '^btn_zoom_([A-Z]{2})([+-])1$', /extract, /subexpr)
    
    pane = tmp[1]
    zoom = tmp[2] eq '+' ? 1 : -1
    
    case pane of
        'RP': begin
            new_zoom = ((*state).zoom + [zoom, 0, 0]) > 1
            view = (*state).ov_RP
            wind = (*state).owin_RP
            ind = 0
        end
        'PS': begin
            new_zoom = ((*state).zoom + [0, zoom, 0]) > 1
            view = (*state).ov_PS
            wind = (*state).owin_PS
            ind = 1
        end
        'RS': begin
            new_zoom = ((*state).zoom + [0, 0, zoom]) > 1
            view = (*state).ov_RS
            wind = (*state).owin_RS
            ind = 2
        end
    endcase
    
    new_dims   = (*state).matrix
    
    ;; Undo whatever the present zoom is and apply the new zoom amount
    (view->getByName('main_model'))->scale, 1.0/(*state).zoom[ind], 1.0/(*state).zoom[ind], 1.
    (view->getByName('main_model'))->scale, new_zoom[ind], new_zoom[ind], 1.

    ;; Save new zoom and update our interal transform matrices
    (*state).zoom = new_zoom 
    ortho_recompute_tx, state
    
    new_dims   = (*state).matrix
    
    ;; update the view obj parameters for new zoom. 
    dims = (invert((*state).tx_RP))[0:1,0:1] ## [new_dims[0:1]]
    dims = fix(abs(dims))
    locs = fix((replicate((*state).maxdim, 3) - dims)/2.) ;; center in viewport    
    
    (*state).ov_RP->setProperty, dimensions=dims, $
                                 viewplane_rect=[0,0, dims[0], dims[1]], $
                                 location=locs

    dims = (invert((*state).tx_PS))[0:1,0:1] ## [new_dims[1],new_dims[2]]
    dims = fix(abs(dims))
    locs = fix((replicate((*state).maxdim, 3) - dims)/2.)

    (*state).ov_PS->setProperty, dimensions=dims, $
                                 viewplane_rect=[0,0, dims[0], dims[1]], $
                                 location=locs

    dims = (invert((*state).tx_RS))[0:1,0:1] ## [new_dims[0],new_dims[2]]
    dims = fix(abs(dims))
    locs = fix((replicate((*state).maxdim, 3) - dims)/2.)
                                 
    (*state).ov_RS->setProperty, dimensions=dims, $
                                 viewplane_rect=[0,0, dims[0], dims[1]], $
                                 location=locs
                                 
    wind->erase, color=(*state).bg_color
    wind->draw, view
    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_rotate_event
;; Created by:
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;     Handles events generated by the rotate buttons. 
;;
;; Editing Information:

pro ortho_rotate_event, event

    name = widget_info(event.id, /uname)
    widget_control, event.top, get_uvalue=state

    new_dims = (*state).matrix * (*state).voxel_dims
   
    case name of 
    
        'btn_rot_RP+90': begin
            ;; move the axis orientation labels appropriately
            AL = shift((*state).axis_lbls_RP, -1)
            widget_control, widget_info((*state).base_main, find_by_uname='lbl_RP_T'), set_value=AL[0]
            widget_control, widget_info((*state).base_main, find_by_uname='lbl_RP_R'), set_value=AL[1]
            widget_control, widget_info((*state).base_main, find_by_uname='lbl_RP_B'), set_value=AL[2]
            widget_control, widget_info((*state).base_main, find_by_uname='lbl_RP_L'), set_value=AL[3]
            (*state).axis_lbls_RP = AL
            ;; ind ~ the panels underconsideration, xy_ind ~ the indices into the dimension array
            ;; for that panel. Used to determine the translation needed to compensate for the rotation.
            ind = 0 & xy_ind=[0,1]
            new_dims *= (*state).zoom[0]
            amt = 90
            view = (*state).ov_RP
            wind = (*state).owin_RP
            modl = view->getByName('main_model')
        end
        
        'btn_rot_PS+90': begin
            AL = shift((*state).axis_lbls_PS, -1)
            widget_control, widget_info((*state).base_main, find_by_uname='lbl_PS_T'), set_value=AL[0]
            widget_control, widget_info((*state).base_main, find_by_uname='lbl_PS_R'), set_value=AL[1]
            widget_control, widget_info((*state).base_main, find_by_uname='lbl_PS_B'), set_value=AL[2]
            widget_control, widget_info((*state).base_main, find_by_uname='lbl_PS_L'), set_value=AL[3]
            (*state).axis_lbls_PS = AL
            ind = 1 & xy_ind=[1,2]
            new_dims *= (*state).zoom[1]
            amt = 90
            view = (*state).ov_PS
            wind = (*state).owin_PS
            modl = view->getByName('main_model')
        end

        'btn_rot_RS+90': begin
            AL = shift((*state).axis_lbls_RS, -1)
            widget_control, widget_info((*state).base_main, find_by_uname='lbl_RS_T'), set_value=AL[0]
            widget_control, widget_info((*state).base_main, find_by_uname='lbl_RS_R'), set_value=AL[1]
            widget_control, widget_info((*state).base_main, find_by_uname='lbl_RS_B'), set_value=AL[2]
            widget_control, widget_info((*state).base_main, find_by_uname='lbl_RS_L'), set_value=AL[3]
            (*state).axis_lbls_RS = AL
            ind = 2 & xy_ind=[0,2]
            new_dims *= (*state).zoom[2]
            amt = 90
            view = (*state).ov_RS
            wind = (*state).owin_RS
            modl = view->getByName('main_model')
        end
        
    endcase
    
    (*state).rotate[ind] = ((*state).rotate[ind] + amt) mod 360
    
    case (*state).rotate[ind] of
        90: begin
            tr_x      = [new_dims[xy_ind[1]], 0]
            view_dims = [new_dims[xy_ind[1]], new_dims[xy_ind[0]]]
            view_vpr  = [0,0,new_dims[xy_ind[1]], new_dims[xy_ind[0]]]
       end
       
       180: begin
            tr_x = [new_dims[xy_ind[0]], 0]
            view_dims = [new_dims[xy_ind[0]], new_dims[xy_ind[1]]]
            view_vpr  = [0,0,new_dims[xy_ind[0]], new_dims[xy_ind[1]]]
       end
       
       270: begin
            tr_x  = [new_dims[xy_ind[1]], 0]
            view_dims = [new_dims[xy_ind[1]], new_dims[xy_ind[0]]]
            view_vpr  = [0,0,new_dims[xy_ind[1]], new_dims[xy_ind[0]]]
       end
       
       0: begin
            tr_x = [new_dims[xy_ind[0]], 0]
            view_dims = [new_dims[xy_ind[0]], new_dims[xy_ind[1]]]
            view_vpr  = [0,0,new_dims[xy_ind[0]], new_dims[xy_ind[1]]]
       end
            
    endcase
    
    new_locs   = fix((replicate((*state).maxdim, 2) - view_dims)/2.) ;; center in viewport    
    
    modl->rotate, [0,0,1], amt
    modl->translate, tr_x[0], tr_x[1], 0
    
    view->setProperty, dimensions=view_dims, $
                        viewplane_rect=view_vpr, $
                        location=new_locs

    ortho_recompute_tx, state
    
    wind->erase, color=(*state).bg_color
    wind->draw, view  

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_flip_event
;; Created by:
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;     Handles events generated by the flip buttons. 
;;
;; Editing Information:

pro ortho_flip_event, event

    widget_control, event.top, get_uvalue=state

    name = widget_info(event.id, /uname)

    dims = (*state).matrix
    
    ;; extract which panel is begin flipped and which direction (H or V)
    tmp = stregex(name, '^btn_flip_([A-Z]{2})_([A-Z])$', $
                  /subexpr, /extract)
    pane      = tmp[1]
    direction = tmp[2]

    ;; prepare to move orientation labels accordingly
    AL_RP = ptr_new((*state).axis_lbls_RP)
    AL_PS = ptr_new((*state).axis_lbls_PS)
    AL_RS = ptr_new((*state).axis_lbls_RS)
    
    case pane of
        'RP': begin
            AL = AL_RP
            new_dims = abs(invert((*state).tx_RP[0:1,0:1]) ## dims[0:1])
            view = (*state).ov_RP
            wind = (*state).owin_RP
            modl = view->getByName('main_model')
         end
        'PS': begin
            AL = AL_PS
            new_dims = abs(invert((*state).tx_PS[0:1,0:1]) ## [dims[1],dims[2]])
            view = (*state).ov_PS
            wind = (*state).owin_PS
            modl = view->getByName('main_model')
         end
        'RS': begin
            AL = AL_RS
            new_dims = abs(invert((*state).tx_RS[0:1,0:1]) ## [dims[0],dims[2]])
            view = (*state).ov_RS
            wind = (*state).owin_RS
            modl = view->getByName('main_model')
         end
     endcase
    
     case direction of
        'H': begin
            tmp = (*AL)[3] & (*AL)[3] = (*AL)[1] & (*AL)[1] = tmp
            widget_control, widget_info((*state).base_main, find_by_uname='lbl_'+pane+'_R'), set_value=(*AL)[1]
            widget_control, widget_info((*state).base_main, find_by_uname='lbl_'+pane+'_L'), set_value=(*AL)[3]
            axis = [0,1,0]
            tran = [new_dims[0], 0, 0]
        end
        'V': begin
            tmp = (*AL)[0] & (*AL)[0] = (*AL)[2] & (*AL)[2] = tmp
            widget_control, widget_info((*state).base_main, find_by_uname='lbl_'+pane+'_T'), set_value=(*AL)[0]
            widget_control, widget_info((*state).base_main, find_by_uname='lbl_'+pane+'_B'), set_value=(*AL)[2]
            axis = [1,0,0]
            tran = [0, new_dims[1], 0]
        end
    endcase

    (*state).axis_lbls_RP = *AL_RP
    (*state).axis_lbls_PS = *AL_PS
    (*state).axis_lbls_RS = *AL_RS
    ptr_free, [AL_RP, AL_PS, AL_RS]
    
    modl->rotate, axis, 180
    modl->translate, tran[0], tran[1], tran[2]
    wind->draw, view
    
    ortho_recompute_tx, state

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_base_event
;; Created by:
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;;     Commented out because it can cause seg faults and very strange
;;     behavior on linux (and possibly other platforms.)
;;
;; Purpose of subroutine:
;;
;;     Handles events generated by the when the main window is
;;     resized. Similar to the zoom feature, this works by passing in
;;     the old geometry structure which is used as a hint to recreate
;;     the GUI with new dimensions.
;;
;; Editing Information:

;pro ortho_base_event, event
;    help, event, /structure
;    return
;end

 pro ortho_base_event, event

     widget_control, event.top, get_uvalue=state

     widget_control, widget_info(event.top, find_by_uname='txt_fq'), get_value=initx
     widget_control, widget_info(event.top, find_by_uname='txt_ph'), get_value=inity
     widget_control, widget_info(event.top, find_by_uname='txt_sl'), get_value=initz

     int = (*state).max_intensity

     tlb = event.top
     dataset = (*state).proj_index
     imagery = (*state).current_imagery

     zoom = (*state).zoom

     geom = widget_info(event.top, /geometry)

     ;;help, geom, /structure
     if (imagery eq 4) then begin
        data_ptr = (*state).data_ptr
     endif

     widget_control, tlb, /destroy

     if (ptr_valid(data_ptr)) then begin
          mas_display_ortho, inittab=1, $
                            initx=fix(initx), inity=fix(inity), initz=fix(initz), $
                            initintensity=int, data_ptr=data_ptr, $
                            initgeometry=geom, initdataset=dataset
     endif else begin
     mas_display_ortho, inittab=1, $
                        initx=fix(initx), inity=fix(inity), initz=fix(initz), $
                        initintensity=int, initimagery=imagery, $
                        initgeometry=geom, initdataset=dataset
    endelse

 end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_data_threshold_event
;; Created by:
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;     Handles events associated with the Min/Max Data Threshold text input boxes
;;
;; Editing Information:

pro ortho_data_threshold_event, event

    ;; data threshold textbox  event
    widget_control, event.top, get_uvalue=state
    widget_control, event.id, get_uvalue = intensity_type
    widget_control, event.id, get_value=value
    if intensity_type eq 0 then begin
        (*state).max_data_thr = float(value[0])
    endif else begin
        (*state).min_data_thr = float(value[0])
    endelse
    
    ortho_refresh_data, state
    ortho_update_imagery, state

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_int_event
;; Created by:
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;     Handles events associated with the Intensity Adjust slider
;;
;; Editing Information:

pro ortho_int_event, event

    ;; intensity slider event
    widget_control, event.top, get_uvalue=state
    widget_control, event.id, get_uvalue = intensity_type
    widget_control, event.id, get_value=value
    if intensity_type eq 0 then begin
        (*state).max_intensity = value; byte(value)
    endif else if intensity_type eq 1 then begin
        (*state).center_intensity = value; byte(value)
    endif else begin
        (*state).min_intensity = value; byte(value)
    endelse
    
    ortho_refresh_data, state
    ortho_update_imagery, state

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_view_event
;; Created by:
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;     Handles events generated as the user moves around the data by
;;     dragging in one of the image planes.
;;
;;      LEFT_CLICK: Normal navigation around adata
;;      RIGHT_CLICK: Pans/Scrolls the data to reveal imagery that may
;;      be outside the viewport boundaries.
;;
;; Editing Information:

pro ortho_view_event, event

    wTarget = event.id
    wWidget = event.top
    
    if (wTarget eq wWidget) then begin
       ortho_base_event, event
       return
    endif

    ;;help, event, /structure

    ;; force event nonnegative
    clk_x = max([0,event.x])
    clk_y = max([0,event.y])

    name = widget_info(wTarget, /uname)

    widget_control, wWidget, get_uvalue=state

    if (event.x eq 0 and event.y eq 0 and event.press eq 0 and event.release eq 0) then begin
        (*state).owin_RP->erase, color=(*state).bg_color
        (*state).owin_RS->erase, color=(*state).bg_color
        (*state).owin_PS->erase, color=(*state).bg_color
        
        (*state).owin_RP->draw, (*state).ov_RP
        (*state).owin_RS->draw, (*state).ov_RS
        (*state).owin_PS->draw, (*state).ov_PS
        return
    endif

    matrix = (*state).matrix

    ;; this takes care of dragging within a view
    if (event.press ne 0) then begin
       (*state).is_pressed = event.press
       (*state).down_x = clk_x
       (*state).down_y = clk_y
    endif
    if (event.release ne 0) then (*state).is_pressed = 0
    
    zoom = float((*state).zoom) * (*state).voxel_dims

    ;; translate the imagery slices to in-data indices
    freq_actual  = (*state).freq
    phase_actual = (*state).phase
    slice_actual = (*state).slice

    ;; has the scan gone away?
    if (not ptr_valid((*state).data_ptr)) then begin
        if (*state).data_available eq 1 then begin
            junk = dialog_message('Scan missing. Some features may no longer be available.', /center)
            (*state).data_available = 0
        endif
    endif

    ;; default intensity displayed
    int = '------'
    dx = ((*state).down_x - clk_x)
    dy = ((*state).down_y - clk_y)

    maxdim = (*state).maxdim

    ;; which image view are they dragging around in?
    case name of 
             
         'draw_RP': begin

             ;; need to translate the viewport (x,y) to in-data (x,y)
             view = (*state).ov_RP
             wind = (*state).owin_RP
             view->getProperty, location=loc, dimensions=matrix
             
             resolved = (*state).tx_RP ## transpose([ clk_x-loc[0], clk_y-loc[1], 1. ])
             
             resolved[0] = min([ max([0, resolved[0]]), (*state).matrix[0]-1 ])
             resolved[1] = min([ max([0, resolved[1]]), (*state).matrix[1]-1 ])
             
             freq_actual  = resolved[0]
             phase_actual = resolved[1]
             
             ;; are they dragging? if so, update the imagery
             if ((*state).is_pressed eq 1 or event.press eq 1) then begin
                 if (resolved[0] ne (*state).freq or resolved[1] ne (*state).phase) then begin
                    ;; only update is they moved to a new voxel location
                    (*state).freq  = resolved[0]
                    (*state).phase = resolved[1]
                    ortho_update_imagery, state, ignore='RP'
                 endif
             endif
         end
         
         ;; last two same as above
         'draw_RS': begin
         
             view = (*state).ov_RS
             wind = (*state).owin_RS
             view->getProperty, location=loc, dimensions=matrix
             
             resolved = (*state).tx_RS ## transpose([ clk_x-loc[0], clk_y-loc[1], 1. ])
             
             resolved[0] = min([ max([0, resolved[0]]), (*state).matrix[0]-1 ])
             resolved[1] = min([ max([0, resolved[1]]), (*state).matrix[2]-1 ])
             
             freq_actual  = resolved[0]
             slice_actual = resolved[1]
             
             if ((*state).is_pressed eq 1 or event.press eq 1) then begin
                 if (resolved[0] ne (*state).freq or resolved[1] ne (*state).slice) then begin
                    (*state).freq  = resolved[0]
                    (*state).slice = resolved[1]
                    ortho_update_imagery, state, ignore='RS'
                 endif
              endif
         end
         
         'draw_PS': begin
             
             view = (*state).ov_PS
             wind = (*state).owin_PS
             view->getProperty, location=loc, dimensions=matrix
      
             resolved = (*state).tx_PS ## transpose([ clk_x-loc[0], clk_y-loc[1], 1. ])
             
             resolved[0] = min([ max([0, resolved[0]]), (*state).matrix[1]-1 ])
             resolved[1] = min([ max([0, resolved[1]]), (*state).matrix[2]-1 ])

             phase_actual = resolved[0]
             slice_actual = resolved[1]

             if ((*state).is_pressed eq 1 or event.press eq 1) then begin
                 if (resolved[0] ne (*state).phase or resolved[1] ne (*state).slice) then begin
                    (*state).phase = resolved[0]
                    (*state).slice = resolved[1]
                    ortho_update_imagery, state, ignore='PS'
                 endif
             endif
          end
         
     endcase

     ;; check if the user is scrolling the image around
     if ( (*state).is_pressed eq 4 or event.press eq 4 ) then begin
         ;; right-{click-drag} => scroll image around in the viewport
         new_loc = [ ( (loc[0]-dx) < 0 ) > (maxdim-matrix[0]) , $
                     ( (loc[1]-dy) < 0 ) > (maxdim-matrix[1]) ]
         if matrix[0] lt maxdim then new_loc[0] = loc[0]
         if matrix[1] lt maxdim then new_loc[1] = loc[1]
       
         view->setProperty, location=new_loc
         (*state).down_x = clk_x
         (*state).down_y = clk_y
         wind->draw, view
     endif
     
     ;; for use with on-image ROI drawing. (not implemented)
     ;;if ( event.modifiers eq 1 and ptr_valid( (*state).loaded_ROI ) ) then begin
     ;;   ( *((*state).loaded_ROI) )[resolved_r, resolved_p, resolved_s] = 1
     ;;endif

     ;; update intensity display
;     if ((*state).data_available) then begin
;       int = (*((*state).data_ptr))[freq_actual, phase_actual, slice_actual, $
;                                    ((*state).current_imagery eq 0) ? (*state).adim : 0]
;       widget_control, widget_info(wWidget, find_by_uname='txt_int'),$
;         set_value=strcompress(string(int),/remove_all)                      
;     endif
;      
;     ;; update the freq/phase/slice textboxes with new values
;     if ((*state).is_pressed or event.press) then begin
;
;         widget_control, widget_info(wWidget, find_by_uname='txt_sl'),$
;           set_value=strcompress(string(fix((*state).slice)), /remove_all)
;         widget_control, widget_info(wWidget, find_by_uname='txt_fq'),$
;           set_value=strcompress(string(fix((*state).freq)), /remove_all)
;         widget_control, widget_info(wWidget, find_by_uname='txt_ph'),$
;           set_value=strcompress(string(fix((*state).phase)), /remove_all)
;
;     endif
     
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_line_data_gui_win_cleanup
;; Created by:
;; Calling Information:
;;
;;     tlb: the top level base of the line profile window
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;     Handles destruction of the line profile plot windows. Either
;;     the case where the user clicks the [X] to close it, or if the
;;     user unchecks the button on the Tools tab.
;;
;; Editing Information:

pro ortho_line_data_gui_win_cleanup, tlb

    widget_control, tlb, get_uvalue=state
    if (ptr_valid(state)) then begin
       drawstate = (*state).user_defined ;; don't free this ptr
                                         ;; it is the main ortho
                                         ;; state ptr!
       if (ptr_valid(drawstate)) then begin
          widget_control, widget_info((*drawstate).base_main,$
                                      find_by_uname='btn_line_data'), $
                          set_button=0
          (*drawstate).line_data_gui_win = -1
          (*drawstate).show_line_profile = 0
       endif
       ;; if the above is not valid then the parent has been closed
       ;; already by the user.
    endif

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_setlink_event
;; Created by:
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;  Handles events related to linking one ortho viewer to another
;;
;; Editing Information:

pro ortho_setlink_event, event

    common scan_data
    
    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        case !ERROR_STATE.name of
            'IDL_M_FILE_CNTOPNDIR': void = dialog_message('Cannot read directory.', /error, /center)
            'IDL_M_FILE_FREADERR':  void = dialog_message('File does not appear to be a valid matrix file.', /error, /center)
            else: void = dialog_message('An unknown error has occurred.', /error, /center)
        endcase
        return
    endif
    
    widget_control, event.top, get_uvalue=state
    name = widget_info(event.id, /uname)
    
    case name of 
    
        'btn_choose_matrix': begin
        
            mat_file = dialog_pickfile(path=project.current_path, /file, /read, filter=['*.mat'])
            if (mat_file eq '') then return
            
            matrix = dblarr(4,4)
            openr, lun, mat_file, /get_lun
            readf, lun, matrix
            close, lun & free_lun, lun
            matrix[*,3] = [0,0,0,1]
            (*state).link_matrix = matrix
            widget_control, widget_info(event.top, find_by_uname='lbl_linkmatrix_file'), set_value=file_basename(mat_file)
            widget_control, widget_info(event.top, find_by_uname='btn_invert_matrix'), set_button=0
            
        end
        
        'dl_set_scanid': begin
        
            if (event.index eq 0) then begin
                (*state).linked_scanid = -1
                return
            endif
            
            id = event.index - 1
            
            if (id eq (*state).proj_index) then begin
                void = dialog_message('Cannot link a scan with itself.', /error, /center)
                widget_control, event.id, set_droplist_select=(*state).linked_scanid+1
                return
            endif
            
            linked_tlb = project.procpramarray[id].ortho_tlb
            
            if (widget_info(linked_tlb, /valid_id) eq 1) then begin
            
                (*state).linked_scanid = id
                
            endif else begin
            
                junk = dialog_message(['An orthogonal viewer must already be open', $
                                        'in order to set up the link.'], /error, /center)
                                        
                (*state).linked_scanid = -1
                (*state).link_matrix = diag_matrix(replicate(1D,4))
                widget_control, event.id, set_droplist_select=0
                                
            endelse
            
        end
        
        'btn_clear_matrix': begin
            ;; the button does not currently exist
            ;; linking is stoped by selecting "None" from the
            ;; droplist
            (*state).link_matrix = diag_matrix(replicate(1D,4))
            ;;(*state).linked_scanid = -1
            widget_control, widget_info(event.top, find_by_uname='lbl_linkmatrix_file'), set_value='<none>'
            widget_control, widget_info(event.top, find_by_uname='btn_invert_matrix'), set_button=0

        end
        
        'btn_invert_matrix': begin
            (*state).link_matrix = invert((*state).link_matrix)
        end
        
    endcase
    
end

pro ortho_set_4d_callback_type_event, event

    widget_control, event.top, get_uvalue=state
 
    if (event.select) then begin
    
        ortho_set_4d_callback_event, event
        
    endif else begin
    
        (*state).navigation_callback = ''
        
    endelse

end

pro ortho_run_4d_callback, event

    widget_control, event.top, get_uvalue=state

    txt_callback_name = widget_info(event.top, find_by_uname='txt_callback_name')
    widget_control, txt_callback_name, get_value=cb_pro_name
    cb_pro_name = cb_pro_name[0]
    
    ;; make sure it looks right
    if (strupcase(strmid(cb_pro_name, 0, 3)) ne 'UCB') then begin
        message, /info, 'Callback procedure name must be prefixed with "UCB" (eg UCB_do_my_stuff)'
        widget_control, event.id, set_button=0
        return
    endif
    
    ;; make sure that the procedure actually exists (ie it is compiled)
    void = where(strmatch(routine_info(), strupcase(cb_pro_name)), nmatch)
    if (nmatch eq 0) then begin
        message, /info, 'Callback procedure not found: '+cb_pro_name
        widget_control, event.id, set_button=0
        return
    endif 
         
    void = where(strmatch(routine_info(), strupcase(cb_pro_name)), nmatch)
    if (nmatch ne 0) then begin
        call_procedure, cb_pro_name, [(*state).freq, (*state).phase, (*state).slice], state
    endif

end


pro ortho_set_4d_callback_event, event

    widget_control, event.top, get_uvalue=state
    
    widget_control, event.id, get_uvalue=cb_pro_name
    
    if (n_elements(cb_pro_name) ne 0) then begin
    
        ;; this code is exclusively for the mas_timeseries_plot button onthe 4D tools tab.
        if (event.select) then begin
            (*state).navigation_callback = cb_pro_name
        endif else begin
            (*state).navigation_callback = ''
        endelse
        
    endif else begin
    
        ;; this code actually works with the "Developer" tab widgets
        if (event.select eq 0) then begin
            (*state).navigation_callback = ''
            return
        endif
        
        ;; get the user's callback name
        txt_callback_name = widget_info(event.top, find_by_uname='txt_callback_name')
        widget_control, txt_callback_name, get_value=cb_pro_name
        cb_pro_name = cb_pro_name[0]
        
        ;; make sure it looks right
        if (strupcase(strmid(cb_pro_name, 0, 3)) ne 'UCB') then begin
            message, /info, 'Callback procedure name must be prefixed with "UCB" (eg UCB_do_my_stuff)'
            widget_control, event.id, set_button=0
            (*state).navigation_callback = ''
            return
        endif
        
        ;; make sure that the procedure actually exists (ie it is compiled)
        void = where(strmatch(routine_info(), strupcase(cb_pro_name)), nmatch)
        if (nmatch eq 0) then begin
            message, /info, 'Callback procedure not found: '+cb_pro_name
            widget_control, event.id, set_button=0
            (*state).navigation_callback = ''
            return
        endif else begin
            message, /info, 'Setting callback procedure: '+cb_pro_name
            (*state).navigation_callback = cb_pro_name
        endelse 
    
    endelse
    
end

pro ortho_run_hnb_ellipsoid, event

    
    widget_control, event.top, get_uvalue=state
    
    x = (*state).freq
    y = (*state).phase
    z = (*state).slice
    
    hnb_ellipsoid_curvefit, x, y, z

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_control_event
;; Created by:
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;; Editing Information:

pro ortho_control_event, event

    common scan_data
    common common_widgets
    forward_function mas_open_line_data_window
    
    ;;help, event, /structure
    
    ;; print, "ortho_control_event: arrived."

    ci = project.ci
    wTarget = event.id
    wWidget = event.top

    widget_control, wWidget, get_uvalue=state
    name = widget_info(wTarget, /uname)
    
    ;; which button was pressed?
    case name of
        'btn_xhair': begin
            (*state).draw_crosshairs = event.select
            ortho_update_imagery, state, ignore=(*state).line_prof_selected_window
        end
        
        'btn_bbox': begin
            ((*state).ov_RP->getByName('main_model/bbox_model'))->setProperty, hide=(1-event.select)
            ((*state).ov_PS->getByName('main_model/bbox_model'))->setProperty, hide=(1-event.select)
            ((*state).ov_RS->getByName('main_model/bbox_model'))->setProperty, hide=(1-event.select)
            (*state).owin_RP->draw, (*state).ov_RP
            (*state).owin_PS->draw, (*state).ov_PS
            (*state).owin_RS->draw, (*state).ov_RS
        end

        'slider_lpmax': begin
           widget_control, event.id, get_value=val
           (*state).line_prof_scale = float(val)/100.
           ortho_update_imagery, state, ignore=(*state).line_prof_selected_window
        end

        'btn_line_data': begin

           (*state).show_line_profile = event.select

           if (event.select) then begin
              line_data_gui_win = mas_open_line_data_window(group_leader=event.top, fov=(*state).fov_fps[0:1])
              widget_control, line_data_gui_win, $
                              kill_notify='ortho_line_data_gui_win_cleanup', $
                              get_uvalue=ld_state
              (*ld_state).user_defined = state
              (*state).line_data_gui_win = line_data_gui_win
           endif else begin
              if (widget_info((*state).line_data_gui_win, /valid_id)) then begin
                 widget_control, (*state).line_data_gui_win, /destroy
                 (*state).line_data_gui_win = 0
              endif
           endelse

           ortho_update_imagery, state, ignore=(*state).line_prof_selected_window

        end

        'btn_update': begin
           ;; someone manually entered freq/phase/slice/array 
           old_freq  = (*state).freq
           old_phase = (*state).phase
           old_slice = (*state).slice
           old_adim  = (*state).adim

           wid_freq  = widget_info(wWidget, find_by_uname='txt_fq')
           wid_phase = widget_info(wWidget, find_by_uname='txt_ph')
           wid_slice = widget_info(wWidget, find_by_uname='txt_sl')
           wid_array = widget_info(wWidget, find_by_uname='txt_ar')

           ;; bounds checking
           widget_control, wid_slice, get_value=new_slice
            if (fix(new_slice) lt 0 or fix(new_slice) gt (*state).matrix[2]-1) then begin
               widget_control, wid_slice, set_value=strcompress(fix(old_slice), /remove_all)
               return
            endif
            
            widget_control, wid_freq, get_value=new_freq
            if (fix(new_freq) lt 0 or fix(new_freq) gt (*state).matrix[0]-1) then begin
               widget_control, wid_freq, set_value=strcompress(fix(old_freq), /remove_all)
               return
            endif

            widget_control, wid_phase, get_value=new_phase
            if (fix(new_phase) lt 0 or fix(new_phase) gt (*state).matrix[1]-1) then begin
               widget_control, wid_phase, set_value=strcompress(fix(old_phase), /remove_all)
               return
            endif

            widget_control, wid_array, get_value=new_adim
            if (fix(new_adim) lt 0 or fix(new_adim) gt project.imndarray[ci].adim-1) then begin
               widget_control, wid_array, set_value=strcompress(fix(old_adim), /remove_all)
               return
            endif

            new_adim = long(new_adim)

            ;; translate to in-data indices
            (*state).freq  = long(new_freq)
            (*state).phase = long(new_phase)
            (*state).slice = long(new_slice)
            
            ;; want new adim?
            if (new_adim eq (*state).adim or $
                new_adim lt 0 or $
                new_adim gt project.imndarray[(*state).proj_index].adim) then begin
                
                widget_control, widget_info(wWidget, find_by_uname='txt_ar'), $
                  set_value=strcompress(string((*state).adim), /remove_all)

                new_adim = (*state).adim

            endif

            ;; set mas's sliders
            project.procpramarray[(*state).proj_index].sdim_start = new_slice
            project.procpramarray[(*state).proj_index].fdim_start = new_freq
            project.procpramarray[(*state).proj_index].pdim_start = new_phase

            ;; update the one pertinent to the current slice axis
            case project.procpramarray[(*state).proj_index].slice_axis of
                0: widget_control, sdim_slider, set_value=strcompress(fix(new_slice), /remove_all)
                1: widget_control, sdim_slider, set_value=strcompress(fix(new_phase), /remove_all)
                2: widget_control, sdim_slider, set_value=strcompress(fix(new_freq), /remove_all)
            endcase
            
            sdim_slider_callback

            ortho_update_imagery, state, want_adim=new_adim

        end
        
        'txt_fq': begin
        
            widget_control, event.id, get_value=val
            
            reg = stregex(val, '[-+=]')
            if (reg ne -1) then begin
                val = strmid(val, reg, 1)
            endif

            if (val eq '+' or val eq '=') then begin
                (*state).freq = ((*state).freq + 1) < ((*state).matrix[0] - 1)
                widget_control, event.id, set_value = strtrim(string((*state).freq),2)
                ortho_update_imagery, state
            endif else if (val eq '-') then begin
                (*state).freq = ((*state).freq - 1) > 0
                widget_control, event.id, set_value = strtrim(string((*state).freq),2)
                ortho_update_imagery, state
            endif
        end
        
        'txt_ph': begin
        
            widget_control, event.id, get_value=val
            
            reg = stregex(val, '[-+=]')
            if (reg ne -1) then begin
                val = strmid(val, reg, 1)
            endif

            if (val eq '+' or val eq '=') then begin
                (*state).phase = ((*state).phase + 1) < ((*state).matrix[1] - 1)
                widget_control, event.id, set_value = strtrim(string((*state).phase),2)
                ortho_update_imagery, state
            endif else if (val eq '-') then begin
                (*state).phase = ((*state).phase - 1) > 0
                widget_control, event.id, set_value = strtrim(string((*state).phase),2)
                ortho_update_imagery, state
            endif
        end
        
        'txt_sl': begin
        
            widget_control, event.id, get_value=val
            
            reg = stregex(val, '[-+=]')
            if (reg ne -1) then begin
                val = strmid(val, reg, 1)
            endif

            if (val eq '+' or val eq '=') then begin
                (*state).slice = ((*state).slice + 1) < ((*state).matrix[2] - 1)
                widget_control, event.id, set_value = strtrim(string((*state).slice),2)
                ortho_update_imagery, state
            endif else if (val eq '-') then begin
                (*state).slice = ((*state).slice - 1) > 0
                widget_control, event.id, set_value = strtrim(string((*state).slice),2)
                ortho_update_imagery, state
            endif
            
        end

        'txt_ar': begin
        
            widget_control, event.id, get_value=val
            
            reg = stregex(val, '[-+=]')
            if (reg ne -1) then begin
                val = strmid(val, reg, 1)
            endif

            if (val eq '+' or val eq '=') then begin
                new_adim = ((*state).adim + 1) mod (project.imndarray[(*state).proj_index].adim)
                widget_control, event.id, set_value = strtrim(string(new_adim),2)
                ortho_update_imagery, state, want_adim=new_adim
            endif else if (val eq '-') then begin
                new_adim = ((*state).adim - 1)
                if (new_adim lt 0) then new_adim = project.imndarray[(*state).proj_index].adim - 1
                widget_control, event.id, set_value = strtrim(string(new_adim),2)
                ortho_update_imagery, state, want_adim=new_adim
            endif
            
        end
        
        else: begin
           print, "ortho_control_event: recv'd event from unknown widget: "+name
        end

    endcase
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_make_panel
;; Created by:
;; Calling Information:
;;
;;  VOXEL_DIMS:  2 element array specifying the voxel dimensions 
;;               as [xdim, ydim]
;;
;;  VIEWPORT_DIMS: 2 element array specifying the image display 
;;                 area within the view panel as [xdim, ydim]
;;
;;  PANEL_DIMS: 2 element array specifying the dimensions of the
;;              panel base. The viewport will be centered in the
;;              panel.
;;              
;;  AXIS_LABELS: 4 elements array specifying the axis labels in 
;;               clockwise order.
;;
;;  DESIGNATION: a text string speficying the "designation" of this panel. 
;;               It is used to uniquely identify this panel and its
;;               sub-elements
;;
;;  PARENT_BASE: the widget ID of the WIDGET_BASE that will enclose this panel
;;
;;  RENDERER: this is passed to IDLgrView indicating whether to use hardware
;;            openGL rendering or software rendering
;;
;;  draw_wid: (OUTPUT) set this to a named variable to receive the widget ID 
;;             of the draw widget for this panel
;;
;;  view_obj: (OUTPUT) set this to a naked variable to receive the IDLgrView object
;;             for the newly created panel
;;              
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;     Handles the creation of a single view panel
;;
;; Editing Information:

pro ortho_make_panel, $
    voxel_dims=voxel_dims, $
    viewport_dims=viewport_dims, $
    panel_dims=panel_dims, $
    axis_labels=axis_labels, $ 
    designation=designation, $
    parent_base=parent_base, $
    renderer=renderer, $
    draw_wid=draw_wid, $
    view_obj=view_obj, stat_view=stat_view
    
    des_base_ctr = 'base_'+designation+'_ctr'
    des_draw     = 'draw_'+designation
    des_label_T  = 'lbl_'+designation+'_T'
    des_label_R  = 'lbl_'+designation+'_R'
    des_label_B  = 'lbl_'+designation+'_B'
    des_label_L  = 'lbl_'+designation+'_L'

    if (not keyword_set(axis_labels)) then begin
        axis_labels = ['','','','']
    endif
    
    if (not keyword_set(renderer)) then begin
        renderer = 1
    endif
    
    ;; creates main view panel with labels for orientation letters
    panel_base     = widget_base(parent_base, /row, /frame, /align_center, xpad=0)
    orient_lbl_L   = widget_label(panel_base, value=axis_labels[3], uname=des_label_L)
    panel_base_ctr = widget_base(panel_base, /column, uname=des_base_ctr)
    orient_lbl_T   = widget_label(panel_base_ctr, value=axis_labels[0], uname=des_label_T)
    base_draw      = widget_base(panel_base_ctr, /align_center)
    draw_wid       = widget_draw(base_draw, uname=des_draw, $
                                 BUTTON_EVENTS=1, $
                                 MOTION_EVENTS=1, $
                                 EXPOSE_EVENTS=1, $
                                 GRAPHICS_LEVEL=2, $
                                 retain=0, $
                                 renderer=1, $
                                 xsize=viewport_dims[0], ysize=viewport_dims[1], $
                                 scr_xsize=panel_dims[0], scr_ysize=panel_dims[1], $
                                 /ALIGN_CENTER)
    orient_lbl_B     = widget_label(panel_base_ctr, value=axis_labels[2], uname=des_label_B)
    orient_lbl_R     = widget_label(panel_base, value=axis_labels[1], uname=des_label_R)
        
    mid_location = (panel_dims - viewport_dims)/2.

    des_main_model   = 'main_model'
    des_IM_model     = 'IM_model'
    des_IM_OVL_model = 'IM_OVL_model'
    des_image        = 'image'
    des_CH_model     = 'CH_model'
    des_stat_model   = 'STAT_model'
    des_crosshair_H  = 'crosshair_H'
    des_crosshair_V  = 'crosshair_V'
    des_bbox_model   = 'bbox_model'
    des_bbox         = 'bbox'
    
    ;; main view
    view_obj = obj_new('idlgrview', $
                       dimensions=viewport_dims, $
                       location=mid_location, $
                       color=[128,128,128], $
                       viewplane_rect=[0,0,viewport_dims[0],viewport_dims[1]])
    ;; main model
    main_model_obj = obj_new('idlgrmodel', name=des_main_model, depth_test_disable=2)
    ;; anisotropic voxel dimension handled here
    main_model_obj->scale, voxel_dims[0], voxel_dims[1], 1.
    
    ;; image model
    IM_model_obj     = obj_new('idlgrmodel', name=des_IM_model, depth_test_disable=2)
    IM_OVL_model_obj = obj_new('idlgrmodel', name=des_IM_OVL_model, depth_test_disable=2)
    xdim = (viewport_dims[0]/voxel_dims[0])>1
    ydim = (viewport_dims[1]/voxel_dims[1])>1
    image = reform(bytarr(xdim,ydim), xdim, ydim) 
    img_obj = obj_new('idlgrimage', name=des_image, transform_mode=1, image, palette=obj_new())
    IM_model_obj->add, img_obj
    
    img_obj_ovl = obj_new('idlgrimage', name=des_image, transform_mode=1, blend_function=[3,1], $
                          image, /hide)
    IM_OVL_model_obj->add, img_obj_ovl
    
    main_model_obj->add, [IM_model_obj, IM_OVL_model_obj]
    
    ;; bounding box
    bbox_model = obj_new('idlgrmodel', name=des_bbox_model)
    bbox_w = viewport_dims[0]/voxel_dims[0] - 0.1
    bbox_h = viewport_dims[1]/voxel_dims[1] - 0.1
    bbox = obj_new('idlgrpolyline', [ [0.1,0.1,0], $
                                      [bbox_w, 0.1, 0], $
                                      [bbox_w, bbox_h, 0], $
                                      [0.1, bbox_h, 0] , [0.1,0.1,0]], $
                   color=[255,255,255], alpha_channel=0.4)
    
    bbox_model->add, bbox
    main_model_obj->add, bbox_model
    
    ;; HEADSUP DISPLAY
    CH_model_obj = obj_new('idlgrmodel', name=des_CH_model)
    CH_model_obj->add, obj_new('idlgrpolyline', color=[255,255,255], alpha=1.0, name=des_crosshair_H)
    CH_model_obj->add, obj_new('idlgrpolyline', color=[255,255,255], alpha=1.0, name=des_crosshair_V)
    main_model_obj->add, CH_model_obj
    
    view_obj->add, main_model_obj
    
    ;; note that on-display stats not implemented
;    stat_view = obj_new('idlgrview', /transparent, viewplane_rect=[0,0,50,50], dimensions=[50,50], name=des_stat_view)
;    stat_model = obj_new('idlgrmodel', name=des_stat_model)
;    case designation of
;        'RP': cir = circle(50,50,0, slice_axis=0)
;        'RS': cir = circle(50,50,1, slice_axis=1)
;        'PS': cir = circle(50,50,2, slice_axis=2)
;    endcase
;    cir_image = obj_new('idlgrimage', cir, interleave=2)
;    stat_model->add, cir_image
;    stat_view->add, stat_model
    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: ortho_recompute_tx
;; Created by:
;; Calling Information:
;;
;;    state: the pointer to the current ortho_display state
;;
;; Bugs or Important Comments to Developers:
;;
;;     This could be eliminated and we could just rely on the
;;     model's transformation matrix.
;;
;; Purpose of subroutine:
;;
;;     Recomputes the internal transformation matrix for each
;;     panel. Based on the model object's internal tx matrix,
;;     but we simplify it a bit. 
;;
;; Editing Information:

pro ortho_recompute_tx, state
 
    tol = 1e-4
    
    ((*state).ov_RP->getByName('main_model'))->getProperty, transform=tx
    toosmall = where(abs(tx) lt tol, n_toosmall)
    if (n_toosmall gt 0) then tx[toosmall] = 0.
    ((*state).ov_RP->getByName('main_model'))->setProperty, transform=tx
    tmp = [ tx[0,0:2], tx[1,0:2], tx[3,0:2] ] & tmp[2,2] = 1.
    tmp = invert(tmp)
    (*state).tx_RP = tmp

    ((*state).ov_PS->getByName('main_model'))->getProperty, transform=tx
    toosmall = where(abs(tx) lt tol, n_toosmall)
    if (n_toosmall gt 0) then tx[toosmall] = 0.
    ((*state).ov_PS->getByName('main_model'))->setProperty, transform=tx
    tmp = [ tx[0,0:2], tx[1,0:2], tx[3,0:2] ] & tmp[2,2] = 1.
    tmp = invert(tmp)
    toosmall = where(abs(tmp) lt tol, n_toosmall)
    if (n_toosmall gt 0) then tmp[toosmall] = 0.
    (*state).tx_PS = tmp

    ((*state).ov_RS->getByName('main_model'))->getProperty, transform=tx
    toosmall = where(abs(tx) lt tol, n_toosmall)
    if (n_toosmall gt 0) then tx[toosmall] = 0.
    ((*state).ov_RS->getByName('main_model'))->setProperty, transform=tx
    tmp = [ tx[0,0:2], tx[1,0:2], tx[3,0:2] ] & tmp[2,2] = 1.
    tmp = invert(tmp)
    toosmall = where(abs(tmp) lt tol, n_toosmall)
    if (n_toosmall gt 0) then tmp[toosmall] = 0.
    (*state).tx_RS = tmp

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Subroutine name: mas_display_ortho
;; Created by:
;; Calling Information:
;;
;;     GROUP_LEDAER: the window group leader. Defaults to MAS's WID_BASE_MAIN
;;     VOXDIMS: a 3 elements floating point array indicating the voxel dimensions
;;              (r, p, s) or (x, y, z). Not currently used -- just a placeholder
;;              for when this becomes stand alone.
;;     INITX: The initial position of the x (READ) location
;;     INITY: The initial position of the y (PHASE) location
;;     INITZ: The initial position of the z (SLICE) location
;;
;;     INITIMAGERY: the initial imagery to show to the user
;;
;;     INITINTENSITY: the initial intensity factor to use
;;
;;     INITIAB: the initially selected tab on the control panel
;;
;;     INITGEOMETRY:  a WIDGET_GEOMETRY structure to define the
;;                    initial window geometry.
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;     Creates the window and widgets, sets up the user interface.
;;
;; Editing Information: 1/13/10 - Added Min, Center, and Max intensity sliders to imagery tab (GA)

pro mas_display_ortho, group_leader=group_leader, voxdims=voxdims, data_ptr=data_ptr, $
                       initx=initx, inity=inity, initz=initz, $; init_zoom=init_zoom, $ ;; zoom is not working properly
                       inittab=inittab, initimagery=initimagery, initintensity=initintensity, $
                       initgeometry=initgeometry, initdataset=initdataset, block=block, state = state

    common scan_data
    common common_widgets
    
    ;;print, "mas_display_ortho: arrived."
    
    if (arg_present(initdataset) and n_elements(initdataset) ne 0) then begin
        ci = initdataset
    endif else begin
        ci = project.ci
    endelse
        
    if (not ptr_valid(data_ptr)) then begin
        if (project.procpramarray[ci].ortho_tlb gt 0) then begin
            if (widget_info(project.procpramarray[ci].ortho_tlb, /valid_id)) then begin
                widget_control, project.procpramarray[ci].ortho_tlb , /show
                return
            endif
        endif
        if (not project.procpramarray[ci].state_1) then mas_load_state_1
        data_ptr = project.dataarray[ci].state1
    endif else begin
        initimagery = 5 ;; 5 indicates that the data displayed is not the original scan data
                        ;; or any of the built-in diffusion/flow data. It is some other derived data.
    endelse

    if ( (size(*data_ptr))[0] lt 3 ) then begin
       void = dialog_message(['Orthogonal viewer only applicable to', $
                              'data sets with 3 spatial dimensions.'], /center, /error)
       return
    endif

    min_data = min(*data_ptr, max=max_data)
    min_data = real_part(min_data)
    max_data = real_part(max_data)
    adim   = (0 < project.procpramarray[ci].adim_start);; < (size(*data_ptr, /dimensions))[4]

    initintensity = (keyword_set(initintensity)) ? fix(initintensity) :  255
    initimagery   = (keyword_set(initimagery))   ? initimagery        :  0

    interp_f = project.procpramarray[ci].freq_interp
    interp_p = project.procpramarray[ci].phase_interp
    interp_s = project.procpramarray[ci].slice_interp

    fov_f = project.imndarray[ci].f_fov
    fov_p = project.imndarray[ci].p_fov
    fov_s = project.imndarray[ci].s_fov

    voxdims_actual = [project.imndarray[ci].f_voxsz/interp_f,$
                      project.imndarray[ci].p_voxsz/interp_p,$
                      project.imndarray[ci].s_voxsz/interp_s]

    voxdims_display = voxdims_actual / min(voxdims_actual)
    voxsz_f = voxdims_display[0]
    voxsz_p = voxdims_display[1]
    voxsz_s = voxdims_display[2]

    initzoom =  [1, 1, 1] 

    title  = project.imndarray[ci].display_name
    title  = string(ci, format='(I0)') + ': ' + title
    
    dims = size(*data_ptr, /dimensions)
    n_dims = n_elements(dims)

    initx = (keyword_set(initx)) ? fix(initx) : dims[0]/2.0
    inity = (keyword_set(inity)) ? fix(inity) : dims[1]/2.0
    initz = (keyword_set(initz)) ? fix(initz) : dims[2]/2.0

    fdim_eff = floor(dims[0] * voxsz_f) * initzoom[0]
    pdim_eff = floor(dims[1] * voxsz_p) * initzoom[1]
    sdim_eff = floor(dims[2] * voxsz_s) * initzoom[2]

    base_geom = create_struct(name='WIDGET_GEOMETRY')
    device, get_screen_size=screen_size
    maxdim = screen_size[0]/3.0 ;; 900/3. - 50

    ;; begin ---
    ;; This works with the ortho_base_event
    ;; procedure. the maxdim = line specifies a size for the new
    ;; widget based on the resized value. The /tbl_size_events flag
    ;; causes events to be generated when the user resizes the window.
    if (keyword_set(initgeometry)) then begin
       geom_info = size(initgeometry, /structure)
       if (geom_info.structure_name eq 'WIDGET_GEOMETRY') then begin
          base_geom = initgeometry
          maxdim = base_geom.scr_xsize/3. - 35
       endif
    endif

    group_l = WID_BASE_MAIN
    if (keyword_set(group_leader)) then begin
       if (widget_info(group_leader, /valid_id)) then begin
          group_l = group_leader
       endif
    endif

    base_main = widget_base(title=title, /align_center, /tlb_size_events, $
                            column=1, group_leader=group_l, $
                            xoffset=base_geom.xoffset-(base_geom.xpad+base_geom.space), $
                            yoffset=base_geom.yoffset-(9*base_geom.ypad))
    ;; --- end

    base = widget_base(base_main, row=1, xpad=0)
    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    axis_labels_RP = mas_orient_get_axis_labels(0, orient_string=os)
    ortho_make_panel, axis_labels=axis_labels_RP, designation='RP', $
                      parent_base=base, $
                      voxel_dims=[voxdims_display[0], voxdims_display[1]], $
                      panel_dims=[maxdim, maxdim], $
                      viewport_dims=[fdim_eff, pdim_eff], $
                      draw_wid=draw_RP, $
                      view_obj=ov_RP;;, stat_view=sv_RP

    axis_labels_PS = mas_orient_get_axis_labels(2, orient_string=os)
    ortho_make_panel, axis_labels=axis_labels_PS, designation='PS', $
                      parent_base=base, $
                      voxel_dims=[voxdims_display[1], voxdims_display[2]], $
                      panel_dims=[maxdim, maxdim], $
                      viewport_dims=[pdim_eff, sdim_eff], $
                      draw_wid=draw_PS, $
                      view_obj=ov_PS;;, stat_view=sv_PS

    axis_labels_RS = mas_orient_get_axis_labels(1, orient_string=os)
    ortho_make_panel, axis_labels=axis_labels_RS, designation='RS', $
                      parent_base=base, $
                      voxel_dims=[voxdims_display[0], voxdims_display[2]], $
                      panel_dims=[maxdim, maxdim], $
                      viewport_dims=[fdim_eff, sdim_eff], $
                      draw_wid=draw_RS, $
                      view_obj=ov_RS;;, stat_view=sv_RS
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    base_control = widget_tab(base_main, event_pro='ortho_null_event')

    info_base = widget_base(base_control, title="Slice Information", $
                            row=1, /align_left, /base_align_center)
    
    base_int = widget_base(info_base, row=1, /align_left, frame=0, /base_align_center)
    lbl_int = widget_label(base_int, value="Scan Intensity:", /align_left)
    txt_int = widget_text(base_int, value='----------', xsize=12, uname='txt_int', /align_left)

    lbl_fq = widget_label(info_base, value='Freq:')
    txt_fq = widget_text(info_base, value='40', xsize=3, uname='txt_fq', /editable, /all_events, event_pro='ortho_control_event')
    lbl_ph = widget_label(info_base, value='Phase:')
    txt_ph = widget_text(info_base, value='40', xsize=3, uname='txt_ph', /editable, /all_events, event_pro='ortho_control_event')
    lbl_sl = widget_label(info_base, value='Slice:')
    txt_sl = widget_text(info_base, value='40', xsize=3, uname='txt_sl', /editable, /all_events, event_pro='ortho_control_event')
    lbl_ar = widget_label(info_base, value='Array:')
    txt_ar = widget_text(info_base, value='40', xsize=3, uname='txt_ar', /editable, /all_events, event_pro='ortho_control_event')
    btn_update = widget_button(info_base, value='Update', uname='btn_update', event_pro='ortho_control_event')

    ;; TOOLS Tab
    info_baser2 = widget_base(base_control, title="Tools", $
                              row=1, /align_left, /base_align_center)

    btn_base = widget_base(info_baser2, /nonexclusive, /align_left, /column)
    btn_xhair = widget_button(btn_base, value='Show Crosshairs', uname='btn_xhair', event_pro='ortho_control_event')
    btn_bbox = widget_button(btn_base, value='Show Bounding Box', uname='btn_bbox', event_pro='ortho_control_event')
    
    widget_control, btn_bbox, /set_button
    widget_control, btn_xhair, /set_button
     
    btn_base = widget_base(info_baser2, /nonexclusive, /align_left, /column)
    btn_line_data = widget_button(btn_base, value='Show Line Data', uname='btn_line_data', event_pro='ortho_control_event')
        
    link_base = widget_base(info_baser2, row=2, xpad=1,  ypad=1, /frame)
    link_base_r1 = widget_base(link_base, /row, xpad=0, ypad=0)
    
    lbl = widget_label(link_base_r1, value="Link with scan #", sensitive=(initimagery eq 4) ? 0 : 1)
    dl_linked_with = widget_droplist(link_base_r1, value=['None', string(indgen(20), format='(I3)')], uname="dl_set_scanid", $
                                        event_pro='ortho_setlink_event', sensitive=(initimagery eq 4) ? 0 : 1)
    link_base_r2 = widget_base(link_base, /row, xpad=0, ypad=0)
    lbl = widget_label(link_base_r2, value="Matrix File:", sensitive=(initimagery eq 4) ? 0 : 1)
    lbl_linkmatrix = widget_label(link_base_r2, value='<none>', /dynamic_resize, uname='lbl_linkmatrix_file')
    btn_choose = widget_button(link_base_r2, value="Choose", uname="btn_choose_matrix", $
                                 event_pro='ortho_setlink_event', sensitive=(initimagery eq 4) ? 0 : 1)
    btn_invert_base = widget_base(link_base_r2, /nonexclusive)
    btn_invert_matrix = widget_button(btn_invert_base, value='Invert', uname='btn_invert_matrix',event_pro='ortho_setlink_event', sensitive=(initimagery eq 4) ? 0 : 1)
    
;    btn_clear = widget_button(link_base_r2, value="Clear", uname="btn_clear_matrix", $
;                                 event_pro='ortho_setlink_event')
    ;; TRANSFORM Tab
    info_baser2 = widget_base(base_control, title="Transform", $
                              row=1, /align_left, /base_align_center)
                              
    base_rotate = widget_base(info_baser2, /row)
    
    panes = [ 'RP', 'PS', 'RS' ] & pos = [ 'Left', 'Center', 'Right' ]
    for p = 0, 2 do begin

        base_tx = widget_base(base_rotate, /column, /frame)
        lbl = widget_label(base_tx, value=pos[p]+" Panel")
        btns = widget_base(base_tx, /row)
        btn_zoom = widget_button(btns, value="Zoom +1", $
                                    uname="btn_zoom_"+panes[p]+"+1", $
                                    event_pro='ortho_zoom_event')
        btn_zoom = widget_button(btns, value="Zoom -1", $
                                    uname="btn_zoom_"+panes[p]+"-1", $
                                    event_pro='ortho_zoom_event')
        btn_rot_90 = widget_button(btns, value="Rot 90", $
                                      uname="btn_rot_"+panes[p]+"+90", $
                                      event_pro='ortho_rotate_event')
        btn_flip_H = widget_button(btns, value="Flip H", $
                                     uname="btn_flip_"+panes[p]+"_H", $
                                     event_pro='ortho_flip_event')
        btn_flip_V = widget_button(btns, value="Flip V", $
                                     uname="btn_flip_"+panes[p]+"_V", $
                                     event_pro='ortho_flip_event')
    endfor

    info_baser2 = widget_base(base_control, title="Imagery", row=1, event_pro='ortho_int_event')
    thr_base = widget_base(info_baser2, /column, /base_align_center, ypad=0)
    lbl = widget_label(thr_base, value='Available Imagery')
    dl_imagery = widget_droplist(thr_base, $
                                 value=['Scan Data', $
                                        'Fractional Anisotropy', $
                                        'Average Diffusivity', $
                                        'FA + EV Orientation', $
                                        'Flow Speed'], $
                                 event_pro='ortho_imagery_event')
    if (initimagery eq 4) then begin
        widget_control, dl_imagery, set_value=['Other / Derived'], sensitive = 0
    endif else begin
        widget_control, dl_imagery, set_droplist_select=initimagery, sensitive=1
    endelse

    ;; Text input boxes that determine min/max data thresholds
    thr_base = widget_base(info_baser2, /column, /base_align_center)
    txt_min_data_thr = widget_text(thr_base, xsize=10, value=strtrim(string(min_data),2), $
                      event_pro='ortho_data_threshold_event', /editable, $
                      uname='txt_min_data_thr', uvalue=1)
    lbl_min_data_thr = widget_label(thr_base, value='Min Data Thresh.')

    thr_base = widget_base(info_baser2, /column, /base_align_center)
    txt_max_data_thr = widget_text(thr_base, xsize=10, value=strtrim(string(max_data),2), $
                      event_pro='ortho_data_threshold_event', /editable, $
                      uname='txt_max_data_thr', uvalue=0)
    lbl_max_data_thr = widget_label(thr_base, value='Max Data Thresh.')

    ;; Sliders to determine post-threshold intensity scaling.
    slider_max_int = widget_slider(info_baser2, min=0, max=255, $
                               value=initintensity,uvalue=0, $
                               uname='slider_max_int', title='Max Intensity', $
                               event_pro='ortho_int_event')
    
    slider__center_int = widget_slider(info_baser2, min=0, max=255, $
                               value=initintensity,uvalue=1, $
                               uname='slider_center_int', title='Center Intensity', $
                               event_pro='ortho_int_event')
    
    slider__min_int = widget_slider(info_baser2, min=0, max=255, $
                               value=0, uvalue=2, $
                               uname='slider_min_int', title='Min Intensity', $
                               event_pro='ortho_int_event')
                               
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    info_baser2 = widget_base(base_control, title="Overlay", /row)
    b = widget_base(info_baser2, /row, /frame)
    lbl = widget_label(b, value='Use open scan:')
    dl_link_scan = widget_droplist(b, value=['None', string(indgen(20), format='(I3)')], uname="dl_set_scan_as_overlay", $
                                        event_pro='ortho_overlay_event', sensitive=1, uvalue=0L)
    b = widget_base(info_baser2, /frame, /row)
    sl_overlay_opacity = widget_slider(b, title='Opacity (%)', $
                                       min=0, max=100, value=100, /drag, $
                                       uname='slider_overlay_opacity', $
                                       event_pro='ortho_overlay_event')
    b1 = widget_base(b, /row)
    l = widget_label(b1, value='Threshold Min:')
    t = widget_text(b1, xsize=6, /editable, sensitive=1, uname='txt_overlay_min', event_pro='ortho_overlay_event')
    b1 = widget_base(b, /row)
    l = widget_label(b1, value='Max:')
    t = widget_text(b1, xsize=6, /editable, sensitive=1, uname='txt_overlay_max', event_pro='ortho_overlay_event')
    btn = widget_button(b1, value="Reset", uname='txt_overlay_reset', event_pro='ortho_overlay_event')
    
    b = widget_base(info_baser2, /row, /frame)
    lbl = widget_label(b, value='Color Scheme:')
    loadct, get_names=ct_names
    dl_overlay_setpalette = widget_droplist(b, value=ct_names, uname='dl_overlay_setpalette', event_pro='ortho_overlay_event')
    

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    info_baser2 = widget_base(base_control, title="Data Markers", /row)
    
    marker_list_base = widget_base(info_baser2, /row)
    marker_list = widget_list(marker_list_base, value='', uname="list_markers", xsize=25, ysize=3)

    btn_base = widget_base(info_baser2, /row)
    btn_add    = widget_button(btn_base, value="Add Marker", uname='btn_add_marker', event_pro='ortho_marker_event')
    btn_remove = widget_button(btn_base, value="Remove Marker", uname='btn_remove_marker', event_pro='ortho_marker_event')
    btn_jump   = widget_button(btn_base, value="Jump to Marker", uname='btn_jump_marker', event_pro='ortho_marker_event')
    
    marker_list_curr = ptr_valid(project.procpramarray[ci].marker_list) ? ptr_new(*project.procpramarray[ci].marker_list) : ptr_new()

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;    

    if (n_elements(dims) gt 3 && dims[3] gt 1) then begin
        info_base2 = widget_base(base_control, title='4D Tools', /row)
        bb = widget_base(info_base2, /column)
        b = widget_base(bb, /row)
        btn_b = widget_base(b, /nonexclusive)
        btn_enable_ts = widget_button(btn_b, value='Enable Timeseries Plot', uvalue='mas_timeseries_plot', event_pro='ortho_set_4d_callback_event')
        btn_run_hb = widget_button(b, value="Run HnB Ellipsoid Fit", event_pro='ortho_run_hnb_ellipsoid')
    endif
    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    info_base2 = widget_base(base_control, title='Developer', /row)
    bb = widget_base(info_base2, /column, xpad=0, ypad=0)
    b = widget_base(bb, /column, xpad=0, ypad=0)
    txt_b = widget_base(b, /row, xpad=0, ypad=0)
    lbl = widget_label(txt_b, value='Name of callback procedure:')
    txt_callback_name = widget_text(txt_b, xsize=30, value="", uname='txt_callback_name', /editable)
    btn_enable_ts = widget_button(txt_b, value='Run now', uname='btn_callback_run_now', event_pro='ortho_run_4d_callback')  
    btn_b = widget_base(b, /nonexclusive, xpad=0, ypad=0)
    btn_enable_ts = widget_button(btn_b, value='Automatically run when crosshairs are moved', uname='btn_enable_callback', event_pro='ortho_set_4d_callback_event')
    
    
    ;btn_b = widget_base(b, /nonexclusive)
    ;btn_enable_ts = widget_button(btn_b, value='Run on mouse event', uname='btn_callback_run_on_mouse_event', event_pro='ortho_set_4d_callback_type_event')
    ;btn_b = widget_base(b)
      


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    widget_control, base_main, /realize
    widget_control, base_control, set_tab_current=(keyword_set(inittab)) ? inittab : 0
    widget_control, draw_RP, get_value=owin_RP
    widget_control, draw_PS, get_value=owin_PS
    widget_control, draw_RS, get_value=owin_RS

    loadct, 0, rgb_table=rgb_table
    opalette = obj_new('idlgrpalette', reform(rgb_table[*,0]), reform(rgb_table[*,1]), reform(rgb_table[*,2]))

    state = ptr_new({$
            base_main: base_main, $
            linked_scanid: -1, $
            link_matrix: diag_matrix(replicate(1D,4)), $
            loaded_ROI: ptr_new(), $ ;;ptr_new(bytarr(dims[0:2])), $
            bg_color: [0,0,0], $
            owin_PS: owin_PS,    $
            owin_RS: owin_RS,    $
            owin_RP: owin_RP,    $
            tx_RP  : diag_matrix([voxdims_display[0], voxdims_display[1], 1.]), $
            tx_RS  : diag_matrix([voxdims_display[0], voxdims_display[2], 1.]), $
            tx_PS  : diag_matrix([voxdims_display[1], voxdims_display[2], 1.]), $
            axis_lbls_RP: axis_labels_RP, $
            axis_lbls_PS: axis_labels_PS, $
            axis_lbls_RS: axis_labels_RS, $
            ov_PS  : ov_PS,      $
            ov_RS  : ov_RS,      $
            ov_RP  : ov_RP,      $
            $;sv_RP  : sv_RP, $
            $;sv_PS  : sv_PS, $
            $;sv_RS  : sv_RS, $
            slice  : fix(initz), $
            freq   : fix(initx), $
            phase  : fix(inity), $
            adim   : -1, $            
            voxel_dims_actual: voxdims_actual, $
            voxel_dims: voxdims_display, $
            fov_fps: voxdims_actual * dims[0:2], $
            matrix : (dims[0:2]), $
            maxdim : maxdim, $
            zoom   : initzoom, $
            rotate : [0,0,0.], $
            datar  : ptr_new(), $
            data_ptr: data_ptr, $ ;;ptr_new(), $
            data_max: 0.0D, $
            overlay_data_ptr: ptr_new(), $
            overlay_opalette: opalette, $
            overlay_threshold: [-1.0,-1.0], $
            proj_index: ci, $
            scan_name: title, $
            data_available: 0, $
            current_imagery: initimagery, $
            is_pressed: 0, $
            down_x: 0, $
            down_y: 0, $
            max_intensity: float(initintensity), $
            center_intensity: 255 , $
            min_intensity: 0.0, $
            min_data_thr: min_data, $
            max_data_thr: max_data, $
            draw_crosshairs: 1,  $
            show_line_profile: 0, $
            line_data_gui_win: 0, $
            line_prof_selected_window: '', $
            time_series_gui_win: 0, $
            axis_units: 'CM', $ ;;  ('CM' or 'MM')
            marker_list: marker_list_curr, $
            marker_list_wid: marker_list, $
            navigation_callback: '', $
            navigation_callback_continuous: 1B $
                    }, /no_copy)
    
    widget_control, txt_sl, set_value=strcompress(string((*state).slice), /remove_all)
    widget_control, txt_fq, set_value=strcompress(string((*state).freq), /remove_all)
    widget_control, txt_ph, set_value=strcompress(string((*state).phase), /remove_all)
    widget_control, txt_ar, set_value=strcompress(string(adim), /remove_all)

    widget_control, base_main, set_uvalue=state
    
    ortho_update_imagery, state, want_adim=0
    ortho_recompute_tx, state
    ortho_update_marker_list, state
    
    if (initimagery ne 4) then begin
        project.procpramarray[ci].ortho_tlb = base_main
    endif
    
    xmanager, 'ortho_view', base_main, cleanup='ortho_view_cleanup', no_block=keyword_set(block) ? 0 : 1

end
