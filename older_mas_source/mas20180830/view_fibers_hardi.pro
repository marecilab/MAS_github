;; $Id$
;; Copyright 2008 University of Florida. All Rights Reserved

function VFH_get_usable_colors

    forward_function FSC_Color
    
    colors = ['Khaki',      'Red',       'Green',        'Lime Green', 'Blue', 'Maroon', 'Cyan',  $
              'Magenta',    'Orange',    'Yellow',       'Orchid',     'Violet', $
              'Deep Pink',  'Tomato',    'Dark Green',   'Olive',$
              'Purple',     'Goldenrod', 'Saddle Brown', 'Coral',$
              'Brown',      'Navy',      'Light Cyan',   'Lavender',   'Green Yellow',$
              'Lawn Green', 'Honeydew']
    act_colors = intarr(3,n_elements(colors))
    for i = 0, n_elements(colors)-1 do begin
        act_colors[*,i] = FSC_Color(colors[i], /triple, /row)
    endfor

    return, act_colors
    
end

function VFH_plane_coeff, p1, p2, p3

    A = p1[1]*(p2[2]-p3[2]) + p2[1]*(p3[2]-p1[2]) + p3[1]*(p1[2] - p2[2])
    B = p1[2]*(p2[0]-p3[0]) + p2[2]*(p3[0]-p1[0]) + p3[2]*(p1[0] - p2[0])
    C = p1[0]*(p2[1]-p3[1]) + p2[0]*(p3[1]-p1[1]) + p3[0]*(p1[1] - p2[1])
    D = p1[0]*(p2[1]*p3[2] - p3[1]*p2[2]) + $
        p2[0]*(p3[1]*p1[2] - p1[1]*p3[2]) + $
        p3[0]*(p1[1]*p2[2] - p2[1]*p1[2])
    
    return, [A, B, C, -D]

end

pro VFH_update_clipping, state, R, P, S

    fib_par = (*state).model->getByName('FIB_parent_model')
    temp = where([R, P, S] ne 0, n_planes)
    if (n_planes eq 0) then begin
        fib_par->setProperty, clip_planes=-1
        xobjview, refresh=(*state).tlb
        return
    endif
    
    clip_planes = lonarr(4,n_planes)
    cp = 0L
    if (R ne 0) then begin
        clip_planes[*,cp++] = VFH_plane_coeff([R,1,1], [R,1,5], [R,5,1])
    endif
    if (P ne 0) then begin
        clip_planes[*,cp++] = VFH_plane_coeff([1,P,1], [1,P,5], [5,P,1])
    endif
    if (S ne 0) then begin
        clip_planes[*,cp++] = VFH_plane_coeff([1,1,S], [5,1,S], [1,5,S])
    endif
    
    fib_par->setProperty, clip_planes=clip_planes
    xobjview, refresh=(*state).tlb

end

;; Subroutine name: VFH_VFH_rebuild_fiber_cache, state
;; Created by: BT 20080421
;; Calling Information:
;;   state - the struct containing the state information. see below
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine: read in the raw fiber data array and cache
;;   information about roi containment, length, etc.
;;
;; Editing Information:
;;
pro VFH_rebuild_fiber_cache, state

    print, state
    
    progressbar = OBJ_NEW('progressbar', Color='red', /fast_loop, $
                          Title='Building fiber cache')
    progressbar->SetProperty, TEXT='Building data cache...'
    progressbar -> Start
    
    ; free the current cache
    if (ptr_valid((*state).fib_cache)) then begin
        old_fib_cache = *((*state).fib_cache)
        for j = 0L, n_elements(old_fib_cache)-1 do begin
            ptr_free, old_fib_cache[j]
            endfor
    endif
    
    fiber_data = (*state).fib_rawdata
    
    oroi_model = (*state).model->getByName('ROI_model')
    
    transform_vox2dat          = diag_matrix(replicate(1.0, 4))
    transform_vox2dat[0:2,0:2] = diag_matrix((*state).voxel_size)
    
    if ((*state).roi_count gt 0) then begin
        roi_work = bytarr([(*state).dimensions, (*state).roi_count])
        fib_arr  = 0
        for r = 0, (*state).roi_count-1 do begin
            roi_name = strcompress('ROI_'+string(r), /remove_all)
            roi = oroi_model->getByName(roi_name)
            roi->getProperty, uvalue=roi_pts
            ;;;if (obj_valid(roi_pts)) then oroi_work[r] = roi_pts
            for ro = 0L, n_elements(roi_pts)/3 - 1 do begin
                pt = round(roi_pts[*,ro])
                roi_work[pt[0],pt[1],pt[2], r] = 1B
            endfor
        endfor
    endif
    
    new_fib_cache = replicate(create_struct(name="VFH_FIBER_CACHE"), n_elements(*fiber_data))
    
    ;; check roi containment and length of fiber
    max_fiber_len = 0L
    ;;;;;;;;;;;;;;;
    terminal_pts = fltarr(3, 2*n_elements(*fiber_data)-1)
    tn=0L
    start_time = systime(1)
    for j = 0L, n_elements(*fiber_data)-1 do begin
        
        progressbar->Update, (float(j)/n_elements(*fiber_data)) * 100.0

        curr_fiber = (*fiber_data)[j]
        
        if not ptr_valid(curr_fiber) then continue

        ;; measured in number of coordinates.
        fsize = n_elements(*curr_fiber)/3
        if (fsize gt max_fiber_len) then max_fiber_len = fsize
        
        tol = 1
        
        terminal = (*curr_fiber)[*,[0, fsize-1]]
        terminal_pts[*,tn:tn+1] = terminal
        tn++
        
        tmp = vert_t3d(*curr_fiber, matrix=transform_vox2dat, /no_divide)

        point_distances = sqrt(total( ( (tmp[*,1:fsize-1] - tmp[*,0:fsize-2]) )^2, 1)) 
        voxlen = total(point_distances);; & print, voxlen ;; in fractions of voxel.
        if (voxlen gt 0) then begin
            mean_stepsize = mean(point_distances[where(point_distances ne 0)]); & print, mean_stepsize
        endif else begin
            mean_stepsize = 0.0
        endelse
        
        if ((*state).roi_count gt 0) then begin
        
            fib_passthru = intarr((*state).roi_count)
            fib_terminal = intarr((*state).roi_count)
            
            ;; TERMINAL INTERSECTIONS ONLY
            for st = 0L, n_elements(terminal)/3-1 do begin
                pt = (floor(terminal[*,st]) > [0,0,0]) < ((*state).dimensions-1)
                for ro = 0, (*state).roi_count-1 do begin
                    if (roi_work[pt[0],pt[1],pt[2],ro] ne 0) then begin
                        fib_terminal[ro] = (st eq 0) ? 1 : fsize-1 
                        fib_passthru[ro] = (st eq 0) ? 1 : fsize-1
                        break
                    endif
                endfor
            endfor
            
            ;; ALL INTERSECTIONS
            temp = lonarr(fsize)
            for ro = 0, (*state).roi_count-1 do begin
                temp[*] = ro
                test = roi_work[round((*curr_fiber)[0,*]), round((*curr_fiber)[1,*]), round((*curr_fiber)[2,*]), temp]
                st = where(test ne 0, n_st)
                if (n_st ne 0) then begin
                    fib_passthru[ro] = (st[0] eq 0) ? 1 : st[0]
                endif
            endfor

        endif else begin
            fib_terminal = 1
            fib_passthru = 1
        endelse

        ;; if fiber doesn't match any, then it belongs in the "default" ROI
        if (total(fib_passthru) eq 0) then begin
            fib_terminal[0] = 1
            fib_passthru[0] = 1
            ;print, j
        endif
        
        new_fib_cache[j].roi_hits_new = ptr_new(fib_passthru, /no_copy)
        new_fib_cache[j].roi_terminal = ptr_new(fib_terminal, /no_copy)
        new_fib_cache[j].roi_passthru = new_fib_cache[j].roi_hits_new
        new_fib_cache[j].mid_point = ( *( (*state).fib_midpoints ) )[j]
        new_fib_cache[j].branch_level = 0
        new_fib_cache[j].voxel_length = voxlen
        new_fib_cache[j].physical_length = voxlen
        new_fib_cache[j].voxel_stepsize = mean_stepsize
        new_fib_cache[j].length = fsize
                                               
        if (progressbar->CheckCancel()) then begin
            progressbar->Destroy
            ptr_free, new_fib_cache.roi_hits_new
            ptr_free, new_fib_cache.roi_terminal
            VFH_free_state, state
            return
        endif

    endfor
    
    terminal_pts = terminal_pts[*,0:tn-1]

    print, "Time: ", systime(1)-start_time
    (*state).fib_cache   = ptr_new(new_fib_cache, /no_copy)
    (*state).fib_nfibers = n_elements(*fiber_data)
    (*state).max_fiber_len = max_fiber_len+2
    (*state).fib_max_len_threshold = max_fiber_len+2
    progressbar->Destroy

end

pro VFH_fiber_cache__define

    struct = { VFH_FIBER_CACHE, $
               roi_hits_new: ptr_new(), $
               roi_terminal: ptr_new(), $
               roi_passthru: ptr_new(), $
               mid_point: 0L, $
               branch_level: 0L, $
               voxel_length: 0.0, $
               physical_length: 0.0, $
               voxel_stepsize: 0.0, $
               length: 0.0 }
end

pro VFH_do_network_export, state, tab_delimited=tab_delimited, rgraph=rgraph, rstatnet=rstatnet



end


pro VFH_export_network_data_event, event

    widget_control, event.top, get_uvalue=state
    
    case event.id of
    
        state.btn_cancel:  begin
            widget_control, event.top, /destroy
        end
        
        state.btn_export: begin
            type = widget_info(state.dl_type, /droplist_select)
            case type of
                0: VFH_do_network_export, state, /TAB_DELIMITED
                1: VFH_do_network_export, state, /RGRAPH
                2: VFH_do_network_export, state, /RSTATNET
                else:
            end
        end
        
        else: begin
        
        end
        
    endcase

end


pro VFH_export_network_data, state

    base = widget_base(title="Export Network Data", group_leader=(*state).gp_tlb, /modal, /column)
    b = widget_base(base, /row)
    lbl = widget_label(b, value="Choose the export type:")
    dl_type = widget_droplist(b, value=['Tab-delimited File', 'R-Graph import script', 'R-statnet import script'])
    b = widget_base(base, /row, /align_center)
    
    btn_cancel = widget_button(b, value="Cancel")
    btn_export = widget_button(b, value="Export...")
    
    state = { btn_export: btn_export, $
              btn_cancel: btn_cancel, $
              dl_type: dl_type, $
              main_state: state }

    widget_control, base, set_uvalue=state, /realize, /no_copy
    xmanager, "VFH_export_network_data", base
    
end


pro VFH_export_network_matrix, state

    common scan_data
    if (not ptr_valid((*state).graph_data)) then begin
        VFH_compute_network, state
    endif
    graph_data = (*state).graph_data
    
    len       = (*graph_data).tot_lengths
    weight    = (*graph_data).weights
    conn      = (*graph_data).conn
    degrees   = (*graph_data).degrees
    roi_names = (*graph_data).roi_names
    
    tmp = *conn
    non_zero = where(tmp ne 0, nnon_zero)
    if (nnon_zero ne 0) then begin
        tmp[non_zero] = 1
    endif
    
    dest = dialog_pickfile(path=project.current_path, /write, /overwrite_prompt, title="Save Adjacency Matrix to...", default_extension='txt')
    if (dest eq '') then return
    
    openw, lun, dest, /get_lun
    for j = 1L, (size(tmp, /dimensions))[0] - 1 do begin
        printf, lun, strjoin(string(reform(tmp[j-1,1:(size(tmp, /dimensions))[1]-1]), format='(I0)'), ",")
    endfor
    close, lun & free_lun, lun
    
    dest = dialog_pickfile(path=project.current_path, title="Save ROI Names to...")
    if (dest eq '') then return
    openw, lun, dest, /get_lun
    printf, lun, transpose(*roi_names)
    close, lun & free_lun, lun
    
end

pro VFH_compute_network, state

    graph_state = (*state).graph_state

    oroi_model = (*state).model->getByName('ROI_model')
    rois = oroi_model->get(ISA='IDLGRPOLYGON', /ALL)
    if (n_elements(rois) eq 1 and total(obj_valid(rois)) eq 0) then begin
        void = dialog_message(['There are no connections, or there are', $
                               'not enough ROIs to define connectivity.'], /error, /center)
        return
    endif
    
    n_rois = n_elements(rois)+1
    terminate = 0
    fib_cache = (*state).fib_cache
    
    conn   = lonarr(n_rois, n_rois)
    len    = fltarr(n_rois, n_rois)
    weight = fltarr(n_rois, n_rois)
    
    for j = 0L, n_elements(*fib_cache)-1 do begin

       if ((*graph_state).passthru_nodes ne 0) then begin
           fib_mask_conn = *((*fib_cache)[j].roi_passthru)
       endif else begin
           fib_mask_conn = *((*fib_cache)[j].roi_terminal)
       endelse

       conn_sorted = sort(fib_mask_conn)

       w = where(fib_mask_conn[1:n_elements(fib_mask_conn)-1] ne 0, count)
       ;if (count lt 2) then begin
       ;     continue
       ;endif
       ;
       ;; first sort the ROIs in the order in which they are
       ;; reached during the tracking
       tmp = fib_mask_conn[w+1]
       tmp_sorted= sort(tmp)
       tmp = tmp[tmp_sorted]
       ;; now arrange ROIs in sorted order. This causes the edges to be 
       ;; created correctly when a fiber "passes through" an ROI on its way
       ;; to another ROI.
       w = w[tmp_sorted]

       if (count gt 1) then begin
          for c = 1, count-1 do begin
              conn[w[c-1], w[c]]++
              conn[w[c], w[c-1]]++
              
              tmplen = (tmp[c] - tmp[c-1]) * (*fib_cache)[j].voxel_stepsize
              len[w[c-1], w[c]] += tmplen
              len[w[c], w[c-1]] += tmplen
              
              weight[w[c-1], w[c]] += 1.0/tmplen
              weight[w[c], w[c-1]] += 1.0/tmplen
          endfor
       endif
          
    endfor
    
    ;; len contains the total length of all of the segments of fibers between ROI A and ROI B
    ;; conn contains the number of fiber segments between ROI A and ROI B
    
    ;;openw, lun, '~/Desktop/G.R', /get_lun
    ;;
    ;;for r = 0, n_elements(rois)-1 do begin
    ;;    roi_name = string(r, format='(I03)')
    ;;    printf, lun, "g <- addNode('"+roi_name+"', g)"
    ;;endfor
    
    degrees = lonarr(n_elements(rois))
    roi_names = strarr(n_elements(rois))
    
    for r = 0, n_elements(rois)-1 do begin

       roi_names[r] = string(r+1, format='(I04)')

       dest = where(conn[*,r] gt 0, count)
       
       degrees[r] = count

       if (count eq 0) then continue

    endfor
    
    ;;close, lun
    if (ptr_valid((*state).graph_data)) then begin
        ptr_free, (*((*state).graph_data)).tot_lengths
        ptr_free, (*((*state).graph_data)).weights
        ptr_free, (*((*state).graph_data)).conn
        ptr_free, (*((*state).graph_data)).degrees
        ptr_free, (*((*state).graph_data)).roi_names
    endif
    
    (*state).graph_data = ptr_new({ tot_lengths: ptr_new(len, /no_copy), $
                                    weights: ptr_new(weight, /no_copy),  $
                                    conn: ptr_new(conn, /no_copy),       $
                                    degrees: ptr_new(degrees, /no_copy), $
                                    roi_names: ptr_new(roi_names, /no_copy) })
end

pro VFH_rebuild_fiber_model_conn, state, passthru_nodes=passthru_nodes
    
    common colors, b_curr, b_orig, g_curr, g_orig, r_curr, r_orig

    if (not ptr_valid((*state).graph_data)) then begin
        print, "WARNING: network not computed, using default values."
        VFH_compute_network, state
    endif
    
    omodel_fibs_par = (*state).model->getByName('FIB_parent_model')
    
    obj_destroy, omodel_fibs_par->getByName('FIB_child_model')
    
    omodel_fibs_chl = obj_new('idlgrmodel', name='FIB_child_model', lighting=0, depth_test_disable=2)
    
    oroi_model = (*state).model->getByName('ROI_model')
    rois = oroi_model->get(ISA='IDLGRPOLYGON', /ALL)
    if (n_elements(rois) eq 1 and total(obj_valid(rois)) eq 0) then begin
        void = dialog_message(['There are no connections, or there are', $
                               'not enough ROIs to define connectivity.'], /error, /center)
        return
    endif
    
    graph_state  = (*state).graph_state
    graph_params = (*state).graph_data
    n_rois = n_elements(rois)+1
    
    conn      = (*((*state).graph_data)).conn
    len       = (*((*state).graph_data)).tot_lengths
    weight    = (*((*state).graph_data)).weights
    roi_names = (*((*state).graph_data)).roi_names
    degrees   = (*((*state).graph_data)).degrees
    
    max_weight = alog(max(*conn))
    
    ;; len contains the total length of all of the segments of fibers between ROI A and ROI B
    ;; conn contains the number of fiber segments between ROI A and ROI B
    mean_len = (*len)/(*conn)
    n = 0
    col = [255,0,0]
;    
;    ;;openw, lun, '~/Desktop/G.R', /get_lun
;    ;;
;    ;;for r = 0, n_elements(rois)-1 do begin
;    ;;    roi_name = string(r, format='(I03)')
;    ;;    printf, lun, "g <- addNode('"+roi_name+"', g)"
;    ;;endfor
;    
;    
    loadct, 1
    for r = 0, n_elements(rois)-1 do begin

       rois[r]->getProperty, data=verts_src, name=roi_name
       
       ;; 0:r instead of * because of symmetry
       dest = where((*conn)[0:r,r] gt (((*graph_state).min_edge-1) > 0), count)
       roi_name = 'ROI: '+string(r+1, format='(I03)') + $
                  ', Node Degree: '+string((*degrees)[r], format='(I0)') + $
                  ', Node Strength: '+string(total((*weight)[*,r]), format='(G0)')
       rois[r]->setProperty, name=roi_name
       
       if (count eq 0) then continue

       for c = 0, count-1 do begin
          
          rois[dest[c]]->getProperty, data=verts_dest
          
          ;;src_pt = total(verts_src, 2)/float(n_elements(verts_src)/3.)
          ;;dst_pt = total(verts_dest, 2)/float(n_elements(verts_dest)/3.)
          
          src_pt = verts_src[*,0]
          dst_pt = verts_dest[*,0]
          
          ;;printf, lun, "g <- addEdge('"+string(r, format='(I03)')+"', '"+string(dest[c], format='(I03)')+"', g, "+string(weight[dest[c],r], format='(G0.5)')+")"
          
          name = 'ROI '+string(r+1, format='(I03)')+ ' <-> ROI '+string(dest[c]+1, format='(I03)')+ $
                 ', N_fibers='+string((*conn)[dest[c],r], format='(I0)')+ $
                 ', Avg_Len='+string(mean_len[dest[c],r], format='(G0)')+'mm'+ $
                 ', Weight='+string((*weight)[dest[c],r], format='(G0)')
          col_pick = fix(alog(float((*conn)[dest[c],r]))/max_weight * 255)
          col = [r_curr[col_pick], g_curr[col_pick], b_curr[col_pick]]
          o = obj_new('idlgrpolyline', $
                      name=name, $
                      [ [ src_pt ], [ dst_pt ] ], $
                      color=col, $
                      alpha=1, $
                      thick=1); (*conn)[dest[c],r])

          lines = (n++ eq 0) ?  o : [ lines, o ]
          
       endfor

    endfor
    loadct, 0
    if (n_elements(lines) ne 0) then begin
        omodel_fibs_chl->add, lines
        omodel_fibs_par->add, omodel_fibs_chl
    endif else begin
        void = dialog_message('No connections found.', /error, /center)
    endelse

    xobjview, refresh=(*state).tlb

end

;; Subroutine name: VFH_rebuild_fiber_model, state
;; Created by: BT 20080421
;; Calling Information:
;;   state - the struct containing the state information. see below
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine: rebuilds fiber graphics model. usually
;;  in response to a GUI event.
;;
;; Editing Information:
;;
pro VFH_rebuild_fiber_model, state

    compile_opt idl2
    print, "VFH_rebuild_fiber_model... rebuilding..."
    
    oroi_model = (*state).model->getByName('ROI_model')
    rois = oroi_model->get(ISA='IDLGRPOLYGON', /ALL)
    fstep = (*state).fib_fstep
    
    terminate = 0
    roi_mask_new   = *((*state).roi_mask_new)
    roi_type_mask  = *((*state).roi_type_mask)
    exclusion_mask = (*((*state).roi_type_mask))[2,*]

    for r = 1L, n_elements(roi_type_mask)/3-1 do begin
        rois[r-1]->setProperty, HIDE=(1-roi_mask_new[r])
    endfor
    
    if ((*state).display_as eq 1B) then begin
        VFH_rebuild_fiber_model_conn, state
        return
    endif
    
    fiber_data = (*state).fib_rawdata

    usable_colors = VFH_get_usable_colors()

    omodel_fibs_par = (*state).model->getByName('FIB_parent_model')
    
    obj_destroy, omodel_fibs_par->getByName('FIB_child_model')
    
    omodel_fibs_chl = obj_new('idlgrmodel', name='FIB_child_model', lighting=0)
    
    ;; set up fibers model
    fibs        = objarr((*state).fib_nfibers)
    n_discarded = 0L
    n_fibs      = 0L
    n_rois      = (*state).roi_count
        
    fib_mask = roi_type_mask[0:1,*]
    fib_cache = (*state).fib_cache
    
    for j = 0L, n_elements(*fiber_data)-1, fstep do begin

        curr_fiber = (*fiber_data)[j]

        ;; this should never happen
        if not ptr_valid(curr_fiber) then continue

        ;; check the length
        fsize = (*fib_cache)[j].length
        if fsize lt (*state).fib_min_len_threshold then continue
        if fsize gt (*state).fib_max_len_threshold then continue

        ;; check whether to draw this fiber based on which
        ;; ROIs the user wants to see
        fib_mask_passthru = *((*fib_cache)[j].roi_passthru) ne 0
        fib_mask_terminal = *((*fib_cache)[j].roi_terminal) ne 0
        
        fib_mask[0,*] = fib_mask_passthru
        fib_mask[1,*] = fib_mask_terminal
        type_mask = (fib_mask and roi_type_mask[0:1,*])

        ;;branch_lvl = 0;;(*((*(*state).fib_cache)[j])).branch_level

        if (total(type_mask) eq 0) then continue
        if (not array_equal(type_mask[0:1,*], roi_type_mask[0:1,*])) then continue
        if (total(exclusion_mask and (fib_mask_passthru or fib_mask_terminal)) ne 0) then continue

        pivot = (*((*state).fib_midpoints))[j]
        fiblen = n_elements(*curr_fiber)/3
        
        low = (pivot - pivot * float((*state).fib_length_pct)/100.0) > 0
        hi  = (pivot + (fiblen) * float((*state).fib_length_pct)/100.0) < fiblen

        fiblen_new = (hi-low) > 1

        if (fiblen_new lt fiblen) then begin
            fib_range = indgen(fiblen_new) + low
            fib_range = (fib_range > 0) < (fiblen - 1)
            fiber = (*curr_fiber)[*,[fib_range]]
        endif else begin
            fiber = *curr_fiber
        endelse
        
        ;; handle coloring of fibers according to user-selected method
        if ((*state).fib_color_style eq 1 and n_elements(fiber)/3 gt 3) then begin ;; directional coloring
            ;; compute difference to get tanget vectors for directional coloring
            d_fib = fiber[*,1:n_elements(fiber)/3-1] - fiber[*,0:n_elements(fiber)/3-2]
            d_fib /= congrid(transpose(sqrt(total(d_fib^2, 1))),3,n_elements(fiber)/3-2,interp=0)
            d_fib = abs(d_fib)*255.0
            fibs[n_fibs++] = obj_new('idlgrpolyline', fiber, $
                vert_colors=d_fib, $
                name='fiber_'+strcompress(string(j), /remove_all), $
                alpha_channel=(*state).fib_alpha,           $
                depth_test_disable=2,      $
                thick=(*state).fib_thick)
        endif else if (*state).fib_color_style eq 2 then begin ;; roi based coloring
            ok_roi = where(fib_mask_passthru ne 0, nok_roi)
            roi_color = reform(usable_colors[*, ((ok_roi+2)) mod (n_elements(usable_colors)/3)])
            roi_color = congrid(roi_color, 3, 5*n_elements(roi_color)/3, interp=0)
            if (nok_roi gt 1) then begin ;; multi-roi case
                fibs[n_fibs++] = obj_new('idlgrpolyline', fiber, $
                    name='fiber_'+strcompress(string(j), /remove_all), $
                    vert_colors = roi_color, $
                    alpha_channel=(*state).fib_alpha,           $
                    depth_test_disable=2,      $
                    thick=(*state).fib_thick)
            endif else begin ;; single roi case
                fibs[n_fibs++] = obj_new('idlgrpolyline', fiber, $
                    name='fiber_'+strcompress(string(j), /remove_all), $
                    color = roi_color, $
                    alpha_channel=(*state).fib_alpha,           $
                    depth_test_disable=2,      $
                    thick=(*state).fib_thick)
            endelse
        end else begin ;; monotone coloring
            fiber_color = [255,0,0]
            fibs[n_fibs++] = obj_new('idlgrpolyline', fiber, $
                name='fiber_'+strcompress(string(j), /remove_all), $
                color=fiber_color, $
                alpha_channel=(*state).fib_alpha,           $
                depth_test_disable=2,      $
                thick=(*state).fib_thick)
        endelse
    
    endfor

    ;; collapse array
    if (n_fibs gt 0) then begin
        ofibs = fibs[0:n_fibs-1]
        fibs  = 0
        omodel_fibs_chl->add, ofibs
    endif

    print, "Number of Fibers matching ROI criteria: ", string(n_elements(ofibs), format="(I0)")
    omodel_fibs_par->add, omodel_fibs_chl

end

pro VFH_rebuild_roi_model, state

    usable_colors = VFH_get_usable_colors()
    
    oroi_model = (*state).model->getByName('ROI_model')
    oroi_lbl_model = (*state).model->getByName('ROI_LBL_model')
    
    old_rois = oroi_model->get( /all, count=ct)
    if (ct ne 0) then begin
        for r = 0, ct-1 do begin
            obj_destroy, old_rois[r]
        endfor
        obj_destroy, oroi_lbl_model    
    endif

    roi_objects = (*state).roi_rawdata
    num_rois = n_elements(*roi_objects)
    
    (*state).roi_count = ++num_rois
    
    orois = objarr(num_rois)
    oroi_lbls = objarr(num_rois)
    
    oroi_lbl_model  = obj_new('idlgrmodel', depth_test_disable=2, $
        name='ROI_LBL_model', lighting=0, hide=0, transform=(*state).transform)
    ofnt_lbl = obj_new('idlgrfont', size=9)
    
    orois[0] = obj_new('idlgrroi', $
        name='ROI_0', $
        depth_test_disable=2, $
        [ [0,0,0] ], $
        style=0, $
        thick=4, $
        type=2, $
        color=[0,0,0], $
        alpha_channel=0.0)
        
    oroi_lbls[0] = obj_new('idlgrtext', $
        'B', $
        locations=[0,0,0], $
        /onglass,    $
        font=ofnt_lbl,$
        color=[255,255,255])
        
    vol = bytarr((*state).dimensions)
    
    for i = 1, num_rois-1 do begin
    
        tmp = (*roi_objects)[i-1]
        
        if (ptr_valid(tmp)) then begin
        
            c = reform((*tmp)[*,0])
            
            for f = 0L, n_elements(*tmp)/3 - 1 do begin
                pt = ((floor((*tmp)[*,f])) > 0) < ((*state).dimensions-1)
                vol[pt[0],pt[1],pt[2]] = 1B
            endfor
            
            roi_color = usable_colors[*,(i+2) mod (n_elements(usable_colors)/3)]
            shade_volume, vol, 0, v, p, /low
            nverts = n_elements(v)/3
            v += [ [replicate(0.5, nverts)], [replicate(0.5, nverts)], [replicate(0.5, nverts)] ] 

            ;;VFH_mask2oroigroup, data_ptr=vol, roigroup=roigroup
            mask_verts = array_indices(vol, where(vol ne 0))
            orois[i] = obj_new('idlgrpolygon', v, polygons=p, $
                color=roi_color, alpha=0.4, $
                depth_test_disable=2, uvalue=mask_verts, $
                name=strcompress('ROI_'+string(i),/remove_all))
            vol[*] = 0B
            
            oroi_lbls[i] = obj_new('idlgrtext', $
                string(i), $
                locations=c, $
                /onglass,    $
                font=ofnt_lbl,$
                color=[255,255,255])
        endif
        
    endfor
    
    oroi_model->add, orois
    oroi_lbl_model->add, oroi_lbls
    
    ((*state).model)->add, oroi_lbl_model
    
    ;; start with all rois visible
    (*state).roi_mask_new = ptr_new(intarr((*state).roi_count)+1)
    (*state).roi_type_mask = ptr_new(intarr(3,(*state).roi_count))
    
end

pro VFH_export_volume_mask, state, dimensions=dimensions, mask=mask, file_name=file_name

    common scan_data
    
    omodel_fibs_par = (*state).model->getByName('FIB_parent_model')
    omodel_fibs_chld = omodel_fibs_par->getByName('FIB_child_model')
    omodel_fibs_par->getProperty, transform=transform
    
    aa = omodel_fibs_chld->get(isa='IDLGRPOLYLINE', /all)
    if (n_elements(aa) eq 1) then begin ;; check to see if -1 was return 
        if (aa eq -1) then return       ;; meaning that no fibers were
    endif
    
    if (not keyword_set(dimensions)) then begin
        dimensions = (*state).dimensions
    endif
    
    mask = lonarr(dimensions)
    prev = [0L,0L,0L]
    
    for a = 0L, n_elements(aa)-1 do begin
        
        if (not obj_valid(aa[a])) then continue
        
        aa[a]->getProperty, data=fibs_native
        uu = vert_t3d(fibs_native, matrix=transform)
        
        for u = 0L, n_elements(uu)/3-1 do begin
            if (uu[0,u] gt dimensions[0]-1 or uu[1,u] gt dimensions[1]-1 or uu[2,u] gt dimensions[2]-1) then continue
            if (uu[0,u] lt 0   or uu[1,u] lt 0   or uu[2,u] lt 0) then continue
            if (round(uu[0,u]) eq prev[0] and round(uu[1,u]) eq prev[1] and round(uu[2,u]) eq prev[2]) then begin
                ;print, "Removing Dup.", uu[*,u], prev
                continue
            endif
            mask[uu[0,u],uu[1,u],uu[2,u]]++
            prev = round([ uu[0,u],uu[1,u],uu[2,u] ])
        endfor
    endfor

    ;;mask[where(mask ne 0)] = 1
    tmp = ptr_new(mask)
    
    if (keyword_set(file_name)) then begin
        if (file_test(file_dirname(file_name), /write)) then begin
            mas_export_nifti, data_ptr=tmp, file_name=file_name
        endif else begin
            void = dialog_message(['Unable to write -- check permissions.', $
                                   file_name], /center, /error)
        endelse
    endif else begin
        file_name = dialog_pickfile(path=project.current_path, /write, /overwrite_prompt, filter='*.nii,*.nii.gz')
        if (file_name ne '') then begin
            mas_export_nifti, data_ptr=tmp, file_name=file_name
        endif
    endelse
    
    ptr_free, tmp
    return
        
end

pro VFH_mask2oroigroup, data_ptr=data_ptr, roigroup=roigroup

    ;if (not ptr_valid(data_ptr)) then return
    
    data_sz = size(data_ptr, /dimensions)
    roigroup = obj_new('idlgrroigroup')
    
    for s = 0, data_sz[2]-1 do begin
    
        tmp = reform((data_ptr)[*,*,s])
        w = where(tmp ne 0, ct)
        if (ct eq 0) then continue
        
        ind3 = fltarr(3, n_elements(w))
        ind3[0:1, *] = array_indices([data_sz[0],data_sz[1]], w, /dimensions)
        ind3[2, *] = s
        
        roigroup->add, obj_new('idlgrroi', ind3)
        
    endfor
    
end

pro VFH_add_roi, state, roi_arr

    if (not ptr_valid(state)) then return
    old_roi_arr = (*state).roi_rawdata
    
    if (ptr_valid(old_roi_arr) && n_elements(*old_roi_arr) ne 0) then begin
        *( (*state).roi_rawdata ) = [ *( (*state).roi_rawdata ), roi_arr ]
    endif else begin
        (*state).roi_rawdata = ptr_new([ roi_arr])
    endelse
    
    VFH_rebuild_roi_model, state
    VFH_rebuild_fiber_cache, state
    VFH_rebuild_fiber_model, state
    VFH_gui_roilist, state
    
    xobjview, refresh=(*state).tlb
    
end

function VFH_create_nodes, state, radius=radius, wm_mask=wm_mask, terminal_pts=terminal_pts

    fiber_data = *((*state).fib_rawdata)
    
    mask = fltarr(size(*wm_mask, /dimensions))
    mask[where(*wm_mask ge 0.4)] = 1.0
    mask = smooth(mask,2) & mask[where(abs(mask) lt 0.5)] = 0 & mask[where(abs(mask) ge 0.5)] = 1

    pp = ptr_new(mask) & mas_display_ortho, data_ptr=pp, /block & ptr_free, pp
    
    ;;bwdist = morph_distance(mask, neighbor_sampling=3)
    ;;mas_display_ortho, data_ptr=ptr_new(bwdist), /block
    
    terminal_pts = fltarr(3, 2*n_elements(fiber_data)-1) 
    tn=0L

    for j = 0L, n_elements(fiber_data)-1 do begin
        
        curr_fiber = fiber_data[j]
        
        if not ptr_valid(curr_fiber) then continue

        ;; measured in number of coordinates.
        fsize = n_elements(*curr_fiber)/3

        terminal = (*curr_fiber)[*,[0, fsize-1]]
        terminal_pts[*,tn:tn+1] = terminal
        tn++
    endfor
    
    max_roi_radius = keyword_set(radius) ? float(radius) : 8.
    test_spheres = ptrarr(max_roi_radius)
    
    for roi_radius = 1, max_roi_radius do begin
        roi_diam = 2.*roi_radius
        volume = dblarr(roi_diam, roi_diam, roi_diam)
        for x=0, roi_diam-1 do begin
            for y=0, roi_diam-1 do begin
                for z=0, roi_diam-1 do begin
                    volume[x,y,z] = -total( (roi_diam/2.-.5-[x,y,z])^2 )
                endfor
            endfor
        endfor
        st=-(roi_diam/4.)^2*3.25
        sph_roi = array_indices(volume, where(volume gt st)) - roi_radius
        test_spheres[roi_radius-1] = ptr_new(sph_roi, /no_copy)
    endfor
    
    lbl = 1L
    roilabels = intarr((*state).dimensions)
    sph_roi = *test_spheres[n_elements(test_spheres)-1]
    erm = intarr((*state).dimensions)
    while (lbl ne 0) do begin

        ;; choose the first terminal point
        pt = terminal_pts[*,0]
        if (mask[pt[0],pt[1],pt[2]] ne 0) then begin
            ;; if this terminal pt lies within the WM mask, disregard it
            if (n_elements(terminal_pts) gt 3) then begin
                terminal_pts = terminal_pts[*,1:n_elements(terminal_pts)/3-1]
                continue
            endif else begin
                break
            endelse
        endif
        
        ;; find all other terminal points with roi_radius of this point
        close_points = where( sqrt( (terminal_pts[0,*] - replicate(pt[0],tn))^2 + $
                                      (terminal_pts[1,*] - replicate(pt[1],tn))^2 + $
                                      (terminal_pts[2,*] - replicate(pt[2],tn))^2) lt roi_radius, $
                              complement=keep, ncomplement=nkeep)
        match = terminal_pts[*,close_points]
;        roilabels[match[0,*], match[1,*], match[2,*]] = lbl
;        
;        roilabels[sph_roi[0,*]+replicate(pt[0], n_elements(sph_roi)/3), $
;                  sph_roi[1,*]+replicate(pt[1], n_elements(sph_roi)/3), $
;                  sph_roi[2,*]+replicate(pt[2], n_elements(sph_roi)/3) ] = lbl++

        ;; 1) turn "on" all voxel within roi_radius sphere of the terminal point
        erm[match[0,*], match[1,*], match[2,*]] = 1
        erm[sph_roi[0,*]+replicate(pt[0], n_elements(sph_roi)/3), $
                  sph_roi[1,*]+replicate(pt[1], n_elements(sph_roi)/3), $
                  sph_roi[2,*]+replicate(pt[2], n_elements(sph_roi)/3) ] = 1
        ;; 2) turn "off" all voxels enabled in step 1 that are in the WM mask
        erm *= 1-mask
        
        ;; 3) label the remaining regions
        ll = label_region(erm, /all_neighbors)
        
        ;; 4) find the label number that contains our original terminal point
        mv = ll[pt[0], pt[1], pt[2]]
        
        ;; 5) get the coordinates for the label contiaining our terminal point
        lblval = where(ll eq mv) & print, n_elements(lblval)
        if (n_elements(lblval) gt 4000) then begin
            ;; if this terminal pt lies within the WM mask, disregard it
            terminal_pts = terminal_pts[*,1:n_elements(terminal_pts)/3-1]
            continue
        endif

        lblind = array_indices(ll, lblval)
        
        ;; 
        roilabels[lblind[0, *], lblind[1,*], lblind[2,*]] = lbl++
        
        erm[*] = 0
        
        if (nkeep eq 0) then begin
            lbl = 0
        endif else begin
            terminal_pts = terminal_pts[*,keep]
            tn = nkeep
        endelse
        
        if (lbl gt 100000) then begin
            print, "over 1000000"
            break
        endif
        
    endwhile

    return, roilabels

end

pro VFH_free_state, state

    obj_destroy, (*state).model

    if (ptr_valid((*state).fib_rawdata)) then begin
        for p = 0L, n_elements(*(*state).fib_rawdata)-1 do begin
            ptr_free, (*((*state).fib_rawdata))[p]
        endfor
        ptr_free, (*state).fib_rawdata
    endif
    
    if (ptr_valid((*state).fib_cache)) then begin
        ptr_free, (*((*state).fib_cache)).roi_hits_new
        ptr_free, (*((*state).fib_cache)).roi_terminal
        ptr_free, (*state).fib_cache
    endif
    
    if (ptr_valid((*state).roi_rawdata)) then begin
        for p = 0L, n_elements(*(*state).roi_rawdata)-1 do begin
            ptr_free, (*(*state).roi_rawdata)[p]
        endfor
        ptr_free, (*state).roi_rawdata
    endif
    
    if (ptr_valid((*state).graph_data)) then begin
        graph_data = (*state).graph_data
        ptr_free, (*graph_data).tot_lengths
        ptr_free, (*graph_data).weights
        ptr_free, (*graph_data).conn
        ptr_free, (*graph_data).degrees
        ptr_free, (*graph_data).roi_names
        ptr_free, graph_data
     endif
     ptr_free, (*state).graph_state
        
    ;; TODO: Note that we need to free the imagery pointer IF it was loaded in by the user
    ;; but NOT if it is part of the active scan in the scan list.
    ptr_free, (*state).roi_mask_new
    ptr_free, (*state).roi_type_mask
    ptr_free, (*state).fib_midpoints
    ptr_free, state

end

;; Subroutine name: view_fibers
;;
;; Created by: BT 20080421
;; Calling Information:
;;   fiber_data    - output from make_fibers
;;   thick         - initial tract thickness (def: 1)
;;   alpha         - initial alpha for tracts(def: .75)
;;   show_axis     - whether or not to show  (def: 1)
;;   len_threshold - tract length threshold  (def: 20)
;;   fsteps        - fibers to skip (def 0)
;;   in_mas        - whether this is being called from within mas or
;;                   from the idl terminal
;;supplemental_omodel - a precreated omodel containing something that
;;                      should be displayed along with the fibers
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine: creates the visualization environment and
;; displays it
;;
;; Editing Information:
;;
pro view_fibers_hardi, fiber_data,     $ 
                   roi_data=roi_data,  $
                   thick=thick,        $
                   alpha=alpha,        $
                   show_axis=show_axis,$
                   len_threshold=len_threshold, $
                   fstep=fstep, $
                   in_mas=in_mas, $
                   working_dir=working_dir, $
                   project_dataset=project_dataset, $
                   supplemental_omodel=supplemental_omodel, $
                   handle=handle
          
    common scan_data, project
    forward_function get_dir_slash

    if (keyword_set(in_mas)) then begin
        in_mas = 1
        rotate_Direction = project.procPramArray[project.ci].rotate_Direction
        flip_Direction   = project.procPramArray[project.ci].flip_Direction
        interp_factors   = double([ project.procPramArray[project.ci].freq_Interp, $
                                     project.procPramArray[project.ci].phase_Interp, $
                                     project.procPramArray[project.ci].slice_Interp ] )
    endif else begin
        in_mas = 0
        rotate_Direction = 0
        flip_Direction   = 0
        interp_factors   = [ 1D,1D,1D ]
    endelse
    
     ;; in mas, we have the imagery available
    if (in_mas) then begin
        image_datapool = project.dataarray[project.ci].frac_ani
    
        if (keyword_set(project_dataset)) then begin
            (*fiber_data).transform = diag_matrix(1D/[interp_factors, 1D]) # (*fiber_data).transform
            (*fiber_data).dimensions = (size(*project.dataarray[project_dataset].state1, /dimensions))[0:2]
            image_datapool = project.dataarray[project_dataset].state1
        endif else begin
            project_dataset = project.ci
        endelse
    endif
    
    dimensions   = (*fiber_data).dimensions
    if (ptr_valid(image_datapool)) then begin
        dimensions   = (size(*image_datapool, /dimensions))[0:2]
    endif
        
    if ptr_valid(roi_data) then begin
    
        num_rois = n_elements((*roi_data).vertices)
        roi_objects = (*roi_data).vertices
    
    endif else num_rois = 0

    if not keyword_set(thick) then thick = 1.0
    if not keyword_set(alpha) then alpha = 0.75
    if not keyword_set(len_threshold) then len_threshold = 20
    if not keyword_set(fstep) then  fstep = 1
    if not keyword_set(working_dir) then working_dir = ''
    
    dep_fun = 4
    ;; Set up the various models
    oview = obj_new('idlexobjview', name='MAIN_view')
    omodel          = obj_new('idlgrmodel', name='MAIN_model', lighting=0)
    omodel_fibs_par = obj_new('idlgrmodel', name='FIB_parent_model', depth_test_disable=2, lighting=0, $
                               transform=(*fiber_data).transform, depth_test_function=dep_fun, clip_planes=-1)
    omodel_fibs_chl = obj_new('idlgrmodel', name='FIB_child_model', depth_test_disable=2, lighting=0)

    omodel_PS       = obj_new('idlgrmodel', name="PS_image", depth_test_disable=2, depth_test_function=dep_fun, lighting=0)
    omodel_FP       = obj_new('idlgrmodel', name="FP_image", depth_test_disable=2, depth_test_function=dep_fun, lighting=0)
    omodel_FS       = obj_new('idlgrmodel', name="FS_image", depth_test_disable=2, depth_test_function=dep_fun, lighting=0)
    omodel_im       = obj_new('idlgrmodel', name='IMG_model', lighting=0, depth_test_function=dep_fun)

    omodel_axes     = obj_new('idlgrmodel', name='AXES_model', lighting=0)

    oroi_model      = obj_new('idlgrmodel', name='ROI_model', depth_test_disable=2, depth_test_function=dep_fun, lighting=2, $
                               transform=(*fiber_data).transform)
    
    ;; order is important
    omodel->add, omodel_fibs_par, position=0
    omodel->add, oroi_model,      position=1
    omodel->add, omodel_im,       position=2    
    omodel->add, omodel_axes,     position=3
    
    omodel_fibs_par->add, omodel_fibs_chl
    omodel_im->add, omodel_PS
    omodel_im->add, omodel_FP
    omodel_im->add, omodel_FS
    
    omodel_FS->rotate, [1,0,0], 90
    omodel_PS->rotate, [1,0,0], 90
    omodel_PS->rotate, [0,0,1], 90
    
    if keyword_set(supplemental_omodel) and obj_valid(supplemental_omodel) then begin
        omodel->add, supplemental_omodel, position=4
    endif

    oview->add, omodel
    
    ;; set up axis
    if keyword_set(show_axis) then begin
    
        oaxfnt = obj_new('idlgrfont', size=8)
        oxaxis = obj_new('idlgraxis', 0, range=[0,dimensions[0]], color=[200,200,200])
        oyaxis = obj_new('idlgraxis', 1, range=[0,dimensions[1]], color=[200,200,200])
        ozaxis = obj_new('idlgraxis', 2, range=[0,dimensions[2]], color=[200,200,200])
        
        oxaxis->getProperty, ticktext=tmptxt
        tmptxt->setProperty, font=oaxfnt
        oyaxis->getProperty, ticktext=tmptxt
        tmptxt->setProperty, font=oaxfnt
        ozaxis->getProperty, ticktext=tmptxt
        tmptxt->setProperty, font=oaxfnt

        omodel_axes->add, oxaxis
        omodel_axes->add, oyaxis
        omodel_axes->add, ozaxis

    endif

;;  model             - the main idlgrmodel object. contains
;;                      everything else
;;  tlb               - the tlb for the main viewport. returned by
;;                      xobjview
;;  dimensions        - the dimensions (from mas) of the environment
;;  roi_on            - 1 ~ user want to see ROIs drawn, 0 ~ otherwise
;;  roi_count         - number of ROIs 
;;  roi_constrain     - 1 ~ user wants to see tracts running through
;;                      every visiable ROI, 0 ~ any ROI
;;  roi_identify      - 1 ~ user wants to see ROI numbers near ROIs, 0
;;                      ~ otherwise
;;  roi_mask          - n-bit vector (n = roi_count), where bit i=1 if
;;                      user wants to see ROI #i, zero otherwise. 
;;                      0<=i<=(roi_count-1)
;;  in_mas            - we're running from mas
;;  active_FP         - active freq-phase slice number
;;  active_FS         - active freq-slice slice number
;;  active_PS         - active phase-slice slice number
;;  image_alpha       - alpha percent for image display
;;  image_contrast    - contrast for image display
;;  image_datapool    - pointer to the data array that contains the
;;                      image data. don't free this pointer!
;;  fibers_front      - show fibers in front of imagery, otherwise
;;                      show imagery in front of fibers
;;  fib_rawdata       - raw fiber data (from make_fibers)
;;  fib_cache         - cache pointer
;;  fib_nfibers       - number of fibers
;;  fib_thick         - fiber thickness
;;  fib_alpha         - fiber alpha
;;  fib_fstep         - display increment
;;  fib_len_threshold - minimum length

    have_vd = where(tag_names(*fiber_data) eq 'VOXEL_DIMS', n_have_vd)
    if (n_have_vd ne 0) then begin
        voxel_dims = (*fiber_data).voxel_dims
    endif else begin
        voxel_dims = [1.0,1.0,1.0]
    endelse
    
    state = ptr_new({ working_directory:working_dir, $
                      view: oview, $
                      model:omodel, $
                      tlb:0L, $
                      roi_tlb: 0L, $
                      gui_tlb: 0L, $
                      gp_tlb: 0L, $
                      clipping_tlb: 0L, $
                      dimensions:dimensions, $
                      voxel_size: voxel_dims, $
                      transform: (*fiber_data).transform, $
                      roi_count:0, $
                      roi_identify:1, $
                      roi_mask_new: ptr_new(), $
                      roi_type_mask: ptr_new(), $
                      in_mas:in_mas, $
                      active_FP:-1, $
                      active_FS:-1, $
                      active_PS:-1, $
                      clipping_rps:[0,0,0L], $
                      min_intensity: long(0), $
                      max_intensity: long(0), $
                      image_alpha: .75, $
                      image_contrast: 0, $
                      image_datapool: ptr_valid(image_datapool) ? image_datapool : ptr_new(), $
                      fibers_front:0, $
                      display_as: 0B, $ ;; 0B = streamlines, 1B = connectivity graphy
                      roi_rawdata: ptr_new(roi_objects), $
                      fib_rawdata: ptr_new((*fiber_data).fibers), $
                      graph_state: ptr_new({passthru_nodes: 0, min_edge: 0, displaying_network: 0}), $
                      graph_data: ptr_new(), $
                      max_fiber_len: 0L, $
                      fib_midpoints: ptr_new((*fiber_data).seed_points), $
                      fib_cache: ptr_new(),$
                      fib_nfibers:(*fiber_data).nfibers, $
                      fib_branch_levels:(*fiber_data).branch_levels, $
                      fib_thick: thick, $
                      fib_alpha: alpha, $
                      fib_fstep: fstep, $
                      fib_length_pct: 100, $
                      fib_color_style:1B, $ ;; 0 = monotone, 1 = directional, 2 = ROI based
                      fib_max_len_threshold: 1000L, $
                      fib_min_len_threshold: len_threshold }, /no_copy)

    if (ptr_valid(image_datapool)) then begin
        (*state).max_intensity = max(*image_datapool, min=min_int)
        (*state).min_intensity = min_int
    endif
        
    ptr_free, fiber_data

    ;; set up ROI model
    if (num_rois gt -1) then begin
        VFH_rebuild_roi_model, state
    
        if (n_elements(ptr_seeds) ne 0) then begin
           ptr_free, ptr_seeds[0]
           ptr_free, ptr_seeds
        endif

    endif else begin
        (*state).roi_mask_new = ptr_new([1])
        (*state).roi_mask = 1
    endelse
    
    ;; build up the stuff
    VFH_rebuild_fiber_cache, state
    if (not ptr_valid(state)) then return ;; must have cancelled
    VFH_rebuild_fiber_model, state

    xobjview, oview, background=[0,0,0], tlb=tlb, xsize=800, ysize=600, renderer=1
    
    (*state).tlb = tlb

    view_fibers_hardi_GUI, state=state

    handle=state

end

