;; $Id$
;;
function VFH_load_image_data, state, slice_type, slice_number

    common scan_data, project
   
    ci = project.ci
    
    ;;pixszs = [ project.imndarray[ci].f_voxsz, project.imndarray[ci].p_voxsz, project.imndarray[ci].s_voxsz ]
    ;;pixszs /= min(pixszs)
    ;;print, pixszs
    
    ;;sized_matrix = (*state).dimensions * pixszs
    ;;print, sized_matrix
    if (*state).in_mas then begin
        
        if (slice_number lt 0) then begin
            case slice_type of
                'FP': data = ptr_new(bytarr((*state).dimensions[0], (*state).dimensions[1]))
                'FS': data = ptr_new(bytarr((*state).dimensions[0], (*state).dimensions[2]))
                'PS': data = ptr_new(bytarr((*state).dimensions[1], (*state).dimensions[2]))
                else: 
            endcase
            
        endif else begin
            
            case slice_type of
                'FP': begin
;;                   data = congrid(reform( (*( (*state).image_datapool ))[*,*,slice_number] ), sized_matrix[0], sized_matrix[1])
                   data = reform( (*( (*state).image_datapool ))[*,*,slice_number] )
                end
                'FS': begin
;;                   data = congrid(reform( (*( (*state).image_datapool ))[*,slice_number,*] ), sized_matrix[0], sized_matrix[2])
                   data = reform( (*( (*state).image_datapool ))[*,slice_number,*] )
                end
                'PS': begin
;;                   data = congrid(reform( (*( (*state).image_datapool ))[slice_number,*,*] ), sized_matrix[1], sized_matrix[2])
                   data = reform( (*( (*state).image_datapool ))[slice_number,*,*] )
                end
                else: 
            endcase
            data = ptr_new(data, /no_copy)
        endelse

    endif else begin

        image_file = strcompress(slice_type+'_'+string(slice_number)+'.tif', /remove_all)
        if not file_test(image_file) then return, -1
        tmp = read_tiff(image_file)
        data = ptr_new(temporary(reverse(tmp,2)), /no_copy)

    endelse
        
    return, data

end

pro VFH_gui_menu_transform_event, event

    common scan_data, project
    forward_function mas_read_transformation_matrix
    
    name = widget_info(event.id, /uname)
    widget_control, event.id, get_uvalue=state
    
    roi_lbl_model = (*state).model->getByName('ROI_LBL_model')
    omodel_fib_parent = (*state).model->getByName('FIB_parent_model')
    oroi_model = (*state).model->getByName('ROI_model')
    
    roi_lbl_model->getProperty, transform=orig_mat
    
    case name of
        'menu_repl_trans_mat' : begin
            do_invert = dialog_message('Invert matrix before using?', /question, /center)
            final_mat = mas_read_transformation_matrix(/ask)
            if (do_invert eq 'Yes') then final_mat = invert(final_mat)
        end
        
        'menu_apply_trans_mat' : begin
            do_invert = dialog_message('Invert matrix before using?', /question, /center)
            mat = mas_read_transformation_matrix(/ask)
            if (do_invert eq 'Yes') then mat = invert(mat)
            final_mat = orig_mat # mat
        end
        
        'menu_reset_trans_mat' : begin
            final_mat = diag_matrix(replicate(1D, 4))
        end
            
    endcase
    
    roi_lbl_model->setProperty, transform=final_mat
    omodel_fib_parent->setProperty, transform=final_mat
    oroi_model->setProperty, transform=final_mat
    xobjview, refresh=(*state).tlb
    
end

pro VFH_gui_graphparams_event, event

    name = widget_info(event.id, /uname)
    widget_control, event.top, get_uvalue=state
    graph_state = (*state).graph_state
    
    case name of 
    
        'btn_export_network': begin
           VFH_export_network_data, state
           
        end
        
        'btn_export_am': begin
            VFH_export_network_matrix, state
        end
                
        'btn_recompute_network': begin
        
            txt_min_fibers_edge = widget_info(event.top, find_by_uname='txt_min_fibers_edge')
            if (not widget_info(txt_min_fibers_edge, /valid_id)) then return
            widget_control, txt_min_fibers_edge, get_value=txtval
            min_edge = long(txtval)
            (*graph_state).min_edge = min_edge

            VFH_compute_network, state
            
            if ((*state).display_as eq 1) then begin
                VFH_rebuild_fiber_model_conn, state
            endif
            
        end
        
        'txt_min_fibers_edge': ;; this is a no-op
        
        'btn_passthru_nodes': begin
            bs = widget_info(event.id, /button_set)
            (*graph_state).passthru_nodes = bs
        end
        
        else: print, name
    
    endcase
    
    
    
end

pro VFH_gui_clipping_event, event

    name = widget_info(event.id, /uname)
    widget_control, event.top, get_uvalue=state
    
    curr_clip = (*state).clipping_rps
    
    case name of
    
        'SL_clip_F': begin
            val = event.value
            (*state).clipping_rps[0] = val
            VFH_update_clipping, state, val, curr_clip[1], curr_clip[2]
        end
    
        'SL_clip_P': begin
            val = event.value
            (*state).clipping_rps[1] = val
            VFH_update_clipping, state, curr_clip[0], val, curr_clip[2]
        end
        
        'SL_clip_S': begin
            val = event.value
            (*state).clipping_rps[2] = val
            VFH_update_clipping, state, curr_clip[0], curr_clip[1], val
        end
        
        'SL_clip_Z': begin
            val = float(event.value)
            (*state).view->setProperty, zclip=[val,-1]
            xobjview, refresh=(*state).tlb
        end
        
    endcase
    
end

pro VFH_gui_menu_event, event

    common scan_data, project
    forward_function mas_read_nifti
    
    name = widget_info(event.id, /uname)
    widget_control, event.id, get_uvalue=state
    
    case name of
        'menu_roi_display' : begin
            roi_list = (*state).roi_tlb
            if (widget_info(roi_list, /valid_id)) then return
            VFH_GUI_ROIlist, state
         end
        'menu_load_roi' : begin
            file = dialog_pickfile(path=(*state).working_directory, /read, filter=['*.nii;*.nii.gz'])
            if (file eq '') then return
            nif = mas_read_nifti(nifti_filename=file, read_status=rs)
            if (rs eq 0) then return
            mask = nif.voxel_data
            mas_hardi_tractography_addROI, state, mask
            ptr_free, mask
        end
        
        'menu_export_fibmask' : begin
            VFH_export_volume_mask, state, dimensions=(*state).dimensions
        end

        'menu_show_clipping': begin
            
            VFH_GUI_ClippingTools, state
        end
        
        else:
    endcase

end

pro VFH_gui_roilist_event, event

    name = widget_info(event.id, /uname)
    widget_control, event.id, get_uvalue=roi_num
    widget_control, event.top, get_uvalue=state
    
    case name of 
    
        'dl_roi_type': begin
            val = widget_info(event.id, /droplist_select)
            case val of
                0: begin
                    (*((*state).roi_type_mask))[0,roi_num] = 0
                    (*((*state).roi_type_mask))[1,roi_num] = 0
                    (*((*state).roi_type_mask))[2,roi_num] = 0
                end
                1: begin
                    (*((*state).roi_type_mask))[0,roi_num] = 1
                    (*((*state).roi_type_mask))[1,roi_num] = 0
                    (*((*state).roi_type_mask))[2,roi_num] = 0
                end
                2: begin
                    (*((*state).roi_type_mask))[0,roi_num] = 0
                    (*((*state).roi_type_mask))[1,roi_num] = 1
                    (*((*state).roi_type_mask))[2,roi_num] = 0              
                end
                3: begin
                    (*((*state).roi_type_mask))[0,roi_num] = 0
                    (*((*state).roi_type_mask))[1,roi_num] = 0
                    (*((*state).roi_type_mask))[2,roi_num] = 1               
                end
            endcase
            
            VFH_rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
            return
            
         end
        'btn_roi_ident':
                
        'btn_roi_show': begin
            widget_control, event.id, get_uvalue=roi
            (*((*state).roi_mask_new))[roi] = event.select
            VFH_rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
            return
        end
        'btn_fibers_show': begin
        end
        
        else:
        
    endcase

    VFH_rebuild_fiber_model, state
    xobjview, refresh=(*state).tlb
end


pro VFH_gui_event, Event

    common scan_data, project
    forward_function mas_read_nifti
    
    wTarget = (widget_info(Event.id,/NAME) eq 'TREE' ?  $
               widget_info(Event.id, /tree_root) : event.id)

    wWidget =  Event.top
    widget_control, wWidget, get_uvalue=state
    name = widget_info(wTarget, /uname)

    case name of
        
        'SL_FreqPhase': begin
            widget_control, wTarget, get_value=val
            pref = 'FP'
            model_name = pref+'_image'
            omodel = (*state).model->getByName('IMG_model/'+model_name)
            omodel->setProperty, hide=(val eq 0) ? 1 : 0

            (*state).active_FP = val-1

            VFH_update_slice_image, omodel, state, val;, data=data
        end
        
        'SL_PhaseSlice': begin
            widget_control, wTarget, get_value=val
            pref = 'PS'
            model_name = pref+'_image'
            omodel = (*state).model->getByName('IMG_model/'+model_name)
            omodel->setProperty, hide=(val eq 0) ? 1 : 0

            (*state).active_PS = val-1
            ;data = VFH_load_image_data(state, pref, val-1)

            VFH_update_slice_image, omodel, state, val;, data=data
        end
        
        'SL_FreqSlice': begin
            widget_control, wTarget, get_value=val
            pref = 'FS'
            model_name = pref+'_image'
            omodel = (*state).model->getByName('IMG_model/'+model_name)
            omodel->setProperty, hide=(val eq 0) ? 1 : 0

            (*state).active_FS = val-1
            ;data = VFH_load_image_data(state, pref, val-1)

            VFH_update_slice_image, omodel, state, -val;, data=data
        end
        
        'ROI_Display': begin
            (*state).roi_on = Event.select
            oroi_model = (*state).model->getByName('ROI_model')
            oroi_model->setProperty, HIDE=(1-(*state).roi_on)
            xobjview, refresh=(*state).tlb
        end
        
        'FIB_displayas': begin
            (*state).display_as = Event.select
            VFH_rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
         end
        
        'ROI_Constrain': begin
            (*state).roi_constrain = Event.select
            VFH_rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
        end
        
        'ROI_Identify': begin
            (*state).roi_identify = Event.select & print, event.select
            roi_lbl_model = (*state).model->getByName('ROI_LBL_model')
            if (obj_valid(roi_lbl_model)) then begin
                roi_lbl_model->setProperty, hide=(1-Event.select)
                xobjview, refresh=(*state).tlb
            endif
        end

        'ROI_HideAll': begin
            roi_model = (*state).model->GetByName('ROI_model')
            roi_model->setProperty, hide=event.select
            xobjview, refresh=(*state).tlb
            return
        end

        'SL_ImageAlpha': begin
            widget_control, wTarget, get_value=val
            widget_control, wTarget, get_uvalue=state
            (*state).image_alpha = float(val)/100.0
            VFH_update_slice_alpha, state
            xobjview, refresh=(*state).tlb
        end

        'SL_ImageContrast': begin
            widget_control, wTarget, get_value=val
            widget_control, wTarget, get_uvalue=state
            (*state).image_contrast = long(val)
            VFH_update_slice_contrast, state
            xobjview, refresh=(*state).tlb
        end

        'SL_MinFiberLenThr': begin
            widget_control, wTarget, get_value=val
            (*state).fib_min_len_threshold = val
            VFH_rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
        end

        'SL_MaxFiberLenThr': begin
            widget_control, wTarget, get_value=val
            (*state).fib_max_len_threshold = val
            VFH_rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
        end

        'SL_FiberLenPct': begin
            widget_control, wTarget, get_value=val
            (*state).fib_length_pct = val > 0
            VFH_rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
        end
        
        'SL_FiberAlpha': begin
            widget_control, wTarget, get_value=val
            (*state).fib_alpha = float(val)/100.0
            VFH_rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
        end

        'SL_FiberThickness': begin
            widget_control, wTarget, get_value=val
            (*state).fib_thick = long(val)
            VFH_rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
        end

        'SL_FiberStepInc': begin
            widget_control, wTarget, get_value=val
            (*state).fib_fstep = long(val)
            VFH_rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
        end

        'BTN_ColorStreamlineMono': begin
            widget_control, wTarget, get_value=val
            (*state).fib_color_style = 0B
            VFH_rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
        end

       'BTN_ColorStreamlineDir': begin
            widget_control, wTarget, get_value=val
            (*state).fib_color_style = 1B
            VFH_rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
        end

       'BTN_ColorStreamlineROI': begin
            widget_control, wTarget, get_value=val
            (*state).fib_color_style = 2B
            VFH_rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
        end

        'SL_ROI_Alpha': begin
            widget_control, wTarget, get_value=val
            oroi_model = (*state).model->getByName('ROI_model')
            rois = oroi_model->get(/ALL, ISA='IDLGRPOLYGON', count=ct)
            for i = 0,ct-1 do begin
                rois[i]->setProperty, alpha=float(val)/100.0
            endfor
            xobjview, refresh=(*state).tlb
        end

        'BTN_FibersFront': begin
            if (not(*state).fibers_front) then return

            omodel = (*state).model
            omodel->move, 0,2
            omodel->move, 0,1
            (*state).fibers_front = 1
            xobjview, refresh=(*state).tlb
        end

        'BTN_ImageryFront': begin
            if ((*state).fibers_front) then return

            omodel = (*state).model
            omodel->move, 0,2
            omodel->move, 0,1
            (*state).fibers_front = 1
            xobjview, refresh=(*state).tlb
        end

        'DL_ImageSource': begin
            widget_control, wTarget, get_value=val
            widget_control, wTarget, get_uvalue=state
            index = Event.index
            case index of
                0: (*state).image_datapool = project.dataarray[project.CI].frac_ani
                1: (*state).image_datapool = project.dataarray[project.CI].avg_dif
                2: (*state).image_datapool = project.dataArray[project.CI].adt
                3: begin ;; load external image
                    nif = mas_read_nifti(read_status=rs, default_dir=(*state).working_directory)
                    if (rs eq 0) then begin
                       return
                    endif else if (not ptr_valid(nif.voxel_data)) then begin
                       print, "NIFTI image did not load."
                       return
                    endif
                    
                    size_data = size(*nif.voxel_data, /dimensions)
                    if (n_elements(size_data) ne 3) then begin
                        void= dialog_message('Image data must have three dimensions.', /error, /center)
                        return
                    endif
                    (*state).image_datapool = nif.voxel_data
                    
                end
                4: begin
                    (*state).image_datapool = ptr_new()
                   end
                
            endcase
            
            if (not ptr_valid((*state).image_datapool)) then begin
                if (index ne 4) then begin
                    void = dialog_message('There is no imagery of that type available for this visualization.', $ 
                                          /error, /center)
                    widget_control, event.id, set_droplist_select=4
                endif
                VHF_gui_refresh_imagery_widgets, state
                return
            endif
            
            (*state).max_intensity = max(*((*state).image_datapool), min=min_int)
            (*state).min_intensity = min_int
            (*state).dimensions = (size(*((*state).image_datapool), /dimensions))[0:2]
                        
            VHF_gui_refresh_imagery_widgets, state
            VFH_update_slice_contrast, state
            xobjview, refresh=(*state).tlb
        end

        'BTN_UpdateFibers': begin
            VFH_rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
        end

        'BTN_network_params': begin
             VFH_GUI_GraphParams, state
         end
        
        else: begin
            print, name
        end
        
    endcase
end

pro VHF_gui_refresh_imagery_widgets, state

    SL_FP = widget_info((*state).gui_tlb, find_by_uname='SL_FreqPhase')
    SL_FS = widget_info((*state).gui_tlb, find_by_uname='SL_FreqSlice')
    SL_PS = widget_info((*state).gui_tlb, find_by_uname='SL_PhaseSlice')
    SL_AL = widget_info((*state).tlb, find_by_uname='SL_ImageAlpha')
    SL_CT = widget_info((*state).tlb, find_by_uname='SL_ImageContrast')
    
    if (ptr_valid((*state).image_datapool)) then begin
    
        widget_control, SL_AL, sensitive=1
        widget_control, SL_CT, sensitive=1
        widget_control, SL_FP, set_slider_max=(*state).dimensions[2]-1, set_value=min([(*state).active_FP, (*state).dimensions[2]-1]), sensitive=1
        widget_control, SL_FS, set_slider_max=(*state).dimensions[1]-1, set_value=min([(*state).active_FS, (*state).dimensions[1]-1]), sensitive=1
        widget_control, SL_PS, set_slider_max=(*state).dimensions[0]-1, set_value=min([(*state).active_PS, (*state).dimensions[0]-1]), sensitive=1
        widget_control, SL_FP, get_value=FP_val
        widget_control, SL_FS, get_value=FS_val
        widget_control, SL_PS, get_value=PS_val
        
        oaxes_model = (*state).model->getByName('AXES_model')
        oaxes = oaxes_model->get(ISA='idlgraxis', /all)
        oaxes[0]->setProperty, range=[0, (*state).dimensions[0]]
        oaxes[1]->setProperty, range=[0, (*state).dimensions[1]]
        oaxes[2]->setProperty, range=[0, (*state).dimensions[2]]
        
        (*state).active_FP = FP_val
        (*state).active_FS = FS_val
        (*state).active_PS = PS_val
        
    endif else begin
    
        widget_control, SL_FP, sensitive=0
        widget_control, SL_FS, sensitive=0
        widget_control, SL_PS, sensitive=0
        widget_control, SL_AL, sensitive=0
        widget_control, SL_CT, sensitive=0
        
    endelse
    
end

pro VFH_update_slice_contrast, state
   
    if (not ptr_valid((*state).image_datapool)) then return
    
    thr = 255.0 - float((*state).image_contrast)

    if ((*state).active_FP ne -1) then begin
        omodel = (*state).model->getByName('IMG_model/FP_image/img')
        if (obj_valid(omodel)) then begin
            img_data = VFH_load_image_data(state, 'FP', (*state).active_FP)
            m3ds_contrastify, img_data, thr, MIN=(*state).min_intensity, MAX=(*state).max_intensity
            omodel->setProperty, data=*img_data, /reset_data
            ptr_free, img_data
        endif
    endif

    if ((*state).active_FS ne -1) then begin
        omodel = (*state).model->getByName('IMG_model/FS_image/img')
        if (obj_valid(omodel)) then begin
            img_data = VFH_load_image_data(state, 'FS', (*state).active_FS)
            m3ds_contrastify, img_data, thr, MIN=(*state).min_intensity, MAX=(*state).max_intensity
            omodel->setProperty, data=*img_data, /reset_data
            ptr_free, img_data
        endif
    endif
    
    if ((*state).active_PS ne -1) then begin
        omodel = (*state).model->getByName('IMG_model/PS_image/img')
        if (obj_valid(omodel)) then begin
            img_data = VFH_load_image_data(state, 'PS', (*state).active_PS)
            m3ds_contrastify, img_data, thr, MIN=(*state).min_intensity, MAX=(*state).max_intensity
            omodel->setProperty, data=*img_data, /reset_data
            ptr_free, img_data
        endif
    endif
    
end

pro VFH_update_slice_alpha, state
    omodel = (*state).model->getByName('IMG_model/FP_image/img')
    if (obj_valid(omodel)) then begin
        omodel->setProperty, alpha_channel=(*state).image_alpha
    endif

    omodel = (*state).model->getByName('IMG_model/FS_image/img')
    if (obj_valid(omodel)) then begin
        omodel->setProperty, alpha_channel=(*state).image_alpha
    endif

    omodel = (*state).model->getByName('IMG_model/PS_image/img')
    if (obj_valid(omodel)) then begin
        omodel->setProperty, alpha_channel=(*state).image_alpha
    endif
end

pro VFH_update_slice_image, omodel, state, val, data=data

    ;;if not file_test(image_name) then return
    if (not ptr_valid((*state).image_datapool)) then return
    offset = 0.5
    if (not obj_valid(omodel)) then begin
        print, 'invalid!'
        return
    endif
    
    old_img = omodel->getByName('img')
    
    if (obj_valid(old_img)) then begin
        old_img->setProperty, location=[ offset, offset, val ], $
                                alpha_channel = (*state).image_alpha
    endif else begin

        oimg = obj_new('idlgrimage', transform_mode     = 1, $
                       alpha_channel      = (*state).image_alpha, $
                       interpolate        = 0, $;1,             $
                       location           = [ offset, offset, val ], $
                       depth_test_disable = 2, $
                       blend_function     = [3,4], $
                       name               = 'img')
        
        omodel->add, oimg
    endelse

    VFH_update_slice_contrast, state
    xobjview, refresh=(*state).tlb

    ;ptr_free, data
    
end

pro VFH_GUI_cleanup, tlb

    widget_control, tlb, get_uvalue=state
    if (ptr_valid(state)) then begin
        VFH_free_state, state
    endif
        
end


pro VFH_GUI_ClippingTools, state

    if (widget_info((*state).clipping_tlb, /valid_id)) then return
    (*state).view->getProperty, zclip=zclip
    
    base = widget_base(title="Fiber Clipping Adjustments", uname='base_clipping', /column, group_leader=(*state).tlb)
    SL_clip_F = widget_slider(base, title="Freq Clipping", uname='SL_clip_F', minimum=0, maximum=(*state).dimensions[0], value=(*state).clipping_rps[0], /drag, event_pro='VFH_gui_clipping_event')
    SL_clip_P = widget_slider(base, title="Phase Clipping", uname='SL_clip_P', minimum=0, maximum=(*state).dimensions[1], value=(*state).clipping_rps[1], /drag, event_pro='VFH_gui_clipping_event')
    SL_clip_S = widget_slider(base, title="Slice Clipping", uname='SL_clip_S', minimum=0, maximum=(*state).dimensions[2], value=(*state).clipping_rps[2], /drag, event_pro='VFH_gui_clipping_event')
    SL_clip_Z = CW_fslider(base, title="Z Clipping", uname='SL_clip_Z', minimum=-1, maximum=1, scroll=0.005, value=zclip[0], double=1, /drag);, event_pro='VFH_gui_clipping_event')
    
    widget_control, base, /realize
    widget_control, base, set_uvalue=state
    xmanager, 'VFH_gui_clipping', base, /no_block
    (*state).clipping_tlb = base
    
end

pro VFH_GUI_GraphParams, state

    gp_tlb = (*state).gp_tlb
    if (widget_info(gp_tlb, /valid_id) ne 0) then begin
        widget_control, gp_tlb, /destroy
    endif
    
    base = widget_base(title="Network Graph Settings", uname="base_graph_params", group_leader=(*state).tlb, /row)
    
    base_left = widget_base(base, /frame, /column)
    
    b = widget_base(base_left, /row)
    lbl = widget_label(b, value='Minimum Fibers for Edge:')
    txt = widget_text(b, value='1', uname='txt_min_fibers_edge', xsize=5, /editable)

    b = widget_base(base_left, /nonexclusive, /row)
    btn = widget_button(b, value='Consider passthru ROIs as nodes', uname='btn_passthru_nodes')

    b = widget_base(base_left, /row)
    btn = widget_button(b, value="Export Adjacency Matrix", uname='btn_export_am')
    btn = widget_button(b, value="Export Network Data", uname='btn_export_network', sensitive=0)
    
    b = widget_base(base_left, /row)
    btn = widget_button(b, value="Recompute Network", uname='btn_recompute_network')
    
    widget_control, base, set_uvalue=state
    widget_control, base, /realize
    
    (*state).gp_tlb = base
    
    xmanager, 'VFH_gui_graphparams', base, /no_block
    
end


pro VFH_GUI_ROIlist, state ;;view_fibers_GUI_ROI_manager, state

    ;existing = where(names eq 'view_fibers_ROIlist', count)
    roi_tlb = (*state).roi_tlb
    if (widget_info(roi_tlb, /valid_id) ne 0) then begin
        widget_control, roi_tlb, /destroy
    endif

        roi_mask_new = *((*state).roi_mask_new)
        roi_type_mask = *((*state).roi_type_mask)
        ;exclusion_mask = (*((*state).roi_type_mask))[2,*]

    base = widget_base(title="ROI List", uname='base_roi_list', /column, group_leader=(*state).tlb, scr_xsize=400, scr_ysize=400, /scroll)

    for j=0,(*state).roi_count - 1 do begin
    
        if (j eq 0) then begin
            roi_title = 'BRN'
        endif else begin
            roi_title = string(j, format='(I03)')
        endelse

        ROI_base = widget_base(base, uname='BASE_ROI_'+roi_title, /row, /frame, /base_align_center)

        b = widget_base(ROI_base, /row, /base_align_center)
        lbl = widget_label(b, value=roi_title)
        
        b = widget_base(ROI_base, /row)
        lbl = widget_label(b, value='ROI Type:')
        com = widget_droplist(b, value=['Disabled', 'Passthrough', 'Terminal', 'Exclusion'], uname='dl_roi_type', uvalue=j)
        active = where(roi_type_mask[*,j] ne 0, n_active)
        if (n_active ne 0) then begin
            widget_control, com, set_droplist_select=active+1
        endif
        
        b = widget_base(ROI_base, /nonexclusive, /row)
        btn = widget_button(b, value="Show ROI", uname='btn_roi_show', uvalue=j, sensitive=(j eq 0) ? 0 : 1)
        widget_control, btn, set_button=roi_mask_new[j]
        
    endfor

    widget_control, base, /realize
    widget_control, base, set_uvalue=state
    xmanager, 'VFH_gui_roilist', base, /no_block
    (*state).roi_tlb = base

end

pro view_fibers_hardi_GUI, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_, state=state

    matrix = (*state).dimensions
    have_imagery = ptr_valid((*state).image_datapool)
    
    base_main = Widget_Base( GROUP_LEADER=(*state).tlb, $
                        UNAME='BASE', $
                        XOFFSET=817, $
                        YOFFSET=0, $
                        TITLE='ImgSelect', $
                        SPACE=3, $
                        XPAD=3, $
                        YPAD=3, $
                        COLUMN=2,$
                        FRAME=0, /base_align_center, /align_center)

    base_tlb = widget_base((*state).tlb, $
                           uname='base_tlb',$
                           row=1,$
                           /align_center, $
                           /base_align_center)
    
    if (*state).in_mas then begin
        v = ['Fract. Anisotropy', 'Avg. Diffusity', 'S0 Image', 'Other (.nii file)...', 'None']
        DL_ImageSource = Widget_Droplist(base_tlb,$
                                         value=v, $
                                         uvalue=state,$
                                         title='Imagery Source:',$
                                         uname='DL_ImageSource',$
                                         event_pro='VFH_gui_event')
        if (not have_imagery) then begin
            widget_control, DL_ImageSource, set_droplist_Select=4
        endif
    endif
    
    BASE_image = widget_base(base_main, /frame, /column, /base_align_center, /align_center)
    lbl = widget_label(BASE_image, value="Image Controls")
    
    SL_FreqPhase = Widget_Slider(BASE_image, $
                                 UNAME='SL_FreqPhase', $
                                 XOFFSET=3, $
                                 minimum=0, $
                                 maximum=matrix[2], $
                                 YOFFSET=3, $ 
                                 TITLE='Freq/Phase', $
                                 scroll=1, $ 
                                 sensitive=have_imagery, $
                                 /drag)
    
    SL_PhaseSlice = Widget_Slider(BASE_image, $
                                  UNAME='SL_PhaseSlice',$
                                  minimum=0, $
                                  maximum=matrix[0],  $
                                  XOFFSET=3, $
                                  YOFFSET=50,$
                                  TITLE='Phase/Slice', $
                                  scroll=1, $
                                  sensitive=have_imagery, $
                                  /drag)
    
    SL_FreqSlice = Widget_Slider(BASE_image, $
                                 UNAME='SL_FreqSlice', $
                                 XOFFSET=3, $
                                 minimum=0, $
                                 maximum=matrix[1], $
                                 YOFFSET=97, $
                                 TITLE='Freq/Slice', $
                                 scroll=1, $
                                 sensitive=have_imagery, $
                                 /drag)

    base2a = widget_base(BASE_image, uname='base2a', /exclusive)

    BTN_FibersFront  = widget_button(base2a, uname='BTN_ImageryFront', value='Imagery in Front')
    BTN_ImageryFront = widget_button(base2a, uname='BTN_FibersFront', value='Fibers in Front')
    widget_control, BTN_FibersFront, /SET_BUTTON

    SL_ImageAlpha = Widget_Slider(base_tlb,$;(*state).tlb,$;BASE, $
                                  uvalue=state,$
                                  UNAME='SL_ImageAlpha', $
                                  XOFFSET=3, $
                                  minimum=0, $
                                  maximum=100,$
                                  VALUE=fix((*state).image_alpha*100), $
                                  YOFFSET=97, $
                                  TITLE='Image Opacity', $
                                  scroll=10, $
                                  sensitive=have_imagery, $
                                  /drag, $
                                  event_pro = 'VFH_gui_event')

    SL_ImageContrast = Widget_Slider(base_tlb,$;(*state).tlb,$;BASE, $
                                     uvalue=state,$
                                     UNAME='SL_ImageContrast', $
                                     XOFFSET=3, $
                                     minimum=0, $
                                     maximum=255, $
                                     value=(*state).image_contrast, $
                                     YOFFSET=97, $
                                     TITLE='Image Contrast', $
                                     scroll=10, $
                                     sensitive=have_imagery, $
                                     /drag, $
                                     event_pro = 'VFH_gui_event')

    BASE_roi = widget_base(base_main, /column, /frame, /base_align_center)
    lbl = widget_label(base_roi, value="Regions of Interest")
    SL_ROI_Alpha = widget_slider(BASE_roi, title="ROI Opacity", $
                                   minimum=0, maximum=100, value=45, uname='SL_ROI_Alpha', /drag)

    base2b = widget_base(BASE_roi, uname='base2b', /nonexclusive)
    ROI_Identify = Widget_Button(BASE2b, $
                                 UNAME='ROI_Identify', $
                                 XOFFSET=3, $
                                 YOFFSET=97, $
                                 VALUE="Identify ROIs")
    ROI_HideAll = Widget_Button(BASE2b, $
                                 UNAME='ROI_HideAll', $
                                 XOFFSET=3, $
                                 YOFFSET=97, $
                                 VALUE="Hide All ROIs")
    FIB_displayas = Widget_Button(BASE2b, $
                                 UNAME='FIB_displayas', $
                                 XOFFSET=3, $
                                 YOFFSET=97, $
                                 VALUE="Display as Network Graph")
    
    widget_control, FIB_displayas, set_button=0
    widget_control, ROI_Identify, /set_button
    
    bb = widget_button(BASE_roi, uname='BTN_network_params', value="Network Graph Settings...")

    base2 = Widget_Base( base_main,$
                         UNAME='BASE2', $
                         SPACE=3, $
                         XPAD=3, $
                         YPAD=3, $
                         COLUMN=1,$
                         FRAME=1)
    lbl = widget_label(BASE2, value="Fiber Controls")
    
    SL_MinFiberLenThr = Widget_Slider(BASE2, $
                                  UNAME='SL_MinFiberLenThr', $
                                  XOFFSET=3, $
                                  minimum=0, $
                                  maximum=(*state).max_fiber_len,$
                                  VALUE=min([(*state).fib_min_len_threshold,(*state).max_fiber_len]), $
                                  YOFFSET=97, $
                                  TITLE='Min. Length Threshold', $
                                  scroll=20, $
                                  /drag)

    SL_MaxFiberLenThr = Widget_Slider(BASE2, $
                                  UNAME='SL_MaxFiberLenThr', $
                                  XOFFSET=3, $
                                  minimum=0, $
                                  maximum=(*state).max_fiber_len,$
                                  VALUE=(*state).fib_max_len_threshold, $
                                  YOFFSET=97, $
                                  TITLE='Max. Length Threshold', $
                                  scroll=20, $
                                  /drag)

    SL_FiberLenPct = Widget_Slider(BASE2, $
                                  UNAME='SL_FiberLenPct', $
                                  XOFFSET=3, $
                                  minimum=1, $
                                  maximum=100,$
                                  VALUE=(*state).fib_length_pct, $
                                  YOFFSET=97, $
                                  TITLE='Length (%)', $
                                  scroll=10, $
                                  /drag)

    SL_FiberStepInc = Widget_Slider(BASE2, $
                                  UNAME='SL_FiberStepInc', $
                                  XOFFSET=3, $
                                  minimum=1, $
                                  maximum=20,$
                                  VALUE=(*state).fib_fstep, $
                                  YOFFSET=97, $
                                  TITLE='Thinning Factor', $
                                  scroll=20)

    SL_FiberAlpha = Widget_Slider(BASE2, $
                                  UNAME='SL_FiberAlpha', $
                                  XOFFSET=3, $
                                  minimum=0, $
                                  maximum=100,$
                                  VALUE=fix((*state).fib_alpha*100), $
                                  YOFFSET=97, $
                                  TITLE='Opacity', $
                                  scroll=10, $
                                  /drag)

    SL_FiberThickness = Widget_Slider(BASE2, $
                                     UNAME='SL_FiberThickness', $
                                     XOFFSET=3, $
                                     minimum=1, $
                                     maximum=10, $
                                     value=(*state).fib_thick, $
                                     YOFFSET=97, $
                                     TITLE='Streamline Thickness', $
                                     scroll=10, $
                                     /drag)
    lbl = widget_label(base2, value='Fiber Coloring:')                                
    fib_col_base = widget_base(base2, /column, /exclusive)
    BTN_ColorStreamlineMono = widget_button(fib_col_base, value='Monotone', uname='BTN_ColorStreamlineMono')
    BTN_ColorStreamlineDir = widget_button(fib_col_base, value='Directional', uname='BTN_ColorStreamlineDir')
    BTN_ColorStreamlineROI = widget_button(fib_col_base, value='ROI Based', uname='BTN_ColorStreamlineROI')
    widget_control, BTN_ColorStreamlineDir, set_button=(*state).fib_color_style eq 1 ? 1 : 0
    widget_control, BTN_ColorStreamlineMono, set_button=(*state).fib_color_style eq 0 ? 1 : 0
    widget_control, BTN_ColorStreamlineROI, set_button=(*state).fib_color_style eq 2 ? 1 : 0
    ;;;
    mbar = widget_info((*state).tlb, find_by_uname='xobjview:mbar')
    if (widget_info(mbar, /valid_id)) then begin
        track_menu = widget_button(mbar, value='Tractography', uname='menu_tractography', event_pro='VFH_gui_menu_event', uvalue=state)
        menu_roi_display = widget_button(track_menu, value='Show ROI List', uname='menu_roi_display', event_pro='VFH_gui_menu_event', uvalue=state)
        menu_load_roi = widget_button(track_menu, value='Load ROI Mask', uname='menu_load_roi', event_pro='VFH_gui_menu_event', uvalue=state)
        menu_transform = widget_button(track_menu, value="Transformation Matrix", /menu)
        menu_repl_trans_mat = widget_button(menu_transform, value='Replace Current Matrix', uname='menu_repl_trans_mat', event_pro='VFH_gui_menu_transform_event', uvalue=state)
        menu_apply_trans_mat = widget_button(menu_transform, value='Apply Matrix to Current', uname='menu_apply_trans_mat', event_pro='VFH_gui_menu_transform_event', uvalue=state)
        menu_reset_trans_mat = widget_button(menu_transform, value='Reset Matrix to Identity', uname='menu_reset_trans_mat', event_pro='VFH_gui_menu_transform_event', uvalue=state)
        menu_export_fibmask = widget_button(track_menu, value='Export Fiber Density Mask', uname='menu_export_fibmask', event_pro='VFH_gui_menu_event', uvalue=state)
        menu_show_clipping = widget_button(track_menu, value='Show Fiber Clipping Sliders', uname='menu_show_clipping', event_pro='VFH_gui_menu_event', uvalue=state)
    endif
    ;;;
    
    widget_control, /REALIZE, base_main
    widget_control, base_main, set_uvalue=state
    xmanager, 'VFH_gui', base_main, cleanup='VFH_GUI_cleanup', /no_block
    (*state).gui_tlb = base_main
    VFH_GUI_ROIlist, state

end

; 
; Empty stub procedure used for autoloading.
; 
;pro xx_view_fibers_GUI, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_, state=state
;  BASE_fibers, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_, state=state
;end
