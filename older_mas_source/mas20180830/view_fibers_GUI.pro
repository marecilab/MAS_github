;; $Id$
;;
function load_image_data_fibers, state, slice_type, slice_number

    common scan_data, project
   
    ci = project.ci
    pixszs = [ project.imndarray[ci].f_voxsz, project.imndarray[ci].p_voxsz, project.imndarray[ci].s_voxsz ]
    pixszs /= min(pixszs)
    ;;print, pixszs
    
    sized_matrix = (*state).matrix * pixszs
    ;;print, sized_matrix
    if (*state).in_mas then begin

        if (slice_number lt 0) then begin
            case slice_type of
                'FP': data = ptr_new(bytarr((*state).matrix[0], (*state).matrix[1]))
                'FS': data = ptr_new(bytarr((*state).matrix[0], (*state).matrix[2]))
                'PS': data = ptr_new(bytarr((*state).matrix[1], (*state).matrix[2]))
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

pro BASE_roi_event_fibers, Event

    wTarget = (widget_info(Event.id,/NAME) eq 'TREE' ?  $
               widget_info(Event.id, /tree_root) : event.id)

    wWidget =  Event.top
    name = widget_info(wTarget, /uname)

    widget_control, wWidget, get_uvalue=state
    widget_control, wTarget, get_uvalue=roi

    mask = (*state).roi_mask
    vec = 2^roi

    if (Event.select) then begin
        mask = mask or vec
    endif else begin
        mask = mask and not vec
    endelse

    (*state).roi_mask = mask
    rebuild_fiber_model, state
    xobjview, refresh=(*state).tlb

end

pro BASE_fibers_event, Event

    common scan_data, project

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
            ;data = load_image_data_fibers(state, pref, val-1)

            update_slice_image_fibers, omodel, state, val;, data=data
        end
        
        'SL_PhaseSlice': begin
            widget_control, wTarget, get_value=val
            pref = 'PS'
            model_name = pref+'_image'
            omodel = (*state).model->getByName('IMG_model/'+model_name)
            omodel->setProperty, hide=(val eq 0) ? 1 : 0

            (*state).active_PS = val-1
            ;data = load_image_data_fibers(state, pref, val-1)

            update_slice_image_fibers, omodel, state, val;, data=data
        end
        
        'SL_FreqSlice': begin
            widget_control, wTarget, get_value=val
            pref = 'FS'
            model_name = pref+'_image'
            omodel = (*state).model->getByName('IMG_model/'+model_name)
            omodel->setProperty, hide=(val eq 0) ? 1 : 0

            (*state).active_FS = val-1
            ;data = load_image_data_fibers(state, pref, val-1)

            update_slice_image_fibers, omodel, state, -val;, data=data
        end
        
        'ROI_Display': begin
            (*state).roi_on = Event.select
            oroi_model = (*state).model->getByName('ROI_model')
            oroi_model->setProperty, HIDE=(1-(*state).roi_on)
            xobjview, refresh=(*state).tlb
        end

        'ROI_Constrain': begin
            (*state).roi_constrain = Event.select
            rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
        end
        
        'ROI_Identify': begin
            (*state).roi_identify = Event.select
            roi_lbl_model = (*state).model->getByName('ROI_model/ROI_LBL_model')
            roi_lbl_model->setProperty, hide=(1-Event.select)
            xobjview, refresh=(*state).tlb
        end

        'SL_ImageAlpha': begin
            widget_control, wTarget, get_value=val
            widget_control, wTarget, get_uvalue=state
            (*state).image_alpha = float(val)/100.0
            update_slice_alpha_fibers, state
            xobjview, refresh=(*state).tlb
        end

        'SL_ImageContrast': begin
            widget_control, wTarget, get_value=val
            widget_control, wTarget, get_uvalue=state
            (*state).image_contrast = long(val)
            update_slice_contrast_fibers, state
            xobjview, refresh=(*state).tlb
        end

        'SL_FiberLenThr': begin
            widget_control, wTarget, get_value=val
            (*state).fib_len_threshold = val
            rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
        end

        'SL_FiberMaxLen': begin
            widget_control, wTarget, get_value=val
            (*state).fib_max_len = val > 0
            rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
        end
        
        'SL_FiberAlpha': begin
            widget_control, wTarget, get_value=val
            (*state).fib_alpha = float(val)/100.0
            rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
        end

        'SL_FiberThickness': begin
            widget_control, wTarget, get_value=val
            (*state).fib_thick = long(val)
            rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
        end

        'SL_FiberStepInc': begin
            widget_control, wTarget, get_value=val
            (*state).fib_fstep = long(val)
            rebuild_fiber_model, state
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
            if ((*state).fibers_front) then return

            omodel = (*state).model
            omodel->move, 0,2
            omodel->move, 0,1
            (*state).fibers_front = 1
            xobjview, refresh=(*state).tlb
        end

        'BTN_ImageryFront': begin
            if (not (*state).fibers_front) then return

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
                0: (*state).image_datapool = project.dataarray[project.ci].frac_ani
                1: (*state).image_datapool = project.dataarray[project.ci].avg_dif
                2: (*state).image_datapool = project.dataArray[project.ci].adt
            endcase
            (*state).max_intensity = max(*((*state).image_datapool), min=min_int)
            (*state).min_intensity = min_int
            update_slice_contrast_fibers, state
            xobjview, refresh=(*state).tlb
        end

        'BTN_UpdateFibers': begin
            rebuild_fiber_model, state
            xobjview, refresh=(*state).tlb
        end

        else: begin
            ;;print, name
        end
        
    endcase
end

pro update_slice_contrast_fibers, state
   
    thr = 255.0 - float((*state).image_contrast)

    if ((*state).active_FP ne -1) then begin
        omodel = (*state).model->getByName('IMG_model/FP_image/img')
        if (obj_valid(omodel)) then begin
            img_data = load_image_data_fibers(state, 'FP', (*state).active_FP)
            m3ds_contrastify, img_data, thr, MIN=(*state).min_intensity, MAX=(*state).max_intensity
            omodel->setProperty, data=*img_data, /reset_data
            ptr_free, img_data
        endif
    endif

    if ((*state).active_FS ne -1) then begin
        omodel = (*state).model->getByName('IMG_model/FS_image/img')
        if (obj_valid(omodel)) then begin
            img_data = load_image_data_fibers(state, 'FS', (*state).active_FS)
            m3ds_contrastify, img_data, thr, MIN=(*state).min_intensity, MAX=(*state).max_intensity
            omodel->setProperty, data=*img_data, /reset_data
            ptr_free, img_data
        endif
    endif
    
    if ((*state).active_PS ne -1) then begin
        omodel = (*state).model->getByName('IMG_model/PS_image/img')
        if (obj_valid(omodel)) then begin
            img_data = load_image_data_fibers(state, 'PS', (*state).active_PS)
            m3ds_contrastify, img_data, thr, MIN=(*state).min_intensity, MAX=(*state).max_intensity
            omodel->setProperty, data=*img_data, /reset_data
            ptr_free, img_data
        endif
    endif
    
end

pro update_slice_alpha_fibers, state
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

pro update_slice_image_fibers, omodel, state, val, data=data

    ;;if not file_test(image_name) then return

    if (not obj_valid(omodel)) then begin
        print, 'invalid!'
        return
    endif
    
    old_img = omodel->getByName('img')
    
    if (obj_valid(old_img)) then begin
        old_img->setProperty, location=[ -0.5,-0.5, val ], $
                                alpha_channel = (*state).image_alpha
    endif else begin

        oimg = obj_new('idlgrimage', transform_mode     = 1,  $
                       alpha_channel      = (*state).image_alpha, $
                       interpolate        = 0, $;1,             $
                       location           = [ -0.5, -0.5, val ], $
                       depth_test_disable = 2,             $
                       blend_function     = [3,4],         $
                       name               = 'img')
        
        omodel->add, oimg
    endelse

    update_slice_contrast_fibers, state
    xobjview, refresh=(*state).tlb

    ;ptr_free, data
    
end

pro view_fibers_cleanup, tlb

    widget_control, tlb, get_uvalue=state

    obj_destroy, (*state).model
;    if keyword_set(show_axis) then obj_destroy, oaxfnt
    ptr_free, (*state).fib_rawdata
    ptr_free, (*state).fib_cache
    ptr_free, state

end

pro BASE_fibers, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_, state=state

    matrix = (*state).matrix

    base = Widget_Base( GROUP_LEADER=(*state).tlb, $
                        UNAME='BASE', $
                        XOFFSET=817, $
                        YOFFSET=0, $
                        TITLE='ImgSelect', $
                        SPACE=3, $
                        XPAD=3, $
                        YPAD=3, $
                        COLUMN=1,$
                        FRAME=1)

    base_tlb = widget_base((*state).tlb, $
                           uname='base_tlb',$
                           row=1,$
                           /align_center, $
                           /base_align_center)
    
    if (*state).in_mas then begin
        v = ['Fract. Anisotropy', 'Avg. Diffusity', 'S0 Image']
        DL_ImageSource = Widget_Droplist(base_tlb,$
                                         value=v, $
                                         uvalue=state,$
                                         title='Image Source:',$
                                         uname='DL_ImageSource',$
                                         event_pro='BASE_fibers_event')
    endif

    SL_FreqPhase = Widget_Slider(BASE, $
                                 UNAME='SL_FreqPhase', $
                                 XOFFSET=3, $
                                 minimum=0, $
                                 maximum=matrix[2], $
                                 YOFFSET=3, $ 
                                 TITLE='Freq/Phase', $
                                 scroll=1, $ 
                                 /drag)
    
    
    SL_PhaseSlice = Widget_Slider(BASE, $
                                  UNAME='SL_PhaseSlice',$
                                  minimum=0, $
                                  maximum=matrix[0],  $
                                  XOFFSET=3, $
                                  YOFFSET=50,$
                                  TITLE='Phase/Slice', $
                                  scroll=1, $
                                  /drag)
    
    
    SL_FreqSlice = Widget_Slider(BASE, $
                                 UNAME='SL_FreqSlice', $
                                 XOFFSET=3, $
                                 minimum=0, $
                                 maximum=matrix[1], $
                                 YOFFSET=97, $
                                 TITLE='Freq/Slice', $
                                 scroll=1, $
                                 /drag)

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
                                  /drag, $
                                  event_pro = 'BASE_fibers_event')

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
                                     /drag, $
                                     event_pro = 'BASE_fibers_event')

    base2 = Widget_Base( base,$
                         UNAME='BASE2', $
                         SPACE=3, $
                         XPAD=3, $
                         YPAD=3, $
                         COLUMN=1,$
                         FRAME=0)

    SL_FiberLenThr = Widget_Slider(BASE2, $
                                  UNAME='SL_FiberLenThr', $
                                  XOFFSET=3, $
                                  minimum=0, $
                                  maximum=5000,$
                                  VALUE=(*state).fib_len_threshold, $
                                  YOFFSET=97, $
                                  TITLE='Fiber Length', $
                                  scroll=20, $
                                  /drag)

    SL_FiberStepInc = Widget_Slider(BASE2, $
                                  UNAME='SL_FiberStepInc', $
                                  XOFFSET=3, $
                                  minimum=1, $
                                  maximum=20,$
                                  VALUE=(*state).fib_fstep, $
                                  YOFFSET=97, $
                                  TITLE='Fiber Step Increment', $
                                  scroll=20, $
                                  /drag)

    SL_FiberMaxLen = Widget_Slider(BASE2, $
                                  UNAME='SL_FiberMaxLen', $
                                  XOFFSET=3, $
                                  minimum=1, $
                                  maximum=100,$
                                  VALUE=(*state).fib_max_len, $
                                  YOFFSET=97, $
                                  TITLE='Max Fiber Length (%)', $
                                  scroll=10, $
                                  /drag)

    SL_FiberAlpha = Widget_Slider(BASE2, $
                                  UNAME='SL_FiberAlpha', $
                                  XOFFSET=3, $
                                  minimum=0, $
                                  maximum=100,$
                                  VALUE=fix((*state).fib_alpha*100), $
                                  YOFFSET=97, $
                                  TITLE='Fiber Opacity', $
                                  scroll=10, $
                                  /drag)

    SL_FiberThickness = Widget_Slider(BASE2, $
                                     UNAME='SL_FiberThickness', $
                                     XOFFSET=3, $
                                     minimum=1, $
                                     maximum=10, $
                                     value=(*state).fib_thick, $
                                     YOFFSET=97, $
                                     TITLE='Fiber Thickness', $
                                     scroll=10, $
                                     /drag)

    SL_ROI_Alpha = widget_slider(Base2, title="ROI Opacity", $
                                   minimum=0, maximum=100, value=45, uname='SL_ROI_Alpha', /drag)

    base2a = widget_base(base2, uname='base2a', /exclusive)

    BTN_FibersFront  = widget_button(base2a, uname='BTN_FibersFront', value='Fibers in Front')
    BTN_ImageryFront = widget_button(base2a, uname='BTN_ImageryFront', value='Imagery in Front')
    widget_control, BTN_ImageryFront, /SET_BUTTON

    base2b = widget_base(base2, uname='base2b', /nonexclusive)
    ROI_Identify = Widget_Button(BASE2b, $
                                 UNAME='ROI_Identify', $
                                 XOFFSET=3, $
                                 YOFFSET=97, $
                                 VALUE="Identify ROIs")
    widget_control, ROI_Identify, /set_button

    ROI_Display = Widget_Button(BASE2b, $
                                UNAME='ROI_Display', $
                                XOFFSET=3, $
                                YOFFSET=97, $
                                VALUE="Show ROI Area")
    widget_control, ROI_Display, /set_button

    ROI_Constrain = Widget_Button(BASE2b, $
                                UNAME='ROI_Constrain', $
                                XOFFSET=3, $
                                YOFFSET=97, $
                                VALUE="Constrain ROIs")
                                   
    tmp = widget_label(base, value='ROIs to display:')

    BASE_roi = widget_base(base, uname='BASE_roi_0', /nonexclusive, /row)
    for j=0,(*state).roi_count - 1 do begin
        if (j eq 0) then begin
            roi_title = 'BR'
        endif else begin
            roi_title = strcompress(string(j), /remove_all)
        endelse

        if (j mod 4 eq 0) then begin
            BASE_roi = widget_base(base, uname='BASE_roi_'+roi_title, /nonexclusive, /row)
        endif
        btn = widget_button(BASE_roi, value=roi_title, uvalue=j, uname=roi_title, event_pro='BASE_roi_event_fibers')
        widget_control, btn, /set_button
    endfor

    widget_control, /REALIZE, base
    widget_control, base, set_uvalue=state
    xmanager, 'BASE_fibers', base, cleanup='view_fibers_cleanup', /no_block

end
; 
; Empty stub procedure used for autoloading.
; 
pro view_fibers_GUI, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_, state=state
  BASE_fibers, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_, state=state
end
