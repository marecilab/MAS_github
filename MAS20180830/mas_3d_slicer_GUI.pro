;; $Id$
;; Copyright 2008 University of Florida. All Rights Reserved

;; Subroutine name: m3ds_load_image_data
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;    state:   The state pointer associated with this mas_3d_slicer session
;;
;;    slice_type: Two letter indicator of the axes currently being loaded. One
;;                of FP (Frequency/Phase), PS (Phase/Slice), FS (Frequency/Slice).
;;
;;    slice_number: The index of the slice along the currently requested axis (slice_type)
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Returns a two-dimensions image of the data at the current slice index and axis.
;;
;; Editing Information:

function m3ds_load_image_data, state, slice_type, slice_number

    common scan_data, project
   
    if (slice_number lt 0) then begin
        case slice_type of
            'FP': data = ptr_new(bytarr((*state).matrix[0], (*state).matrix[1]))
            'FS': data = ptr_new(bytarr((*state).matrix[0], (*state).matrix[2]))
            'PS': data = ptr_new(bytarr((*state).matrix[1], (*state).matrix[2]))
            else: 
        endcase

    endif else begin

        case slice_type of
            'FP': data = reform( (*( (*state).image_datapool ))[*,*,slice_number] )
            'FS': data = reform( (*( (*state).image_datapool ))[*,slice_number,*] )
            'PS': data = reform( (*( (*state).image_datapool ))[slice_number,*,*] )
            else: 
        endcase
        data = ptr_new(data, /no_copy)
    endelse

    return, data

end

;; Subroutine name: mas_3d_slicer_GUI_cleanup
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Cleans up data structures associated with this session
;;
;; Editing Information:

pro mas_3d_slicer_GUI_cleanup, tlb

  print, 'mas_3d_slicer_GUI_cleanup: ',tlb

   widget_control, tlb, get_uvalue=state
   if (ptr_valid(state)) then begin
      if (ptr_valid( (*state).marker_list )) then begin
         nmarkers = n_elements( *((*state).marker_list) )
         for m = 0, nmarkers-1 do begin
            obj_destroy, (*((*state).marker_list))[m]
         endfor
         ptr_free, (*state).marker_list
      endif

      obj_destroy, (*state).model

      ptr_free, state
      
   endif

end

;; Subroutine name: mas_3d_slicer_GUI_event
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Event handler for the widgets on the control window
;;
;; Editing Information:

pro mas_3d_slicer_GUI_event, Event

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
            (*state).active_FP = val
            omodel->setProperty, hide=(val lt 0) ? 1 : 0
            data = m3ds_load_image_data(state, pref, val)
            val = val*(*state).s_scale + (*state).s_scale/2.0 ;; voxel size adjustment
            m3ds_update_slice_image, omodel, state, val, data=data;, scale=(*state).s_scale
        end
        
        'SL_PhaseSlice': begin
            widget_control, wTarget, get_value=val
            pref = 'PS'
            model_name = pref+'_image'
            omodel = (*state).model->getByName('IMG_model/'+model_name)
            (*state).active_PS = val
            omodel->setProperty, hide=(val lt 0) ? 1 : 0
            data = m3ds_load_image_data(state, pref, val)
            val = val*(*state).f_scale + (*state).f_scale/2.0
            m3ds_update_slice_image, omodel, state, val, data=data;, scale=(*state).f_scale
        end
        
        'SL_FreqSlice': begin
            widget_control, wTarget, get_value=val
            pref = 'FS'
            model_name = pref+'_image'
            omodel = (*state).model->getByName('IMG_model/'+model_name)

            (*state).active_FS = val
            omodel->setProperty, hide=(val lt 0) ? 1 : 0
            data = m3ds_load_image_data(state, pref, val)
            val = val*(*state).p_scale + (*state).p_scale/2.0
            m3ds_update_slice_image, omodel, state, -val, data=data;, scale=(*state).p_scale
        end
        
        'SL_ImageAlpha': begin
            widget_control, wTarget, get_value=val
            widget_control, wTarget, get_uvalue=state
            (*state).image_alpha = float(val)/100.0
            m3ds_update_slice_alpha, state
            xobjview, refresh=(*state).tlb
        end

        'SL_ImageContrast': begin
            widget_control, wTarget, get_value=val
            widget_control, wTarget, get_uvalue=state
            (*state).image_contrast = long(val)
            m3ds_update_slice_contrast, state
            xobjview, refresh=(*state).tlb
        end

;; TODO: Implement this
;;        'DL_ImageSource': begin
;;            widget_control, wTarget, get_value=val
;;            widget_control, wTarget, get_uvalue=state
;;            index = Event.index
;;            case index of
;;                0: (*state).image_datapool = project.dataarray[project.ci].frac_ani
;;                1: (*state).image_datapool = project.dataarray[project.ci].avg_dif
;;                2: (*state).image_datapool = project.dataArray[project.CI].adt
;;            endcase
;;            m3ds_update_slice_contrast, state
;;            xobjview, refresh=(*state).tlb
;;        end

        'btn_add_marker': begin
            m3ds_add_marker, state
         end

        'btn_rem_marker': begin
           sel = widget_info((*state).marker_list_wid, /LIST_SELECT)
           m3ds_remove_marker, state, sel
        end

        'marker_list': ;; null case to prevent messages printed from
                       ;; the default "else:" case

        else: begin
            print, 'mas_3d_slicer_GUI_event: recv;d event from unknown wid: '+name
        end
        
    endcase
end

;; Subroutine name: m3ds_contrastify
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;     pdata:  pointer to data to be contrast-adjusted
;;
;;     thr: thresholding value
;;
;;     MIN: user requested minimum value. default is the data's min
;;
;;     MAX: user requested maximum value. default is the data's max
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Adjust the contrast of an image to display.
;;
;; Editing Information:

pro m3ds_contrastify, pdata, thr, MIN=min, MAX=max

    if (not keyword_set(min)) then begin
        min = min(*pdata)
    endif
    
    if (not keyword_set(max)) then begin
        max = max(*pdata)
    endif

    ;; first bytscl the image like normal, then adjust to requested contrast
    *pdata = temporary(bytscl($
                               bytscl(*pdata, min=min, max=max), $
                             min=0, max=thr))

end

;; Subroutine name: m3ds_update_slice_contrast
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;   state:  the state pointer for this session
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Adjust the contrast of an image to display.
;;
;; Editing Information:

pro m3ds_update_slice_contrast, state
   
    thr = 255 - float((*state).image_contrast)

    ;; only update if there is something to display
    if ((*state).active_FP ne -1) then begin
        omodel = (*state).model->getByName('IMG_model/FP_image/img')
        if (obj_valid(omodel)) then begin
            img_data = m3ds_load_image_data(state, 'FP', (*state).active_FP)
            m3ds_contrastify, img_data, thr, MIN=(*state).min_intensity, MAX=(*state).max_intensity
            omodel->setProperty, data=*img_data, /reset_data
            ptr_free, img_data
        endif
    endif

    if ((*state).active_FS ne -1) then begin
        omodel = (*state).model->getByName('IMG_model/FS_image/img')
        if (obj_valid(omodel)) then begin
            img_data = m3ds_load_image_data(state, 'FS', (*state).active_FS)
            m3ds_contrastify, img_data, thr, MIN=(*state).min_intensity, MAX=(*state).max_intensity
            omodel->setProperty, data=*img_data, /reset_data
            ptr_free, img_data
        endif
    endif
    
    if ((*state).active_PS ne -1) then begin
        omodel = (*state).model->getByName('IMG_model/PS_image/img')
        if (obj_valid(omodel)) then begin
            img_data = m3ds_load_image_data(state, 'PS', (*state).active_PS)
            m3ds_contrastify, img_data, thr, MIN=(*state).min_intensity, MAX=(*state).max_intensity
            omodel->setProperty, data=*img_data, /reset_data
            ptr_free, img_data
        endif
    endif
    
end

;; Subroutine name: m3ds_update_slice_alpha
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;   state:  the state pointer for this session
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Adjust the alpha channel of an image to display.
;;
;; Editing Information:

pro m3ds_update_slice_alpha, state
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

;; Subroutine name: m3ds_remove_marker
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;   state:  the state pointer for this session
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Revmoes a marker from the marker list and the 3D space
;;
;; Editing Information:

pro m3ds_remove_marker, state, marker_num

    if (not ptr_valid(state)) then return
  
    markers = (*state).marker_list
    
    if (not ptr_valid(markers)) then return
    
    n_markers = n_elements(*markers)
    
    if (marker_num lt 0 or marker_num gt n_markers) then return
    
    obj_destroy, (*markers)[marker_num]
    
    valid_markers = where(obj_valid(*markers), n_valid)
    if (n_valid gt 0) then begin
       temp = (*markers)[valid_markers]
       (*state).marker_list = ptr_new(temp)
    endif else begin
       ptr_free, (*state).marker_list
       (*state).marker_list = ptr_new()
    endelse
    
    m3ds_update_marker_list, state

    ptr_free, markers  
    xobjview, refresh=(*state).tlb
    
end

;; Subroutine name: m3ds_add_marker
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;   state:  the state pointer for this session
;;
;;   coord:  the [r,p,s] location in space to add the marker
;;
;;   style:  the style of marker. default is a 3D "+". This keyword is
;;            passed directly to obj_new('idlgrsymbol'). See that
;;            entry in IDL's help for more information.
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Adds a marker in the 3D space located at <coord>. If coord is not
;;   specified, then it is added at the current intersection point of
;;   the three displayed image planes.
;;
;; Editing Information:

pro m3ds_add_marker, state, coord, style=style

    if (not ptr_valid(state)) then return
    
    marker_model = (*state).model->getByName('MARKER_model')
    
    if (not obj_valid(marker_model)) then return
    
    style = n_elements(style) ne 0 ? fix(style) : 1
    if (style gt 9 or style lt 0) then style = 1
    
    sym = obj_new('idlgrsymbol', style, color=[255,0,0])
    
    if (n_elements(coord) eq 0) then begin
       coord = [(*state).active_PS, (*state).active_FS, (*state).active_FP]
    endif

    name = 'Marker: ['+strcompress(string(coord[0])+','+ $
                                   string(coord[1])+','+ $
                                   string(coord[2]), /remove_all)+']'
    ;; voxel XYZ * voxel dimension + 1/2 voxel dim to center it in the destination voxel.
    ;; assumes voxel location is bottom-left corner of 3D cube.
    coord = [coord[0] * (*state).f_scale + (*state).f_scale/2., $
             coord[1] * (*state).p_scale + (*state).p_scale/2., $
             coord[2] * (*state).s_scale + (*state).s_scale/2.]

    omarker = obj_new('idlgrpolyline', coord, color=[255,0,0], symbol=sym, name=name)
    
    marker_model->add, omarker
    
    if (not ptr_valid((*state).marker_list)) then begin
       markers = [ omarker ]
       (*state).marker_list = ptr_new(markers)
    endif else begin
       *((*state).marker_list) = [ *((*state).marker_list), omarker ]
    endelse
    
    m3ds_update_marker_list, state

    xobjview, refresh=(*state).tlb

end

;; Subroutine name: m3ds_update_marker_list
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;   state:  the state pointer for this session
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Ensures that the list of markers displayed to the user is
;;   up-to-date with the marker list backing store.
;;
;; Editing Information:

pro m3ds_update_marker_list, state

    marker_list = (*state).marker_list

    if (not ptr_valid(marker_list)) then begin
       widget_control, (*state).marker_list_wid, set_value=''
       return
    endif

    n_markers = n_elements(*marker_list)

    marker_str = strarr(n_markers)

    for i = 0, n_markers-1 do begin

       (*marker_list)[i]->getProperty, name=marker_name

       ;; MARKER: [xx,yy,zz]
       ;; 012345678
       marker_name = strmid(marker_name, 8)
       marker_str[i] = marker_name

    endfor
   
    widget_control, (*state).marker_list_wid, set_value=marker_str

end

;; Subroutine name: m3ds_update_slice_image
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;  omodel:  the IDLgrModel containing this image for an axis
;;
;;   state:  the state pointer for this session
;;
;;  val: the location along the axis to display the image
;;
;;  scale: the amount of offset to add based on the voxel dimension scaling
;;
;;  data: the image data to display at the location/scale in the view area
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Adjust the alpha channel of an image to display.
;;
;; Editing Information:

pro m3ds_update_slice_image, omodel, state, val, data=data, scale=scale

    ;;if not file_test(image_name) then return

    if (not obj_valid(omodel)) then begin
        print, 'invalid!'
        return
    endif
    
    old_img = omodel->getByName('img')
    
    if (obj_valid(old_img)) then begin
        obj_destroy, old_img
    endif
 
    if (keyword_set(data)) then begin
        if (n_elements(*data) eq 0) then return
    endif else begin
        return
    endelse

    m3ds_contrastify, data, 255.0 - float((*state).image_contrast), $
      min=(*state).min_intensity, max=(*state).max_intensity

    oimg = obj_new('idlgrimage', *data,                $
                   transform_mode      = 1,             $
                   alpha_channel       = (*state).image_alpha, $
                   interpolate         = 0,             $
                   location            = [ 0, 0, val*(keyword_set(scale) ? scale : 1) ], $
                   depth_test_disable  = 2,             $
                   blend_function      = [3,4],         $
                   name                = 'img')
    
    ptr_free, data
    omodel->add, oimg
    xobjview, refresh=(*state).tlb
    
end

;; Subroutine name: mas_3d_slicer_GUI
;; Created by: BT, 2008-06
;; Calling Information:
;;
;;   state:  the state pointer for this session
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Set up and display the GUI control window for this session
;;
;; Editing Information:

pro mas_3d_slicer_GUI, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_, state=state

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

;    base_tlb = widget_base((*state).tlb, $
;                           uname='base_tlb',$
;                           row=1,$
;                           /align_center, $
;                           /base_align_center)
    
    SL_FreqPhase = Widget_Slider(BASE, $
                                 UNAME='SL_FreqPhase', $
                                 XOFFSET=3, $
                                 minimum=-1, $
                                 maximum=matrix[2]-1, $
                                 YOFFSET=3, $ 
                                 TITLE='Frequency/Phase', $
                                 scroll=1, $ 
                                 /drag)
    
    
    SL_PhaseSlice = Widget_Slider(BASE, $
                                  UNAME='SL_PhaseSlice',$
                                  minimum=-1, $
                                  maximum=matrix[0]-1,  $
                                  XOFFSET=3, $
                                  YOFFSET=50,$
                                  TITLE='Phase/Slice', $
                                  scroll=1, $
                                  /drag)
    
    
    SL_FreqSlice = Widget_Slider(BASE, $
                                 UNAME='SL_FreqSlice', $
                                 XOFFSET=3, $
                                 minimum=-1, $
                                 maximum=matrix[1]-1, $
                                 YOFFSET=97, $
                                 TITLE='Frequency/Slice', $
                                 scroll=1, $
                                 /drag)

    SL_ImageAlpha = Widget_Slider(BASE, $
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
                                  event_pro = 'mas_3d_slicer_GUI_event')

    SL_ImageContrast = Widget_Slider(BASE,$
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
                                     event_pro = 'mas_3d_slicer_GUI_event')

    marker_base = widget_base(base, /frame, /column)
    lbl = widget_label(marker_base, value="Markers:")

    li_markers = Widget_List(marker_BASE, uname="marker_list", ysize=10)
    b = widget_base(marker_base, /row)
    btn_add = widget_button(b, uname="btn_add_marker", value="Add...")
    btn_rem = widget_button(b, uname="btn_rem_marker", value="Remove")

    (*state).marker_list_wid = li_markers
    widget_control, /REALIZE, base
    print, "mas_3d_slicer: tlb is", base
    widget_control, base, set_uvalue=state
    xmanager, 'mas_3d_slicer_GUI', base, cleanup='mas_3d_slicer_GUI_cleanup', /no_block

end
; 
; Empty stub procedure used for autoloading.
; 
;pro mas_3d_slicer_GUI, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_, state=state
;  BASE, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_, state=state
;end
