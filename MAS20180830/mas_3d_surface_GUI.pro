;; Subroutine name: mas_3d_surface_GUI_cleanup
;; Created by: KM, 2015-02
;; Calling Information:  
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Cleans up data structures associated with this session
;;
;; Editing Information:

pro mas_3d_surface_GUI_cleanup, tlb

  print, 'mas_3d_surface_GUI_cleanup: ',tlb

   widget_control, tlb, get_uvalue=state
   if (ptr_valid(state)) then begin

      obj_destroy, (*state).model

      ptr_free, state
      
   endif

end

;; Subroutine name: mas_3d_surface_update
;; Created by: KM, 2015-02
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Updates the object's surface rendering using new iso-values specified by the user
;;
;; Editing Information:

pro mas_3d_surface_update, omodel, state

    if (not obj_valid(omodel)) then begin
        print, 'invalid!'
        return
    endif
    
; Destroy old image object    
    old_img = omodel->getByName('Image_object')
    
    if (obj_valid(old_img)) then begin
        obj_destroy, old_img
    endif

; Create new surface rendering with user-specified isovalues
    INTERVAL_VOLUME, *(*state).image_datapool, (*state).min_isoval, (*state).max_isoval, verts, conn
    conn = TETRA_SURFACE(verts, conn)
    omodel_surface = obj_new('IdlgrPolygon', verts, POLYGONS=conn, $
          COLOR=[200,200,200], SHADING=1, name='Surface_image')

; Add surface image model to main model and scale to appropriate size, then display    
    omodel_image = obj_new('idlgrmodel', name='Image_object')    
    omodel_image->add, omodel_surface
    omodel_image->scale, (*state).f_scale, (*state).p_scale, (*state).s_scale
    omodel->add, omodel_image
 
    xobjview, refresh=(*state).tlb

end 

;; Subroutine name: mas_3d_slicer_GUI_event
;; Created by: KM, 2015-02
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;
;;   Event handler for the widgets on the control window
;;
;; Editing Information:

pro mas_3d_surface_GUI_event, Event

    common scan_data, project

    wTarget = (widget_info(Event.id,/NAME) eq 'TREE' ?  $
               widget_info(Event.id, /tree_root) : event.id)

    wWidget =  Event.top
    widget_control, wWidget, get_uvalue=state
    name = widget_info(wTarget, /uname)
    
    case name of
    
        'Min_isovalue': begin
            widget_control, wTarget, get_value=val
            (*state).min_isoval = val
        end
        
        'Max_isovalue': begin
            widget_control, wTarget, get_value=val
            (*state).max_isoval = val
        end

; 'Recalc' button hides the current image object, deletes it, then adds recalculated image object and displays        
        'Recalc': begin
            start_mem = MEMORY(/CURRENT) 
            omodel = (*state).model->getByName('IMG_model/Image_object') 
            omodel->setProperty, hide=1 
            omodel = (*state).model->getByName('IMG_model')
            mas_3d_surface_update, omodel, state   
            mem_required = (MEMORY(/HIGHWATER) - start_mem)/(1024.0)^3    
            PRINT, 'Memory required: ',  mem_required
            if mem_required gt 0.01 then return  
              
        end
        
        else: begin
            print, 'mas_3d_surface_GUI_event: recv;d event from unknown wid: '+name
        end
        
    endcase
end


;; Subroutine name: mas_3d_surface_GUI
;; Created by: Kevin Montes, 2015-02
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

pro mas_3d_surface_GUI, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_, state=state

    matrix = (*state).matrix
    min_isovalue = min(*(*state).image_datapool)
    max_isovalue = max(*(*state).image_datapool)
    
;Create top level base for surface rendering GUI    
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

;Draw widget for voxel density histogram                     
    Plot = WIDGET_DRAW(BASE, UNAME='plot', XSIZE = 320, YSIZE = 320)

;Add buttons for isovalue sliders    
    base_isovalues = Widget_Base(BASE, /ROW, /ALIGN_CENTER)
    
    Min_isovalue = CW_FSlider(base_isovalues, $
                                 UNAME='Min_isovalue', $                                 
                                 maximum=max_isovalue, $
                                 TITLE='Minimum Iso-value', $
                                 Value=(max_isovalue/4),/EDIT)
                                 
    Max_isovalue = CW_FSlider(base_isovalues, $
                                 UNAME='Max_isovalue', $
                                 maximum=max_isovalue, $ 
                                 TITLE='Maximum Iso-value', $
                                 Value=3*(max_isovalue/4),/EDIT)
                                 

;Add button to recalculate the surface rendering using the user-input iso-values                                 
    button_base = Widget_Base(BASE, /ALIGN_CENTER, /ROW)
    
    Recalc = Widget_Button(button_base, VALUE='Recalculate', UVALUE='Recalc', UNAME='Recalc', $
                          /ALIGN_CENTER, TOOLTIP='Display surface rendering with new iso-values', $
                          event_pro = 'mas_3d_surface_GUI_event')

;Realize widget and display histogram
    widget_control, /REALIZE, base
    print, "mas_3d_slicer: tlb is", base
    widget_control, base, set_uvalue=state
    
    widget_control, plot, GET_VALUE = plotid
    wset, plotid
    plot, histogram(*(*state).image_datapool), $
    TITLE='Voxel Density for Selected Image', $
    XTITLE='Voxel Intensity', $
    YTITLE='Number of Voxels'
    
    xmanager, 'mas_3d_surface_GUI', base, cleanup='mas_3d_surface_GUI_cleanup', /no_block
    
end