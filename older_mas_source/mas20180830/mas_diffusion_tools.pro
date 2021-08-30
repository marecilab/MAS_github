;; $Id$
;;

;; Subroutine name: mdt_get_b_matrix
;; Created by: BT, Jan 2009
;; Calling Information:
;;   dir_index    - the "adim" index for the b-matrix
;;
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine: returns a properly formatted symmetric
;; b-matrix for the requested adim (diffusion weighted volume) index
;;
;; Editing Information:

function mdt_get_b_matrix, dir_index

    common scan_data

    imnd_b = reform((*project.imndarray[project.ci].b_matrix)[*,dir_index])

    ;; these entries are multiplied by two 
    ;; when the matrix is read by mas_open
    imnd_b[3:5] /= 2.
    
    B = [ [imnd_b[0], imnd_b[3], imnd_b[4] ], $
          [imnd_b[3], imnd_b[1], imnd_b[5] ], $
          [imnd_b[4], imnd_b[5] ,imnd_b[2] ] ]
    
    return, B
 end

;; Subroutine name: mdt_make_default_state
;; Created by: BT, Jan 2009
;; Calling Information:
;;    < no paramteters >
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine: returns a state pointer that is used to
;; store paramters for the diffusion tools
;;
;; Editing Information:

function mdt_make_default_state

    return, ptr_new({ DIFFUSION_TOOLS_STATE, $
                      adt_base: 0L, $;; tlb for adt tab
                      dot_base: 0L, $
                      mow_base: 0L, $
                      odf_base: 0L, $
                      glyph_size: 12,  $
                      glyph_size_mult: 0.5, $
                      glyph_shiny: 4,  $
                      l0_intensity: 60,$
                      l1_intensity: 25,$
                      l2_intensity: 0, $
                      use_adt_color: 0, $
                      glyph_opacity: 100, $
                      image_opacity: 60, $
                      disp_3d:0, $
                      mask_roi: 0, $
                      crop_roi: 0, $
                      apply_acq_matrix: 1}, /no_copy)
end

;; Subroutine name: as_diffusion_tools_GUI
;; Created by: BT, Jan 2009
;; Calling Information:
;;    < no paramteters >
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine: sets up the GUI for the mas_diffusion_tools
;; main window. 
;;
;; Editing Information:

pro mas_diffusion_tools_GUI

    common scan_data
    common common_widgets

    ci = project.ci

    if (not ptr_valid(project.imndarray[ci].bval_array)) then begin
       void = dialog_message('Not a diffusion scan.', /error, /center)
       return
    endif

    mas_load_state_1

    ;; make sure the window isn't already open; if so redraw it.
    if (project.difftools_tlb gt 0) then begin
        mdt_redraw        
        widget_control, project.difftools_tlb, /show
        return
    endif

    ;; check that this scan has a valid state pointer, otherwise
    ;; create one.
    if (ptr_valid(project.procpramarray[ci].difftools_state)) then begin
        difftools_state = project.procpramarray[ci].difftools_state
    endif else begin
        difftools_state = mdt_make_default_state()
        project.procpramarray[ci].difftools_state = difftools_state
    endelse
 
    main = widget_base(group_leader=WID_BASE_MAIN, $
                       title='MAS Diffusion Visualization Tools', $
                       row=2);, /grid)

    ;; in many of these widgets, the uvalue (one of
    ;; TAB/RENDER/DISPLAY) is used to indicate what purpose this
    ;; widget serves. Mainly for the event handler.
    tab_base = widget_tab(main, uvalue='TAB')

    ;; the gui code for the respective tabs are in the respective .pro files
    adt_tab_base = widget_base(tab_base, title=' ADT ', uname='ADT_tab', xsize=400)
    mas_adt_sp_gui_make, adt_tab_base, 0
    (*difftools_state).adt_base = adt_tab_base
    
    dot_tab_base = widget_base(tab_base, title=' DOT ', uname='DOT_tab')
    mas_dot_gui_make, dot_tab_base, 0
    (*difftools_state).dot_base = dot_tab_base

    mow_tab_base = widget_base(tab_base, title=' MOW ', uname='MOW_tab')
    mas_mow_gui_make, mow_tab_base, 0
    (*difftools_state).mow_base = mow_tab_base

    odf_tab_base = widget_base(tab_base, title=' ODF ', uname='ODF_tab')
    mas_odf_gui_make, odf_tab_base, 0
    (*difftools_state).odf_base = odf_tab_base

    render_base = widget_base(main, xpad=2, ypad=2, /row)

    hiq_base = widget_base(render_base, column=2, /grid, /frame)
 
    ;; "size" of each glyph
    glyph_size = widget_slider(hiq_base, $
                               uname   = 'glyph_size', $
                               title   = 'Glyph Size',  $
                               minimum = 8, $
                               uvalue  = 'RENDER', $
                               maximum = 128, $
                               value   = (*difftools_state).glyph_size, $
                               scroll  = 2)

    glyph_shiny = widget_slider(hiq_base, $
                                minimum = 0, $
                                maximum = 25, $
                                value   = (*difftools_state).glyph_shiny, $
                                scroll  = 1, $
                                uvalue  = 'RENDER', $
                                uname   = 'glyph_shiny', $
                                title   = "Shininess")
    

    l0_intensity = widget_slider(hiq_base, $
                                 minimum = 0, $
                                 maximum = 100, $
                                 value   = (*difftools_state).l0_intensity, $
                                 scroll  = 1, $
                                 uvalue  = 'RENDER', $
                                 uname   = 'l0_intensity', $
                                 title   = "Overhead Light")

    l1_intensity = widget_slider(hiq_base, $
                                 minimum = 0, $
                                 maximum = 100, $
                                 value   = (*difftools_state).l1_intensity, $
                                 scroll  = 1, $
                                 uvalue  = 'RENDER', $
                                 uname   = 'l1_intensity', $
                                 title   = "Headlight")

    l2_intensity = widget_slider(hiq_base, $
                                 minimum = 0, $
                                 maximum = 100, $
                                 scroll  = 1, $
                                 uvalue  = 'RENDER', $
                                 value   = (*difftools_state).l2_intensity, $
                                 uname   = 'l2_intensity', $
                                 title   = "Rear Light")
    
    glyph_opacity = widget_slider(hiq_base, $
                                  minimum = 0, $
                                  maximum = 100, $
                                  scroll  = 1, $
                                  uvalue  = 'RENDER', $
                                  value   = (*difftools_state).glyph_opacity, $
                                  uname   = 'glyph_opacity', $
                                  title   = 'Glyph Opacity')
    
    image_opacity = widget_slider(hiq_base, $
                                  minimum = 0, $
                                  maximum = 100, $
                                  scroll  = 1, $
                                  uvalue  = 'RENDER', $
                                  value   = (*difftools_state).image_opacity, $
                                  uname   = 'image_opacity', $
                                  title   = 'Image Opacity')

    main_c2 = widget_base(render_base, /column)
         
    non_ex_base = widget_base(main_c2, /align_left, /column, /nonexclusive)
    btn_use_adt_color = widget_button(non_ex_base, $
                                      value="ADT Orientation Color", $
                                      uname='btn_use_adt_color', $
                                      uvalue='RENDER', $
                                      sensitive=project.procpramarray[ci].adt_proccess_flag)

    widget_control, btn_use_adt_color, set_button=(*difftools_state).use_adt_color
    
    btn_apply_roi = widget_button(non_ex_base, $
                                  value="Apply Active ROI", $
                                  uvalue='RENDER', $
                                  uname='btn_apply_roi')
    widget_control, btn_apply_roi, set_button=(*difftools_state).mask_roi
    
    btn_crop_roi = widget_button(non_ex_base, $
                                 value="Crop to ROI", $
                                 uvalue='RENDER', $
                                 uname='btn_crop_roi', xoffset=5)
    widget_control, btn_crop_roi, set_button=(*difftools_state).crop_roi, sensitive=(*difftools_state).mask_roi

    btn_use_3d = widget_button(non_ex_base, $
                               value="Display in 3D", $
                               uvalue='RENDER', $
                               uname='btn_use_3d', xoffset=5)
    widget_control, btn_use_3d, set_button=(*difftools_state).disp_3d, sensitive=1

    btn_use_acq_matrix = widget_button(non_ex_base, $
                                       value="Apply Acquisiton Matrix", $
                                       uvalue='RENDER', $
                                       uname='btn_use_acq_matrix', xoffset=5)
    widget_control, btn_use_acq_matrix, set_button=(*difftools_state).apply_acq_matrix, sensitive=1

    btn_base = widget_base(main_c2, /column)
    btn_create = widget_button(btn_base, $
                               value='Create Image', $
                               event_pro='mas_adt_sp_configure_event', $
                               uvalue='DISPLAY', $
                               uname='btn_create')
    btn_batch  = widget_button(btn_base, $
                               value='Batch Images', $
                               event_pro='mas_adt_sp_configure_event', $
                               uvalue='DISPLAY', $
                               uname='btn_batch')
    
;;    active_method: Currently selected tab. ADT/MOW/DOT/ODF
;;    adt_base: tlb for adt tab gui
;;    dot_base: tlb for dot tab gui
;;    mow_base: tlb for mow tab gui
;;    odf_base: tlb for odf tab gui
;;    btn_create: widget id for "create" button
;;    btn_batch:  widget if for "batch" button
;;    glyph_size:  the rest are widget ids for their resp. widget
;;    glyph_shiny: 
;;    l0_intensity:
;;    l1_intensity: 
;;    l2_intensity: 
;;    glyph_opacity:
;;    image_opacity:
;;    btn_crop_roi: 
;;    btn_mask_roi:
;;    btn_use_3d: 
;;    btn_use_adt_color:
    
    widget_state = ptr_new({$
                             active_method: 'ADT', $
                             adt_base: adt_tab_base, $
                             dot_base: dot_tab_base, $
                             mow_base: mow_tab_base, $
                             odf_base: odf_tab_base, $
                             btn_create: btn_create, $
                             btn_batch:  btn_batch, $
                             glyph_size: glyph_size,  $
                             glyph_shiny: glyph_shiny,  $
                             l0_intensity: l0_intensity,  $
                             l1_intensity: l1_intensity,  $
                             l2_intensity: l2_intensity,  $
                             glyph_opacity: glyph_opacity,  $
                             image_opacity: image_opacity,  $
                             btn_crop_roi: btn_crop_roi, $
                             btn_mask_roi: btn_apply_roi, $
                             btn_use_3d: btn_use_3d, $
                             btn_use_adt_color: btn_use_adt_color }, /no_copy)

    widget_control, main, set_uvalue=widget_state
    project.difftools_tlb = main

    widget_control, main, /realize
    
    xmanager, 'mas_diffusion_tools', main, /no_block, cleanup='mdt_cleanup'

end

;; Subroutine name: mdt_redraw
;; Created by: BT, Jan 2009
;; Calling Information:
;;    < no paramteters >
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine: redraws the diffusion tools window
;;
;; Editing Information:

pro mdt_redraw

    common scan_data

    ci = project.ci

    if (project.difftools_tlb lt 0) then begin
;        print, 'mdt_redraw: nothing to do'
        return
    endif else if not widget_info(project.difftools_tlb, /valid_id) then begin
        print, 'mdt_redraw: inconsistent tlb ld'
        project.difftools_tlb = -1
        return
    endif else begin
        widget_control, project.difftools_tlb, get_uvalue=widget_state
        if not ptr_valid(widget_state) then begin
            print, 'mdt_redraw: tlb has no state information'
            return
        endif
    endelse

    difftools_state = project.procpramarray[ci].difftools_state
    sen = 1
    free_state_when_done = 0

    if (not ptr_valid(difftools_state)) then begin

        print, 'mdt_redraw: no state ptr for this scan'
        difftools_state = mdt_make_default_state()

        if (project.ni eq 0) then begin
            sen = 0
            free_state_when_done = 1
        endif else begin
            project.procpramarray[ci].difftools_state = difftools_state
        endelse
    endif

    adt_processed = project.procpramarray[ci].adt_proccess_flag

    if (adt_processed eq 0) then adt_sen = 0 else adt_sen = 1

    widget_control, (*widget_state).btn_create, sensitive=sen 
    widget_control, (*widget_state).btn_batch, sensitive=sen
    widget_control, (*widget_state).glyph_size, sensitive=sen, set_value=(*difftools_state).glyph_size
    widget_control, (*widget_state).glyph_shiny, sensitive=sen, set_value=(*difftools_state).glyph_shiny
    widget_control, (*widget_state).l0_intensity, sensitive=sen, set_value=(*difftools_state).l0_intensity
    widget_control, (*widget_state).l1_intensity, sensitive=sen, set_value=(*difftools_state).l1_intensity
    widget_control, (*widget_state).l2_intensity, sensitive=sen, set_value=(*difftools_state).l2_intensity
    widget_control, (*widget_state).glyph_opacity, sensitive=sen, set_value=(*difftools_state).glyph_opacity
    widget_control, (*widget_state).image_opacity, sensitive=sen, set_value=(*difftools_state).image_opacity
    widget_control, (*widget_state).btn_mask_roi, sensitive=sen, set_button=(*difftools_state).mask_roi
    widget_control, (*widget_state).btn_use_adt_color, sensitive=adt_sen*sen, set_button=(*difftools_state).use_adt_color
    widget_control, (*widget_state).btn_crop_roi, sensitive=sen*(*difftools_state).mask_roi, set_button=(*difftools_state).crop_roi
   
    if (free_state_when_done) then ptr_free, difftools_state

    if (ptr_valid(project.imndarray[ci].bval_array)) then begin
       mas_odf_redraw, (*widget_state).odf_base
       mas_adt_sp_redraw, (*widget_state).adt_base
       mas_mow_redraw, (*widget_state).mow_base
       mas_dot_redraw, (*widget_state).dot_base
    endif

end

;; Subroutine name: mdt_cleanup
;; Created by: BT, Jan 2009
;; Calling Information:
;;    < no paramteters >
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine: frees pointers and such when tlb is destroyed
;;
;; Editing Information:

pro mdt_cleanup, tlb

    common scan_data
    print, 'mdt_cleanup: cleaning up:', tlb

    widget_control, tlb, get_uvalue=widget_state 
    if (ptr_valid(widget_state)) then begin
        ptr_free, widget_state
    endif

    ;; if, for some reason, the project_tlb != the tlb to be
    ;; destroyed, then we just print a message. This has never
    ;; happened, but if it does, then there is a bug somewhere.
    if (project.difftools_tlb ne tlb) then begin
        print, "difftools tlb not consistent!"
    endif

    project.difftools_tlb = -1

end

;; Subroutine name: mas_diffusion_tools_EVENT
;; Created by: BT, Jan 2009
;; Calling Information:
;;    < no paramteters >
;;
;; Bugs or Important Comments to Developers: Note that in this
;; subroutine, events in tabs are re-routed to the event handlers
;; within the .pro files for those visualization methods.
;;
;; the widget_state.active method holds which tab is currently active.
;;
;; Purpose of subroutine: event handler for diffusion tools events.
;;
;; Editing Information:

pro mas_diffusion_tools_EVENT, event

    common scan_data

    ci = project.ci
    difftools_state = project.procpramarray[ci].difftools_state

    if (not ptr_valid(difftools_state)) then begin
        print, "No Diff Tools State"
        return
    endif

    ;; get widget state pointer
    widget_control, event.top, get_uvalue=widget_state
    ;; get "purpose" class of widget
    widget_control, event.id, get_uvalue=uvalue
    uname = widget_info(event.id, /uname)

    ;; must be a visualization specific event... rethrow event.
    if (n_elements(uvalue) eq 0) then begin

        case (*widget_state).active_method of
            'MOW': mas_mow_gui_event, event
            'ADT': mas_adt_sp_configure_event, event
            'DOT': mas_dot_gui_event, event
            'ODF': mas_odf_gui_event, event
        endcase
        ;; must be a visualization-specific widget
        ;; we will update the event.top and re-throw the event to 
        ;; the visualization-specific handler

        ;;  ...

        return

     endif else if (uvalue eq 'TAB') then begin
        ;; a tab focus has been changed. need to update the event_pros
        ;; to route events to the respective visualization handlers
        current_tab = event.tab
        children = widget_info(event.id, /all_children)
        uname = widget_info(children[current_tab], /uname)
        case uname of 
            'ADT_tab': begin
                (*widget_state).active_method = 'ADT'
                widget_control, (*widget_state).btn_use_3d, sensitive=1
                widget_control, (*widget_state).btn_create, event_pro='mas_adt_sp_configure_event'
                widget_control, (*widget_state).btn_batch,  event_pro='mas_adt_sp_configure_event'
            end

            'MOW_tab': begin
                (*widget_state).active_method = 'MOW'
                widget_control, (*widget_state).btn_use_3d, sensitive=1
                widget_control, (*widget_state).btn_create, event_pro='mas_mow_gui_event'
                widget_control, (*widget_state).btn_batch,  event_pro='mas_mow_gui_event'
            end

            'DOT_tab': begin
                (*widget_state).active_method = 'DOT'
                widget_control, (*widget_state).btn_use_3d, sensitive=1
                widget_control, (*widget_state).btn_create, event_pro='mas_dot_gui_event'
                widget_control, (*widget_state).btn_batch,  event_pro='mas_dot_gui_event'
            end

            'ODF_tab': begin
                (*widget_state).active_method = 'ODF'
                widget_control, (*widget_state).btn_use_3d, sensitive=1
                widget_control, (*widget_state).btn_create, event_pro='mas_odf_gui_event'
                widget_control, (*widget_state).btn_batch,  event_pro='mas_odf_gui_event'
            end

        endcase

        return

    endif else if (uvalue eq 'RENDER') then begin
        ;; rendering options have been changed. these affect all visualizations
        widget_control, event.id, get_value=val
        case uname of 
            'glyph_size'     : (*difftools_state).glyph_size  = long(val)
            'glyph_shiny'    : (*difftools_state).glyph_shiny = float(val)
            'glyph_opacity'  : (*difftools_state).glyph_opacity = float(val)
            'image_opacity'  : (*difftools_state).image_opacity = float(val)
            'l0_intensity'   : (*difftools_state).l0_intensity = float(val)
            'l1_intensity'   : (*difftools_state).l1_intensity = float(val)
            'l2_intensity'   : (*difftools_state).l2_intensity = float(val)
            'btn_use_adt_color': (*difftools_state).use_adt_color = event.select
            'btn_use_3d'     : (*difftools_state).disp_3d = event.select
            'btn_use_acq_matrix': (*difftools_state).apply_acq_matrix = event.select
            'btn_apply_roi'  : begin
                (*difftools_state).mask_roi    = event.select
                widget_control, (*widget_state).btn_crop_roi, sensitive=event.select
            end
            
            'btn_crop_roi'   : (*difftools_state).crop_roi    = event.select
            else: print, 'mdt_event: unknown RENDER event from: ', uname

        endcase

    endif else if (uvalue eq 'DISPLAY') then begin

        case uname of 
            'btn_create':
            'btn_batch':
            else: print, 'mdt_event: unknown DISPLAY event from: ', uname
        endcase

    endif else begin
        
        print, "UNKNOWN EVENT"
        
    endelse
        
end

;; Subroutine name: mdt_get_adt_background_image
;; Created by: BT, Jan 2009
;; Calling Information:
;;    
;;    WANT_SLICE: Slice requested by the user
;;    GOT_SLICE:  Slice returned. Used when WANT_SLICE may not be
;;                specified
;;    BG_IMAGE:   Set to a named variable that will receive the image
;;                data
;;    USE_MY_IMAGE: Set this keyword to a variable that contains image
;;                  data supplied by the user. The image data will be
;;                  cropped and masked according to the other parameters
;;                  and returned as BG_IMAGE
;;    MASK_ROI:     Set this keyword to have the image masked by the
;;                  currently selected ROI in the MAS main window
;;    DATA_MASK:    Set this to a named variable to have the resulting
;;                  mask returned.
;;    CROP_ROI:     Set this keyword to have the image cropped to the
;;                  currently selected ROI. Outer-most bounds are used
;;                  for cropping.
;;    {X,Y,Z}DIM:   Set these keyword to a named variable that will
;;                  contain the dimensions of the resulting image
;;    {X,Y}START:   Contains the starting (x,y) coordinate in the
;;                  original data that corresponds to the (0,0)
;;                  coordinate of the cropped data .
;;    SLICE_DATA:   Set this to a named keyword that will contain the
;;                  slice data (original voxel data) for the requested
;;                  slice.
;;    FA:           Set this to a named variable that will receive an
;;                  FA image of the current slice (for FA thresh.)
;;
;; Bugs or Important Comments to Developers: 
;;
;;  This really serves a very specific purpose and could probably be
;;  integrated into the calling procedure.
;;
;; Purpose of subroutine: To acquire a background image upon which to
;;                        overlay diffusion visualization glyphs.
;;
;; Editing Information:
;;

pro mdt_get_adt_background_image, want_slice=want_slice, $
                                  got_slice=got_slice, $
                                  bg_image=bg_image, $
                                  use_my_image=use_my_image, $
                                  mask_roi=mask_roi, $
                                  fa=fa, $
                                  slice_data=slice_data, $
                                  xdim=xdim, ydim=ydim, zdim=zdim, $
                                  xstart=xstart, ystart=ystart, $
                                  data_mask=data_mask, crop_roi=crop_roi
  
    common scan_data

    ci = project.ci
    sz_img = size(*project.dataarray[ci].state1, /dimensions)
    have_fa = ptr_valid(project.dataarray[ci].frac_ani)

    case project.procPramArray[project.ci].slice_axis of
        0: begin                ; freq_phase
            xdim = sz_img[0]
            ydim = sz_img[1]
            zdim = sz_img[2]
            xfov = project.imndarray[ci].f_fov
            yfov = project.imndarray[ci].p_fov
            slice      = keyword_set(want_slice) ? want_slice : project.procpramArray[ci].sdim_start
            slice_data = reform((*project.dataarray[ci].state1)[*,*,slice,*])
            if (have_fa) then fa = reform((*project.dataarray[ci].frac_ani)[*,*,slice])
        end
        
        1: begin                ; freq_slice
            xdim = sz_img[0]
            ydim = sz_img[2]
            zdim = sz_img[1]
            xfov = project.imndarray[ci].f_fov
            yfov = project.imndarray[ci].s_fov
            slice      = keyword_set(want_slice) ? want_slice : project.procpramArray[ci].pdim_start
            slice_data = reform((*project.dataarray[ci].state1)[*,slice,*,*])
            if (have_fa) then fa = reform((*project.dataarray[ci].frac_ani)[*,slice,*])
        end
        
        2: begin                ; phase_slice
            xdim = sz_img[1]
            ydim = sz_img[2]
            zdim = sz_img[0]
            xfov = project.imndarray[ci].p_fov
            yfov = project.imndarray[ci].s_fov
            slice      = keyword_set(want_slice) ? want_slice : project.procpramArray[ci].fdim_start
            slice_data = reform((*project.dataarray[ci].state1)[slice,*,*,*])
            if (have_fa) then fa = reform((*project.dataarray[ci].frac_ani)[slice,*,*])
        end
    endcase

    xstart = 0
    ystart = 0

    if (keyword_set(use_my_image)) then begin
       help, use_my_image
       tmp = bytscl(use_my_image)
    endif else begin
       tmp = adt_make_image(slice,0)
    endelse

    case ((size(tmp))[0]) of
       2: begin
          ;; adt_make_image returned a b/w image, use that
          bg_image = temporary(tmp)
          
          if (project.procpramarray[ci].adt_display_type eq 4) then begin
             ;; adt_make_image is inconsistent when it comes to the
             ;; diffusion-weighted image, so we un-de-rotate it
             if (project.procpramarray[ci].flip_direction gt 0) then begin
                bg_image = reverse(bg_image, project.procpramarray[ci].flip_direction)
             endif
             
             ;; inverse rotation
             bg_image = rotate(bg_image, (4-(project.procpramarray[ci].rotate_direction)) mod 4)
          endif
       end
       3: begin
          bg_image = bytscl(temporary(tmp))
       end
       else: begin
          junk = dialog_message(['Could not acquire ADT image.', $
                                 'Please check ADT display selection'], $
                                /center, /error)
          use_derived_anis_image = 1
       end
    endcase

    ;; Apply ROI Masking and cropping, recompute dimensions
    if (keyword_set(mask_roi)) then begin
        if (keyword_set(crop_roi)) then begin
            data_mask = mas_roi_get_current_mask([xdim, ydim], crop=crop_dims, /no_transform)
        endif else begin
            data_mask = mas_roi_get_current_mask([xdim, ydim], /no_transform)
        endelse

        if (not ptr_valid(data_mask)) then begin
            mask_roi = 0
            if (n_elements(crop_dims) ne 0) then begin
                junk = temporary(crop_dims)
            endif
        endif
    endif

    if (n_elements(crop_dims) ne 0) then begin
        crop = temporary(bg_image[ crop_dims[0]:crop_dims[2]-1, $
                                   crop_dims[1]:crop_dims[3]-1, $
                                   * ])
        bg_image = temporary(crop)

        crop = (*data_mask)[ crop_dims[0]:crop_dims[2]-1, $
                        crop_dims[1]:crop_dims[3]-1 ]
        ptr_free, data_mask
        data_mask = ptr_new(crop, /no_copy)
        
        old_xdim = xdim
        old_ydim = ydim

        xdim = crop_dims[2]-crop_dims[0]
        ydim = crop_dims[3]-crop_dims[1]
        
        xfov *= float(xdim)/float(old_xdim)
        yfov *= float(ydim)/float(old_ydim)

        slice_data_dims = size(slice_data, /dimensions)
        
        sd = fltarr(xdim, ydim, slice_data_dims[2])
        sd[*,*,*] = slice_data[crop_dims[0]:crop_dims[2]-1, $
                               crop_dims[1]:crop_dims[3]-1, * ]
        
        slice_data = temporary(sd)

        xstart = crop_dims[0] ;; used for indexing into the FA, Eval, Evec array
        ystart = crop_dims[1] ;; when glyph coloring is enabled

    endif
        
    got_slice = slice

end

;; Subroutine name: mdt_make_glyphimage
;; Created by: BT, Jan 2009
;;
;; Calling Information:
;;    
;;  RECO_OBJ:  A diffusion bulk reconstructor object, such as: 
;;             MAS_MOW_BULK_RECONSTRUCTOR,
;;             MAS_ODF_BULK_RECONSTRUCTOR,
;;             MAS_DOT_BULK_RECONSTRUCTOR.
;;
;;  VDIM:      A number between 8 and 64 that indicates the glyph
;;             size. In particular, this size cooresponds to the 1/2
;;             the radius of an undistorted sphere centered at the
;;             origin.     
;;  
;;  VMULT:     A multiplier that is used to simplify the spherical
;;             mesh. A mesh is created using VDIM/VMULT size, then the
;;             resulting surface size is increased.
;;  
;;  ANIS_THRESHOLD: A cutoff value of FA. Voxels having FA below this
;;                  cutoff will be ignored
;;
;;  THRESHOLD: A cutoff value for the voxel data itself.
;;
;;  VIS_SLICE: the requested slice in the currently chosen slice axis
;;             in the mas main window. If not present, the currently
;;             selected slice in the mas main window will be used.
;;
;;  IMAGE_IN: Image data to be used as a background for the glyph
;;            image. If not present, the currently selected ADT image
;;            will be used.
;;
;;  IMAGE_OUT: Set this to a named variable that will contain the
;;             resulting image.
;;
;;  DISPLAY:  Set this keyword to display the image using MAS's
;;            tabbed image display routines
;;
;; Bugs or Important Comments to Developers: 
;;
;;
;;
;; Purpose of subroutine: To acquire a background image upon which to
;;                        overlay diffusion visualization glyphs.
;;
;; Editing Information:
;;

pro mdt_make_glyphimage, $
   reco_obj, $
   vdim=vdim,  $
   vmult=vmult, $
   vis_slice=vis_slice, $
   anis_threshold=anis_threshold, $
   threshold=threshold, $
   display=display, $
   image_in=image_in, $
   image_out=image_out

    common scan_data
    forward_function mas_dot_color_mapper
    forward_function get_orient_rgb_axis

    ci = project.ci

    mas_load_state_1

    if (not ptr_valid(project.imndArray[ci].bval_Array)) then begin
       void = dialog_message(['B values are not present for this scan.', $
                              'Please make sure that this is a diffusion-',$
                              'weighted scan.'], /error, /center)
       return
    endif

    if n_elements(vdim) eq 0      then vdim = 10
    if n_elements(vmult) eq 0     then vmult = 1.0
    if n_elements(threshold) eq 0 then begin
       threshold = 0.0
    endif 
    if n_elements(anis_threshold) eq 0 then begin
       anis_thr = 1.0
    endif else begin
       anis_thr = anis_threshold
    endelse

    ;; this forces it to always display the image.
    if n_elements(display) eq 0  then display = 1

    if (ptr_valid(project.procpramarray[ci].difftools_state)) then begin
       difftools_state = project.procpramarray[ci].difftools_state
    endif else begin
       difftools_state = mdt_make_default_state()
    endelse

    use_derived_anis_image = 0
    crop_roi     = (*difftools_state).crop_roi
    mask_roi     = (*difftools_state).mask_roi
    vdim         = (*difftools_state).glyph_size
    vmult        = (*difftools_state).glyph_size_mult
    l0_intensity = (*difftools_state).l0_intensity / 100.0
    l1_intensity = (*difftools_state).l1_intensity / 100.0
    l2_intensity = (*difftools_state).l2_intensity / 100.0
    image_opac   = (*difftools_state).image_opacity   / 100.0
    glyph_opac   = (*difftools_state).glyph_opacity   / 100.0
    glyph_shiny  = (*difftools_state).glyph_shiny
    use_adt_color = (*difftools_state).use_adt_color
    apply_acq_matrix = (*difftools_state).apply_acq_matrix
    
    ;; incorporate freq, phase, slice viewing direction
    sz_img = size(*project.dataarray[ci].state1, /dimensions)

    ;; used for fa thresholding
    have_fa = project.procpramarray[ci].adt_proccess_flag
    
    ;; acquire the background image from ADT, mask, crop, etc
    mdt_get_adt_background_image, $
       use_my_image=image_in, $
       want_slice=vis_slice, $
       mask_roi=keyword_set(mask_roi), $
       crop_roi=keyword_set(crop_roi), $
       bg_image=bg_image, fa=fa, $
       slice_data=slice_data, $
       xdim=xdim, ydim=ydim, zdim=zdim, $
       got_slice=sl, data_mask=mask, $
       xstart=xstart, ystart=ystart

    if (n_elements(bg_image) eq 0) then begin
        ;; user id notified by above procedure.
        return
    endif
    
    ;; recompute optimal vdim based on image dimensions
    vdim_opt = mas_glyphscene_renderer_getOptimalVdim(xdim, ydim, vdim)

    if (vdim_opt ne vdim) then begin
       vdim = vdim_opt
    endif

    ;; at least 20 vertices are required for a decent shape
    if (vdim lt 24) then vmult = 1.0

    gs = obj_new('mas_glyphscene_renderer', bg_image=bg_image, $
                 vmult=vmult, $
                 vdim=vdim, $
                 light0_int=l0_intensity, $
                 light1_int=l1_intensity, $
                 light2_int=l2_intensity, $
                 bg_img_opac=image_opac, $
                 glyph_opac=glyph_opac, $
                 /skip_bufsize_warnings)

    vdim_effective = vmult*vdim

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
    ;mas_tessellator_make_sphere, vdim=vdim_effective, $
    ;                            vertlist=vertlist, $
    ;                            polylist=polylist, /normalize
    mas_tessellator_make_icos, level=3, vertlist=vertlist, polylist=polylist
     
    junk = get_orient_rgb_axis(tx_matrix=tx)
    

    reco_obj->setReconstructionDirections, vertlist
    ;;if (0 and obj_class(reco_obj) ne 'MAS_MOW_BULK_RECONSTRUCTOR') then begin
    if (apply_acq_matrix ne 0) then begin
       ;; The mow is known to use the b-matrix, which is rotated by
       ;; Bruker and/or MAS's PAR/REC,Philips DICOM readers on import
       tx = *project.imndarray[project.ci].acq_matrix # tx
       
    endif
    vertlist = vert_t3d(vertlist, matrix=tx, /no_copy, /no_divide, /double)
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    progressBar = obj_new('progressbar', title='Computing Probabilities', $
                          text='Computing Probabilities...', /fast_loop)
    progressBar->start

    print, "Creating surfaces..."
    syst = systime(1)

    ;; will contain the rendered image
    prendered_img = ptr_new(bytarr(xdim*vdim, ydim*vdim, 3), /no_copy)

    ;; compute thresholding
    max_signal = (project.procpramarray[ci].threshold)[1]
    pct_threshold = max_signal * threshold/1000.0
    print, "Threshold value = ", pct_threshold

    ;; holds the surfaces indexed by their (x,y) location
    osurfs = objarr(xdim,ydim)
    iter = 0

    ;; default coloring
    glyph_color    = [50,50,100] ;fltarr(3)
    glyph_diffuse  = [90,90,110]
    glyph_specular = [150,150,170]
    glyph_scaling  = 'min-max'
    shift_in_x = 0

    ;; load the RGB color table for anisotropy-color scaled glyphs
    loadct, 25, rgb_table=glyph_coltab
    nvertices = n_elements(vertlist)/3
    vert_colors = transpose([ [replicate(1,nvertices)], $
                              [replicate(1,nvertices)], $
                              [replicate(1,nvertices)] ])
        
;;;;;;;;;;;;;;;;;;;;;;;;;

    for x = 0,xdim-1 do begin
       for y = 0,ydim-1 do begin

          if (keyword_set(mask_roi) and ptr_valid(mask)) then begin
             if ((*mask)[x,y] eq 0) then continue
          endif
          
          if (slice_data[x,y,0] lt pct_threshold) then begin
             continue
          endif

          if (have_fa) then begin
             if (fa[x+xstart,y+ystart] lt anis_thr) then continue
          endif

          Pr_orig = reco_obj->getDisplacementProbability(data=reform(slice_data[x,y,*]))

          if ( total(finite(Pr_orig, /nan) + total(finite(Pr_orig, /inf)) ) gt 0) then begin
             continue
          endif

          ;; minmax adjust to accentuate directionality
          max_Pr = max(Pr_orig, min=min_Pr)
          Pr_orig_over_max_pr = Pr_orig/max_Pr
          
          if (glyph_scaling eq 'min-max') then begin
              Pr_scl = (Pr_orig - min_Pr)/(max_Pr-min_Pr) * (vdim_effective * 0.48)
          endif else begin
              Pr_scl = (vdim_effective * 0.48) * (Pr_orig_over_max_pr) * fa[x+xstart,y+ystart]
          endelse
          
          verts = vertlist
          verts[0,*] *= Pr_scl & verts[1,*] *= Pr_scl & verts[2,*] *= Pr_scl 

          ;; note that shift_in_x is zero since the glyph renderer
          ;; creates the image line-by-line in the y-direction:
          verts[0,*] += ((shift_in_x*x)*vdim_effective + vdim_effective/2.)
          verts[1,*] += (y*vdim_effective + vdim_effective/2.)

          if (keyword_set(use_adt_color) and $
              project.procpramarray[ci].adt_proccess_flag) then begin
              glyph_color = adt_get_orientation_color(x+xstart, y+ystart, sl) 
              vert_colors = transpose([ [replicate(glyph_color[0],nvertices)], $
                                        [replicate(glyph_color[1],nvertices)], $
                                        [replicate(glyph_color[2],nvertices)] ])
          endif else begin
              vert_colors[0,*] = glyph_coltab[fix(Pr_orig_over_max_pr * 255), 0]
              vert_colors[1,*] = glyph_coltab[fix(Pr_orig_over_max_pr * 255), 1]
              vert_colors[2,*] = glyph_coltab[fix(Pr_orig_over_max_pr * 255), 2]
          endelse
          
          tmp = obj_new('idlgrpolygon', verts, $
                        polygons=polylist,  $
                        shading=1,  $
                        style=2, $
                        diffuse=glyph_diffuse, $
                        vert_colors=vert_colors, $
                        ;color=glyph_color, $
                        $;specular=glyph_specular,$
                        alpha_channel=glyph_opac, $
                        reject=0, $
                        shininess=glyph_shiny, $
                        depth_test_disable=2)

          gs->addGlyph, tmp, x, y, /replace
          osurfs[x,y] = tmp

       endfor

       progressbar->update, (100.0) * x/xdim
       if (progressBar->checkcancel()) then begin
          progressBar->destroy
          obj_destroy, gs
          if (ptr_valid(mask)) then ptr_free, mask
          obj_destroy, osurfs
          obj_destroy, reco_obj
          return
       endif
       
    endfor
    
    progressBar->destroy

    print, '... Done. (Time = '+string(systime(1)-syst, format='(F0.5)')+')'
    syst = systime(1)

    print, 'Rendering...'

    gs->render, rendered_img=image_out, use_this_ptr=prendered_img

    if (not ptr_valid(prendered_img)) then return

    image_out = *prendered_img
    ptr_free, prendered_img

    obj_destroy, gs
    obj_destroy, osurfs

    if (ptr_valid(mask)) then ptr_free, mask

    print, '... Done. (Time = '+string(systime(1)-syst, format='(F0.5)')+')'

    if (n_elements(image_out) eq 0) then return

    pimage_out = ptr_new(image_out, /no_copy)
    mas_rotate_flip, pimage_out
    image_out = temporary(*pimage_out) ;; gets "returned"
    ptr_free, pimage_out

    if ( use_adt_color ) then begin
       orient_circle = circle(150,150,1)
       pcircle = ptr_new(orient_circle, /no_copy)
       mas_rotate_flip, pcircle
    endif else begin
       pcircle = ptr_new()
    endelse

    if (keyword_set(display)) then begin
       mas_display_multi, $
          image_out, $
          tab_title=('DIFF Scan '+strcompress(string(ci+1), /remove_all)+ $
                     '/Slice '  +strcompress(string(sl), /remove_all)), $
          fov_x=xfov, fov_y=yfov, fov_units=1, $
          extra_thumbnail=(ptr_valid(pcircle)) ? temporary(*pcircle) : 0
    endif

end

;; Subroutine name: mdt_make_glyph3d
;; Created by: BT, 2008-06
;; Calling Information:
;;
;; Bugs or Important Comments to Developers:
;;
;;
;; Purpose of subroutine:
;;
;; Creates (a) glyph image plane(s) that can be rotated in xobjview
;;
;; Editing Information:

pro mdt_make_glyph3d, reco_obj, slice_range, $
                      min_fa=min_fa, max_fa=max_fa, $
                      threshold=threshold

    common scan_data
    forward_function mas_roi_get_current_mask
    ci = project.ci

    mas_load_state_1

;    catch, error_status
;    if (error_status ne 0) then begin
;        print, "Error!"
;    endif

    if (not ptr_valid(project.imndArray[ci].bval_Array)) then begin
        void = dialog_message(['B values are not present for this scan.', $
                               'Please make sure that this is a diffusion-',$
                               'weighted scan.'], /error, /center)
        return
    endif

    if not keyword_set(min_fa) then min_fa = 0.0
    if not keyword_set(max_fa) then max_fa = 1.0
    if not keyword_set(threshold) then threshold = 3.

    if (ptr_valid(project.procpramarray[ci].difftools_state)) then begin
        difftools_state = project.procpramarray[ci].difftools_state
    endif else begin
        difftools_state = mdt_make_default_state()
    endelse

    crop_roi     = (*difftools_state).crop_roi
    mask_roi     = (*difftools_state).mask_roi
    vdim         = (*difftools_state).glyph_size
    vmult        = (*difftools_state).glyph_size_mult
    l0_intensity = (*difftools_state).l0_intensity / 100.0
    l1_intensity = (*difftools_state).l1_intensity / 100.0
    l2_intensity = (*difftools_state).l2_intensity / 100.0
    image_opac   = (*difftools_state).image_opacity   / 100.0
    glyph_opac   = (*difftools_state).glyph_opacity   / 100.0
    glyph_shiny  = (*difftools_state).glyph_shiny
    use_adt_color= (*difftools_state).use_adt_color
    ;; incorporate freq, phase, slice viewing direction
    sz_img = size(*project.dataarray[ci].state1, /dimensions)
    adim = project.imndarray[ci].adim

    ;; Set up dimensions and grab the slice data
    case project.procPramArray[project.ci].slice_axis of
        0: begin                ; freq_phase
            xdim = sz_img[0]
            ydim = sz_img[1]
            zdim = sz_img[2]
            xfov = project.imndarray[ci].f_fov
            yfov = project.imndarray[ci].p_fov
            slice_data = fltarr(xdim, ydim, zdim, adim)
            for sl = 0, n_elements(slice_range)-1 do begin
                slice_data[*,*,sl,*] = reform((*project.dataarray[ci].state1)[*,*,slice_range[sl],*])
            endfor
        end
        
        1: begin                ; freq_slice
            xdim = sz_img[0]
            ydim = sz_img[2]
            zdim = sz_img[1]
            xfov = project.imndarray[ci].f_fov
            yfov = project.imndarray[ci].s_fov
            slice_data = fltarr(xdim, ydim, zdim, adim)
            for sl = 0, n_elements(slice_range)-1 do begin
                slice_data[*,*,sl,*] = reform((*project.dataarray[ci].state1)[*,slice_range[sl],*,*])
            endfor
        end
        
        2: begin                ; phase_slice
            xdim = sz_img[1]
            ydim = sz_img[2]
            zdim = sz_img[0]
            xfov = project.imndarray[ci].p_fov
            yfov = project.imndarray[ci].s_fov
            slice_data = fltarr(xdim, ydim, zdim, adim)
            for sl = 0, n_elements(slice_range)-1 do begin
                slice_data[*,*,sl,*] = reform((*project.dataarray[ci].state1)[slice_range[sl],*,*,*])
            endfor
        end
    endcase
    
    xstart = 0
    ystart = 0

    glyph_color    = [50,50,100] 
    glyph_diffuse  = [90,90,110]
    glyph_specular = [150,150,170]

    vdim_effective = vmult*vdim
    
    max_signal = max(*project.dataarray[ci].state1) 

    pct_threshold = max_signal * threshold/1000.0
    print, "threshold value: ", pct_threshold

    osurfs = objarr(xdim,ydim, n_elements(slice_range))
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    mas_tessellator_make_sphere, vdim=vdim_effective, /normalize, $
                                vertlist=vertlist, polylist=polylist
    reco_obj->setReconstructionDirections, vertlist

    progressBar = obj_new('progressbar', title='Computing Probabilities', $
                          text='Computing Probabilities...', /fast_loop)
    progressBar->start
    iter = 0
    cols = indgen(255)
    junk = get_orient_rgb_axis(tx_matrix=tx)
    if (obj_class(reco_obj) ne 'MAS_MOW_BULK_RECONSTRUCTOR') then begin
       tx = *project.imndarray[project.ci].acq_matrix # tx
    endif
    vertlist = vert_t3d(vertlist, matrix=tx)

    syst = systime(1)
    
    if (mask_roi) then begin
       mask = mas_roi_get_current_mask([xdim, ydim], crop=crop_dims, /no_transform)
    endif else begin
       mask = ptr_new()
    endelse

    for sl = 0, n_elements(slice_range)-1 do begin
        
        glyph_color=abs([255-8*cols[2*sl], sl*20, 10*cols[3*((sl mod 2)+1)]])
        for x = 0, xdim-1 do begin

            for y = 0, ydim-1 do begin

                if (ptr_valid(mask)) then begin
                    if ((*mask)[x,y] eq 0) then continue
                endif

                voxel_data = reform(slice_data[x,y,sl,*])

                if (max(voxel_data) lt pct_threshold) then continue

                fit = reco_obj->getDisplacementProbability(data=voxel_data)
                if ( total(finite(fit, /nan) + total(finite(fit, /inf)) ) gt 0) then continue

                max_fit = max(fit, min=min_fit)
                fit = (fit - min_fit)/(max_fit-min_fit)
                fit = (vdim_effective * 0.48) * fit

                verts = vertlist
                verts[0,*] = verts[0,*]*fit + (x*vdim_effective  + vdim_effective/2.)
                verts[1,*] = verts[1,*]*fit + (y*vdim_effective  + vdim_effective/2.)
                verts[2,*] = verts[2,*]*fit + (slice_range[sl]*vdim_effective + vdim_effective/2.)
                
                if (keyword_set(use_adt_color) and project.procpramarray[ci].adt_proccess_flag) then begin
                    glyph_color = adt_get_orientation_color(x+xstart, y+ystart, slice_range[sl]) 
                endif
                
                tmp = obj_new('idlgrpolygon', verts, $
                              polygons=polylist,  $
                              shading=1,  $
                              diffuse=glyph_diffuse, $
                              color=glyph_color, $
                              specular=glyph_specular,$
                              alpha_channel=glyph_opac, $
                              reject=0, $
                              shininess=glyph_shiny, $
                              depth_test_disable=2)
                
                osurfs[x,y,sl] = tmp
                
            endfor
            
            progressbar->update, (100.0) * x/xdim
            if (progressBar->checkcancel()) then begin
                progressBar->destroy
                if (ptr_valid(mask)) then ptr_free, mask
                obj_destroy, osurfs
                return
            endif

        endfor
        
    endfor

    progressBar->destroy

    print, '... Done. (Time = '+string(systime(1)-syst, format='(F0.5)')+')'
    syst = systime(1)
    
    sf = obj_valid(osurfs)
    ptr_free, mask

    xobjview, osurfs[where(sf eq 1)], /block, background=[0,0,0]
    obj_destroy, osurfs
    return

end

pro mas_diffusion_tools

end
