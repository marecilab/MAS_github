;; $id$
;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; Info will be put here that describes how the display
; handler works.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Subroutine name: mmd_init
; Created by: BT, 2008-06
; Calling Information:
; simple       - IN - if this keyword is set, the resulting window
;                will be a standalone, untabbed window index
; tlb_title    - IN - the title for the tlb's title bar
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; RETURNS: the widget id of the top level base
; Editing Information:

function mmd_init, $
                   simple=simple,  $
                   xsize=xsize,    $
                   ysize=ysize,    $
                   tlb_title=tlb_title, $
                   show_line_data_btn=show_line_data_btn, $
                   show_mini_window=show_mini_window

    compile_opt hidden
    common scan_data
    common common_widgets

    if (keyword_set(simple)) then simple = 1 else simple = 0
    
    if (not keyword_set(tlb_title)) then tlb_title='MAS Multi Display'

    tlb = widget_base(title=tlb_title, $
                      group_leader=WID_BASE_MAIN, $
                      /row, $
                      /tlb_size_events, $
                      xpad=0, $
                      ypad=0, event_pro='mas_multi_display_event')

    ; Create the base for the context menu.
    contextBase = widget_base(tlb, /context_menu, UNAME = 'drawContext')

    ; Create the buttons for the context menu.
    loadCTButton = widget_button(contextBase, value='Load Color Table...', $
                                 uname='ctm_load_color_table', $
                                 event_pro='mmd_draw_context_menu_event')
    loadCTButton = widget_button(contextBase, value='Show/Update Color Bar...', $
                                 uname='ctm_show_color_bar', $
                                 event_pro='mmd_draw_context_menu_event')

    toolbar_base  = widget_base(tlb, /toolbar, /column,  xpad=0, ypad=0)

    if (not simple) then begin

        action_base   = widget_base(toolbar_base, /column, frame=0,$
                                    xpad=0, ypad=0, /align_center)
        btn_save      = widget_button(action_base, value="Save Tab", $
                                      event_pro='mmd_save_current_tab')
        btn_close     = widget_button(action_base, value="Close Tab", $
                                      event_pro='mmd_delete_current_tab')
        btn_rename    = widget_button(action_base, value="Rename Tab",$
                                      event_pro='mmd_rename_current_tab')
        btn_detach    = widget_button(action_base, value="Detach Tab",$
                                      event_pro='mmd_detach_current_tab')
        btn_animate   = widget_button(action_base, value="Animate Tabs",$
                                      event_pro='mmd_animate')
        
        tool_bar      = widget_base(toolbar_base, /column, /toolbar, /exclusive, $
                                    /frame, /base_align_center, /align_center,$
                                    xpad=0, ypad=0)
        btn_zoom      = widget_button(tool_bar, value="Zoom", uname='btn_zoom',$
                                      event_pro='mmd_tool_event', $
                                      tooltip='Left Btn: In; Right Btn: Out')
        btn_scroll    = widget_button(tool_bar, value="Scroll", uname='btn_scroll',$
                                      event_pro='mmd_tool_event')
        btn_intensity = widget_button(tool_bar, value="Intensity", uname='btn_intensity',$
                                      event_pro='mmd_tool_event', $
                                      tooltip='Mouse: Up - Increase, Dn - Decrease')
        btn_distance  = widget_button(tool_bar, value="Distance", uname='btn_distance',$
                                      event_pro='mmd_tool_event')

        readout_base  = widget_base(toolbar_base, /column, /frame, xpad=0, ypad=0, /align_left)
        lbl_color     = widget_label(readout_base, value='FP: xxxxxxxxxxxx', /align_left)
        lbl_coord     = widget_label(readout_base, value='Pos: (xxxx,yyyyy)', /align_left)
        lbl_zoom      = widget_label(readout_base, value='Zoom: 100%', /align_left)
        lbl_dist      = widget_label(readout_base, value='Dist: ------- px', /align_left)

        mini_base     = widget_base(toolbar_base, /column, xpad=0, ypad=0, /align_center)
        draw_mini     = widget_draw(toolbar_base, /align_center, $
                                    uname='draw_mini', $
                                    graphics_level=2, $
                                    renderer=1, $
                                    event_pro='mmd_mini_window_event', $
                                    /button_events, /motion_events)

        extra_base    = widget_base(toolbar_base, /column, xpad=0, ypad=0, /align_center)
        
        tab_base      = widget_tab(tlb)
        widget_control, btn_scroll, /set_button
    
    endif else begin

        action_base   = widget_base(toolbar_base, /column, frame=0, xpad=0, ypad=0, /align_center)
        btn_save      = widget_button(action_base, value="Save Image",$
                                      event_pro='mmd_save_current_tab')

        tool_bar      = widget_base(toolbar_base, /column, /toolbar, /exclusive,$
                                    frame=1, /align_center, /base_align_center, xpad=0, ypad=0)

        btn_zoom      = widget_button(tool_bar, value="Zoom", uname='btn_zoom', $
                                      event_pro='mmd_tool_event', $
                                      tooltip='Left Btn: In; Right Btn: Out')
        btn_scroll    = widget_button(tool_bar, value="Scroll", uname='btn_scroll', $
                                      event_pro='mmd_tool_event')
        btn_intensity = widget_button(tool_bar, value="Intensity", uname='btn_intensity',$
                                      event_pro='mmd_tool_event')
        btn_distance  = widget_button(tool_bar, value="Distance", uname='btn_distance',$
                                      event_pro='mmd_tool_event')
                                      
        if (keyword_set(show_line_data_btn)) then begin
            btn_line_data = widget_button(tool_bar, value="Line Data", uname='btn_line_data',$
                                          event_pro='mmd_tool_event')
        endif
        
        readout_base  = widget_base(toolbar_base, /column, /frame, xpad=0, ypad=0, /align_left)
        lbl_color     = widget_label(readout_base, value='FP: xxxxxxxxxxxx', /align_left)
        lbl_coord     = widget_label(readout_base, value='Pos: (xxxx,yyyyy)', /align_left)
        lbl_zoom      = widget_label(readout_base, value='Zoom: 100%', /align_left)
        lbl_dist      = widget_label(readout_base, value='Dist: ------- px', /align_left)
        
        if (keyword_set(show_mini_window)) then begin
            mini_base     = widget_base(toolbar_base, /column, xpad=0, ypad=0, /align_center)
            draw_mini     = widget_draw(toolbar_base, /align_center, $
                                    uname='draw_mini', $
                                    graphics_level=2, $
                                    renderer=1, $
                                    event_pro='mmd_mini_window_event', $
                                    /button_events, /motion_events)
        endif else begin
            draw_mini     = -1
        endelse

        extra_base    = widget_base(toolbar_base, /column, xpad=0, ypad=0, /align_center)

        tab_base      = widget_base(tlb, xpad=0, ypad=0)
        widget_control, btn_scroll, /set_button
        widget_control, toolbar_base

    endelse

;            ----------------
    str = [ '     .....      ', $ ;; 1 | Poor attempt at
            '    .#####.     ', $ ;; 2 | creating a magnifying
            '   .##...##.    ', $ ;; 3 | glass style cursor
            '  .##.   .##.   ', $ ;; 4 
            ' .##.     .##.  ', $ ;; 5 
            ' .##.      .##. ', $ ;; 6 
            '  .##. $    .##.', $ ;; 7 
            '   .###.    .##.', $ ;; 8
            '  .#####.  .##. ', $ ;; 9
            ' .###..##..##.  ', $ ;; 10
            '.###.  .####.   ', $ ;; 11
            '.##.    ....    ', $ ;; 12
            '##.             ', $ ;; 13
            '..              ', $ ;; 14
            '                ', $ ;; 15
            '                ' ]  ;; 16
;            1234567890123456

    cur = create_cursor(str, hotspot=hotspot, mask=mask)
    register_cursor, 'ZOOM', cur, hotspot=hotspot, mask=mask

    if (not keyword_set(xsize)) then xsize = 0
    if (not keyword_set(ysize)) then ysize = 0

    topstate = ptr_new({ tlb:tlb, $
                         tab_base:tab_base, $
                         tabbed:(simple) ? 0 : 1, $
                         extra_base:extra_base, $
                         extra_base_filled:0, $
                         tab_counter: 0, $
                         lbl_color:lbl_color, $
                         lbl_coord:lbl_coord, $
                         lbl_zoom:lbl_zoom, $
                         lbl_dist:lbl_dist, $
                         draw_mini:(widget_info(draw_mini, /valid_id)) ? draw_mini : -1, $
                         draw_mini_click_hold: 0, $
                         mini_window:obj_new(), $
                         cursor:'MOVE', $
                         view_xsize: xsize, $
                         view_ysize: ysize, $
                         toolbar_xsize: 0L, $
                         current_tool:'SCROLL'}, /no_copy)
    
    widget_control, tlb, set_uvalue=topstate
    widget_control, tlb, /realize

    tb_geom = widget_info(toolbar_base, /geometry)
    (*topstate).toolbar_xsize = tb_geom.xsize + tb_geom.space
    
    xmanager, 'mas_multi_display', tlb, cleanup='mmd_cleanup', /no_block
    
    return, tlb

end

; Subroutine name: mmd_add_tab
; Created by: BT, 2008-06
; Calling Information:
;
; topstate     - IN - the state structure associated with the top
;                     level base
;
; image_data   - IN - the image to be displayed in the draw tab.
; title        - IN - the tab title
; fov_in       - IN - field of view: [xfov, yfov]
; fov_units    - IN - units of the fov: 0 = px, 1 = cm, 2 = mm
; fp_vals      - IN - pointer to array containing the floating point
;                     values to display at (x,y)
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Creates a new tab in a tabbed display widget. Sets up all parameters.
; Editing Information:

pro mmd_add_tab, $
                 topstate,     $
                 image_data,   $
                 title=title,  $
                 fov_in=fov_in, $
                 fov_units=fov_units, $
                 fp_vals=fp_vals, $
                 fp_range=fp_range, $
                 show_axis_labels=show_axis_labels, $
                 axis_labels=axis_labels, $
                 drawstate=drawstate

    compile_opt hidden
    scr_size = get_screen_size()

    imsz = size(image_data, /dimensions)
    
    if (n_elements(fp_range) ne 2) then begin
        float_range = ptr_new()
    endif else begin
        float_range = ptr_new(fp_range)
    endelse
    
    if (not ptr_valid(fp_vals)) then begin
        fp_vals = ptr_new()
        float_range = ptr_new([0, 255])
    endif

    if ((*topstate).view_xsize ne 0 and (*topstate).view_ysize ne 0) then begin
        
        xsize = (*topstate).view_xsize
        ysize = (*topstate).view_ysize

    endif else begin

        xsize = min([imsz[1], scr_size[0]*0.8])  > 200
        ysize = min([imsz[2], scr_size[1]*0.75]) > 220

        (*topstate).view_xsize = xsize
        (*topstate).view_ysize = ysize

    endelse

    tab_base  = (*topstate).tab_base

    ;; keep track of number of tabs.
    (*topstate).tab_counter++

    if (not keyword_set(title)) then begin
        title = strcompress("Untitled "+ string((*topstate).tab_counter))
    endif

    draw_base = widget_base(tab_base, title=title, xpad=0, ypad=0)

    main_draw = widget_draw(draw_base, $
                            retain=0,  $
                            xsize=xsize, ysize=ysize, $
                            /button_events, $
                            /expose_events, $
                            /motion_events, $
                            keyboard_events=2,$
                            renderer=1, $
                            graphics_level=2, $
                            kill_notify='mmd_tab_cleanup', $
                            event_pro='mmd_image_event')

    widget_control, main_draw, /realize
    widget_control, main_draw, get_value=owindow

    ;; start with the image centered in the viewplane
    viewplane_rect = float([imsz[1]/2 - xsize/2, imsz[2]/2 - ysize/2, xsize, ysize])
    oview  = obj_new('idlgrview', viewplane_rect=viewplane_rect, color=[128,128,128])
    omodel = obj_new('idlgrmodel')

    show_axis_labels = keyword_set(show_axis_labels) ? 1 : 0
    if (n_elements(axis_labels) ne 4) then begin
       
       axis_labels = [' ',' ',' ',' ']

    endif
    
    oaxis_labels = obj_new('idlgrtext', strings=axis_labels, color=[255,0,0], /onglass, $
                           hide = (1-show_axis_labels), $
                           locations=[ [imsz[1]/2,imsz[2]-10], $
                                       [imsz[1]-10, imsz[2]/2], $
                                       [imsz[1]/2, 1], $
                                       [1, imsz[2]/2] ] )
    
    oimage = obj_new('idlgrimage', image_data, interleave=0)

    ;; used to draw distance measure line and zoom box
    odist      = obj_new('idlgrpolyline', color=[0,255,0], alpha_channel=.75, thick=2)
    ozoom_box  = obj_new('idlgrpolyline', color=[255,0,0], alpha_channel=.75, thick=2)

    ;; used to draw the orthogonal crosshairs for line data profile plotting
    oline_data_h = obj_new('idlgrpolyline', color=[255,0,0])
    oline_data_v = obj_new('idlgrpolyline', color=[0,0,255])

    omodel->add, oimage    
    omodel->add, [odist, oline_data_h, oline_data_v, ozoom_box, oaxis_labels]
    
    oview->add, omodel

    owindow->setCurrentCursor, (*topstate).cursor

    ;; Used to draw the mini, comprehensive window.
    draw_mini = (*topstate).draw_mini
    
    ;; make sure mini-window is present
    if (widget_info(draw_mini, /valid_id) ne 0) then begin

        ;; resize mini window for aspect ratio of image
        geom = widget_info(draw_mini, /geometry)
        
        if (imsz[1] gt imsz[2]) then begin

            resz = float(imsz[2])/float(imsz[1])
            ns = fix(100.*resz)
            widget_control, draw_mini, draw_xsize=100, draw_ysize=100*resz

        endif else begin

            resz = float(imsz[1])/float(imsz[2])
            ns = fix(100.*resz)
            widget_control, draw_mini, draw_ysize=100, draw_xsize=100*resz

        endelse

        ;; set up the mini window
        if not obj_valid( (*topstate).mini_window ) then begin
            widget_control, draw_mini, get_value=mini_window
            if (not obj_valid(mini_window)) then begin
                message, 'Can''t create mini window!'
            endif
            (*topstate).mini_window = mini_window
        endif else begin
            mini_window = (*topstate).mini_window
        endelse

        mini_view  = obj_new('idlgrview', viewplane_rect=[0.,0., imsz[1], imsz[2]])
        mini_model = obj_new('idlgrmodel')
        mini_model->add, oimage, /alias
        bl = [viewplane_rect[0], viewplane_rect[1]]
        tr = [viewplane_rect[0]+viewplane_rect[2], viewplane_rect[1]+viewplane_rect[3]]
 
        box = [ [ bl[0], bl[1] ], $
                [ tr[0], bl[1] ], $
                [ tr[0], tr[1] ], $
                [ bl[0], tr[1] ], $
                [ bl[0], bl[1] ] ]

        mini_box = obj_new('idlgrpolyline', box, color=[255,0,0])
        mini_model->add, mini_box
        mini_view->add, mini_model
        mini_window->draw, mini_view

    endif else begin

        mini_view   = obj_new()
        mini_model  = obj_new()
        mini_box    = obj_new()
        mini_window = obj_new()

    endelse

    ;; we will store an extra copy for the purposes of returning the
    ;; floating point value
    tmp = image_data
    mx = max(tmp, min=mn)

    loadct, 0, /silent, rgb_table=color_table
    
    drawstate = ptr_new({ oview:oview,   $
                          omodel:omodel, $
                          oimage:oimage, $
                          owindow:owindow, $
                          $
                          odist: odist, $
                          ozoom_box:ozoom_box, $
                          oaxis_labels: oaxis_labels, $
                          show_axis_labels:show_axis_labels, $
                          line_data_gui_win: -1, $
                          oline_data_h: oline_data_h, $
                          oline_data_v: oline_data_v, $
                          $
                          tab_title:title, $
                          vpr:viewplane_rect, $
                          mini_view: mini_view, $
                          mini_model: mini_model, $
                          mini_box: mini_box, $
                          mini_window: mini_window, $
                          zoom_factor:.1, $ ;; see mmd_image event
                          zoom_level:0,  $  ;; see mmd_image event
                          zoom_amount:1.0, $
                          draw_dims:[xsize, ysize], $
                          img_effective_size:[imsz[1],imsz[2]], $
                          img_actual_size:[imsz[1],imsz[2]], $
                          fov:fov_in,       $
                          fov_units:fov_units, $
                          x_pos: 0,      $  ;; for scrolling
                          y_pos: 0,      $
                          down_x:0,      $ ;; down_{x,y} used for click-hold
                          down_y:0,      $  ;; dragging events
                          is_down:0,     $ ;; track if user is dragging
                          float_vals:fp_vals,$
                          float_range:float_range, $
                          image_data_max:long(mx), $
                          image_it_max:long(mx), $ 
                          image_it_min:long(mn), $ 
                          image_it_top:long(255), $
                          locked_with:ptr_new(), $
                          extra_thumbnail:ptr_new(), $
                          image_color_table:ptr_new(color_table, /no_copy), $
                          image_data:ptr_new(tmp, /no_copy) }, $
                        /no_copy)
    widget_control, $
      main_draw,    $
      set_uvalue = drawstate

    mmd_update_axis_label_position, drawstate
    owindow->draw, oview

    ;; make this new tab auto-active
    widget_control, tab_base, set_tab_current=widget_info(tab_base, /tab_number)-1

end


; Subroutine name: mmd_get_tab_id
; Created by: BT, 2008-06
; Calling Information:
; topstate     - IN - the state structure for the top display base
;  active      - IN - set this keyword to return the id of the
;                     currently active tab_id
; tab_index    - IN - set this to get the id of a specific tab by its
;                     index
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; RETURNS: the widget id of the tab.
; Editing Information:

function mmd_get_tab_id, topstate, active=active, tab_index=tab_index

    compile_opt hidden

    tab_base = (*topstate).tab_base

    if ((*topstate).tabbed eq 0) then begin
        return, widget_info(tab_base, /child)
    endif

    if (keyword_set(active) or not arg_present(tab_index)) then begin
        tab_index = widget_info(tab_base, /tab_current)
    endif

    all_tabs = widget_info(tab_base, /all_children)

    if (n_elements(all_tabs) le tab_index) then return, -1
    
    return, all_tabs[tab_index]
    
end

; Subroutine name: mmd_get_pixel_value
; Created by: BT, 2008-06
; Calling Information:
;
;   drawstate - IN - the state structure associates with a draw widget
;   coord     - IN - the (x,y) location relative to the
;                    current window view. Generally an
;                    event.{x,y} 
;  translated - OUT - set this to a named variable to receive the
;                     translated coordinates relative to the actual image
;                     being displayed. Accounts for zoom.
;  float_val  - OUT - set this to a named variable to receive the
;                     floating point value at a particular image
;                     coordinate. If no floating point values were
;                     set when this draw widget was initted, -1 is
;                     returned
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; RETURNS: image pixel value at a coordinate
; Editing Information:

function mmd_get_pixel_value, drawstate, coord, translated=translated, float_val=float_val

    zoom_f = (*drawstate).zoom_amount

    vpr = (*drawstate).vpr
    imsize = size((*(*drawstate).image_data), /dimensions)

    new_x = round(vpr[0] + (coord[0])*(1.0/zoom_f)) < (imsize[1]-1)
    new_y = round(vpr[1] + (coord[1])*(1.0/zoom_f)) < (imsize[2]-1)
    new_x = max([0, new_x])
    new_y = max([0, new_y])

    translated = [new_x, new_y]
    
    if (ptr_valid( (*drawstate).float_vals )) then begin
        float_val = double( (*(*drawstate).float_vals)[translated[0], translated[1]])
    endif else begin
        float_val = -1
    endelse

    return, long((*(*drawstate).image_data)[*,new_x, new_y])

end

; Subroutine name: mmd_update_axis_label_position
; Created by: BT, 2008-06
; Calling Information:
;
;   drawstate - IN - the state structure associates with a draw widget

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; 
;  Updates the position of the red orientation letters when the image
;  is resized, scrolled, or zoomed
;
; Editing Information:

pro mmd_update_axis_label_position, drawstate

    if (not ptr_valid(drawstate)) then return
    if (*drawstate).show_axis_labels eq 0 then return

    vpr  = (*drawstate).vpr
    zoom = (*drawstate).zoom_amount
    imsz = (*drawstate).img_actual_size

    top   = [ vpr[0]+vpr[2]/2. , vpr[1]+vpr[3]-12./zoom ]
    right = [ vpr[0]+vpr[2]-12./zoom , vpr[1]+vpr[3]/2. ]
    bot   = [ vpr[0]+vpr[2]/2.,  vpr[1] ]
    left  = [ vpr[0]+3./zoom, vpr[1]+vpr[3]/2. ]
    
    loc = [ [ top ], [ right ], [ bot ], [ left ] ]

    (*drawstate).oaxis_labels->setProperty, locations = loc

end

; Subroutine name: mmd_update_mini_window
; Created by: BT, 2008-06
; Calling Information:
;
;   drawstate  - IN - state structure for top-level base
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;  Updates the mini window with viewplane rect and image
;
; Editing Information:

pro mmd_update_mini_window, drawstate

    if (not ptr_valid(drawstate)) then return
    if (not obj_valid((*drawstate).mini_window)) then return
    
    vpr = (*drawstate).vpr
    
    ;; bl = bottom-left, tr = top-right
    bl = [vpr[0], vpr[1]]
    tr = [vpr[0]+vpr[2], vpr[1]+vpr[3]]
    
    box = [ [ bl[0], bl[1] ], $
            [ tr[0], bl[1] ], $
            [ tr[0], tr[1] ], $
            [ bl[0], tr[1] ], $
            [ bl[0], bl[1] ] ]
    
    (*drawstate).mini_box->setProperty, data=box
    (*drawstate).mini_window->draw, (*drawstate).mini_view

end

; Subroutine name: mmd_update_readout
; Created by: BT, 2008-06
; Calling Information:
;
;   topstate  - IN - state structure for top-level base
;      zoom   - IN - zoom percent
;     fpval   - IN - floating point value at pixel location
;    dist     - IN - string - distance with units
;   {x,y}pos  - IN - pixel location within image
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;  Updates the information label widgets with info specified in
;  arguments
;
; Editing Information:

pro mmd_update_readout, topstate, zoom=zoom, dist=dist, fpval=fpval, xpos=xpos, ypos=ypos

    compile_opt hidden

    if (n_elements(xpos) and n_elements(ypos)) then begin
        widget_control, (*topstate).lbl_coord, $
          set_value='Pos: ('+strcompress(string(xpos)+','+string(ypos)+')', /remove_all)
    endif

    if (arg_present(fpval)) then begin
        if (n_elements(fpval) ne 0) then begin
            fp_string = string(fpval, format='(F0.5)')
        endif else begin
            fp_string = 'N/A'
        endelse
        widget_control, (*topstate).lbl_color, $
          set_value='FP: '+strcompress(fp_string, /remove_all)
    endif

    if (keyword_set(zoom)) then begin
        widget_control, (*topstate).lbl_zoom, $
          set_value='Zoom: '+strcompress(string(zoom), /remove_all)+'%'
    endif

    if (keyword_set(dist)) then begin
        widget_control, (*topstate).lbl_dist, $
          set_value=strcompress('Dist: '+string(dist))
    endif

end

; Subroutine name: mmd_do_scroll
; Created by: BT, 2008-06
; Calling Information:
;
;   drawstate - IN - the state structure associates with a draw widget
;      event  - IN - the original event (re-thrown to us)
;   center_on - IN - set this keyword to have the image centered in
;                    the viewport
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;  Handles the fine details of scrolling or moving an image around in
;  the viewport. Generally this procedure will be called from from
;  other event handlers to do the actuall scrolling.
;
; Editing Information:

pro mmd_do_scroll, drawstate, event, center_on=center_on

    ;; scroll sensitivity based on zoom level 
    scroll_sense = (*drawstate).zoom_amount

    dx = (*drawstate).down_x - event.x
    dy = (*drawstate).down_y - event.y
    
    vpr = (*drawstate).vpr

    if (keyword_set(center_on)) then begin
        vpr[0] = (event.x - (0.5 * vpr[2])) > 0.
        vpr[1] = (event.y - (0.5 * vpr[3])) > 0.
    endif else begin
        ;; bound the scrolling to the image dimensions
        ;; bottom-left corner
        vpr[0] = (vpr[0] + dx/scroll_sense) > 0.
        vpr[1] = (vpr[1] + dy/scroll_sense) > 0.
    endelse

    ;; top-right corner
    if vpr[2] gt (*drawstate).img_actual_size[0] then begin
        vpr[0] = ((*drawstate).img_actual_size[0] - vpr[2])/2.
    endif else if vpr[0] + vpr[2] gt (*drawstate).img_actual_size[0] then begin
        vpr[0] = ((*drawstate).img_actual_size[0] - vpr[2])
    endif
    
    if vpr[3] gt (*drawstate).img_actual_size[1] then begin
        vpr[1] = ((*drawstate).img_actual_size[1] - vpr[3])/2.
    endif else if vpr[1] + vpr[3] gt (*drawstate).img_actual_size[1] then begin
        vpr[1] = ((*drawstate).img_actual_size[1] - vpr[3])
    endif
    
    ;; update the new values, clear any distance lines, update
    ;; the state structure with the current 'down' location.
    (*drawstate).vpr = vpr
    (*drawstate).oview->setProperty, viewplane_rect=(*drawstate).vpr
    ;(*drawstate).odist->setProperty, data=[[0,0,0],[0,0,0]], hide=1
    ;(*drawstate).oline_data_h->setProperty, data=[[0,0],[0,0]], hide=1
    ;(*drawstate).oline_data_v->setProperty, data=[[0,0],[0,0]], hide=1
    mmd_update_axis_label_position, drawstate

    (*drawstate).owindow->draw, (*drawstate).oview
    (*drawstate).down_x = event.x
    (*drawstate).down_y = event.y

    mmd_update_mini_window, drawstate

end

; Subroutine name: mmd_do_zoom
; Created by: BT, 2008-06
; Calling Information:
;
;   drawstate - IN - the state structure associates with a draw widget
;      event  - IN - the original event (re-thrown to us)
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
;  Handles the fine details of rezizing an image in the viewport
;
; Editing Information:

pro mmd_do_zoom, drawstate, event

    ;; double-click
    if (event.clicks eq 2) then begin
        if (event.press ne 1) then return
        vpr = (*drawstate).vpr
        vpr[2] = (*drawstate).draw_dims[0]
        vpr[3] = (*drawstate).draw_dims[1]
        if vpr[0] + vpr[2] gt (*drawstate).img_actual_size[0] then begin
            vpr[0] = ((*drawstate).img_actual_size[0] - vpr[2])/2.
        endif
        
        if vpr[1] + vpr[3] gt (*drawstate).img_actual_size[1] then begin
            vpr[1] = ((*drawstate).img_actual_size[1] - vpr[3])/2.
        endif
        
        (*drawstate).zoom_amount = 1.0
        (*drawstate).vpr = vpr
        
    endif else if (*drawstate).is_down then begin
        
        if abs(((*drawstate).down_x-event.x)*((*drawstate).down_y-event.y)) lt 9 then return
        
        junk = mmd_get_pixel_value(drawstate, [(*drawstate).down_x, (*drawstate).down_y], translated=tr_down)
        junk = mmd_get_pixel_value(drawstate, [event.x, event.y], translated=tr_event)
        
        ;; bl = bottom-left, tr = top-right
        bl = [ min([tr_event[0], tr_down[0]]) , min([tr_event[1], tr_down[1]]) ]
        tr = [ max([tr_event[0], tr_down[0]]) , max([tr_event[1], tr_down[1]]) ]
        
        sq = [ [ bl[0], bl[1] ], $
               [ tr[0], bl[1] ], $
               [ tr[0], tr[1] ], $
               [ bl[0], tr[1] ], $
               [ bl[0], bl[1] ] ]
        
        (*drawstate).ozoom_box->setProperty, data=sq, hide=0
        (*drawstate).owindow->draw, (*drawstate).oview
        mmd_update_mini_window, drawstate
        return
        
    endif else if (event.press eq 4) then begin
        
        ;; Zoom Out
        zfact = (*drawstate).zoom_factor
        
        vpr = (*drawstate).vpr
        
        ;; changing the viewplane rect is sufficient.
        vpr[2:3] *= (1 + zfact)
        (*drawstate).img_effective_size /= (1 + zfact)
        
        ;; we want this thing centered if the window is larger
        ;; than the zoomed image.
        
        vpr[0] = max([0,vpr[0] + event.x - (*drawstate).draw_dims[0]/2])
        vpr[1] = max([0,vpr[1] + event.y - (*drawstate).draw_dims[1]/2])
        
        if vpr[0] + vpr[2] gt (*drawstate).img_actual_size[0] then begin
            vpr[0] = ((*drawstate).img_actual_size[0] - vpr[2])/2.
        endif
        
        if vpr[1] + vpr[3] gt (*drawstate).img_actual_size[1] then begin
            vpr[1] = ((*drawstate).img_actual_size[1] - vpr[3])/2.
        endif
        
        ;; update data structure and clear distance lines (if any)
        (*drawstate).odist->setProperty, hide=1
        (*drawstate).vpr = vpr
        (*drawstate).zoom_amount = float((*drawstate).draw_dims[1]/vpr[3])
        
    endif else begin
        
        (*drawstate).ozoom_box->getProperty, data=sq, hide=hi
        
        if (n_elements(sq) eq 0) then return
        if (hi ne 0)             then return
        
        vpr = (*drawstate).vpr
        new_width = [ sq[4]-sq[0], sq[5]-sq[1] ]
        
        new_width[1] = new_width[0]*vpr[3]/vpr[2]
        
        zoom_amt = float((*drawstate).draw_dims[1]/new_width[1])
        
        (*drawstate).zoom_amount = zoom_amt
        (*drawstate).vpr = [ sq[0], sq[1], new_width[0], new_width[1] ]
        
    endelse   
    
    (*drawstate).oview->setProperty, viewplane_rect=(*drawstate).vpr
    (*drawstate).ozoom_box->setProperty, hide=1
    mmd_update_axis_label_position, drawstate
    (*drawstate).owindow->draw, (*drawstate).oview
    mmd_update_mini_window, drawstate
end

; Subroutine name: mmd_redraw_tab
; Created by: BT, 2008-06
; Calling Information:
;
;  topstate  - pointer to top level base state structure
;  tab_index - index of tab to redraw
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; redraws the tab at index tab_index, updates the "mini window" if it
; exists 
;
; Editing Information:

pro mmd_redraw_tab, topstate, tab_index

    compile_opt hidden
    common scan_data
    
    tab_id = mmd_get_tab_id(topstate, tab_index=tab_index)

    if (not widget_info(tab_id, /valid_id)) then return
    
    widget_control, widget_info(tab_id,/child), get_uvalue=drawstate

    imsz = (*drawstate).img_actual_size

    draw_mini = (*topstate).draw_mini
        
    if (imsz[0] gt imsz[1]) then begin
        
        resz = float(imsz[1])/float(imsz[0])
        ns = fix(100.*resz)
        widget_control, draw_mini, draw_xsize=100, draw_ysize=100*resz
        
    endif else begin
        
        resz = float(imsz[0])/float(imsz[1])
        ns = fix(100.*resz)
        widget_control, draw_mini, draw_ysize=100, draw_xsize=100*resz
        
    endelse
    
    if (not ptr_valid(drawstate)) then return

    (*drawstate).owindow->draw, (*drawstate).oview

    mmd_update_readout, topstate, $
      zoom=(*drawstate).zoom_amount*100., $
      dist='---'

    if (obj_valid((*drawstate).mini_window)) then begin
        (*drawstate).mini_window->draw, (*drawstate).mini_view
    endif
    
    if (ptr_valid( (*drawstate).extra_thumbnail ) ) then begin
        widget_control, (*topstate).extra_base, /show
        wd = widget_info ((*topstate).extra_base, /child)
        if (widget_info(wd, /valid_id)) then begin
            widget_control, wd, get_value=wind, /show
            wset, wind
            tv, *((*drawstate).extra_thumbnail), true=3
        endif
    endif else begin
        wd = widget_info ((*topstate).extra_base, /child)
        if (widget_info(wd, /valid_id)) then begin
            widget_control, wd, get_value=wind, /show
            wset, wind
            erase
            
        endif
    endelse

end

; Subroutine name: mmd_detach_cuttent_tab
; Created by: BT, 2008-06
; Calling Information:
;
; event - the event structure
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Event handler for the "Detach Tab" button. This will take the
; information in the current tab and redisplay it in a "standalone"
; window. The current tab will be destroyed
;
; Editing Information:

pro mmd_detach_current_tab, event

    compile_opt hidden

    widget_control, event.top, get_uvalue=topstate
    tab_base = (*topstate).tab_base
    
    active_tab_id = mmd_get_tab_id(topstate, /active)
    widget_control, widget_info(active_tab_id, /child), get_uvalue=drawstate

    new_tlb = mmd_init(/simple, xsize=600, ysize=400,  $
                       tlb_title=(*drawstate).tab_title,$
                       /show_mini_window)
    widget_control, new_tlb, get_uvalue=new_topstate
    mmd_add_tab, new_topstate,  *((*drawstate).image_data), $
      title=(*drawstate).tab_title, $
      fov_in=(*drawstate).fov, $
      fov_units=(*drawstate).fov_units, $
      fp_vals=(*drawstate).float_vals, $
      drawstate=new_drawstate

    if ptr_valid((*drawstate).extra_thumbnail) then begin

        ex = (*(*drawstate).extra_thumbnail)
        wd = widget_info((*new_topstate).extra_base, /child)
        
        if (not widget_info(wd, /valid_id)) then begin
            wd = widget_draw((*new_topstate).extra_base, /align_center, $
                             graphics_level=1, retain=2)
            widget_control, (*new_topstate).tlb, /realize
        endif
        
        geom = widget_info(wd, /geometry)
        
        widget_control, wd, get_value=wind
        wset, wind
        
        tv, ex, true=3 

        (*new_drawstate).extra_thumbnail = ptr_new(ex , /no_copy)
        
    endif

    mmd_delete_current_tab, event

end

; Subroutine name: mmd_rename_current_tab
; Created by: BT, 2008-06
; Calling Information:
;
;  event - event to be handled
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; handles the event of "Rename Tab" button
;
; Editing Information:

pro mmd_rename_current_tab, event

    compile_opt hidden
    widget_control, event.top, get_uvalue=topstate
    tab_base = (*topstate).tab_base
    
    active_tab_id = mmd_get_tab_id(topstate, /active)

    widget_control, widget_info(active_tab_id, /child), get_uvalue=drawstate

    if (not ptr_valid(drawstate)) then return

    newname = mas_get_text_dialog(event.top, $
                                  initial_text=(*drawstate).tab_title, $
                                  message_text='Please enter a new name for this tab:', $
                                  text_size=25)
    
    if (newname eq '') then return

    (*drawstate).tab_title = newname
    widget_control, active_tab_id, base_set_title=newname

end

; Subroutine name: mmd_delete_current_tab
; Created by: BT, 2008-06
; Calling Information:
;
;  event - event to be handled
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; handles the event of "Close Tab" button
;
; Editing Information:

pro mmd_delete_current_tab, event

    compile_opt hidden
    widget_control, event.top, get_uvalue=topstate
    tab_base = (*topstate).tab_base
    
    active_tab_id = mmd_get_tab_id(topstate, /active)
    widget_control, active_tab_id, /destroy
    
    if (widget_info(tab_base, /tab_number) eq 0) then begin
        widget_control, event.top, /destroy
    endif else begin
        curr_tab_index = widget_info(tab_base, /tab_current)
        mmd_redraw_tab, topstate, curr_tab_index
    endelse
    
end

; Subroutine name: mmd_save_current_tab
; Created by: BT, 2008-06
; Calling Information:
;
;  event - event to be handled
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; handles the event of "Save Tab" button
;
; Editing Information:

pro mmd_save_current_tab, event

    compile_opt hidden
    common scan_data

    widget_control, event.top, get_uvalue=topstate

    active_tabid = mmd_get_tab_id(topstate, /active)
    active_drawid = widget_info(active_tabid, /child)

    widget_control, active_drawid, get_uvalue=drawstate

    if (not ptr_valid(drawstate)) then return

    ;; create a buffer to render the image 
    obuf = obj_new('idlgrbuffer', dimensions=(size(*((*drawstate).image_data)))[2:3]) 
    (*drawstate).oview->setProperty, $
      viewplane_rect=[0,0,(size(*((*drawstate).image_data)))[2:3]]

    obuf->draw, (*drawstate).oview
    obuf->getProperty, image_data=image

    (*drawstate).oview->setProperty, viewplane_rect=(*drawstate).vpr
    obj_destroy, obuf
    
    ;; resize to effective size (accounting for zoom level)
    image = reverse(congrid(image, $
                            3, $
                            (*drawstate).img_effective_size[0], $ 
                            (*drawstate).img_effective_size[1], $
                            cubic=-0.6, /center), 3)
    
    status = dialog_write_image(image, $
                                path=project.current_path, $
                                /warn_exist, $
                                /fix_type, $
                                type='TIFF', $
                                title='Save current tab as...')

end

; Subroutine name: mmd_mini_window_event
; Created by: BT, 2008-06
; Calling Information:
;
;  event - event to be handled
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; handles mouse click within mini-window (left tool bar) by
; centering the scroll rect on the mouse click's coordinates
;
; Editing Information:

pro mmd_mini_window_event, event

    widget_control, event.top, get_uvalue=topstate

    if (event.press) then begin
       (*topstate).draw_mini_click_hold = 1
    endif else if (event.release) then begin
       (*topstate).draw_mini_click_hold = 0
    endif

    if (event.press ne 4 and (*topstate).draw_mini_click_hold eq 0) then return

    ex = event.x
    ey = event.y

    geom = widget_info(event.id, /geometry)

    ;; not implemented
    if (event.press eq 4) then begin

         tab_base = (*topstate).tab_base
         all_tabs = widget_info(tab_base, /all_children)
         curr_tab = widget_info(tab_base, /tab_current)

         for i = 0, n_elements(all_tabs)-1 do begin
             
             draw_id = widget_info(all_tabs[i], /child)
             widget_control, draw_id, get_uvalue=drawstate

             if (i eq curr_tab) then final_drawstate = drawstate

             scale_factor = float((*drawstate).img_actual_size[0])/float(geom.draw_xsize)
             event.x = ex * scale_factor
             event.y = ey * scale_factor
             
             ;; we re-throw this event to the scroll handler
             mmd_do_scroll, drawstate, event, /center_on
             
         endfor
         
         mmd_update_mini_window, final_drawstate
         
     endif else begin
         
         tab_id = mmd_get_tab_id(topstate, /active)
         widget_control, widget_info(tab_id,/child), get_uvalue=drawstate
         
         scale_factor = float((*drawstate).img_actual_size[0])/float(geom.draw_xsize)
         event.x *= scale_factor
         event.y *= scale_factor
         
         mmd_do_scroll, drawstate, event, /center_on
         
     endelse

end

; Subroutine name: mmd_image_event
; Created by: BT, 2008-06
; Calling Information:
;
;  event  - event to be handled
;
; Bugs or Important Comments to Developers:
;
;  It is important to understand IDL's view objects "viewplane_rect"
;  property. The viewplane_rect 'moves' over the image and projects
;  the image data within its boundaries into the view window. changing
;  the size of the viewplane_rect is equivalent to zooming in and out
;  (since a smaller area (within the vpr) is projected into a larger
;  area (the window's view) while changing the location (x,y) is
;  equivalent to scrolling the image.
;
; Purpose of subroutine:
;
; Takes care of user clicks/drags based on currently selected tool.
;
; Editing Information:

pro mmd_image_event, event

    compile_opt hidden
    widget_control, event.id, get_uvalue=drawstate
    widget_control, event.top, get_uvalue=topstate

    current_tool = (*topstate).current_tool
    ;; if we are running in standalone mode, modifiers will determine
    ;; the tool
    if (event.modifiers eq 1) then current_tool = 'ZOOM'

    if (event.press eq 4) then begin
        ctm_id = widget_info(event.top, find_by_uname='drawContext')
        if (widget_info(ctm_id, /valid_id)) then begin
            widget_displaycontextmenu, event.id, event.x, event.y, ctm_id 
        endif
    endif else if (event.press and event.type eq 0) then begin
        ;; mouse click
        (*drawstate).is_down = 1
        (*drawstate).down_x = event.x
        (*drawstate).down_y = event.y
    endif else if (event.release and event.type eq 1) then begin
        ;; mouse release
        (*drawstate).is_down = 0
    endif else if (event.type eq 4) then begin
        ;; "expose" event - just dedraw
        (*drawstate).owindow->draw, (*drawstate).oview
        return
    endif else if (event.type eq 2) then begin
        ;; mouse-over event - update the float vals and position
        vals = long(mmd_get_pixel_value(drawstate, [event.x, event.y], $
                                        translated=tr_coord, float_val=fp))

        mmd_update_readout, topstate, xpos=tr_coord[0], ypos=tr_coord[1], fpval=fp

    endif

    case current_tool of
        
        'SCROLL': begin
            ;; we only care if mouse is being held
            if not (*drawstate).is_down then return
            mmd_do_scroll, drawstate, event
        end

        'ZOOM': begin
            mmd_do_zoom, drawstate, event
            mmd_update_readout, topstate, zoom=round(((*drawstate).zoom_amount*100.))
        end
        
        'DISTANCE': begin
            
            ;; we only care if they're click-holding (drag)
            if (not (*drawstate).is_down) then return
            
            ;; we want to draw a line on the actual image so we need
            ;; the translated coordinates.
            vals = long(mmd_get_pixel_value(drawstate,  $
                                            [(*drawstate).down_x, (*drawstate).down_y], $
                                            translated=tr_down))
            vals = long(mmd_get_pixel_value(drawstate, $
                                            [event.x, event.y], translated=tr_event))
            vec = [ [ tr_down[0] , tr_down[1] , 1 ], $
                    [ tr_event[0], tr_event[1], 1 ] ]

            ;; compute straight-line distance
            dx = 1.0 * (tr_down[0]-tr_event[0]-1)/(*drawstate).img_actual_size[0] * (*drawstate).fov[0]
            dy = 1.0 * (tr_down[1]-tr_event[1]-1)/(*drawstate).img_actual_size[1] * (*drawstate).fov[1]
            dist = sqrt(dx^2 + dy^2)

            ;; update the readout and show the distance line
            case (*drawstate).fov_units of
                0: un = ' px'
                1: un = ' cm'
                2: un = ' mm'
                else:
            endcase
            mmd_update_readout, topstate, dist=string(dist, format='(F0.3)')+un
            (*drawstate).odist->setProperty, data=vec, hide=0
            (*drawstate).owindow->draw, (*drawstate).oview
        end
        
        'INTENSITY': begin
            ;; we only care if they're click-holding (drag)
            if (event.clicks eq 2) then begin
                new_data = *((*drawstate).image_data)
                ;(*drawstate).oimage->setProperty, data=new_data
                ;(*drawstate).owindow->draw, (*drawstate).oview
                (*drawstate).image_it_max = max(new_data)
                (*drawstate).image_it_min = min(new_data)
                (*drawstate).image_it_top = 255 ; min(new_data)
                goto, INTENSITY_REDRAW
                return
            endif
            if (not (*drawstate).is_down) then return
            
            dx = ((*drawstate).down_x - event.x)
            dy = ((*drawstate).down_y - event.y)

            if (dy gt 1) then begin
                increment_y = 1
            endif else if (dy lt -1) then begin
                increment_y = -1
            endif else increment_y = 0

            if (dx gt 1) then begin
                increment_x = 1
            endif else if (dx lt -1) then begin
                increment_x = -1
            endif else increment_x = 0

            (*drawstate).image_it_max = ( ((*drawstate).image_it_max + 1.5*increment_y) > 1 ) < 255
            (*drawstate).image_it_min = ( ((*drawstate).image_it_min + 1.5*increment_x) > 1 ) < 255
            (*drawstate).image_it_top = ( ((*drawstate).image_it_top + 1.5*increment_x) > 1 ) < 255

            INTENSITY_REDRAW:
            
            new_data = bytscl( *((*drawstate).image_data),    $
                               min=(*drawstate).image_it_min, $
                               max=(*drawstate).image_it_max, $
                               top=(*drawstate).image_it_top)

            sz_img = size(new_data, /dimensions)
            ct = (*drawstate).image_color_table
            new_data[0,*,*] = reform(((*ct)[*,0])[(new_data[0,*,*])[*]], sz_img[1], sz_img[2]) 
            new_data[1,*,*] = reform(((*ct)[*,1])[(new_data[1,*,*])[*]], sz_img[1], sz_img[2]) 
            new_data[2,*,*] = reform(((*ct)[*,2])[(new_data[2,*,*])[*]], sz_img[1], sz_img[2]) 
            
            (*drawstate).oimage->setProperty, data=new_data
            (*drawstate).owindow->draw, (*drawstate).oview
            (*drawstate).down_y = event.y
            (*drawstate).down_x = event.x

            print, strcompress('it_max: '+string((*drawstate).image_it_max)+  $
                               ', it_min: '+string((*drawstate).image_it_min)+  $
                               ', it_top: '+string((*drawstate).image_it_top), /remove_all)
            
            mmd_draw_color_bar, drawstate, /update_only
        end
        
        'LINE_DATA': begin
              
           ;; this is where the data values come from
           if (not ptr_valid((*drawstate).float_vals)) then return

           ;; we only care if they're click-holding (drag)
           if (not (*drawstate).is_down) then return
           
           ;; we want to draw a line on the actual image so we need
           ;; the translated coordinates.
           vals = long(mmd_get_pixel_value(drawstate, $
                                           [event.x, event.y], translated=tr_event))

           ;; draw the lines
           line_pos = [ tr_event[0] + 0.5 , tr_event[1] + 0.5 ]
           (*drawstate).oline_data_h->setProperty, hide=0, $
              data=[ [line_pos[0],0],[line_pos[0],((*drawstate).img_actual_size)[1]] ]
           (*drawstate).oline_data_v->setProperty, hide=0, $
              data=[ [0,line_pos[1]],[((*drawstate).img_actual_size)[0],line_pos[1]] ]
           (*drawstate).owindow->draw, (*drawstate).oview
           
           looping = 0
           catch, error_state
           if (error_state ne 0) then begin
              WINDOW_NONEXISTENT:
              catch, /cancel ;; this line _may_ need to be before the label
              base = mas_open_line_data_window(group_leader=event.top, fov=(*drawstate).fov)
              (*drawstate).line_data_gui_win = base
              looping = 1
           endif

           if ((*drawstate).line_data_gui_win eq -1) then begin
              if (not looping) then goto, WINDOW_NONEXISTENT else return
           endif else begin
              mas_update_line_data_window, window_id=(*drawstate).line_data_gui_win, $
                                           horiz_data=(*((*drawstate).float_vals))[*,tr_event[1]], $
                                           vert_data=(*((*drawstate).float_vals))[tr_event[0],*]
           endelse

           catch, /cancel

        end

        else: print, "Tool not implemented: " + current_tool

    endcase

end

; Subroutine name: mmd_tool_event
; Created by: BT, 2008-06
; Calling Information:
;
;  event - event to be handled
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; handles the switching of tools in the tool bar
;
; Editing Information:

pro mmd_tool_event, event

    compile_opt hidden
    name = widget_info(event.id, /uname)
    widget_control, event.top, get_uvalue=topstate

    tab_id = mmd_get_tab_id(topstate, /active)
    widget_control, widget_info(tab_id,/child), get_uvalue=drawstate
    (*drawstate).odist->setProperty, data=[[0,0,0],[0,0,0]], hide=1
    (*drawstate).oline_data_h->setProperty, data=[[0,0],[0,0]], hide=1
    (*drawstate).oline_data_v->setProperty, data=[[0,0],[0,0]], hide=1
    (*drawstate).owindow->draw, (*drawstate).oview

    ;; sets the tool identifier and the cursor name for the tool
    case name of 
        'btn_zoom': begin
            (*topstate).current_tool='ZOOM'
            (*topstate).cursor = 'ZOOM'
        end

        'btn_scroll': begin
            (*topstate).current_tool='SCROLL'
            (*topstate).cursor = 'MOVE'
        end

        'btn_intensity': begin
            (*topstate).current_tool='INTENSITY'
            (*topstate).cursor = 'SIZE_NS'
        end

        'btn_distance' : begin
            (*topstate).current_tool='DISTANCE'
            (*topstate).cursor = 'CROSSHAIR'
        end
        
        'btn_line_data': begin
           (*topstate).current_tool='LINE_DATA'
           (*topstate).cursor = 'CROSSHAIR'

           ;; the window is created once the user clicks with this
           ;; tool. see CASE: 'LINE_DATA' in mmd_image_event
           
           if (event.select eq 0) then begin
              if (widget_info((*drawstate).line_data_gui_win, /valid_id)) then begin
                 widget_control, (*drawstate).line_data_gui_win, /destroy
                 (*drawstate).line_data_gui_win = -1
              endif
           endif

        end

        else: begin
            print, "Unknown Tool: "+name
            return
        end
    endcase

    tab_base = (*topstate).tab_base
    all_tabs = widget_info(tab_base, /all_children)
    
    for i=0, n_elements(all_tabs)-1 do begin
        draw_id = widget_info(all_tabs[i], /child)
        widget_control, draw_id, get_uvalue=drawstate
        if (ptr_valid(drawstate)) then begin
            (*drawstate).owindow->setCurrentCursor, (*topstate).cursor
        endif
    endfor

end

; Subroutine name: mmd_draw_context_menu_event
; Created by: BT, 2011-09
; Calling Information:
;
;  event - event to be handled
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; handles the loading of a color table based on result of context menu.
;
; Editing Information:

pro mmd_draw_context_menu_event, event
    
    widget_control, event.top, get_uvalue=topstate
    active_tabid = mmd_get_tab_id(topstate, /active)
    active_drawid = widget_info(active_tabid, /child)
    widget_control, active_drawid, get_uvalue=drawstate

    uname = widget_info(event.id, /uname)
    
    case uname of 
    
        'ctm_load_color_table': begin
            mmd_change_color_table, drawstate
        end        
        
        'ctm_show_color_bar': begin
            mmd_draw_color_bar, drawstate
        end
        
        else:
        
    endcase

end

;+
; :Description:
;    Opens or updates a window containing a color bar representation
;    of the image dynamic range of the current image window or tab
;
; :Params:
;    drawstate - the internal state pointer
;
; :Keywords:
;    update_only - set this keyword to have the color bar update only
;                  if there is a window already open
;
; :Author: wtriplett
;-
pro mmd_draw_color_bar, drawstate, update_only=update_only

    cbar_xsize = 150 ;; size of color bar window
    cbar_ysize = 250
    
    ;; get the window open state. 
    device, window_state=w
    if (w[12] ne 1) then begin
        if (keyword_set(update_only) eq 1) then return
        window, 12, xsize=cbar_xsize, ysize=cbar_ysize
    endif
    
    p_ct = (*drawstate).image_color_table

    if (not ptr_valid(p_ct)) then return
    
    ;; create a gradient and scale it to the current min/max intensity scale
    ind = long(findgen(cbar_ysize)/cbar_ysize * 255)
    ind = bytscl(ind, min=(*drawstate).image_it_min, $
                      max=(*drawstate).image_it_max, $
                      top=(*drawstate).image_it_top)
    
    ;; resize to window dimensions and colorize using color table
    gradient = reform(rebin(transpose(ind), cbar_xsize, cbar_ysize, 3))
    for c = 0, 2 do $
        gradient[*,*,c] = reform( ((*p_ct)[*,c])[(gradient[*,*,c])[*]], $
                                  cbar_xsize, cbar_ysize )

    ;; get the float values for the current window.
    p_fp = (*drawstate).float_vals

    ;; if window has float vals we can print those as a legend
    if (ptr_valid(p_fp)) then begin
        gradient[0:70,*,*] = 0
        if (ptr_valid( (*drawstate).float_range )) then begin
            p_fp_range = (*drawstate).float_range
            min_fp = (*p_fp_range)[0]
            min_fp = (*p_fp_range)[0]
            max_fp = (*p_fp_range)[1]
            mid_fp = (min_fp+max_fp)/2
        endif else begin
        min_fp = min(*p_fp, max=max_fp)
        mid_fp = (min_fp+max_fp)/2
        endelse
        
        tv, gradient, true=3

        n_scale_pts = 8
        scale_pts = (findgen(n_scale_pts)/(n_scale_pts-1)) * (max_fp-min_fp) + min_fp
        txt_ypos = (findgen(n_scale_pts)/(n_scale_pts-1)) * (0.925) + 0.025
        min_format='(%"%0.3e")'
        for p = 0, n_elements(scale_pts)-1 do begin
            print, 1.0 - (float(p)/(n_elements(scale_pts)-1)), scale_pts[p]
            xyouts, 0.05, txt_ypos[p], $
                    string(scale_pts[p], format=min_format), /normal
        endfor
        
    endif else begin
        tv, gradient, true=3
    endelse
    
end

;+
; :Description:
;    This procedure allows for the replacement of image data by
;    some external process, provided that the process has access
;    to the drawstate pointer. The drawstate pointer is returned
;    by the procedure that creates the display window.
;
; :Params:
;    drawstate - the internal state pointer
;    new_img_in - the new image data
;
; :Keywords:
;    new_fpvals - unscaled floating point values represented by the image
;
; :Author: wtriplett
;-
pro mmd_replace_image_data, drawstate, new_img_in, new_fpvals=new_fpvals, new_fprange=new_fprange

    imsz = size(new_img_in, /dimensions)
    
    mx = max(new_img_in, min=mn)
    
    if (n_elements(imsz) ne 3) then begin
        ;; black & white image -> three channel RGB
        image_data = bytarr(3,imsz[0], imsz[1])
        image_data[0,*,*] = new_img_in
        image_data[1,*,*] = new_img_in
        image_data[2,*,*] = new_img_in
        imsz = size(image_data, /dimensions)
    endif else if imsz[2] eq 3 then begin
        ;; three channel RGB -> swich channel order 
        image_data = bytarr(3,imsz[0], imsz[1])
        image_data[0,*,*] = new_img_in[*,*,0]
        image_data[1,*,*] = new_img_in[*,*,1]
        image_data[2,*,*] = new_img_in[*,*,2]
        imsz = size(image_data, /dimensions)
    endif else if imsz[0] eq 3 then begin
        ;; image is OK as-is
        image_data = new_img_in
    endif else begin
        void = dialog_message(["Image must be 2 channel greyscale or 3 channel RGB.", $
                               "Image Dimensions: ("+strcompress(string(imsz[0])+','+$
                                                                 string(imsz[1])+','+$
                                                                 string(imsz[2])+')',$
                                                                 /remove_all)],  $
                              /error, /center)
        return
    endelse
    

    
    (*drawstate).image_it_min = long(mn)
    (*drawstate).image_it_max = long(mx)
    (*drawstate).image_it_top = 255
    
    image_data = bytscl(image_data, min=(*drawstate).image_it_min, $
                                    max=(*drawstate).image_it_max, $
                                    top=(*drawstate).image_it_top)

    ct = (*drawstate).image_color_table
    image_data[0,*,*] = reform(((*ct)[*,0])[(image_data[0,*,*])[*]], imsz[1], imsz[2]) 
    image_data[1,*,*] = reform(((*ct)[*,1])[(image_data[1,*,*])[*]], imsz[1], imsz[2]) 
    image_data[2,*,*] = reform(((*ct)[*,2])[(image_data[2,*,*])[*]], imsz[1], imsz[2]) 
    
    (*drawstate).oimage->setProperty, data=image_data
    (*drawstate).owindow->draw, (*drawstate).oview
    
    if (n_elements(new_fpvals) ne 0) then begin 
        ptr_free, (*drawstate).float_vals
        (*drawstate).float_vals = ptr_new(new_fpvals)
    endif

    if (n_elements(new_fprange) ne 0) then begin
        ptr_free, (*drawstate).float_range
        (*drawstate).float_range = ptr_new(new_fprange)
    endif

    mmd_draw_color_bar, drawstate, /update_only
    
end

;+
; :Description:
;    Applies one of IDL's built-in color tables to the
;    displayed image.
;
; :Params:
;    drawstate - the state pointer
;
; :Author: wtriplett
;-
pro mmd_change_color_table, drawstate

    ;; IDL's common block for color tabling.
    common colors, R_CURR, R_ORIG, B_CURR, B_ORIG, G_CURR, G_ORIG
    
    ;; user selects color table
    xloadct, silent=0, /block
    new_img = *((*drawstate).image_data)
    sz_img = size(new_img, /dimensions)
    
    
    ;; now recolor the image with new color table
    (*drawstate).image_color_table = ptr_new(transpose([transpose(r_curr), $
                                                        transpose(g_curr), $
                                                        transpose(b_curr)]))

    new_img = bytscl(new_img, min=(*drawstate).image_it_min, $
                              max=(*drawstate).image_it_max, $
                              top=(*drawstate).image_it_top)
                              
    new_img[0,*,*] = reform(r_curr[(new_img[0,*,*])[*]], sz_img[1], sz_img[2]) 
    new_img[1,*,*] = reform(g_curr[(new_img[1,*,*])[*]], sz_img[1], sz_img[2]) 
    new_img[2,*,*] = reform(b_curr[(new_img[2,*,*])[*]], sz_img[1], sz_img[2]) 
    
    ;; redraw image
    (*drawstate).oimage->setProperty, data=new_img
    (*drawstate).owindow->draw, (*drawstate).oview
    
    loadct, 0, /silent

    mmd_draw_color_bar, drawstate, /update_only

end

; Subroutine name: mas_multi_display_event
; Created by: BT, 2008-06
; Calling Information:
;
; event - the event to be handled
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; handles resize event and tab switches. Takes care of updating the
; tabs' information with new sizes and also updates the 'mini window'
; if it exists.
;
; Editing Information:

pro mas_multi_display_event, event

    compile_opt hidden

    widget_control, event.top, get_uvalue=topstate
    
     if (tag_names(event, /structure_name) eq 'WIDGET_BASE') then begin
         ;; adjust for left toolbar size. IDL will add that much to
         ;; whatever window size the user ended up choosing.
         event.x -= (*topstate).toolbar_xsize
         tab_base = (*topstate).tab_base
         all_tabs = widget_info(tab_base, /all_children)
         curr_tab = widget_info(tab_base, /tab_current)

        for i=0, n_elements(all_tabs)-1 do begin
             draw_id = widget_info(all_tabs[i], /child)

             dx = (*topstate).view_xsize - event.x
             dy = (*topstate).view_ysize - event.y

             widget_control, draw_id, get_uvalue=drawstate
             if (not ptr_valid(drawstate)) then continue

             vpr = (*drawstate).vpr
             scroll_sense = (*drawstate).zoom_amount
             vpr[2] -= dx/scroll_sense
             vpr[3] -= dy/scroll_sense

             ;; bound the scrolling to the image dimensions
             ;; bottom-left corner
             vpr[0] = (vpr[0] + dx/scroll_sense) > 0.
             vpr[1] = (vpr[1] + dy/scroll_sense) > 0.
             
             ;; top-right corner
             if vpr[2] gt (*drawstate).img_actual_size[0] then begin
                 vpr[0] = ((*drawstate).img_actual_size[0] - vpr[2])/2.
             endif else if vpr[0] + vpr[2] gt (*drawstate).img_actual_size[0] then begin
                 vpr[0] = ((*drawstate).img_actual_size[0] - vpr[2])
             endif
             
             if vpr[3] gt (*drawstate).img_actual_size[1] then begin
                 vpr[1] = ((*drawstate).img_actual_size[1] - vpr[3])/2.
             endif else if vpr[1] + vpr[3] gt (*drawstate).img_actual_size[1] then begin
                 vpr[1] = ((*drawstate).img_actual_size[1] - vpr[3])
             endif
             
             (*drawstate).vpr = vpr
             (*drawstate).oview->setProperty, viewplane_rect=(*drawstate).vpr
             (*drawstate).draw_dims = [event.x, event.y]
             widget_control, draw_id, draw_xsize=event.x, draw_ysize=event.y

             mmd_update_mini_window, drawstate
             mmd_update_axis_label_position, drawstate

         endfor

         (*topstate).view_xsize = event.x
         (*topstate).view_ysize = event.y

     endif else begin
         mmd_redraw_tab, topstate, event.tab
     endelse

end

; Subroutine name: mmd_animate
; Created by: BT, 2008-06
; Calling Information:
;
; event - the event to be handled
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
;  Handles the "animate tabs" button click. Creates an animation from
;  all of the tabs currently in the window.
;
; Editing Information:

pro mmd_animate, event

    common scan_data
    ci = project.ci

    if (project.multi_display_tlb lt 0) then return

    widget_control, project.multi_display_tlb, get_uvalue=topstate
    if (not ptr_valid(topstate)) then begin
        return
    endif

    tab_base = (*topstate).tab_base
    all_tabs = widget_info(tab_base, /all_children)
    
    n_tabs = n_elements(all_tabs)
    if (n_tabs le 1) then return

    pimages = ptrarr(1)
    n_image = 0
    
    for i=0, n_tabs-1 do begin
        
        draw_id = widget_info(all_tabs[i], /child)
        widget_control, draw_id, get_uvalue=drawstate
        if (not ptr_valid(drawstate)) then continue
        
        if (n_image eq 0) then begin
            pimages[0] = (*drawstate).image_data
            im_sz = size(*((*drawstate).image_data), /dimensions)
            
            base = WIDGET_BASE(TITLE = project.imndarray[ci].display_name)
            animate = CW_ANIMATE(base, im_sz[1], im_sz[2], n_tabs, /NO_KILL )
        endif else begin
            pimages = [ pimages, (*drawstate).image_data ]
        endelse

        CW_ANIMATE_LOAD, animate, FRAME=n_image, IMAGE=*((*drawstate).image_data)

        n_image++
        
    endfor

    WIDGET_CONTROL, /REALIZE, base
    CW_ANIMATE_GETP, animate, pixmap_vect
    CW_ANIMATE_RUN, animate, 25
    XMANAGER, 'CW_ANIMATE Demo', base, EVENT_HANDLER = 'EHANDLER', /NO_BLOCK

end

; Subroutine name: mmd_tab_cleanup
; Created by: BT, 2008-06
; Calling Information:
;
; draw_id - the id of the draw widget to be cleaned up
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; destroys everything associated with a tab's draw widget before
; destroying the tab
;
; Editing Information:

pro mmd_tab_cleanup, draw_id

    compile_opt hidden
    widget_control, draw_id, get_uvalue=drawstate

    if (not ptr_valid(drawstate)) then return

    obj_destroy, [(*drawstate).oimage, (*drawstate).omodel, $
                  (*drawstate).oview, (*drawstate).owindow, $
                  (*drawstate).odist, (*drawstate).ozoom_box ]

    obj_destroy, [(*drawstate).mini_model, (*drawstate).mini_view, $
                  (*drawstate).mini_box ]

    ptr_free, (*drawstate).image_data
    ptr_free, (*drawstate).image_color_table
    ptr_free, (*drawstate).float_vals
    ptr_free, (*drawstate).extra_thumbnail
    ptr_free, drawstate

end

; Subroutine name: mmd_cleanup
; Created by: BT, 2008-06
; Calling Information:
;
; tlb -  the id of top level base to be cleaned up
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; destroys everything associated with a top level base
;
; Editing Information:

pro mmd_cleanup, tlb

    common scan_data
    compile_opt hidden

    widget_control, tlb, get_uvalue=state
    if (*state).tabbed eq 1 then begin
        project.multi_display_tlb = -1    
    endif

    obj_destroy, (*state).mini_window
    ptr_free, state

end

; Subroutine name: mas_display_multi
; Created by: BT, 2008-06
; Calling Information:
;
; image_data_in  -  IN - the image data to be displayed.
; tab_title      -  IN - the title of the tab or title bar
; standalone     -  IN - set this keyword to display the image using
;                        the standalone format. 
; fov_x          -  IN - the field of view in the x direction
; fov_y          -  IN - the field of view in the y direction
; fov_units      -  IN - the units for the field of view: 1 = cm, 2 =
;                        mm, 0 = pixels
; fp_vals        -  IN - pointer to array of values to be displayed
;                        in the "FP VAL" field
; fp_range       -  IN - the min/max range of the floating point values.These will be used
;                        to scale the displayed image. If not provided, min/max scaling will be used.
; drawstate_handle - OUT - contains a pointer to the internal state info
;                          you can really mess stuff up by changing things, 
;                          but it is required to have if you want to
;                          externally modify the image data without having
;                          to open a new display window.
;                          
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Displays an image
;
; Editing Information:

pro mas_display_multi, $
                       image_data_in,        $
                       tab_title=tab_title,  $
                       standalone=standalone,$
                       fov_x=fov_x,          $
                       fov_y=fov_y,          $
                       fov_units=fov_units,  $
                       fp_vals=fp_vals,      $
                       fp_range=fp_range,    $
                       show_axis_labels=show_axis_labels, $
                       axis_labels=axis_labels, $
                       extra_thumbnail=extra_thumbnail, $
                       drawstate_handle=drawstate

    common scan_data
    ci = project.ci

    show_axis_labels = project.procpramarray[ci].have_orientation

    if n_elements(axis_labels) ne 4 then begin
       axis_labels = mas_orient_get_axis_labels(/apply_flip_rotate)
    endif

    imsz = size(image_data_in, /dimensions)
    if (n_elements(imsz) ne 3) then begin
        ;; black & white image -> three channel RGB
        image_data = bytarr(3,imsz[0], imsz[1])
        image_data[0,*,*] = image_data_in
        image_data[1,*,*] = image_data_in
        image_data[2,*,*] = image_data_in
        imsz = size(image_data, /dimensions)
    endif else if imsz[2] eq 3 then begin
        ;; three channel RGB -> swich channel order 
        image_data = bytarr(3,imsz[0], imsz[1])
        image_data[0,*,*] = image_data_in[*,*,0]
        image_data[1,*,*] = image_data_in[*,*,1]
        image_data[2,*,*] = image_data_in[*,*,2]
        imsz = size(image_data, /dimensions)
    endif else if imsz[0] eq 3 then begin
        ;; image is OK as-is
        image_data = image_data_in
    endif else begin
        void = dialog_message(["Image must be 2 channel greyscale or 3 channel RGB.", $
                               "Image Dimensions: ("+strcompress(string(imsz[0])+','+$
                                                                 string(imsz[1])+','+$
                                                                 string(imsz[2])+')',$
                                                                 /remove_all)],  $
                              /error, /center)
        return
    endelse

    if (not keyword_set(fov_units)) then fov_units = 0 else begin
        ;; make sure units are ok
    endelse
    
    if (not keyword_set(fov_x)) then fov_x = imsz[1] else begin
        ;; make sure units are ok
    endelse

    if (not keyword_set(fov_y)) then fov_y = imsz[2] else begin
        ;; make sure units are ok
    endelse

    if (n_elements(fp_vals) ne 0) then begin
        fp_vals_in = ptr_new(fp_vals)
    endif

    fov_w = [fov_x, fov_y]

    if keyword_set(standalone) then begin

        if (imsz[2] gt 600) then smw=1

        tlb = mmd_init(/simple, show_mini_window=smw, $
                       show_line_data_btn=ptr_valid(fp_vals_in))
                       
        widget_control, tlb, get_uvalue=topstate

        if (keyword_set(tab_title)) then begin
            widget_control, tlb, tlb_set_title=tab_title
        endif

        mmd_add_tab, $
           topstate,  $
           image_data,$
           fov_in=fov_w,$
           fov_units=fov_units,$
           fp_vals=fp_vals_in, $
           fp_range=fp_range, $
           show_axis_labels=show_axis_labels, $
           axis_labels=axis_labels, $
           drawstate=drawstate

    endif else begin

        if (project.multi_display_tlb eq -1) then begin
            tlb = mmd_init()
            project.multi_display_tlb = tlb
        endif else begin
            tlb = project.multi_display_tlb 
        endelse

        widget_control, tlb, get_uvalue=topstate

        mmd_add_tab, $
           topstate,  $
           image_data,$
           fov_in=fov_w,$
           fov_units=fov_units,$
           title=tab_title, $
           fp_vals=fp_vals_in, $
           fp_range=fp_range, $
           show_axis_labels=show_axis_labels, $
           axis_labels=axis_labels, $
           drawstate=drawstate

    endelse

    if (keyword_set(extra_thumbnail)) then begin
        
        wd = widget_info((*topstate).extra_base, /child)
        
        if (not widget_info(wd, /valid_id)) then begin
            wd = widget_draw((*topstate).extra_base, /align_center, $
                             graphics_level=1, retain=2)
            widget_control, (*topstate).tlb, /realize
        endif
        
        geom = widget_info(wd, /geometry)
        ex = ptr_new(congrid(extra_thumbnail, geom.draw_xsize, $
                             geom.draw_ysize, 3), /no_copy)
        
        widget_control, wd, get_value=wind
        wset, wind
        
        tv, *ex, true=3
        
        (*topstate).extra_base_filled = 1
        (*drawstate).extra_thumbnail = ex
        
    endif

end
