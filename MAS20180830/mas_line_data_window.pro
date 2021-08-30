; Subroutine name: mas_line_data_option_window_event
; Created by: BT, 2009-03
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:

pro mas_line_data_option_window_event, event

    name = widget_info(event.id, /uname)
    
    widget_control, event.top, get_uvalue=state

    case name of 
       
       'btn_ok': begin
          plot_state = (*state).plot_state

          widget_control, (*state).txt_save_resolution, get_value=val
          if (strcompress(val, /remove_all) ne '') then begin
             (*plot_state).save_resolution = float(val)
          endif
          widget_control, (*state).txt_save_width, get_value=val
          if (strcompress(val, /remove_all) ne '') then begin
             (*plot_state).save_width = float(val)
          endif
          widget_control, (*state).txt_save_height, get_value=val
          if (strcompress(val, /remove_all) ne '') then begin
             (*plot_state).save_height = float(val)
          endif
          
          widget_control, (*state).sldr_yscale, get_value=val
          (*plot_state).yscale = float(val)
          widget_control, (*state).sldr_plotline_thk, get_value=val
          (*plot_state).plotline_thick = fix(val)
          
          widget_control, (*state).txt_vert_title, get_value=val
          (*plot_state).vert_title = val
          widget_control, (*state).txt_vert_xtitle, get_value=val
          (*plot_state).vert_xtitle = val
          widget_control, (*state).txt_vert_ytitle, get_value=val
          (*plot_state).vert_ytitle = val
          
          widget_control, (*state).txt_horiz_title, get_value=val
          (*plot_state).horiz_title = val
          widget_control, (*state).txt_horiz_xtitle, get_value=val
          (*plot_state).horiz_xtitle = val
          widget_control, (*state).txt_horiz_ytitle, get_value=val
          (*plot_state).horiz_ytitle = val

          (*plot_state).save_in_color = widget_info((*state).btn_save_color, /button_set)

         mas_update_line_data_window, window_id=(*plot_state).window_id, $
                                       /refresh_only
 
          ptr_free, state
          widget_control, event.top, /destroy
       end

       'btn_save_color': begin ;; no op
       end

       'sldr_yscale': begin
          widget_control, (*state).sldr_yscale, get_value=val
          plot_state = (*state).plot_state
          (*plot_state).yscale = val
          
          mas_update_line_data_window, window_id=(*plot_state).window_id, $
                                       /refresh_only
       end
       'sldr_plotline_thk': begin
          widget_control, (*state).sldr_plotline_thk, get_value=val
          plot_state = (*state).plot_state
          (*plot_state).plotline_thick = val
          mas_update_line_data_window, window_id=(*plot_state).window_id, $
                                       /refresh_only
       end
       
       'btn_xaxis_pixels': begin
          plot_state = (*state).plot_state
          (*plot_state).xaxis_units = 'PIXELS'
          widget_control, (*state).txt_vert_xtitle, set_value='Pixel Location'
          widget_control, (*state).txt_horiz_xtitle, set_value='Pixel Location'
          (*plot_state).vert_xtitle = 'Pixel Location'
          (*plot_state).horiz_xtitle = 'Pixel Location'
          
          mas_update_line_data_window, window_id=(*plot_state).window_id, $
                                       /refresh_only
       end
       
       'btn_xaxis_cm': begin
          plot_state = (*state).plot_state
          if ((*plot_state).have_fov eq 0) then return
          (*plot_state).xaxis_units = 'CM'
          widget_control, (*state).txt_vert_xtitle, set_value='CM'
          widget_control, (*state).txt_horiz_xtitle, set_value='CM'
          (*plot_state).vert_xtitle = 'CM'
          (*plot_state).horiz_xtitle = 'CM'
          mas_update_line_data_window, window_id=(*plot_state).window_id, $
                                       /refresh_only
       end
       
       'btn_cancel': begin
          plot_state = (*state).plot_state
          (*plot_state).yscale = (*state).orig_yscale
          (*plot_state).plotline_thick = (*state).orig_plotline_thick
          mas_update_line_data_window, window_id=(*plot_state).window_id, $
                                       /refresh_only
          widget_control, event.top, /destroy
       end

       else: print, "recv'd event from unknown widget: "+name

    endcase

end

; Subroutine name: mas_line_data_window_event
; Created by: BT, 2009-03
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:

pro mas_line_data_window_event, event

    name = widget_info(event.id, /uname)
    
    case name of

       'save': mas_line_data_save_tiff, event
       'opts': mas_line_data_option_gui, event

    endcase

end

; Subroutine name: mas_line_data_option_gui
; Created by: BT, 2009-03
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:

pro mas_line_data_option_gui, event

    widget_control, event.top, get_uvalue=state
    
    base = widget_base(group_leader=event.top, /modal, $
                       title='Plotting Options', $
                       /column, /align_left)

    ;; export to tiff options
    void = widget_base(base, /column, /frame, /align_left)
    lbl = widget_label(void, value="Export to TIFF options:", /align_left)

    B0 = widget_base(void, /row)

    B1 = widget_base(B0, /row)
    lbl = widget_label(b1, value="Width:")
    txt_save_width = widget_text(b1, xsize=4, /editable, $
                                 uname="txt_save_width", $
                                 value=strcompress((*state).save_width, /remove_all))
    lbl = widget_label(b1, value="in.")

    B1 = widget_base(B0, /row)
    lbl = widget_label(b1, value="Height:")
    txt_save_height = widget_text(b1, xsize=4, /editable, $
                                  uname="txt_save_height", $
                                  value=strcompress((*state).save_height, /remove_all))
    lbl = widget_label(b1, value="in.")

    B1 = widget_base(B0, /row)
    lbl = widget_label(b1, value="Resolution:")
    txt_save_resolution = widget_text(b1, xsize=4, /editable, $
                                      uname="txt_save_resolution", $
                                      value=strcompress((*state).save_resolution, /remove_all))
    lbl = widget_label(b1, value="dpi")

    B2 = widget_base(void, /nonexclusive)
    btn_save_color = widget_button(b2, value="Save in color", $
                                   uname='btn_save_color')
    widget_control, B2, set_button=(*state).save_in_color


    ;; plot titles options
    BBB = widget_base(base, /row, /grid)
    
    B1 = widget_base(BBB, /column, /frame)
    lbl = widget_label(B1, value="Vertical Plot Labels:", /align_left)

    BB = widget_base(b1, /row)
    lbl = widget_label(bb, value="Title :")
    txt_vert_title = widget_text(bb, xsize=20, /editable, uname="txt_vert_title", $
                                 value=(*state).vert_title)

    BB = widget_base(b1, /row)
    lbl = widget_label(bb, value="X Axis:")
    txt_vert_xtitle = widget_text(bb, xsize=20, /editable, uname="txt_vert_xtitle", $
                                  value=(*state).vert_xtitle)

    BB = widget_base(b1, /row)
    lbl = widget_label(bb, value="Y Axis:")
    txt_vert_ytitle = widget_text(bb, xsize=20, /editable, uname="txt_vert_ytitle", $
                                  value=(*state).vert_ytitle)

    B2 = widget_base(BBB, /column, /frame)
    lbl = widget_label(B2, value="Horizontal Plot Labels:", /align_left)

    BB = widget_base(B2, /row)
    lbl = widget_label(bb, value="Title :")
    txt_horiz_title = widget_text(bb, xsize=20, /editable, uname="txt_horiz_title", $
                                  value=(*state).horiz_title)
    BB = widget_base(B2, /row)
    lbl = widget_label(bb, value="X Axis:")
    txt_horiz_xtitle = widget_text(bb, xsize=20, /editable, uname="txt_horiz_xtitle", $
                                  value=(*state).horiz_xtitle)

    BB = widget_base(B2, /row)
    lbl = widget_label(bb, value="Y Axis:")
    txt_horiz_ytitle = widget_text(bb, xsize=20, /editable, uname="txt_horiz_ytitle", $
                                  value=(*state).horiz_ytitle)
    


    BBB = widget_base(base, /row, /grid, /align_center)

    void = widget_base(BBB, /column, /frame, /align_left)
    lbl = widget_label(void, value="X-axis Units:")
    BOM = widget_base(void, /column, /exclusive)
    btn_pixels = widget_button(BOM, value='Pixels', uname='btn_xaxis_pixels')
    btn_cm     = widget_button(BOM, value='CM', uname='btn_xaxis_cm', sensitive=(*state).have_fov)
    widget_control, btn_pixels, set_button=((*state).xaxis_units eq 'PIXEL')
    widget_control, btn_cm    , set_button=((*state).xaxis_units eq 'CM')
    
    ;; axis scaling option
    void = widget_base(BBB, /column, /frame, /align_left)
    lbl = widget_label(void, value="Axis Scaling:")
    sldr_yscale  = widget_slider(void, title="Y-Axis Scaling", uname="sldr_yscale", $
                                 min=1, max=200, /drag, value=(*state).yscale)
    
    ;; plotline thickness
    void = widget_base(BBB, /column, /frame, /align_left)
    lbl = widget_label(void, value="Plotline Thickness")
    sldr_plotline_thk = widget_slider(void, title="Plotline Thickness", uname="sldr_plotline_thk", $
                                      min=1, max=5, /drag, value=(*state).plotline_thick)
                                      
    ;; ok/cancel buttons
    btnbase = widget_base(base, /row, /align_center)
    btn_cancel = widget_button(btnbase, value="Cancel", uname="btn_cancel")
    btn_ok = widget_button(btnbase, value="Save Settings", uname="btn_ok")

    widget_control, base, set_uvalue=ptr_new({$
                    btn_ok: btn_ok, $
                    btn_cancel: btn_cancel, $
                    btn_save_color: btn_save_color, $
                    txt_save_resolution: txt_save_resolution, $
                    txt_save_height: txt_save_height, $
                    txt_save_width: txt_save_width, $
                    txt_vert_title: txt_vert_title, $
                    txt_vert_xtitle: txt_vert_xtitle, $
                    txt_vert_ytitle: txt_vert_ytitle, $
                    txt_horiz_title: txt_horiz_title, $
                    txt_horiz_xtitle: txt_horiz_xtitle, $
                    txt_horiz_ytitle: txt_horiz_ytitle, $
                    sldr_yscale: sldr_yscale, $
                    sldr_plotline_thk: sldr_plotline_thk, $
                    orig_yscale: (*state).yscale, $
                    orig_plotline_thick: (*state).plotline_thick, $
                    plot_state: state})
    
    widget_control, base, /realize
    xmanager, 'mas_line_data_option_window', base

end

; Subroutine name: mas_line_data_save_tiff
; Created by: BT, 2009-03
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:

pro mas_line_data_save_tiff, event

    widget_control, event.top, get_uvalue=state

    export_w = (*state).save_width  * (*state).save_resolution
    export_h = (*state).save_height * (*state).save_resolution
    export_dims = [export_w, 2*export_h]
    
    old_plot = !D.name
    old_p_multi = !P.multi
    set_plot, 'Z'
    !P.multi=[0,1,2] ;; set device to 2-up vertical
    if ( (*state).save_in_color ) then begin
       device, set_pixel_depth=24
       data = bytarr(3, export_w, 2*export_h)
       colorv = '0000FF'x
       colorh = 'FF0000'x
    endif else begin
       device, set_pixel_depth=8
       data = bytarr(export_w, 2*export_h)
       colorv = '000000'x
       colorh = colorv
    endelse

    device, set_resolution=export_dims, decomposed=1

    plot, *((*state).vert_data), /device, $
          yrange=[(*state).yrange_min,(*state).yrange_max], $
          background='FFFFFF'x, color=colorv, xstyle=1, $
          title=(*state).vert_title, thick=(*state).plotline_thick, $
          xtitle=(*state).vert_xtitle, ytitle=(*state).vert_ytitle

    plot, *((*state).horiz_data), /device, $
          yrange=[(*state).yrange_min,(*state).yrange_max], $
          background='FFFFFF'x, color=colorh, xstyle=1, $
          title=(*state).horiz_title, thick=(*state).plotline_thick, $
          xtitle=(*state).horiz_xtitle, ytitle=(*state).horiz_ytitle

    if ( (*state).save_in_color  ) then begin
       data = reverse(tvrd(true=1), 3)
    endif else begin
       data = reverse(tvrd(),2)
    endelse
    
    !P.multi = old_p_multi
    set_plot, old_plot
    
    filename = dialog_pickfile(dialog_parent=event.top, $
                               filter=['*.tif','*.tiff'], $
                               file='plot_output.tif', $
                               /fix_filter, $
                               /write)

    if (filename ne '') then begin
       write_tiff, filename, data, units=2, $
                   xresol=(*state).save_resolution, $
                   yresol=(*state).save_resolution
    endif

end

; Subroutine name: mas_open_line_data_window
; Created by: BT, 2009-03
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:

function mas_open_line_data_window, group_leader=group_leader, fov=fov

    GL = -1
    if (keyword_set(group_leader)) then begin
       if (widget_info(group_leader, /valid_id)) then begin
          GL = group_leader
       endif
    endif
    
    if (n_elements(fov) ne 2) then begin
        fov_xy = [0.0,0.0]
        have_fov = 0B
        fov_units = 'PIXEL'
    endif else begin
        fov_xy = float(fov)
        have_fov = 1B
        fov_units = 'CM'
    endelse
    
    if (GL ne -1) then begin
       base = widget_base(title="Line Profile", group_leader=GL, mbar=mbar)
    endif else begin
       base = widget_base(title="Line Profile", mbar=mbar)
    endelse

    file = widget_button(mbar, value="Options")
    save = widget_button(file, uname='save', value="Save as TIFF...")
    opts = widget_button(file, uname='opts', value="Plotting Options...")

    draw = widget_draw(base, xsize=300, ysize=400, retain=2)
    widget_control, base, /realize

    widget_control, draw, get_value=window
    
    case strlowcase(!VERSION.os_family) of 
       'unix': plot_device = 'x'
       'windows': plot_device = 'win'
       else: plot_device = 'x' ;; This assumes that Mac = unix, which
                               ;; it is in mac os x.
    endcase

    if (fov_units eq 'CM') then begin
        horiz_xtitle = 'CM'
        vert_xtitle  = 'CM'
    endif else begin
        horiz_xtitle = 'Pixel Location'
        vert_xtitle  = 'Pixel Location'
    endelse
    
    state = ptr_new({ LINE_DATA_STATE, $ 
                      window_id: base, $
                      plot_device: plot_device, $
                      draw: draw, $
                      window: window, $
                      yscale: 100, $
                      yrange_max: -1., $
                      yrange_min: 0., $
                      plotline_thick: 1, $
                      save_width: 3, $ ;; inches
                      save_height: 2, $;; inches
                      save_resolution: 300, $ ;; dpi
                      save_in_color: 0, $
                      xaxis_units: (have_fov eq 1) ? 'CM' : 'PIXEL', $
                      fov_xy: fov_xy, $
                      have_fov: have_fov, $
                      horiz_data: ptr_new(), $
                      vert_data: ptr_new(), $
                      horiz_title: "Horizontal Line Plot", $
                      horiz_xtitle: horiz_xtitle, $
                      horiz_ytitle: "Pixel Value", $
                      vert_title: "Vertical Line Plot", $
                      vert_xtitle: vert_xtitle, $
                      vert_ytitle: "Pixel Value", $
                      user_defined: ptr_new() })

    widget_control, base, set_uvalue=state
    xmanager, 'mas_line_data_window', base, /no_block
    return, base

end

; Subroutine name: mas_update_line_data_window
; Created by: BT, 2009-03
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; Editing Information:

pro mas_update_line_data_window, window_id=window_id, $
                                 refresh_only=refresh_only, $
                                 horiz_data=horiz_data, $ 
                                 vert_data=vert_data, $
                                 yrange_max=yrange_max, $
                                 yrange_min=yrange_min, $
                                 background=background, fov=fov
    
    widget_control, window_id, get_uvalue=state

    if (keyword_set(refresh_only)) then begin
       if ( not ptr_valid((*state).vert_data) or $
            not ptr_valid((*state).horiz_data) ) then begin
          return
       endif
    endif else begin
       ptr_free, [ (*state).vert_data, (*state).horiz_data ]
       (*state).vert_data  = ptr_new(vert_data)
       (*state).horiz_data = ptr_new(horiz_data)
    endelse

    if (keyword_set(yrange_max)) then begin
       (*state).yrange_max = yrange_max
    endif else begin
       (*state).yrange_max = max([ max(*((*state).vert_data)), $
                                   max(*((*state).horiz_data)) ])
    endelse

    if (keyword_set(yrange_min)) then begin
       (*state).yrange_min = yrange_min
    endif else begin
       (*state).yrange_min = min([ min(*((*state).vert_data)), $
                                   min(*((*state).horiz_data)) ])
    endelse
    
    if (n_elements(fov) eq 2) then begin
        (*state).fov_xy = float(fov)
        (*state).have_fov = 1B
    endif
    
    nv = n_elements(*((*state).vert_data))
    nh = n_elements(*((*state).horiz_data))
    if ((*state).xaxis_units eq 'CM') then begin
        xaxis_vert  = findgen(nv)/(nv-1) * (*state).fov_xy[1] - ((*state).fov_xy[1]/2) 
        xaxis_horiz = findgen(nh)/(nh-1) * (*state).fov_xy[0] - ((*state).fov_xy[0]/2)
    endif else begin
        xaxis_vert  = findgen(nv)
        xaxis_horiz = findgen(nh)
    endelse
    
    yrange_max = (*state).yrange_max * float((*state).yscale)/100.
    yrange_min = (*state).yrange_min * float((*state).yscale)/100.

    ;; Get and save original plot and device values before changing
    ;; them.

    old_plot = !D.name
    set_plot, (*state).plot_device

    old_p_multi = !P.multi
    !P.multi = [0,1,2]

    device, get_decomposed=old_decomp
    device, decomposed=1

    ;; set active the our line data window
    wset, (*state).window
   
    plot, xaxis_vert,  *((*state).vert_data),  yrange=[yrange_min, yrange_max], $
          background='FFFFFF'x, color='0000FF'x, xstyle=1, $
          title=(*state).vert_title, thick=(*state).plotline_thick, $
          xtitle=(*state).vert_xtitle, ytitle=(*state).vert_ytitle

    plot, xaxis_horiz, *((*state).horiz_data), yrange=[yrange_min, yrange_max], $
          background='FFFFFF'x, color='FF0000'x, xstyle=1, $
          title=(*state).horiz_title, thick=(*state).plotline_thick, $
          xtitle=(*state).horiz_xtitle, ytitle=(*state).horiz_ytitle

    ;; restore plot and device values before returning
    device, decomposed=old_decomp
    !P.multi = old_p_multi
    set_plot, old_plot

end

pro mas_histogram_barplot, coord, scan_list_id

    common scan_data

    values = reform((*project.dataarray[scan_list_id].state1)[coord[0], $
                                                              coord[1], $
                                                              coord[2], *])
    hist = histogram(values, nbins=15, locations=loc)

    bar_plot, hist, $
              barnames=string(loc, format='(F0.3)'), $
              barwidth=0.75, $
              colors=indgen(n_elements(hist))+128
end
