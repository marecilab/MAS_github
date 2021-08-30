; IDL routine to display 2D maximum intensity projection (MIP) of a volume

; Subroutine name: mip_gui
; Created by: Magdoom Kulam
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Main MIP GUI
;
; Editing Information:

pro mip_gui

  common scan_data
  common common_widgets
  common mip_params, zhi, zlo, thres
  
  if (xregistered('mip_gui') ne 0) then return

  CI = project.ci
  slice_axis = project.procPramArray[CI].slice_axis     ; Slice axis (Read/Phase,Read/Slice or Phase/Slice)
  is_3d = project.imndarray[CI].dimensions eq 3
  zmin = 0
  case slice_axis of
   0 : zmax = project.imndArray[CI].sdim*project.procpramArray[CI].slice_interp*(1+is_3d*project.procpramArray[CI].zpad_flag)-1
   1 : zmax = project.imndArray[CI].pdim*project.procpramArray[CI].phase_interp*(1+is_3d*project.procpramArray[CI].zpad_flag)-1
   2 : zmax = project.imndArray[CI].fdim*project.procpramArray[CI].freq_interp*(1+is_3d*project.procpramArray[CI].zpad_flag)-1  
  endcase
   
  topbase = widget_base(title = 'MIP : Scan #' + strtrim(project.ci+1,2) ,tlb_frame_attr=1,/grid_layout,/kbrd_focus_events)
  mainbase = widget_base(topbase,row=3)

  zbase = Widget_Base(mainbase, row = 2, /ALIGN_CENTER)
  zmin_MIP = widget_slider(zbase,title = 'Min Z', uvalue = 'minZ', value = 0, maximum = zmax, minimum = 0, xsize = 150)
  zmax_MIP = widget_slider(zbase,title = 'Max Z', uvalue = 'maxZ', value = zmax, maximum = zmax, minimum = 0, xsize = 150)
  
  Ithres_MIP = widget_slider(mainbase, title = 'Image Threshold (%)', uvalue = 'Image Threshold',xsize = 160)

  disp_button_base = widget_base(mainbase,column = 2)
  disp_button_RE = widget_button(disp_button_base, value = 'Display', uvalue = 'Display',/no_release, xsize = 80)
  movie_button_RE = widget_button(disp_button_base, value = 'Movie', uvalue = 'Movie',/no_release, xsize = 80)
  
  thres = 0
  zlo = zmin
  zhi = zmax
  
  widget_control, topbase, /realize
  xmanager, 'mip_gui', topbase, cleanup='mip_gui_cleanup',/NO_BLOCK

end

; Subroutine name: region_extract_gui_event
; Created by: Magdoom Kulam
; Calling Information:
;
;  event
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Event handler for the main eddy current gui
;
; Editing Information:

pro mip_gui_event, event

  common scan_data
  common common_widgets
  common mip_params

  ; To catch error ahead and prevent crashing
  catch,error_status
  IF (error_status NE 0) THEN BEGIN
     help, calls= trace_back
     dummy = DIALOG_MESSAGE(['Please send a copy of this error message to the developer of this program.', $
     !ERROR_STATE.msg,trace_back], /ERROR, TITLE='Error in WID_BASE_MAIN_event')
     widget_control, event.top, /destroy
     RETURN
  ENDIF

  widget_control, event.id, get_uvalue = widget
  CI = project.ci
  
  if event.top eq event.id then begin 
    widget_control, event.top, tlb_set_title = 'MIP : Scan #' + strtrim(CI+1,2)
    slice_axis = project.procPramArray[CI].slice_axis     ; Slice axis (Read/Phase,Read/Slice or Phase/Slice)
    is_3d = project.imndarray[CI].dimensions eq 3
    zmin = 0
    case slice_axis of
      0 : zmax = project.imndArray[CI].sdim*project.procpramArray[CI].slice_interp*(1+is_3d*project.procpramArray[CI].zpad_flag)-1
      1 : zmax = project.imndArray[CI].pdim*project.procpramArray[CI].phase_interp*(1+is_3d*project.procpramArray[CI].zpad_flag)-1
      2 : zmax = project.imndArray[CI].fdim*project.procpramArray[CI].freq_interp*(1+is_3d*project.procpramArray[CI].zpad_flag)-1  
    endcase
    widget_control, zmin_MIP, set_slider_min = zmin
    widget_control, zmin_MIP, set_slider_max = zmax
    widget_control, zmax_MIP, set_slider_min = zmin
    widget_control, zmax_MIP, set_slider_max = zmax
    
    if zhi gt zmax then zhi = zmax
    if zlo gt zmax then zlo = zmax-1
         
  endif else begin
    case widget of

      'minZ'            : zlo = event.value

      'maxZ'            : zhi = event.value
      
      'Image Threshold' : thres = event.value

      'Display'         : begin
                          if project.procPramArray[CI].single_Multi_flag eq 0 then single_Multi_flag_toggle
                          mas_load_state_2
                          data = *project.dataarray[CI].state2
                          N = size(data,/Dimension)
                          zlim = [zlo,zhi]
                          if n_elements(N) ne 3 then h = dialog_message('Please select a multi-slice or 3D dataset',/error) $
                          else if zlo lt zhi then mas_mip,data, zlim = zlim, threshold = thres, 0
                          end   
      
      'Movie'           : begin
                          if project.procPramArray[CI].single_Multi_flag eq 0 then single_Multi_flag_toggle
                          mas_load_state_2
                          data = *project.dataarray[CI].state2
                          N = size(data,/Dimension)
                          zlim = [zlo,zhi]
                          if n_elements(N) ne 3 then h = dialog_message('Please select a multi-slice or 3D dataset',/error) $
                          else if zlo lt zhi then mas_mip,data, zlim = zlim, threshold = thres, 1
                          end   

      else              : return

    endcase  
  endelse
  
end

; Subroutine name: mas_mip
; Created by: Magdoom Kulam
; Calling Information:
; 
; data      - 3D array
; z1        - Lower clipping limit index in the 3rd dimension
; z2        - Upper clipping limit index in the 3rd dimension
; threshold - Threshold level (%) of maximum of clipped 3D array
; type      - 2D(0) or cine(1)
; 
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Displays 2D or cine MIP of the given 3D array limited in 3rd index by [z1,z2] with a threshold.
;
; Editing Information:

pro mas_mip,data, zlim = zlim, threshold = threshold, type

  common scan_data
  
  CI = project.ci
  if keyword_set(zlim) then Isub = data[*,*,zlim[0]:zlim[1]] else Isub = data
  if keyword_set(threshold) then begin
    Imask= abs(Isub) ge max(abs(Isub))*float(threshold)/100.0
    Isub*= Imask
  endif
  N = size(Isub,/Dimension)
    
  if type eq 0 then begin
     slice_axis = project.procPramArray[CI].slice_axis     ; Slice axis (Read/Phase,Read/Slice or Phase/Slice)
     case slice_axis of
        0 : islice = project.procPramArray[CI].sdim_start
        1 : islice = project.procPramArray[CI].pdim_start
        2 : islice = project.procPramArray[CI].fdim_start
     endcase
     iimage, data[*,*,islice], layout = [2,1,1], aspect_ratio = N[1]/N[0], title = 'Reference Image', font_size = 15
     iimage, max(Isub,dimension = 3), layout = [2,1,2],/current, aspect_ratio = N[1]/N[0], title = 'MIP', font_size = 15
   
  endif else begin
    
     top_base = widget_base(title = 'Cine MIP',/grid_layout)
     anim = cw_animate(top_base,N[0],N[1],19, /NO_KILL)
      
     progressbar = Obj_New('progressbar', Color='red', Text='Loading frames', /NOCANCEL)
     progressbar -> Start
     ; Load the animation
     FOR j=0,18 DO BEGIN
         rotData = Transform_Volume(Isub,Rotation=[0,(j*10) MOD 360,0], Missing=0) 
         cw_animate_load, anim, frame = j, Image = Reform(Max(bytscl(rotData), DIMENSION=3))
         progressBar -> Update, (float(j+1)/float(19))*100.0
     ENDFOR
     progressbar -> Destroy

     widget_control, /realize, top_base
     cw_animate_run, anim,1
     XMANAGER, 'CW_ANIMATE Demo', top_base, EVENT_HANDLER = 'anim_ehandler'
    endelse
    
end

; Subroutine name: region_extract_gui_cleanup
; Created by: Magdoom Kulam
; Calling Information:
;
; top_base
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Cleanup procedure before exiting the eddy current GUI
;
; Editing Information:

pro mip_gui_cleanup, topbase

  widget_control, topbase, get_uvalue=state
  if (ptr_valid(state)) then ptr_free, state
  return

end
