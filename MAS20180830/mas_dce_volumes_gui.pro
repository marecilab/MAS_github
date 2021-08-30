; Subroutine name: dce_volumes_ROI_edit
; Created by: Garrett Astary
; Calling Information:
;
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Selects ROI for each slice
;
; Editing Information:

pro dce_volumes_ROI_edit, info, event
  COMPILE_OPT IDL2
  COMMON scan_data
  COMMON common_widgets
  CI = project.CI
  
  selected_slice = info.selected_slice
  main_window_slice = project.procPramArray[CI].sdim_start 
  project.procPramArray[CI].sdim_start = selected_slice
  project.procPramArray[ci].state_2 = 0
  sort_dir  =  project.procPramArray[CI].sort_dir
  project.procPramArray[CI].sort_dir = 0
  mas_load_state_2
  image = (*project.dataArray[ci].state2)[*,*]
  mas_windowing, image
  pimage = ptr_new(image,/no_copy)
   
  if (not ptr_valid(pImage)) then return
  
  img_dims = (size(*pImage, /dimensions))[0:1]
  
  if not (ptr_valid(project.roi.pROIs_volume[selected_slice])) then begin
    
    xroi, *pimage , regions_out=oRoi, /block, /modal ,group = event.top
    project.roi.xroi_running = 0

    if (project.procpramarray[ci].no_transform_roi eq 0) then begin
        mas_roi_transform, oRoi, img_dims, /native
    endif
  endif else begin
    oRoi =  *project.roi.pROIs_volume[selected_slice]  ;change to pROIs_volume
    if (project.procpramarray[ci].no_transform_roi eq 1) then begin
    
        xroi, *pImage, regions_in = oRoi, regions_out=oRoi, /block, /modal ,group = event.top
        project.roi.xroi_running = 0
        
    endif else begin
    
        mas_roi_transform, oRoi, img_dims, /current

        xroi, *pImage, regions_in = oRoi, regions_out=oRoi, /block, /modal ,group = event.top
        project.roi.xroi_running = 0

        mas_roi_transform, oRoi, img_dims, /native
    
    endelse
   endelse
    
  project.roi.pROIs_volume[selected_slice] = ptr_new(oRoi)
    
  project.procPramArray[ci].state_2 = 0
  project.procPramArray[CI].sdim_start = main_window_slice
  project.procPramArray[CI].sort_dir = sort_dir
  widget_control, dce_volumes_base, set_uvalue=info
end

; Subroutine name: dce_volumes_SE
; Created by: Garrett Astary
; Calling Information:
;
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Calculates the signal enhancement for the entire volume
;
; Editing Information:
; Edited by Magdoom (12/1/2014)
; Added a selective curve fitting functionality

function dce_volumes_SE
  COMPILE_OPT IDL2
  COMMON scan_data
  COMMON common_widgets
  CI = project.CI
  
  if ptr_valid(project.procpramArray[ci].array_select) then s = *project.procpramArray[ci].array_select $
  else s = indgen(project.imndarray[ci].adim)  
  
  n_Pre = project.imndArray[CI].n_Pre
  total_adim = (size(s,/Dimension))[0]   
  n_Pst = total_adim - n_Pre
  fdim = project.ImndArray[ci].fdim 
  pdim = project.ImndArray[ci].pdim 
  sdim = project.ImndArray[ci].sdim
  
  pre_avg = fltarr(fdim, pdim, sdim)
  for i=0, n_Pre-1 do begin
    pre_avg = pre_avg + (*project.dataArray[ci].state1)[*,*,*,i]
  endfor
  pre_avg = pre_avg/n_Pre
  
  SE_volume = fltarr(fdim, pdim, sdim, n_pst)
  for j=0, n_pst-1 do begin
    SE_volume[*,*,*,j] = ((*project.dataArray[ci].state1)[*,*,*,s[j+n_pre]]-pre_avg)/pre_avg*100.
  endfor
  
return, SE_volume
end

; Subroutine name: dce_volumes_segment
; Created by: Garrett Astary
; Calling Information:
;
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Calculates the tracer distribution volume of a given ROI for all chosen array points
;
; Editing Information:
; Edited by Magdoom (12/1/2014)
; 1) Added a selective curve fitting functionality
; 2) Fixed Bug  : Looping over all slices (slices without ROIs drawn resulted in error) 
;    Solution : Slice loop removed

function dce_volumes_segment, SE_volume, threshold, islice
  COMPILE_OPT IDL2
  COMMON scan_data
  COMMON common_widgets
  CI = project.CI
  
  if ptr_valid(project.procpramArray[ci].array_select) then s = *project.procpramArray[ci].array_select $
  else s = indgen(project.imndarray[ci].adim)  
  
  n_Pre = project.imndArray[CI].n_Pre
  total_adim = (size(s,/Dimension))[0]     
  n_Pst = total_adim - n_Pre
  fdim = project.ImndArray[ci].fdim 
  pdim = project.ImndArray[ci].pdim 
  sdim = project.ImndArray[ci].sdim
  ffov = project.ImndArray[ci].f_fov 
  pfov = project.ImndArray[ci].p_fov
  sfov = project.ImndArray[ci].s_fov
  
  Vd = fltarr(n_pst+1)
  for post=1, n_pst do begin
     regions = (*project.roi.pROIs_volume[islice])
     mask = regions[0] ->ComputeMask(dimensions = [fdim, pdim], MASK_RULE = 2)
     mask = mask/255.
     data = reform(SE_volume[*,*,islice,s[post-1]])
     masked_data = data*mask
     segment_index = where(masked_data ge threshold*100.)
     sz_segment_index = size(segment_index)
       if sz_segment_index[0] ne 0 then begin
        Vd[post] = Vd[post] + sz_segment_index[1]
       endif else begin
        Vd[post] = Vd[post] + 0.0
       endelse 
  endfor
  
  rres = ffov*10./fdim
  pres = pfov*10./pdim
  sres = sfov*10./sdim
  Vd = Vd*rres*pres*sres  
  return, Vd
end

; Subroutine name: dce_volumes_display
; Created by: Garrett Astary
; Calling Information:
;
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Plots tracer distribution volume in a ROI as a function of time
;
; Editing Information:
; Edited by Magdoom (12/1/2014)
; Added a selective curve fitting functionality

pro dce_volumes_display, info
  COMPILE_OPT IDL2
  COMMON scan_data
  COMMON common_widgets
  CI = project.CI
  threshold = info.threshold
  if ptr_valid(project.procpramArray[ci].array_select) then s = *project.procpramArray[ci].array_select $
  else s = indgen(project.imndarray[ci].adim)  
  
  signal_enhancement_volumes = dce_volumes_SE()
  Vd_vector = dce_volumes_segment(signal_enhancement_volumes, threshold, info.selected_slice)
  
  time_array = (*project.imndArray[CI].time.array)[s]
  n_Pre = project.imndArray[CI].n_Pre
  total_adim = (size(s,/Dimension))[0]   
  n_Pst = total_adim - n_Pre
  time = fltarr(n_pst +1)
  final_pre = time_array[n_Pre-1]
  for i=0, n_pst do begin
    time[i] = time_array[i+2] - final_pre
  endfor
  
  IPLOT, time, Vd_vector , $
                title = 'Distribution Volume', $
                IDENTIFIER=6, $
                NAME='Vd', $
                SYM_size = 2, $
                SYM_INDEX = 1, $
                XTitle = 'Time (sec)', $
                YTitle = 'Vd (mm^3)', $
                thick = 2.0,$
                COLOR = FSC_Color('Black',/Triple, /Row)  
  
  
end

; Subroutine name: dce_volumes_SE
; Created by: Garrett Astary
; Calling Information:
;
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Event handler for the DCE volumes GUI
;
; Editing Information:
; Edited by Magdoom (11/30/2014)
; Fixed Bug : Widget variables not passed through the GUI widget. 
; Solution : Neccessary variables are now defined in common_widgets in mas.pro 

pro mas_dce_volumes_gui_event, event
  COMPILE_OPT IDL2
  COMMON scan_data
  COMMON common_widgets
  
  widget_control, event.top, get_uvalue=info, /no_copy
  CI = project.CI
  
  CASE event.id OF
    
    dce_volumes_threshold_text: BEGIN
      widget_control, dce_volumes_threshold_text, get_value=thr_temp
      info.threshold = thr_temp*.01
    END
    
    dce_volumes_ROI_slider: BEGIN
      widget_control, dce_volumes_ROI_slider, get_value = slice
      info.selected_slice = slice
    END
    
    dce_volumes_ROI_edit_button: BEGIN
      widget_control, dce_volumes_ROI_slider, get_value = slice
      info.selected_slice = slice
      dce_volumes_ROI_edit, info, event
    END
    
    dce_volumes_display_button: BEGIN
      widget_control, dce_volumes_threshold_text, get_value=thr_temp
      info.threshold = thr_temp*.01
      dce_volumes_display, info
    END
  
  ENDCASE
widget_control, event.top, set_uvalue=info 
end

; Subroutine name: mas_dce_volumes_gui
; Created by: Garrett Astary
; Calling Information:
;
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: GUI for DCE volumes which plots the tracer distribution volume (mm^3) in a ROI drawn on a slice over time. 
;
; Editing Information:

pro mas_dce_volumes_gui

  HEAP_GC
  COMPILE_OPT IDL2
  COMMON scan_data
  COMMON common_widgets
  CI = project.CI
 
  dce_volumes_base = widget_base(column=1, TITLE='Calculate DCE-MRI Vd', xoffset=420, /base_align_center)
    
  dce_volumes_threshold_base = widget_base(dce_volumes_base, column = 2)
  dce_volumes_threshold_label = widget_label(dce_volumes_threshold_base, value = 'Threshold (Signal Enhancement %):')
  dce_volumes_threshold_text = widget_text(dce_volumes_threshold_base, value='10', xsize=5, /editable, sensitive = 1)   
    
  dce_volumes_ROI_label_base = widget_base(dce_volumes_base)
  dce_volumes_ROI_label = widget_label(dce_volumes_ROI_label_base, value = 'Segmentation ROI')
  
  n_slices = project.imndArray[CI].slices
  dce_volumes_ROI_base = widget_base(dce_volumes_base, column=2)
  dce_volumes_ROI_slider_base = widget_base(dce_volumes_ROI_base, row = 2)
  dce_volumes_ROI_slider = widget_slider(dce_volumes_ROI_slider_base, minimum = 0, maximum = n_slices-1)
  dce_volumes_ROI_slider_label = widget_label(dce_volumes_ROI_slider_base, value = 'Slice')
  dce_volumes_ROI_edit_base = widget_base(dce_volumes_ROI_base)
  dce_volumes_ROI_edit_button = widget_button(dce_volumes_ROI_edit_base, value = 'Edit ROI') 
  
  dce_volumes_action_base = widget_base(dce_volumes_base, column = 2)
  dce_volumes_display_button = widget_button(dce_volumes_action_base, value = 'Display Vd vs. t')
 
  widget_control, dce_volumes_base, /realize
  info = {threshold: 0.2, selected_slice:0, roi_flag: fltarr(n_slices)}
 
  widget_control, dce_volumes_base, set_uvalue=info, /no_copy

  xmanager, 'mas_dce_volumes_gui', dce_volumes_base

end