; IDL routine to extract connected regions in an image


; Subroutine name: region_extract_gui
; Created by: Magdoom Kulam
; Calling Information:
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Main region extract GUI
;
; Editing Information:

pro region_extract_gui

  common scan_data
  common params, ROInum, SDM
  
  if (xregistered('region_extract_gui') ne 0) then return
  topbase = widget_base(title = 'Extract Region : Scan #' + strtrim(project.ci+1,2) ,tlb_frame_attr=1,/grid_layout)
  mainbase = widget_base(topbase,row=3)
  
  ROI_base_RE = widget_base(mainbase, column = 2)
  ROI_label_RE = widget_label(ROI_base_RE, value='Select ROI')
  
  ROIlist_val = make_array(n_elements(*project.roi.pROIs[1]),/string)
  for i = 0, n_elements(*project.roi.pROIs[1])-1 do begin
    (*project.roi.pROIs[1])[i]->Getproperty,name = temp
    ROIlist_val[i] = temp
  endfor
  ROIdroplist_RE = widget_droplist(ROI_base_RE,value= ROIlist_val, uvalue= 'ROIregion', /dynamic_resize)
  
  SD_base = Widget_Base(mainbase, /ROW, /ALIGN_CENTER)
  SD_RE = CW_FSlider(SD_base,title = 'Standard Deviation Multiplier', uvalue = 'SDM', value=0.0, /EDIT)
    
  disp_button_base = widget_base(mainbase)
  disp_button_RE = widget_button(disp_button_base, value = 'Display', uvalue = 'Display',/no_release, xsize = 82)
  
  ROInum = 0
  SDM = 1.0
   
  widget_control, topbase, /realize
  xmanager, 'region_extract_gui', topbase, cleanup='region_extract_gui_cleanup',/NO_BLOCK

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

pro region_extract_gui_event, event

common scan_data
common params

;; To catch error ahead and prevent crashing
;catch,error_status
;IF (error_status NE 0) THEN BEGIN
;   help, calls= trace_back
;   dummy = DIALOG_MESSAGE(['Please send a copy of this error message to the developer of this program.', $
;   !ERROR_STATE.msg,trace_back], /ERROR, TITLE='Error in WID_BASE_MAIN_event')
;   widget_control, event.top, /destroy
;   RETURN
;ENDIF

widget_control, event.top, get_uvalue=info
widget_control, event.id, get_uvalue = widget
CI = project.ci

case widget of

'ROIregion': ROInum = event.index
             
'SDM'      : SDM = event.value
             
'Display'  : region_extract,ROInum,SDM

else       : return

endcase
end

pro region_extract,region_num, SD_mult
common scan_data

I = *project.dataarray[project.ci].state2
N = size(I,/Dimension)
ROImask = mas_roi_get_current_mask(region_num = region_num)
roiPixels = where(*ROImask ne 0)

if n_elements(N) eq 3 then begin
  for j = 0,N[2]-1 do begin
    newROIPixels = REGION_GROW(I[*,*,j], roiPixels, STDDEV_MULTIPLIER = SD_mult,/nan)
    ; Color the grown region, and draw the new image.
    img1 = I[*,*,j]
    img1[newROIPixels] = 255b
    imgTrue = REBIN(I[*,*,j], N[0], N[1], 3)
    imgTrue[*,*,1] = img1
    iimage,imgTrue, layout = [sqrt(N[2]),1+sqrt(N[2]),j+1], aspect_ratio = N[1]/N[0],/current
  endfor
   
endif else begin
  
  newROIPixels = REGION_GROW(I, roiPixels, STDDEV_MULTIPLIER = SD_mult,/nan)

  iimage, *ROImask,layout = [2,1,1], aspect_ratio = N[1]/N[0], title = 'Seed region'
  iimage, I, /current, layout = [2,1,1], transparency = 50, aspect_ratio = N[1]/N[0]

  ; Color the grown region, and draw the new image.
  img1 = I
  img1[newROIPixels] = 255b
  imgTrue = REBIN(I, N[0], N[1], 3)
  imgTrue[*,*,1] = img1
  iimage,imgTrue, /current, layout = [2,1,2], aspect_ratio = N[1]/N[0], title = 'Grown region'
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

pro region_extract_gui_cleanup, topbase

  widget_control, topbase, get_uvalue=state
  if (ptr_valid(state)) then ptr_free, state
  return

end