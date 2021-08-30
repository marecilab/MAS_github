; Subroutine name: select_data_dce
; Created by: Magdoom Kulam
; Calling Information:
;
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Data selection GUI which allows which arrays to fit for DCE analysis for the chosen scan
;
; Editing Information:

pro select_data_dce
common scan_data
select_data_gui,index = project.ci
end

;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved

; *****************************************************************
; *                                                               *
; * procedure tofts_fix, t, p, f, [pder]                          *
; *                                                               *
; * description:                                                  *
; *   triexponential function for curvefitting of signal          *
; *   enhancement data, using fixed plasma curve parameters from  *
; *   common block 'plasma_params'                                *
; *                                                               *
; * input via common block 'plasma_params':                       *
; *   A, B, m1, m2                                                *
; *                                                               *
; * input arguments:                                              *
; *   t    : vector of values for the independent variable        *
; *   p    : vector of function parameters                        *
; *   f    : reference for output of dependent variable vector    *
; *   pder : array of partial derivatives with respect to         *
; *          parameters in a (optional)                           *
; *                                                               *
; * created: April 24, 2001,                                      *
; *          by Gustav Persson                                    *
; *                                                               *
; * modified: 7-22-03 Ty Black hard coded plasma params.          *
; * HS - 20061004. Fixed spelling mistakes and commenting.        *
; *****************************************************************

pro tofts_fix, t, p, f, pder
  ;COMMON plasma_params ; A, B, m1, m2
  print, 'tofts_fix'
  A=15
  B=3.8
  m1 = 0.016
  m2 = 0.00033

  k = p[0]
  v = p[1]

  kv = k/v

  ekv = exp(-kv * t)
  em1 = exp(-m1 * t)
  em2 = exp(-m2 * t)

  f = k * ( A * (em1 - ekv) / (kv - m1) + $
            B * (em2 - ekv) / (kv - m2) )

; Calculate partial derivative array
; if variable for storage is supplied
;---------------------------------------->

  if N_PARAMS() ge 4 then begin
    dfdk = A * ( 1.0 / (kv - m1) - kv / (kv - m1)^2 ) * $
               (em1 - ekv) + $
           B * ( 1.0 / (kv - m2) - kv / (kv - m2)^2 ) * $
               (em2 - ekv) + $
           ( A / (kv - m1) + B / (kv - m2) ) * kv * t * ekv
    dfdv = kv^2 * ( A * (em1 - ekv) / (kv - m1)^2 + $
                    B * (em2 - ekv) / (kv - m2)^2 - $
                    ( A / (kv - m1) + B / (kv - m2) ) * t * ekv )

    pder = [[dfdk],[dfdv]]

  endif

;<----------------------------------------

end

pro DCE_ROI_CURVE_FIT_GUI
  common scan_data

  if xregistered('DCE_ROI_CURVE_FIT_GUI') ge 1 then return
  top_base = widget_base(title = 'DCE ROI Curves' ,tlb_frame_attr=1,/grid_layout)
  main_base = widget_base(top_base,row = 5)
  button_val = ['Standard','Normalized']
  dat_button_grp = cw_bgroup(main_base,column = 2,button_val, uvalue = 'Disp_type', /exclusive, button_uvalue = button_val, frame = 2, set_value=0,/no_release)
  roi_set_base = widget_base(main_base, row=2,sensitive = 0)
  roi_set_label = widget_label(roi_set_base, value = 'Select ROI Set')
  roi_set_droplist = widget_droplist(roi_set_base,value = *project.roi.pdisplay_names, uvalue='roi_set_droplist')
  roi_label = widget_label(roi_set_base, value='ROI')
  roi_droplist = widget_droplist(roi_set_base,value='      ', uvalue='roi_droplist', /dynamic_resize)
  smooth_base = widget_base(main_base, row =1,/nonexclusive)
  smooth_enable = widget_button(smooth_base, uvalue = 'Smooth', value='Smooth time series data')
  misc_buttons_txt = ['Display', 'Close']
  misc_button_grp = cw_bgroup(main_base,column = 2, misc_buttons_txt, uvalue = 'misc_buttons', button_uvalue = misc_buttons_txt, /no_release, space = 50)
  
  infoptr = ptr_new({dat_button_grp:dat_button_grp, misc_button_grp:misc_button_grp,roi_set_base:roi_set_base,roi_set_droplist:roi_set_droplist,roi_droplist:roi_droplist,smooth_flag:0})
  widget_control, top_base, set_uvalue=infoptr
  widget_control, top_base, /realize
  xmanager, 'DCE_ROI_CURVE_FIT_GUI', top_base, cleanup='DCE_ROI_CURVE_FIT_GUI_CLEANUP',/no_block

end

; Subroutine name: select_data_gui_event
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Event handler for the data selection GUI
;
; Editing Information:

pro DCE_ROI_CURVE_FIT_GUI_EVENT, event

  common scan_data

  widget_control, event.top, get_uvalue = infoptr
  info = *infoptr
  widget_control, event.id, get_uvalue = widget
 
 case widget of
  'misc_buttons' : begin
                   case event.value of
                   'Display' : begin
                               widget_control, info.dat_button_grp, get_value = disp_type
                               roi_set_index = widget_info(info.roi_set_droplist, /DROPLIST_SELECT)
                               roi_index = widget_info(info.roi_droplist, /DROPLIST_SELECT)
                               if disp_type eq 0 then DCE_ROI_CURVE_FIT, smooth_flag = info.smooth_flag else DCE_ROI_CURVE_FIT, smooth_flag = info.smooth_flag, type = 2, roi_num = roi_index, roi_set_num = roi_set_index
                               end
                    'Close'   : widget_control, event.top,/destroy
                     endcase
                     end
                     
  'Disp_type'      : begin
                     widget_control, event.id, get_value = disp_type
                     if disp_type eq 0 then widget_control,info.roi_set_base, sensitive = 0 else widget_control,info.roi_set_base, sensitive = 1                                     
                     end
                     
 'roi_set_droplist': begin
                     roi_set_index = widget_info(info.roi_set_droplist, /DROPLIST_SELECT)
                     if ptr_valid(project.roi.pRois[roi_set_index]) then rois = (*project.roi.pRois[roi_set_index]) else return
                     if (size(rois))[0] EQ 0 THEN num_roi = 1 ELSE num_roi =  (size(rois))[1]
                        roi_names = strarr(num_roi)
                        FOR temp=0 , num_roi-1 DO BEGIN
                          rois[temp] -> GETPROPERTY, name= name
                          roi_names[temp] = name
                        END
                        widget_control, info.roi_droplist, set_value=roi_names
                      end 
   'Smooth'        : (*infoptr).smooth_flag = event.select                   
   else            : return
   
   endcase
end

; Subroutine name: select_data_gui_cleanup
; Created by: Magdoom Kulam
; Calling Information:
;
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: Cleanup procedure before exiting the data selection GUI
;
; Editing Information:

pro DCE_ROI_CURVE_FIT_GUI_CLEANUP, top_base

  widget_control, top_base, get_uvalue=state
  if ptr_valid(state) then ptr_free, state
  return

end

; Subroutine name: DCE_ROI_CURVE_FIT
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; This procedure will do calculations based upoin the roi's
; drawn by the user.
; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro DCE_ROI_CURVE_FIT, smooth_flag = smooth_flag, type = type, roi_num = roi_num, roi_set_num = roi_set_num
    COMPILE_OPT IDL2

    ;bring in the appropriate global varables
    COMMON scan_data
    CI = project.ci
    sdimStart = project.procPramArray[CI].sdim_Start
    adimStart = project.procPramArray[CI].adim_Start
    if ptr_valid(project.procpramArray[ci].array_select) then adim = (size(*project.procpramArray[ci].array_select,/Dimension))[0] else adim = project.imndarray[ci].adim
    n_Pre = project.imndArray[CI].n_Pre
    time_Start = project.imndArray[ci].time.min

    ;check to see it the data set has a valid time array.
    IF PTR_VALID(project.imndArray[ci].time.array) EQ 0    THEN RETURN
    time_End = project.imndArray[ci].time.max
    if ptr_valid(project.procpramArray[ci].array_select) then time_Array = (*project.imndArray[ci].time.array)[*project.procpramArray[ci].array_select] $
    else time_Array = (*project.imndArray[ci].time.array) 
    time_length= (*project.imndArray[ci].time.length)

    time_Zero =  time_Array[n_pre] - time_length[n_pre]/2
    ;move the zero time
    time_Array -= time_zero

    if time_array[0] gt 0.0 then begin
       void = dialog_message(['There are no times less than zero.','pre image times have the same time stamp'],/error)
       return
    end

    multi_ROI_graphs = project.procPramArray[project.ci].multi_ROI_graphs

    ;check to make sure that the user didn't select the first ROI from the ROI list
    ;which has special meaning.
    if project.roi.ci eq 0 then begin
        void= dialog_message(['Please select a different item in the ROI list',$
                        'The first item in the ROI list is for creating a new ROI'],/info)
        return
    end

    ;get the data object from the roi program.
    if not(ptr_valid(project.roi.pROIs[project.roi.ci])) then begin
       void= dialog_message(['No ROI selected','Please select a different ROI from ROI list'],/info)
       return
    end
    ;now that we know that the pointer is valid try and extract the object
    ;at then end of the pointer
    regions = (*project.roi.pROIs[project.roi.ci])

    ;check to make sure the object regions is valid
    if not(total(obj_valid(regions))) then begin
       ;update_status_bar, 'No region returned from ROI tool'
       return
    end

    ;make sure the image is in state1
    mas_load_state_1

    ;copy the image locally. this image will have the roi drawn on it
    ;process the image the way the user has set up.
    p_image = ptr_new((*project.dataArray[ci].state1)[*,*,sdimStart,adimStart])

;; FIXME ROI
    if (project.procpramarray[ci].no_transform_roi eq 1) then begin
        mas_zoom , p_image
        mas_rotate_flip , p_image
    endif

    ;clear the status bar
    update_status_bar, ''

    ;figure out how many roi's there are
    if     (size(regions))[0] eq 0 then numRoi = 1 $
    else numRoi =  (size(regions))[1]
    ;print, 'numRoi=',numRoi

    ;we need an array to hold the roi names.
    name_array = strarr(numRoi)

    ;create the colors to loop through when over plotting
	;color_names = ['Red','Yellow','Blue','Green','Brown','Cyan','Purple','Maroon','Gray','Orange']
    color_names = ['Red','Blue','Forest Green','Orange','Gray','Maroon','Brown','Purple','Dark Gray','Hot Pink']
    length_color_names = (size(color_names))[1]


    ; make a pointer array to hold all the masks
    mask = ptrarr(numRoi)
    for temp=0 , numRoi-1 do begin
       mask[temp] = ptr_new(regions[temp] -> ComputeMask( dimensions =[ (size((*p_image)))[1], (size((*p_image)))[2] ] , MASK_RULE=2))
       regions[temp] -> GETPROPERTY, name= name
       name_array[temp] = name
    end

    ;after we know how many roi's there are we can now make an array to hold all the info so that we can
    ;display it later on after the data is calculated.the array is 1 larger than what it needs to be
    ;so that the first entry will be the independent points
    ROI_Data = fltArr(adim,numRoi)

    ;now we create an array to hold the fit parameters and chisquare for each ROI
    ; Fit_Data = fltArr(numRoi,4)

    ;create an string array that will hold all the info to display to the screen
    ;1 for the time array and 2*numRoi, 1 line for roi name and 1 for roi data.
    sData = strarr(2+numRoi*2)
    sData[0] = string('Time')
    sData[1] = strjoin(strTrim(time_array,2),string(09b))
    sCounter=2
    ;make an array to hold the average intensity for each adim
    avgIntensityInROI = fltarr(adim)

    maskCounter = 0

; Adding counters for cycling colors and markers.
    sym_counter = 2
    marker_counter = 0
    line_style_counter = 2

    for maskCounter= 0 , numRoi-1 do begin

       ;loop through the adim and where the roi mask is
       ; average that data and save it into the avgIntensityInROI matrix
       for ii=0, adim-1 do begin
         ;we have to bring state1 back in b/c we byte scaled it earlier
         ; and we need to extract data from the ffted image and not the bytescalled
         ; also the mask was drawn on the zoomed image so we have to zoom the image also
          
         if ptr_valid(project.procpramArray[ci].array_select) then begin
         c = (*project.procpramarray[ci].array_select)[ii]
         p_image = ptr_new((*project.dataArray[ci].state1)[*,*,sdimStart,c]) 
         endif else p_image = ptr_new((*project.dataArray[ci].state1)[*,*,sdimStart,ii]) 

;; FIXME ROI
        if (project.procpramarray[ci].no_transform_roi eq 1) then begin
            mas_zoom , p_image
            mas_rotate_flip , p_image
        endif

         ;image_statistics will calculte the mean for me automatically
         ; then i just save that to the avg intensity roi
         IMAGE_STATISTICS , (*p_image), mask=*(mask[maskCounter]) $
                          , SUM_OF_SQUARES=image_sum_of_squares $
                          , VARIANCE=image_var $
                          , Count=image_count $
                          , Mean = d
       avgIntensityInROI[ii] = d                   
       end

       ROI_Data[*,maskCounter] = avgIntensityInROI
      
       ;calculate the average intensity for the pre image's roi
       average_intensity_pre = mean(ROI_Data[where(time_array lt 0.0),maskCounter])

       ;calculate the average change vs the pre images
       for jj=0, adim-1 do ROI_Data[jj,maskCounter] = ((ROI_Data[jj,maskCounter] - average_intensity_pre)/average_intensity_pre)*100
       if smooth_flag eq 1 then begin
         filter_size = 4
         ROI_Data_smooth  = dblarr(adim)
         ROI_Data_smooth[0:filter_size-2] = ROI_Data[0:filter_size-2,maskCounter]
         FOR t = filter_size-1,adim-1 DO ROI_Data_smooth[t] = mean(ROI_Data[t-filter_size+1:t,maskCounter],dimension = 1)  ; Moving average filter
         ROI_Data[*,maskCounter] = ROI_Data_smooth
       endif
   end
   for maskCounter= 0 , numRoi-1 do begin  
       ;plot single or multiple roi graphs

       iPlot_identifier = 6
       if keyword_set(type) then begin
        (*project.roi.pROIs[roi_set_num])[roi_num]->Getproperty,name = ROIname
        plot_title = 'DCE_ROI data normalized by ' + ROIname
       endif else plot_title = 'DCE_ROI data'
       if smooth_flag eq 1 then plot_title+=' smoothed by moving average filter' 
       
       if maskCounter eq 0 then begin
         if keyword_set(type) then begin
           if type eq 2 then plot_data = 100*ROI_Data[*,maskCounter]/max(abs(ROI_Data[*,roi_num])) 
         endif else  plot_data = ROI_Data[*,maskCounter]
               h = PLOT(time_array, plot_data , $
                title = plot_title, $
                IDENTIFIER=iPlot_identifier, $
                NAME=name_array[maskCounter], $
                $;/SCATTER, $
                SYM_size = 0.3, $
                ;SYM_INDEX = 1, $
                XTitle = 'Time (sec)', $
                YTitle = '% enhancement', $
                thick = 2.0,$
                COLOR = FSC_Color('Black',/Triple, /Row))

       end else begin
          if keyword_set(type) then begin
            if type eq 2 then plot_data = 100*ROI_Data[*,maskCounter]/max(abs(ROI_Data[*,roi_num])) 
          endif else  plot_data = ROI_Data[*,maskCounter]
           h = PLOT(time_array ,plot_data , $
                OVERPLOT=iPlot_identifier, $
                $;/SCATTER, $
                name= name_array[maskCounter], $
                SYM_size = 0.3, $
                ;SYM_INDEX = sym_counter, $
                LINESTYLE = line_style_counter,$
                COLOR = FSC_Color(color_names[(maskCounter-1+marker_counter) mod 10], /Triple, /Row), $
                thick = 2.0)

			;if (sym_counter eq 8) then sym_counter = 2 else sym_counter = sym_counter + 1
			if (maskCounter eq length_color_names) then marker_counter = marker_counter + 2
			if (line_style_counter eq 5) then line_style_counter = 1 else line_style_counter = line_style_counter + 1

       end
           
       ;store the name of the roi and the values for it.
       sData[sCounter] = 'Data Values for '+name_array[maskCounter]
       sCounter++
       sData[sCounter] =  strjoin(strTrim(plot_data,2),string(09b))
       sCounter++

    endfor

    display_stats, sData, 'Measure data for each Dynamics ROI'

    if multi_ROI_graphs eq 1 then print,'multi graphs not supported yet'

end


; Subroutine name: mas_dce_roi_init
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
;
; This procedure will allow the user to specifiy an roi and then allow them
; look at the enhancement curves associated with there specied roi's
;
; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_dce_roi_init
    COMPILE_OPT IDL2

    COMMON common_widgets
    COMMON scan_data
    ci = project.ci

    ;bring in the common data
    sdimStart = project.procPramArray[CI].sdim_Start
    adimStart = project.procPramArray[CI].adim_Start
    adim = project.imndArray[CI].adim

    ;make sure the image is in state1
    mas_load_state_1


    ;copy the image locally. this image will have the roi drawn on it
    ;process the image the way the user has set up.
    p_image = ptr_new((*project.dataArray[ci].state1)[*,*,sdimStart,adimStart])
    mas_windowing , *p_image
    mas_zoom , p_image
    mas_rotate_flip , p_image

    ;tell the user to draw ROI's and then exit the ROI
    update_status_bar, 'Draw ROIs then click Process'

    tlb=widget_base(/column)
    ;mas_xroi, (*p_image), outgroup=og,apptlb=tlb, group=tlb

    b=widget_button(tlb,value='Process', event_pro='DCE_ROI_CURVE_FIT', uvalue=og)

    tlb2=widget_base(tlb,/row)
    GRAPH_BASE = Widget_Base(tlb2, UNAME='GRAPH_BASE  ' ,COLUMN=1 ,/NONEXCLUSIVE)

    b= Widget_Button(GRAPH_BASE, UNAME='multi_Graphs' ,/ALIGN_LEFT ,VALUE='Multiple Graphs' $
       , event_pro='multi_graphs_toggle_event'   )
    if project.procPramArray[CI].multi_ROI_graphs eq 1 then WIDGET_CONTROL, /SET_BUTTON, b


end



; Subroutine name: mas_signal_enhancement_event
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_signal_enhancement_event , Event
    COMPILE_OPT IDL2

    COMMON common_widgets
    COMMON scan_data

    catch,error_status
    if (error_status ne 0) then begin
       help, calls= trace_back

       dummy = DIALOG_MESSAGE(['Please send a copy of this error message to the developer of this program.', $
         !ERROR_STATE.msg,trace_back], /ERROR, $
           TITLE='Error in WID_BASE_MAIN_event')


       RETURN
    endif

    ci = project.ci
    wWidget =  Event.top

    case Event.id of
      SE_WINDOW_DISPLAY: begin
           if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then begin

          numPre = project.imndArray[ci].n_Pre
          minRef =  project.procPramArray[ci].min_Ref
          intensity_max_se = project.procPramArray[ci].intensity_max_se

          ;there is still a chance that the user could send in data with out
          ;any pre images so if there is no preimages then quit.
          if numPre eq 0 then begin
            update_status_bar,'scan does not contain any PRE images'
            return
          end

          ;bring in the data that we need.
          sdim_start = project.procPramArray[CI].sdim_start
          if ptr_valid(project.procpramArray[ci].array_select) then begin
          ind = *project.procpramArray[ci].array_select
          image = reform((*project.dataArray[ci].state1)[*,*,sdim_start,ind])
          endif else image = reform((*project.dataArray[ci].state1)[*,*,sdim_start,*])
          
         ;now we need to generate a mask for our data set.
         ;sence we have multiple slices we have to generate it from the begining
         ;and loop through then freq then phase then slice. that way we don't
         ;have to recreate them unncessarly.
         IF project.procPramArray[CI].remove_background_active EQ 1 THEN BEGIN
          dims = size(image,/DIMENSIONS )

          mask = bytarr(dims[0], dims[1])
          key_adim = project.procPramArray[CI].key_adim

          p_image_mask = PTR_NEW(REFORM((*project.dataArray[CI].state1)[*,*,sdim_start,key_adim]))
          mask[*,*]= *(generate_background_mask( p_image_mask))

          FOR ii=0, dims[2]-1 DO $
              image[*,*,ii] *= mask
         END

          sizeOfImage = size(image)

         ;use the roi mask
         FOR ii=0, sizeOfImage[3]-1 DO BEGIN
          p_image_temp = ptr_new(image[*,*,ii])
          mas_roi_mask, p_image_temp
          image[*,*,ii] = *p_image_temp
         END

          fdim = sizeOfImage[1]
          pdim = sizeOfImage[2]
          sdim = sizeOfImage[3]

          SEData = FltArr(fdim, pdim, sdim)
          mean_pre_flt = FltArr(fdim, pdim)

          update_status_bar,'Calculating mean pre-injection signal'

          ; Calculate mean pre-injection signal
          ;========================================>

          For rr = 0, fdim-1 Do $
              For pp = 0, pdim-1 Do begin
                  mean_pre_flt[rr,pp] =Mean(image[rr,pp,0:numPre-1])
              end
          ;<========================================



          minRef = minRef * max(mean_pre_flt)
          print, 'minref after', minRef
          update_status_bar,'Calculating signal enhancement data'
          For ii = 0, fdim-1 Do $
              For jj = 0, pdim-1 Do $
                  if mean_pre_flt[ii,jj] GT minRef then $
                      SEData[ii,jj,*] = ( (image)[ii,jj,*] / mean_pre_flt[ii,jj] ) - 1 $
                  else $
                      SEData[ii,jj,*] = 0


          update_status_bar, 'Intenisfying signal enhancement image...Start'

          pSEData=ptr_new(bytscl(SEData, max=intensity_max_se))
          mas_windowing , image

          pImage = ptr_new(image)
          
          signal_Enhancement = ptr_new(bytarr((size(SEData))[1]*2,(size(SEData))[2],sdim))

          (*signal_Enhancement)[0:(size(SEData))[1]-1,*,*] = *pimage
          (*signal_Enhancement)[fdim:*  ,*,*] = *pSEData

          project.dataArray[ci].signal_Enhancement = signal_Enhancement

          mas_movie_se

         endif
        end
       SE_MAX_REBIN: begin
           if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_SLIDER' ) then begin
               WIDGET_CONTROL , SE_MAX_REBIN , GET_VALUE=temp
          print, temp
               project.procPramArray[ci].intensity_max_se = temp


           endif

        end

        SE_MINIMUM_REF: begin
           if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_SLIDER' ) then begin
               WIDGET_CONTROL , SE_MINIMUM_REF , GET_VALUE=temp
          ;print, temp
               project.procPramArray[ci].min_Ref = float(float(temp)/100.0)


           endif

        end

    endcase
end


; Subroutine name: mas_signal_enhancement
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_signal_enhancement
    COMPILE_OPT IDL2
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets

    ci = project.ci

    mas_load_state_1

    update_status_bar,'Set signal enhancement parameters.'

    SE_WINDOW_BASE = widget_base(/COLUMN, TITLE='Signal Enhancement', UVALUE = SE_WINDOW_BASE ,SCR_XSIZE=200)

    SE_MINIMUM_REF = widget_slider(SE_WINDOW_BASE, MINIMUM=0, FRAME=0, MAXIMUM=100, $
                       TITLE = 'Set fit threshold' ,   $
                       UNAME='SE_MINIMUM_REF', $
                       VALUE=FIX(project.procPramArray[ci].min_Ref*100)  )
    W = Widget_Label(SE_WINDOW_BASE, VALUE='')

    SE_MAX_REBIN = widget_slider(SE_WINDOW_BASE, MINIMUM=0, FRAME=0, MAXIMUM=20 ,$
                       TITLE = 'Set max fractional enhancement' ,   $
                       UNAME='SE_MAX_REBIN',$
                       VALUE= project.procPramArray[ci].intensity_max_se )

    W = Widget_Label(SE_WINDOW_BASE, VALUE='')
                        
    SE_WINDOW_DISPLAY= widget_button(SE_WINDOW_BASE,value='Display'$
              ,UNAME = 'SE_WINDOW_DISPLAY')

    SE_LABEL = Widget_Label(SE_WINDOW_BASE, UNAME='SE_LABEL' ,SCR_XSIZE=150, $
                   VALUE='          ')

    widget_control, SE_WINDOW_BASE, /realize

    xmanager, 'mas_signal_enhancement',SE_WINDOW_BASE,/no_block

end



pro mas_sig_enhance
end
