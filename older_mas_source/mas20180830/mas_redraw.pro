;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved

; HS 20061002. All "dimenshions" renamed to dimensions.


; Subroutine name: mas_redraw_GUI
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:
; This procedure will redraw almost all the base widgets and take values from the current
; index and use the common data structure to poplulate the widgets value.

; Editing Information:
    ;Edited by HS 2006/10/02.
    ;Fix spelling mistakes and commenting.


pro mas_redraw_GUI, event
COMPILE_OPT IDL2


COMMON scan_data
COMMON common_widgets

ii = project.ci
jj = project.scan_Open_Flag

if (project.imndarray[ii].image_type eq 99) then begin
    is_imaging = 0
    is_spect   = 1
endif else begin
    is_imaging = 1
    is_spect   = 0
endelse

;first i will update the scan selection dialog
if project.ni gt 0 then begin
    newListStringArray = strarr(project.ni)
    newListStringArray = project.scan_list[0:project.ni-1]
    Widget_control, WID_SCAN_LIST, SET_LIST_SELECT=ii, SET_VALUE= newListStringArray ,SENSITIVE=jj
end else Widget_control, WID_SCAN_LIST , SET_VALUE= '' ,SENSITIVE=jj

;next i will update the imnd labels

widget_control, W_MENU_Trace_Phase, SENSITIVE=jj*is_imaging
widget_control, W_MENU_FILTER,      SENSITIVE=jj*is_imaging

widget_control, widget_info(WID_BASE_MAIN, find_by_uname='LABEL_SF')    , set_value="SF(MHz):"
widget_control, LABEL_SF_DATA , SET_VALUE= STRTRIM(project.imndArray[ii].spect_bf1,1) ,SENSITIVE=jj

if (is_spect eq 1) then begin

    widget_control, widget_info(WID_BASE_MAIN, find_by_uname='LABEL_RECOVTIME'), set_value="AT(ms):"
    widget_control, LABEL_RECOVTIME_DATA , SET_VALUE= STRTRIM(project.imndArray[ii].spect_acq_time,1) ,SENSITIVE=jj
    
    widget_control, widget_info(WID_BASE_MAIN, find_by_uname='LABEL_ECHOTIME'), set_value="NUCLEUS:"
    widget_control, LABEL_ECHOTIME_DATA , SET_VALUE= STRTRIM(project.imndArray[ii].spect_nucleus,1) ,SENSITIVE=jj
     
    widget_control, widget_info(WID_BASE_MAIN, find_by_uname='LABEL_RARE')    , set_value="TR(ms):"
    widget_control, LABEL_RARE_DATA , SET_VALUE = STRTRIM(project.imndArray[ii].spect_rep_time,1) ,SENSITIVE=jj
    
    widget_control, widget_info(WID_BASE_MAIN, find_by_uname='LABEL_FFOV')    , set_value="SW(Hz):"
    widget_control, widget_info(WID_BASE_MAIN, find_by_uname='LABEL_PFOV')    , set_value="PL(1e-6 s):"
    ;widget_control, widget_info(WID_BASE_MAIN, find_by_uname='LABEL_SFOV')    , set_value="NA:"
    widget_control, LABEL_FFOV_DATA , SET_VALUE = STRTRIM(project.imndArray[ii].spect_spectral_width,1) ,SENSITIVE=jj
    widget_control, LABEL_PFOV_DATA , SET_VALUE = STRTRIM(project.imndArray[ii].spect_pulse_length,1) ,SENSITIVE=jj
    widget_control, LABEL_NAVG_DATA , SET_VALUE = STRTRIM(project.imndArray[ii].spect_num_avg,1) ,SENSITIVE=jj
    widget_control, LABEL_SFOV_DATA , SENSITIVE=jj*is_imaging
    
    widget_control, widget_info(WID_BASE_MAIN, find_by_uname='LABEL_SLICE')   , set_value="Matrix:"
    ;; mas stores the total number of elements.
    widget_control, LABEL_SLICE_DATA , SET_VALUE = STRTRIM(project.imndArray[ii].spect_acq_size,1) ,SENSITIVE=jj
    
endif else begin

    widget_control, widget_info(WID_BASE_MAIN, find_by_uname='LABEL_RECOVTIME'), set_value="TR(ms):"
    widget_control, LABEL_RECOVTIME_DATA , SET_VALUE= STRTRIM(project.imndArray[ii].recov_time,1), SENSITIVE=jj
   
    widget_control, widget_info(WID_BASE_MAIN, find_by_uname='LABEL_ECHOTIME'), set_value="TE(ms):"
    widget_control, LABEL_ECHOTIME_DATA , SET_VALUE= STRTRIM(project.imndArray[ii].echo_time,1) ,SENSITIVE=jj
    
    widget_control, widget_info(WID_BASE_MAIN, find_by_uname='LABEL_NAVG')    , set_value="NA:"
    widget_control, LABEL_NAVG_DATA , SET_VALUE= STRTRIM(project.imndArray[ii].n_avg,1) ,SENSITIVE=jj
    
    widget_control, widget_info(WID_BASE_MAIN, find_by_uname='LABEL_RARE')    , set_value="ETL:"
    widget_control, LABEL_RARE_DATA , SET_VALUE = STRTRIM(project.imndArray[ii].rare,1) ,SENSITIVE=jj
    
    widget_control, widget_info(WID_BASE_MAIN, find_by_uname='LABEL_FFOV')    , set_value="FFOV(mm):"
    widget_control, widget_info(WID_BASE_MAIN, find_by_uname='LABEL_PFOV')    , set_value="PFOV(mm):"
    widget_control, widget_info(WID_BASE_MAIN, find_by_uname='LABEL_SFOV')    , set_value="SFOV(mm):"
    widget_control, LABEL_FFOV_DATA , SET_VALUE = STRTRIM(project.imndArray[ii].f_fov*10,1) ,SENSITIVE=jj
    widget_control, LABEL_PFOV_DATA , SET_VALUE = STRTRIM(project.imndArray[ii].p_fov*10,1) ,SENSITIVE=jj
    widget_control, LABEL_SFOV_DATA , SET_VALUE = STRTRIM(project.imndArray[ii].s_fov*10,1) ,SENSITIVE=jj

    widget_control, widget_info(WID_BASE_MAIN, find_by_uname='LABEL_SLICE')   , set_value="SLICE:"
    widget_control, LABEL_SLICE_DATA , SET_VALUE = STRTRIM(project.imndArray[ii].slices,1) ,SENSITIVE=jj
    
endelse


widget_control, LABEL_DIMENSIONS_DATA , SET_VALUE= STRTRIM(project.imndArray[ii].dimensions,1) ,SENSITIVE=jj
widget_control, LABEL_SLICE_ORIENT_DATA  , SET_VALUE=project.imndArray[ii].orientation[0] ,SENSITIVE=jj*is_imaging

if project.imndArray[ii].orientation[0] NE 'Arbitrary' then begin
    temp = project.imndArray[ii].orientation[1:2]
    temp = strjoin (temp,/single,' to ')
end else begin
    temp = project.imndArray[ii].orientation[1]
endelse

;widget_control, LABEL_SLICE_DIR_DATA  , SET_VALUE=temp ,SENSITIVE=jj*is_imaging

widget_control, LABEL_THICK_DATA , SET_VALUE = STRTRIM(project.imndArray[ii].thick,1) ,SENSITIVE=jj*is_imaging
widget_control, LABEL_FDIM_DATA , SET_VALUE = STRTRIM(project.imndArray[ii].fdim,1) ,SENSITIVE=jj
widget_control, LABEL_PDIM_DATA , SET_VALUE = STRTRIM(project.imndArray[ii].pdim,1) ,SENSITIVE=jj*is_imaging
widget_control, LABEL_SDIM_DATA , SET_VALUE = STRTRIM(project.imndArray[ii].sdim,1) ,SENSITIVE=jj*is_imaging

widget_control, LABEL_SCANDATE_DATA , SET_VALUE=project.imndArray[ii].scan_date ,SENSITIVE=jj
widget_control, LABEL_SCANNAME_DATA  , SET_VALUE=project.imndArray[ii].scan_name ,SENSITIVE=jj
widget_control, LABEL_ADIM_DATA , SET_VALUE=STRTRIM(project.imndArray[ii].adim,1) ,SENSITIVE=jj*is_imaging
widget_control, LABEL_NECHO_DATA , SET_VALUE = STRTRIM(project.imndArray[ii].n_echo,1),SENSITIVE=jj*is_imaging

CASE project.imndArray[ii].fat_sup OF
	1: text = 'Yes'
	2: text = 'SPIR'
	ELSE: text = 'No'
ENDCASE

widget_control, LABEL_FATSUP_DATA , SET_VALUE=text ,SENSITIVE=jj*is_imaging

widget_control, SLIDER_DISPLAY_CONTROL ,SENSITIVE=jj

widget_control, SLIDER_DIMENSION_CONTROL ,SENSITIVE=jj

widget_control, DISPLAY_MOVIE ,SENSITIVE=jj*is_imaging
if project.procPramArray[ii].single_Multi_flag eq 0 then $
    widget_control, DISPLAY_MOVIE ,SENSITIVE=0


widget_control, W_MENU_SAVE_TIFF ,SENSITIVE=jj*is_imaging

;check to see if they are running in the virtual machine
IF LMGR(/VM) GT 0 THEN widget_control, W_MENU_SAVE_MPEG ,SENSITIVE=0 $
ELSE widget_control, W_MENU_SAVE_MPEG ,SENSITIVE=jj*is_imaging



widget_control, W_MENU_SAVE_FLT ,SENSITIVE=jj*is_imaging

widget_control, W_MENU_SAVE_SUBSET ,SENSITIVE=jj*is_imaging

widget_control, W_MENU_SAVE_SLICER ,SENSITIVE=jj*is_imaging

widget_control, W_MENU_SAVE_NIFTI ,SENSITIVE=jj*is_imaging

widget_control, W_MENU_SAVE_DICOM, SENSITIVE=jj*is_imaging

widget_control, menu_save_raw_data ,SENSITIVE=jj*is_imaging

;widget_control, W_MENU_SAVE_FID ,SENSITIVE=jj

;widget_control, W_MENU_SAVE_PROJECT ,SENSITIVE=jj

widget_control, W_MENU_REMOVE_SCAN ,SENSITIVE=jj

widget_control, DISPLAY_GO ,SENSITIVE=jj

widget_control, freq_interp_slider ,SENSITIVE=jj*is_imaging , SET_VALUE = project.procPramArray[ii].freq_interp
widget_control, phase_interp_slider ,SENSITIVE=jj*is_imaging , SET_VALUE = project.procPramArray[ii].Phase_interp
widget_control, slice_interp_slider ,SENSITIVE=jj*is_imaging , SET_VALUE = project.procPramArray[ii].Slice_interp
WIDGET_CONTROL ,lock_interp_chk_box ,SET_BUTTON=project.procPramArray[ii].lock_interp, SENSITIVE=jj*is_imaging
widget_control, down_sample_chk_box, SET_BUTTON=project.procPramArray[ii].down_sample,SENSITIVE=jj*is_imaging
widget_control, zpad_chk_box, SET_BUTTON=project.procPramArray[ii].zpad_flag,SENSITIVE=jj*is_imaging


widget_control, x_zoom_slider ,SENSITIVE=jj*is_imaging , SET_VALUE = project.procPramArray[ii].x_zoom
widget_control, y_zoom_slider ,SENSITIVE=jj*is_imaging , SET_VALUE = project.procPramArray[ii].y_zoom


widget_control, SLIDER_PREIMAGES ,SENSITIVE=jj*is_imaging , SET_VALUE = project.imndArray[ii].n_pre $
       ,SET_SLIDER_MAX=project.imndArray[ii].adim-1

widget_control, n_avg_flag_button, SET_BUTTON = project.procPramArray[ii].n_avg_flag ,SENSITIVE=jj
widget_control, phi_unwrap_flag_button, SET_BUTTON = project.procPramArray[ii].phi_unwrap_flag ,SENSITIVE=jj

widget_control, W_MENU_DCE_ROI_CURVES ,SENSITIVE=jj*is_imaging
widget_control, W_MENU_DCE_SCC ,SENSITIVE=jj*is_imaging
widget_control, W_MENU_DCE_VOL ,SENSITIVE=jj*is_imaging

widget_control, W_MENU_SLICER, SENSITIVE=jj*is_imaging

widget_control, W_MENU_volume, SENSITIVE=jj*is_imaging

widget_control, W_MENU_surface, SENSITIVE=jj*is_imaging

widget_control, W_MENU_MIP, SENSITIVE=jj*is_imaging

widget_control, W_MENU_orthoview, SENSITIVE=jj*is_imaging



widget_control, W_MENU_VOLUMETRICS, SENSITIVE=jj*is_imaging

;; Note the curvefit routines were replaced by mas_curvefit_2. I'm leaving
;; this here just in case someone would like to re-enable it later on.
;widget_control, W_MENU_T1SR_SE_ROI , SENSITIVE=jj*is_imaging
;widget_control, W_MENU_T1SR_SE_IMAGE , SENSITIVE=jj*is_imaging
;widget_control, W_MENU_T1SR_GE_ROI , SENSITIVE=jj*is_imaging
;widget_control, W_MENU_T1SR_GE_IMAGE , SENSITIVE=jj*is_imaging
;widget_control, W_MENU_T1IR_SE_ROI , SENSITIVE=jj*is_imaging
;widget_control, W_MENU_T1IR_SE_IMAGE , SENSITIVE=jj*is_imaging
;widget_control, W_MENU_T1IR_GE_ROI , SENSITIVE=jj*is_imaging
;widget_control, W_MENU_T1IR_GE_IMAGE , SENSITIVE=jj*is_imaging
;
;widget_control, W_MENU_T2MSE_ROI , SENSITIVE=jj*is_imaging
;widget_control, W_MENU_T2MSE_IMAGE , SENSITIVE=jj*is_imaging

widget_control, W_MENU_calculate_delta_b_least_squares, SENSITIVE=jj*is_imaging
widget_control, W_MENU_CURVEFIT_2, SENSITIVE=jj*is_imaging

;; Note the curvefit routines were replaced by mas_curvefit_2. I'm leaving
;; this here just in case someone would like to re-enable it later on.
;if  project.imndArray[ii].image_type eq 3 then begin
;    widget_control, W_MENU_ADC_ROI , SENSITIVE=0
;    widget_control, W_MENU_ADC_IMAGE , SENSITIVE=0
;    widget_control, W_MENU_ADC_STRETCH_ROI , SENSITIVE=0
;    widget_control, W_MENU_ADC_STRETCH_IMAGE , SENSITIVE=0
;end else begin
;    widget_control, W_MENU_ADC_ROI , SENSITIVE=jj*is_imaging
;    widget_control, W_MENU_ADC_IMAGE , SENSITIVE=jj*is_imaging
;    widget_control, W_MENU_ADC_STRETCH_ROI , SENSITIVE=jj*is_imaging
;    widget_control, W_MENU_ADC_STRETCH_IMAGE , SENSITIVE=jj*is_imaging
;end

widget_control, W_MENU_IMAGE_STATS , SENSITIVE=jj*is_imaging
widget_control, W_MENU_IMAGE_HISTOGRAM , SENSITIVE=jj*is_imaging
widget_control, W_MENU_SIGNAL_ENHANCEMENT , SENSITIVE=jj*is_imaging
widget_control, W_MENU_ADT_REGRESS , SENSITIVE=jj*is_imaging
widget_control, W_MENU_ALG , SENSITIVE=jj*is_imaging
widget_control, W_MENU_calculate_delta_B_algebra , SENSITIVE=jj*is_imaging
widget_control, W_MENU_CALC_PROTON_DENSITY , SENSITIVE=jj*is_imaging
; FILTERING MENU
widget_control, W_MENU_FREQUENCY_FILTER, SENSITIVE=jj*is_imaging
widget_control, W_MENU_FREQUENCY_FILTER_LOWPASS, SENSITIVE=jj*is_imaging
widget_control, W_MENU_FREQUENCY_FILTER_HIGHPASS, SENSITIVE=jj*is_imaging
widget_control, W_MENU_IMAGE_DOMAIN, SENSITIVE=jj*is_imaging

;if !d.name eq 'WIN' then widget_control, W_MENU_SMOOTH , SENSITIVE=0 $
;else widget_control, W_MENU_SMOOTH , SENSITIVE=jj


widget_control, smooth_none_button, SENSITIVE=jj*is_imaging
if project.procPramArray[ii].smooth_Direction eq 0 then $
    widget_control, smooth_none_button, /SET_BUTTON

widget_control, smooth_button, SENSITIVE=jj*is_imaging
if project.procPramArray[ii].smooth_Direction eq 1 then $
    widget_control, smooth_button, /SET_BUTTON

widget_control, median_button, SENSITIVE=jj*is_imaging
if project.procPramArray[ii].smooth_Direction eq 2 then $
    widget_control, median_button, /SET_BUTTON



widget_control, flip_none_button, SENSITIVE=jj*is_imaging
if project.procPramArray[ii].flip_Direction eq 0 then $
    widget_control, flip_none_button, /SET_BUTTON

widget_control, flip_horz_button, SENSITIVE=jj*is_imaging
if project.procPramArray[ii].flip_Direction eq 1 then $
    widget_control, flip_horz_button, /SET_BUTTON

widget_control, flip_vert_button, SENSITIVE=jj*is_imaging
if project.procPramArray[ii].flip_Direction eq 2 then $
    widget_control, flip_vert_button, /SET_BUTTON


widget_control, deg_0_button , SENSITIVE=jj*is_imaging
if project.procPramArray[ii].rotate_Direction eq 0 then $
    widget_control, deg_0_button , /SET_BUTTON

widget_control, deg_90_button , SENSITIVE=jj*is_imaging
if project.procPramArray[ii].rotate_Direction eq 1 then $
    widget_control, deg_90_button , /SET_BUTTON

widget_control, deg_180_button , SENSITIVE=jj*is_imaging
if project.procPramArray[ii].rotate_Direction eq 2 then $
    widget_control, deg_180_button , /SET_BUTTON

widget_control, deg_270_button , SENSITIVE=jj*is_imaging
if project.procPramArray[ii].rotate_Direction eq 3 then $
    widget_control, deg_270_button , /SET_BUTTON



if project.procPramArray[ii].slice_axis eq 0 then $
    widget_control, freq_phase_button , /SET_BUTTON


if project.procPramArray[ii].slice_axis eq 1 then $
    widget_control, freq_slice_button , /SET_BUTTON


if project.procPramArray[ii].slice_axis eq 2 then $
    widget_control, phase_slice_button , /SET_BUTTON


widget_control , freq_phase_button  ,SENSITIVE=jj*is_imaging
widget_control , freq_slice_button  ,SENSITIVE=jj*is_imaging
widget_control , phase_slice_button ,SENSITIVE=jj*is_imaging


IF project.procPramArray[ii].mc_enable EQ 1 THEN BEGIN
    widget_control, W_MENU_SUB_MOTION_CORRECT_MI_ON, SENSITIVE=1
    widget_control, W_MENU_SUB_MOTION_CORRECT_ED_ON, SENSITIVE=0
    widget_control, W_MENU_SUB_MOTION_CORRECT_OFF,   SENSITIVE=1
ENDIF ELSE IF project.procPramArray[ii].mc_enable eq 2 then begin
    widget_control, W_MENU_SUB_MOTION_CORRECT_MI_ON, SENSITIVE=0
    widget_control, W_MENU_SUB_MOTION_CORRECT_ED_ON, SENSITIVE=1
    widget_control, W_MENU_SUB_MOTION_CORRECT_OFF,   SENSITIVE=1
ENDIF ELSE BEGIN
    widget_control, W_MENU_SUB_MOTION_CORRECT_MI_ON, SENSITIVE=1
    widget_control, W_MENU_SUB_MOTION_CORRECT_ED_ON, SENSITIVE=1
    widget_control, W_MENU_SUB_MOTION_CORRECT_OFF,   SENSITIVE=0
ENDELSE

if (project.imndArray[ii].image_type eq 5 or $
    project.imndarray[ii].image_type eq 16 or $
    project.imndarray[ii].image_type eq 18) then begin

    widget_control, W_MENU_MAGNITUDE_TYPE, SENSITIVE=0
    if project.procPramArray[ii].signal_type eq 0 then $
        widget_control, W_MENU_MAGNITUDE_TYPE, SENSITIVE=0

    widget_control, W_MENU_REAL_TYPE , SENSITIVE=0
    if project.procPramArray[ii].signal_type eq 1 then $
        widget_control, W_MENU_REAL_TYPE , SENSITIVE=0

    widget_control, W_MENU_IMAGNARY_TYPE , SENSITIVE=0
    if project.procPramArray[ii].signal_type eq 2 then $
        widget_control, W_MENU_IMAGNARY_TYPE , SENSITIVE=0

    widget_control, W_MENU_INTENSITY_TYPE , SENSITIVE=0
    if project.procPramArray[ii].signal_type eq 3 then $
        widget_control, W_MENU_INTENSITY_TYPE , SENSITIVE=0

    widget_control, W_MENU_PHASE , SENSITIVE=0
    if project.procPramArray[ii].signal_type eq 4 then $
        widget_control, W_MENU_PHASE , SENSITIVE=0

	widget_control, W_MENU_RAWM , SENSITIVE=0
    if project.procPramArray[ii].signal_type eq 5 then $
        widget_control, W_MENU_RAWM , SENSITIVE=0

	widget_control, W_MENU_RAWR , SENSITIVE=0
    if project.procPramArray[ii].signal_type eq 6 then $
        widget_control, W_MENU_RAWR , SENSITIVE=0

	widget_control, W_MENU_RAWI , SENSITIVE=0
    if project.procPramArray[ii].signal_type eq 7 then $
        widget_control, W_MENU_RAWI , SENSITIVE=0
  
  widget_control, W_MENU_COMPLEX_TYPE, SENSITIVE=0
    if project.procPramArray[ii].signal_type eq 9 then $
        widget_control, W_MENU_COMPLEX_TYPE, SENSITIVE=0


end else begin

    widget_control, W_MENU_MAGNITUDE_TYPE, SENSITIVE=jj
    if project.procPramArray[ii].signal_type eq 0 then $
        widget_control, W_MENU_MAGNITUDE_TYPE, SENSITIVE=0

    widget_control, W_MENU_REAL_TYPE , SENSITIVE=jj
    if project.procPramArray[ii].signal_type eq 1 then $
        widget_control, W_MENU_REAL_TYPE , SENSITIVE=0

    widget_control, W_MENU_IMAGNARY_TYPE , SENSITIVE=jj
    if project.procPramArray[ii].signal_type eq 2 then $
        widget_control, W_MENU_IMAGNARY_TYPE , SENSITIVE=0

    widget_control, W_MENU_INTENSITY_TYPE , SENSITIVE=jj
    if project.procPramArray[ii].signal_type eq 3 then $
        widget_control, W_MENU_INTENSITY_TYPE , SENSITIVE=0

    widget_control, W_MENU_PHASE , SENSITIVE=jj
    if project.procPramArray[ii].signal_type eq 4 then $
        widget_control, W_MENU_PHASE , SENSITIVE=0

	widget_control, W_MENU_RAWM , SENSITIVE=jj
    if project.procPramArray[ii].signal_type eq 5 then $
        widget_control, W_MENU_RAWM , SENSITIVE=0

	widget_control, W_MENU_RAWR , SENSITIVE=jj
    if project.procPramArray[ii].signal_type eq 6 then $
        widget_control, W_MENU_RAWR , SENSITIVE=0

	widget_control, W_MENU_RAWI , SENSITIVE=jj
    if project.procPramArray[ii].signal_type eq 7 then $
        widget_control, W_MENU_RAWI , SENSITIVE=0
  
   widget_control, W_MENU_COMPLEX_TYPE , SENSITIVE=jj
    if project.procPramArray[ii].signal_type eq 9 then $
        widget_control, W_MENU_COMPLEX_TYPE , SENSITIVE=0

end

widget_control , INTENSITY_SLIDER_MAX , SET_VALUE = project.procPramArray[ii].intensity_max
widget_control , INTENSITY_SLIDER     , SET_VALUE = project.procPramArray[ii].intensity_cen
widget_control , INTENSITY_SLIDER_MIN , SET_VALUE = project.procPramArray[ii].intensity_min

widget_control , INTENSITY_SLIDER_MAX ,SENSITIVE=1*is_imaging
widget_control , INTENSITY_SLIDER     ,SENSITIVE=1*is_imaging
widget_control , INTENSITY_SLIDER_MIN ,SENSITIVE=1*is_imaging

if (project.procpramarray[ii].state_1 eq 1) then begin
    widget_control , DISPLAY_SCALE_RESET_BUTTON, sensitive=1*is_imaging
    widget_control , DISPLAY_DATA_MAX, sensitive=1*is_imaging, $
                     SET_VALUE=strtrim(string(project.procPramArray[ii].max_display_thresh),2)
    widget_control , DISPLAY_DATA_MIN, sensitive=1*is_imaging, $
                     SET_VALUE=strtrim(string(project.procPramArray[ii].min_display_thresh),2)
endif else begin
    widget_control , DISPLAY_SCALE_RESET_BUTTON, sensitive=0
    widget_control , DISPLAY_DATA_MAX, sensitive=0, $
                     SET_VALUE=strtrim(string(project.procPramArray[ii].max_display_thresh),2)
    widget_control , DISPLAY_DATA_MIN, sensitive=0, $
                     SET_VALUE=strtrim(string(project.procPramArray[ii].min_display_thresh),2)
endelse


temp = project.imndArray[ii].adim
if (temp eq 1) or (temp eq 0) then begin ; have to set the slider to not-sensitive
    widget_control, ADIM_SLIDER ,SENSITIVE=0 , SET_SLIDER_MAX=temp
    widget_control, W_MENU_multi_slice_movie, SENSITIVE=0
endif  else begin
    widget_control, ADIM_SLIDER ,SENSITIVE=1 , SET_SLIDER_MAX=temp
    widget_control, W_MENU_multi_slice_movie, SENSITIVE=1
end

widget_control , ADIM_SLIDER , SET_VALUE = project.procPramArray[ii].adim_start $
                   , SET_SLIDER_MAX=project.imndArray[ii].adim -1

sdim_slider_callback

if project.procPramArray[ii].single_Multi_flag eq 1 then begin
    widget_control, SLIDER_DISPLAY_CONTROL , SET_VALUE='Multiple' ,SENSITIVE=1*is_imaging
endif else begin
    widget_control, SLIDER_DISPLAY_CONTROL , SET_VALUE='Single' ,SENSITIVE=1*is_imaging
endelse


if project.procPramArray[ii].sort_dir eq 0 then BEGIN
    CASE project.procPramArray[ii].slice_axis OF
       0:widget_control, SLIDER_DIMENSION_CONTROL , SET_VALUE='Slice dir',SENSITIVE=1*is_imaging
       1:widget_control, SLIDER_DIMENSION_CONTROL , SET_VALUE='Phase dir',SENSITIVE=1*is_imaging
       2:widget_control, SLIDER_DIMENSION_CONTROL , SET_VALUE='Freq  dir',SENSITIVE=1*is_imaging
    ENDCASE

END else widget_control, SLIDER_DIMENSION_CONTROL , SET_VALUE='Array dir',SENSITIVE=1*is_imaging



if project.imndArray[ii].image_type eq 5 then sensitivity = 0 $
else sensitivity = 1

WIDGET_CONTROL,TIME_ZERO_SLIDER, SENSITIVE=jj*is_imaging


widget_control    ,p_dim_shift,SENSITIVE=jj*is_imaging   $
            ,SET_SLIDER_Min=-project.imndArray[ii].pdim $
            ,SET_SLIDER_MAX=project.imndArray[ii].pdim $
            ,SET_VALUE=project.imndArray[ii].pdim_shift

widget_control    ,s_dim_shift ,SENSITIVE=jj*is_imaging $
         ,SET_SLIDER_Min=-project.imndArray[ii].sdim $
         ,SET_SLIDER_MAX=project.imndArray[ii].sdim $
                 ,SET_VALUE=project.imndArray[ii].sdim_shift

widget_control    ,f_dim_shift ,SENSITIVE=jj*is_imaging $
         ,SET_SLIDER_Min=-project.imndArray[ii].fdim $
         ,SET_SLIDER_MAX=project.imndArray[ii].fdim $
                 ,SET_VALUE=project.imndArray[ii].fdim_shift

widget_control  ,TIME_ZERO_SLIDER $
         ,SET_SLIDER_Min = project.imndArray[ii].time.min $
         ,SET_VALUE = project.imndArray[ii].time.zero $
         ,SET_SLIDER_MAX = project.imndArray[ii].time.max

widget_control    ,kspace_pdim_shift,SENSITIVE=jj*is_imaging   $
            ,SET_SLIDER_Min=-project.imndArray[ii].pdim $
            ,SET_SLIDER_MAX=project.imndArray[ii].pdim $
            ,SET_VALUE=project.imndArray[ii].k_pdim_shift

widget_control    ,kspace_sdim_shift ,SENSITIVE=jj*is_imaging $
         ,SET_SLIDER_Min=-project.imndArray[ii].sdim $
         ,SET_SLIDER_MAX=project.imndArray[ii].sdim $
                 ,SET_VALUE=project.imndArray[ii].k_sdim_shift

widget_control    ,kspace_fdim_shift ,SENSITIVE=jj*is_imaging $
         ,SET_SLIDER_Min=-project.imndArray[ii].fdim $
         ,SET_SLIDER_MAX=project.imndArray[ii].fdim $
                 ,SET_VALUE=project.imndArray[ii].k_fdim_shift

widget_control    ,kspace_fdim_subsamp ,SENSITIVE=jj*is_imaging $
         ,SET_SLIDER_Min=min([1, project.imndarray[ii].k_fdim_span_max]) $
         ,SET_SLIDER_MAX=project.imndarray[ii].k_fdim_span_max $
         ,SET_VALUE=project.imndArray[ii].k_fdim_span

widget_control    ,kspace_pdim_subsamp ,SENSITIVE=jj*is_imaging $
         ,SET_SLIDER_Min=min([1, project.imndarray[ii].k_pdim_span_max]) $
         ,SET_SLIDER_MAX=project.imndarray[ii].k_pdim_span_max $
         ,SET_VALUE=project.imndArray[ii].k_pdim_span

if (project.imndarray[ii].k_sdim_span_max eq 1) then enable_slice = 0 else enable_slice = 1
widget_control    ,kspace_sdim_subsamp ,SENSITIVE=enable_slice*jj*is_imaging $
         ,SET_SLIDER_Min=min([1, project.imndarray[ii].k_sdim_span_max]) $
         ,SET_SLIDER_MAX=project.imndarray[ii].k_sdim_span_max $
         ,SET_VALUE=project.imndArray[ii].k_sdim_span

roi_mask_display = (*project.roi.pdisplay_names)
roi_mask_display[0,0] = 'None'

IF project.roi.mask GE (SIZE(roi_mask_display))[2]-1 THEN $
      project.roi.mask = (SIZE(roi_mask_display))[2]-1

widget_control, w_roi_mask, set_value=roi_mask_display,  SET_DROPLIST_SELECT= project.roi.mask
widget_control, w_roi_xform, set_button=project.procpramarray[project.ci].no_transform_roi
widget_control, roi_graph_toggle, set_button=project.procpramarray[project.ci].curvefit_graph

; ====== This part changes the sensitivity of the Fiber Track Mapping button  depending on the available
; processed data from MAS ===========
; It is commented because user might want to load a track file directly from the fibers modal window.
; HS-20061204 (at 35,000 feet somewhere over craphole Texas)
; BTW I am keeping this code commented out here as a future reference on how to change the sensitivity of a widget.

IF ptr_valid(project.dataArray[project.CI].adt) AND  $
  ptr_valid(project.dataArray[project.CI].eign_val) AND $
  ptr_valid(project.dataArray[project.CI].eign_Vec) AND  $
  ptr_valid(project.dataArray[project.CI].frac_Ani) AND $
  ptr_valid(project.dataArray[project.CI].Avg_Dif) THEN BEGIN
    widget_control, W_MENU_FIBERS, sensitive = 1*is_imaging
    widget_control, W_MENU_FIBERS_HARDI, sensitive=1*is_imaging ;;PROB, sensitive = 1
    WIDGET_CONTROL, W_MENU_DIFFTOOLS, sensitive=1
ENDIF else begin
    widget_control, W_MENU_FIBERS, sensitive = 0
    widget_control, W_MENU_FIBERS_HARDI, sensitive=1*is_imaging;; PROB, sensitive = 0
    WIDGET_CONTROL, W_MENU_DIFFTOOLS, sensitive=0
endelse

; Make the IMND option sensitive
WIDGET_CONTROL,W_MENU_display_IMND, SENSITIVE=jj
WIDGET_CONTROL,W_MENU_display_PROCPAR, SENSITIVE=jj
WIDGET_CONTROL,W_MENU_display_PROCPAR_pretty, SENSITIVE=jj
WIDGET_CONTROL,W_MENU_display_ACQP, SENSITIVE=jj
WIDGET_CONTROL,W_MENU_display_CONFIG, SENSITIVE=jj
WIDGET_CONTROL,W_MENU_display_PAR, SENSITIVE=jj

mdt_redraw

redraw_ADT
mas_background_filter_redraw
mas_cfit2_gui_redraw, project.procpramarray[project.ci].curvefit2_state

end


; Subroutine name: mas_redraw
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/02
    ;Fix spelling mistakes and commenting.

PRO mas_redraw
COMPILE_OPT IDL2

END
