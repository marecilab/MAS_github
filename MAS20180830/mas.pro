;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved

pro WID_BASE_MAIN_realized, event
    print, 'WID_BASE_MAIN_realized: called.'
end


;****************************************************************************************************
;
; NAME:
;   WID_BASE_MAIN_event
;
; PURPOSE: Main eveny handler for MAS's main window.
;
;
; ARGUMENTS: Event - event structure
;
;
; MODIFICATION HISTORY:
;
;
;***************************************************************************************************
pro WID_BASE_MAIN_event, Event
  HEAP_GC
  COMPILE_OPT IDL2

;   ;error handling.
;   catch,error_status
;   if (error_status ne 0) then begin
;     help, calls= trace_back
;
;     dummy = DIALOG_MESSAGE(['Please send a copy of this error message to the developer of this program.', $
;      !ERROR_STATE.msg,trace_back], /ERROR, $
;         TITLE='Error in WID_BASE_MAIN_event')
;
;
;     RETURN
;   endif


  wWidget =  Event.top
 ; help, event, /structure

  case Event.id of

    Widget_Info(wWidget, FIND_BY_UNAME='do_nothing'): begin

    end

    ; =================[ pull down menu widgets ]=========================

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_MAIN_EXIT'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then  $
        mas_exit, Event
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_IMND'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        mas_open, 0
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_ONEPULSE'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        mas_open_onepulse_spect
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_SUMMARY_AGILENT'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        mas_varian_summarize_experiment_dir
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_AGILENT_SINGLE'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        mas_varian_open_fid_dir
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_AGILENT_FDF'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        mas_varian_open_fdf_dir
    end
       Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_PROJECT'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        mas_open_project, Event
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_DYNAMICS'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         open_dynamics_multi
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_T2MULTI'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         mas_open_T2_Multi
   end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_ADCMULTI'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         mas_open_ADC_Multi
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_DICOM_GENERIC'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         mas_dicom_opendir
    end
    ;; replacede by generic DICOM reader
    ;;Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_DICOM_PHILIPS'): begin
    ;;  if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
    ;;     mas_open_multi_dicom
    ;;end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_PHILIPS_DICOM_SPECT'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         read_philips_dicom_spect_dir
    end
    
	Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_PARREC'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         readPARdata
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_SPARSDAT'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         mas_read_philips_sparsdat_spect
    end
    
	Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_NIFTI_DATA'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         mas_open_nifti
    end
    
  Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_NIFTI_MAPS'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         trial_open_qmri_maps_dir
    end
    
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_ADT'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         mas_open_adt
    end
    ;; Adding menu option to open ADT single acquisitions. This just
    ;; maps to mas_open, since mas_open determines the correct opener
    ;; based on contents.
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_ADT_SINGLE'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         mas_open
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_phased_array'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         mas_open_phase_array
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_bruker_image'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         mas_open_bruker_image
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_bruker_multi_image'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         mas_open_bruker_multi_image
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_multi_slice_movie'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         open_multi_slice_movie
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_FLT'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         open_flt_file
    end
	Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPEN_RAW'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         open_raw_file
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_SAVE_PROJECT'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        save_project, Event
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_SAVE_TIFF'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        save_tiff
  end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_SAVE_MPEG'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        MAS_MPEG_WIDGET;save_mpeg
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_SAVE_FLT'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        save_flt
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_SAVE_SUBSET'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        save_subset
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_SAVE_SLICER'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        write_slicer_image
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_SAVE_NIFTI'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         mas_export_nifti_gui
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_SAVE_DICOM'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         mas_export_dicom_gui
    end

    Widget_Info(wWidget, FIND_BY_UNAME='menu_save_raw_data'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        save_raw
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_REMOVE_SCAN'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mas_remove_scan
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_REMOVE_SCAN_ALL'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mas_remove_scan_all
    end

;    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_OPTIONS'): begin
;      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
;        mas_options ,event
;    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_COLOR_TABLE'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        xloadct
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_SLICER'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mas_slicer_3d
  end

  Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_orthoview'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mas_display_ortho
    end
    
   Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_multi_slice_movie'): begin
    if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
      mas_multi_movie_gui
  end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_volume'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mas_volume_slicer; mas_ivolume
    end
    
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_surface'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mas_surface_render
    end
    
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_MIP'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mip_gui
    end
    
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_CURVEFIT_2'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mas_cfit2_gui
    end
    
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_QMRI_PD_BTN'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        trial_make_proton_map 
    end
    
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_QMRI_TOOL_BTN'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        trial_qmri_gui 
    end
    
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_MREIT'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mas_mreit
    end
    
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_IMAGE_STATS'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        do_image_statistics, event ;mas_image_statistics
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_IMAGE_HISTOGRAM'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        do_image_histogram, event
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_DCE_DATA_SELECT'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        select_data_dce
    end
    
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_SIGNAL_ENHANCEMENT'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mas_signal_enhancement
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_DCE_ROI_CURVES'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        DCE_ROI_CURVE_FIT_GUI   ;mas_dce_roi_init
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_ALG'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mas_algebra
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_ADT_GRADIENT'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        adt_show_gradient_profile
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_ADT_REGRESS'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mas_adt_regress
  end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_DIFFTOOLS'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
       mas_diffusion_tools_gui
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_FIBERS'): begin
       if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        FibersGUI ;fibers
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_FIBERS_HARDI'): begin
        if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
          mas_hardi_tractography_gui ;;mas_probtrack_gui
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_SMOOTH'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mas_smooth
    end
    
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_FLOW_GRADIENT'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
         flow_show_gradient_profile 
    end
    
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_FLOW_EXTRACT'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $ 
       flow_gui
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_GPE'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        Bez_gui
    end
    
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_VOLUMETRICS'): begin
      mas_measure_volume
   end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_Trace_Phase'): begin
      mas_Trace, /phase
    end
    
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_calculate_delta_b_algebra'): begin
      mas_calc_delta_b_algebra_gui
    end
    
     Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_calculate_delta_b_least_squares'): begin
      mas_calc_delta_b_least_squares_gui
    end

     Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_CALC_PROTON_DENSITY'): begin
      mas_curvefit_pd_init
    end
    
       Widget_Info(wWidget, FIND_BY_UNAME='deg_0_button'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_rotate_Direction,0
    end

       Widget_Info(wWidget, FIND_BY_UNAME='deg_90_button'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_rotate_Direction,1
    end

       Widget_Info(wWidget, FIND_BY_UNAME='deg_180_button'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_rotate_Direction,2
    end

       Widget_Info(wWidget, FIND_BY_UNAME='deg_270_button'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_rotate_Direction,3
    end

       Widget_Info(wWidget, FIND_BY_UNAME='flip_none_button'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_flip_Direction, 0
    end
       Widget_Info(wWidget, FIND_BY_UNAME='flip_horz_button'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_flip_Direction, 1
    end
       Widget_Info(wWidget, FIND_BY_UNAME='flip_vert_button'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_flip_Direction, 2
    end


; Ami's smoothing buttons

 	Widget_Info(wWidget, FIND_BY_UNAME='smooth_none_button'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_smooth_Direction, 0
    end
       Widget_Info(wWidget, FIND_BY_UNAME='smooth_button'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_smooth_Direction, 1
    end
       Widget_Info(wWidget, FIND_BY_UNAME='median_button'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_smooth_Direction, 2
    end

    Widget_Info(wWidget, FIND_BY_UNAME='freq_phase_button'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_slice_axis, 0
    end
    Widget_Info(wWidget, FIND_BY_UNAME='freq_slice_button'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_slice_axis, 1
    end
    Widget_Info(wWidget, FIND_BY_UNAME='phase_slice_button'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_slice_axis, 2
    end

       Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_MAGNITUDE_TYPE'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_signal_type, 0
    end

       Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_REAL_TYPE'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_signal_type, 1
    end

       Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_IMAGNARY_TYPE'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_signal_type, 2
    end

       Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_INTENSITY_TYPE'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_signal_type, 3
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_PHASE'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_signal_type, 4
    end

	Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_RAWM'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_signal_type, 5
    end

	Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_RAWR'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_signal_type, 6
    end

	Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_RAWI'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        set_signal_type, 7
    end
    
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_COMPLEX_TYPE'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
      set_signal_type, 9
      end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_IMAGE_PHASE'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then begin
            mas_image_phase
      endif
    end
      
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_SUB_MOTION_CORRECT_OFF'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        toggle_motion_correction, TYPE=0
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_SUB_MOTION_CORRECT_MI_ON'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        toggle_motion_correction, TYPE=2
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_SUB_MOTION_CORRECT_ED_ON'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        toggle_motion_correction, TYPE=1
    end

    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_SUB_MOTION_CORRECT_PCA'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mas_mc_pca_gui
    end

    ;===================[ scan list selection ]============================

    Widget_Info(wWidget, FIND_BY_UNAME='WID_SCAN_LIST'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_LIST' ) then $
        update_scan_selection_event, Event
    end

    ;=====================[ zoom control ]===========================
    Widget_Info(wWidget, FIND_BY_UNAME='freq_interp_slider'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_SLIDER' ) then $
        mas_zoom_tab_callback, event
    end
    Widget_Info(wWidget, FIND_BY_UNAME='phase_interp_slider'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_SLIDER' ) then $
        mas_zoom_tab_callback, event
    end
    Widget_Info(wWidget, FIND_BY_UNAME='slice_interp_slider'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_SLIDER' ) then $
        mas_zoom_tab_callback, event
    end
    Widget_Info(wWidget, FIND_BY_UNAME='lock_interp_chk_box'): $
        mas_zoom_tab_callback, event

	Widget_Info(wWidget, FIND_BY_UNAME='down_sample_chk_box'): $
        mas_zoom_tab_callback, event
  
  Widget_Info(wWidget, FIND_BY_UNAME='zpad_chk_box'): $
        mas_zoom_tab_callback, event        

    Widget_Info(wWidget, FIND_BY_UNAME='x_zoom_slider'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_SLIDER' ) then $
        mas_zoom_tab_callback, event
    end
    Widget_Info(wWidget, FIND_BY_UNAME='y_zoom_slider'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_SLIDER' ) then $
        mas_zoom_tab_callback, event
    end
    
    ;=====================[ intensity control ]=========================
    Widget_Info(wWidget, FIND_BY_UNAME='INTENSITY_SLIDER'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_SLIDER' ) then $
        intensity_slider_callback
    end
    Widget_Info(wWidget, FIND_BY_UNAME='INTENSITY_SLIDER_MIN'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_SLIDER' ) then $
        intensity_slider_callback
    end
    Widget_Info(wWidget, FIND_BY_UNAME='INTENSITY_SLIDER_MAX'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_SLIDER' ) then $
        intensity_slider_callback
    end

    Widget_Info(wWidget, FIND_BY_UNAME='DISPLAY_SCALE_RESET_BUTTON'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        display_threshold_callback, event
    end

    Widget_Info(wWidget, FIND_BY_UNAME='DISPLAY_DATA_MAX'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_TEXT_CH' ) then $
        display_threshold_callback, event
    end

    Widget_Info(wWidget, FIND_BY_UNAME='DISPLAY_DATA_MIN'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_TEXT_CH' ) then $
        display_threshold_callback, event
    end


    ;=====================[ setting control ]============================
    Widget_Info(wWidget, FIND_BY_UNAME='f_dim_shift'): $
        mas_settings_event, event

    Widget_Info(wWidget, FIND_BY_UNAME='p_dim_shift'): $
        mas_settings_event, event

    Widget_Info(wWidget, FIND_BY_UNAME='s_dim_shift'): $
        mas_settings_event, event

    Widget_Info(wWidget, FIND_BY_UNAME='TIME_ZERO_SLIDER'): $
        mas_settings_event, event

    Widget_Info(wWidget, FIND_BY_UNAME='SLIDER_PREIMAGES'): $
        mas_settings_event, event

    Widget_Info(wWidget, FIND_BY_UNAME='n_avg_flag_button'): $
        mas_settings_event, event
    
    Widget_Info(wWidget, FIND_BY_UNAME='phi_unwrap_flag_button'): $
        mas_settings_event, event
    
    Widget_Info(wWidget, FIND_BY_UNAME='kspace_fdim_shift'): $
        mas_settings_event, event

    Widget_Info(wWidget, FIND_BY_UNAME='kspace_pdim_shift'): $
        mas_settings_event, event    
    
    Widget_Info(wWidget, FIND_BY_UNAME='kspace_sdim_shift'): $
        mas_settings_event, event  

    Widget_Info(wWidget, FIND_BY_UNAME='kspace_fdim_span'): $
        mas_settings_event, event

    Widget_Info(wWidget, FIND_BY_UNAME='kspace_pdim_span'): $
        mas_settings_event, event    
    
    Widget_Info(wWidget, FIND_BY_UNAME='kspace_sdim_span'): $
        mas_settings_event, event  


    ;=====================[ ROI control ]======================


    Widget_Info(wWidget, FIND_BY_UNAME='w_roi_list'): $
        mas_roi_event, event


    Widget_Info(wWidget, FIND_BY_UNAME='w_roi_name'): $
        mas_roi_event, event


    Widget_Info(wWidget, FIND_BY_UNAME='w_roi_load'): $
        mas_roi_event, event

    Widget_Info(wWidget, FIND_BY_UNAME='w_roi_new_edit'): $
        mas_roi_event, event

    Widget_Info(wWidget, FIND_BY_UNAME='w_roi_save'): $
        mas_roi_event, event

    Widget_Info(wWidget, FIND_BY_UNAME='w_roi_delete'): $
        mas_roi_event, event

    Widget_Info(wWidget, FIND_BY_UNAME='w_roi_mask'): $
        mas_roi_event, event

    Widget_Info(wWidget, FIND_BY_UNAME='w_roi_xform'): $
        mas_roi_event, event

    Widget_Info(wWidget, FIND_BY_UNAME='roi_graph_toggle'): $
        mas_roi_event, event
    ;=====================[ display control ]======================


    Widget_Info(wWidget, FIND_BY_UNAME='SLIDER_DISPLAY_CONTROL'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )  then $
        single_Multi_flag_toggle
    end
    Widget_Info(wWidget, FIND_BY_UNAME='SLIDER_DIMENSION_CONTROL'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        sort_dir_toggle
    end
     Widget_Info(wWidget, FIND_BY_UNAME='DISPLAY_GO'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mas_display
    end
    Widget_Info(wWidget, FIND_BY_UNAME='DISPLAY_MOVIE'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mas_movie
    end
    Widget_Info(wWidget, FIND_BY_UNAME='iImage_flag_button'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mas_settings_event, event
    end
    Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_AUTO_DISPLAY'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        auto_display_toggle, event
    end
    Widget_Info(wWidget, FIND_BY_UNAME='SDIM_SLIDER'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_SLIDER' ) then $
        sdim_slider_callback
    end
    Widget_Info(wWidget, FIND_BY_UNAME='ADIM_SLIDER'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_SLIDER' ) then $
        adim_slider_callback
    end

; HS June 12 2007
; Display IMND routine

 Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_display_IMND'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        show_text_file, 'imnd'
      end

 Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_display_ACQP'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        show_text_file, 'acqp'
      end

 Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_display_CONFIG'): begin
        if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
          show_text_file, 'config'
      end
      
; CD August 20 2007
 Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_display_PAR'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        show_text_file, 'PAR'
      end

 Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_display_PROCPAR'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        show_text_file, 'procpar'
      end
      
 Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_display_PROCPAR_Pretty'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        show_text_file, 'procpar-pretty'
      end

; HS August 05, 2007
; Frequency domain filtering

 Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_FREQUENCY_FILTER_LOWPASS'): begin
	if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        WID_BASE_Frequency_Domain_lowpass, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_
    end

 Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_FREQUENCY_FILTER_HIGHPASS'): begin
	if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        WID_BASE_Frequency_Domain_highpass, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_
    end


; HS August 10, 2007
; Image domain filtering

Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_IMAGE_DOMAIN'): begin
	if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        WID_BASE_IMAGE_DOMAIN_FILTERING, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_
    end

; BT Jan 2009
; data orientation routines
Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_ORIENT_DATA'): begin
	if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        mas_orient_data
    end
;GA July 2009
;DCE Processing GUI
Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_DCE_SCC'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mas_dce_processing   ;mas_dce_StoC
    end
;GA Oct 2010
;DCE Volume Processing GUI
Widget_Info(wWidget, FIND_BY_UNAME='W_MENU_DCE_VOL'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' ) then $
        mas_dce_volumes_gui  
    end

  endcase


end


;   Makes a common block "global variables" of widget id's that can be called later on.
;   Creates the main widget on which every other widget hangs on to.
;   This procedure is divided up into location of widgets.


;*****************************************************************************************************
;
; NAME:
;   WID_BASE_MAIN
;
; PURPOSE:
;   Makes a common block of widget ID's that can be called later on.
;   Creates the main widget on which every other widget hangs on to.
;   This procedure was created by the idlde and then later edited.
;
; ARGUMENTS:
;
;
; MODIFICATION HISTORY:
;
;	HS - 20061207
;	No more mas_lite flag verification.
; Added B0 and removed slice dir (Magdoom, 05/28/2016)
;*****************************************************************************************************
pro WID_BASE_MAIN, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_ , mas_lite
  COMPILE_OPT IDL2

  ;note: FRAME_ON is a variable used to help place the label widgets on the screen b/c one can see
  ; the boundaries for all the labels
  FRAME_ON =0 ;

  common scan_data, project
  ;note: by making the widgets that i change alot common to all procedures this is very handy
  ;if you are going to change the VALUE of one
  
  ;; BT: 2013-01-28 - This common block contains the widget variables for pretty much every widget
  ;;                  related to MAS. As new features have been added, their GUI specific widgets
  ;;                  have been put into a compact struct which is attached to the top-level widget
  ;;                  base. This is an IDL programming practice which makes things a lot cleaner when
  ;;                  you have a lot of widgets.
  common common_widgets, $
    WID_BASE_MAIN, WID_SCAN_LIST, $
    LABEL_RARE_DATA, LABEL_MAX_DATA,LABEL_SF_DATA, $
    LABEL_FDIM_DATA, LABEL_PDIM_DATA, LABEL_SDIM_DATA, LABEL_SLICE_DIR_DATA, LABEL_SLICE_ORIENT_DATA, $
    LABEL_ADIM_DATA, LABEL_SCANNAME_DATA, LABEL_RECOVTIME_DATA ,LABEL_ECHOTIME_DATA, $
    LABEL_NECHO_DATA, LABEL_NAVG_DATA, LABEL_FFOV_DATA ,LABEL_PFOV_DATA, LABEL_SCANDATE_DATA, $
    LABEL_SFOV_DATA, LABEL_SLICE_DATA ,LABEL_THICK_DATA, LABEL_FATSUP_DATA ,LABEL_DIMENSIONS_DATA, $
    $
    W_MENU_COLOR_TABLE, W_MENU_SAVE_TIFF, W_MENU_SAVE_PROJECT, W_MENU_REMOVE_SCAN_MAIN, W_MENU_REMOVE_SCAN, $
    W_MENU_REMOVE_SCAN_ALL, W_MENU_AUTO_DISPLAY, W_MENU_SLICER, W_MENU_orthoview, W_MENU_multi_slice_movie, $
    W_MENU_OPTIONS, W_MENU_OPEN_DYNAMICS, W_MENU_SIGNAL_ENHANCEMENT, W_MENU_OPEN_DICOM_GENERIC, $
    W_MENU_OPEN_DICOM_PHILIPS, W_MENU_ADT_REGRESS, W_MENU_FREQ_PHASE, W_MENU_FREQ_SLICE, W_MENU_PHASE_SLICE, $
    W_MENU_DCE_ROI_CURVES, W_MENU_DCE_SCC, W_MENU_DCE_VOL, W_MENU_MAGNITUDE_TYPE, W_MENU_REAL_TYPE, W_MENU_IMAGNARY_TYPE, $
    W_MENU_INTENSITY_TYPE, W_MENU_COMPLEX_TYPE, W_MENU_PHASE, $
    W_MENU_CURVEFIT_2, $
    $;; Commented-out original Curve Fit common block members
    $;;W_MENU_T1SR_GE_ROI, W_MENU_T1SR_GE_IMAGE, $
    $;;W_MENU_T1SR_SE_ROI, W_MENU_T1SR_SE_IMAGE, W_MENU_T1IR_SE_ROI, W_MENU_T1IR_SE_IMAGE, $
    $;;W_MENU_T1IR_GE_ROI, W_MENU_T1IR_GE_IMAGE, $
    $;;W_MENU_T2MSE_ROI, W_MENU_T2MSE_IMAGE, W_MENU_ADC_ROI, W_MENU_ADC_IMAGE, W_MENU_ADC_STRETCH_ROI, W_MENU_ADC_STRETCH_IMAGE, $
    W_MENU_OPEN_ADT, W_MENU_SAVE_MPEG, W_MENU_OPEN_ADCMULTI, W_MENU_OPEN_T2MULTI, W_MENU_IMAGE_STATS, $
    W_MENU_CALC_PROTON_DENSITY, $
    $
    SLIDER_DISPLAY_CONTROL, SLIDER_DIMENSION_CONTROL, ZOOM_SLIDER, SDIM_SLIDER, ADIM_SLIDER, $
    INTENSITY_SLIDER_MIN, INTENSITY_SLIDER, INTENSITY_SLIDER_MAX, $
    DISPLAY_SCALE_RESET_BUTTON, DISPLAY_DATA_MAX, DISPLAY_DATA_MIN, $
    DISPLAY_GO, REFFT_IMAGE,DISPLAY_MOVIE,  OPTIONS_WINDOW_BASE, SLIDER_PREIMAGES, $
    OPTIONS_WINDOW_DONE, SE_WINDOW_BASE, SE_MINIMUM_REF, $
    SE_WINDOW_DISPLAY, SE_LABEL, LABEL_STATUS_BAR, SE_MAX_REBIN, SLIDER_ADTThreshold, $
    BVALUES_WINDOW_BASE, SLIDER_BVALUES, BVALUES_WINDOW_DONE, BVALUES_WINDOW_NEXT, BVALUES_WINDOW_BASE1, $
    BVALUES_TABLE, SDIM_SLIDER_TITLE, ZOOM_3D_TOGGLE, SDIM_SLIDER_MAX, $
    S0_min_slider, S0_max_slider, T_min_slider, T_max_slider, threshold_slider, image_fit_process, $
    image_fit_done, image_fit_base, image_fit_auto, image_fit_batch, $
    p_dim_shift, s_dim_shift, fit_zoomed_button,ADT_SMOOTH,ADT_FIT,ADT_DISPLAY,ADT_MULTI,ADT_SELECT_DATA, $
    ADT_FLT, ADT_TIF, ADT_Type0, ADT_Type4, ADT_Type2, ADT_Type3, ADT_type5, ADT_type6, Slider_ADT_Start, $
    Slider_ADT_Start_label, ADT_ROI, ADT_MULTI_display, wid_ADT_MOVIE, ADT_WINDOW_BASE, $
    bXX,bXY,bE1,bYY,bXZ,bE2,bZZ,bYZ,bE3,bS0,bFA,bAD,adt_zpad_stat, $
    ADT_scale_avg_d_max, ADT_scale_avg_d_min, ADT_scale_off_dia, ADT_scale_dia_max, $
    ADT_scale_dia_min, ADT_scale_fa_max, wid_ADT_MPG, TIME_ZERO_SLIDER, f_dim_shift, $
    kspace_fdim_shift, kspace_pdim_shift, kspace_sdim_shift, $
    kspace_fdim_subsamp, kspace_pdim_subsamp, kspace_sdim_subsamp, $
    alg_window_base , W_MENU_ALG,alg_load1,alg_load2,alg_add, alg_sub,   alg_mult, alg_div, alg_abs, $
    alg_inv, alg_text,W_MENU_volume, W_MENU_surface,W_MENU_MIP, iImage_flag_button, W_MENU_SAVE_FLT, $
    W_MENU_SMOOTH,W_MENU_SAVE_SUBSET, W_MENU_OPEN_phased_array, W_MENU_OPEN_bruker_image, $
    W_MENU_OPEN_bruker_multi_image, $
    freq_phase_button, freq_slice_button, phase_slice_button, $
    deg_0_button, deg_90_button, deg_180_button, deg_270_button, flip_none_button, flip_horz_button, $
    flip_vert_button, smooth_none_button, smooth_button, median_button, $
    freq_interp_slider, phase_interp_slider, slice_interp_slider, lock_interp_chk_box, down_sample_chk_box, $
    zpad_chk_box, x_zoom_slider, y_zoom_slider, menu_save_raw_data, W_MENU_OPEN_multi_slice_movie, $
    w_roi_list, w_roi_name, w_roi_load, w_roi_new_edit, w_roi_mask, w_roi_xform, roi_graph_toggle, $
    W_MENU_IMAGE_HISTOGRAM, ADT_ROI_WINDOW_BASE, $
    vol_window_base, vol_start, vol_stop, vol_draw,  W_MENU_VOLUMETRICS, n_avg_flag_button,phi_unwrap_flag_button, $
    ADT_MPEG_WINDOW_BASE, MAS_MPEG_WINDOW_BASE, REMOVE_BACKGROUND_BASE, W_MENU_FILTER, ADT_SMOOTH_WINDOW_BASE, $
    W_MENU_Trace_Phase, W_MENU_calculate_delta_b, W_MENU_calculate_delta_b_algebra, W_MENU_calculate_delta_b_least_squares, W_MENU_FIBERS, W_MENU_FIBERS_HARDI, ADT_WINDOW_BASE3, W_MENU_ADT_GRADIENT, W_MENU_DIFFTOOLS, $
    W_MENU_display_IMND, W_MENU_display_ACQP, W_MENU_display_CONFIG, W_MENU_Bruker_Files, W_MENU_display_procpar, W_MENU_display_procpar_pretty, W_MENU_SAVE_SLICER, wid_ADT_SLICER, $
    W_MENU_RAWM, W_MENU_RAWR, W_MENU_RAWI, W_MENU_IMAGE_PHASE, W_MENU_OPEN_PARREC, W_MENU_OPEN_NIFTI, W_MENU_SAVE_NIFTI, W_MENU_SAVE_DICOM, $
    W_MENU_FREQUENCY_FILTER, W_MENU_FREQUENCY_FILTER_LOWPASS,$
    W_MENU_FREQUENCY_FILTER_HIGHPASS, W_MENU_IMAGE_DOMAIN, W_MENU_Philips_Files, W_MENU_display_PAR, $
    W_MENU_MOTION_CORRECT, W_MENU_SUB_MOTION_CORRECT_ED_ON, W_MENU_SUB_MOTION_CORRECT_PCA, $
    W_MENU_SUB_MOTION_CORRECT_MI_ON, W_MENU_SUB_MOTION_CORRECT_OFF, $
    ADT_REG_STANDARD, ADT_REG_MOW,dce_window_base,  dce_t1_enable,dce_porosity_enable, dce_r1_text, dce_r2_text, dce_aif_base, aifroi_toggle,aifroi_caudate_droplist, aifroi_set_caudate_droplist, aifroi_set_droplist, aifroi_droplist, aifroi_set_base, aifroi_label, aifroi_set_label, kep_m_base, kep_a_base,aifroi_calc, aifroi_calc_data $
    ,dce_t_one_text, dce_t_two_text, dce_phi_text,dce_kie_text,dce_t2_enable,dce_sel_slice_slider, dce_sel_adim_slider,smooth_ktrans_enable, plroi_set_droplist, plroi_droplist,plot_droplist, aifroi_caudate_slice_droplist, aifroi_slice_droplist, plroi_display,sl_display_button, $
    t1_fit_lm_button, t1_fit_mp_button, t1_display_button, t1_select_button, conc_nm_button, conc_am_button,pl_display_button,conc_button_grp,vessroi,caudateroi, $
    conc_sm_button, cc_display_button, cc_export_button, kep_fit_options_base, kep_fit_lm_button, conc_2sx_button, $
    kep_fit_mp_button, kep_droplist_dce, kep_droplist_mask, mask_droplist, vtr_droplist,kep_export_button, kep_display_button, kep_droplist_vtr,  $
    nroi_set_droplist, nroi_droplist, kep_m1_text, kep_m2_text, kep_a1_text, kep_a2_text,ktransroi_set_droplist,ktransroi_droplist,kep_button_grp,$
    kep_dose_text, max_t1_text, fit_thresh_text,str_crit_slider, num_pst_slider, torigin_slider, conc_curves_f_text $
    ,conc_curves_p_text, conc_curves_s_text,cc_curves_button,filter_vessels_enable, kep_normal_toggle, nroi_set_base, nroi_toggle $
    , msm_slice_1_droplist, msm_slice_2_droplist, msm_slice_3_droplist, msm_slice_4_droplist, msm_slice_5_droplist $
    , msm_slice_6_droplist, msm_play_button, msm_export_button, cdb_ge_button, cdb_fp_button, cdb_loadx_button, cdb_loady_button, cdb_cd_button, cdb_cd3d_button, cdb_cd_display_droplist, cdb_cdb_gr_text $
    , cdb_display_button, cdb_loady_text, cdb_loadx_text,cdb_cdb_TE_label, cdb_cdb_TE_text, cdb_cd_display_button $
    , cdb_cd3D_display_button, cdb_3D_display_button, cdb_filter_toggle, cdb_threshold_text, cdb_unwrap_toggle, cdb_hipass_toggle $
    , dce_volumes_threshold_text, dce_volumes_ROI_slider, dce_volumes_ROI_edit_button, dce_volumes_display_button, dce_volumes_base $
    $ ;flow widget related variables 
    , vel_button, vroi_button_stat, vroi_button, vel_button_grp,vis_button_grp,movie_button $
    , disp_button,gfield,tfield1,tfield2, shpfield,schfield, fcfield,vencfield, ncyclfield, parbase, dispbase, pc_button, ioverlay_button, ecc_button $
    , Ithres,Vthres, ROIsetlist1_ecc, ROIlist1_ecc, ROIsetlist2_ecc, ROIlist2_ecc, Fitfn_button_grp, ecvel_button_grp, ROI_button_grp $
    , Bezbase, parbase_Bez, gaxis, gfield_Bez, tfield_Bez, schfield_Bez, Axlist_Bez, shpfield_Bez, Bencfield, pc_button_Bez,ap_button_Bez, Bez_constraint, $
     Ithres_Bez,sm_button_Bez, Tlist_Bez, BW_Bez, vis_button_grp_Bez,model_button_grp_Bez, movie_button_Bez,visbase_Bez, $
     Xlist, Ylist, Zlist,zmax_MIP,zmin_MIP, zf_button_Bez, Nfit_Bez
    
  mainXSixz = 450
  mainYSize = 700


  title = 'MAS Release #'+project.mas_version


    ;this is the main widget where every thing is atached to
    WID_BASE_MAIN = Widget_Base( GROUP_LEADER=wGroup, $
     UNAME='WID_BASE_MAIN' ,SCR_XSIZE=mainXSixz  $
        ,/ALIGN_CENTER ,/BASE_ALIGN_CENTER ,TITLE=title $
        ,SPACE=3 ,XPAD=3 ,YPAD=3 ,COLUMN=1  $
        ,MBAR=WID_BASE_MAIN_MBAR)


    fileXSize = mainXSixz-14
    fileYSize = 120

    ;this is the list widget where all the names to the scan are placed in a
    ;  string array.
      WID_SCAN_LIST = Widget_List(WID_BASE_MAIN, UNAME='WID_SCAN_LIST'  $
       ,FRAME=1 , SENSITIVE=0 $
        ,SCR_XSIZE=fileXSize ,SCR_YSIZE=fileYSize ,XSIZE=11 ,YSIZE=0)


;============================[holder info]=============================
 hXsize = fileXSize
 hYsize = 190

 BASE_HEADER = Widget_Base(WID_BASE_MAIN, UNAME='BASE_HEADER' $
      ,FRAME=0 $
      ,SCR_XSIZE=hXsize );,SCR_YSIZE=hYsize ,SPACE=0)

    ;this is the size of each column used to hold the label widget
 col1Xsize = hXsize/4-36
 col2Xsize = hXsize/4+36
 col3Xsize = hXsize/4-15
 col4Xsize = hXsize/4+15
 colYsize = hYsize-2

    ;this is the 4 columns that hold all the scan labels
 BASE_COL1 = Widget_Base(BASE_HEADER, UNAME='BASE_COL1' $
      ,FRAME=1 ,XOFFSET=0 ,YOFFSET=0 $
      ,SCR_XSIZE=col1Xsize ,SPACE=3 $;,SCR_YSIZE=colYsize $
      ,XPAD=0 ,YPAD=3 ,COLUMN=1 ,/ALIGN_RIGHT)

 BASE_COL2 = Widget_Base(BASE_HEADER, UNAME='BASE_COL2' $
      ,FRAME=1 ,XOFFSET=col1Xsize ,YOFFSET=0 $
      ,SCR_XSIZE=col2Xsize ,SPACE=3 $;,SCR_YSIZE=colYsize $
      ,XPAD=0 ,YPAD=3 ,COLUMN=1 ,/ALIGN_RIGHT)

 BASE_COL3 = Widget_Base(BASE_HEADER, UNAME='BASE_COL3' $
      ,FRAME=1 ,XOFFSET= col1Xsize+col2Xsize,YOFFSET=0 $
      ,SCR_XSIZE=col3Xsize  ,SPACE=3 $;,SCR_YSIZE=colYsize $
      ,XPAD=0 ,YPAD=3 ,COLUMN=1 ,/ALIGN_RIGHT)

 BASE_COL4 = Widget_Base(BASE_HEADER, UNAME='BASE_COL4' $
      ,FRAME=1 ,XOFFSET=col1Xsize+col2Xsize+col3Xsize ,YOFFSET=0 $
      ,SCR_XSIZE=col4Xsize  ,SPACE=3 $;,SCR_YSIZE=colYsize $
      ,XPAD=0 ,YPAD=3 ,COLUMN=1 ,/ALIGN_RIGHT)


;============================[header info]=============================

labYsize = 13

 ;SCANNAME
 LABEL_SCANNAME = Widget_Label(BASE_COL1, UNAME='LABEL_SCANNAME'  $
    ,SCR_XSIZE=col1Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='METHOD:')

 LABEL_SCANNAME_DATA = Widget_Label(BASE_COL2, UNAME='LABEL_SCANNAME_DATA'  $
    ,SCR_XSIZE=col2Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)

 ;SCANDATE
 LABEL_SCANDATE = Widget_Label(BASE_COL1, UNAME='LABEL_SCANDATE'  $
    ,SCR_XSIZE=col1Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='SCANDATE:')

 LABEL_SCANDATE_DATA = Widget_Label(BASE_COL2, UNAME='LABEL_SCANDATE_DATA'  $
    ,SCR_XSIZE=col2Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)

  ;FIELD STRENGTH
  LABEL_SF = Widget_Label(BASE_COL1, UNAME='LABEL_SF'  $
    ,SCR_XSIZE=col1Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='SF(MHz):')

  LABEL_SF_DATA = Widget_Label(BASE_COL2, UNAME='LABEL_SF_DATA'  $
    ,SCR_XSIZE=col2Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)
    
 ;RECOVTIME
 LABEL_RECOVTIME = Widget_Label(BASE_COL1, UNAME='LABEL_RECOVTIME'  $
    ,SCR_XSIZE=col1Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='TR(ms):')

 LABEL_RECOVTIME_DATA = Widget_Label(BASE_COL2, UNAME='LABEL_RECOVTIME_DATA'  $
    ,SCR_XSIZE=col2Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)

;ECHOTIME
 LABEL_ECHOTIME = Widget_Label(BASE_COL1, UNAME='LABEL_ECHOTIME'  $
    ,SCR_XSIZE=col1Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='TE(ms):')

 LABEL_ECHOTIME_DATA = Widget_Label(BASE_COL2, UNAME='LABEL_ECHOTIME_DATA'  $
    ,SCR_XSIZE=col2Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)

;NAVG
 LABEL_NAVG = Widget_Label(BASE_COL1, UNAME='LABEL_NAVG'  $
    ,SCR_XSIZE=col1Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='NA:')

 LABEL_NAVG_DATA = Widget_Label(BASE_COL2, UNAME='LABEL_NAVG_DATA'  $
    ,SCR_XSIZE=col2Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)

 ;RARE
 LABEL_RARE = Widget_Label(BASE_COL1, UNAME='LABEL_RARE'  $
    ,SCR_XSIZE=col1Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='ETL:')

 LABEL_RARE_DATA = Widget_Label(BASE_COL2, UNAME='LABEL_RARE_DATA'  $
    ,SCR_XSIZE=col2Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)

 ;FFOV
 LABEL_FFOV = Widget_Label(BASE_COL1, UNAME='LABEL_FFOV'  $
    ,SCR_XSIZE=col1Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='FFOV(cm):')

 LABEL_FFOV_DATA = Widget_Label(BASE_COL2, UNAME='LABEL_FFOV_DATA'  $
    ,SCR_XSIZE=col2Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)

 ;PFOV
 LABEL_PFOV = Widget_Label(BASE_COL1, UNAME='LABEL_PFOV'  $
    ,SCR_XSIZE=col1Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='PFOV(cm):')

 LABEL_PFOV_DATA = Widget_Label(BASE_COL2, UNAME='LABEL_PFOV_DATA'  $
    ,SCR_XSIZE=col2Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)

 ;SFOV
 LABEL_SFOV = Widget_Label(BASE_COL1, UNAME='LABEL_SFOV'  $
    ,SCR_XSIZE=col1Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='SFOV(cm):')

 LABEL_SFOV_DATA = Widget_Label(BASE_COL2, UNAME='LABEL_SFOV_DATA'  $
    ,SCR_XSIZE=col2Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)

 
;=================[COLUMN 3 AND 4 NOW]=============

;SLICE
LABEL_SLICE = Widget_Label(BASE_COL3, UNAME='LABEL_SLICE'  $
  ,SCR_XSIZE=col1Xsize ,SCR_YSIZE=labYsize $
  ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='SLICE:')

LABEL_SLICE_DATA = Widget_Label(BASE_COL4, UNAME='LABEL_SLICE_DATA'  $
  ,SCR_XSIZE=col2Xsize ,SCR_YSIZE=labYsize $
  ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)

 ;SLICE_ORIENT
 LABEL_SLICE_ORIENT = Widget_Label(BASE_COL3, UNAME='LABEL_SLICE_ORIENT'  $
    ,SCR_XSIZE=col3Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='SLICE ORIENT:')

 LABEL_SLICE_ORIENT_DATA = Widget_Label(BASE_COL4, UNAME='LABEL_SLICE_ORIENT_DATA'  $
    ,SCR_XSIZE=col4Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)

; ;SLICE_DIR
; LABEL_SLICE_DIR = Widget_Label(BASE_COL3, UNAME='LABEL_SLICE_DIR'  $
;    ,SCR_XSIZE=col3Xsize ,SCR_YSIZE=labYsize $
;    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='SLICE DIR:')
;
; LABEL_SLICE_DIR_DATA = Widget_Label(BASE_COL4, UNAME='LABEL_SLICE_DIR_DATA'  $
;    ,SCR_XSIZE=col4Xsize ,SCR_YSIZE=labYsize $
;    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)

 ;ADIM
 LABEL_ADIM = Widget_Label(BASE_COL3, UNAME='LABEL_ADIM'  $
    ,SCR_XSIZE=col3Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='ASIZE:')

 LABEL_ADIM_DATA = Widget_Label(BASE_COL4, UNAME='LABEL_ADIM_DATA'  $
    ,SCR_XSIZE=col4Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)


 ;NECHO
 LABEL_NECHO = Widget_Label(BASE_COL3, UNAME='LABEL_NECHO'  $
    ,SCR_XSIZE=col3Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='NECHO:')

 LABEL_NECHO_DATA = Widget_Label(BASE_COL4, UNAME='LABEL_NECHO_DATA'  $
    ,SCR_XSIZE=col4Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)

 ;FATSUP
 LABEL_FATSUP = Widget_Label(BASE_COL3, UNAME='LABEL_FATSUP'  $
    ,SCR_XSIZE=col3Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='FAT SUP:')

 LABEL_FATSUP_DATA = Widget_Label(BASE_COL4, UNAME='LABEL_FATSUP_DATA'  $
    ,SCR_XSIZE=col4Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)

;DIMENSIONS
  LABEL_DIMENSIONS = Widget_Label(BASE_COL3, UNAME='LABEL_DIMENSIONS'  $
    ,SCR_XSIZE=col3Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='DIMENSIONS:')

 LABEL_DIMENSIONS_DATA = Widget_Label(BASE_COL4, UNAME='LABEL_DIMENSIONS_DATA'  $
    ,SCR_XSIZE=col4Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)

;FDIM
 LABEL_FDIM = Widget_Label(BASE_COL3, UNAME='LABEL_FDIM'  $
    ,SCR_XSIZE=col3Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='FSIZE:')

 LABEL_FDIM_DATA = Widget_Label(BASE_COL4, UNAME='LABEL_FDIM_DATA'  $
    ,SCR_XSIZE=col4Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)

 ;PDIM
 LABEL_PDIM = Widget_Label(BASE_COL3, UNAME='LABEL_PDIM'  $
    ,SCR_XSIZE=col3Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='PSIZE:')

 LABEL_PDIM_DATA = Widget_Label(BASE_COL4, UNAME='LABEL_PDIM_DATA'  $
    ,SCR_XSIZE=col4Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)

 ;SDIM
 LABEL_SDIM = Widget_Label(BASE_COL3, UNAME='LABEL_SDIM'  $
    ,SCR_XSIZE=col3Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='SSIZE:')

 LABEL_SDIM_DATA = Widget_Label(BASE_COL4, UNAME='LABEL_SDIM_DATA'  $
    ,SCR_XSIZE=col4Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)


 ;THICK
 LABEL_THICK = Widget_Label(BASE_COL3, UNAME='LABEL_THICK'  $
    ,SCR_XSIZE=col3Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_RIGHT ,VALUE='THICK(mm):')

 LABEL_THICK_DATA = Widget_Label(BASE_COL4, UNAME='LABEL_THICK_DATA'  $
    ,SCR_XSIZE=col4Xsize ,SCR_YSIZE=labYsize $
    ,FRAME=FRAME_ON ,/ALIGN_LEFT ,VALUE='            ', SENSITIVE=0)


;============================[tab bases]=====================
hXsize = hXsize-8
tab_main_base = widget_tab(WID_BASE_MAIN ,UNAME='do_nothing',xsize=hXsize )

zoom_tab = widget_base(tab_main_base,title='Zoom',UNAME='do_nothing')
intensity_tab = widget_base(tab_main_base,title='Intensity',UNAME='do_nothing')
settings_tab = widget_base(tab_main_base,title='Settings',UNAME='do_nothing')
roi_tab = widget_base(tab_main_base,title='ROI',UNAME='do_nothing', column=2)
smoothing_tab = widget_base(tab_main_base,title='Smoothing',UNAME='do_nothing')
display_tab = widget_base(tab_main_base,title='Display',UNAME='do_nothing')

; Tab #5 (display tab) will be the default
WIDGET_CONTROL, SET_TAB_CURRENT=5, tab_main_base

DimensionXSize=(hXsize/3)

;============================[display control widgets]=====================

w = widget_base(display_tab, XOFFSET=280, row=2, /GRID_LAYOUT)

;this is use to display which dimension to traverse when displaying images
SLIDER_DIMENSION_CONTROL = Widget_Button (w ,UNAME='SLIDER_DIMENSION_CONTROL'  $
    ,/ALIGN_CENTER ,FRAME=0  , VALUE='Array dim', SENSITIVE=0  )

;controls the number of images to display ether multiple or single slice display
SLIDER_DISPLAY_CONTROL = Widget_Button (w ,UNAME='SLIDER_DISPLAY_CONTROL'  $
    ,/ALIGN_CENTER  ,FRAME=0 , VALUE='Single  ', SENSITIVE=0)

;display_go forces the display of the widget regardless of the state of auto display
DISPLAY_GO = Widget_Button (w ,UNAME='DISPLAY_GO'  $
    ,/ALIGN_CENTER ,FRAME=0 , VALUE='Display', SENSITIVE=0  )

;this display will try to display the image as a movie.
DISPLAY_MOVIE = Widget_Button (w ,UNAME='DISPLAY_MOVIE'  $
    ,/ALIGN_CENTER  ,FRAME=0 , VALUE='Movie', SENSITIVE=0  )

;selects the desired sdim index
SDIM_SLIDER = Widget_Slider(display_tab, UNAME='SDIM_SLIDER' , MAXIMUM=8  ,SCR_XSIZE=DimensionXSize $
    , MINIMUM = 0, SENSITIVE=0 ,XOFFSET=DimensionXSize*2 ,YOFFSET=53 )
xy_offset = [DimensionXSize*2,83]
if !d.name eq 'X' then xy_offset = [DimensionXSize*2,86]

SDIM_SLIDER_TITLE = Widget_Label(display_tab, UNAME='SDIM_SLIDER_TITLE'  $
    ,XOFFSET=xy_offset[0] , YOFFSET=xy_offset[1] , xsize = 40  $
    ,FRAME=0 ,/ALIGN_LEFT ,VALUE='Slice', SENSITIVE=0)

xy_offset = [DimensionXSize*2+97,83]
if !d.name eq 'X' then xy_offset = [DimensionXSize*2+97,86]

SDIM_SLIDER_MAX = Widget_Label(display_tab, UNAME='SDIM_SLIDER_MAX'  $
    ,XOFFSET=xy_offset[0] , YOFFSET=xy_offset[1], xsize = 40 $
    ,FRAME=0 ,/ALIGN_RIGHT ,VALUE='     ', SENSITIVE=0)



;selects the desired adim index
ADIM_SLIDER = Widget_Slider(display_tab, UNAME='ADIM_SLIDER' , MAXIMUM=8  ,SCR_XSIZE=DimensionXSize $
    , MINIMUM = 0, title= 'Array', SENSITIVE=0 ,XOFFSET=DimensionXSize*2 ,YOFFSET=98);90;,SCR_XSIZE=72 )



w = Widget_Base(display_tab,COLUMN=1 ,/NONEXCLUSIVE,YOFFSET=0)

iImage_flag_button = Widget_Button(w,/ALIGN_LEFT ,VALUE='iImage' $
                    ,UNAME='iImage_flag_button' )


;W_MENU_AUTO_DISPLAY = Widget_Button (w ,UNAME='W_MENU_AUTO_DISPLAY'  $
;    ,/ALIGN_CENTER ,FRAME=0 , VALUE='Auto Display', SENSITIVE=1 )

;WIDGET_CONTROL, /SET_BUTTON, iImage_flag_button

xframe = 1
if !d.name eq 'X' then xframe = 0

slice_base =Widget_Base(display_tab,COLUMN=1, XOFFSET=0, YOFFSET=28, frame=xframe )
void = widget_label(slice_base,value='Slice Image')
w = Widget_Base(slice_base,COLUMN=1 ,/EXCLUSIVE )
freq_phase_button = Widget_Button(w,/ALIGN_LEFT ,VALUE='Freq Phase', SENSITIVE=0 $
                    ,UNAME='freq_phase_button' )
freq_slice_button = Widget_Button(w,/ALIGN_LEFT ,VALUE='Freq Slice', SENSITIVE=0 $
                    ,UNAME='freq_slice_button' )
phase_slice_button = Widget_Button(w,/ALIGN_LEFT ,VALUE='Phase Slice', SENSITIVE=0 $
                    ,UNAME='phase_slice_button' )
WIDGET_CONTROL, freq_phase_button, /SET_BUTTON


rotate_base =Widget_Base(display_tab,COLUMN=1, XOFFSET=100, YOFFSET=4, frame=xframe )
void = widget_label(rotate_base,value='Rotate Image')
w = Widget_Base(rotate_base,COLUMN=1 ,/EXCLUSIVE )
deg_0_button =  Widget_Button(w,/ALIGN_LEFT ,VALUE='0',   SENSITIVE=0 $
                    ,UNAME='deg_0_button' )
deg_90_button = Widget_Button(w,/ALIGN_LEFT ,VALUE='90',  SENSITIVE=0 $
                    ,UNAME='deg_90_button' )
deg_180_button =Widget_Button(w,/ALIGN_LEFT ,VALUE='180', SENSITIVE=0 $
                    ,UNAME='deg_180_button' )
deg_270_button =Widget_Button(w,/ALIGN_LEFT ,VALUE='270', SENSITIVE=0 $
                    ,UNAME='deg_270_button' )
WIDGET_CONTROL, deg_0_button, /SET_BUTTON


flip_base =Widget_Base(display_tab,COLUMN=1, XOFFSET=170, YOFFSET=4, frame=xframe )
void = widget_label(flip_base,value='Flip Image')
w = Widget_Base(flip_base,COLUMN=1 ,/EXCLUSIVE)
flip_none_button = Widget_Button(w,/ALIGN_LEFT ,VALUE='none', SENSITIVE=0 $
                    ,UNAME='flip_none_button' )
flip_horz_button = Widget_Button(w,/ALIGN_LEFT ,VALUE='Horizontal', SENSITIVE=0 $
                    ,UNAME='flip_horz_button' )
flip_vert_button = Widget_Button(w,/ALIGN_LEFT ,VALUE='Vertical', SENSITIVE=0 $
                    ,UNAME='flip_vert_button' )
WIDGET_CONTROL, freq_phase_button, /SET_BUTTON




;============================[zoom control]============================


;zoom slider controls the zoom level of the image
freq_interp_slider  = Widget_Slider(zoom_tab, UNAME='freq_interp_slider', MAXIMUM=8 ,SCR_XSIZE=DimensionXSize $
    , MINIMUM = 1, title= 'Freq Interpolation', SENSITIVE=0  ,XOFFSET=0 ,YOFFSET=0 )

phase_interp_slider = Widget_Slider(zoom_tab, UNAME='phase_interp_slider', MAXIMUM=8 ,SCR_XSIZE=DimensionXSize $
    , MINIMUM = 1, title= 'Phase Interpolation', SENSITIVE=0  ,XOFFSET=DimensionXSize ,YOFFSET=0 )

slice_interp_slider = Widget_Slider(zoom_tab, UNAME='slice_interp_slider', MAXIMUM=8 ,SCR_XSIZE=DimensionXSize $
    , MINIMUM = 1, title= 'Slice Interpolation', SENSITIVE=0  ,XOFFSET=DimensionXSize*2 ,YOFFSET=0 )

w = Widget_Base(zoom_tab, row=1 ,/NONEXCLUSIVE, YOFFSET=50, frame=0)


lock_interp_chk_box = Widget_Button(w, UNAME='lock_interp_chk_box'  $
                   ,/ALIGN_LEFT,VALUE='Lock Interpolation', SENSITIVE=0)

;GA 2008 09 25
down_sample_chk_box = Widget_Button(w, UNAME = 'down_sample_chk_box' $
                                    ,VALUE = 'Down Sample', SENSITIVE=0)

zpad_chk_box = Widget_Button(w, UNAME='zpad_chk_box'  $
                   ,VALUE='Zerofill (x2)', SENSITIVE=0)
                                                       
;WIDGET_CONTROL, /SET_BUTTON, lock_interp_chk_box


x_zoom_slider = Widget_Slider(zoom_tab, UNAME='x_zoom_slider', MAXIMUM=15 ,SCR_XSIZE=DimensionXSize $
    , MINIMUM = 1, title= 'X Zoom', SENSITIVE=0  ,XOFFSET=0 ,YOFFSET=80 )

y_zoom_slider = Widget_Slider(zoom_tab, UNAME='y_zoom_slider', MAXIMUM=15 ,SCR_XSIZE=DimensionXSize $
    , MINIMUM = 1, title= 'Y Zoom', SENSITIVE=0  ,XOFFSET=DimensionXSize ,YOFFSET = 80 )
      

;============================[scaling image intensity]=================




;controls the min of the intensity window
INTENSITY_SLIDER_MIN = Widget_Slider(intensity_tab, UNAME='INTENSITY_SLIDER_MIN'  $
    ,SCR_XSIZE=DimensionXSize, MINIMUM = 0, MAXIMUM=255 ,  title= 'MIN Intensity',VALUE = 0, SENSITIVE=0 $
    ,XOFFSET=0 ,YOFFSET=0)

;sets the intensity of the image
INTENSITY_SLIDER = Widget_Slider(intensity_tab, UNAME='INTENSITY_SLIDER'  $
      ,SCR_XSIZE=DimensionXSize, MINIMUM = 0, MAXIMUM=255 , title= 'Center Intensity',VALUE = 255, SENSITIVE=0 $
      ,XOFFSET=DimensionXSize,YOFFSET=0 )

;sets the max of the intensity window
INTENSITY_SLIDER_MAX = Widget_Slider(intensity_tab, UNAME='INTENSITY_SLIDER_MAX'  $
      ,SCR_XSIZE=DimensionXSize , MINIMUM = 0, MAXIMUM=255 , title= 'MAX Intensity',VALUE = 255 , SENSITIVE=0 $
      ,XOFFSET=DimensionXSize*2 ,YOFFSET=0)

;controls over all scaling based on data value threholds.
min_base = widget_base(intensity_tab, XOFFSET=0 ,YOFFSET=64, /column)
lbl = widget_label(min_base, value="Min Display Value")
DISPLAY_DATA_MIN = widget_text(min_base, uname='DISPLAY_DATA_MIN', /editable, /all_events $
      ,VALUE = '', SENSITIVE=0, xsize=10)

max_base = widget_base(intensity_tab, XOFFSET=DimensionXSize ,YOFFSET=64, /column)
lbl = widget_label(max_base, value="Max Display Value")
DISPLAY_DATA_MAX = widget_text(max_base, uname='DISPLAY_DATA_MAX', /editable, /all_events $
      ,VALUE = '', SENSITIVE=0, xsize=10)

DISPLAY_SCALE_RESET_BUTTON = widget_button(intensity_tab, XOFFSET=2*DimensionXSize, sensitive=1 $
     ,YOFFSET=64+64*0.33, value='Reset thresholds', uname='DISPLAY_SCALE_RESET_BUTTON')


;============================[settings]=================================
f_dim_shift = widget_slider(settings_tab, FRAME=0 ,SENSITIVE=0 $
                        ,MINIMUM=0 ,MAXIMUM=1 ,XOFFSET=0 ,YOFFSET=0 $
                        ,TITLE = 'Frequency Shift' $
                        ,UNAME='f_dim_shift' $
                        ,VALUE=0  ,SCR_XSIZE=DimensionXSize)

p_dim_shift = widget_slider(settings_tab, FRAME=0 ,SENSITIVE=0 $
                        ,MINIMUM=0 ,MAXIMUM=1 ,XOFFSET=DimensionXSize ,YOFFSET=0 $
                        ,TITLE = 'Phase Shift' $
                        ,UNAME='p_dim_shift' $
                        ,VALUE=0 ,SCR_XSIZE=DimensionXSize )

s_dim_shift = widget_slider(settings_tab, FRAME=0 ,SENSITIVE=0 $
                        ,MINIMUM=0 ,MAXIMUM=1 ,XOFFSET=DimensionXSize*2 ,YOFFSET=0 $
                        ,TITLE = 'Slice Shift' $
                        ,UNAME='s_dim_shift' $
                        ,VALUE= 0 ,SCR_XSIZE=DimensionXSize)

w = widget_base(settings_tab,XOFFSET=0 ,YOFFSET=50)

TIME_ZERO_SLIDER = CW_FSLIDER(w , FRAME=0, UNAME='TIME_ZERO_SLIDER' $
                 ,MINIMUM=0.0 ,MAXIMUM=1.0  $
                 ,VALUE=0  $
                 ,TITLE ='Zero Time', /EDIT )

WIDGET_CONTROL,TIME_ZERO_SLIDER, SENSITIVE=0

SLIDER_PREIMAGES = widget_slider(settings_tab, FRAME=0 ,SENSITIVE=0 $
                        ,MINIMUM=1 ,MAXIMUM=4 ,XOFFSET=DimensionXSize ,YOFFSET=64 $
                        ,TITLE = 'Pre Images' $
                        ,UNAME='SLIDER_PREIMAGES' $
                        ,VALUE= 3 ,SCR_XSIZE=DimensionXSize)


w = Widget_Base(settings_tab,COLUMN=1 ,/NONEXCLUSIVE $
          ,XOFFSET=2*DimensionXSize, YOFFSET=64 )

n_avg_flag_button = Widget_Button(w, /ALIGN_LEFT , VALUE='raw / # averages' $
              ,UNAME='n_avg_flag_button',SENSITIVE=0 )

phi_unwrap_flag_button = Widget_Button(w, /ALIGN_LEFT , VALUE='Unwrap Phase' $
              ,UNAME='phi_unwrap_flag_button',SENSITIVE=0 )
              
kspace_settings = Widget_Base(settings_tab, ROW=1, YOFFSET = 125, scr_xsize=3*DimensionXSize)

kspace_fdim_shift = widget_slider(kspace_settings, FRAME=0 ,SENSITIVE=0 $
                        ,MINIMUM=0 ,MAXIMUM=1 ,XOFFSET=0 ,YOFFSET=0 $
                        ,TITLE = 'k-space Frequency Shift' $
                        ,UNAME='kspace_fdim_shift' $
                        ,VALUE=0  ,SCR_XSIZE=DimensionXSize*0.99)

kspace_pdim_shift = widget_slider(kspace_settings, FRAME=0 ,SENSITIVE=0 $
                        ,MINIMUM=0 ,MAXIMUM=1 ,XOFFSET=DimensionXSize ,YOFFSET=0 $
                        ,TITLE = 'k-space Phase Shift' $
                        ,UNAME='kspace_pdim_shift' $
                        ,VALUE=0 ,SCR_XSIZE=DimensionXSize*0.99 )

kspace_sdim_shift = widget_slider(kspace_settings, FRAME=0 ,SENSITIVE=0 $
                        ,MINIMUM=0 ,MAXIMUM=1 ,XOFFSET=DimensionXSize ,YOFFSET=0 $
                        ,TITLE = 'k-space Slice Shift' $
                        ,UNAME='kspace_sdim_shift' $
                        ,VALUE=0 ,SCR_XSIZE=DimensionXSize*0.99 )

kspace_subsamp = Widget_Base(settings_tab,ROW=1, YOFFSET = 180, scr_xsize=3*DimensionXSize)

kspace_fdim_subsamp = widget_slider(kspace_subsamp, FRAME=0 ,SENSITIVE=0 $
                        ,MINIMUM=0 ,MAXIMUM=1 ,XOFFSET=0 ,YOFFSET=0 $
                        ,TITLE = 'k-space Frequency Span' $
                        ,UNAME='kspace_fdim_span' $
                        ,VALUE=0  ,SCR_XSIZE=DimensionXSize*0.99)

kspace_pdim_subsamp = widget_slider(kspace_subsamp, FRAME=0 ,SENSITIVE=0 $
                        ,MINIMUM=0 ,MAXIMUM=1 ,XOFFSET=DimensionXSize ,YOFFSET=0 $
                        ,TITLE = 'k-space Phase Span' $
                        ,UNAME='kspace_pdim_span' $
                        ,VALUE=0 ,SCR_XSIZE=DimensionXSize*0.99 )

kspace_sdim_subsamp = widget_slider(kspace_subsamp, FRAME=0 ,SENSITIVE=0 $
                        ,MINIMUM=0 ,MAXIMUM=1 ,XOFFSET=DimensionXSize ,YOFFSET=0 $
                        ,TITLE = 'k-space Slice Span' $
                        ,UNAME='kspace_sdim_span' $
                        ,VALUE=0 ,SCR_XSIZE=DimensionXSize*0.99 )

;============================[Smoothing Tab]=================================


smooth_base =Widget_Base(smoothing_tab,COLUMN=1, XOFFSET=0, YOFFSET=4, frame=xframe )
void = widget_label(smooth_base,value='Smooth Image')
w = Widget_Base(smooth_base,COLUMN=1 ,/EXCLUSIVE)
smooth_none_button = Widget_Button(w,/ALIGN_LEFT ,VALUE='none', SENSITIVE=0 $
                    ,UNAME='smooth_none_button' )
smooth_button = Widget_Button(w,/ALIGN_LEFT ,VALUE='Smooth', SENSITIVE=0 $
                    ,UNAME='smooth_button' )
median_button = Widget_Button(w,/ALIGN_LEFT ,VALUE='Median', SENSITIVE=0 $
                    ,UNAME='median_button' )
WIDGET_CONTROL, freq_phase_button, /SET_BUTTON


;===========================[ ROI tab ]===================================
;roi_tab

 roi_tab_colum1 = widget_base(roi_tab, frame=0, /column, ypad=0 )
 roi_tab_colum2 = widget_base(roi_tab, frame=0, SCR_XSIZE=hXsize/2-5, column=1)
 w_roi_list = Widget_List(roi_tab_colum1, UNAME='w_roi_list'  $
       ,FRAME=0 , SENSITIVE=1 $
       ,SCR_XSIZE=hXsize/2, SCR_YSIZE=128, Value=['New ROI'] )
       
 roi_xform_base= widget_base(roi_tab_colum1, /nonexclusive)
 w_roi_xform = widget_button(roi_xform_base, value="Don't transform ROIs", uname='w_roi_xform')
 roi_graph_toggle = widget_button(roi_xform_base, value='Display Graph with Curve fit', uname='roi_graph_toggle')
widget_control, roi_graph_toggle, set_button=1

void = Widget_Label(roi_tab_colum2  $
    ,FRAME=0 ,/ALIGN_LEFT ,VALUE='Name', SENSITIVE=1)

w_roi_name = Widget_text(roi_tab_colum2, UNAME='w_roi_name'  $
    ,FRAME=0 ,/ALIGN_LEFT ,VALUE='', SCR_XSIZE=hXsize/3,/EDITABLE)

table = widget_base(roi_tab_colum2, frame=0, column=2)

w_roi_load = Widget_Button(table, UNAME='w_roi_load' , SCR_XSIZE=hXsize/6 $
                  ,/ALIGN_LEFT ,VALUE='Load ', SENSITIVE=1)

w_roi_new_edit = Widget_Button(table, UNAME='w_roi_new_edit' , SCR_XSIZE=hXsize/6 $
                  ,/ALIGN_LEFT ,VALUE='New  ', SENSITIVE=1)

w_roi_save = Widget_Button(table, UNAME='w_roi_save' , SCR_XSIZE=hXsize/6 $
                  ,/ALIGN_LEFT ,VALUE='Save', SENSITIVE=1)

w_roi_delete = Widget_Button(table, UNAME='w_roi_delete' , SCR_XSIZE=hXsize/6 $
                  ,/ALIGN_LEFT ,VALUE='Delete', SENSITIVE=1)

w_roi_mask = WIDGET_DROPLIST( roi_tab_colum2 , UNAME='w_roi_mask', /DYNAMIC_RESIZE  $
         ,value=['None'] ,title='Mask with')

;roi_graph_toggle_base = widget_base(roi_tab_colum2, /nonexclusive)

;===========================[status bar]===================================

;This is where the status of the system is displayed.
;If the version that is running is an older version then don't make the frame sunken
;Else make the frame sunken
if !version.release gt 5.5 then begin
    LABEL_STATUS_BAR = Widget_Label(WID_BASE_MAIN, UNAME='LABEL_STATUS_BAR'  $
    ,SCR_XSIZE=hXsize-5 ,SCR_YSIZE=labYsize , /SUNKEN_FRAME $
    ,FRAME=0 ,/ALIGN_LEFT ,VALUE='Please open a scan')
endif else begin
    LABEL_STATUS_BAR = Widget_Label(WID_BASE_MAIN, UNAME='LABEL_STATUS_BAR'  $
    ,SCR_XSIZE=hXsize-5 ,SCR_YSIZE=labYsize $
    ,FRAME=0 ,/ALIGN_LEFT ,VALUE='Please open a scan')
end


;===========================[menu buttons]===================================
;these are the menu items that are to be displayed.

    W_MENU_FILE = Widget_Button(WID_BASE_MAIN_MBAR, UNAME='W_MENU_FILE'  $
      ,XOFFSET=3 ,YOFFSET=3 ,/MENU ,VALUE='File')

    W_MENU_OPEN = Widget_Button(W_MENU_FILE,  $
      UNAME='W_MENU_OPEN' ,/MENU,VALUE='Open')

    W_MENU_OPEN_BRUKER = Widget_button(W_MENU_OPEN, uname="W_MENU_OPEN_BRUKER", value="Bruker", /menu)
    W_MENU_OPEN_AGILENT = Widget_button(W_MENU_OPEN, uname="W_MENU_OPEN_AGILENT", value="Agilent", /menu)

    W_MENU_SUMMARY_AGILENT = Widget_button(W_MENU_OPEN_AGILENT, uname="W_MENU_SUMMARY_AGILENT", value='Browse Experiment Directory (.fid)')
    W_MENU_OPEN_AGILENT_SINGLE = Widget_button(W_MENU_OPEN_AGILENT, uname="W_MENU_OPEN_AGILENT_SINGLE", value="Single FID Directory (.fid)")
    W_MENU_OPEN_AGILENT_FDF = Widget_button(W_MENU_OPEN_AGILENT, uname="W_MENU_OPEN_AGILENT_FDF", value="FDF Directory (.img)")
    
    W_MENU_OPEN_IMND = Widget_Button(W_MENU_OPEN_BRUKER, UNAME='W_MENU_OPEN_IMND'  $
      ,VALUE='Single Scan')

    W_MENU_OPEN_ONEPULSE = Widget_Button(W_MENU_OPEN_BRUKER, UNAME='W_MENU_OPEN_ONEPULSE'  $
      ,VALUE='OnePulse Spectroscopy')

   	W_MENU_DYNAMICS = Widget_Button(W_MENU_OPEN_BRUKER, $
     UNAME='W_MENU_DYNAMICS' ,/MENU,VALUE='Dynamics')

	; HS - 20061102.
	; Mistake on my part. Will delete it as soon as I can verify it did not
	; break anything.
   	;W_MENU_DIFFUSION = Widget_Button(W_MENU_OPEN, $
     ;UNAME='W_MENU_DIFFUSION' ,/MENU,VALUE='DIFFUSION')

   	W_MENU_OPEN_DYNAMICS = Widget_Button(W_MENU_DYNAMICS, UNAME='W_MENU_OPEN_DYNAMICS'  $
     ,VALUE='DCE MRI')


    W_MENU_OPEN_T2MULTI = Widget_Button(W_MENU_OPEN_BRUKER, UNAME='W_MENU_OPEN_T2MULTI'  $
      ,VALUE='T2 Multi')

    ;; Dec. 2008 - BT - Diffusion expriment using variable b-values
    ;;             over a fixed direction. The scans are stored in a
    ;;             multi-folder format.
    W_MENU_OPEN_ADC_MULTI = Widget_Button(W_MENU_OPEN_BRUKER, UNAME='W_MENU_OPEN_ADCMULTI'  $
      ,VALUE='ADC Multi')

    ;if mas_lite eq 0 then $
       W_MENU_BRUKER_ADT = Widget_button(W_MENU_OPEN_BRUKER, uname="W_MENU_BRUKER_OPEN_ADT", value="ADT Data", /menu)
       W_MENU_OPEN_ADT        = Widget_Button(W_MENU_BRUKER_ADT, UNAME='W_MENU_OPEN_ADT' ,VALUE='DWI Multi Acquisition')
       W_MENU_OPEN_ADT_SINGLE = Widget_Button(W_MENU_BRUKER_ADT, UNAME='W_MENU_OPEN_ADT_SINGLE' ,VALUE='DWI Single Acquisition')

;    if mas_lite eq 0 then $
;        W_MENU_OPEN_PROJECT = Widget_Button(W_MENU_OPEN, UNAME='W_MENU_OPEN_PROJECT'  $
;          ,VALUE='Project')



    ;if mas_lite eq 0 then $
        W_MENU_OPEN_phased_array = Widget_Button(W_MENU_OPEN_BRUKER, UNAME='W_MENU_OPEN_phased_array'  $
          ,VALUE='Phased Array RMS')

    ;if mas_lite eq 0 then $
        W_MENU_OPEN_bruker_image = Widget_Button(W_MENU_OPEN_BRUKER, UNAME='W_MENU_OPEN_bruker_image'  $
          ,VALUE='Single pData')

    ;if mas_lite eq 0 then $
        W_MENU_OPEN_bruker_multi_image = Widget_Button(W_MENU_OPEN_BRUKER, UNAME='W_MENU_OPEN_bruker_multi_image'  $
          ,VALUE='Multi pData')

    ;if mas_lite eq 0 then $
        W_MENU_OPEN_multi_slice_movie = Widget_Button(W_MENU_OPEN_BRUKER, UNAME='W_MENU_OPEN_multi_slice_movie'  $
          ,VALUE='Time Series Scans')

    ;; replaced by generic DICOM reader
    ;if mas_lite eq 0 then $
    ;;    W_MENU_OPEN_DICOM_PHILIPS = Widget_Button(W_MENU_OPEN, UNAME='W_MENU_OPEN_DICOM_PHILIPS'  $
    ;;      ,VALUE='Philips Diffusion DICOM')

    ;if mas_lite eq 0 then $
        W_MENU_OPEN_PHILIPS = Widget_Button(W_MENU_OPEN, value='Philips', UNAME='W_MENU_OPEN_PHILIPS', /menu)
        W_MENU_OPEN_PARREC = Widget_Button(W_MENU_OPEN_PHILIPS, UNAME='W_MENU_OPEN_PARREC'  $
          ,VALUE='Philips PARREC')

        W_MENU_OPEN_SPARSDAT = Widget_Button(W_MENU_OPEN_PHILIPS, UNAME='W_MENU_OPEN_SPARSDAT'  $
          ,VALUE='Philips SPAR/SDAT Spect.')

        W_MENU_OPEN_PHILIPS_DICOM_SPECT = Widget_Button(W_MENU_OPEN_PHILIPS, UNAME='W_MENU_OPEN_PHILIPS_DICOM_SPECT'  $
          ,VALUE='Philips DICOM Spect.')

    ;if mas_lite eq 0 then $
        W_MENU_OPEN_DICOM_GENERIC = Widget_Button(W_MENU_OPEN, UNAME='W_MENU_OPEN_DICOM_GENERIC'  $
          ,VALUE='Generic DICOM Imaging Import')

    ;; Jan 2008 - BT - Nifti 1.1 data import.
        W_MENU_OPEN_NIFTI = Widget_Button(W_MENU_OPEN, UNAME='W_MENU_OPEN_NIFTI'  $
          ,VALUE='NIFTI',/MENU)
        
        W_MENU_OPEN_NIFTI_DATA = Widget_Button(W_MENU_OPEN_NIFTI, UNAME = 'W_MENU_OPEN_NIFTI_DATA',$
          VALUE='Scan Data')
          
        W_MENU_OPEN_NIFTI_MAPS = Widget_Button(W_MENU_OPEN_NIFTI, UNAME = 'W_MENU_OPEN_NIFTI_MAPS',$
          VALUE='QMRI maps')

    void = Widget_Button(W_MENU_OPEN, UNAME='W_MENU_OPEN_FLT'  $
          ,VALUE='*.flt file',/SEPARATOR)

    void = Widget_Button(W_MENU_OPEN, UNAME='W_MENU_OPEN_RAW'  $
          ,VALUE='*.raw file')


    W_MENU_REMOVE_SCAN_MAIN = Widget_Button(W_MENU_FILE, UNAME='W_MENU_REMOVE_SCAN_MAIN'  $
          ,VALUE='Remove Scan', SENSITIVE=1, /MENU)

	W_MENU_REMOVE_SCAN = Widget_Button(W_MENU_REMOVE_SCAN_MAIN, UNAME='W_MENU_REMOVE_SCAN'  $
          ,VALUE='Remove Selected', SENSITIVE=1)

    W_MENU_REMOVE_SCAN_ALL = Widget_Button(W_MENU_REMOVE_SCAN_MAIN, UNAME='W_MENU_REMOVE_SCAN_ALL'  $
          ,VALUE='Remove All', SENSITIVE=1)



    W_MENU_save = Widget_Button(W_MENU_FILE,  $
      UNAME='W_MENU_save' ,/MENU,VALUE='Save/Export')

;    if mas_lite eq 0 then $
;    W_MENU_SAVE_PROJECT = Widget_Button(W_MENU_save, UNAME='W_MENU_SAVE_PROJECT'  $
;      ,VALUE='Save Project', SENSITIVE=0)

    W_MENU_SAVE_TIFF = Widget_Button(W_MENU_save, UNAME='W_MENU_SAVE_TIFF'  $
      ,VALUE='Export to TIFF', SENSITIVE=0)

    W_MENU_SAVE_MPEG = Widget_Button(W_MENU_save, UNAME='W_MENU_SAVE_MPEG'  $
      ,VALUE='Export to MPEG', SENSITIVE=1)

    W_MENU_SAVE_FLT = Widget_Button(W_MENU_save, UNAME='W_MENU_SAVE_FLT'  $
      ,VALUE='Export Preprocessed Images', SENSITIVE=0)

    W_MENU_SAVE_SUBSET = Widget_Button(W_MENU_save, UNAME='W_MENU_SAVE_SUBSET'  $
      ,VALUE='Export Subset Preprocessed Images', SENSITIVE=0)

    menu_save_raw_data = Widget_Button(W_MENU_save, UNAME='menu_save_raw_data'  $
      ,VALUE='Export Raw K-space', SENSITIVE=0)

    W_MENU_SAVE_SLICER = Widget_Button(W_MENU_save, UNAME='W_MENU_SAVE_SLICER'  $
                                       ,VALUE='Export to SLICER 3D', SENSITIVE=0)

    W_MENU_SAVE_NIFTI = Widget_Button(W_MENU_save, UNAME='W_MENU_SAVE_NIFTI'  $
                                      ,VALUE='Export to NIFTI 1.1', SENSITIVE=0)
                                      
    W_MENU_SAVE_DICOM = Widget_Button(W_MENU_save, uname='W_MENU_SAVE_DICOM', $
                                       VALUE='Export to DICOM', SENSITIVE=0)
                                       
;
;    W_MENU_OPTIONS = Widget_Button(W_MENU_FILE, UNAME='W_MENU_OPTIONS'  $
;      ,VALUE='Options')

    W_MENU_MAIN_UPDATE = Widget_Button(W_MENU_FILE, $
      UNAME='W_MENU_MAIN_UPDATE', value="Check for Updates", event_pro='mas_check_update')
      
    W_MENU_MAIN_EXIT = Widget_Button(W_MENU_FILE,  $
      UNAME='W_MENU_MAIN_EXIT' , VALUE='Exit')


    ;======================= analyze menu ============================

    W_MENU_ANALYZE = Widget_Button(WID_BASE_MAIN_MBAR,  $
      UNAME='W_MENU_ANALYZE' ,/MENU ,VALUE='Analyze')

;; Note the curvefit routines were replaced by mas_curvefit_2. I'm leaving
;; this here just in case someone would like to re-enable it later on.
;    W_MENU_CURVEFIT = Widget_Button (W_MENU_ANALYZE ,UNAME='W_MENU_CURVEFIT'  $
;       , VALUE='Curve fit', SENSITIVE=1 ,/MENU)
;    
;    W_MENU_T1_SAT_RECOV_SE = Widget_Button (W_MENU_CURVEFIT ,UNAME='W_MENU_T1_SAT_RECOV_SE'  $
;       , VALUE='T1 Saturation Recovery SE', SENSITIVE=1 ,/MENU )
;
;    W_MENU_T1SR_SE_ROI = Widget_Button (W_MENU_T1_SAT_RECOV_SE ,UNAME='W_MENU_T1SR_SE_ROI'  $
;         , VALUE='ROI', SENSITIVE=0 )   
;     
;    W_MENU_T1SR_SE_IMAGE = Widget_Button (W_MENU_T1_SAT_RECOV_SE ,UNAME='W_MENU_T1SR_SE_IMAGE'  $
;         , VALUE='Image', SENSITIVE=0 )  
;    
;    W_MENU_T1_SAT_RECOV_GE = Widget_Button (W_MENU_CURVEFIT ,UNAME='W_MENU_T1_SAT_RECOV_GE'  $
;       , VALUE='T1 Saturation Recovery GE', SENSITIVE=1 ,/MENU )
;
;      W_MENU_T1SR_GE_ROI = Widget_Button (W_MENU_T1_SAT_RECOV_GE ,UNAME='W_MENU_T1SR_GE_ROI'  $
;         , VALUE='ROI', SENSITIVE=0 )
;
;      W_MENU_T1SR_GE_IMAGE = Widget_Button (W_MENU_T1_SAT_RECOV_GE ,UNAME='W_MENU_T1SR_GE_IMAGE'  $
;         , VALUE='Image', SENSITIVE=0 )
;
;
;    W_MENU_T1_INV_RECOV_SE = Widget_Button (W_MENU_CURVEFIT ,UNAME='W_MENU_T1_INV_RECOV_SE'  $
;       , VALUE='T1 Inversion Recovery SE', SENSITIVE=1 ,/MENU )
;
;    W_MENU_T1IR_SE_ROI = Widget_Button (W_MENU_T1_INV_RECOV_SE ,UNAME='W_MENU_T1IR_SE_ROI'  $
;          , VALUE='ROI', SENSITIVE=0 )
;    
;    W_MENU_T1IR_SE_IMAGE = Widget_Button (W_MENU_T1_INV_RECOV_SE ,UNAME='W_MENU_T1IR_SE_IMAGE'  $
;          , VALUE='Image', SENSITIVE=0 )
;
;    W_MENU_T1_INV_RECOV_GE = Widget_Button (W_MENU_CURVEFIT ,UNAME='W_MENU_T1_INV_RECOV_GE'  $
;       , VALUE='T1 Inversion Recovery GE', SENSITIVE=1 ,/MENU )
;
;    W_MENU_T1IR_GE_ROI = Widget_Button (W_MENU_T1_INV_RECOV_GE ,UNAME='W_MENU_T1IR_GE_ROI'  $
;          , VALUE='ROI', SENSITIVE=0 )
;    
;    W_MENU_T1IR_GE_IMAGE = Widget_Button (W_MENU_T1_INV_RECOV_GE ,UNAME='W_MENU_T1IR_GE_IMAGE'  $
;          , VALUE='Image', SENSITIVE=0 )
;    W_MENU_T2_MULTI_SPIN_ECHO = Widget_Button (W_MENU_CURVEFIT ,UNAME='W_MENU_T2_MULTI_SPIN_ECHO'  $
;       , VALUE='T2 Multi Spin Echo', SENSITIVE=1 ,/MENU )
;
;      W_MENU_T2MSE_ROI = Widget_Button (W_MENU_T2_MULTI_SPIN_ECHO ,UNAME='W_MENU_T2MSE_ROI'  $
;         , VALUE='ROI', SENSITIVE=0 )
;
;      W_MENU_T2MSE_IMAGE = Widget_Button (W_MENU_T2_MULTI_SPIN_ECHO ,UNAME='W_MENU_T2MSE_IMAGE'  $
;         , VALUE='Image', SENSITIVE=0 )
;
;    ;if mas_lite eq 0 then begin
;      W_MENU_ADC = Widget_Button (W_MENU_CURVEFIT ,UNAME='W_MENU_ADC'  $
;        , VALUE='ADC ', SENSITIVE=1 ,/MENU )
;
;      W_MENU_ADC_ROI = Widget_Button (W_MENU_ADC ,UNAME='W_MENU_ADC_ROI'  $
;         , VALUE='ROI (mono-exponent)', SENSITIVE=0 )
;
;      W_MENU_ADC_IMAGE = Widget_Button (W_MENU_ADC ,UNAME='W_MENU_ADC_IMAGE'  $
;         , VALUE='Image (mono-exponent)', SENSITIVE=0 )
;
;      W_MENU_ADC_STRETCH_ROI = Widget_Button (W_MENU_ADC ,UNAME='W_MENU_ADC_STRETCH_ROI'  $
;         , VALUE='ROI (stretch-exponent)', SENSITIVE=0 )
;
;      W_MENU_ADC_STRETCH_IMAGE = Widget_Button (W_MENU_ADC ,UNAME='W_MENU_ADC_STRETCH_IMAGE'  $
;         , VALUE='Image (stretch-exponent)', SENSITIVE=0 )

    ;end
       
    W_MENU_CURVEFIT_2 = Widget_Button (W_MENU_ANALYZE ,UNAME='W_MENU_CURVEFIT_2'  $
       , VALUE='Curve fit 2', SENSITIVE=0)

    W_MENU_IMAGE_STATS = Widget_Button (W_MENU_ANALYZE ,UNAME='W_MENU_IMAGE_STATS'  $
       , VALUE='Image Statistics', SENSITIVE=0 )

    W_MENU_IMAGE_HISTOGRAM = Widget_Button (W_MENU_ANALYZE ,UNAME='W_MENU_IMAGE_HISTOGRAM'  $
       , VALUE='Histogram', SENSITIVE=0 )

    W_MENU_QMRI = Widget_Button (W_MENU_ANALYZE ,UNAME='W_MENU_QMRI'  $
       , VALUE='Quantitative Imaging', SENSITIVE=1, /MENU )
;    W_MENU_QMRI = Widget_Button (W_MENU_ANALYZE ,UNAME='W_MENU_QMRI'  $
;       , VALUE='Quantitative Imaging', SENSITIVE=1 )
           
    W_MENU_QMRI_PD_BTN = WIDGET_BUTTON(W_MENU_QMRI,UNAME='W_MENU_QMRI_PD_BTN'$
       , VALUE = 'Compute PD map', SENSITIVE =1)
       
    W_MENU_QMRI_TOOL_BTN = WIDGET_BUTTON(W_MENU_QMRI, UNAME='W_MENU_QMRI_TOOL_BTN'$
       , VALUE = 'QMRI', SENSITIVE = 1)
       
       W_MENU_MREIT = Widget_Button(W_MENU_ANALYZE, $
       UNAME='W_MENU_MREIT' ,VALUE='MREIT')
       
   W_MENU_DYNAMICS2 = Widget_Button(W_MENU_ANALYZE, $
     UNAME='W_MENU_DYNAMICS' ,/MENU,VALUE='Dynamics')

   W_MENU_DCE_DATA_SELECT = Widget_Button (W_MENU_DYNAMICS2 ,UNAME='W_MENU_DCE_DATA_SELECT'  $
     , VALUE='Select Data', SENSITIVE=1 )
     
   W_MENU_SIGNAL_ENHANCEMENT = Widget_Button (W_MENU_DYNAMICS2 ,UNAME='W_MENU_SIGNAL_ENHANCEMENT'  $
     , VALUE='DCE Images', SENSITIVE=0 )

   W_MENU_DCE_ROI_CURVES = Widget_Button (W_MENU_DYNAMICS2 ,UNAME='W_MENU_DCE_ROI_CURVES'  $
     , VALUE='DCE ROI Curves', SENSITIVE=0 )
   
   W_MENU_DCE_SCC = Widget_Button (W_MENU_DYNAMICS2 ,UNAME='W_MENU_DCE_SCC'  $
     , VALUE='DCE Processing', SENSITIVE=0 )
     
   W_MENU_DCE_VOL = Widget_Button (W_MENU_DYNAMICS2 ,UNAME='W_MENU_DCE_VOL'  $
     , VALUE='DCE Volumes', SENSITIVE=0 )  

   W_MENU_ALG = Widget_Button (W_MENU_ANALYZE ,UNAME='W_MENU_ALG'  $
     , VALUE='Image Algebra', SENSITIVE=0 )

   W_MENU_DIFFUSION2 = Widget_Button(W_MENU_ANALYZE, $
     UNAME='W_MENU_DIFFUSION' ,/MENU,VALUE='Diffusion')

   W_MENU_ADT_GRADIENT = widget_button(W_MENU_DIFFUSION2, uname='W_MENU_ADT_GRADIENT', $
                                       value='Show Gradient Profile')

   W_MENU_ADT_REGRESS = Widget_Button (W_MENU_DIFFUSION2 ,UNAME='W_MENU_ADT_REGRESS'  $
     ,VALUE='ADT Regression', SENSITIVE=0 )

   W_MENU_DIFFTOOLS = Widget_Button(W_MENU_DIFFUSION2, UNAME='W_MENU_DIFFTOOLS', $
      VALUE='Diffusion Visualization Tools', SENSITIVE=0)

   W_MENU_FIBERS = Widget_Button (W_MENU_DIFFUSION2 ,UNAME='W_MENU_FIBERS'  $
     ,VALUE='ADT Fiber Track Mapping', SENSITIVE=0 ) ; Change to 1 if I ever implement an option
     ; to read previously saved trackings from the fibers GUI

   W_MENU_FIBERS_HARDI =  Widget_Button (W_MENU_DIFFUSION2 ,UNAME='W_MENU_FIBERS_HARDI'  $
     ,VALUE='HARDI Fiber Track Mapping', SENSITIVE=1 )

;   W_MENU_SMOOTH = Widget_Button (W_MENU_ANALYZE ,UNAME='W_MENU_SMOOTH'  $
;     ,VALUE='ADT Smoothing & Regression', SENSITIVE=0 )

     W_MENU_FLOW = Widget_Button(W_MENU_ANALYZE, $
     UNAME='W_MENU_FLOW' ,/MENU ,VALUE='Flow')
     
     W_MENU_FLOW_GRADIENT = widget_button(W_MENU_FLOW, uname='W_MENU_FLOW_GRADIENT', $
     value='Show Gradient Profile')

     W_MENU_FLOW_EXTRACT = Widget_Button (W_MENU_FLOW ,UNAME='W_MENU_FLOW_EXTRACT'  $
     ,VALUE='Velocity Extraction' )
     
     W_MENU_GPE = Widget_Button(W_MENU_ANALYZE, $
     UNAME='W_MENU_GPE', VALUE='Gradient Pre-emphasis')
     
    W_MENU_VOLUMETRICS = Widget_Button (W_MENU_ANALYZE ,UNAME='W_MENU_VOLUMETRICS'  $
     ,VALUE='Measure Volume', SENSITIVE=0 )


    W_MENU_Trace_Phase = Widget_Button (W_MENU_ANALYZE ,UNAME='W_MENU_Trace_Phase'  $
     ,VALUE='Trace Phase', SENSITIVE=0 )

     W_MENU_calculate_delta_b = Widget_Button (W_MENU_ANALYZE ,UNAME='W_MENU_calculate_delta_b'  $
     ,VALUE='Calculate Delta B', /MENU )
     
      W_MENU_calculate_delta_b_algebra = Widget_Button (W_MENU_calculate_delta_b ,UNAME='W_MENU_calculate_delta_b_algebra'  $
     ,VALUE='Image Algebra', SENSITIVE=0 )
     
      W_MENU_calculate_delta_b_least_squares = Widget_Button (W_MENU_calculate_delta_b ,UNAME='W_MENU_calculate_delta_b_least_squares'  $
     ,VALUE='Least Squares', SENSITIVE=0 )
     
     W_MENU_CALC_PROTON_DENSITY = widget_button(W_MENU_ANALYZE, uname='W_MENU_CALC_PROTON_DENSITY' $
      ,VALUE='Calculate Proton Density (T1/T2)', sensitive=0)
      
    ;======================= display menu ============================

    W_MENU_DISPLAY = Widget_Button(WID_BASE_MAIN_MBAR,  $
      UNAME='W_MENU_DISPLAY' ,VALUE='Display')

;    W_MENU_CHANGE_LOG = Widget_Button(W_MENU_DISPLAY, UNAME='W_MENU_CHANGE_LOG'  $
;      ,VALUE='MAS Version History', SENSITIVE=1)

    ;; BT Jan 2009 - Data orientation routine
    W_MENU_ORIENT_DATA = Widget_Button(W_MENU_DISPLAY, UNAME='W_MENU_ORIENT_DATA'  $
      ,VALUE='Attach Orientation', SENSITIVE=1)

    W_MENU_COLOR_TABLE = Widget_Button(W_MENU_DISPLAY, UNAME='W_MENU_COLOR_TABLE'  $
      ,VALUE='Load Color Table')

    W_MENU_SLICER = Widget_Button(W_MENU_DISPLAY, UNAME='W_MENU_SLICER'  $
      ,VALUE='Slicer 3d', SENSITIVE=0 )

    W_MENU_volume= Widget_Button(W_MENU_DISPLAY, UNAME='W_MENU_volume'  $
      ,VALUE='Volume', SENSITIVE=0 )

    W_MENU_surface = Widget_Button(W_MENU_DISPLAY, UNAME='W_MENU_surface' $
      ,VALUE='Surface Rendering', SENSITIVE=0 )
    
    W_MENU_MIP = Widget_Button(W_MENU_DISPLAY, UNAME='W_MENU_MIP'  $
      ,VALUE='Maximum Intensity Projection', SENSITIVE=0 )
        
    W_MENU_orthoview = Widget_Button(W_MENU_DISPLAY, UNAME='W_MENU_orthoview'  $
      ,VALUE='Orthogonal Slice View', SENSITIVE=0 )
      
    ;; GA Oct 2009 - Movies in Array dimension for multiple slices at once
    W_MENU_multi_slice_movie = Widget_Button(W_MENU_DISPLAY, UNAME='W_MENU_multi_slice_movie'  $
      ,VALUE='Multi Slice Movie', SENSITIVE=0 )  

    W_MENU_Bruker_Files = Widget_Button(W_MENU_DISPLAY, UNAME='W_MENU_Bruker_Files' $
		,VALUE='Bruker Parameter Files',/MENU, SENSITIVE = 1)

     W_MENU_display_IMND = Widget_Button(W_MENU_Bruker_Files, UNAME='W_MENU_display_IMND'  $
      ,VALUE='Display IMND or Method File', SENSITIVE=0 )

	W_MENU_display_ACQP = Widget_Button(W_MENU_Bruker_Files, UNAME='W_MENU_display_ACQP'  $
      ,VALUE='Display ACQP File', SENSITIVE=0 )
      
  W_MENU_display_CONFIG = Widget_Button(W_MENU_Bruker_Files, UNAME='W_MENU_display_CONFIG'  $
      ,VALUE='Display CONFIG File', SENSITIVE=0 )

    W_MENU_Philips_Files = Widget_Button(W_MENU_DISPLAY, UNAME='W_MENU_Philips_Files' $
		,VALUE='Philips Parameter Files',/MENU, SENSITIVE = 1)

     W_MENU_display_PAR = Widget_Button(W_MENU_Philips_Files, UNAME='W_MENU_display_PAR'  $
      ,VALUE='Display PAR File', SENSITIVE=0 )

    W_MENU_AGILENT_Files = Widget_Button(W_MENU_DISPLAY, UNAME='W_MENU_AGILENT_Files' $
        ,VALUE='Agilent Parameter Files',/MENU, SENSITIVE = 1)

     W_MENU_display_PROCPAR = Widget_Button(W_MENU_AGILENT_Files, UNAME='W_MENU_display_PROCPAR'  $
      ,VALUE='Display PROCPAR File', SENSITIVE=0 )

     W_MENU_display_PROCPAR_Pretty = Widget_Button(W_MENU_AGILENT_Files, UNAME='W_MENU_display_PROCPAR_Pretty'  $
      ,VALUE='Display PROCPAR File (Pretty)', SENSITIVE=0 )

    ;======================= image menu ============================

     W_MENU_IMAGE = Widget_Button(WID_BASE_MAIN_MBAR,  $
      UNAME='W_MENU_DISPLAY' ,VALUE='Image')

    W_MENU_IMAGE_TYPE = Widget_Button(W_MENU_IMAGE , $
      UNAME='W_MENU_IMAGE_TYPE' ,/MENU,VALUE='Type')

    W_MENU_MAGNITUDE_TYPE = Widget_Button(W_MENU_IMAGE_TYPE, UNAME='W_MENU_MAGNITUDE_TYPE'  $
    ,VALUE='Magnitude', SENSITIVE=0)
    W_MENU_REAL_TYPE = Widget_Button(W_MENU_IMAGE_TYPE, UNAME='W_MENU_REAL_TYPE'  $
    ,VALUE='Real', SENSITIVE=0)
    W_MENU_IMAGNARY_TYPE = Widget_Button(W_MENU_IMAGE_TYPE, UNAME='W_MENU_IMAGNARY_TYPE'  $
    ,VALUE='Imaginary', SENSITIVE=0)
    W_MENU_INTENSITY_TYPE = Widget_Button(W_MENU_IMAGE_TYPE, UNAME='W_MENU_INTENSITY_TYPE'  $
    ,VALUE='Intensity', SENSITIVE=0)
    W_MENU_PHASE = Widget_Button(W_MENU_IMAGE_TYPE, UNAME='W_MENU_PHASE'  $
    ,VALUE='Phase', SENSITIVE=0)
    W_MENU_COMPLEX_TYPE = Widget_Button(W_MENU_IMAGE_TYPE, UNAME='W_MENU_COMPLEX_TYPE'  $
    ,VALUE='Complex', SENSITIVE=0)
    W_MENU_RAWM = Widget_Button(W_MENU_IMAGE_TYPE, UNAME='W_MENU_RAWM'  $
    ,VALUE='Raw - Magnitude', SENSITIVE=0)
    W_MENU_RAWR = Widget_Button(W_MENU_IMAGE_TYPE, UNAME='W_MENU_RAWR'  $
    ,VALUE='Raw - Real', SENSITIVE=0)
    W_MENU_RAWI = Widget_Button(W_MENU_IMAGE_TYPE, UNAME='W_MENU_RAWI'  $
    ,VALUE='Raw - Imaginary', SENSITIVE=0)

   W_MENU_IMAGE_PHASE = Widget_Button(W_MENU_IMAGE , $
      UNAME='W_MENU_IMAGE_PHASE', VALUE='Adjust Image Phase')

	;======================= motion correction menu ============================

	; use turn on / turn off menu items



    W_MENU_MOTION_CORRECT = Widget_Button(WID_BASE_MAIN_MBAR,  $
                                          UNAME='W_MENU_MOTION_CORRECT', $
                                          VALUE='Motion Correction')

    W_MENU_SUB_MOTION_CORRECT_OFF = Widget_Button(W_MENU_MOTION_CORRECT, $
                                                 UNAME='W_MENU_SUB_MOTION_CORRECT_OFF',  $
                                                 VALUE='Turn Motion Correction OFF', $
                                                 SENSITIVE=0)

    W_MENU_SUB_MOTION_CORRECT_MI_ON = Widget_Button(W_MENU_MOTION_CORRECT, $
                                                 UNAME='W_MENU_SUB_MOTION_CORRECT_MI_ON',  $
                                                 VALUE='Turn MI Motion Correction ON', $
                                                 SENSITIVE=1)

    W_MENU_SUB_MOTION_CORRECT_ED_ON = Widget_Button(W_MENU_MOTION_CORRECT, $
                                                  UNAME='W_MENU_SUB_MOTION_CORRECT_ED_ON',  $
                                                  VALUE='Turn Edge MC ON', $
                                                  SENSITIVE=1)

    W_MENU_SUB_MOTION_CORRECT_PCA = Widget_Button(W_MENU_MOTION_CORRECT, $
                                                  UNAME='W_MENU_SUB_MOTION_CORRECT_PCA',  $
                                                  VALUE='Principal Component Analysis', $
                                                  SENSITIVE=1)

    ;==================== FILTER menu ===========================

    W_MENU_FILTER = Widget_Button(WID_BASE_MAIN_MBAR,  $
      UNAME='W_MENU_FILTER' ,VALUE='Filter', SENSITIVE=0)

    VOID = Widget_Button(W_MENU_FILTER, UNAME='W_MENU_BACKGROUND_FILTER'  $
        ,VALUE='Background Masking' , EVENT_PRO='mas_background_filter')

    VOID = Widget_Button(W_MENU_FILTER, UNAME='W_MENU_BACKGROUND_MASK_LOAD'  $
        ,VALUE='Load Background Mask' , EVENT_PRO='mas_background_filter_load')

;Frequency Filtering
	W_MENU_FREQUENCY_FILTER = Widget_Button(W_MENU_FILTER, /MENU, UNAME='W_MENU_FREQUENCY_FILTER' $
		,VALUE='Frequency Filtering' ,SENSITIVE=0 )

	W_MENU_FREQUENCY_FILTER_LOWPASS = Widget_Button(W_MENU_FREQUENCY_FILTER, UNAME='W_MENU_FREQUENCY_FILTER_LOWPASS' $
		,VALUE='Low Pass Filter' ,SENSITIVE=0 )

	W_MENU_FREQUENCY_FILTER_HIGHPASS = Widget_Button(W_MENU_FREQUENCY_FILTER, UNAME='W_MENU_FREQUENCY_FILTER_HIGHPASS' $
		,VALUE='High Pass Filter' ,SENSITIVE=0 )

;Image Space filtering

	W_MENU_IMAGE_DOMAIN = Widget_Button(W_MENU_FILTER, UNAME='W_MENU_IMAGE_DOMAIN' $
		,VALUE='Image Filtering' ,SENSITIVE=0 )



; ======================= END OF MAIN MENU DEFINITION ==========


    Widget_Control, /REALIZE, WID_BASE_MAIN

    XManager, 'WID_BASE_MAIN', WID_BASE_MAIN , /NO_BLOCK , CLEANUP='mas_exit';, /JUST_REG



end

;+
; :Description:
;    Procedure "opens" a WWW url using the operating system's
;    mechanisms for such practices. Generally, the OS will open
;    the default browser pointed at the url.
;
; :Params:
;    url - the URL to open
;
;
; :Author: wtriplett
;-
pro mas_open_url, url

    if (!VERSION.OS_FAMILY eq 'Windows') then begin
        spawn, 'start '+url, /hide, /nowait
        ;oJavaDesktop = OBJ_NEW('IDLJavaObject$Static$JAVA_AWT_DESKTOP', 'java.awt.Desktop')
        ;oJavaURI = OBJ_NEW('IDLJavaObject$JAVA_NET_URI', 'java.net.URI', url)
        ;if (oJavaDesktop->isDesktopSupported()) then begin 
        ;    oBrowser = oJavaDesktop->getDesktop(); 
        ;    oBrowser->browse,oJavaURI 
        ;    OBJ_DESTROY, oBrowser 
        ;endif
        ;Remove Objects 
        ;OBJ_DESTROY, oJavaURI, oJavaDesktop 
     endif else if (!VERSION.OS eq 'darwin') then begin
        spawn, 'open '+url
     endif else if (!VERSION.OS eq 'linux') then begin
        spawn, 'xdg-open '+url
     endif

end


;+
; :Description:
;    Checks the MAS download server to see if a newer version
;    exists. This will only work if the release of MAS is performed
;    in a consistent manner. 
;    
;    If the release practices ever change, this will have to be changed also..
;
; :Params:
;    event - even thrown from main menu.
;
;
;
; :Author: wtriplett
;-
pro mas_check_update, event

    common scan_data

    catch, error_code
    if (error_code ne 0) then begin
        catch, /cancel
        if (obj_valid(onet)) then begin
            onet->getProperty, response_code=rc
            if (rc eq 404) then begin
                void = dialog_message('Could not determine current version from web site.', /error, /center)
                obj_destroy, onet
                return
            endif
        endif
        void = dialog_message(['There was a problem connecting to the server', $
                               'or opening a web browser.'], /error, /center)
        if (obj_valid(onet)) then obj_destroy, onet
        return
    endif

    url = project.mas_update_url
    base_dir = strmid(url, 0, strpos(url, 'current_version.txt'))
    
    onet = obj_new('idlneturl', timeout=60)
    current_ver = onet->get(url=url, /string_array)
    onet->GetProperty, response_code=rc
    if (rc eq 404) then begin
    endif else if (rc ne 200) then begin
        void = dialog_message('Could not properly connect to update server.', /error, /center)
        obj_destroy, onet
        return
    endif
    
    if (long(current_ver) gt long(project.mas_version)) then begin
        yesno = dialog_message(["A newer version is available: "+current_ver, $
                                "Would you like to visit the MAS web site?"], /question, /center)
        if (yesno eq 'Yes') then begin
            mas_open_url, base_dir
        endif
        
    endif else begin
        void = dialog_message("Your MAS is up to date.", /center, /info, title="You're up to date.")
    endelse
    
    obj_destroy, onet
    
end


;*****************************************************************************************************
;+
; NAME:
;   mas
;
; PURPOSE:
;   This is the entry point for mas. It was generated by the idl gui creater initally
;   This procedure sets up all the data structures for the whole program.
;   Then calls the procedure to create the widgets
;
; ARGUMENTS:
;
;
; MODIFICATION HISTORY:
;
;-
;*****************************************************************************************************
pro mas, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_
    COMPILE_OPT IDL2
    ;mas_lite is a flag that is used to limit what the user can to do to keep them from getting
    ;into trouble while using mas
    ; 0 is for the full version
    ; 1 is for the lite user limited version
    mas_lite = 0

	;!EXCEPT = 2

;    mas_callback
;    mas_redraw
;    MAS_SAVE
;    mas_sig_enhance
;    mas_extract_files
;    mas_volume
;    mas_background_filter_file
;    mas_roi
;    FSC_Color_file
  ;;  RESOLVE_ROUTINE, 'cw_animate', /EITHER

 	;ITRESOLVE

;    RESOLVE_ALL , /CONTINUE_ON_ERROR


    ;make this optional on window machines.
    device, retain = 2

    ;this is where the data structure for mas lies
    ; i shall call it minnie project.
    COMMON scan_data, project
    COMMON display_data, display_baseID, display_mouseStatus, display_labelID
    display_baseID = 0
    display_mouseStatus = 0

    ;zero  is the time that is going to be used as the zero point in the array
    ;   this may depend on what type of scan the user has defined, for dynamics it is the start of the first
    ;   preimage
    ;min   is the time the acquisition started for the whole multi data set
    ;max   is the time the acquisition finished for the whole multi data set
    ;length    is the actual time each scan took to complete the scan length
    ;array is the actual time units of interest when did each scan start from the middle of the scan.
    time = { zero:0.0, min:0.0, max:0.1, length:ptr_new(), array:ptr_new() }

;    B0             0       B0 field strength (Tesla)
;    fdim           0       dimension for the freq direction
;    pdim           1       dimension for the phase direction
;    sdim           1       dimension for the slice direction
;    adim           1       dimension for the array direction
;   file_adim     ptr_new    is the adim for each file
;
;    sdim_shift     0.0     the shift in the slice dimension
;    pdim_shift     0.0     the shift in the phase dimension
;    fdim_shift     0.0     the shift in the frequency dimensions
;    k_fdim_shift   0.0     the shift in frequency dimension in kspace (pre-fft)
;    k_pdim_shift   0.0     the shift in phase dimension in kspace (pre-fft)
;    k_sdim_shift   0.0     the shift in slice dimension in kspace (pre-fft)
;    k_fdim_span    0L      the effective span in frequency dimension in kspace used for resolution reduction sliders
;    k_pdim_span    0L      the effective span in phase dimension in kspace used for resolution reduction sliders
;    k_sdim_span    0L      the effective span in slice dimension in kspace used for resolution reduction sliders
;    k_fdim_span_max        the max (per-acquisition) span in frequency dimension in kspace
;    k_pdim_span_max        the max (per-acquisition) span in phase dimension in kspace
;    k_sdim_span_max        the max (per-acquisition) span in slice dimension in kspace
    
;    scandate       ''      the date when the scan was taken
;    scanname       ''      the type of pulse sequence used to collect the data
;
;    recovtime      0.0     the recovery time
;    echotime       0.0     the echo time
;    necho          0.0     the number of echos
;    navg           0.0     the number of averages
;    pn_avg      ptr_new() this points to an array where the average for each scan
;                opened is stored. some adt scans have changing number of averages.
;    n_reps                 number of repetitions
;    pn_reps      ptr_new() arrayed number of repetitions
;    rare           0.0     the rare factor
;    f_fov           0.0     the freq field of view
;    p_fov           0.0     the phase field of view
;    s_fov           0.0     ths slice field of view
;    f_voxsz        0.0     voxel dimensions r
;    p_voxsz        0.0     voxel dimensions p
;    s_voxsz        0.0     voxel dimensions s
;    slices         0       the number of slices in the scan
;    thick          0.0     the thickness of a slice
;    fat_sup         0       whether or no fat is supressed 1 = true, 0 = false
;    dimensions     0       the number of dimensions in the scan
;    slice_scheme_ptr  ptr_new()    the way in which the slices were collected
;    file_Path      ''          the directory path to get the slices data
;    displayName    ''          the name to display in the scan list widget
;    orientation    STRARR(3)   the string array of how the slice cuts the sample
;    n_Pre          0           the number of pre-injection slices in diffusion scans
;    flow           'y'         Flow encoding flag ('y' - ON, 'n' - OFF)
;    gflow          0.0         strength of the flow encoding gradient
;    tflow          0.0         half-duration of the flow encoding gradient
;    bptype         't'         shape of the flow encoding gradient ('t' - trapezoid, 's'-sine, 'c'-shifted cosine,'g' - Gaussian)
;    scheme         'h'         Flow encoding scheme
;    fro            ptr_new()   Readout component of flow encoding gradient direction 
;    fpe            ptr_new()   Phase component of flow encoding gradient direction 
;    fsl            ptr_new()   Slice component of flow encoding gradient direction 
;    fc             'y'         flow compensation flag ('y' - ON, 'n' - OFF)
;    bvalArray      ptr_new()   b-values from the scan
;    b_matrix       ptr_new()   matrix contating the b-matrix.
;    repTime_ptr    ptr_new()   repetition time array
;    echoTime_ptr   ptr_new()   echo time array
;    inversionTime_ptr ptr_new() inversion time array
;    epi_num_ref_scans 0     number of reference scans in an EPI experiment (Agilent)
;    
;    imndFile       ''          contains the entire imnd file as 1 string
;    image_type          this flag represents tell what type of scan this is
;                    0   2d,3d
;                    1   dynamic contrast enhancement
;                    2   T2
;                    3   diffusion weighted image
;                    3   apparent diffusion tensor
;                    3   Eigenvalue
;                    3   Eigenvectors
;                    3   fractional anisotropy
;                    4   multi phased array
;                    5   single processed scan
;                    6   multi processed scan
;                    7   spectra
;                    8   chemical shift image
;                    9   variable T1
;                    10  fiber tracking
;                    11  dicom images
;                    12  T2 single scan
;                    13  single slice movie scan
;                    14  multi slice movie scan
;                    15  floating point file
;                    16  Nifti import.
;                    17  17T "DW_Tensor" single scan
;                    18  Generic DICOM
;    multi_scan_file_array  array that contains all the files that are associated with this data set.
;    time                   see structure above
;    acq_matrix   4x4 matrix indicating the transformation from the image space to the logical magnet axis.
;    angle_theta ptr_new()  array of diffusion gradient unit vectors in spherical coordinates
;    angle_phi   ptr_new()  array of diffusion gradient unit vectors in spherical coordinates
;    DICOMtype         Set in mas_dicom_reader.pro
;    DICOMslicelayout  mas_dicom_reader.pro: point to an arrch containing the number of slices per volume
;                      in a multi-volume dicom data set
;    DICOMfiles        mas_dicom_reader.pro: An ptr to array of paths to the actual dicom files for the data set
;    big_delta   ptr_new()  an array containing the "big delta" values for a diffusion experiment
;    small_delta ptr_new()  an array containing the "small delta" values for a diffusion experiment
;    REC_rs_int   Philips par/rec "rescale intercept" parameter used to map stored pixel values to FP or displayed
;                 value ranges
;    REC_rs_slp   "rescale slope" 
;    REC_sc_slp   "scale slope"
;    sform_matrix  sform matrix - when reading a nifti data set, the incoming SFORM matrix is stored here
;    qform_bcd     qform b, c and d - when reading a nifti data set, the incoming QFORM quaternions are stored here
;    qoffset_xyz   The qoffset x, y, z parameters from a nifti data set
;    qfac          the qfac parameter from a nifti data set.
;    Philips_init_fdim - the initial fdim as read from PAR/REC .PAR file, sometimes it is changed when state1 is loaded.
;    Philips_init_pdim - the initial pdim as read from PAR/REC .PAR file, sometimes it is changed when state1 is loaded.
;; Spectroscopy ------------
;    spect_nucleus: string containing the nucleus 
;    spect_spectral_width: The spectral width parameter in Hz
;    spect_acq_size:       The acquisition size in complex points
;    spect_rep_time:       The repetition time in ms
;    spect_pulse_length:   The pulse length in 1e-6 seconds
;    spect_num_avg:        Number of averages
;    spect_bf1:            The "base frequency" in MHz
;    spect_acq_time:       Acquisition time in ms
;    spect_phase_corr      phase correction parameter array: (slope, intercept, 1st order pivot point)
;; State 1 Loader -- This is the name of a procedure which gets called when
;;                   "state1" is loaded.
;   state1_load_procedure  

    imnd = { B0:0.0, gcoil:'undefined', max_slew_rate:0.4615, fdim:0, pdim:1, sdim:1, adim:1, file_adim:ptr_new(), theta:ptr_new(), phi:ptr_new(), psi:ptr_new(), $
             sdim_shift:0.0, pdim_shift:0.0 , fdim_shift:0.0, k_fdim_shift:0.0, k_pdim_shift:0.0, k_sdim_shift:0.0, $
             k_fdim_span: 1L, k_pdim_span: 1L, k_sdim_span: 1L, $
             k_fdim_span_max: 1L, k_pdim_span_max: 1L, k_sdim_span_max: 1L, $
             scan_date:'', scan_name:'', $
             recov_time:0.0, echo_time:0.0, n_echo:0.0, $
             n_avg:0.0, pn_avg:ptr_new(), alpha:90.0, $
             n_reps:0, pn_reps:ptr_new(), $
             rare:1.0, f_fov:0.0, p_fov:0.0, s_fov:0.0, f_fov_offset:0.0, p_fov_offset:0.0, s_fov_offset:0.0, slices:0, thick:0.0,  $
             f_voxsz:0.0, p_voxsz:0.0, s_voxsz:0.0, $
             fat_sup:0, dimensions:0, slice_scheme_ptr:ptr_new(), file_Path:'', $
             display_Name:'',orientation:STRARR(3), n_Pre:3, $
             diff:'n', bval_Array:ptr_new() , b_matrix:ptr_new() , n_bvals:0, $
             flow:'n',gflow:0.0, tflow:0.0, bptype:'t', fc:'y', scheme:'h', ncycl:1, fro:ptr_new(), fpe:ptr_new(), fsl:ptr_new(), $
             tgrdelay:ptr_new(), tgraxis:'s', tgrshape:'t', tgramp:ptr_new(), tgrtime:0.0, tasc:0.0, tdesc:0.0, fpdelay:ptr_new(), $
             BWFresnel:0.0, betaFresnel:0.0, BWSinc:0.0, betaSinc:0.0, etaSinc:0.0, scaleSinc:0.0, $
             rep_Time_ptr:ptr_new(), echo_Time_ptr:ptr_new() , inversion_Time_ptr: ptr_new(), $
             epi_num_ref_scans: 0L, $
             imnd_File:'', image_type:0, multi_scan_file_array:ptr_new(), time:time,$
             acq_matrix:ptr_new(diag_matrix([1.,1.,1.,1.])), $
             angle_theta:ptr_new(), angle_phi:ptr_new(), DICOMtype:0, DICOMslicelayout: ptr_new(), $
             DICOMfiles:ptr_new(), big_delta:ptr_new(), small_delta:ptr_new(), $
             REC_rs_int:PTR_NEW(), REC_rs_slp:PTR_NEW(), REC_sc_slp:PTR_NEW(), $
             sform_matrix:ptr_new(), $
             qform_bcd:ptr_new(), qoffset_xyz:ptr_new(), qfac: 1.0, $
             PARRECtype:0, Philips_init_fdim:0, Philips_init_pdim:0, $
             spect_nucleus: '', spect_spectral_width: 0.0, spect_acq_size: 0L, $
             spect_rep_time: 0.0, spect_pulse_length: 0.0, spect_num_avg: 1L, $
             spect_bf1: 0.0, spect_acq_time: 0.0, spect_phase_corr: [0.0,0,0,0.0], $ ;; [slope,intercept,pivot]
             state1_load_procedure: 'mas_load_state_1_orig',acqtime:0, bandwidth:0, Gfe:0 }


; data array - note these are all pointers unless otherwise indicated.
;   state1  -  loaded "raw" data from whatever source. contents depend on signal selection Image > Type menu
;   state2  -  scaled, rotated, zoomed, windowed image. This is what is displayed to the user
;   signal_enhancement - pointer to signal enhancement image
;   adt     -  apparent diffusion tensor data [*,*,<S0,XX,YY,ZZ,XY,XZ,YZ>]
;   eign_vec - eigenvectors collected from ADT images
;   eign_val - eigenvalues from ADT images
;   frac_ani - fractional anisoptry data collected from ADT fitting
;   avg_dif  - average diffusion
;   flow     - processed flow data <U,V,W,Vres,Image>[*,*,*]
;   flow_ecc_fit - Eddy current correction fit for the three components of velocity vector
;   fit_residuals
;   imgfit1   T1,T2,ADC...etc. Image Fit data produced by fitting routines
;   smooth
;   fibers_adt       Copies of the previous parameters to be used by FIBERS
;   fibers_eign_Vec  The reason for copying is that when mas_rotate_flip is called
;   fibers_frac_Ani  from MAS it acts on the data. When we load Fibers this data would appear
;   fibers_eign_val  rotated.
;   fibers_Avg_Dif
;   dotcoeffs  
;   latt_ani
;   ofmri_coherence
;   curvefit2_volume_fit - Contains the volume fit data from MAS curvefit 2.

    data = { state1:ptr_new(), state2:ptr_new(), $
             signal_Enhancement:ptr_new(), dce_model_fit:ptr_new(), adt:ptr_new(), eign_Vec:ptr_new(), eign_val:ptr_new(), $
             frac_Ani:ptr_new(), Avg_Dif:ptr_new(), flow: ptr_new(), Bez: ptr_new(), fit_residuals:ptr_new(), $
             imgFit1:ptr_new(), smooth:ptr_new(), fibers_adt:ptr_new(), fibers_eign_Vec:ptr_new(), $
             fibers_eign_val:ptr_new(), fibers_frac_Ani:ptr_new(), fibers_Avg_Dif:ptr_new(), $
             dotcoeffs:ptr_new(), latt_Ani:ptr_new(), ofmri_coherence:ptr_new(), $
             curvefit2_volume_fit: ptr_new() }

    ;======================= button parameters ==========================
    ;
    ;sort_dir       0   sort data by which direction
    ;                 [0|1] 0=sort by sdim or equivalent 1=sort by adim
    ;single_Multi_flag 1     display a single image or a multi set of data
    ;                 [0|1] 0=single ,1=multiple
    ;n_avg_flag     0  flag that represents whether or not to divide the raw signal by the number of averages
    ;======================= slider parameters ==========================
    ;
    ;Freq_zoom      1   this is the zoom level in the freq direction
    ;Phase_zoom     1   this is the zoom level in the phase direction
    ;Slice_zoom     1   this is the zoom level in the slice direction
    ;intensity_min     0   all the values of a scan are scaled within the range of [0-255]
    ;                 any value less than this variable will be set to 0 or black
    ;intensity_cen     255       all the values that are being scaled will be set to a max intensity of
    ;                 this variable.
    ;intensity_Max     255       this is the top cut off range for the image any value above this
    ;                 variable will be set to 255 or white.
    ;fdim_start         0       start appropriate actions pertaining to the fdim at this variable
    ;pdim_start         0       start appropriate actions pertaining to the pdim at this variable
    ;sdim_start         0       start appropriate actions pertaining to the sdim at this variable
    ;adim_start         0       start appropriate actions pertaining to the adim at this variable
    ;
    ;======================= menu parameters ============================

    ;flip_Direction       0     direction images are to be flipped over
    ;                 [0|1|2]    0=none, 1=horzontal, 2=vertical
    ;rotate_Direction  0       angle to rotate image to
    ;                 [0|1|2|3] 0=0, 1=90, 2=180 , 3=270
    ;slice_axis         0       axis that will be rotated in the position of x and y
    ;                 [0|1|2]  0=freq phase, 1=freq slice 2=phase slice
    ;signal_type     0       the type of signal that is to processed
    ;                 [0|1|2|3|4] 0=magnitude, 1=real, 2=imagnary, 3=intensity, 4=phase
    ;
    ;======================= display state parameters ====================
    ;
    ;state_2      0     flag that represents whether or not a parameter has change in state 2
    ;                 sence the last display
    ;state_1      0     flag that represnets whether or not a parameter has change in state 1
    ;                 sence the last display
    ;
    ;======================= smooth parameters============================
    ;
    ;smooth_Direction	 0		direction images are to be smoothed
    ;					[0|1|2]	   0=none,	1=smooth function,	2=median function
    ;
    ;======================= option window parameters' ===================
    ;
    ;lock_zoom_TOGGLE   1   flag that tells wether or not to zoom all zooms have the same value
    ;pdim_shift_2          0.0     the % of the fov to shift the image a second time in the phase dim
    ;sdim_shift_2          0.0     the % of the fov to shift the image a second time in the slice dim
    ;
    ;======================= signal enhancement ==========================
    ;
    ;min_Ref      0.0        set the minimum reference for the signal enhancement
    ;intensity_max_se  5.0   controls the intensity max for the signal enhancement data.

    ;======================= curve fitting ==========================
    ;
    ;array_select          Holds the array numbers which are selected for fitting
    ;ROI_fit_type        0      flag that controls which type of processing to do on the ROI selected
    ;                           [0|1|2] 0=T1 progresive saturation, 1=T1 Inversion Recovery
    ;                           ,2=T2 Multi Spin Echo , 3=Apperant Diffusion Coefficient
    ;multi_ROI_graphs    0       tell whether or not to make separate curves for each ROI drawn
    ;                           [0|1] 0=single graph with multiple ROI, 1=separate graphs per ROI
    ;image_fit_type      0     flag that controls which type of processing to do on an image
    ;                           see ROI_fit_type
    ;image_fit_auto      0     automatically calculates and displays the T* image fit when
    ;                           an event is produced by the image fit main base
    ;image_fit_batch     0      toggle for batch export of T2 image fit for all slices in data set
    ;S0_min              0.0     where to clip S0, values below this will be set to 0
    ;S0_max              10.0  where to clip S0, values above this will be set to 0
    ;T_min               0.0       where to clip T[1|2], values below this will be set to 0
    ;T_max               10.0   where to clip T[1|2], values above this will be set to 0
    ;image_fit_threshold 0.20  % threshold of max image intensity,
    ;                                 any pixel above this will be processed, else pixel will be set to zero
    ;image_fit_flag       0     has the current scan been processed.
    ;
    ;====================== Flow processing ==========================
    ;flow_eccmask           Pointer containing the eddy current correction ROI masks
    ;flow_proccess_flag 0   has the whole data set been processed
    ;
    ;====================== ADT fitting ==========================
    ;adt_proccess_flag 0          has the whole data set been fit.
    ;adtThreshold          0.05     percent threshold for apparent diffusion imaging
    ;b_correct_flag      0          whether or not the b-correction is needed for this scan
    ;adt_multi_slice    0        whether or not to process a single slice or a the whole data set.
    ;regresion_type     0          what type of regresion to to
    ;                               [0|1] 0=linear, 1=nonLinear
    ;adt_display_type   0        tell what type of adt image the user selected.
    ;                               [0|1|2|3|4] 0=ADT tensor, 1=eignvectors, 2=frac Ani, 3=average Diffusion
    ;                               4=difusion-weighted
    ;ADT_Start           0          what slice in the ADT matrix the user wants.
    ;adt_display_multi 0        whether the user want to output a single image or the sheet of 11 images
    ;EigenVec_Start         0           what slice in the eigen vector matrix the user wants
    ;adt_img_stat_flg_array 0    array of flags that represent whether or not the user wants to calculate
    ;                                   image stats for a specified ADT image
    ;                                   bS0,bXX,bYY,bZZ,bXY,bXZ,bYZ,bAD,bFA,bE1,bE2,bE3
    ;avg_d_max
    ;avg_d_min
    ;off_dia
    ;dia_max
    ;dia_min
    ;fa_max
    ;
    ;==================== volumetrics ===============================
    ;vol_start_slice   1    the slice to start the volumetrics roi from
    ;vol_stop_slice       10    the plce to stop the volumetrics roi at.
    ;
    ;==================== remove background =========================
    ;slice_range   [0,sdim] array to specify the start and stop slice for not removal
    ;threshold     [0,0] array to specify the lowest and higest value to consider
    ;set_thres_to  [0,0] array to specify what to set values outside of the threshold range
    ;disp_slice       0        slice to display for preview
    ;disp_adim     0      adim to display for preview
    ;region_selection -1   what region the users would like to use as the mask
    ;             [-1|0:*] -1 means to do auto selection 0:* choose that region
    ;preview     0    whether or not to make a preview for the user.
    ;remove_background_active 0 flag that tell whether or not to activate the display.
    ;remove_background_flag    0  flag that tell whether or not to use the background removial.
    ;remove_background_masks   ptr_new pointer to hold the masks for the current image.

                ;button parameters
;nr -    procPram = {sort_dir:0 ,single_Multi_flag:0, n_avg_flag:0 ,$
    procPram = {sort_dir:0 ,single_Multi_flag:0, n_avg_flag:1 ,phi_unwrap_flag:0, no_transform_roi: 0, curvefit_graph: 1, $
                $ ;slider parameters
                intensity_min:0  ,intensity_cen:255   ,intensity_max:255, $
                min_display_thresh: 0.0, max_display_thresh: 1.0, $
                fdim_start:0 ,pdim_start:0 ,sdim_start:0 ,adim_start:0,  $
                $ ;menu parameters
                flip_Direction:0 ,rotate_Direction:0 ,slice_axis:0 ,signal_type:0, $
                $ ; orientation data
                orientation_Imin:4, orientation_Jmin:3, orientation_Kmin:2, have_orientation: 0 ,$
                $ ;display state parameters
                state_1:0 ,state_2:0, state1_max:0.0, $
                $ ;image_var
                image_var:0.000000, $
                 $ ;smooth parameters
                smooth_Direction:0, $
                $;Frequency Filtering
                filterflag:'',subtype:'',design:'',fstart:0.0,fstop:0.0,order:1,$
                $               ; Image Domain Filtering
                image_domain_filterflag: '', image_domain_subtype:'',image_domain_parameter:1,$
                $ ;zoom tab parameters
                lock_interp:0 ,down_sample:0,Freq_interp:1 ,Phase_interp:1, Slice_interp:1, $
                zpad_flag:0, x_zoom:1, y_zoom:1, $
                $ ;signal enhancement
                min_Ref:0.0, intensity_max_se:5.0, $
                $ ;curve Fitting
                ROI_fit_type:0 ,multi_ROI_graphs:0 ,S0_min:0.0  , S0_max:10.0, T_min:0.0, T_max:10.0, $
                image_fit_threshold:0.20, image_fit_flag:0 ,image_fit_type:0, image_fit_auto:0, image_fit_batch:0,  $
                $ ;Eddy current mapping
                Bez_proccess_flag:0, Bez_up_stat:0, Bez_smooth:0, Bez_Ithreshold:0,Bez_time_index:0,Bez_alpha_max:2.5e-2,Bez_axis:0,Bez_apodize:0,Bez_zf_stat:0,Bez_Nfit:4, $
                $ ;flow processing
                flow_proccess_flag:0,flow_up_stat:0, flow_ecc_stat:0, flow_Ithreshold:0, flow_Vthreshold:0,flow_ecc_roi1_set_index:0,flow_ecc_roi1_index:0,flow_ecc_roi2_set_index:0,flow_ecc_roi2_index:0,flow_ecc_fit_fn:replicate(0,6),flow_ecc_vel:replicate(0,3),$ 
                $ ;adt parameters                
                adt_proccess_flag:0, adt_multi_slice:1, regression_type:-1, $
                adt_display_type:3, ADT_Start:0 ,adt_display_multi:0, $
                EigenVec_Start:0,adt_img_stat_flg_array:bytArr(12) ,adtThreshold:0.02, $
                avg_d_max:0.002 ,avg_d_min:0.0000 ,off_dia:0.0002 ,dia_max:0.002 ,dia_min:0.0 ,fa_max:1.0, $
                alpha:1.0, beta:1.0, iterations:1, $
                adt_regress_method:1, $
                $ ; DIFFUSION TOOLS state pointer
                difftools_state:ptr_new(), $
                $ ; DOT state pointer (see mas_DOT.pro)
                dot_state:ptr_new(), $
                $ ; ODF state pointer (see mas_hardi_odf.pro)
                odf_state:ptr_new(), $
                $ ; MOW state pointer (see mas_mow_compute.pro)
                mow_state:ptr_new(), $
                $ ; ADT glyph state pointer (see mas_adt_superquadric.pro)
                adt_sp_state: ptr_new(), $
                $; ProbTrack state pointer (see mas_probtrack.pro)
                probtrack_state: ptr_new(), $
                ortho_tlb: -1, $ ; tlb for ortho display (-1 => not open, >1 => tlb of window base)
                marker_list: ptr_new(), $                
                curvefit2_state: ptr_new(), $
                array_select: ptr_new(), $  ; selects which arrays to analyze (T1,T2,ADC,DCE, etc)
                $ ;volumetrics
                vol_start_slice:1, vol_stop_slice:10, $
                $ ; remove background
                threshold:[0.0,0.0], set_thres_to:[0.0,0.0], key_slice:0, $
                key_adim:0, region_selection:-1, remove_background_active:0, $
                remove_background_thres_image:ptr_new(), $
                latt_slice_status:PTR_NEW(), latt_slice_index:-1, use_3D_lattice:0, latt_cur_slice:-1, $
                $ ; motion correction params
                mc_enable:0, $
                $ ; Mutual Infomrmation MC
                mc_mi_rng_low:-1,   $
                mc_mi_rng_hi:1,     $
 ;;               mc_mi_rng_inc:.125, $
                mc_mi_slice:-1,     $
                mc_mi_refidx:0,     $
                mc_mi_binsize:3,    $
                mc_mi_start_mesh:.25,$
                mc_mi_n_passes:20,   $
                mc_mi_search_tolerance: -5, $
                mc_mi_use_translation: 1, $
                mc_mi_use_rotation: 0, $
                mc_mi_use_shear: 1, $
                mc_mi_use_dilation: 1, $
                mc_mi_use_dilation_h: 1, $
                mc_mi_use_dilation_v: 1, $
                mc_mi_correct_dir_s:1, $
                mc_mi_correct_dir_r:0, $
                mc_mi_correct_dir_p:0, $
                mc_mi_min_method: 1,   $ ;; 0 = Refinement, 1 = Powells
                mc_mi_normalized: 0, $
                $ ; Edge MC
                mc_ed_vxl_chk:2, mc_ed_spv:5, mc_ed_save_fldr:'', mc_ed_refidx:0 }
                $

    max_num_roi = 30
    sroi = { ci:0, ni:1, mask:0, xroi_running:0, $
           pdisplay_names:ptr_new(['New ROI']),$
           pROIs:ptrarr(max_num_roi), $
           max_num_roi:max_num_roi, pROIs_volume:ptrarr(max_num_roi)}

    ;CI                 :0                 Current Index in the project to the current scan selected in the list widget
    ;NI                 :0                 Next Index in the project to insert the next scan
    ;mas_update_url  string URL containing the current version as a text file
    ;mas_version     the string version of this version of MAS
    ;scan_list     :STRARR(50)         string array that is diplayed in the list widget
    ;imndArray          :REPLICATE(imnd, 50)     imnd array structure
    ;dataArray          :REPLICATE(data , 50)    data array structure
    ;procPramArray        :REPLICATE(procPram,50) processing parameter array structure
    ;autoDisplay      :1             flag when a scan is selected or opened as to whether or to make a display for the scan
    ;sliderControlDisplay  :1               flag (1 multislice / 0 single slice ) to be displayed
    ;sliderControlDemension:1                  flag (1 adim / 0 sdim) which way to travers the data to be displayed
    ;scanOpenFlag       :0                  flag whether or not a scan is open (1 true / 0 false)
    ;currentPath      :''               file path of last scan opened
    ;max_num_scans        max_num_scans            maximum number of scans that can be opened
    ;algebra1      :ptr_new()  this is used to do image algebra between 2 images this is the first image
    ;algebra2      :ptr_new()  this is used to do image algebra between 2 images this is the second image
    ;alg_k          :0.0   this is the scaler applied to algebra image 2
    ;regions         :ptr_new()   this is the region that is transfered between subsquent calls to mas_xroi
    ;mas_lite        :mas_lite  this is a flag as to whether or not to run mas lite which is a user limited version
    ;iImage_flag       :1         flag that tell whether or not use iImage as the main display routine.
    ;big_endian         : whether platform is big or little endian
    ;multi_display_tlb: tlb used for global multi display window
    ;difftools_tld: tlb used to display diffusion tools GUI
    ;{odf,dot,mow}_tlb: tlb used to display gui for respective visualization method
    ;debug_mode         : not implemented

    max_num_scans = 20
    project = { ci:0,  $
                ni:0,  $
                mas_update_url: 'http://marecilab.mbi.ufl.edu/software/MAS/current_version.txt', $
                mas_version: '20180830', $
                scan_list:STRARR(max_num_scans+1), $
                imndArray:REPLICATE(imnd, max_num_scans+1), $
                dataArray: REPLICATE(data , max_num_scans+1), $
                procPramArray:REPLICATE(procPram,max_num_scans+1), $
                auto_Display:0,  $
                scan_Open_Flag:0, $
                current_Path:'',  $
                max_num_scans:max_num_scans, $
                algebra1:ptr_new(), $
                algebra2:ptr_new(), $
                algebra1_3D:ptr_new(), $
                algebra2_3D:ptr_new(), $
                algebra_result:ptr_new(),$
                algebra_result_3D:ptr_new(), $
                alg_k:0.0, $
                roi:sroi, $
                mas_lite:mas_lite, $
                iImage_flag:0, $
                use_motion_correction:0,$
                multi_display_tlb:-1, $
                difftools_tlb:-1, $
                curvefit2_gui: ptr_new(), $
                dot_tlb:-1, $
                odf_tlb:-1, $
                mow_tlb:-1, $
                big_endian:0, $
                debug_mode:0, $
                deltaB_result:ptr_new(),$
                deltaB_mask:ptr_new(),$
                deltaB_r2:ptr_new(),$
                deltaB_sigma:ptr_new() }

    case !VERSION.ARCH of

        'x86_64': project.big_endian = 0
        'x86_32': project.big_endian = 0
        'x86'   : project.big_endian = 0
        'i386'  : project.big_endian = 0
        'ppc'   : project.big_endian = 1
        'ppc_32': project.big_endian = 1
        'ppc_64': project.big_endian = 1
        'PPC'   : project.big_endian = 1
        'PPC_32': project.big_endian = 1
        'PPC_64': project.big_endian = 1
	else : begin
		project.big_endian = 1 - (byte(1,0,1))[0]
	end

    endcase

    WID_BASE_MAIN, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_, mas_lite
end
