;; $Id$
;; 
; IDL Widget Interface Procedures. This Code is automatically 
;     generated and should not be modified.
;
; Note, I modified the automatically generated code -BT
pro mas_motion_correct_mi_GUI ;;MC_WINDOW_BASE, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_

    COMMON common_widgets
    COMMON MC_MI_GUI_DATA, MC_WINDOW_BASE, $
                           MC_btn_start, MC_btn_cancel, $
                           MC_refidx, MC_rnglow, MC_rnginc, MC_rnghi, $
                           MC_binsz, MC_slice_sel, MC_slice_all, $
                           MC_startmesh, MC_n_passes, LBL_HISTVAL, $
                           MC_use_powells, MC_use_refinement, LBL_TOLERANCE, $
                           MC_search_tolerance, MC_use_translation, MC_use_rotation, $
                           MC_use_shear, MC_use_dilation, MC_use_dilation_h, MC_use_dilation_v, $
                           MC_normalize_mi, MC_correct_s, MC_correct_p, MC_correct_r
    COMMON scan_data

    parent_dir = project.imndArray[project.ci].file_path
    data_file = parent_dir+'.motion_correction.dat'
    exists = file_test(data_file, /READ)

    if (exists ne 0) then begin
        
        ;; Notify user that motion correction parameters exist
        result = dialog_message(['Motion Correction values exist for this scan.', $
                                 'Do you want to apply them?'], $
                                /QUESTION, /center)
        
        if (result eq 'Yes') then begin
            ;; Apply them if they want, and we're done.
            mc_load_saved_values, data_file
            return
        endif
        
    endif

    title  = STRING('MAS Motion Correction Parameters (MI)')
    params = ptr_new(project.procPramArray[project.ci], /NO_COPY)

  ;; making this modal prevents user from changing slice/aidx in main
  ;; mas window - not a good idea
  MC_WINDOW_BASE = Widget_Base(GROUP_LEADER=WID_BASE_MAIN,  $
      UNAME='MC_WINDOW_BASE' ,XOFFSET=425 ,YOFFSET=5 ,/ALIGN_CENTER  $
      ,TITLE=title ,SPACE=10 ,XPAD=3 ,YPAD=3 ,ROW=4 $;;, /MODAL  $
      ,TLB_FRAME_ATTR=9)
  
  BASE_ROW1 = Widget_Base(MC_WINDOW_BASE, UNAME='BASE_ROW1' ,FRAME=1  $
      ,XOFFSET=3 ,YOFFSET=3 ,TITLE='IDL' ,SPACE=1 ,XPAD=3 ,YPAD=3  $
      ,COLUMN=1, /align_center)

;  LBL_SEARCH_RANGE = Widget_Label(BASE_ROW1, UNAME='LBL_SEARCH_RANGE'  $
;      ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT ,VALUE='Search Range')

;  BASE_R1_C2 = Widget_Base(BASE_ROW1, UNAME='BASE_R1_C2' ,XOFFSET=3  $
;      ,YOFFSET=19 ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3 ,COLUMN=4)

;  LBL_LOW = Widget_Label(BASE_R1_C2, UNAME='LBL_LOW' ,XOFFSET=3  $
;      ,YOFFSET=3 ,/ALIGN_LEFT ,VALUE='Low:')
  
;  MC_rnglow = Widget_Text(BASE_R1_C2, UNAME='MC_rnglow' ,XOFFSET=3  $
;      ,YOFFSET=21 ,/EDITABLE ,VALUE=STRCOMPRESS(STRING((*params).mc_mi_rng_low),/REMOVE_ALL) ,XSIZE=3 ,YSIZE=1)
  
;  LBL_HIGH = Widget_Label(BASE_R1_C2, UNAME='LBL_HIGH' ,XOFFSET=31  $
;      ,YOFFSET=3 ,/ALIGN_LEFT ,VALUE='High:')
  
;  MC_rnghi = Widget_Text(BASE_R1_C2, UNAME='MC_rnghi' ,XOFFSET=31  $
;      ,YOFFSET=21 ,/EDITABLE ,VALUE=STRCOMPRESS(STRING((*params).mc_mi_rng_hi),/REMOVE_ALL) ,XSIZE=3 ,YSIZE=1)
  
;  LBL_STEP = Widget_Label(BASE_R1_C2, UNAME='LBL_STEP' ,XOFFSET=61  $
;      ,YOFFSET=3 ,/ALIGN_LEFT ,VALUE='Step:')

;  MC_startmesh = Widget_Text(BASE_R1_C2, UNAME='MC_startmesh' ,XOFFSET=61  $
;      ,YOFFSET=21 ,/EDITABLE ,VALUE=STRCOMPRESS(STRING((*params).mc_mi_start_mesh),/REMOVE_ALL) ,XSIZE=4 ,YSIZE=1)
  
;  LBL_PASSES = Widget_Label(BASE_R1_C2, UNAME='LBL_PASSES'  $
;      ,XOFFSET=92 ,YOFFSET=3 ,/ALIGN_LEFT ,VALUE='Passes:')

;  MC_n_passes = Widget_Text(BASE_R1_C2, UNAME='MC_n_passes' ,XOFFSET=92  $
;      ,YOFFSET=21 ,/EDITABLE ,VALUE=STRCOMPRESS(STRING((*params).mc_mi_n_passes),/REMOVE_ALL) ,XSIZE=4 ,YSIZE=1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  BASE_R1_C3 = Widget_Base(BASE_ROW1, UNAME='BASE_R1_C3' ,XOFFSET=3  $
      ,YOFFSET=65 ,SCR_XSIZE=144 ,SCR_YSIZE=21 ,TITLE='IDL' ,SPACE=3  $
      ,XPAD=3 ,YPAD=3 ,COLUMN=1)
  
  LBL_INCLUDE_TX = Widget_Label(BASE_R1_C3, UNAME='LBL_INCLUDE_TX'  $
      ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT ,VALUE='Include'+ $
      ' Transformations:')
  
  BASE_R1_C4 = Widget_Base(BASE_ROW1, UNAME='BASE_R1_C4' ,XOFFSET=3  $
      ,YOFFSET=87 ,TITLE='IDL' ,SPACE=1 ,XPAD=0 ,YPAD=0 ,ROW=3)

  
  WID_BASE_2 = Widget_Base(BASE_R1_C4, UNAME='WID_BASE_2' ,XOFFSET=3  $
      ,YOFFSET=3 ,TITLE='IDL' ,COLUMN=1 ,/NONEXCLUSIVE)

    MC_use_translation =  Widget_Button(WID_BASE_2,  $
      UNAME='MC_use_translation' ,SENSITIVE=1 , /ALIGN_LEFT  $
      ,VALUE='Translations')

  MC_use_rotation = Widget_Button(WID_BASE_2, UNAME='MC_use_rotation'  $
      ,SENSITIVE=1 ,/ALIGN_LEFT ,VALUE='Rotations')

  MC_use_shear = Widget_Button(WID_BASE_2, UNAME='MC_use_shear'  $
      ,SENSITIVE=1 ,/ALIGN_LEFT ,VALUE='Shears')

  MC_use_dilation = Widget_Button(WID_BASE_2, UNAME='MC_use_dilation'  $
      ,SENSITIVE=1 ,/ALIGN_LEFT ,VALUE='Dilations')

  BASE_DIL = widget_base(BASE_R1_C4, uname='BASE_DIL', XPAD=18, YPAD=1, /nonexclusive)
  MC_use_dilation_h = Widget_Button(BASE_DIL, UNAME='MC_use_dilation_h'  $
      ,SENSITIVE=project.procpramarray[project.ci].mc_mi_use_dilation ,/ALIGN_LEFT, VALUE='Horizontal')
  MC_use_dilation_v = Widget_Button(BASE_DIL, UNAME='MC_use_dilation_v'  $
      ,SENSITIVE=project.procpramarray[project.ci].mc_mi_use_dilation ,/ALIGN_LEFT, VALUE='Vertical')

  WIDGET_CONTROL, MC_use_translation, SET_BUTTON=project.procpramarray[project.ci].mc_mi_use_translation
  WIDGET_CONTROL, MC_use_rotation,    SET_BUTTON=project.procpramarray[project.ci].mc_mi_use_rotation
  WIDGET_CONTROL, MC_use_shear,       SET_BUTTON=project.procpramarray[project.ci].mc_mi_use_shear
  WIDGET_CONTROL, MC_use_dilation,    SET_BUTTON=project.procpramarray[project.ci].mc_mi_use_dilation
  WIDGET_CONTROL, MC_use_dilation_h,  SET_BUTTON=project.procpramarray[project.ci].mc_mi_use_dilation_h
  WIDGET_CONTROL, MC_use_dilation_v,  SET_BUTTON=project.procpramarray[project.ci].mc_mi_use_dilation_v

;;  LBL_MIN_METHOD = Widget_Label(BASE_ROW1, UNAME='LBL_MIN_METHOD'  $
;;      ,XOFFSET=3 ,YOFFSET=160 ,/ALIGN_LEFT ,VALUE='Optimization'+ $
;;      ' Method:')
  
  BASE_R1_C5 = Widget_Base(BASE_ROW1, UNAME='BASE_R1_C5' ,XOFFSET=3  $
      ,YOFFSET=176 ,TITLE='IDL' ,XPAD=3 ,YPAD=3 ,COLUMN=1, /align_center)

;;  WID_BASE_0 = Widget_Base(BASE_R1_C5, UNAME='WID_BASE_0' ,XOFFSET=3  $
;;      ,YOFFSET=3 ,TITLE='IDL' ,COLUMN=1 ,/EXCLUSIVE)
;;
;;  MC_use_powells = Widget_Button(WID_BASE_0, UNAME='MC_use_powells'  $
;;      ,/ALIGN_LEFT ,VALUE='Powells')
;;
;;  MC_use_refinement = Widget_Button(WID_BASE_0,  sensitive=0, $
;;      UNAME='MC_use_refinement' ,/ALIGN_LEFT ,VALUE='Refinement')

  ;; set default button
;;  selected = ((*params).mc_mi_min_method eq 0) ? MC_use_refinement : MC_use_powells
;;  widget_control, selected, /SET_BUTTON

  ;; update text fields (disabled for powells, enabled for refinement)
;  widget_control, MC_rnghi,     SENSITIVE=(1-(*params).mc_mi_min_method)
;  widget_control, MC_rnglow,    SENSITIVE=(1-(*params).mc_mi_min_method)
;  widget_control, MC_startmesh, SENSITIVE=(1-(*params).mc_mi_min_method)
;  widget_control, MC_n_passes,  SENSITIVE=(1-(*params).mc_mi_min_method)

;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;  WID_BASE_1 = Widget_Base(BASE_R1_C5, UNAME='WID_BASE_1' ,XOFFSET=3  $
;      ,YOFFSET=25 ,TITLE='IDL' ,COLUMN=1 ,/EXCLUSIVE)
  
  LBL_TOLERANCE = Widget_Label(BASE_ROW1, UNAME='LBL_TOLERANCE'  $
      ,XOFFSET=3 ,YOFFSET=227 ,/ALIGN_CENTER, XSIZE=100  $
      ,VALUE=STRING(1/10^(float(-(*params).mc_mi_search_tolerance)), FORMAT='(D0.10)') )

  MC_search_tolerance = Widget_Slider(BASE_ROW1,  $
      UNAME='MC_search_tolerance' ,XOFFSET=3 ,YOFFSET=243 ,MINIMUM=2 $
      ,MAXIMUM=9, TITLE='Search Tolerance', /SUPPRESS_VALUE $
      ,SCROLL=1, /DRAG $
      ,VALUE=STRCOMPRESS(STRING(-(*params).mc_mi_search_tolerance),/REMOVE_ALL) )

;
;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  BASE_ROW2 = Widget_Base(MC_WINDOW_BASE, UNAME='BASE_ROW2' ,FRAME=1  $
      ,XOFFSET=3 ,YOFFSET=80 $;SCR_XSIZE=138 ,SCR_YSIZE=66  $
      ,TITLE='IDL' ,SPACE=1 ,XPAD=3 ,YPAD=3 ,COLUMN=1, /align_center)

  LBL_HISTVAL = Widget_Label(BASE_ROW2, UNAME='LBL_HISTVAL' ,XOFFSET=3  $
      ,YOFFSET=3 ,/ALIGN_CENTER ,XSIZE=30, VALUE=STRCOMPRESS(STRING(2^(*params).mc_mi_binsize),/REMOVE_ALL))

  MC_binsz = Widget_Slider(BASE_ROW2, UNAME='MC_binsz' ,XOFFSET=3,/align_center  $
      ,YOFFSET=19 ,SCR_XSIZE=130 ,SCR_YSIZE=44 ,TITLE='Histogram Bin'+ $
      ' Size' ,MINIMUM=0 ,MAXIMUM=7, SCROLL=1, /SUPPRESS_VALUE, /DRAG,$
      VALUE=STRCOMPRESS(STRING((*params).mc_mi_binsize),/REMOVE_ALL))

  BASE_NORM = widget_base(BASE_ROW2, uname='BASE_NORM', /nonexclusive, /align_center)

  MC_normalize_mi = Widget_Button(BASE_NORM, sensitive=1,  $
      UNAME='MC_normalize_mi' ,/ALIGN_LEFT ,VALUE='Use normalized MI')

  WIDGET_CONTROL, MC_normalize_mi,  SET_BUTTON=project.procpramarray[project.ci].mc_mi_normalized

  BASE_ROW3 = Widget_Base(MC_WINDOW_BASE, UNAME='BASE_ROW3' ,FRAME=1  $
      ,XOFFSET=3 ,YOFFSET=156 ,TITLE='IDL' ,SPACE=1 ,XPAD=3 ,YPAD=3  $
      ,ROW=3)
  
  WID_BASE_45 = Widget_Base(BASE_ROW3, UNAME='WID_BASE_45' ,XOFFSET=3  $
      ,YOFFSET=19 ,TITLE='IDL' ,ROW=1 ,/NONEXCLUSIVE)

  MC_correct_s = Widget_Button(WID_BASE_45, uname='MC_correct_s', value="Slice")
  MC_correct_r = Widget_Button(WID_BASE_45, uname='MC_correct_r', value="Read")
  MC_correct_p = Widget_Button(WID_BASE_45, uname='MC_correct_p', value="Phase")
  widget_control, MC_correct_s, set_button=(*params).mc_mi_correct_dir_s
  widget_control, MC_correct_r, set_button=(*params).mc_mi_correct_dir_r
  widget_control, MC_correct_p, set_button=(*params).mc_mi_correct_dir_p
  
  LBL_SLICE_SEL = Widget_Label(BASE_ROW3, UNAME='LBL_SLICE_SEL'  $
      ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_LEFT ,VALUE='Slice Selection')

  WID_BASE_5 = Widget_Base(BASE_ROW3, UNAME='WID_BASE_5' ,XOFFSET=3  $
      ,YOFFSET=19 ,TITLE='IDL' ,COLUMN=1 ,/EXCLUSIVE)
  
  MC_slice_all = Widget_Button(WID_BASE_5, UNAME='MC_slice_all'  $
      ,/ALIGN_LEFT ,VALUE='Correct all slices')

  MC_slice_sel = Widget_Button(WID_BASE_5, UNAME='MC_slice_sel'  $
      ,/ALIGN_LEFT ,VALUE='Correct selected slice')

  ;; set default button
  selected = ((*params).mc_mi_slice ne -1) ? MC_slice_sel : MC_slice_all
  widget_control, selected, /SET_BUTTON

  BASE_ROW4 = Widget_Base(MC_WINDOW_BASE, UNAME='BASE_ROW4'  $
      ,XOFFSET=3 ,YOFFSET=233 ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3  $
      ,COLUMN=2)
  
  MC_btn_cancel = Widget_Button(BASE_ROW4, UNAME='MC_btn_cancel'  $
      ,XOFFSET=3 ,YOFFSET=3 ,/ALIGN_CENTER ,VALUE='Cancel')
  
  MC_btn_start = Widget_Button(BASE_ROW4, UNAME='MC_btn_start'  $
      ,XOFFSET=54 ,YOFFSET=3 ,/ALIGN_CENTER ,VALUE='Start')

  Widget_Control, /REALIZE, MC_WINDOW_BASE

  xmanager, 'mas_motion_correct_mi_GUI', MC_WINDOW_BASE, $
    Event_Handler='mas_motion_correct_mi_GUI_event'

end

PRO mas_motion_correct_mi_GUI_event, Event

    COMMON MC_MI_GUI_DATA
    COMMON common_widgets
    COMMON scan_data

    ci = project.ci

    CASE Event.id OF

        MC_binsz: BEGIN
            WIDGET_CONTROL, MC_binsz, GET_VALUE=TEMP
            project.procPramArray[ci].mc_mi_binsize = TEMP
            WIDGET_CONTROL, LBL_HISTVAL, SET_VALUE=STRCOMPRESS(STRING(2^TEMP),/REMOVE_ALL)

        END

        MC_search_tolerance: BEGIN
            WIDGET_CONTROL, MC_search_tolerance, GET_VALUE=TEMP
            project.procPramArray[ci].mc_mi_search_tolerance = -TEMP
            WIDGET_CONTROL, LBL_TOLERANCE, SET_VALUE=STRING(1*(10^float(-TEMP)), FORMAT='(D0.10)')
        END
        
        MC_slice_sel: BEGIN
            
            WIDGET_CONTROL, SDIM_SLIDER, GET_VALUE=SLICE
            project.procPramArray[ci].mc_mi_slice = SLICE
            
        END
        
        MC_slice_all: BEGIN
            
            project.procPramArray[ci].mc_mi_slice = -1

        END

        MC_correct_s: project.procPramArray[ci].mc_mi_correct_dir_s = Event.select
        MC_correct_p: project.procPramArray[ci].mc_mi_correct_dir_p = Event.select
        MC_correct_r: project.procPramArray[ci].mc_mi_correct_dir_r = Event.select

;        MC_use_powells: BEGIN
;
;            project.procPramArray[ci].mc_mi_min_method = 1
;            WIDGET_CONTROL, MC_rnglow, SENSITIVE=0
;            WIDGET_CONTROL, MC_rnghi, SENSITIVE=0
;            WIDGET_CONTROL, MC_startmesh, SENSITIVE=0
;            WIDGET_CONTROL, MC_n_passes, SENSITIVE=0
;
;        END
;
;        MC_use_refinement: BEGIN
;
;            project.procPramArray[ci].mc_mi_min_method = 0
;            WIDGET_CONTROL, MC_rnglow, SENSITIVE=1
;            WIDGET_CONTROL, MC_rnghi, SENSITIVE=1
;            WIDGET_CONTROL, MC_startmesh, SENSITIVE=1
;            WIDGET_CONTROL, MC_n_passes, SENSITIVE=1
;           
;        END

        MC_use_translation: BEGIN

            project.procPramArray[ci].mc_mi_use_translation = $
              (Event.select) ? 1 : 0
            
        END
        
        MC_use_rotation: BEGIN

            project.procPramArray[ci].mc_mi_use_rotation = $
              (Event.select) ? 1 : 0

        END
        
        MC_use_shear: BEGIN

            project.procPramArray[ci].mc_mi_use_shear = $
              (Event.select) ? 1 : 0

        END

        MC_use_dilation: BEGIN

            project.procPramArray[ci].mc_mi_use_dilation = $
              (Event.select) ? 1 : 0
            WIDGET_CONTROL, MC_use_dilation_h, SENSITIVE=Event.select
            WIDGET_CONTROL, MC_use_dilation_v, SENSITIVE=Event.select

        END

        MC_use_dilation_h: BEGIN

            project.procPramArray[ci].mc_mi_use_dilation_h = $
              (Event.select) ? 1 : 0

        END

        MC_use_dilation_v: BEGIN

            project.procPramArray[ci].mc_mi_use_dilation_v = $
              (Event.select) ? 1 : 0

        END

        MC_normalize_mi: BEGIN

            project.procPramArray[ci].mc_mi_normalized = $
              (Event.select) ? 1 : 0

        END

        MC_btn_start: BEGIN

;            WIDGET_CONTROL, MC_rnglow, GET_VALUE=TEMP
;            RNG_LOW = float(TEMP[0])
;            project.procPramArray[ci].mc_mi_rng_low = RNG_LOW

;            WIDGET_CONTROL, MC_startmesh, GET_VALUE=TEMP
;            START_MESH = float(TEMP[0])
;            project.procPramArray[ci].mc_mi_start_mesh = START_MESH

;            WIDGET_CONTROL, MC_n_passes, GET_VALUE=TEMP
;            N_PASSES = float(TEMP[0])
;            project.procPramArray[ci].mc_mi_n_passes = N_PASSES

;            WIDGET_CONTROL, MC_rnghi, GET_VALUE=TEMP
;            RNG_HI = float(TEMP[0])
;            project.procPramArray[ci].mc_mi_rng_hi  = RNG_HI

            WIDGET_CONTROL, MC_binsz, GET_VALUE=TEMP
            BINSZ = float(TEMP[0])
            project.procPramArray[ci].mc_mi_binsize = BINSZ

            WIDGET_CONTROL, ADIM_SLIDER, GET_VALUE=REF_IDX
            
            if (project.procPramArray[ci].mc_mi_slice ne -1) then begin
                WIDGET_CONTROL, SDIM_SLIDER, GET_VALUE=SLICE
                project.procPramArray[ci].mc_mi_slice = SLICE
            endif

            ;;if (REF_IDX gt 0) then REF_IDX = REF_IDX-1
            REF_IDX = fix(REF_IDX) > 0

            mas_motion_correct_mi,  $
              REF_IDX,              $
              -1    , $;RNG_LOW,              $
              1     , $;RNG_HI,               $
              0.25  , $;START_MESH,           $
              1     , $;N_PASSES,             $
              2^BINSZ,              $
              project.procPramArray[ci].mc_mi_slice, $
              TOL=1.0*10^(float(project.procPramArray[ci].mc_mi_search_tolerance))
            
            WIDGET_CONTROL, MC_WINDOW_BASE, /DESTROY
        END

        MC_btn_cancel: BEGIN

            WIDGET_CONTROL, MC_WINDOW_BASE, /DESTROY
 
        END

        ELSE:
        
    ENDCASE
END

