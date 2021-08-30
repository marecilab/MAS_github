
PRO trial_call_osv_event, event
  
  WIDGET_CONTROL, event.top, get_uvalue=info_ptr
  widget_control, (*info_ptr).inten_slider, get_value = inten
  
  image = (*info_ptr).final_image
  fin_img_ptr = ptr_new(image)
  
  disp_id = WIDGET_INFO(event.top, find_by_uname='osv_disp_btn')
  upd_id = WIDGET_INFO(event.top, find_by_uname='osv_upd_btn')
  
  mas_display_ortho, data_ptr = fin_img_ptr, state = orthostate
  (*orthostate).max_data_thr = inten
  (*info_ptr).ortho_view_state = orthostate

  
END

PRO trial_call_osv_update, event
  
  WIDGET_CONTROL, event.top, get_uvalue=info_ptr
   
  if (ptr_valid((*info_ptr).ortho_view_state)) then begin
      image = (*info_ptr).final_image
      widget_control, (*info_ptr).inten_slider, get_value = inten
      disp_id = WIDGET_INFO(event.top, find_by_uname='osv_disp_btn')
      upd_id = WIDGET_INFO(event.top, find_by_uname='osv_upd_btn')
      orthostate = (*info_ptr).ortho_view_state
      *((*orthostate).data_ptr) = (*info_ptr).final_image
      (*orthostate).max_data_thr = inten
      ortho_refresh_data, orthostate
      ortho_update_imagery, orthostate
  endif
  
END

;  CASE event.id OF 
;  
;  disp_id: BEGIN
;             mas_display_ortho, data_ptr = fin_img_ptr, state = state
;             END
;  upd_id : BEGIN
;             ortho_update_imagery, state
;             ortho_refresh_data, state
;             END
;     ELSE: ;;
;  ENDCASE
             
  
  
  

;END

;
;
; EVENT PROCEDURE TO SAVE THE FINAL IMAGE AFTER SELECTING THE TR,TE & OTHER PARMS'.
;
;

pro save_data_rat_event, event
      
      widget_control, event.top, get_uvalue=info_ptr
      button_id = widget_info((*info_ptr).button_id,find_by_uname='Save_button')
      image = (*info_ptr).final_image
      fin_img_ptr = ptr_new(image)
      filename = 'Final_image'+'.nii.gz' 
      mas_export_nifti, data_ptr = fin_img_ptr, $
                        file_name = filename
      cd, current=dest_dir
      mesg = dialog_message(['The image has been saved at:',dest_dir],/information, /center)
      
      psequence = widget_info((*info_ptr).pulseseq,/droplist_select)
      seqname = (*info_ptr).seqname
      
      WIDGET_CONTROL, (*info_ptr).TRslider, get_value = TRvalue 
      WIDGET_CONTROL, (*info_ptr).TEslider, get_value = TEvalue
      WIDGET_CONTROL, (*info_ptr).bval_slider, get_value = bvalue   
      WIDGET_CONTROL, (*info_ptr).angle_slider, get_value = tip_angle
      WIDGET_CONTROL, (*info_ptr).InvRslider, get_value = invtime
      WIDGET_CONTROL, (*info_ptr).inten_slider, get_value = max_inten
      
      param_array = STRARR(10)
      param_array[0] = ('***** Imaging Parameters ***** :')
      param_array[1] = ('Pulse Sequence = ') + strtrim(seqname[psequence],2)
      param_array[2] = 'TR (sec) = ' + strtrim(string(TRvalue),2)
      param_array[3] = 'TE (sec) = ' + strtrim(string(TEvalue),2)
      param_array[4] = 'Bvalue (sec/mm2)= ' + strtrim(string(bvalue),2)
      param_array[5] = 'Tip Angle (deg) = ' + strtrim(string(tip_angle),2)
      param_array[6] = 'Inversion Time (sec)= ' + strtrim(string(invtime),2)
      param_array[7] = 'Max Intensity = ' + strtrim(string(max_inten),2)
      param_array[8] = 'Min Intensity = ' + strtrim('0',2)
      param_array[9] = string(10B)
      
      fpath = DIALOG_PICKFILE( Title = 'Select Directory')
      OpenW, lun,fpath + '.txt', /get_lun
      PrintF, lun,format='(a50)', param_array
      close, lun
      free_lun, lun
      
end

;
;
;--------------------FUNCTION TO CLEAN THE DATA, REMOVE INFs and NANs IN THE DATASETS.
;
;

function clean_image, image
      
      inf_check = where(finite(image, /infinity))
      nan_check = where(finite(image, /nan))
      neg_nan_check = where(finite(image,/nan, sign=-1))
        if (inf_check eq -1) then begin
            print, 'No infinite values'
        endif else begin
            image[inf_check] = 0
            print, 'Image cleaned'
        endelse

        if(nan_check eq -1) then begin
            print, 'No NAN values'
        endif else begin
            image[nan_check]=0
            print, 'Image Cleaned'
        endelse
        
        if(neg_nan_check eq -1) then begin
            print, 'No -NAN values'
        endif else begin
            image[nan_check]=0
            print, 'Image Cleaned'
        endelse
   return, image
   
end

;
;
; EVENT PROCEDURE TO HANDLE THE EVENTS OCCURING WITH SLIDERS.
;
;      

pro trial_qmri_tool_event, event
    common scan_data
    
    widget_control, event.top, get_uvalue=info_ptr
    
    osv_flag = WIDGET_INFO(event.top, find_by_uname='osv_btn')   
    value1 = WIDGET_INFO((*info_ptr).reslice_btn1, /BUTTON_SET)
    value2 = WIDGET_INFO((*info_ptr).reslice_btn2, /BUTTON_SET)
    value3 = WIDGET_INFO((*info_ptr).reslice_btn3, /BUTTON_SET)

    case 1 of 
      value1: widget_control, (*info_ptr).slice_slider , SET_SLIDER_MAX=(*info_ptr).sdim-1     
      value2: widget_control, (*info_ptr).slice_slider , SET_SLIDER_MAX=(*info_ptr).pdim-1
      value3: widget_control, (*info_ptr).slice_slider , SET_SLIDER_MAX=(*info_ptr).fdim-1
      else: return
    endcase
    
    widget_control, (*info_ptr).TRslider, get_value = TR_value
    widget_control, (*info_ptr).TEslider, get_value = TE_value
    widget_control, (*info_ptr).bval_slider, get_value = b_value
    widget_control, (*info_ptr).slice_slider, get_value = slice_pos
    widget_control, (*info_ptr).inten_slider, get_value = inten
    widget_control, (*info_ptr).reslice_btn1, get_value = reslice_btn1
    widget_control, (*info_ptr).reslice_btn2, get_value = reslice_btn2
    widget_control, (*info_ptr).reslice_btn3, get_value = reslice_btn3
       
    seqindex = widget_info((*info_ptr).pulseseq,/droplist_select)
    
    ;-------------------------------------------RE-PARAMETERIZATION---------------------------------------------
    T1 = (*info_ptr).T1
    T2 = (*info_ptr).T2
    rho = (*info_ptr).Proton
    D = (*info_ptr).AD
    diff_fac = exp(-b_value*D)
    
    TR = TR_value + make_array((*info_ptr).fdim,(*info_ptr).pdim,(*info_ptr).sdim)
    TE = TE_value + make_array((*info_ptr).fdim,(*info_ptr).pdim,(*info_ptr).sdim)
    
    ;----------------------------------------------SPIN ECHO PULSE SEQUENCE-------------------------------------
    if(seqindex eq 1) then begin
        WIDGET_CONTROL, (*info_ptr).angle_slider, sensitive=0
        A = (1.0 - exp(-TR/T1))
        B = exp(-TE/T2)
        S = rho*A*B*diff_fac
        (*info_ptr).final_image = S
        CI = project.ci
        rotation = project.procpramarray[CI].rotate_direction
        Sscaled=bytscl(S,min=0,max=inten)
        (*info_ptr).fin_imagesc = Sscaled
        new_img_ptr = ptr_new(Sscaled)
        IF(value1 eq 1) THEN BEGIN
                ERASE
                Sfinal = clean_image(reform(Sscaled[*,*,slice_pos]))
;                tv, rotate(Sfinal,rotation);, NOINTERPOLATION, KEEP_ASPECT_RATIO
                tv, Sfinal
        ENDIF
        IF(value2 eq 1) THEN BEGIN
                ERASE
                Sfinal = clean_image(reform(Sscaled[*,slice_pos,*]))
                tv, rotate(Sfinal,rotation);, NOINTERPOLATION, KEEP_ASPECT_RATIO
        ENDIF
        IF(value3 eq 1) THEN BEGIN
                ERASE
                Sfinal = clean_image(reform(Sscaled[slice_pos,*,*]))
                tv, rotate(Sfinal,rotation);, NOINTERPOLATION, KEEP_ASPECT_RATIO
        ENDIF
    endif
    
;-----------------------------------------------SPIN ECHO - INVERSION RECOVERY-------------------------------------
     if(seqindex eq 3) then begin
        widget_control, (*info_ptr).InvRslider, get_value = InvRTime
        TI = InvRTime + make_array((*info_ptr).fdim,(*info_ptr).pdim,(*info_ptr).sdim)
        A = abs((1 - 2*exp(-TI/T1)+ exp(-TR/T1)))
        B = exp(-TE/T2)
        S = rho*A*B*diff_fac
        (*info_ptr).final_image = S
        CI = project.ci
        rotation = project.procpramarray[CI].rotate_direction
        Sscaled=bytscl(S,min=0,max=inten)
        (*info_ptr).fin_imagesc = Sscaled
        IF(value1 eq 1) THEN BEGIN
                ERASE
                Sfinal = clean_image(reform(Sscaled[*,*,slice_pos]))
;                tv, rotate(Sfinal,rotation);, NOINTERPOLATION, KEEP_ASPECT_RATIO
                tv, Sfinal
        ENDIF
        IF(value2 eq 1) THEN BEGIN
                ERASE
                Sfinal = clean_image(reform(Sscaled[*,slice_pos,*]))
                tv, rotate(Sfinal,rotation);, NOINTERPOLATION, KEEP_ASPECT_RATIO
        ENDIF
        IF(value3 eq 1) THEN BEGIN
                ERASE
                Sfinal = clean_image(reform(Sscaled[slice_pos,*,*]))
                tv, rotate(Sfinal,rotation);, NOINTERPOLATION, KEEP_ASPECT_RATIO
        ENDIF
     endif 
;-------------------------------------------------GRADIENT ECHO SEQUENCE----------------------------------------------
    if (seqindex eq 2) then begin
        widget_control, (*info_ptr).angle_slider, get_value = angle
        pi = 3.1416
        theta = (angle*pi)/180 + make_array((*info_ptr).fdim,(*info_ptr).pdim,(*info_ptr).sdim)
        E1 = exp(-TR/T1)
        A = rho*sin(theta)
        Z = 1 - (E1*cos(theta))
        B = (1 - E1)/Z
        C = exp(-TE/T2)
        S = A*B*C*diff_fac
        (*info_ptr).final_image = S
        Sscaled = bytscl(S,min=0,max=inten)
        (*info_ptr).fin_imagesc = Sscaled
        new_img_ptr = ptr_new(Sscaled)
        WIDGET_CONTROL, (*info_ptr).angle_slider, sensitive=1
;        (*info_ptr).fin_img_ptr = new_img_ptr
        CI = project.ci
        rotation = project.procpramarray[CI].rotate_direction
        IF(value1 eq 1) THEN BEGIN
                ERASE
                Sfinal = clean_image(reform(Sscaled[*,*,slice_pos]))
;                tv, rotate(Sfinal,rotation);, NOINTERPOLATION, KEEP_ASPECT_RATIO
                tv, Sfinal
        ENDIF
        IF(value2 eq 1) THEN BEGIN
                ERASE
                Sfinal = clean_image(reform(Sscaled[*,slice_pos,*]))
                tv, rotate(Sfinal,rotation);, NOINTERPOLATION, KEEP_ASPECT_RATIO
        ENDIF
        IF(value3 eq 1) THEN BEGIN
                ERASE
                Sfinal = clean_image(reform(Sscaled[slice_pos,*,*]))
                tv, rotate(Sfinal,rotation);, NOINTERPOLATION, KEEP_ASPECT_RATIO
        ENDIF
    endif

    

;------------------------------RESLICE OPTIONS---------------------------------
;  
;    IF (event.id eq (*info_ptr).reslice_btn) THEN BEGIN
;        img = (*info_ptr).final_image
;        CASE event.value OF
;            
;            0:  BEGIN
;                img_sc = bytscl(reform(img[*,*,slice_pos]),min=0,max=inten)
;                CI = project.ci
;                rotation = project.procpramarray[CI].rotate_direction
;                erase
;                tv, clean_image(rotate(img_sc,rotation))
;                END
;            1:  BEGIN
;                img_sc = bytscl(reform(img[*,slice_pos,*]),min=0,max=inten)
;                CI = project.ci
;                rotation = project.procpramarray[CI].rotate_direction
;                erase
;                tv, clean_image(rotate(img_sc,rotation))
;                END
;            2:  BEGIN
;                img_sc = bytscl(reform(img[slice_pos,*,*]),min=0,max=inten)
;                CI = project.ci
;                rotation = project.procpramarray[CI].rotate_direction
;                erase
;                tv, clean_image(rotate(img_sc,rotation))
;                END
;         ELSE:  ;;
;        ENDCASE
;    ENDIF

    trial_call_osv_update, event
    
end

;
;
; EVENT PROCEDURE FOR HANDLING MOUSEOVER & BUTTON EVENTS FROM THE DRAW WINDOW
;
;

pro mouseover_draw_event, event

    WIDGET_CONTROL, event.top, get_uvalue=info_ptr
    WIDGET_CONTROL, (*info_ptr).x_text, get_value= x_loc
    WIDGET_CONTROL, (*info_ptr).y_text, get_value= y_loc
    WIDGET_CONTROL, (*info_ptr).inten_text, get_value = pix_int
    WIDGET_CONTROL, (*info_ptr).slice_slider, get_value = slice_val

    value1 = WIDGET_INFO((*info_ptr).reslice_btn1, /BUTTON_SET)
    value2 = WIDGET_INFO((*info_ptr).reslice_btn2, /BUTTON_SET)
    value3 = WIDGET_INFO((*info_ptr).reslice_btn3, /BUTTON_SET)

    case 1 of 
      value1: begin
          x_loc = event.x 
          y_loc = event.y 
          z_loc = slice_val
      end
      
      value2: begin
          x_loc = event.x
          z_loc = event.y 
          y_loc = slice_val
      end
      
      value3: begin
          y_loc = event.x
          z_loc = event.y
          x_loc = slice_val
      end
      else: return
      
    endcase
   
    image = (*info_ptr).final_image 
    
;    x_loc = strtrim(event.x)
;    y_loc = strtrim(event.y)

    if (x_loc ge (*info_ptr).fdim or y_loc ge (*info_ptr).pdim or z_loc ge (*info_ptr).sdim) then begin
        pix_int = '---'
    endif else begin
        pix_int = strtrim(image[x_loc,y_loc,z_loc])
    endelse
    
    print, x_loc,y_loc,z_loc,pix_int
    WIDGET_CONTROL, (*info_ptr).x_text, set_value= strtrim(x_loc,2)
    WIDGET_CONTROL, (*info_ptr).y_text, set_value= strtrim(y_loc,2)
    WIDGET_CONTROL, (*info_ptr).inten_text, set_value = pix_int
    WIDGET_CONTROL, event.top, set_uvalue=info_ptr
    
end


pro zoom_on_click_event,event

;    COMMON scan_data, project
;    CI = project.ci
;    
;    WIDGET_CONTROL, event.top, get_uvalue=info_ptr
;    WIDGET_CONTROL, (*info_ptr).slice_slider, get_value = z_value
;    fdim = (*info_ptr).fdim
;    pdim = (*info_ptr).pdim
;    image = (*info_ptr).final_image
;    rotation = project.procpramarray[CI].rotate_direction
;    array = rotate(reform(image[*,*,z_value]),rotation)
;;    test_widzoom, array, info_ptr
;;    im_ptr = ptr_new(array)
;;    tv, array

end

pro draw_widget_event,event

    WIDGET_CONTROL, event.top, get_uvalue = info_ptr
    
    CASE event.type OF 
    
    0 : zoom_on_click_event,  event
    2 : mouseover_draw_event, event
    ELSE : ;;
    
    ENDCASE

end

;---------------------------------------------------CLEAN UP POINTERS---------------------------------------------
pro trial_qmri_tool_cleanup, id

  print, 'Cleaning Up'
  widget_control, id, get_uvalue=info_ptr
  ptr_free, info_ptr
  print, 'Done'

end
;
;
;-------------------------------------------------------MAIN-------------------------------------------------------
;
;
pro trial_qmri_tool,T1,T2,Proton,AD,S0,fdim,pdim,sdim,infoptr

  COMMON scan_data, project
       
       T1 = T1
       T2 = T2
       AD = AD
       S0 = S0
       Proton = Proton
       fdim = fdim
       pdim = pdim
       sdim = sdim 
  
;-------------------------------Creating the Widgets--------------------------

   topbase = widget_base(title='Welcome to the World of Quantitative MRI',row=1)
              
 ; - Creating the display window
        draw_base = widget_base(topbase, /frame, column=1)
        draw_window = widget_draw(draw_base, scr_xsize=fdim, scr_ysize=pdim,$
;        draw_window = widget_draw(draw_base, scr_xsize=300, scr_ysize=300,$
                      /motion_events,/button_events,event_pro='draw_widget_event')
        widget_control, draw_window, get_value=draw_id
;        wset, draw_id
 
 ; - Creating the Controls window
        tab_base = widget_base(topbase,/column)
        tab_widget = widget_tab(tab_base,uname='Tab Widget')
        control_base = widget_base(tab_widget,/frame,/column,title='Controls')
       
 ; - Creating the dropdown to hold the names of various pulse sequences
        seqname = strarr(4)
        seqname = ['--------', 'Spin Echo', 'Gradient echo', 'SE - Inv Recovery']
        seq_base = widget_base(control_base, /frame, column=1)
        pulseseq = widget_droplist(seq_base,title = ' Choose a Pulse Sequence', value=seqname,$
                            /dynamic_resize)
        widget_control, pulseseq, get_value=sequence
        
        angle_base = widget_base(seq_base,/row)
        angle_slider = widget_slider(angle_base,title='Tip Angle (degrees)',/drag,/sensitive, $
                        minimum=0,value=90,maximum = 180)
        
        coord_base = widget_base(draw_base,/column,/frame, /base_align_right)
        label_base = widget_base(coord_base,/row)
        inten_label = widget_label(label_base,xoffset=10,value='Intensity:')
        inten_text = widget_text(label_base,xoffset=40, value='0')
        x_base = widget_base(coord_base,/row)
        x_label = widget_label(x_base,xoffset=10, value = 'X value:')
        x_text = widget_text(x_base,xoffset=40,uname='x_text_widget',value='0')
        y_base = widget_base(coord_base,/row)
        y_label = widget_label(y_base,xoffset=10,value = 'Y value:')
        y_text = widget_text(y_base,xoffset=40,uname='y_text_widget',value='0')
        
        TR_base = widget_base(control_base,/row)
        TRslider = cw_fslider(TR_base, title = 'TR(sec)', /drag, /edit, minimum = 0, $
                      value=0, maximum = 10)
        TE_base = widget_base(TR_base)
        TEslider = cw_fslider(TE_base, title = 'TE(sec)', /drag, /edit, minimum = 0, $
                     value=0,maximum = 0.5)

        bval_base = widget_base(control_base,/row)
;        bval_slider = cw_fslider(bval_base,title='bvalue',minimum=0, /edit,$)
;              maximum=5000,scroll=100,value=0,/drag)
              
        bval_slider = widget_slider(bval_base,title='bvalue',minimum=0, $)
              maximum=5000,scroll=100,value=0,/drag, /sensitive)
        
        
              
        InvR_base = widget_base(bval_base)
        InvRslider = cw_fslider(InvR_base, title = 'InvTime(sec)', /drag, /edit, minimum = 0, $
                      value=0, maximum = 10)
 
        slice_base = widget_base(control_base,/row)
        slice_slider = widget_slider(slice_base, title='Slice number',/drag,$
                       value=4,minimum = 0, maximum = sdim-1, scroll=1)        
                        
        inten_base = WIDGET_BASE(tab_widget,title='Display',/column)
        inten_slider = CW_FSLIDER(inten_base,title='Max intensity',/drag,/edit,minimum =0, $
                        value=50, maximum=5000)
        
;        overlay_lbl_base = WIDGET_BASE(inten_base,/column)
;        overlay_lbl = WIDGET_LABEL(overlay_lbl_base,value='OVERLAY IMAGERY')
;        overlay_base = WIDGET_BASE(overlay_lbl_base,/column,/frame)
;        overlay_btn_base = WIDGET_BASE(overlay_base)
;        load_img = WIDGET_BUTTON(overlay_btn_base, value='Load Overlay', event_pro='trial_overlay',$
;                                  uname='Load Overlay',xoffset=10)
;        rem_img = WIDGET_BUTTON(overlay_btn_base, value = 'Remove Overlay',event_pro='trial_overlay',$
;                                uname='Remove Overlay',xoffset = 175)
;        trans_slider = CW_FSLIDER(overlay_base, title='Transparency',/drag,/edit,minimum=0,$
;                                  value=1,maximum=1)

        
        reslice_lbl_base = WIDGET_BASE(inten_base,/column)
        reslice_lbl = WIDGET_LABEL(reslice_lbl_base,value='RESLICE OPTIONS')
        reslice_base = WIDGET_BASE(reslice_lbl_base,/column,/frame,/EXCLUSIVE)
;        reslice_base = widget_base(reslice_lbl_base,/column)
;        btn_values = ['Freq Phase','Freq Slice','Phase Slice']
;        reslice_btn_group = CW_BGROUP(reslice_base,btn_values,/exclusive,set_value=0)
        reslice_btn1 = WIDGET_BUTTON(reslice_base,value= 'Freq Phase',uname='Freq Phase')
        reslice_btn2 = WIDGET_BUTTON(reslice_base,value= 'Freq Slice',uname='Freq Slice')
        reslice_btn3 = WIDGET_BUTTON(reslice_base,value= 'Phase Slice',uname='Phase Slice')
        
        osv_base = WIDGET_BASE(inten_base,/row)
        osv_disp_btn = WIDGET_BUTTON(osv_base,value='Display in OrthoSliceViewer', uname='osv_disp_btn',$
                                /dynamic_resize, event_pro='trial_call_osv_event')
;        osv_upd_btn = WIDGET_BUTTON(osv_base, value='Update OrthoSliceViewer',uname='osv_upd_btn',$
;                                /dynamic_resize, event_pro='trial_call_osv_event')
        export_base = WIDGET_BASE(tab_widget,title='Export',/column)
        button_id = WIDGET_BUTTON(export_base, /dynamic_resize,yoffset = 10, $
                          uname='Save_button',event_pro='save_data_rat_event', value='Save Image')       
        
                                      
                      
; --------------------------------- REALISING THE WIDGETS-----------------------------------------
 
        widget_control , topbase, /realize
        IF (((*infoptr).AD_index eq 0) or ((*infoptr).S0_index eq 0)) THEN BEGIN
            WIDGET_CONTROL, bval_slider, sensitive=0
        ENDIF
        info_ptr = ptr_new({draw_window: draw_window, $
                 draw_id: draw_id, $
                 TRslider: TRslider, $
                 TEslider: TEslider, $
                 InvRslider : InvRslider, $
                 slice_slider: slice_slider,$
                 angle_slider: angle_slider, $
                 inten_label: inten_label, $
                 inten_text: inten_text,$
                 x_text:x_text, y_text:y_text, $
                 inten_slider: inten_slider, $
                 ortho_view_state: ptr_new(), $
                 seqname: seqname, $
                 pulseseq: pulseseq, $
                 button_id: button_id, $
                 bval_slider:bval_slider, $
                 sequence: sequence, $
                 reslice_btn1:reslice_btn1,$
                 reslice_btn2:reslice_btn2, $
                 reslice_btn3:reslice_btn3,$
                 fdim: fdim, pdim: pdim, sdim: sdim, $
                 T1:T1, T2:T2, Proton:Proton, AD:AD, $
                 fin_imagesc:fltarr(fdim,pdim,sdim), $
                 final_image:fltarr(fdim,pdim,sdim)})

        widget_control, topbase, set_uvalue=info_ptr
        xmanager, 'trial_qmri_tool',topbase, cleanup='trial_qmri_tool_cleanup', /no_block
end

