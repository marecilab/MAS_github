;+
; :Description:
;    This procedure was written to calculate the magnetic field Bz produced by a current, I, passing
;    through a wire.
;
;    It uses the phase images that are computed from the measurements.
;
;    NOTE:
;    Currently does not work with MultiEcho ICNE acquisitions and SingleEcho Single Slice acquisitions. 
;    
; :Author: Aditya Kumar
;-

PRO mask_gen, infoptr,output1=output1, output2=output2
    
    COMMON scan_data, project
    p_data_x = ptr_new(abs((*infoptr).image1))
    p_data_y = ptr_new(abs((*infoptr).image2))
    size_x_y = size(*p_data_x)
    mask_x = dblarr(size_x_y[1], size_x_y[2], size_x_y[3]) + 1.0
    mask_y = dblarr(size_x_y[1], size_x_y[2], size_x_y[3]) + 1.0
    max_x = max(*p_data_x)
    max_y = max(*p_data_y)
    threshold_x = where(*p_data_x lt 0.05*max_x)
    threshold_y = where(*p_data_y lt 0.05*max_y)
    mask_x[threshold_x] = 0
    mask_y[threshold_y] = 0
    output1 = mask_x
    output2 = mask_y
    (*infoptr).mask_flag=1
    
END

PRO mas_mreit_event, event
      
      WIDGET_CONTROL, event.top, get_uvalue=infoptr
      
      COMMON scan_data
      CI = project.ci
      
      rdim = project.imndarray[CI].fdim
      pdim = project.imndarray[CI].pdim
      sdim = project.imndarray[CI].sdim
      slicepos = project.procpramarray[CI].sdim_start
      adim = project.imndarray[CI].adim
      
      CASE event.id OF
      
      (*infoptr).Load_btn1: BEGIN
                      mas_load_state_1
                      test_data = *project.dataarray[CI].state1
                      IF ((size(test_data))[0] eq 4) THEN BEGIN
                          (*infoptr).image1 = test_data(*,*,*,(project.procpramarray[CI].adim_start))
                      ENDIF ELSE BEGIN
                          (*infoptr).image1 = *project.dataarray[CI].state1
                      END
                      (*infoptr).scan1 = CI
;                      help, (*infoptr).scan1, project.imndarray[CI].dimensions
                      scan_name_full = project.scan_list[CI]
                      scan_name_split = strsplit(scan_name_full,PATH_SEP(),/EXTRACT)
                      sz_split = size(scan_name_split)
                      scan_name = PATH_SEP() + scan_name_split[sz_split[1]-2] + PATH_SEP() + $
                                  scan_name_split[sz_split[1]-1] + PATH_SEP()
                      WIDGET_CONTROL, (*infoptr).Load_btn1_text, set_value=scan_name
                      help, scan_name, (*infoptr).scan1
                      END
                       
      (*infoptr).Load_btn2: BEGIN
                      mas_load_state_1
                      test_data = *project.dataarray[CI].state1
                      IF ((size(test_data))[0] eq 4) THEN BEGIN
                          (*infoptr).image2 = test_data(*,*,*,(project.procpramarray[CI].adim_start))
                      ENDIF ELSE BEGIN
                          (*infoptr).image2 = *project.dataarray[CI].state1
                      END
                      (*infoptr).scan2 = CI
;                      help, (*infoptr).scan2, project.imndarray[CI].dimensions
                      scan_name_full = project.scan_list[CI]
                      scan_name_split = strsplit(scan_name_full,PATH_SEP(),/EXTRACT)
                      sz_split = size(scan_name_split)
                      scan_name = PATH_SEP() + scan_name_split[sz_split[1]-2] + PATH_SEP() + $
                                  scan_name_split[sz_split[1]-1] + PATH_SEP()
                      WIDGET_CONTROL, (*infoptr).Load_btn2_text, set_value=scan_name
                      help, scan_name, (*infoptr).scan2
                      END
                      
      (*infoptr).Load_btn3: BEGIN
                      mas_load_state_1
                      test_data = *project.dataarray[CI].state1
                      IF ((size(test_data))[0] eq 4) THEN BEGIN
                          (*infoptr).image3 = test_data(*,*,*,(project.procpramarray[CI].adim_start))
                      ENDIF ELSE BEGIN
                          (*infoptr).image3 = *project.dataarray[CI].state1
                      END
                      (*infoptr).scan3 = CI
;                      help, (*infoptr).scan3, project.imndarray[CI].dimensions
                      scan_name_full = project.scan_list[CI]
                      scan_name_split = strsplit(scan_name_full,PATH_SEP(),/EXTRACT)
                      sz_split = size(scan_name_split)
                      scan_name = PATH_SEP() + scan_name_split[sz_split[1]-2] + PATH_SEP() + $
                                  scan_name_split[sz_split[1]-1] + PATH_SEP()
                      WIDGET_CONTROL, (*infoptr).Load_btn3_text, set_value=scan_name
                      WIDGET_CONTROL, (*infoptr).Ic_text, editable=1
                      help, scan_name, (*infoptr).scan3
                      END
                      
      (*infoptr).Ic_text:   BEGIN
                      WIDGET_CONTROL, (*infoptr).Ic_text, get_value = Ic_value
                      (*infoptr).Ic_value = FLOAT(Ic_value)
                      help, (*infoptr).Ic_value
                      WIDGET_CONTROL, (*infoptr).Tc_text, editable=1      
                      END
      
      (*infoptr).Tc_text:   BEGIN
                      WIDGET_CONTROL, (*infoptr).Tc_text, get_value = Tc_value
                      (*infoptr).Tc_value = FLOAT(Tc_value)
                      help, (*infoptr).Tc_value
                      WIDGET_CONTROL, (*infoptr).icne_btn, sensitive = 1
                      WIDGET_CONTROL, (*infoptr).sems_btn, sensitive = 1
                      WIDGET_CONTROL, (*infoptr).mems_btn, sensitive = 1                    
                      WIDGET_CONTROL, (*infoptr).CD_btn, sensitive = 1
                      END
                      
      (*infoptr).sems_btn: BEGIN
                      sems_flag = WIDGET_INFO((*infoptr).sems_btn, /BUTTON_SET)
                      (*infoptr).sems_flag = sems_flag
                      IF ((*infoptr).sems_flag eq 1) THEN BEGIN
                          WIDGET_CONTROL, (*infoptr).mems_btn, sensitive = 0
                      ENDIF ELSE BEGIN
                          WIDGET_CONTROL,(*infoptr).mems_btn, sensitive = 1
                      END
                      END
                      
      (*infoptr).mems_btn: BEGIN
                      mems_flag = WIDGET_INFO((*infoptr).mems_btn, /BUTTON_SET)
                      (*infoptr).mems_flag = mems_flag
                      IF ((*infoptr).mems_flag eq 1) THEN BEGIN
                          WIDGET_CONTROL, (*infoptr).sems_btn, sensitive = 0
                      ENDIF ELSE BEGIN
                          WIDGET_CONTROL, (*infoptr).sems_btn, sensitive = 1
                      END
                      END
                      
      (*infoptr).icne_btn:  BEGIN
                      icne_flag = WIDGET_INFO((*infoptr).icne_btn, /BUTTON_SET)
                      (*infoptr).icne_flag = icne_flag
                      help, (*infoptr).icne_flag
;                      IF((*infoptr).icne_flag eq 1) THEN BEGIN
;                           project.procpramarray[(*infoptr).scan1].signal_type = 0
                          mesg = dialog_message('Please change the no current dataset image type to MAGNITUDE.', /information, /center)
;                      ENDIF ELSE BEGIN
;                           project.procrpramarray[(*infoptr).scan1].signal_type = 9
;                      END 
                      END
                     
      (*infoptr).CD_btn:    BEGIN
                       IF ((*infoptr).mems_flag eq 0 ) THEN BEGIN                     ;;; If Single echo acquisition
                          IF((*infoptr).icne_flag eq 0) THEN BEGIN
                            (*infoptr).cd_image = ((*infoptr).image2)/((*infoptr).image3)
                            WIDGET_CONTROl, (*infoptr).Phase_btn, sensitive=1
                            WIDGET_CONTROL, (*infoptr).mask_btn, sensitive=1
                            WIDGET_CONTROL, (*infoptr).phaseunwrap_btn, sensitive=1
                          ENDIF ELSE BEGIN                                            ;;; If Single echo ICNE
                            mas_mreit_icne, infoptr, output1 = image          
                            (*infoptr).Bc = image
                            WIDGET_CONTROL, (*infoptr).Field_btn, sensitive=1
                            WIDGET_CONTROL, (*infoptr).mask_btn, sensitive=1
                            WIDGET_CONTROL, (*infoptr).phaseunwrap_btn, sensitive=1
                          END
                        ENDIF ELSE BEGIN                                              ;;; If MultiEcho acquisition
                          mas_mreit_mems, infoptr, output1 = image
                          (*infoptr).phaseimage = image
                          WIDGET_CONTROl, (*infoptr).Phase_btn, sensitive=1
                          WIDGET_CONTROL, (*infoptr).mask_btn, sensitive=1
                          WIDGET_CONTROL, (*infoptr).phaseunwrap_btn, sensitive=1
                      END
                    END
                          
      (*infoptr).mask_btn: BEGIN
                      mask_gen, infoptr, output1 = mask_1, output2 = mask_2
                      (*infoptr).mask_image = (*infoptr).mask_image*mask_1*mask_2
                      END
                      
      (*infoptr).Phase_btn: BEGIN
                      IF ((*infoptr).mems_flag eq 0)THEN BEGIN
                          (*infoptr).phaseimage = atan((*infoptr).cd_image, /phase)
                      END
                      IF((size((*infoptr).phaseimage))[0] eq 2) THEN BEGIN
                        phase2dimage = (*infoptr).phaseimage[*,*]
;                        iimage, phase2dimage
                      ENDIF ELSE BEGIN
                        phase_data_ptr = ptr_new((*infoptr).phaseimage)
                        mas_display_ortho, data_ptr = phase_data_ptr
                      END
                      WIDGET_CONTROL, (*infoptr).Field_btn, sensitive=1
                      END
      
      (*infoptr).Field_btn: BEGIN
                      IF((*infoptr).icne_flag eq 0) THEN BEGIN
                          (*infoptr).Bc = ((*infoptr).phaseimage)/(2*((*infoptr).Gamma_value)*((*infoptr).Tc_value))
                      ENDIF ELSE BEGIN
;                          (*infoptr).Bc = (*infoptr).cd_image   ; IF ICNE, cd_image is field image.
                          IF((*infoptr).mask_flag eq 1) THEN BEGIN
                                (*infoptr).Bc = (*infoptr).Bc*(*infoptr).mask_image
                          END
                      END               
                      Bc_data_ptr = ptr_new((*infoptr).Bc)
                      IF((size((*infoptr).Bc))[0] eq 2) THEN BEGIN
;                        iimage, (*infoptr).Bc
                      ENDIF ELSE BEGIN
                        mas_display_ortho, data_ptr = Bc_data_ptr
;                        iimage, (*infoptr).Bc[*,*,slicepos]
                      END
                      WIDGET_CONTROL, (*infoptr).export_MRBz_btn, sensitive=1
                      WIDGET_CONTROL, (*infoptr).export_nifti_btn, sensitive=1
                      END
                      
      (*infoptr).export_MRBz_btn: BEGIN 
                      data_ptr = ptr_new((*infoptr).Bc)
                      mag_data_ptr = ptr_new(abs((*infoptr).image1))
                      mas_mreit_save_to_coreha, data_ptr, mag_data_ptr
                      END    
      
      (*infoptr).export_nifti_btn: BEGIN
                      field_data_ptr = ptr_new((*infoptr).Bc)
                      mag_data_ptr = ptr_new(abs((*infoptr).image1))
                      dirname = dialog_pickfile(/directory)
                      mag_fname = dirname+'MR.nii.gz'
                      mas_export_nifti, data_ptr = mag_data_ptr, file_name = mag_fname
                      bz_fname = dirname+'Bz.nii.gz'
                      mas_export_nifti, data_ptr = field_data_ptr, file_name = bz_fname
                      mesg = dialog_message(['The images have been saved at:',dirname],/information, /center)
                      END
      (*infoptr).raw_data_btn: BEGIN
                      raw_data_flag = WIDGET_INFO((*infoptr).icne_btn, /BUTTON_SET)
                      (*infoptr).raw_data_flag = raw_data_flag
                      END
      ENDCASE                

END

PRO mas_mreit_cleanup, id

  print, 'Cleaning Up'
  widget_control, id, get_uvalue=infoptr
  ptr_free, infoptr
  print, 'Done'

END

PRO mas_mreit

      COMMON scan_data, project
      CI = project.ci
      rdim = project.imndarray[CI].fdim
      pdim = project.imndarray[CI].pdim
      sdim = project.imndarray[CI].sdim
      adim = project.imndarray[CI].adim
      
      bw = project.imndarray[CI].bandwidth
      slicepos = project.procpramarray[CI].sdim_start
      
      Bc = FLTARR(rdim,pdim,sdim)
      phaseimage = FLTARR(rdim,pdim,sdim)
      mask_image = FLTARR(rdim,pdim,sdim)+1.0
      cd_image = COMPLEX(FLTARR(rdim,pdim,sdim))
      
      Gamma_value = 267.5*1e6                    ; Gyromagnetic Ratio, (rad s-1 T-1)

;
;      
;      ********************* CREATING THE GUI *********************
;      
;

     topbase = WIDGET_BASE(TITLE='MREIT_TOOLBOX',/column)
     label_base = WIDGET_BASE(topbase,/column)
     label = WIDGET_LABEL(label_base, value='PLEASE CHOOSE THE DATASETS TO BEGIN',$
                          /align_center)
     Btn_base = WIDGET_BASE(label_base,/column)
     
     Image1_base = WIDGET_BASE(Btn_base,/row)
     Load_btn1 = WIDGET_BUTTON(Image1_base,value="Load 1st image",uname="Load Image1")
     Load_btn1_label = WIDGET_LABEL(Image1_base,value="No Current injection:")
     Load_btn1_text = WIDGET_TEXT(Image1_base, value = 'Empty', xsize = 74, editable = 0)
     
     Image2_base = WIDGET_BASE(Btn_base,/row)
     Load_btn2 = WIDGET_BUTTON(Image2_base,value="Load 2nd image",uname="Load Image2")
     Load_btn2_label = WIDGET_LABEL(Image2_base,value="Positive Current injection:")
     Load_btn2_text = WIDGET_TEXT(Image2_base, value = 'Empty', xsize = 49, editable = 0)
     
     Image3_base = WIDGET_BASE(Btn_base,/row)
     Load_btn3 = WIDGET_BUTTON(Image3_base,value="Load 3rd image",uname="Load Image1")
     Load_btn3_label = WIDGET_LABEL(Image3_base,value="Negative Current Injection:")
     Load_btn3_text = WIDGET_TEXT(Image3_base, value = 'Empty', xsize = 49, editable = 0)
     
     Props_base = WIDGET_BASE(Btn_base,/row)
     Ic_base = WIDGET_BASE(Props_base,/row)
     Ic_label = WIDGET_LABEL(Ic_base,value='Injected Current (mA):')
     Ic_text = WIDGET_TEXT(Ic_base,value='0',editable=0)     
     Tc_base = WIDGET_BASE(Props_base,/row)
     Tc_label = WIDGET_LABEL(Tc_base,value = 'Current Duration Tc, (s)):')
     Tc_text = WIDGET_TEXT(Tc_base,value='0',editable=0)
     
     Magnet_type_base = WIDGET_BASE(Btn_base,/row) ;;;; Added by AK2 for the 3T data processing ability. WIP.
     Philips_base = WIDGET_BASE(Magnet_type_base,/nonexclusive)
     Philips_flag = WIDGET_BUTTON(Philips_base,value='HumanMagnet',uname='check_if_humanmagnet',sensitive=0)
     Agilent_base = WIDGET_BASE(Magnet_type_base,/nonexclusive)
     Agilent_flag = WIDGET_BUTTON(Agilent_base,value='HorizontalMagnet',uname='check_if_horizontalmagnet',sensitive=0)
     Bruker_base = WIDGET_BASE(Magnet_type_base,/nonexclusive)
     Bruker_flag = WIDGET_BUTTON(Bruker_base, value='VerticalMagnet',uname='check_if_verticalmagnet',sensitive=0)
     
     Datatype_base = WIDGET_BASE(Btn_base, /row)      ;;;; Added by AK2 for the 3T data processing ability. WIP.
     Rawdata_base = WIDGET_BASE(Magnet_type_base,/nonexclusive)
     Rawdata_flag = WIDGET_BUTTON(rawdata_base, value='',uname='',sensitive=0)
     Recondata_base = WIDGET_BASE(Magnet_type_base,/nonexclusive)
     Recondata_flag = WIDGET_BUTTON(Recondata_base, value='',uname='',sensitive=0)
     
     Acqstyle_base = WIDGET_BASE(Btn_base,/row)
     sems_base = WIDGET_BASE(Acqstyle_base,/nonexclusive)
     sems_btn = WIDGET_BUTTON(sems_base,value='SingleEchoMultislice', uname='check_if_sems',sensitive=0)
     mems_base = WIDGET_BASE(Acqstyle_base,/nonexclusive)
     mems_btn = WIDGET_BUTTON(mems_base,value='MultiEchoMultiSlice', uname='check_if_mems',sensitive=0)
     icne_base = WIDGET_BASE(Acqstyle_base,/nonexclusive)
     icne_btn = WIDGET_BUTTON(icne_base,value='ICNE', uname='check_if_ICNE',sensitive=0)
         
     operations_base = WIDGET_BASE(Btn_base,/row)
     CD_base = WIDGET_BASE(operations_base,/column)
     CD_btn = WIDGET_BUTTON(CD_base,value="Complex Divide",uname="CD_button",sensitive=0)
     
     mask_base = WIDGET_BASE(operations_base,/nonexclusive)
     mask_btn = WIDGET_BUTTON(mask_base,value='Mask', uname='Calc_MaskImage',sensitive=0)
     
     phaseunwrap_base = WIDGET_BASE(operations_base,/nonexclusive)
     phaseunwrap_btn = WIDGET_BUTTON(phaseunwrap_base,value='Phase Unwrap', uname='UnwrapPhase',sensitive=0)
 
     Phase_base = WIDGET_BASE(operations_base,/column)
     Phase_btn = WIDGET_BUTTON(Phase_base,value="Phase Image",uname="Disp_PhaseImage",sensitive=0)
     
     Field_base = WIDGET_BASE(operations_base,/column)
     Field_btn = WIDGET_BUTTON(Field_base,value="Magnetic Field Image",uname="Disp_FieldImage",sensitive=0)
     
     export_MRBz_base = WIDGET_BASE(Btn_base,/column)
     export_MRBz_btn = WIDGET_BUTTON(export_MRBz_base,value='Save MR & Bz images for CoReHa',uname = 'MR_Bz_for_CoReHA',sensitive=0)
     
     export_nifti_base = WIDGET_BASE(Btn_base,/column)
     export_nifti_btn = WIDGET_BUTTON(export_nifti_base,value='Export MR & Bz images in Nifti',uname = 'MR_Bz_for_Nifti',sensitive=0)
     
     
     WIDGET_CONTROL, topbase, /REALIZE
     
;     INFOPTR STRUCTURE - variables description --------------------------------------------------------     
;     scan1 , scan2, scan3    =  flags for the no current and injected current datasets.
;     image1, image2, image3  = complex arrays to store the corresponding image datasets.
;     cd_image                = complex divided image dataset.
;     Gamma_value             = gyromagnetic ratio
;     Bc                      = induced magnetic field image.
;     phaseimage              = induced phase image.
;     mask_image              = mask image for filtering the noise.
;     Ic_value                = injected current value.
;     Tc_value                = total injected current duration per TR.
;     icne_flag               = flag to check whether icne type of acquisition.
;     sems/mems_flag          = flags to check whether it is a single echo or a multi echo acquisition.
;     mask_flag               = flag to check if mask exists. Prevents recomputation of mask.
;     phaseunwrap             = flag to check for phase unwrapping.
;     threshold               = harcoded to be 5%. Can be changed in the code. Will be allocated in the GUI soon. 
;     rdim,pdim,sdim,adim     = image dataset dimensions consistent with selected MAS scan.      
;     The remaining variables are the appropriate widget ids and handlers.
     
     infoptr = PTR_NEW({scan1:-1,scan2:-1,scan3:-1,image1:COMPLEX(FLTARR(rdim,pdim,sdim)),image2:COMPLEX(FLTARR(rdim,pdim,sdim)),cd_image:cd_image,$
                        image3:COMPLEX(FLTARR(rdim,pdim,sdim)),Gamma_value:Gamma_value,Bc:Bc,phaseimage:phaseimage,mask_image:mask_image, $
                        Load_btn1:Load_btn1,Load_btn1_text:Load_btn1_text, Load_btn3:Load_btn3, Ic_text:Ic_text, Tc_text:Tc_text,$
                        Load_btn2:Load_btn2,Load_btn2_text:Load_btn2_text, Load_btn3_text:Load_btn3_text,bw:bw, Ic_value:0.0, Tc_value:0.0,$
                        CD_btn:CD_btn, Phase_btn:Phase_btn, Field_btn:Field_btn, mask_btn:mask_btn, export_nifti_btn:export_nifti_btn, $
                        export_MRBz_btn:export_MRBz_btn, phaseunwrap_btn:phaseunwrap_btn,icne_btn:icne_btn,sems_btn:sems_btn,mems_btn:mems_btn,$
                        rdim:rdim,pdim:pdim,sdim:sdim,adim:adim,threshold:0.05,mask_flag:0,phaseunwrap:0,icne_flag:0,sems_flag:0,mems_flag:0})
     
     WIDGET_CONTROL, topbase, set_uvalue=infoptr, /no_copy
     XMANAGER, 'mas_mreit', topbase, /no_block
     
END