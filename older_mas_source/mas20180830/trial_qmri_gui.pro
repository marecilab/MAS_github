

FUNCTION prep_wid_list, scanlist
    COMMON scan_data
    NI = project.ni
    maps_list = STRARR(NI)
    fnames = STRARR(NI)
    FOR i = 0,NI-1 DO BEGIN
        map = scanlist[i]
        map_split = STRSPLIT(map, PATH_SEP(),/EXTRACT)
        sz = SIZE(map_split)
        map_temp = PATH_SEP() + map_split[sz[1]-2] + PATH_SEP() + map_split[sz[1]-1]
        maps_list[i] = map_temp
        fnames[i] = map_split[sz[1]-1]
    ENDFOR
    map_list = ['................................................................', $
                maps_list]
    list_struct = {map_list:map_list,fnames:fnames}
    RETURN, list_struct
END

PRO trial_call_qmri_tool, event
  
        COMMON scan_data, project
        WIDGET_CONTROL, event.top, get_uvalue = infoptr
        
        flag = (*infoptr).flag_array
        T1_index = (*infoptr).T1_index
        T2_index = (*infoptr).T2_index
        PD_index = (*infoptr).PD_index
        AD_index = (*infoptr).AD_index
        S0_index = (*infoptr).S0_index
        
          IF(T1_index gt 0) THEN BEGIN 
          flag[0] = 1
          ENDIF
          IF(T2_index gt 0) THEN BEGIN
          flag[1] = 1
          ENDIF
          IF(PD_index gt 0) THEN BEGIN 
          flag[2] = 1
          ENDIF 
          IF(AD_index gt 0) THEN BEGIN 
          flag[3] = 1
          ENDIF
          IF(S0_index gt 0) THEN BEGIN 
          flag[4] = 1
          ENDIF
          
          IF ((flag[0] eq 0) or (flag[1] eq 0)) THEN BEGIN
                res = DIALOG_MESSAGE(['Please select at least the T1, T2 and PD maps.','The tool cannot be opened without them.'])
                RETURN
          ENDIF
          
          IF (flag[2] eq 0) THEN BEGIN
                res = DIALOG_MESSAGE(['Please select at least the T1, T2 and PD maps.','The tool cannot be opened without them.'])
                RETURN
          ENDIF
          
          IF ((flag[3] eq 0) or (flag[4]) eq 0) THEN BEGIN
              T1_index = (*infoptr).T1_index
              T2_index = (*infoptr).T2_index
              PD_index = (*infoptr).PD_index
              
              T1 = *project.dataarray[T1_index-1].state1
              T2 = *project.dataarray[T2_index-1].state1
              Proton = *project.dataarray[PD_index-1].state1
              
              fdim = project.imndarray[T1_index-1].fdim
              pdim = project.imndarray[T1_index-1].pdim
              sdim = project.imndarray[T1_index-1].sdim     
              AD = MAKE_ARRAY(fdim,pdim,sdim)
              S0 = MAKE_ARRAY(fdim,pdim,sdim)
              
              trial_qmri_tool, T1,T2,Proton,AD,S0,fdim,pdim,sdim,infoptr
          ENDIF ELSE BEGIN
        
              T1_index = (*infoptr).T1_index
              T2_index = (*infoptr).T2_index
              PD_index = (*infoptr).PD_index
              AD_index = (*infoptr).AD_index
              S0_index = (*infoptr).S0_index
         
              T1 = *project.dataarray[T1_index-1].state1
              T2 = *project.dataarray[T2_index-1].state1
              Proton = *project.dataarray[PD_index-1].state1
              AD = *project.dataarray[AD_index-1].state1
              S0 = *project.dataarray[S0_index-1].state1
             
              fdim = project.imndarray[T1_index-1].fdim
              pdim = project.imndarray[T1_index-1].pdim
              sdim = project.imndarray[T1_index-1].sdim
        
              trial_qmri_tool, T1,T2,Proton,AD,S0,fdim,pdim,sdim,infoptr
          ENDELSE

END


PRO trial_qmri_gui_event, event

      COMMON scan_data
      WIDGET_CONTROL, event.top, get_uvalue=infoptr 
      flag_array = (*infoptr).flag_array
      map_list = (*infoptr).map_list
;      IF (event.id eq (*infoptr).upd_btn) THEN BEGIN
;            new_list = project.scan_list
;            req_maplist = ['T1_map.nii.gz','T2_map.nii.gz','ProtonDensity_map.nii.gz','AD_map.nii.gz','S0_map.nii.gz']
;            scanlist = (*infoptr).scanlist
;            sz_scanlist = (*infoptr).sz
;            PRINT, scanlist
;            PRINT, new_list
;            PRINT, sz_scanlist
;            FOR i = 0, sz_scanlist[1]-1 DO BEGIN
;                IF (scanlist[i] eq '') THEN BEGIN
;                    CONTINUE
;                ENDIF
;                fname = get_filename(scanlist[i])
;                PRINT, fname
;                test_res = WHERE(STRMATCH(req_maplist,fname))
;                PRINT, 'test_res =', test_res
;                new_index = WHERE(STRMATCH(new_list,scanlist[i]))
;                PRINT, 'new_index',new_index
;                    IF (new_index eq -1) THEN BEGIN
;                          IF(test_res ne -1) THEN BEGIN
;                              CASE test_res OF
;                                  0:  BEGIN
;                                      note=DIALOG_MESSAGE(['Required data deleted. Please reload the nifti maps from the folder.',$
;                                      'Please ensure atleast the T1,T2 and PD maps exist.'])
;                                      RETURN
;                                      END
;                                  1:  BEGIN
;                                      note=DIALOG_MESSAGE(['Required data deleted. Please reload the nifti maps from the folder.',$
;                                      'Please ensure atleast the T1,T2 and PD maps exist.'])
;                                      RETURN
;                                      END
;                                  2:  BEGIN
;                                      note=DIALOG_MESSAGE(['Required data deleted. Please reload the nifti maps from the folder.',$
;                                      'Please ensure atleast the T1,T2 and PD maps exist.'])
;                                      RETURN
;                                      END
;                                  3:  BEGIN
;                                      flag_array[test_res] = 0
;                                      END
;                                  4:  BEGIN
;                                      flag_array[test_res] = 0
;                                      END
;
;                                ENDCASE
;                           ENDIF ELSE BEGIN
;                           PRINT, 'File Removed.'
;                           ENDELSE
;                    ENDIF ELSE BEGIN
;                          IF(test_res ne -1) THEN BEGIN
;                              CASE test_res OF
;                                  0:  BEGIN
;                                      (*infoptr).T1_index = new_index + 1
;;                                      flag_array[0] = 1
;                                      END
;                                  1:  BEGIN
;                                      (*infoptr).T2_index = new_index + 1
;;                                      flag_array[1] = 1
;                                      END
;                                  2:  BEGIN
;                                      (*infoptr).PD_index = new_index + 1
;;                                      flag_array[2] = 1
;                                      END
;                                  3:  BEGIN
;                                      (*infoptr).AD_index = new_index + 1
;;                                      flag_array[3] = 1
;                                      END
;                                  4:  BEGIN
;                                      (*infoptr).S0_index = new_index + 1
;;                                      flag_array[4] = 1
;                                      END
;                              ENDCASE
;                           ENDIF ELSE BEGIN
;                           PRINT, 'Complete'
;                           ENDELSE
;                    ENDELSE
;            ENDFOR
;            PRINT, (*infoptr).T1_index, (*infoptr).T2_index , (*infoptr).PD_index ,(*infoptr).AD_index ,(*infoptr).S0_index 
;      ENDIF ELSE BEGIN
;      PRINT, (*infoptr).T1_index, (*infoptr).T2_index , (*infoptr).PD_index ,(*infoptr).AD_index ,(*infoptr).S0_index        
      
      CASE event.id OF
      
        (*infoptr).list1: BEGIN
                          (*infoptr).T1_index = event.index
                          IF(event.index ne 0) THEN flag_array[0] = 1
                          PRINT, (*infoptr).T1_index, flag_array
                          END
                          
        (*infoptr).list2: BEGIN
                          (*infoptr).T2_index = event.index
                          IF(event.index ne 0) THEN flag_array[1] = 1
                          PRINT, (*infoptr).T2_index, flag_array
                          END
                          
        (*infoptr).list3: BEGIN
                          (*infoptr).PD_index = event.index
                          IF(event.index ne 0) THEN flag_array[2] = 1
                          PRINT, (*infoptr).PD_index, flag_array
                          END
                          
        (*infoptr).list4: BEGIN
                          (*infoptr).AD_index = event.index
                          IF(event.index ne 0) THEN flag_array[3] = 1
                          PRINT, (*infoptr).AD_index, flag_array
                          END
        
        (*infoptr).list5: BEGIN
                          (*infoptr).S0_index = event.index
                          IF(event.index ne 0) THEN flag_array[4] = 1
                          PRINT, (*infoptr).S0_index, flag_array
                          END

      ENDCASE
;      ENDELSE
END

PRO trial_qmri_gui

    COMMON scan_data
    NI = project.ni
    CI = project.ci
   
    scanlist = project.scan_list
    sz = SIZE(scanlist)
    flag_array = MAKE_ARRAY(5)
;    FOR i =0,4 DO BEGIN
;        flag_array[i] = 1
;    ENDFOR
    
    list_struct = prep_wid_list(scanlist)
    map_list = list_struct.map_list
    fname_list = list_struct.fnames
                
    PRINT, map_list, fname_list                
    topbase = WIDGET_BASE(title='Select the appropriate maps',/column)
    label_base = WIDGET_BASE(topbase,/base_align_center)
    label = WIDGET_LABEL(label_base, value='PLEASE CHOOSE THE DATASETS TO BEGIN',$
                          /align_center)
    list_base = WIDGET_BASE(topbase,/column)
    
    list1_base = WIDGET_BASE(list_base,/column)
    list1 = WIDGET_DROPLIST(list1_base,title = 'Choose T1 map', $
                            value=map_list,/dynamic_resize)
    list2_base = WIDGET_BASE(list_base,/column)
    list2 = WIDGET_DROPLIST(list2_base,title = 'Choose T2 map', $
                            value=map_list,/dynamic_resize)
    list3_base = WIDGET_BASE(list_base,/column)
    list3 = WIDGET_DROPLIST(list3_base,title = 'Choose PD map', $
                            value=map_list,/dynamic_resize)
    list4_base = WIDGET_BASE(list_base,/column)
    list4 = WIDGET_DROPLIST(list4_base,title = 'Choose AD map', $
                            value=map_list,/dynamic_resize)
    list5_base = WIDGET_BASE(list_base,/column)
    list5 = WIDGET_DROPLIST(list5_base,title = 'Choose S0 map', $
                            value=map_list,/dynamic_resize)
    
    btn_base = WIDGET_BASE(topbase, /row)
    updt_base = WIDGET_BASE(btn_base)
;    upd_btn = WIDGET_BUTTON(updt_base,value='Update List',uname='UpdateButton',$
;                            xoffset=15, xsize=100, event_pro='trial_update_scanlist.pro')
;    upd_btn = WIDGET_BUTTON(updt_base,value='Update List',uname='UpdateButton',$
;                            xoffset=15, xsize=100)
    gen_base = WIDGET_BASE(btn_base,/base_align_center)
;    gen_btn = WIDGET_BUTTON(gen_base,value='Open Contrast Manipulator',$
;                            uname='OpenConMan',xoffset=50, event_pro='trial_opencm.pro')
    gen_btn = WIDGET_BUTTON(gen_base,value='Generate Tool',$
                            uname='OpenConMan',xoffset=125, event_pro='trial_call_qmri_tool')
    WIDGET_CONTROL, topbase, /REALIZE
    

    infoptr = PTR_NEW({list1:list1,list2:list2,list3:list3,list4:list4,list5:list5,$
                       T1_index:0,T2_index:0,PD_index:0,AD_index:0,S0_index:0,$
;                       upd_btn:upd_btn,gen_btn:gen_btn,$
                       gen_btn:gen_btn, $
                       scanlist:scanlist, map_list:map_list,sz:sz,$
                       flag_array:flag_array})

    
    WIDGET_CONTROL, topbase, set_uvalue=infoptr 
    XMANAGER, 'trial_qmri_gui', topbase, /no_block
 
END