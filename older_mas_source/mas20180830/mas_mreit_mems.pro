;+
;
; This code is written to perform the MREIT data processing for multi echo acquisitions. It is called from the program "mas_mreit.pro."
; It stores the ids of the datasets from the main MAS scanlist window and accesses the multiple echoes through the adim parameter.
; 
; OUTPUT: 
; It returns the average phase image from the acquired datasets to generate the reduced noise phase/field map.
; 
; NOTE: 
; The number of echoes in the positive and negative current injection SHOULD BE the same.
; 
;  
; :Author: aditya
;-


PRO mas_mreit_mems, infoptr, output1=image

  COMMON scan_data
  forward_function get_dir_slash
  dir_sep = (get_dir_slash())[1]
  
  scan_list_array = [(*infoptr).scan1, (*infoptr).scan2, (*infoptr).scan3]
;  dataptr = [ptr_new(), ptr_new(), ptr_new()]
  rdim = project.imndarray[scan_list_array(1)].fdim
  pdim = project.imndarray[scan_list_array(1)].pdim
  sdim = project.imndarray[scan_list_array(1)].sdim
  adim = project.imndarray[scan_list_array(1)].adim
  slicepos = project.procpramarray[scan_list_array(1)].sdim_start
  
  image1 = COMPLEX(FLTARR(rdim,pdim,sdim))
  phase = FLTARR(rdim,pdim,sdim)
  phase_add = FLTARR(rdim,pdim,sdim)
  
  data_no_current = *project.dataarray[scan_list_array[0]].state1       ; No current arrayed data.
  data_pos = *project.dataarray[scan_list_array(1)].state1              ; Positive current arrayed data.
  data_neg = *project.dataarray[scan_list_array(2)].state1              ; Negative current arrayed data.
  
  IF((*infoptr).icne_flag eq 0) THEN BEGIN
    FOR j = 0, adim-1 DO BEGIN
        FOR i = 0, sdim-1 DO BEGIN
          image1[*,*,i] = data_pos[*,*,i,j]/data_neg[*,*,i,j]
        END
        phase = atan(image1,/phase)
;        iimage, phase[*,*,slicepos]                       ; test look at phase data from each echo.
        IF ((size(phase))[0] eq 2) THEN BEGIN
              iimage, phase                              ; For single slice acquisitions
        END
        IF (j eq 1) THEN BEGIN
            phase_add = phase_add - phase
            iimage, phase_add[*,*,slicepos]               ; cumulative phase until the jth echo
        ENDIF ELSE BEGIN
            phase_add = phase_add + phase
            iimage, phase_add[*,*,slicepos]               ; cumulative phase until the jth echo
        END
    END 
    image = phase_add/adim
    iimage, image[*,*,slicepos]                           ; average phase image.
  ENDIF ELSE BEGIN                                                     ; For multi echo ICNE Acquisition.
      mesg = DIALOG_MESSAGE('MultiEcho ICNE still not available in this program')
  END
  
END