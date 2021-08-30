
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


PRO trial_make_proton_map_event, event

  WIDGET_CONTROL, event.top, get_uvalue = infoptr
  WIDGET_CONTROL, (*infoptr).tr_text, get_value = tr_value
        
        tr_string = strtrim(tr_value)
        reptime = double(tr_string)
        T1 = (*infoptr).T1
        Mo = (*infoptr).Mo
        M = (*infoptr).M
        fac = (*infoptr).fac
        factor = (*infoptr).factor
        rho = (*infoptr).rho
        rhoclean = (*infoptr).rhoclean
        fdim = (*infoptr).fdim
        pdim=(*infoptr).pdim
        sdim =(*infoptr).sdim
        TR = (*infoptr).TR
        help, TR, reptime
        
        M = clean_image(Mo)
        TR = make_array(fdim, pdim, sdim, /double, value=reptime)
        help, M , TR
        
        fac = (1 - exp( -TR/T1))
        factor = clean_image(fac)
        help , fac, factor
        
        rho = M/factor
        rhoclean = clean_image(rho)
        rhoptr = ptr_new(rhoclean)
        help, rho, rhoclean
        
        mas_export_nifti, data_ptr = rhoptr, file_name = 'ProtonDensity_map.nii.gz'
        cd, current=thisDir
        res = DIALOG_MESSAGE(['The image has been saved at',thisDir])

END
;
;
;---------------------------------CREATE PROTON DENSITY MAP FROM THE GIVEN T1 & T2 DATASETS--------------------------
;
;

;--- T1 dataset : T1 map generated from a variable TR or a variable TI dataset.
;--- T2 dataset : T2 map generated from a variable TE dataset.
;--- Mo dataset : Mo map generated from the T2 fit used to obtain the T2 map.

PRO trial_make_proton_map 
  
  COMMON scan_data, project
  
    T1 = *project.dataarray[0].state1 
    T2 = *project.dataarray[1].state1
    Mo = *project.dataarray[2].state1
    fdim = project.imndarray[0].fdim
    pdim = project.imndarray[0].pdim
    sdim = project.imndarray[0].sdim
    
    TR = dblarr(fdim,pdim,sdim)
    M = dblarr(fdim,pdim,sdim)
    fac = dblarr(fdim,pdim,sdim)
    PD = dblarr(fdim,pdim,sdim)
    factor = dblarr(fdim,pdim,sdim)
    rho = dblarr(fdim,pdim,sdim)
    rhoclean = dblarr(fdim,pdim,sdim)
     
        ; CREATING GUI
        
        topbase = WIDGET_BASE(title='Compute Proton Density',/row)
        tr_label = WIDGET_LABEL(topbase,value='Enter TR value (secs):')
        tr_text = WIDGET_TEXT(topbase,/editable)
        
        WIDGET_CONTROL, topbase, /realize
        
        infoptr = PTR_NEW({T1:T1,T2:T2,Mo:Mo,$
                           fdim:fdim, pdim:pdim, sdim:sdim, $
                           PD:PD, fac:fac, factor:factor, rho:rho, rhoclean:rhoclean, M:M,$
                           TR:TR,$
                           tr_label:tr_label, tr_text:tr_text})
                           
        WIDGET_CONTROL, topbase, set_uvalue = infoptr
        XMANAGER, 'trial_make_proton_map', topbase, /no_block
        
END
