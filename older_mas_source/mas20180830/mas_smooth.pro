;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved

; Subroutine name: rflt
; Created by: Originally written by E. Bossart.
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; Read in the *.flt data file, file_name, with the header then return the
; data to the calling program without the header

; Editing Information:
; Modified by T. Mareci, 17 August 2000
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

Function rflt, file_name

; Set up the array to determine the dimensionality of the data (up to 5)
d_num = BytArr(4)

; Read in the byte number defining the dimensionality of the data
OpenR, lun, file_name, /Get_LUN, /SWAP_IF_LITTLE_ENDIAN
ReadU, lun, d_num

; Determine the dimensionality
d = 256^3*d_num(0) + 256^2*d_num(1) + 256*d_num(2) + d_num(3)

If d EQ 2 Then hdr_sz = 16
If d EQ 3 Then hdr_sz = 20
If d EQ 4 Then hdr_sz = 24
If d EQ 5 Then hdr_sz = 28

hdr = BytArr(hdr_sz)
ReadU, lun, hdr

;get the number of blocks, slices, points from the file header
nfreq = bytarr(4)   ;number of complex freq points
nphase = bytarr(4)  ;number of phase encode points
nslice = bytarr(4)  ;number of slices
narray = bytarr(4)  ;number of array values
nvec = bytarr(4)    ;number of vector values

For i=0,3 Do Begin
    nphase(i) = hdr(i)
    nfreq(i) = hdr(i+4)
EndFor

If d EQ 2 Then Begin
    nslice = [0, 0, 0, 1]
EndIf Else Begin
    For i=0,3 Do nslice(i) = hdr(i+8)
EndElse

If d EQ 4 Then For i=0,3 Do Begin
    narray(i) = hdr(i+12)
    adim = 256^3*narray(0) + 256^2*narray(1) + 256*narray(2) + narray(3)
EndFor

If d EQ 5 Then For i=0,3 Do Begin
    narray(i) = hdr(i+12)
    adim = 256^3*narray(0) + 256^2*narray(1) + 256*narray(2) + narray(3)
    nvec(i) = hdr(i+16)
    vdim = 256^3*nvec(0) + 256^2*nvec(1) + 256*nvec(2) + nvec(3)
EndFor

; Size of freq, phase, and slice dimensions
fdim = 256^3*nfreq(0) + 256^2*nfreq(1) + 256*nfreq(2) + nfreq(3)
pdim = 256^3*nphase(0) + 256^2*nphase(1) + 256*nphase(2) + nphase(3)
sdim = 256^3*nslice(0) + 256^2*nslice(1) + 256*nslice(2) + nslice(3)

If d EQ 2 Then data_in = Make_Array(pdim,fdim,sdim)
If d EQ 3 Then data_in = Make_Array(pdim,fdim,sdim)
If d EQ 4 Then data_in = Make_Array(pdim,fdim,sdim,adim)
If d EQ 5 Then data_in = Make_Array(pdim,fdim,sdim,adim,vdim)

ReadU, lun, data_in
Free_LUN, lun

Return, data_in

End


; Subroutine name: adt_smooth_event
; Created by
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting


pro adt_smooth_event, ev
    compile_opt idl2
    common scan_data
    CI = project.ci

    widget_control, ev.id, get_uvalue=uvalue

    case uvalue of
       'linear': project.procPramArray[CI].regression_type = -1
       'linear_smooth': project.procPramArray[CI].regression_type = 0
       'Non-Linear_smooth': project.procPramArray[CI].regression_type = 2
       'Unified_smooth': project.procPramArray[CI].regression_type = 14
       'Alpha':project.procPramArray[CI].alpha = float(ev.VALUE)
       'Beta' :project.procPramArray[CI].beta = float(ev.VALUE)
       'Iterations':project.procPramArray[CI].iterations = fix(ev.VALUE)
       else:print, 'uvalue not recgnized'
    endcase
    project.procpramArray[CI].adt_proccess_flag = 0
end


; Subroutine name: mas_smooth
; Created by
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_smooth
    compile_opt idl2
    common scan_data
    COMMON common_widgets
    CI = project.ci

    ;If the ADT ROI window is already on the screen then dont make another one.
    IF N_ELEMENTS(ADT_ROI_WINDOW_BASE) EQ 1 THEN $
        IF     widget_info( ADT_ROI_WINDOW_BASE,  /VALID_ID ) eq 1 THEN return

    ;if the current scan does not have a good b-matrix then make every thing insensitive
    SENSITIVE = 0
    if !d.name EQ 'X' then begin
       SENSITIVE = 1
    end

    title = string ('Choose Fit Model', strtrim(CI,2))
    ADT_SMOOTH_WINDOW_BASE = widget_base(TITLE=title, $
            UVALUE = 'ADT_SMOOTH_WINDOW_BASE', $
            XOFFSET=420 ,YOFFSET=637, /column  )


    w=widget_label(ADT_SMOOTH_WINDOW_BASE, value='Fit Models')

    tlb=widget_base(ADT_SMOOTH_WINDOW_BASE, /EXCLUSIVE)

    b0=widget_button(tlb,value='Linear Model',uvalue='linear', SENSITIVE=SENSITIVE )
    b1=widget_button(tlb,value='Linear Model with Smoothing',uvalue='linear_smooth', SENSITIVE=SENSITIVE )
    b2=widget_button(tlb,value='Non-Linear Model wiht Smoothing',uvalue='Non-Linear_smooth' , SENSITIVE=SENSITIVE )
    b3=widget_button(tlb,value='Unified Non-Linear with smoothing',  uvalue='Unified_smooth' , SENSITIVE=SENSITIVE )

    if project.procPramArray[CI].regression_type eq -1 then WIDGET_CONTROL, /SET_BUTTON, b0
    if project.procPramArray[CI].regression_type eq 0  then WIDGET_CONTROL, /SET_BUTTON, b1
    if project.procPramArray[CI].regression_type eq 2  then WIDGET_CONTROL, /SET_BUTTON, b2
    if project.procPramArray[CI].regression_type eq 14 then WIDGET_CONTROL, /SET_BUTTON, b3

    w=cw_field( ADT_SMOOTH_WINDOW_BASE, /floating, /RETURN_EVENTS,  $
       value = project.procPramArray[CI].alpha ,title=' S0 smooth', uvalue='Alpha')

    w=cw_field( ADT_SMOOTH_WINDOW_BASE, /floating, /RETURN_EVENTS, $
       value = project.procPramArray[CI].beta  ,title='ADT smooth', uvalue='Beta')

    w=cw_field( ADT_SMOOTH_WINDOW_BASE, /INTEGER , /RETURN_EVENTS, $
       value = project.procPramArray[CI].iterations ,title='Iterations', uvalue='Iterations')

;    if project.procPramArray[CI].adt_img_stat_flg_array[1] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bXX
;    if project.procPramArray[CI].adt_img_stat_flg_array[2] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bYY
;    if project.procPramArray[CI].adt_img_stat_flg_array[3] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bZZ
;    if project.procPramArray[CI].adt_img_stat_flg_array[4] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bXY
;    if project.procPramArray[CI].adt_img_stat_flg_array[5] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bXZ
;    if project.procPramArray[CI].adt_img_stat_flg_array[6] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bYZ
;    if project.procPramArray[CI].adt_img_stat_flg_array[7] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bAD
;    if project.procPramArray[CI].adt_img_stat_flg_array[8] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bFA
;    if project.procPramArray[CI].adt_img_stat_flg_array[9] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bE1
;    if project.procPramArray[CI].adt_img_stat_flg_array[10] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bE2
;    if project.procPramArray[CI].adt_img_stat_flg_array[11] eq 1 then WIDGET_CONTROL, /SET_BUTTON, bE3


    widget_control, ADT_SMOOTH_WINDOW_BASE, /realize

    xmanager, 'adt_smooth',ADT_SMOOTH_WINDOW_BASE,/no_block, GROUP_LEADER=ADT_WINDOW_BASE


end


; Subroutine name: mas_do_smooth
; Created by
; Calling Information:

; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:

; Editing Information:
    ;Edited by HS 2006/10/04.
    ;Fix spelling mistakes and commenting

pro mas_do_smooth
    common scan_data
    CI = project.ci

    if !d.name eq 'WIN' then begin
       void= dialog_message('Smoothing can only be run on faraday')
       return
    end

    if project.procPramArray[CI].regression_type eq -1 then return

    ;send the user a message to make sure they want to run this program.
    ;if they don't then let them exit quietly.
;    user_message = [ $
;             'The Smoothing program can only be run on Fraraday', $
;             'Running elsewhere could make you computer act funny',$
;             '', $
;             'This process takes a long time to complete please be patient', $
;             '',$
;             'Do you wish to continue?']
;    if dialog_message(user_message, /question) eq 'No' then return


    ;since the users said that they want to continue then we will make
    ;a lot of assumptions first of all that we are running on faraday.
    ;this means that certan directories exist and there path is static

    ;in order for smooth to work we have to know the exact location of the
    ;executables will be at
    ;/export/faraday/thmlab/smooth
    ;and the exact location of all the work that is to be done will be at
    ;/export/faraday/thmlab/smooth/buffer

    buffer_pointer = '/export/faraday/thmlab/smooth/buffer/'

    ;first we will have to populate a b-matrix file.
    ;so we will need to bring in the b-matrix file and
    ;then start the writing process.

    ;we have to make sure that the data was extracted from the files
    mas_load_state_1

    ;bring in the b-value array
    bvals = *project.imndArray[CI].b_matrix

    ;create the file array that will be written to a file
    ;each collumn is seperated by a tab caracter
    sz_bvals = size(bvals)
    file_string = strarr(sz_bvals[2])
    for ii=0, sz_bvals[2]-1 do $
        file_string[ii] = $
           strtrim(bvals[0,ii],2)+string([09b])+ $
           strtrim(bvals[1,ii],2)+string([09b])+ $
           strtrim(bvals[2,ii],2)+string([09b])+ $
           strtrim(bvals[3,ii],2)+string([09b])+ $
           strtrim(bvals[4,ii],2)+string([09b])+ $
           strtrim(bvals[5,ii],2)

    ;now we write the array out to a file.

    buffer_directory = '/export/faraday/thmlab/smooth/buffer/'


    ;this file is the b_value file
    OpenW, lun0, buffer_directory+'b_val.txt', /SWAP_IF_LITTLE_ENDIAN, /get_lun
    printf, lun0,file_string
    Free_LUN, lun0


    ;now we create a dwi.flt(diffusion weighted image.floating point file) file
    sz_dwi = size(*project.dataArray[project.CI].state1)
    OpenW, lun0, buffer_directory+'dwi.flt', /SWAP_IF_LITTLE_ENDIAN, /get_lun
    WriteU, lun0, sz_dwi, *project.dataArray[project.CI].state1
    Free_LUN, lun0


    ;now that the data and b_val are written out we now write there locations
    ;to a data config file
    file_string = buffer_directory+'dwi.flt'+' '+buffer_directory+'b_val.txt'
    OpenW, lun0, buffer_directory+'data.cfg', /SWAP_IF_LITTLE_ENDIAN, /get_lun
    printf, lun0,file_string
    Free_LUN, lun0


    ;now we whave to make the config file for the whole process
    file_string = strarr(6)
    file_string[0] = buffer_directory+'data.cfg'
    file_string[1] = buffer_directory+'smoothed_dwi.flt'
    file_string[2] = strtrim(project.procPramArray[CI].regression_type,2)
    file_string[3] = strtrim(project.procPramArray[CI].alpha,2)+' '+strtrim(project.procPramArray[CI].beta,2)
    file_string[4] = '-1 -1 -1 -1 -1 -1 -1 -1'
    file_string[5] = strtrim(project.procPramArray[CI].iterations,2)

    OpenW, lun0, buffer_directory+'smooth.cfg', /SWAP_IF_LITTLE_ENDIAN, /get_lun
    printf, lun0,format= '(%"%s\n%s\n%s\n%s\n%s")', $
              file_string[0],file_string[1],file_string[2],file_string[3],file_string[4],file_string[5]
    Free_LUN, lun0

    update_status_bar,'Running. Please be patient, this could take hours'


    ;from this point all we have to do is spawn a new process and tell it where smooth.cfg is


    cd, '/export/faraday/thmlab/smooth/'

    spawn, './GenDT_MRI buffer/smooth.cfg'


    break = 0
    file_size = 0

    while break eq 0 do begin

       wait, 2
       ;get the parameters for the file that is going to be created
       ;and if the file exists then measure the size
       ;if the size is not equil the file_size varable then set file_size to the
       ;current size. else if the size equils file_size then no more data has
       ;been written to the file.
       fileParameters = FILE_INFO(buffer_directory+'smoothed_dwi.flt' )


       if fileParameters.exists eq 1 then $
         if fileParameters.SIZE ne file_size then file_size = fileParameters.SIZE $
         else break = 1
    end

    update_status_bar, 'files done'

;   project.dataArray[ci].smooth =

    ADT = rflt(buffer_directory+'smoothed_dwi.flt')
    ADT[*,*,*,1:*] /= 1000
    ;scale the image so that it ranges from 0 - 100.00 floating point
;   temp = max(adt)
;    adt = adt*100.0/temp

    fdim = project.imndArray[CI].fdim
    pdim = project.imndArray[CI].pdim
    sdim = project.imndArray[CI].sdim
    adim = project.imndArray[CI].adim

    project.procpramArray[CI].adt_proccess_flag = 1

     ; Make arrays for the alp, ADT and sigma matrices as well as 'Basser' result
;       ln_sig_1d = Make_Array(1, adim)
;       fit_result = Make_Array(1, 7)
;       ADT = Make_Array(fdim, pdim, sdim, 7, value = 0.0)
;       sigma = Make_Array(fdim, pdim, sdim, 6, value = 0.0)
    A_matrix = Make_array(3, 3)
    evals = Make_array(fdim, pdim, sdim, 3)
    evecs = Make_array(fdim, pdim, sdim, 3, 3)
    bra_D_ket = Make_array(fdim, pdim, sdim)
    Fract_Anisot = Make_Array(fdim,pdim,sdim)
    ;Rel_Anisot = Make_Array(fdim,pdim,sdim)

    nan_count = 0


    progressbar = Obj_New('progressbar', Color='red', Text='Calculating Other ADT Parameters',/NOCANCEL)
    progressbar -> Start

       For i = 0, fdim-1 Do Begin
         For j = 0, pdim-1 Do Begin
           For k = 0, sdim-1 Do Begin
;             If sig[i, j, k, 0] GT thr Then Begin


;                 ln_sig_1d = alog(sig[i, j, k, *])
;                   ln_sarr = Reform(ln_sig_1d, adim)

;          ;after we take the log we check for the presence of NAN Not A Number
;          ;if we find any skip to the next iteration
;                   void = where( FINITE(ln_sarr) ne 1, n_nan)
;                   if n_nan gt 0 then continue
;
;                fit_result = Regress(-bvals,ln_sarr,CONST=a0,SIGMA=sigmafit )

;                ADT[i, j, k, 0] = exp(a0)
;                ADT[i, j, k, 1:6] = fit_result
;                sigma[i, j, k, *] = sigmafit

                A_matrix = [[ADT[i,j,k,1], ADT[i,j,k,4], ADT[i,j,k,5]], $
                    [ADt[i,j,k,4], ADT[i,j,k,2], ADT[i,j,k,6]], $
                    [ADT[i,j,k,5], ADT[i,j,k,6], ADT[i,j,k,3]]]
                ev = 1    ; initialize eigenvectors (see IDL function EIGENQL)
                eigenval = EIGENQL(A_matrix, EIGENVECTORS = ev)
                evals[i, j, k, *] = eigenval
                evecs[i, j, k, *, *] = ev
                I1=ADT[i,j,k,1]+ADT[i,j,k,2]+ADT[i,j,k,3]
                I2=ADT[i,j,k,1]*ADT[i,j,k,2]+ADT[i,j,k,1]*ADT[i,j,k,3]+$
                    ADT[i,j,k,2]*ADT[i,j,k,3]-ADT[i,j,k,4]^2- $
                    ADT[i,j,k,5]^2-ADT[i,j,k,6]^2
                Bra_D_ket[i,j,k]=abs(I1/3.)

                temp = sqrt(1.-I2/(I1^2-2.*I2)) < 1.
          temp >= 0.
          void = where( FINITE(temp) ne 1, n_nan)
          nan_count += n_nan
                if n_nan gt 0 then temp *= 0

          Fract_Anisot[i,j,k]= temp


;             EndIf
          EndFor
         EndFor
         ;update_status_bar, string ('ADT fitting ', strtrim(i,2),'/',strtrim(fdim-1,2))
         progressBar -> Update, (float(i)/float(fdim-1))*100.0
       EndFor
       project.procpramArray[project.CI].adt_proccess_flag = 1
       progressbar -> Destroy

    ;****************************


;   print, 'nan count', nan_count
;
;
;   print, 'min and max of adt',min(adt, max=temp)
;   print, '                  ',temp
;   print, ''
;   print,  'min and max of eign_val',min(evals, max = temp)
;   print, '                        ',temp
;   print, ''
;   print, 'min and max of evecs',min(evecs, max=temp)
;   print, '                  ',temp
;   print, ''
;   print, 'min and max of Fract_Anisot',min(Fract_Anisot, max = temp)
;   print, '                       ',temp
;   print, ''
;   print, 'min and max of Bra_D_ket',min(Bra_D_ket, max = temp)
;   print, '                       ',temp


    update_status_bar, 'Saving structures to memory'
    project.dataArray[project.CI].adt = ptr_new(ADT,/no_copy)
    project.dataArray[project.CI].eign_val = ptr_new(evals,/no_copy)
    project.dataArray[project.CI].eign_Vec = ptr_new(evecs,/no_copy)
    project.dataArray[project.CI].frac_Ani = ptr_new(Fract_Anisot,/no_copy)
    project.dataArray[project.CI].Avg_Dif = ptr_new(Bra_D_ket,/no_copy)
    update_status_bar, ''


    heap_gc

end
