;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved
;
;This procedure will handle fiber tracking from mas. it will create a gui
;after which the user will enter the parameters they want specify an roi
;and the go from there. it will also let the user load in flt files.
;

; HS - 20070914
; NAME: seeds_to_3Dspace
; Called by: rd_slicer.pro
; Parameters: seedptr, fdim, pdim, sdim
; This function returns a 3D array of the size determined by fdim,pdim,sdim that
; contains the representation of the ROI in space.

function seeds_to_3Dspace, seedsptr, roi_fdim, roi_pdim, roi_sdim

	dataspace = bytarr(roi_fdim,roi_pdim,roi_sdim)

	num_rois = size(seedsptr)
		num_rois = num_rois[1]

	for i=0,num_rois-1 do begin

		seeds =  *seedsptr[i]

		xpoints = seeds[*,0]
		ypoints = seeds[*,1]
		zpoints = seeds[*,2]

		xpoints = fix(temporary(xpoints))
		ypoints = fix(temporary(ypoints))
		zpoints = fix(temporary(zpoints))

		xsize = size(xpoints)
		xsize = xsize[1]
		ysize = size(ypoints)
		ysize = ysize[1]
		zsize = size(zpoints)
		zsize = zsize[1]

		;print, seeds



		FOR k=0,zsize-1 DO BEGIN
			x_index = xpoints[k]
			y_index = ypoints[k]
			z_index = zpoints[k]
			dataspace[x_index,y_index,z_index] = 255

		ENDFOR
end
	;print, dataspace[*,*,7]
	return, dataspace

END


; Subroutine name: Build_fibers
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:
; Must include certain paths that contain subroutine calls not in the Fibers directory.

; Purpose of subroutine:
; Builds the .SAV file that runs Fibers with Virtual Machine.

; Editing Information:
    ;Edited by HS 2006/10/12
    ;Fix spelling mistakes and commenting.

;PRO Build_fibers
;    compile_opt idl2
;
;    RESOLVE_ALL , /CONTINUE_ON_ERROR
;
;    file = dialog_pickfile(/directory)
;    print, file
;    save, /routines, filename = file+'fibers.sav'
;    ;save, /routines, filename = 'Z:\thmlab\idl\mas\mas.sav'
;
;END

; HS 20061013 try to see event related info. not working.

;pro fibers2_event, event
;
;print, 'Event Detected'
;help, event, /structure
;
;End


; ------------------------------------------------------------------------------
; Fiber-Tract Mapping Software v2.0
; See User Manual for using Instructions
;
; Purpose of subroutine:
; Subroutine that generates the main form of the program and launches the appropriate
; subroutines.
;
; Written by E. Ozarslan on October 18, 2001
; Last modified by E. Ozarslan on August 31, 2003
; Edited by HS 2006/10/12
; Fix commenting.
; HS - 20061031
; Integrating into MAS
; ------------------------------------------------------------------------------


pro fibers2, defaults

    compile_opt idl2, hidden
    common scan_data
    ;;common share_fibers_main, pTHR, pFA, pS_0, pAvD, pDIR, loadedthr, loadeddir
    

    ;; Describes the components and values for the main form in Fibers.
    desc = [ $
             '1, base, , column, ', $
             '0, label, Project Directory, center', $
             '1, base, , row, frame', $
             '0, text,'+defaults.directory+',width=59, tag=path', $
             '1, base, , row', $
             '2, button, Change Directory, tag=direct, quit', $
             '2, label, , center',$
             '2, label, , center',$
             ;;
    '1, base, , column, ', $
      '0, label, Resolution , center', $
      '1, base, , row', $
      '1, base, , row, frame', $
      '0, float, '+defaults.resx+', label_left=  X:, width=10, tag=resx', $
      '0, float, '+defaults.resy+', label_left=             Y:, width=10, tag=resy', $
      '2, float, '+defaults.resz+', label_left=             Z:, width=10, tag=resz', $
      '2, label, , center',$
      '2, label, , center',$
      
    '1, base, , column, ', $
      '0, label, Parameters For Tracking, center', $
      '1, base, , column, frame, center', $
      '1, base, , row, ', $
      '1, base, , row, frame, center', $
      '0, button, Fractional Anisotropy|Average Diffusion|S0|Choose File,' $
      + 'frame, exclusive, label_top=ROI Selection File, set_value='+defaults.tfile+', tag=gfile', $
      '0, label, , center',$
      '0, label, , center',$
      '0, label, , center',$
      '1, base, ,column, ', $
      '0, label, Interpolation Scheme in the ROI window, center', $
      '1, base, ,column, frame', $
      '2, button,Nearest Neighbor|Linear|Congrid, exclusive,'$
      + 'set_value='+defaults.interpol+', tag=interpol', $
      '2, label, , center',$
      '2, label, , center',$
      '2, label, , center',$
      
    '1, base, , row, ', $
      '1, base, ,column, frame', $
      '0, droplist, Fractional Anisotropy|Directional Coherence|FA of coherence,' $
      + 'label_left=Thresholding File, set_value='+defaults.tfile+', tag=tfile', $
      '0, droplist, Diffusion Directions|Coherence Directions,' $
      + 'label_left=Directions File, set_value='+defaults.dfile+', tag=dfile', $
      '0, text, '+defaults.fname+', label_left=Output File Name:, width=17, tag=fname', $
      '0, float, '+defaults.ss+', label_left=Step Size    :, width=7, tag=ss', $
      '0, float, '+defaults.thr+', label_left=Threshold    :, width=7, tag=thr', $
      '2, float, '+defaults.sd+', label_left=Seed Density :, width=7, tag=sd', $
      '1, base, , column, ', $
      '0, label, Visualization Files, center', $
      '1, base, ,column, frame', $
      '2, button,From MAS|Tracts|Fractional Anisotropy|Directional Coherence|S0|Read,'$
      + 'column, set_value='+defaults.rfiles+', tag=rfiles', $
      '0, button, Repeat Tracking, tag=repeated, quit', $
      '2, button, Specify new ROI, tag=roi, quit', $
      '2, label, , center',$
      '2, label, , center',$
      '2, label, , center',$
      
    '1, base, , column, ', $
      '0, label, Visualization Tool, center', $
      '1, base, ,row, frame', $
      '0, integer, '+defaults.rfx+', label_left=  Visualization Factor, width=10, tag=rfx', $
      '0, float, '+defaults.ll+', label_left=   Length Limit, width=10, tag=ll', $
      '1, base, , row', $
      '2, button,  Visualize, tag=vis, quit', $
      '2, label, , center',$
      '2, label, , center',$
      
    '1, base, ,row', $
      '0, button, Set Defaults, tag=set_defs, ', $
      '2, button, Exit, tag=cancel, quit' $
      ]
    
    ;; Generates the main form. IDL specific call, see page 260 of IDL Reference Guide A-M.
    b=cw_form(desc, /column, title='Fiber Tract Mapping Software v2.0')
    
    ;; User presses the "Change Directory" button.
    if b.direct then begin
        
        thisdir=dialog_pickfile(/must_exist, /directory,$
                                Get_Path=path, Title='Choose the project folder')
        
        cd, thisdir
        
        defaults.directory=strjoin(strsplit(thisdir, '\', /extract), '\\')
        
        if ptr_valid(pTHR) then ptr_free, pTHR, pFA, pS_0, pAvD, pDIR
        
        fibers2, defaults
        
        if ptr_valid(pTHR) then ptr_free, pTHR, pFA, pS_0, pAvD, pDIR
        
        return
        
    endif
    
; User presses the "Set Defaults" button.
    if b.set_defs then begin
        defaults.gfile=b.gfile
        defaults.interpol=b.interpol
        defaults.tfile=b.tfile
        defaults.dfile=b.dfile
        defaults.fname=b.fname
        defaults.resx=b.resx
        defaults.resy=b.resy
        defaults.resz=b.resz
        defaults.rfx=b.rfx
        defaults.ss=b.ss
        defaults.thr=b.thr
        defaults.ll=b.ll
        defaults.sd=b.sd
        rfile0=strtrim(string((b.rfiles)[0]),2)
        rfile1=strtrim(string((b.rfiles)[1]),2)
        rfile2=strtrim(string((b.rfiles)[2]),2)
        rfile3=strtrim(string((b.rfiles)[3]),2)
        rfile4=strtrim(string((b.rfiles)[4]),2)
        rfile5=strtrim(string((b.rfiles)[5]),2)
        ;; These are the values of the checkboxes for what data do you want to use for drawing the ROIs and tracking
        defaults.rfiles='['+rfile0+'\,'+rfile1+'\,'+rfile2+'\,'+rfile3+'\,'+rfile4+'\'+rfile5+']'
    endif
    
; User presses the "Repeat Tracking" or "Specify new ROI" button.
    if b.repeated or b.roi then begin
        if b.tfile eq 0 then tfile='Fract_Anisotropy'
        if b.tfile eq 1 then tfile='dir_coherence'
        if b.tfile eq 2 then tfile='fran_coh'
        if b.dfile then dfile='coherence_dir' else dfile='Eigenvectors'
        
; User presses the "Repeat Tracking" button.
        if b.repeated then begin
            tractfile=dialog_pickfile(/must_exist, filter='*.trs', Title='Read *.trs file.')
            ft_track, b, pTHR, pDIR, ptr_seeds, seeds_from=tractfile
        endif
        
        
; User presses the "Specify new ROI" button.
        if b.roi then begin
            if b.gfile eq 0 then gfile='Fract_Anisotropy'
            
            if b.gfile eq 1 then gfile='Aver_D'
            
            if b.gfile eq 2 then gfile='S0'
            
            if b.gfile eq 3 then begin
                
                gfile=dialog_pickfile(/must_exist, filter='*.flt',$
                                      Get_Path=path, Title='Choose *.flt file.')
                slashpos=strpos(gfile, '\', /reverse_search)
                gfile=strmid(gfile, slashpos+1, strlen(gfile)-1)
                dotpos=strpos(gfile, '.', /reverse_search)
                gfile=strmid(gfile, 0, dotpos)
                
            endif
            
; ==================== THIS IS WHERE DATA FROM MAS IS LOADED =================
            
;; HS - 20061101
;; If the pointers to the DWI data are valid AND the user selects to
;; read from MAS then we will not require the selection of files.

	if strtrim(string((b.rfiles)[0]),2) eq 1            AND $
          ptr_valid(project.dataArray[project.CI].adt)      AND $
          ptr_valid(project.dataArray[project.CI].eign_val) AND $
          ptr_valid(project.dataArray[project.CI].eign_Vec) AND $
          ptr_valid(project.dataArray[project.CI].frac_Ani) AND $
          ptr_valid(project.dataArray[project.CI].Avg_Dif) THEN BEGIN
            
            print, 'Reading in available data from MAS'
            
;
;;	**** This is how the data is stored by the MAS_adt_regress.pro subroutine. ****
            
;;	 project.dataArray[project.CI].adt = ptr_new(ADT)
;;   project.dataArray[project.CI].eign_val = ptr_new(evals)
;;   project.dataArray[project.CI].eign_Vec = ptr_new(evecs)
;;   project.dataArray[project.CI].frac_Ani = ptr_new(Fract_Anisot)
;;   project.dataArray[project.CI].Avg_Dif = ptr_new(Bra_D_ket)

;;	 project.dataArray[project.CI].fibers_adt = ptr_new(ADT,/no_copy)
;;   project.dataArray[project.CI].fibers_eign_val = ptr_new(evals,/no_copy)
;;   project.dataArray[project.CI].fibers_eign_Vec = ptr_new(evecs,/no_copy)
;;   project.dataArray[project.CI].fibers_frac_Ani = ptr_new(Fract_Anisot,/no_copy)
;;   project.dataArray[project.CI].fibers_Avg_Dif = ptr_new(Bra_D_ket,/no_copy)

;;	**** End Example of Code ****

            ;; The next two lines are because of how Evren reads the data.
            ;; Originally I thought the pTHR was the "adt" matrix but the dimensions do not
            ;; agree. Further in the program the pTHR and pFA have the same data so that's
            ;; what I will be using here.
            
            pTHR=project.dataArray[project.CI].fibers_frac_Ani
            pFA=project.dataArray[project.CI].fibers_frac_Ani
            ;; pFA and pTHR are the SAME pointer!
            
            ;; HS - 20061105
            ;; s0 is stored like this by Evren
            ;; s0 = ptr_new((*project.dataArray[project.CI].adt)[*,*,*,0])
            ;; so I will read it in.
            pS_0_temp = (*project.dataArray[project.CI].fibers_adt)[*,*,*,0]
            pS_0= ptr_new(pS_0_temp)
            
            pAvD=project.dataArray[project.CI].fibers_Avg_Dif
            
            ;; Only the first 5th dimension is read in the eigenvector matrix. This
            ;; is how Evren did it for the files stored in the directory.
            pDIR_TEMP = reform((*project.dataArray[project.CI].fibers_eign_Vec)[*,*,*,*,0])
            pDIR=ptr_new(pDIR_TEMP)
            
            ;; We only need to rotate pTHR OR pFA because they are the same pointer
            mas_rotate_flip, pTHR
            ;;mas_rotate_flip, pFA
            
            ;; Rotate the other images
            mas_rotate_flip, pS_0
            mas_rotate_flip, pAvD
            mas_rotate_flip, pDIR
            
            ;; Following exactly what Evren did for the ROI selection
            ;; ptr_seeds=ft_roi(b.resx, b.resy, b.resz, b.interpol, file=gfile)
            ptr_seeds=ft_roi(b.resx, b.resy, b.resz, b.interpol)
            
            if total(ptr_valid(ptr_seeds)) then ft_track, b, pTHR, pDIR, ptr_seeds
            
            ;; Now I have to restore the original values of the ADT data in case the user closes Fibers and
            ;; restarts it later (because we need to reset the rotation of images).
            project.dataArray[project.CI].fibers_adt = ptr_new(*project.dataArray[project.CI].adt)
            project.dataArray[project.CI].fibers_eign_val = ptr_new(*project.dataArray[project.CI].eign_val)
            project.dataArray[project.CI].fibers_eign_Vec = ptr_new(*project.dataArray[project.CI].eign_Vec)
            project.dataArray[project.CI].fibers_frac_Ani = ptr_new(*project.dataArray[project.CI].frac_Ani)
            project.dataArray[project.CI].fibers_Avg_Dif = ptr_new(*project.dataArray[project.CI].Avg_Dif)
            
            
; ================ THIS IS WHERE ALL THE FILES ARE READ IN  (NOT FROM MAS) ===============
            
            ;; This happens when you are loading data from .flt files
        endif else begin
            
            if gfile ne '' and file_test(gfile+'.flt') then begin
                
                if ptr_valid(pS_0) eq 0 then begin
                    
                    if  file_test(tfile+'.flt') and file_test('Fract_Anisotropy.flt')$
                      and file_test('S0.flt') and file_test('Aver_D.flt') and $
                      file_test(dfile+'.flt') then begin
                        
                        pTHR   = ptr_new(rflt(tfile+'.flt'))
                        pFA    = ptr_new(rflt('Fract_Anisotropy.flt'))
                        pS_0   = ptr_new(rflt('S0.flt'))
                        pAvD   = ptr_new(rflt('Aver_D.flt'))
                        pDIR   = b.dfile ? ptr_new(rflt(dfile+'.flt')) : $
                          ptr_new(reform((rflt('Eigenvectors.flt'))[*,*,*,*,0]))
                        loadedthr = b.tfile
                        loadeddir = b.dfile
                        
                    endif else begin
                        ;; User has no valid files
                        ok=dialog_message('Please set the project directory.')
                        fibers2, defaults
                        
                    endelse
                    
                endif else begin
                    
                    
                    if b.tfile ne loadedthr then begin
                        if tfile eq 'Fract_Anisotropy' then pTHR=pFA else $
                          pTHR=ptr_new(rflt(tfile+'.flt'))
                        loadedthr=b.tfile
                    endif
                    
                    if b.dfile ne loadeddir then begin
                        pDIR=b.dfile ? ptr_new(rflt(dfile+'.flt')) : $
                          ptr_new(reform((rflt('Eigenvectors.flt'))[*,*,*,*,0]))
                        loadeddir=b.dfile
                    endif
                    
                endelse
                
                ptr_seeds=ft_roi(b.resx, b.resy, b.resz, b.interpol, file=gfile)
                if total(ptr_valid(ptr_seeds)) then ft_track, b, pTHR, pDIR, ptr_seeds
                
            endif
            
            ;; HS 20061101
            ;; This one closes the pointer check (reading data directly from MAS)
        endelse
        
    endif



endif









; User presses the "Visualize" button
if b.vis then rd_slicer, aspect_ratio=float(b.resz)/float(b.resx), $
     rf=[b.rfx, b.rfx, b.rfx], length_limit=float(b.ll)/float(b.ss)


; User presses the "Cancel" button
  if b.cancel ne 1 then fibers2, defaults
  ; We should not close these pointers in case the user relaunches Fibers from MAS
  ;if ptr_valid(pTHR) then ptr_free, pTHR, pFA, pS_0, pAvD, pDIR
end


; Subroutine name: fibers
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; Main program that creates the default values and then launches Fibers 2 which actually
; creates the main form.

; Editing Information:
    ;Edited by HS 2006/10/12
    ;Fix spelling mistakes and commenting.


pro fibers

  compile_opt idl2

  common scan_data

  window, /pixmap & wdelete
  device, bypass_translation=0
  device, retain=2

  cd, current=thisdir

  cd, thisdir

; ==== Going to send different defaults if there is data from MAS =====
if ptr_valid(project.dataArray[project.CI].adt) AND ptr_valid(project.dataArray[project.CI].eign_val) AND $
ptr_valid(project.dataArray[project.CI].eign_Vec) AND ptr_valid(project.dataArray[project.CI].frac_Ani) $
	AND ptr_valid(project.dataArray[project.CI].Avg_Dif) THEN BEGIN


  defaults={defaults, directory:strjoin(strsplit(thisdir, '\', /extract), '\\'), $
                gfile:'0', $
                tfile:'0', $
                dfile:'0', $
                fname:'tracts', $
                interpol:'1', $
                resx:'1', $
                resy:'1', $
                resz:'1', $
                rfx:'2', $
                ss:'0.1', $
                thr:'0.175', $
                ll:'10.0', $
                sd:'1.0', $
                rfiles:'[1\,0\,0\,0\,0\,0]' $
                }

endif else begin

defaults={defaults, directory:strjoin(strsplit(thisdir, '\', /extract), '\\'), $
                gfile:'0', $
                tfile:'0', $
                dfile:'0', $
                fname:'tracts', $
                interpol:'1', $
                resx:'1', $
                resy:'1', $
                resz:'1', $
                rfx:'2', $
                ss:'0.1', $
                thr:'0.175', $
                ll:'10.0', $
                sd:'1.0', $
                rfiles:'[0\,1\,1\,0\,0\,0]' $
                }


endelse

  fibers2, defaults
end

pro WID_BASE_Fibers_Main, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_

     common scan_Data
     common common_widgets

;    catch, error_state
;    if (error_state ne 0) then begin
;        catch, /cancel
;        return
;    endif

;;; Create the FIBERS data pointers now, instead of after ADT regressions, since
;;; If we are here, then we will actually need them. First free them, in case they're already
;;; defined, then recreate them.
    ptr_free, project.dataArray[project.CI].fibers_adt
    ptr_free, project.dataArray[project.CI].fibers_eign_val
    ptr_free, project.dataArray[project.CI].fibers_eign_Vec
    ptr_free, project.dataArray[project.CI].fibers_frac_Ani
    ptr_free, project.dataArray[project.CI].fibers_Avg_Dif

; These are the ones that will be edited in MAS and Fibers
    project.dataArray[project.CI].fibers_adt = ptr_new(*project.dataArray[project.CI].adt)
    project.dataArray[project.CI].fibers_eign_val = ptr_new(*project.dataArray[project.CI].eign_val)
    project.dataArray[project.CI].fibers_eign_Vec = ptr_new(*project.dataArray[project.CI].eign_Vec)
    project.dataArray[project.CI].fibers_frac_Ani = ptr_new(*project.dataArray[project.CI].frac_Ani)
    project.dataArray[project.CI].fibers_Avg_Dif = ptr_new(*project.dataArray[project.CI].Avg_Dif)


    
   WID_BASE_Fibers_Main = WIDGET_BASE(UNAME='WID_BASE_Fibers_Main', TITLE='Fiber Track Mapping', $
                                         Notify_realize='WID_BASE_Main_Realized', /column)
   
   main_rows = widget_base(WID_BASE_Fibers_Main, /column)
   top_portion = widget_base(main_rows, /row)
   base_left = widget_base(top_portion, /column, /frame)
   base_right = widget_base(top_portion, /column, /frame)
   
  ;; BUTTONS
  
  base_btns = widget_base(base_left, row=4, /grid)

  WID_BUTTON_SpecifyROI = Widget_Button(base_btns,  $
      UNAME='WID_BUTTON_SpecifyROI' ,EVENT_PRO='WID_BUTTON_SpecifyROI_Pressed',/ALIGN_CENTER , $
      VALUE='Specify ROI')

  WID_BUTTON_RepeatTrack = Widget_Button(base_btns,  $
      UNAME='WID_BUTTON_RepeatTrack' ,EVENT_PRO='WID_BUTTON_RepeatTrack_Pressed', /ALIGN_CENTER, $
      VALUE='Repeat Track')

  WID_BUTTON_SpecifyDirectory = Widget_Button(base_btns,  $
      UNAME='WID_BUTTON_SpecifyDirectory' ,EVENT_PRO='WID_BUTTON_SpecifyDirectory_Pressed', /ALIGN_CENTER, $
      VALUE='Select Directory')

  WID_BUTTON_RevisualizeTracts = Widget_Button(base_btns,  $
      UNAME='WID_BUTTON_Revisualize' ,EVENT_PRO='WID_BUTTON_Revisualize_Pressed', /ALIGN_CENTER, $
      VALUE='Revisualize')

 ;; RESOLUTION
  base_res = widget_base(base_left, /column, /frame)
  
  WID_LABEL_ResolutionTitle =  $
      Widget_Label(base_res,  $
      UNAME='WID_LABEL_ResolutionTitle',/ALIGN_CENTER ,VALUE='Resolution')

   base_val = widget_base(base_res, /row)
   
  WID_LABEL_Res_X = Widget_Label(base_val,  $
      UNAME='WID_LABEL_Res_X',/ALIGN_LEFT  $
      ,VALUE='X:')
  WID_Text_Res_X_Change = Widget_Text(base_val,  $
      UNAME='WID_Text_Res_X_Change' $
      ,/ALIGN_LEFT ,VALUE=' 1', xsize=8)

  base_val = widget_base(base_res, /row)

  WID_LABEL_Res_Y = Widget_Label(base_val,  $
      UNAME='WID_LABEL_Res_Y',/ALIGN_LEFT  $
      ,VALUE='Y:')

  WID_Text_Res_Y_Change = Widget_Text(base_val,  $
      UNAME='WID_Text_Res_Y_Change' $
      ,/ALIGN_LEFT ,VALUE=' 1', xsize=8)

 base_val = widget_base(base_res, /row)
 
 WID_LABEL_0 = Widget_Label(base_val,  $
      UNAME='WID_LABEL_0' ,/ALIGN_LEFT  $
      ,VALUE='Z:')

  WID_Text_Res_Z_Change = Widget_Text(base_val,  $
      UNAME='WID_Text_Res_Z_Change' $
      ,/ALIGN_LEFT ,VALUE=' 1', xsize=8)

  ;; CONTROLS
  base_controls = widget_base(base_right, column=2, /frame)

  temp_base = widget_base(base_controls, row=2)
  WID_LABEL_InterpolationScheme = Widget_Label(temp_base,  $
      UNAME='WID_LABEL_InterpolationScheme', /align_left $
      ,VALUE='Interpolation Scheme:')
  
  WID_DROPLIST_InterpolationScheme = Widget_Droplist(temp_base,  $
      UNAME='WID_DROPLIST_InterpolationScheme' ,EVENT_PRO='WID_DROPLIST_InterpolationScheme_Change'  $
      ,VALUE=[ 'Linear', 'Nearest Neighbor', 'Congrid' ], /align_left )


  temp_base = widget_base(base_controls, row=2)
  WID_LABEL_ROI_Selection = Widget_Label(temp_base,  $
      UNAME='WID_LABEL_ROI_Selection'  $
      ,VALUE='ROI Selection Image:', /align_left )

  WID_DROPLIST_ROI_Selection = Widget_Droplist(temp_base,  $
      UNAME='WID_DROPLIST_ROI_Selection' ,EVENT_PRO='WID_DROPLIST_ROI_Selection_Change' $
      ,VALUE=[ 'Fractional Anisotropy', 'Average Diffusion', 'S0', 'Orientation' ], /align_left )

  WID_LABEL_ThresholdingFile = Widget_Label(base_controls,  $
      UNAME='WID_LABEL_ThresholdingFile' $
      ,VALUE='Thresholding by:', /align_left )

  WID_DROPLIST_Thresholding = Widget_Droplist(base_controls,  $
      UNAME='WID_DROPLIST_Thresholding' ,EVENT_PRO='WID_DROPLIST_Thresholding_Change'  $
      ,VALUE=[ 'Fractional Anisotropy', 'Generalized Anisotropy' ], /align_left )

  ;; OTHER CONTROLS
  
  WID_LABEL_Parameters = Widget_Label(base_controls,  $
      UNAME='WID_LABEL_Parameters' ,XOFFSET=135 ,YOFFSET=110  $
      ,/ALIGN_LEFT ,VALUE='Parameters')
  
  row = widget_base(base_controls, column=3, /grid)
  
  temp_base = widget_base(row, column=1, /frame)
  
  WID_LABEL_StepSize = Widget_Label(temp_base,  $
      UNAME='WID_LABEL_StepSize' ,/ALIGN_LEFT  $
      ,VALUE='Step Size:')

  WID_TEXT_StepSizeValue = Widget_Text(temp_base,  $
      UNAME='WID_TEXT_StepSizeValue' ,/ALL_EVENTS, EVENT_PRO='WID_TEXT_StepSizeValue_Event' $
      ,/EDITABLE ,VALUE=[ ' 0.100000' ]  $
      ,XSIZE=10)

  WID_SLIDER_StepSize = Widget_Slider(temp_base,  $
      UNAME='WID_SLIDER_StepSize' ,EVENT_PRO='WID_SLIDER_StepSize_ChangeValue' $
      ,/SUPPRESS_VALUE ,MAXIMUM=20 ,VALUE=1)

  temp_base = widget_base(row, column=1, /frame)
  
  WID_LABEL_ThresholdParameter = Widget_Label(temp_base,  $
      UNAME='WID_LABEL_ThresholdParameter' $
      ,/ALIGN_LEFT ,VALUE='Threshold:')

  WID_TEXT_ThresholdValue = Widget_Text(temp_base,  $
      UNAME='WID_TEXT_ThresholdValue',/ALL_EVENTS, EVENT_PRO='WID_TEXT_ThresholdValue_Event'  $
      ,/EDITABLE ,VALUE=[ ' 0.175000' ]  $
      ,XSIZE=10 )
  
  WID_SLIDER_Threshold = Widget_Slider(temp_base,  $
      UNAME='WID_SLIDER_Threshold' ,EVENT_PRO='WID_SLIDER_Threshold_ChangeValue'  $
      ,/SUPPRESS_VALUE ,MAXIMUM=200 ,VALUE=35)

  temp_base = widget_base(row, column=1, /frame)

  WID_LABEL_SeedDensity = Widget_Label(temp_base,  $
      UNAME='WID_LABEL_SeedDensity'  $
      ,/ALIGN_LEFT ,VALUE='Seed Density:')
  
  WID_TEXT_SeedDensityValue = Widget_Text(temp_base,  $
      UNAME='WID_TEXT_SeedDensityValue',/ALL_EVENTS, EVENT_PRO='WID_TEXT_SeedDensityValue_Event'  $
      ,SCR_XSIZE=45 , /EDITABLE ,VALUE=[ ' 1' ] ,XSIZE=20)

  WID_SLIDER_SeedDensity = Widget_Slider(temp_base,  $
      UNAME='WID_SLIDER_SeedDensity' ,EVENT_PRO='WID_SLIDER_SeedDensity_ChangeValue' $
      ,/SUPPRESS_VALUE ,MAXIMUM=5 ,VALUE=1)

  row = widget_base(base_controls, column=3, /grid)
  temp_base = widget_base(row, column=1, /frame)
  
  WID_LABEL_LengthLimit = Widget_Label(temp_base,  $
      UNAME='WID_LABEL_LengthLimit' $
      ,/ALIGN_LEFT ,VALUE='Length Limit:')
  
  WID_TEXT_LengthLimitValue = Widget_Text(temp_base,  $
      UNAME='WID_TEXT_LengthLimitValue' ,/ALL_EVENTS ,EVENT_PRO='WID_TEXT_LengthLimitValue_Event' $
      ,SCR_XSIZE=65 , /EDITABLE ,VALUE=[ ' 10.0000' ]  $
      ,XSIZE=45)

  WID_SLIDER_LengthLimit = Widget_Slider(temp_base,  $
      UNAME='WID_SLIDER_LengthLimit',EVENT_PRO='WID_SLIDER_LengthLimit_ChangeValue'  $
      ,/SUPPRESS_VALUE ,MINIMUM=1 ,MAXIMUM=200 ,VALUE=100)

  temp_base = widget_base(row, column=1, /frame)
  
  WID_LABEL_Visualization = Widget_Label(temp_base,  $
      UNAME='WID_LABEL_Visualization' $
      ,/ALIGN_LEFT ,VALUE='Vis. Factor:')
  
  WID_TEXT_VisualizationFactorValue = Widget_Text(temp_base,  $
      UNAME='WID_TEXT_VisualizationFactorValue',/ALL_EVENTS ,EVENT_PRO='WID_TEXT_VisualizationFactorValue_Event' $
       ,SCR_XSIZE=45 , /EDITABLE ,VALUE=[ ' 2' ] ,XSIZE=20)

  WID_SLIDER_VisualizationFactor = Widget_Slider(temp_base,  $
      UNAME='WID_SLIDER_VisualizationFactor' ,EVENT_PRO='WID_SLIDER_VisualizationFactor_ChangeValue' $
      ,/SUPPRESS_VALUE ,MINIMUM=1 ,MAXIMUM=4 ,VALUE=2)

  temp_base = widget_base(row, column=1, /frame)
  
  WID_LABEL_AngleLimit = Widget_Label(temp_base,  $
      UNAME='WID_LABEL_AngleLimit' $
      ,/ALIGN_LEFT ,VALUE='Angle Cutoff:')
  
  WID_TEXT_AngleLimitValue = Widget_Text(temp_base,  $
      UNAME='WID_TEXT_AngleLimitValue',/ALL_EVENTS ,EVENT_PRO='WID_TEXT_AngleLimitValue_Event' $
       ,SCR_XSIZE=45 , /EDITABLE ,VALUE=[ ' 75' ] ,XSIZE=20)

  WID_SLIDER_AngleLimit = Widget_Slider(temp_base,  $
      UNAME='WID_SLIDER_AngleLimit',EVENT_PRO='WID_SLIDER_AngleLimit_ChangeValue' $
      ,/SUPPRESS_VALUE ,MINIMUM=0 ,MAXIMUM=90 ,VALUE=75)

  ;; PATH & FILES
  
  base = widget_base(base_right, /column)
  
  WID_LABEL_OutputFile = Widget_Label(base,  $
      UNAME='WID_LABEL_OutputFile'  $
      ,/ALIGN_LEFT ,VALUE='Output Filename:')

  WID_TEXT_Output_Filename = Widget_Text(base,  $
      UNAME='WID_TEXT_Output_Filename' $
      ,/EDITABLE ,VALUE=[ 'tracts' ] ,XSIZE=70 )

  WID_LABEL_Transform_File = Widget_Label(base,  $
      UNAME='WID_LABEL_Transform_File'  $
      ,/ALIGN_LEFT ,VALUE='Transformation Matrix File:')

  goof = widget_base(base, /row)
  WID_TEXT_Transform_File = Widget_Text(goof,  $
      UNAME='WID_TEXT_Transform_File' $
      ,/EDITABLE ,VALUE=[ 'transform.matrix' ] ,XSIZE=70 )
  WID_BUTTON_Transform_File = widget_button(goof, value="Choose...", $
      event_pro='WID_BUTTON_Transform_File_Pressed')
  
  WID_LABEL_Project_Dataset = Widget_Label(base, $
      UNAME='WID_LABEL_Project_Dataset', $
      /ALIGN_LEFT, VALUE='Project onto dataset number (from scan list):')
  
  WID_TEXT_Project_Dataset = Widget_Text(base, $
      UNAME='WID_TEXT_Project_Dataset', $
      /EDITABLE, VALUE=[ string(project.ci, format='(I2)') ])
      
  widget_control, wid_base_fibers_main, /realize
  XManager, 'WID_BASE_Fibers_Main', WID_BASE_Fibers_Main, /NO_BLOCK, GROUP_LEADER=WID_BASE_MAIN
  
end

pro xWID_BASE_Fibers_Main, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_

	HEAP_GC
    COMPILE_OPT IDL2
    COMMON scan_data
    COMMON common_widgets
    CI = project.CI

;;; Create the FIBERS data pointers now, instead of after ADT regressions, since
;;; If we are here, then we will actually need them. First free them, in case they're already
;;; defined, then recreate them.
    ptr_free, project.dataArray[project.CI].fibers_adt
    ptr_free, project.dataArray[project.CI].fibers_eign_val
    ptr_free, project.dataArray[project.CI].fibers_eign_Vec
    ptr_free, project.dataArray[project.CI].fibers_frac_Ani
    ptr_free, project.dataArray[project.CI].fibers_Avg_Dif

; These are the ones that will be edited in MAS and Fibers
    project.dataArray[project.CI].fibers_adt = ptr_new(project.dataArray[project.CI].adt)
    project.dataArray[project.CI].fibers_eign_val = ptr_new(project.dataArray[project.CI].eign_val)
    project.dataArray[project.CI].fibers_eign_Vec = ptr_new(project.dataArray[project.CI].eign_Vec)
    project.dataArray[project.CI].fibers_frac_Ani = ptr_new(project.dataArray[project.CI].frac_Ani)
    project.dataArray[project.CI].fibers_Avg_Dif = ptr_new(project.dataArray[project.CI].Avg_Dif)

  WID_BASE_Fibers_Main = Widget_Base(UNAME='WID_BASE_Fibers_Main' ,XOFFSET=420 ,YOFFSET=380  $
      ,SCR_XSIZE=500 ,SCR_YSIZE=400  $
      ,NOTIFY_REALIZE='WID_BASE_Main_Realized' ,TITLE='Fiber Track'+ $
      ' Mapping' ,SPACE=3 ,XPAD=3 ,YPAD=3)

  WID_BASE_Fibers_DisplayStats = Widget_Base(WID_BASE_Fibers_Main,  $
      UNAME='WID_BASE_Fibers_DisplayStats' ,FRAME=1 ,XOFFSET=10  $
      ,YOFFSET=260 ,SCR_XSIZE=90 ,SCR_YSIZE=100 ,TITLE='IDL' ,SPACE=3  $
      ,XPAD=3 ,YPAD=3)

  WID_LABEL_ResolutionTitle =  $
      Widget_Label(WID_BASE_Fibers_DisplayStats,  $
      UNAME='WID_LABEL_ResolutionTitle' ,XOFFSET=20 ,YOFFSET=5  $
      ,SCR_XSIZE=50 ,/ALIGN_CENTER ,VALUE='Resolution')


  WID_LABEL_Res_X = Widget_Label(WID_BASE_Fibers_DisplayStats,  $
      UNAME='WID_LABEL_Res_X' ,XOFFSET=10 ,YOFFSET=30 ,/ALIGN_LEFT  $
      ,VALUE='X:')


  WID_LABEL_Res_Y = Widget_Label(WID_BASE_Fibers_DisplayStats,  $
      UNAME='WID_LABEL_Res_Y' ,XOFFSET=10 ,YOFFSET=55 ,/ALIGN_LEFT  $
      ,VALUE='Y:')


  WID_LABEL_0 = Widget_Label(WID_BASE_Fibers_DisplayStats,  $
      UNAME='WID_LABEL_0' ,XOFFSET=10 ,YOFFSET=80 ,/ALIGN_LEFT  $
      ,VALUE='Z:')


  WID_Text_Res_X_Change = Widget_Text(WID_BASE_Fibers_DisplayStats,  $
      UNAME='WID_Text_Res_X_Change' ,XOFFSET=25 ,YOFFSET=25  $
      ,/ALIGN_LEFT ,VALUE=' 1',SCR_XSIZE=60 )


  WID_Text_Res_Y_Change = Widget_Text(WID_BASE_Fibers_DisplayStats,  $
      UNAME='WID_Text_Res_Y_Change' ,XOFFSET=25 ,YOFFSET=50  $
      ,/ALIGN_LEFT ,VALUE=' 1',SCR_XSIZE=60)


  WID_Text_Res_Z_Change = Widget_Text(WID_BASE_Fibers_DisplayStats,  $
      UNAME='WID_Text_Res_Z_Change' ,XOFFSET=25 ,YOFFSET=75  $
      ,/ALIGN_LEFT ,VALUE=' 1',SCR_XSIZE=60)

  WID_BUTTON_SpecifyROI = Widget_Button(WID_BASE_Fibers_Main,  $
      UNAME='WID_BUTTON_SpecifyROI' ,EVENT_PRO='WID_BUTTON_SpecifyROI_Pressed',XOFFSET=10 ,YOFFSET=10  $
      ,SCR_XSIZE=90 ,SCR_YSIZE=30 ,/ALIGN_CENTER ,VALUE='Specify'+ $
      ' ROI')

  WID_BUTTON_RepeatTrack = Widget_Button(WID_BASE_Fibers_Main,  $
      UNAME='WID_BUTTON_RepeatTrack' ,EVENT_PRO='WID_BUTTON_RepeatTrack_Pressed',XOFFSET=10 ,YOFFSET=45  $
      ,SCR_XSIZE=90 ,SCR_YSIZE=30 ,/ALIGN_CENTER ,VALUE='Repeat'+ $
      ' Track')

  WID_BUTTON_SpecifyDirectory = Widget_Button(WID_BASE_Fibers_Main,  $
  UNAME='WID_BUTTON_SpecifyDirectory' ,EVENT_PRO='WID_BUTTON_SpecifyDirectory_Pressed',XOFFSET=10 ,YOFFSET=80  $
      ,SCR_XSIZE=90 ,SCR_YSIZE=30 ,/ALIGN_CENTER ,VALUE='Select Directory')

  WID_BUTTON_RevisualizeTracts = Widget_Button(WID_BASE_Fibers_Main,  $
  UNAME='WID_BUTTON_Revisualize' ,EVENT_PRO='WID_BUTTON_Revisualize_Pressed',XOFFSET=10 ,YOFFSET=115  $
      ,SCR_XSIZE=90 ,SCR_YSIZE=30 ,/ALIGN_CENTER ,VALUE='Revisualize')

  WID_BASE_0 = Widget_Base(WID_BASE_Fibers_Main, UNAME='WID_BASE_0'  $
      ,FRAME=1 ,XOFFSET=115 ,YOFFSET=10 ,SCR_XSIZE=370 ,SCR_YSIZE=350  $
      ,TITLE='IDL' ,SPACE=3 ,XPAD=3 ,YPAD=3)


  WID_DROPLIST_InterpolationScheme = Widget_Droplist(WID_BASE_0,  $
      UNAME='WID_DROPLIST_InterpolationScheme' ,EVENT_PRO='WID_DROPLIST_InterpolationScheme_Change' ,XOFFSET=5 ,YOFFSET=25  $
      ,VALUE=[ 'Linear', 'Nearest Neighbor', 'Congrid' ])


  WID_LABEL_InterpolationScheme = Widget_Label(WID_BASE_0,  $
      UNAME='WID_LABEL_InterpolationScheme' ,XOFFSET=5 ,YOFFSET=5  $
      ,/ALIGN_LEFT ,VALUE='Interpolation Scheme')


  WID_LABEL_ROI_Selection = Widget_Label(WID_BASE_0,  $
      UNAME='WID_LABEL_ROI_Selection' ,XOFFSET=5 ,YOFFSET=60  $
      ,/ALIGN_LEFT ,VALUE='ROI Selection Image')


  WID_DROPLIST_ROI_Selection = Widget_Droplist(WID_BASE_0,  $
      UNAME='WID_DROPLIST_ROI_Selection' ,EVENT_PRO='WID_DROPLIST_ROI_Selection_Change',XOFFSET=5 ,YOFFSET=80  $
      ,VALUE=[ 'Fractional Anisotropy', 'Average Diffusion', 'S0', 'Orientation' ])


  WID_DROPLIST_Thresholding = Widget_Droplist(WID_BASE_0,  $
      UNAME='WID_DROPLIST_Thresholding' ,EVENT_PRO='WID_DROPLIST_Thresholding_Change',XOFFSET=200 ,YOFFSET=80  $
      ,VALUE=[ 'Fractional Anisotropy', 'Generalized Anisotropy' ])


  WID_LABEL_ThresholdingFile = Widget_Label(WID_BASE_0,  $
      UNAME='WID_LABEL_ThresholdingFile' ,XOFFSET=200 ,YOFFSET=60  $
      ,/ALIGN_CENTER ,VALUE='Thresholding by:')


  WID_LABEL_OutputFile = Widget_Label(WID_BASE_0,  $
      UNAME='WID_LABEL_OutputFile' ,XOFFSET=200 ,YOFFSET=5  $
      ,/ALIGN_LEFT ,VALUE='Output Filename:')


  WID_TEXT_Output_Filename = Widget_Text(WID_BASE_0,  $
      UNAME='WID_TEXT_Output_Filename' ,XOFFSET=200 ,YOFFSET=25  $
      ,/EDITABLE ,VALUE=[ 'tracts' ] ,XSIZE=20 )


  WID_LABEL_Parameters = Widget_Label(WID_BASE_0,  $
      UNAME='WID_LABEL_Parameters' ,XOFFSET=135 ,YOFFSET=110  $
      ,/ALIGN_LEFT ,VALUE='Parameters')


  WID_SLIDER_StepSize = Widget_Slider(WID_BASE_0,  $
      UNAME='WID_SLIDER_StepSize' ,EVENT_PRO='WID_SLIDER_StepSize_ChangeValue',XOFFSET=5 ,YOFFSET=160  $
      ,/SUPPRESS_VALUE ,MAXIMUM=20 ,VALUE=1)


  WID_LABEL_StepSize = Widget_Label(WID_BASE_0,  $
      UNAME='WID_LABEL_StepSize' ,XOFFSET=5 ,YOFFSET=135 ,/ALIGN_LEFT  $
      ,VALUE='Step Size:')


  WID_TEXT_StepSizeValue = Widget_Text(WID_BASE_0,  $
      UNAME='WID_TEXT_StepSizeValue' ,/ALL_EVENTS, EVENT_PRO='WID_TEXT_StepSizeValue_Event',XOFFSET=125 ,YOFFSET=130  $
      ,SCR_XSIZE=65 ,/EDITABLE ,VALUE=[ ' 0.100000' ]  $
      ,XSIZE=20)


  WID_LABEL_ThresholdParameter = Widget_Label(WID_BASE_0,  $
      UNAME='WID_LABEL_ThresholdParameter' ,XOFFSET=200 ,YOFFSET=135  $
      ,/ALIGN_LEFT ,VALUE='Threshold:')


  WID_TEXT_ThresholdValue = Widget_Text(WID_BASE_0,  $
      UNAME='WID_TEXT_ThresholdValue',/ALL_EVENTS, EVENT_PRO='WID_TEXT_ThresholdValue_Event' ,XOFFSET=285 ,YOFFSET=130  $
      ,SCR_XSIZE=65 ,/EDITABLE ,VALUE=[ ' 0.175000' ]  $
      ,XSIZE=20 )


  WID_SLIDER_Threshold = Widget_Slider(WID_BASE_0,  $
      UNAME='WID_SLIDER_Threshold' ,EVENT_PRO='WID_SLIDER_Threshold_ChangeValue',XOFFSET=200 ,YOFFSET=160  $
      ,/SUPPRESS_VALUE ,MAXIMUM=200 ,VALUE=35)


  WID_LABEL_SeedDensity = Widget_Label(WID_BASE_0,  $
      UNAME='WID_LABEL_SeedDensity' ,XOFFSET=5 ,YOFFSET=210  $
      ,/ALIGN_LEFT ,VALUE='Seed Density:')


  WID_TEXT_SeedDensityValue = Widget_Text(WID_BASE_0,  $
      UNAME='WID_TEXT_SeedDensityValue',/ALL_EVENTS, EVENT_PRO='WID_TEXT_SeedDensityValue_Event',XOFFSET=125 ,YOFFSET=205  $
      ,SCR_XSIZE=45 , /EDITABLE ,VALUE=[ ' 1' ] ,XSIZE=20)


  WID_SLIDER_SeedDensity = Widget_Slider(WID_BASE_0,  $
      UNAME='WID_SLIDER_SeedDensity' ,EVENT_PRO='WID_SLIDER_SeedDensity_ChangeValue',XOFFSET=5 ,YOFFSET=235  $
      ,/SUPPRESS_VALUE ,MAXIMUM=5 ,VALUE=1)


  WID_LABEL_LengthLimit = Widget_Label(WID_BASE_0,  $
      UNAME='WID_LABEL_LengthLimit' ,XOFFSET=200 ,YOFFSET=210  $
      ,/ALIGN_LEFT ,VALUE='Length Limit:')


  WID_TEXT_LengthLimitValue = Widget_Text(WID_BASE_0,  $
      UNAME='WID_TEXT_LengthLimitValue' ,/ALL_EVENTS ,EVENT_PRO='WID_TEXT_LengthLimitValue_Event',XOFFSET=280 ,YOFFSET=205  $
      ,SCR_XSIZE=65 , /EDITABLE ,VALUE=[ ' 10.0000' ]  $
      ,XSIZE=45)


  WID_SLIDER_LengthLimit = Widget_Slider(WID_BASE_0,  $
      UNAME='WID_SLIDER_LengthLimit',EVENT_PRO='WID_SLIDER_LengthLimit_ChangeValue',XOFFSET=200 ,YOFFSET=235  $
      ,/SUPPRESS_VALUE ,MINIMUM=1 ,MAXIMUM=200 ,VALUE=100)


  WID_LABEL_Visualization = Widget_Label(WID_BASE_0,  $
      UNAME='WID_LABEL_Visualization' ,XOFFSET=5 ,YOFFSET=285  $
      ,/ALIGN_LEFT ,VALUE='Visualization Factor:')


  WID_TEXT_VisualizationFactorValue = Widget_Text(WID_BASE_0,  $
      UNAME='WID_TEXT_VisualizationFactorValue',/ALL_EVENTS ,EVENT_PRO='WID_TEXT_VisualizationFactorValue_Event',XOFFSET=140  $
      ,YOFFSET=280 ,SCR_XSIZE=45 , /EDITABLE ,VALUE=[  $
      ' 2' ] ,XSIZE=20)


  WID_SLIDER_VisualizationFactor = Widget_Slider(WID_BASE_0,  $
      UNAME='WID_SLIDER_VisualizationFactor' ,EVENT_PRO='WID_SLIDER_VisualizationFactor_ChangeValue',XOFFSET=5 ,YOFFSET=310  $
      ,/SUPPRESS_VALUE ,MINIMUM=1 ,MAXIMUM=4 ,VALUE=2)


  WID_SLIDER_AngleLimit = Widget_Slider(WID_BASE_0, $
                                        UNAME='WID_SLIDER_AngleLimit', $
                                        xoffset=200, $
                                        yoffset=285, $
                                        title='Angle Cutoff (deg, keep <)', $
                                        minimum=0, $
                                        value=80,$
                                        maximum=90);, $
  ;;EVENT_PRO='WID_SLIDER_AngleLimit_Event')

  Widget_Control, /REALIZE, WID_BASE_Fibers_Main

  XManager, 'WID_BASE_Fibers_Main', WID_BASE_Fibers_Main, /NO_BLOCK, GROUP_LEADER=WID_BASE_MAIN

end

;
; Empty stub procedure used for autoloading.
; (Created by the GUI editor)
;
pro FibersGUI, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_
  WID_BASE_Fibers_Main, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_
end


pro mas_fibers



end
