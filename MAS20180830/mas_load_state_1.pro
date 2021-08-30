;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved

; Subroutine name: mas_kspace_shift
; Created by: Garrett Astary
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:
; Shift data in kspace.

; Editing Information:
 

pro mas_kspace_shift, p_data, ISDICOM=isdicom

    COMPILE_OPT IDL2
    COMMON scan_data

    CI = project.ci
    dim  =  project.imndArray[CI].dimensions
    k_fdim_shift = project.imndArray[CI].k_fdim_shift
    k_pdim_shift = -project.imndArray[CI].k_pdim_shift
    k_sdim_shift =  project.imndArray[CI].k_sdim_shift
    
    if project.procPramArray[ci].zpad_flag eq 1 then begin
        k_fdim_shift*=2
        k_pdim_shift*=2
        if dim eq 3 then k_sdim_shift*=2
    endif
    
    if (k_fdim_shift eq 0 and k_pdim_shift eq 0 and k_sdim_shift eq 0) then begin
        return ;; don't shift if theres no need
    endif
    
    sz_data = size(*p_data)
    fdim = sz_data[1]
    pdim = sz_data[2]
    sdim = sz_data[3]
    adim = sz_data[4]

    ;;shift the data according to it's dimensions
    ;;if the data is 2 dimensional then shift the data.
    if dim eq 2 then begin
        ;; Apply spatial shift in all dimensions

        If sz_data[0] EQ 4 then (*p_DATA) = temporary(shift((*p_DATA),k_fdim_shift,k_pdim_shift,k_sdim_shift,0))
        If sz_data[0] EQ 3 then (*p_DATA) = temporary(shift((*p_DATA),k_fdim_shift,k_pdim_shift,k_sdim_shift))
        If sz_data[0] EQ 2 then (*p_DATA) = temporary(shift((*p_DATA),k_fdim_shift,k_pdim_shift))

    endif

    ;;if the data is 3 dimensional then we have to shift this data differently
    if dim eq 3 then begin


        progressbar = Obj_New('progressbar', Color='red', Text='Shift Acquisition correction 1/2',/NOCANCEL)
        progressbar -> Start

        for ii=0, sdim-1 do begin
            (*p_DATA)[*,*,ii] = shift((*p_DATA)[*,*,ii],k_fdim_shift,k_pdim_shift)
            if ii mod 10 eq 0 then progressBar->Update, (float(ii)/float(sdim-1))*100.0
        endfor

        progressbar -> Destroy

        progressbar = Obj_New('progressbar', Color='red', Text='Shift Acquisition correction 2/2',/NOCANCEL)
        progressbar -> Start

        for ii=0, fdim-1 do begin
            for jj=0, pdim-1 do begin
                (*p_DATA)[ii,jj,*] = shift((*p_DATA)[ii,jj,*],k_sdim_shift)
            end
            if ii mod 10 eq 0 then progressBar->Update, (float(ii)/float(fdim-1))*100.0
        endfor

        progressbar -> Destroy

    endif

    update_status_bar,''

end

; Subroutine name: mas_kspace_subsample
; Created by: Bill Triplett
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; Subsamples kspace by truncating in the read/phase directions.

; Editing Information:
; Bug fixes (Magdoom, 09/20/15)

pro mas_kspace_subsample, pkspace

    common scan_data, project
    
    ci = project.ci
   
    ;; check if the scan is a true 3D scan
    is_3d = project.imndarray[ci].dimensions eq 3
    
    ;; the selected span
    new_xres = project.imndarray[ci].k_fdim_span
    new_yres = project.imndarray[ci].k_pdim_span
    new_zres = project.imndarray[ci].k_sdim_span

    ;; the max available span
    max_xres = project.imndarray[ci].k_fdim_span_max
    max_yres = project.imndarray[ci].k_pdim_span_max
    max_zres = project.imndarray[ci].k_sdim_span_max
    
    ;; Matrix dimensions
    fdim = project.imndarray[ci].fdim
    pdim = project.imndarray[ci].pdim
    sdim = project.imndarray[ci].sdim
    adim = project.imndarray[ci].adim
   
    if project.procPramArray[ci].zpad_flag eq 1 then begin
      new_xres*=2
      new_yres*=2
      if is_3d eq 1 then new_zres*=2
      
      max_xres*=2
      max_yres*=2
      if is_3d eq 1 then max_zres*=2 
      
      fdim*=2
      pdim*=2
      if is_3d eq 1 then sdim*=2
    endif
    
    *pkspace = reform(*pkspace, [fdim, pdim, sdim, adim])
    size_k = size(*pkspace, /dimensions)
    
    ;; check to make sure that the new values are at least in range.
    if (new_xres gt max_xres) then begin
        void = dialog_message('New X resolution cannot be larger than existing X resolution.', $
                              /error, /center)
        return
    endif else if (new_yres gt max_yres) then begin
        void = dialog_message('New Y resolution cannot be larger than existing Y resolution.', $
                              /error, /center)
        return    
    endif else if (new_zres gt max_zres) then begin
        void = dialog_message('New Z resolution cannot be larger than existing Z resolution.', $
                              /error, /center)
        return    
    endif
    
    ;; Create a new kspace
    if (is_3d eq 1) then begin
        new_k = complexarr(new_xres, new_yres, new_zres, adim)
    endif else begin
        new_k = complexarr(new_xres, new_yres, size_k[2], adim)
    endelse
    
    ;; loop for each volume
    for a = 0, adim-1 do begin
    
        ;; True 3D data sets are handled differently.
        if (is_3d eq 1) then begin
            
            new_k[*,*,*,a] = (*pkspace)[(size_k[0]-new_xres)/2:(size_k[0]+new_xres)/2-1, $
                                        (size_k[1]-new_yres)/2:(size_k[1]+new_yres)/2-1, $
                                        (size_k[2]-new_zres)/2:(size_k[2]+new_zres)/2-1, a]
            
        endif else begin
        
        ;; 2D Multislice data sets need a bit more overhead.
        ;; First we can resample kspace over the in-plane dimensions, then
        ;; after imagery we can reduce the resolution by interpolating in
        ;; the image domain. (Not implemented yet)
            for z = 0, size_k[2]-1 do begin
                new_k[*,*,z,a] = (*pkspace)[(size_k[0]-new_xres)/2:(size_k[0]+new_xres)/2-1, $
                                            (size_k[1]-new_yres)/2:(size_k[1]+new_yres)/2-1, $
                                            z, a]
            endfor
            
        endelse

    endfor

    *pkspace = temporary(new_k)

    ci = project.ci
    
;    ;; Set the new dimensions into the parameter array
;    project.imndarray[ci].fdim = new_xres
;    project.imndarray[ci].pdim = new_yres
;    
;    project.procpramarray[ci].fdim_start = project.procpramarray[ci].fdim_start < (new_xres-1)
;    project.procpramarray[ci].pdim_start = project.procpramarray[ci].pdim_start < (new_yres-1)
    
;    if (is_3d eq 1) then begin
;        project.imndarray[ci].sdim = new_zres
;        project.procpramarray[ci].sdim_start = project.procpramarray[ci].sdim_start < (new_zres-1)
;    endif
    
    mas_redraw_gui
    
end

; Subroutine name: mas_kspace_subsample_multislice
; Created by: Bill Triplett
; Calling Information:

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; Interpolated in the slice direction for resolution reduction
; in multi-slice data.

; Editing Information:
; Bug fixes (Magdoom, 09/20/15)

pro mas_kspace_subsample_multislice, pimagery

    common scan_data
    
    ci = project.ci
    
    ;; the selected span
    new_zres = project.imndarray[ci].k_sdim_span
    
    ;; the max available span
    max_zres = project.imndarray[ci].k_sdim_span_max
    
    ;; check if the scan is a true 3D scan
    is_3d = project.imndarray[ci].dimensions eq 3
    if (is_3d) then return
    
    new_xres = project.imndarray[ci].k_fdim_span
    new_yres = project.imndarray[ci].k_pdim_span
    fdim = project.imndarray[ci].fdim
    pdim = project.imndarray[ci].pdim
    sdim = project.imndarray[ci].sdim
    adim = project.imndarray[ci].adim
    
    if project.procPramArray[ci].zpad_flag eq 1 then begin
      new_xres*=2
      new_yres*=2
      fdim*=2
      pdim*=2
    endif 
    
    *pimagery = reform(*pimagery, new_xres, new_yres, sdim, adim)
    size_im = size(*pimagery, /dimensions)
    
    ;; Create a new imagery
    signal_type = project.procPramArray[CI].signal_type
    if (signal_type eq 9) then begin
        new_im = complexarr(new_xres, new_yres, new_zres, adim)
    endif else begin
        new_im = fltarr(new_xres, new_yres, new_zres, adim)
    endelse
    
    slice_interp = findgen(new_zres)/(new_zres) * size_im[2]
    ;; loop for each volume
    for a = 0, adim-1 do begin
    
        ;; 2D Multislice data sets need a bit more overhead.
        ;; First we can resample kspace over the in-plane dimensions, then
        ;; after imagery we can reduce the resolution by interpolating in
        ;; the image domain. (Not implemented yet)
        for x = 0, new_xres-1 do begin
            for y = 0, new_yres-1 do begin
            
                temp = reform((*pimagery)[x,y,*,a])
                new_im[x,y,*,a] = interpolate(temp, slice_interp, cubic=-0.5, missing=0)

            endfor
            
        endfor

    endfor

    *pimagery = temporary(new_im)

    ci = project.ci
    
;    ;; Set the new dimensions into the parameter array
;    project.imndarray[ci].fdim = new_xres
;    project.imndarray[ci].pdim = new_yres
;    project.imndarray[ci].sdim = new_zres
;    
;    project.procpramarray[ci].fdim_start = project.procpramarray[ci].fdim_start < (new_xres-1)
;    project.procpramarray[ci].pdim_start = project.procpramarray[ci].pdim_start < (new_yres-1)
;    project.procpramarray[ci].sdim_start = project.procpramarray[ci].sdim_start < (new_zres-1)
    
    mas_redraw_gui

end


; Subroutine name: mas_fft
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:



; Purpose of subroutine:
; Program FFTs the data depending on the dimension of the Data matrix.

; Editing Information:
; Edited by HS 2006/10/02
; Fix spelling mistakes and commenting.
; Optimized 3D data FFT (Magdoom, 09/11/15 and 09/23/16)
    
pro mas_fft , p_DATA
    COMPILE_OPT IDL2

    COMMON scan_data

    CI = project.ci

    ;check to see if the data needs to be ffted.
    image_type = project.imndArray[CI].image_type
    if image_type eq 5 or image_type eq 6 or image_type eq 10 or image_type eq 11 then return

    ;bring in the data.
    dim  = project.imndArray[CI].dimensions 
    fdim = project.imndarray[ci].k_fdim_span
    pdim = project.imndarray[ci].k_pdim_span
    sdim = project.imndarray[ci].k_sdim_span
;    pdim = project.imndArray[CI].pdim
;    sdim = project.imndArray[CI].sdim
    
    if project.procPramArray[CI].zpad_flag eq 1 then begin
        fdim*=2
        pdim*=2
        if dim eq 3 then sdim*=2
    endif 
    adim = project.imndArray[CI].adim
    signal_type = project.procPramArray[CI].signal_type

    ;if the data is 1 dimensional it must be spectroscopy
    if dim eq 1 then begin

       *p_DATA = temporary(fft(*p_DATA))
       *p_DATA = temporary(abs(*p_DATA))
       
       iPlot, *p_DATA

    end
    ;if the data is 2 dimensional then
    if dim eq 2 then begin
    ;2d fft

       ;update_status_bar , 'FFT'
       progressbar = Obj_New('progressbar', Color='red', Text='2D FFT', /NOCANCEL)
       progressbar -> Start

       For ii=0,adim-1 Do Begin
         For jj=0,sdim-1 Do Begin
            (*p_DATA)[*,*,jj,ii] = temporary(fft((*p_DATA)[*,*,jj,ii]))
            ;We need to account for FFT being performed from 0 -> N when kspace is ordered from -k -> k
            for x=0,fdim -1 do begin
                for y=0,pdim-1 do begin
                  if x mod 2 ne 0 and y mod 2 eq 0 then begin
                    (*p_data)[x,y,jj,ii] = -(*p_data)[x,y,jj,ii]
                  endif
                  if x mod 2 eq 0 and y mod 2 ne 0 then begin
                    (*p_data)[x,y,jj,ii] = -(*p_data)[x,y,jj,ii]
                  endif
               endfor
             endfor
         EndFor
         ;update_status_bar , string('FFT ',strtrim(ii,2),'/',strtrim(adim-1,2))
         progressBar -> Update, (float(ii)/float(adim-1))*100.0

       EndFor

       progressbar -> Destroy

    end
    if dim eq 3 then begin
        ;3d fft

        ;update_status_bar,'3d FFT'
        progressbar = Obj_New('progressbar', Color='red', Text='3D FFT',/NOCANCEL)
        progressbar -> Start
        for aa = 0, adim-1 do begin
          (*p_DATA)[*,*,*, aa] = fft((*p_DATA)[*,*,*,aa])
          progressBar -> Update, (float(aa)/float(adim))*100.0
        endfor
        progressbar -> Destroy
        
;        progressbar = Obj_New('progressbar', Color='red', Text='3D FFT Part 1/2',/NOCANCEL)
;        progressbar -> Start
;        for aa = 0, adim-1 do begin
;            For ii=0,sdim-1 Do Begin
;    
;                (*p_DATA)[*,*,ii, aa] = temporary(fft((*p_DATA)[*,*,ii,aa]))
;                 
;                if ii mod 10 eq 0 then begin
;                    ;update_status_bar,string('3d FFT ',strtrim(ii,2),'/',strtrim(sdim-1,2))
;                    progressBar -> Update, (float(ii)/float(sdim-1))*100.0
;                end
;            end
;        endfor
;        progressbar -> Destroy
;
;        ;update_status_bar,'3d FFT part 2.'
;        progressbar = Obj_New('progressbar', Color='red', Text='3D FFT Part 2/2',/NOCANCEL)
;        progressbar -> Start
;        for aa = 0, adim-1 do begin
;            For ii=0,fdim-1 Do begin
;                for jj=0,pdim-1 Do $
;                    (*p_DATA)[ii,jj,*,aa] = temporary(fft((*p_DATA)[ii,jj,*,aa]))
;                   
;                
;                if ii mod 10 eq 0 then begin
;                    ;update_status_bar,string('3d FFT ',strtrim(ii,2),'/',strtrim(fdim-1,2))
;                    progressBar -> Update, (float(ii)/float(fdim-1))*100.0
;                end
;            end
;        endfor
;        progressbar -> Destroy
        
        progressbar = Obj_New('progressbar', Color='red', $
                              Text='Correcting for FFT offset...', $
                              /NOCANCEL)
        progressbar->Start
  
        for ii=0, adim-1 do begin  
           (*p_data)[0:*:2,0:*:2,1:*:2,ii]=-(*p_data)[0:*:2,0:*:2,1:*:2,ii]
           (*p_data)[0:*:2,1:*:2,0:*:2,ii]=-(*p_data)[0:*:2,1:*:2,0:*:2,ii]
           (*p_data)[1:*:2,0:*:2,0:*:2,ii]=-(*p_data)[1:*:2,0:*:2,0:*:2,ii]
           (*p_data)[1:*:2,1:*:2,1:*:2,ii]=-(*p_data)[1:*:2,1:*:2,1:*:2,ii]     
          progressbar->update, float(ii)/float(adim) * 100
        endfor
        progressbar->Destroy
        
;        for ii=0, adim-1 do begin
;            for z =0, sdim-1 do begin
;              for x=0,fdim-1  do begin
;                for y=0, pdim -1 do begin
;                  if x mod 2 eq 0 and y mod 2 eq 0 and z mod 2 ne 0 then begin
;                        (*p_data)[x,y,z,ii] = -(*p_data)[x,y,z,ii]
;                  endif
;                  
;                  if x mod 2 eq 0 and y mod 2 ne 0 and z mod 2 eq 0 then begin
;                        (*p_data)[x,y,z,ii] = -(*p_data)[x,y,z,ii]
;                  endif
;                  
;                  if x mod 2 ne 0 and y mod 2 eq 0 and z mod 2 eq 0 then begin
;                        (*p_data)[x,y,z,ii] = -(*p_data)[x,y,z,ii]
;                  endif
;                  
;                  if x mod 2 ne 0 and y mod 2 ne 0 and z mod 2 ne 0 then begin
;                        (*p_data)[x,y,z,ii] = -(*p_data)[x,y,z,ii]
;                  endif
;                endfor
;               endfor
;               progressbar->update, float(z)/float(sdim) * 100
;              endfor  
;          endfor
;          progressbar->Destroy
    end
end


; Subroutine name: mas_shift
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:
; Shift data.

; Editing Information:
    ;Edited by HS 2006/10/02
    ;Fix spelling mistakes and commenting.

pro mas_shift, p_data, ISDICOM=isdicom

    COMPILE_OPT IDL2
    COMMON scan_data

    CI = project.ci
    dim  =  project.imndArray[CI].dimensions
    fdim_shift = project.imndArray[CI].fdim_shift
    pdim_shift = -project.imndArray[CI].pdim_shift
    sdim_shift =  project.imndArray[CI].sdim_shift
       
    if project.procPramArray[CI].zpad_flag eq 1 then begin
      fdim_shift*=2
      pdim_shift*=2
      if dim eq 3 then sdim_shift*=2
    endif 
 
    ;; In the event of a user-defined kspace span for resolution reduction, this
    ;; is a correction to the read-in phase shifts to take into account the
    ;; lower resolution.
    rshift_mult = float(project.imndarray[CI].k_fdim_span)/project.imndarray[ci].k_fdim_span_max
    pshift_mult = float(project.imndarray[CI].k_pdim_span)/project.imndarray[ci].k_pdim_span_max
    sshift_mult = float(project.imndarray[CI].k_sdim_span)/project.imndarray[ci].k_sdim_span_max

    fdim_shift *= rshift_mult
    pdim_shift *= pshift_mult
    if (project.imndarray[ci].dimensions eq 3) then sdim_shift *= sshift_mult
    
    sz_data = size(*p_data)
    fdim = sz_data[1]
    pdim = sz_data[2]
    sdim = sz_data[0] ge 3 ? sz_data[3] : 1
    adim = sz_data[0] ge 4 ? sz_data[4] : 1
    sz_data = size(*p_data)
    
    ;;shift the data according to it's dimensions
    ;;if the data is 2 dimensional then shift the data.
    if dim eq 2 then begin

        IF KEYWORD_SET(ISDICOM) NE 1 THEN BEGIN
            case sz_data[0] of
                4: begin;;(*p_DATA) = temporary(shift((*p_DATA),fdim/2,pdim/2,0,0))
                    for aa = 0, adim-1 do begin
                        for ss = 0, sdim-1 do begin
                            (*p_DATA)[*,*,ss,aa] = shift((*p_DATA)[*,*,ss,aa],fdim/2, pdim/2)
                        endfor
                    endfor
                end
                3: begin;;(*p_DATA) = temporary(shift((*p_DATA),fdim/2,pdim/2,0))
                    for ss = 0, sdim-1 do begin
                        (*p_DATA)[*,*,ss,0] = shift((*p_DATA)[*,*,ss,0],fdim/2, pdim/2)
                    endfor
                end
                2: (*p_DATA) = temporary(shift((*p_DATA),fdim/2,pdim/2))
            endcase
        ENDIF
        
        ;; Apply spatial shift in all dimensions
        If sz_data[0] EQ 4 then begin ;;(*p_DATA) = temporary(shift((*p_DATA),fdim_shift,pdim_shift,sdim_shift,0))
            for aa = 0, adim-1 do begin
                if (fdim_shift ne 0 or pdim_shift ne 0) then begin
                    for ss = 0, sdim-1 do begin
                        (*p_DATA)[*,*,ss,aa] = shift((*p_DATA)[*,*,ss,aa], fdim_shift, pdim_shift)
                    endfor
                endif
                if (sdim_shift ne 0) then begin
                    for rr = 0, fdim-1 do begin
                        for pp = 0, pdim-1 do begin
                            (*p_DATA)[rr, pp, *, aa] = shift((*p_DATA)[rr,pp,*,aa], sdim_shift)
                        endfor
                    endfor
                endif
            endfor
            
        endif

 ;       If sz_data[0] EQ 4 then begin
 ;           (*p_DATA) = temporary(shift((*p_DATA),fdim_shift,pdim_shift,sdim_shift, 0))
 ;       endif
        
        If sz_data[0] EQ 3 then begin
            (*p_DATA) = temporary(shift((*p_DATA),fdim_shift,pdim_shift,sdim_shift))
        endif
        
        If sz_data[0] EQ 2 then begin
            (*p_DATA) = temporary(shift((*p_DATA),fdim_shift,pdim_shift))
        endif

    endif

    ;;if the data is 3 dimensional then we have to shift this data differently
    if dim eq 3 then begin

        update_status_bar,'Shifting'
;        progressbar = Obj_New('progressbar', Color='red', Text='Shift FFT correction',/NOCANCEL)
;        progressbar -> Start
        for aa = 0, adim-1 do (*p_DATA)[*,*,*,aa] = shift((*p_DATA)[*,*,*,aa],fdim/2,pdim/2,0)
    
;        progressbar -> Destroy

        ;; Apply spatial shift in phase encode direction

       ;Phase_shift = Fix[RND{(IMND_slicepack_position/(sfov*10))*sdim)}]
       ;Slice_shift = Fix[RND{((IMND_phase1_offset/(pfov*10))+0.5)*pdim)}]
       ;
       ;Where does this come from?
       ;
       ;1) In the MAS program (3D version), phase and slice are reversed with
       ;respect to their definition in Paravision. Ultimately, this reversal doesn't
       ;matter because both directions are phase encoded in the acquisition. The
       ;effect is the same. But you may want to change it to correspond directly
       ;with Bruker assignments. This is the reason the slice offsets are used to
       ;calculate the phase shifts, etc.
       ;2) The method of phase encode acquisition defines the formula used. As
       ;defined in MAS, the PHASE direction is acquired using Paravision slice
       ;encoding. The PE of the Paravision slice direction is linear and sequential
       ;starting at 1 and progressing by integers to n, where n is the number of
       ;slices. Therefore, the middle of array is defined as the correct RECO_rotate
       ;parameter multiplied by the Paravision slice direction. Eg.: Given a 3D
       ;acquisition of 256x64x64 over a 8x2x2 mm FOV with offsets 0.5x0.05x-1.5mm,
       ;the RECO_rotate parameter in the slice direction is RECO_rotate(s)=(slice
       ;offset)/sFOV = -1.5/2 = -0.75 and the MAS PHASE_SHIFT = RECO_rotate(s)*sdim
       ;= -0.75*64 = -48.
       ;3) Again, the method of phase encode acquisition defines the formula used.
       ;As defined in MAS, the SLICE direction is acquired using Paravision phase
       ;encoding. The PE of the Paravision slice direction is linear and interlaced
       ;starting at -1 and progressing by equally spaced fractions to approximately
       ;1. Therefore, the middle of array is defined as (the correct RECO_rotate
       ;parameter plus 0.5 to find the center of the array) multiplied by the
       ;Paravision phase direction. Eg.: Given a 3D acquisition of 256x64x64 over a
       ;8x2x2 mm FOV with offsets 0.5x0.05x-1.5mm, the RECO_rotate parameter in the
       ;phase direction is RECO_rotate(p)=(phase1 offset)/pFOV = -0.05/2 = 0.025 and
       ;the MAS PHASE_SHIFT = (RECO_rotate(p)+0.5)*pdim = (0.025+0.5)*64 = -34.
       ;4) Paravision rounds it's numbers to the nearest array/pixel value. I don't
       ;know the IDL command but have inserted it as RND.

;      slencode_shift = Fix((sdim_shift/(sfov*10))*sdim)
;      phencode_shift = Fix(((pdim_shift/(pfov*10))+0.5)*pdim)

        for aa = 0, adim-1 do (*p_DATA)[*,*,*,aa] = shift((*p_DATA)[*,*,*,aa],fdim_shift,pdim_shift,sdim_shift)
;        progressbar = Obj_New('progressbar', Color='red', Text='Shift Acquisition correction 2/2',/NOCANCEL)
;        progressbar -> Start
;        for aa = 0, adim-1 do begin
;            for ii=0, fdim-1 do begin
;                for jj=0, pdim-1 do begin
;                    (*p_DATA)[ii,jj,*,aa] = shift((*p_DATA)[ii,jj,*,aa],sdim_shift)
;                end
;                if ii mod 10 eq 0 then progressBar->Update, (float(ii)/float(fdim-1))*100.0
;            endfor
;        endfor
;        progressbar -> Destroy

    endif

    update_status_bar,''

end


; Subroutine name: mas_signal_selection
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:
; This procedure will allow the user to specify what type of signal to pass on.

; Editing Information:
    ;Edited by HS 2006/10/02
    ;Fix spelling mistakes and commenting.
    ;Edited by GA 2010/01/13 - Fixed Phase image display
    ;IDL's arctan function (atan) results in the arctan mapped from -pi/2 to pi/2 if a ratio is supplied as an argument.
    ;The arctan result is mapped from -pi to pi if the imaginary and real parts are supplied as separate entities rather than
    ;a ratio. MR phase images should be mapped from -pi to pi. I simply replaced the division symbol with a comma to fix this error. 
    ; See the IDL help on the atan function for more information. 

pro mas_signal_selection , p_data
    COMPILE_OPT IDL2
    COMMON scan_data

    CI = project.ci
    
    dim  = project.imndArray[CI].dimensions
    fdim = project.imndArray[CI].fdim
    pdim = project.imndArray[CI].pdim
    sdim = project.imndArray[CI].sdim
    adim = project.imndArray[CI].adim
    if project.procPramArray[CI].zpad_flag eq 1 then begin
      fdim*=2
      pdim*=2
      if dim eq 3 then sdim*=2
    endif 
    
    signal_type = project.procPramArray[CI].signal_type

    update_status_bar,'Signal Selection'
    progressbar = Obj_New('progressbar', Color='red', Text='Signal Selection',/NOCANCEL)
    progressbar -> Start

	; allow for magnitude/real/imaginary selection of fft or non-fft(raw)data
	; 0 = fft magnitude		5 = raw magnitude
	; 1 = fft real			6 = raw real
	; 2 = fft imaginary		7 = raw imaginary
	; 3 = fft intensity     8 = fixed, magnitude -- already reconstructed.
	; 4 = fft phase          9 = fft complex

    if signal_type eq 0 or signal_type eq 5 then begin

                                ; magnatude part
        (*p_data) = temporary(abs(*p_data))

    endif else if signal_type eq 1 or signal_type eq 6 then  begin

                                ; real part
        (*p_DATA) = REAL_PART((*p_DATA))
        ;Need to account for FFT starting at 0 when kspace starts at -k 
;        if signal_type eq 1 then begin  
;          IF dim EQ 1 THEN BEGIN
;             
;             for i = 0,fdim-1 do begin
;                if i mod 2 ne 0 then begin
;                  (*p_data)[i] = -(*p_data)[i]
;                endif  
;             endfor
;          
;          ENDIF ELSE IF dim EQ 2 THEN BEGIN   
;              
;              for i=0,fdim -1 do begin
;                for j=0,pdim-1 do begin
;                  if i mod 2 ne 0 and j mod 2 eq 0 then begin
;                    (*p_data)[i,j] = -(*p_data)[i,j]
;                  endif
;                  if i mod 2 eq 0 and j mod 2 ne 0 then begin
;                    (*p_data)[i,j] = -(*p_data)[i,j]
;                  endif
;                endfor
;              endfor
;          
;          ENDIF ELSE IF dim EQ 3 THEN BEGIN
;          
;            for k=0, sdim-1 do begin  
;               if k mod 2 ne 0 then begin
;                  (*p_data)[*,*,k] = -(*p_data)[*,*,k]
;                endif 
;              for i=0,fdim -1 do begin
;                for j=0,pdim-1 do begin
;                  if i mod 2 ne 0 and j mod 2 eq 0 then begin
;                    (*p_data)[i,j] = -(*p_data)[i,j,k]
;                  endif
;                  if i mod 2 eq 0 and j mod 2 ne 0 then begin
;                    (*p_data)[i,j] = -(*p_data)[i,j,k]
;                  endif
;                endfor
;              endfor
;            endfor
;          
;          ENDIF
;         endif 

    endif else if signal_type eq 2 or signal_type eq 7 then begin

                                ; imag part
        (*p_DATA) = IMAGINARY((*p_DATA))
     
    endif else if signal_type eq 3 then begin

                                ; intensity type display
        (*p_DATA) = REAL_PART((*p_DATA)) / cos(atan(IMAGINARY((*p_DATA)) / REAL_PART((*p_DATA))))
        (*p_DATA) = abs((*p_DATA))

    endif else if signal_type eq 4 then begin

                                ; phase type display, edited by GA
        ;(*p_DATA) = atan(abs(IMAGINARY((*p_DATA))) , abs(REAL_PART((*p_DATA))))
    
        if project.procPramArray[CI].phi_unwrap_flag eq 1 then (*p_DATA) = mas_phase_unwrap(project.imndarray[CI].dimensions,(*p_DATA)) $
          else (*p_DATA) = atan((*p_DATA), /PHASE)
            
    endif else if signal_type eq 8 then begin
    
        ;; Nothing to do .. this implies already reconstructed imagery, as in
        ;; DICOM or NIFTI or PAR/REC, etc
        
    endif
    
    progressbar -> Destroy
    update_status_bar,''

end

;************************************

; Subroutine name: mas_smoothing
; Created by: AMV
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:
; This procedure will actually do all the smoothing of the images.

; Editing Information:
     ;Edited by AMV 2006/12/05
     ;Wrote the program

pro mas_smoothing , p_data
    common scan_data

    CI = project.ci

    smooth_Direction=project.procPramArray[CI].smooth_Direction

    sz_data = size(*p_data) ;Transfering the siza of the pointer p_data to sz_data

    if smooth_Direction eq 0 then begin ;smooth_Direction is equal zero if the user does not want to smooth the image
        return                  ;then end the program
    end

    if smooth_Direction eq 1 then begin ;if smooth_Direction is equal one the user wants to smooth the image using the
                                          ;smooth function

        if sz_data[0] eq 2 then begin ;check what type of image it is, depending on sz_data[0], then use the smooth function
;            scan = fltarr(sz_data[1],sz_data[2]) ;Creating an array
;            scan = *p_data     ;Transfering the data in p_data to scan
            (*p_data)[*,*] = smooth((*p_data)[*,*], 127 , /EDGE_TRUNCATE ,/NAN) ;Use the smooth function
        end


        if sz_data[0] eq 3 then begin ;here the image type is 3
;           scan = fltarr(sz_data[1],sz_data[2],sz_data[3]) ;Creating an array
;           scan = *p_data     ;Transfering the data in p_data to scan
            for i=0,sz_data[3]-1  DO BEGIN ;using a for statement to wary the value of sdim
                (*p_data)[*,*,i] = smooth((*p_data)[*,*,i], 3 , /EDGE_TRUNCATE ,/NAN) ;Use the smooth function
            ENDfor
        end

        if sz_data[0] eq 4 then begin ;here the image type is 4
;            scan = fltarr(sz_data[1],sz_data[2],sz_data[3],sz_data[4]) ;Creating an array
;            scan = *p_data     ;Transfering the data in p_data to scan
            for i=0,sz_data[3]-1  DO BEGIN ;using a for statement to wary the value of sdim
                for j=0,sz_data[4]-1 DO BEGIN ;using a for statement to wary the value of n_scans
                    (*p_data)[*,*,i,j] = smooth((*p_data)[*,*,i,j], 3 , /EDGE_TRUNCATE ,/NAN) ;Use the smooth function
                ENDfor
            ENDfor
        end

    end

    if smooth_Direction eq 2 then begin ;if smooth_Direction is equal two the user wants to smooth the image using the
                                ;median function
        if sz_data[0] eq 2 then begin ;check what type of image it is, depending on sz_data[0], then use the median function
;            scan = fltarr(sz_data[1],sz_data[2]) ;Creating an array
;            scan = *p_data     ;Transfering the data in p_data to scan
            (*p_data)[*,*] = median((*p_data)[*,*], 127) ;Use the median function
        end

        if sz_data[0] eq 3 then begin ;the image type is 3
;            scan = fltarr(sz_data[1],sz_data[2],sz_data[3]) ;Creating an array
;            scan = *p_data     ;Transfering the data in p_data to scan
            for i=0,sz_data[3]-1  DO BEGIN ;use a for statement to wary the value of sdim
                (*p_data)[*,*,i] = median((*p_data)[*,*,i], 3) ;Use the median function
            ENDfor
        end

        if sz_data[0] eq 4 then begin ;the image type is 4
;            scan = fltarr(sz_data[1],sz_data[2],sz_data[3],sz_data[4]) ;Creating an array
;            scan = *p_data     ;Transfering the data in p_data to scan
            for i=0,sz_data[3]-1  DO BEGIN ;use a for statement to wary the value of sdim
                for j=0,sz_data[4]-1 DO BEGIN ;use a for statement to wary the value of n_scans
                    (*p_data)[*,*,i,j] = median((*p_data)[*,*,i,j], 3) ;Use the median function
                ENDfor
            ENDfor
        endif
    endif

end


;***********************************

; Subroutine name: mas_phase_array
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:
; This procedure will actually do all the combining of the images using square root of sum of the square algorithm.

; Editing Information:
    ;Edited by HS 2006/10/02
    ;Fix spelling mistakes and commenting.

    ;Edited by AMV 2006/10/15
    ;Wrote the algorithm for the combining.

function  mas_phase_array, p_data
    common scan_data

    sz_data = size(*p_data)  ;Transfering the siza of the pointer p_data to sz_data

    ;first make an array to hold all the data.

    rms_data2 = ptr_new(fltarr(sz_data[1],sz_data[2],sz_data[3],sz_data[4]))
    ;temp = fltarr(sz_data[1],sz_data[1])

    image = fltarr(sz_data[1],sz_data[2],sz_data[3],sz_data[4])  ;Creating an array

    image = *p_data  ;Transfering the data in p_data to image



  ;Creating three for statements, because the pixels in x,y and z coordinate must multiple element by element
  ;(using square root of sum of the square algorithm) to get proper data

  for xx=0, sz_data[1]-1 DO BEGIN
    for yy=0, sz_data[2]-1 DO BEGIN
       for zz=0, sz_data[3]-1 DO BEGIN

         (*rms_data2)[xx,yy,zz,0] = sqrt((image[xx,yy,zz,0]^2)+(image[xx,yy,zz,1]^2)+(image[xx,yy,zz,2]^2)+(image[xx,yy,zz,3]^2))

       ENDfor
    ENDfor
  ENDfor



;    for ii=0, sz_data[3]-1 do begin
;       for jj=0, sz_data[4]-1 do begin
;
;          temp += ((*data)[*,*,ii,jj])^2
;       end
;       (*rms_data)[*,*,ii] = temp
;       ;(*rms_data)[*,*,ii] = ((*data)[*,*,ii,0])^2+((*data)[*,*,ii,1])^2 +((*data)[*,*,ii,2])^2 +((*data)[*,*,ii,3])^2
;
;       (*rms_data)[*,*,ii] = sqrt((*rms_data)[*,*,ii])
;
;    end

;    data = rms_data
;    ptr_free, rms_data

    project.imndArray[project.ci].adim = 1

    rms_data = ptr_new(*rms_data2, /no_copy)
    ;obj_destroy , image[*,*,*,*]
    heap_free , image
    ptr_free , rms_data2
    return, rms_data

end

pro mas_load_state_1_flt

    common scan_data, project
    
    mas_shift, project.dataArray[CI].state1
    mas_signal_selection , project.dataArray[CI].state1

end

pro mas_load_state_1_bruker_onepulse

    common scan_data, project
    forward_function mas_extract_data
    
    ci = project.ci
    project.dataArray[CI].state1 = mas_extract_data()

    if (ptr_valid(project.dataarray[ci].state1)) then begin
        project.procpramarray[CI].state_1 = 1
    endif
    
end

;***********************************

; Subroutine name: mas_load_state_1_bruker
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:
; This procedure will process state1 raw data

; Editing Information:
; Edited by Magdoom 2015/04/12
; Incorporated zeropadding

pro mas_load_state_1_bruker

    common scan_data, project
    forward_function mas_extract_data
    
    ci = project.ci
    
    project.dataArray[CI].state1 = mas_extract_data()
    
    if project.procPramArray[ci].zpad_flag eq 1 then mas_zeropad, project.dataarray[CI].state1
    
    mas_kspace_subsample, project.dataarray[ci].state1
    
    mas_kspace_shift, project.dataArray[CI].state1
    
    IF (project.procPramArray[ci].filterflag eq 'Frequency') THEN BEGIN
        frequency_filter
    ENDIF
    
    signal_type = project.procPramArray[CI].signal_type
    
    ;; signal types 5,6,7 are raw (non-fft) data
    ;; see mas_signal_selection for signal_type definitions
    
    IF signal_type LT 5 OR signal_type GT 7 THEN BEGIN
    
        ;; data that needs FFT, if motion correct, apply
        ;; now
        mas_fft,   project.dataArray[CI].state1
        mas_shift, project.dataArray[CI].state1
        
    ENDIF
    
    mas_signal_selection, project.dataArray[CI].state1
    
    IF signal_type LT 5 OR signal_type GT 7 THEN BEGIN
    
        ;; slice-direction part for 2D multislice "kspace" span reduction.
        mas_kspace_subsample_multislice, project.dataarray[CI].state1
    
        IF project.procPramArray[CI].mc_enable NE 0 THEN BEGIN
            mas_motion_correct
        ENDIF
        
        mas_smoothing , project.dataArray[CI].state1
        
        ;; HS August 10, 2007
        ;; Image domain filtering
        
        IF (project.procPramArray[ci].image_domain_filterflag eq 'Image Domain') THEN BEGIN
            image_domain_filter
        ENDIF
        
    ENDIF
    
    ;;if the scan is multi phase array then join them.
    if project.imndArray[ci].image_type eq 4 then begin
        print, 'RMS phase array'
        project.dataArray[CI].state1 = mas_phase_array( project.dataArray[CI].state1 )
    endif
    
    project.dataArray[CI].state1 = temporary(mas_interpolate( project.dataArray[CI].state1 ))

    project.procpramarray[CI].state_1 = 1
    
end

; Subroutine name: mas_load_state_1
; Created by: Unknown
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:
; Prepare to load state 1 and dispatch loading task to the proper procedure
; for whatever data source (Varian, Bruker, Philips, etc.)

; Editing Information:
; Magdoom (05/14/2015)
; - Added an optional keyword to choose the project index to load

pro mas_load_state_1, ci = ci

    common scan_data, project
    
    if ~keyword_set(ci) then ci = project.ci
    
    if (project.procpramarray[ci].state_1 eq 1) then return
    
    if (project.procpramarray[ci].have_orientation eq 0) then begin
        orient_file = project.imndarray[ci].file_path+'.orientation.dat'
        if (file_test(orient_file, /read)) then begin
            orientation = intarr(3)
            openr, lun, orient_file, /get_lun
            readf, lun, orientation
            close, lun
            free_lun, lun
            yesno = dialog_message(['Orientation hint file found: '+ $
                strjoin($
                    string($
                        mas_orient_from_matrix(mas_orient_get_matrix(orientation), /lettercode))), $
                'Do you want to apply it?'], /question, /center)
            if (yesno eq 'Yes') then begin
                project.procpramarray[ci].orientation_Imin = orientation[0]
                project.procpramarray[ci].orientation_Jmin = orientation[1]
                project.procpramarray[ci].orientation_Kmin = orientation[2]
                project.procpramarray[ci].have_orientation = 1
            endif
        endif
    endif
    
    ;; adt should be refit and leftover adt data pointers should
    ;; be free'd since the data may not be valid any longer
    mas_free_data_pointers
    
    project.imndarray[ci].fdim = max([project.imndarray[ci].fdim, project.imndarray[ci].k_fdim_span_max])
    project.imndarray[ci].pdim = max([project.imndarray[ci].pdim, project.imndarray[ci].k_pdim_span_max])
    project.imndarray[ci].sdim = max([project.imndarray[ci].sdim, project.imndarray[ci].k_sdim_span_max])
    
    project.procpramarray[ci].adt_proccess_flag = 0
    
    state1_loader = project.imndarray[ci].state1_load_procedure
    
    call_procedure, state1_loader

    if (project.procpramarray[ci].state_1 ne 1) then begin
        update_status_bar, "An error occurred while loading data."
    endif else begin
        update_status_bar , ' '
    
        ;;after processing state 1 set the flag back to true
        ;;project.procPramArray[ci].state_1 = 1
    
        ;;set the threshold and max threshold
        state1_min = min((*project.dataArray[CI].state1),   max = state1_max )
        project.procPramArray[ci].threshold = [state1_min,state1_max]
        project.procPramArray[ci].min_display_thresh = state1_min
        project.procPramArray[ci].max_display_thresh = state1_max
   endelse
    
   mas_redraw_gui
    
   HEAP_GC

end

pro mas_free_data_pointers

    common scan_data, project
    
    ci = project.ci
    
    if (ptr_valid(project.dataArray[ci].adt)) then ptr_free, project.dataArray[ci].adt
    if (ptr_valid(project.dataArray[ci].eign_val)) then ptr_free, project.dataArray[ci].eign_val
    if (ptr_valid(project.dataArray[ci].eign_Vec)) then ptr_free, project.dataArray[ci].eign_Vec
    if (ptr_valid(project.dataArray[ci].frac_Ani)) then ptr_free, project.dataArray[ci].frac_Ani
    if (ptr_valid(project.dataArray[ci].Avg_Dif)) then ptr_free, project.dataArray[ci].Avg_Dif
    
    if (ptr_valid(project.dataArray[ci].fibers_adt)) then ptr_free, project.dataArray[ci].fibers_adt
    if (ptr_valid(project.dataArray[ci].fibers_eign_val)) then ptr_free, project.dataArray[ci].fibers_eign_val
    if (ptr_valid(project.dataArray[ci].fibers_eign_Vec)) then ptr_free, project.dataArray[ci].fibers_eign_Vec
    if (ptr_valid(project.dataArray[ci].fibers_frac_Ani)) then ptr_free, project.dataArray[ci].fibers_frac_Ani
    if (ptr_valid(project.dataArray[ci].fibers_Avg_Dif)) then ptr_free, project.dataArray[ci].fibers_Avg_Dif
    
    ;; get rid of the old data, if any
    if (ptr_valid(project.dataArray[ci].state1)) then ptr_free, project.dataarray[ci].state1
    if (ptr_valid(project.dataArray[ci].state2)) then ptr_free, project.dataarray[ci].state2
    
    heap_gc
    
end

;*****************************************************************************************************
;
; NAME:
;   mas_load_state_1
;
; PURPOSE:
;   load images from files into memory for processing.
;   then fft the data if necessary
;   then interpolate if necessary
;   then filter which signal type the user wants
;   then store the data into state1
;
; ARGUMENTS:
;
;
; MODIFICATION HISTORY:
;   20080701 - Now frees any existing data pointers before reloading
;              Support for rudimentary NIFTI reader to overload state1 data
;
;
;
;*****************************************************************************************************
pro mas_load_state_1_orig
    COMPILE_OPT IDL2
    common scan_data
    
    CI = project.ci
    
    if project.procPramArray[ci].state_1 ne 0 then return

    ;; if the image doesn't need to be extracted then don't. ie dicom images
    if project.imndArray[ci].image_type ne 5 $
        and project.imndArray[ci].image_type ne 6 $
        and project.imndarray[ci].image_type ne 16 $
        and project.imndArray[ci].image_type ne 15 $
        and project.imndArray[ci].image_type ne 18 then begin
        
        ;; if the image type is 11 (dicom) then check to see if extraction
        ;; necessary, skip all other steps below.
        IF project.imndArray[ci].image_type EQ 11 THEN BEGIN

            ;; PAR/REC and Old Philips DICOM went here.
            
        ENDIF ELSE BEGIN

            ;;; Bruker went here
            
        ENDELSE
        
    end else if project.imndArray[ci].image_type eq 15 then begin
    
        ;; "flt" file state 1 loading went here
        
    endif else if project.imndarray[ci].image_type eq 16 then begin
    
        ;; Nifti state1 loading went here
        
    endif else if project.imndarray[ci].image_type eq 18 then begin
    
        ;; generic DICOM loading went here
        
    endif
    
;    update_status_bar , ' '
;    
;    ;;after processing state 1 set the flag back to true
;    project.procPramArray[ci].state_1 = 1
;    
;    ;;set the threshold and max threshold
;    state1_min = min((*project.dataArray[CI].state1),   max = state1_max )
;    project.procPramArray[ci].threshold = [state1_min,state1_max]
;    
;    
;    mas_redraw_gui
;    
;    HEAP_GC
    
end
