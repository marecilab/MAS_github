;; $Id$
;;
; HS 20070805
; Collection of filtering techniques
; Called by: MAS.PRO

; =============================== Custom Windows ================================

; NAME:
;   kaiser_bessel
;
; PURPOSE:
;   Window function for Fourier Transform filtering.
;
; CATEGORY:
;   Signal, image processing.
;
; CALLING SEQUENCE:
;
;   Result = kaiser_bessel(N1, N2, alpha)
;
; INPUTS:
;   N1: The number of columns of the result.
;
;   N2: The number of rows of the result.
;
; Keyword Parameters:
;
; Alpha is the steepness of the filter
;
;
; OUTPUTS:
;   Result(i) = BESELI(pi*alpha*sqrt(1-((2*index)/n1-1)^2),0)/BESELI(pi*alpha,0)
;
;   For two dimensions, the result is the same except that "i" is replaced
;   with "i*j", where i and j are the row and column subscripts.

function kaiser_bessel, n1In, n2In, alpha

    compile_opt IDL2

    on_error,2                              ;Return to caller if an error occurs
    n1 = n1In
    n2 = n2In
    one = 1
    pi = !PI

    index = FINDGEN(n1)

	temp = pi*alpha*sqrt(1-((2*index)/n1-1)^2)

    row = BESELI(temp,0)/BESELI(pi*alpha,0)


    index = FINDGEN(n2)

	temp = pi*alpha*sqrt(1-((2*index)/n1-1)^2)

    col = BESELI(temp,0)/BESELI(pi*alpha,0)

	; Generates a symmetrical matrix with the shape of the filter
    RETURN,(row # col)

 end

; Blackman window
function blackman, n1In, n2In, alphazero, alphaone

    compile_opt IDL2

    on_error,2                              ;Return to caller if an error occurs
    n1 = n1In
    n2 = n2In
    one = 1
    pi = !PI
    alpha0 = alphazero
	alpha1 = alphaone
	alpha2 = 0.08

    index = FINDGEN(n1)

	row = alpha0 - alpha1*cos(2*pi*index/(n1-1)) + alpha2*cos(4*pi*index/(n1-1))

    index = FINDGEN(n2)

    col = alpha0 - alpha1*cos(2*pi*index/(n1-1)) + alpha2*cos(4*pi*index/(n1-1))

	RETURN,(row # col)

 end




pro frequency_filter
    COMPILE_OPT IDL2, HIDDEN
    COMMON scan_data
    CI = project.ci

    print, 'Inside of frequency filter'
    
    p_DATA = *project.dataArray[CI].state1
    
    if project.procPramArray[ci].zpad_flag eq 1 then begin
      fdim = 2*project.imndArray[CI].fdim
      pdim = 2*project.imndArray[CI].pdim
    endif else begin
      fdim = project.imndArray[CI].fdim
      pdim = project.imndArray[CI].pdim
    endelse
      
    sdim = project.imndArray[CI].sdim
    adim = project.imndArray[CI].adim
    
    subtype = project.procPramArray[ci].subtype
    design  = project.procPramArray[ci].design
    fstop   = project.procPramArray[ci].fstop
    fstart  = project.procPramArray[ci].fstart
    order   = project.procPramArray[ci].order


    if (subtype eq 'Lowpass') then begin
        
        if (design eq 'Butterworth') then begin
            
            if (fstop eq 0) then fstop = long(fdim/2)
            
            print, 'Butterworth Lowpass Filter'
            filter = 1.0 / ( 1.0d + (DIST(fdim,pdim)/fstop)^2*order )
            Display_DrawResS_Create, bytscl(filter) , 'Butterworth Filter', fdim, pdim

        endif else if (design eq 'Hanning') then begin
            
            print, 'Hanning Window'
            filter = HANNING(fdim,pdim, ALPHA=fstop)
            Display_DrawResS_Create, bytscl(filter) , 'Hanning Window Filter', fdim, pdim
            
        endif else if (design eq 'Hamming') then begin
            
            print, 'Hamming Window'
            ;;Change equation
            filter = HANNING(fdim,pdim,  ALPHA=0.56)
            Display_DrawResS_Create, bytscl(filter) , 'Hamming Window Filter', fdim, pdim
            
            
        endif else if (design eq 'Kaiser-Bessel' ) then begin
            if (fstop eq 0) then fstop = 2
            print, 'Kaiser-Bessel Window'
            filter = kaiser_bessel( fdim,pdim ,fstop)
            Display_DrawResS_Create, bytscl(filter) , 'Kaiser-Bessel Window Filter', fdim, pdim
            
        endif else if (design eq 'Blackman' ) then begin
            if (fstop eq 0) then fstop = 0.42
            if (order eq 0) then order = 0.5
            print, 'Blackman Window'
            filter = blackman( fdim,pdim ,fstop,order)
            Display_DrawResS_Create, bytscl(filter) , 'Kaiser-Bessel Window Filter', fdim, pdim
        endif
        
    endif else if (subtype eq 'Highpass') then begin
        
        if (design eq 'Butterworth') then begin
            
            print, 'Butterworth Highpass Filter'
            filter = 1.0 / ( 1.0d + (fstart/DIST(fdim,pdim))^2*order )
            Display_DrawResS_Create, bytscl(filter) , 'Butterworth Filter', fdim, pdim
            

        endif else if (design eq 'Distance') then begin
            
            print, 'Distance Window Filter'
                                ;Change equation
            filter = DIST(fdim,pdim)
            Display_DrawResS_Create, bytscl(filter) , 'Distance Window Filter', fdim, pdim
            
            
        endif
        
    endif


    FOR i=0,sdim-1 DO BEGIN
        
        FOR j=0,adim-1 do begin

            p_DATA[*,*,i,j] = temporary(p_DATA[*,*,i,j]) * filter

        ENDFOR
        
    ENDFOR

    project.dataArray[CI].state1 = ptr_new(p_DATA)
    
end


PRO image_domain_filter

    COMPILE_OPT IDL2, HIDDEN
    COMMON scan_data
    CI = project.ci

    print, 'Inside of Image Domain Filter'
    
    p_DATA = *project.dataArray[CI].state1

    fdim = project.imndArray[CI].fdim
    pdim = project.imndArray[CI].pdim
    sdim = project.imndArray[CI].sdim
    adim = project.imndArray[CI].adim

    subtype   = project.procPramArray[ci].image_domain_subtype
    parameter = project.procPramArray[ci].image_domain_parameter

    cancel_status = 0
    
    if (subtype eq 'Lee Adaptive Filter') then begin
        
        progressbar = Obj_New('progressbar', Color='red', Text='Calculating Lee Filter',/fast_loop)
        progressbar -> Start
        if (adim eq 1) then begin
            FOR i=0,sdim-1 DO BEGIN
                p_DATA[*,*,i] = LEEFILT(p_DATA[*,*,i], parameter, /DOUBLE, /EXACT)
                progressBar -> Update, (float(i)/float(sdim-1))*100.0
                if (progressBar->checkCancel()) then begin
                    cancel_status = 1
                    break
                endif
            ENDFOR
        endif else begin
            iter = 0
            for j = 0, adim-1 do begin
                FOR i=0,sdim-1 DO BEGIN
                    iter++
                    p_DATA[*,*,i,j] = LEEFILT(p_DATA[*,*,i,j], parameter, /DOUBLE, /EXACT)
                    progressBar -> Update, (float(iter)/float((adim-1)*(sdim-1)))*100.0
                    if (progressBar->checkCancel()) then begin
                        cancel_status = 1
                        j = adim - 1
                        break
                    endif
                ENDFOR
            endfor
        endelse
        progressBar -> Destroy
        
    endif else if (subtype eq 'Wavelet Transform') then begin

        progressbar = Obj_New('progressbar', Color='red', Text='Calculating Wavelet Filter',/fast_loop)
        progressbar -> Start
        processed_data = fltarr(fdim,pdim,sdim)
        if (adim eq 1) then begin
            FOR i=0,sdim-1 DO BEGIN
                temp = p_DATA[*,*,i]
                temp2 = WV_DENOISE(temp, 'Daubechies', 2, PERCENT=parameter);, DENOISE_STATE=denoise_state)
                processed_data[*,*,i] = temp2[0:fdim-1, 0:pdim-1]
                p_data[*,*,i]  = processed_data[*,*,i]
                progressBar -> Update, (float(i)/float(sdim-1))*100.0
                if (progressBar->checkCancel()) then begin
                    cancel_status = 1
                    break
                endif
            ENDFOR
        endif else begin
            iter=0
            FOR j = 0, adim-1 do begin
                FOR i=0, sdim-1 DO BEGIN
                    iter++
                    temp = p_DATA[*,*,i,j]
                    temp2 = WV_DENOISE(reform(temp), 'Daubechies', 2, PERCENT=parameter);, DENOISE_STATE=denoise_state)
                    p_data[*,*,i,j] = temp2[0:fdim-1, 0:pdim-1]
                    progressBar -> Update, (float(iter)/float((adim-1)*(sdim-1))*100.0)
                    if (progressBar->checkCancel()) then begin
                        cancel_status = 1
                        j = adim - 1
                        break
                    endif
                ENDFOR
            ENDFOR
        endelse

        progressBar -> Destroy
;        print, 'Percent of power retained: ', denoise_state.percent

    endif else if (subtype eq 'Sharpening') then begin
        
        if (parameter eq 1) then begin
            kernelSize = [3, 3]
            kernel = REPLICATE(-1./9., kernelSize[0], kernelSize[1])
            kernel[1, 1] = 1.
        endif
        
        progressbar = Obj_New('progressbar', Color='red', Text='Sharpening',/fast_loop)
        progressbar -> Start
        
        temp = p_DATA

        if (adim eq 1) then begin
            FOR i=0,sdim-1 DO BEGIN
                p_DATA[*,*,i] = CONVOL(FLOAT(temporary(p_DATA[*,*,i])), kernel, /CENTER, /EDGE_TRUNCATE)
                progressBar -> Update, (float(i)/float(sdim-1))*100.0
                if (progressBar->checkCancel()) then begin
                    cancel_status = 1
                    break
                endif
            ENDFOR
        endif else begin
            iter = 0
            for j = 0, adim-1 do begin
                FOR i=0,sdim-1 DO BEGIN
                    iter++
                    p_DATA[*,*,i,j] = CONVOL(FLOAT(temporary(p_DATA[*,*,i,j])), kernel, /CENTER, /EDGE_TRUNCATE)
                    progressBar -> Update, (float(iter)/float((adim-1)*(sdim-1)))*100.0
                    if (progressBar->checkCancel()) then begin
                        cancel_status = 1
                        j = adim - 1
                        break
                    endif
                ENDFOR
            endfor
        endelse
        p_DATA = temporary(p_DATA)+temp
        progressBar -> Destroy
        
    endif else if (subtype eq 'Low Pass Filter') then begin
        
        if (parameter eq 1) then begin
            kernelSize = [3, 3]
            kernel = REPLICATE((1./(kernelSize[0]*kernelSize[1])), kernelSize[0], kernelSize[1])
        endif
        
        progressbar = Obj_New('progressbar', Color='red', Text='Low Pass Filter', /fast_loop)
        progressbar -> Start
        if (adim eq 1) then begin
            FOR i=0,sdim-1 DO BEGIN
                p_DATA[*,*,i] = CONVOL(FLOAT(temporary(p_DATA[*,*,i])), kernel, /CENTER, /EDGE_TRUNCATE)
                progressBar -> Update, (float(i)/float(sdim-1))*100.0
                if (progressBar->checkCancel()) then begin
                    cancel_status = 1
                    break
                endif
            ENDFOR
        endif else begin
            iter = 0
            for j = 0, adim-1 do begin
                FOR i=0,sdim-1 DO BEGIN
                    iter++
                    p_DATA[*,*,i,j] = CONVOL(FLOAT(temporary(p_DATA[*,*,i,j])), kernel, /CENTER, /EDGE_TRUNCATE)
                    progressBar -> Update, (float(iter)/float((adim-1)*(sdim-1)))*100.0
                    if (progressBar->checkCancel()) then begin
                        cancel_status = 1
                        j = adim - 1
                        break
                    endif
                ENDFOR
            ENDFOR
        endelse
        progressBar -> Destroy
        
        
    endif else if (subtype eq 'High Pass Filter') then begin
        
        if (parameter eq 1) then begin
            kernelSize = [3, 3]
            kernel = REPLICATE(-1./9., kernelSize[0], kernelSize[1])
            kernel[1, 1] = 8/9.
        endif
        
        progressbar = Obj_New('progressbar', Color='red', Text='Sharpening')
        progressbar -> Start
        
        if adim eq 1 then begin
            FOR i=0,sdim-1 DO BEGIN
                p_DATA[*,*,i] = CONVOL(FLOAT(temporary(p_DATA[*,*,i])), kernel, /CENTER, /EDGE_TRUNCATE)
                progressBar -> Update, (float(i)/float(sdim-1))*100.0
                if (progressBar->checkCancel()) then begin
                    cancel_status = 1
                    break
                endif
            ENDFOR
        endif else begin
            iter = 0
            for j = 0, adim-1 do begin
                FOR i=0,sdim-1 DO BEGIN
                    p_DATA[*,*,i,j] = CONVOL(FLOAT(temporary(p_DATA[*,*,i,j])), kernel, /CENTER, /EDGE_TRUNCATE)
                    progressBar -> Update, (float(iter)/float((adim-1)*(sdim-1)))*100.0
                    if (progressBar->checkCancel()) then begin
                        cancel_status = 1
                        j = adim - 1
                        break
                    endif
                ENDFOR
            endfor
        endelse
        progressBar -> Destroy
        
    endif
    
    if (cancel_status eq 0) then begin
        project.dataArray[CI].state1 = ptr_new(p_DATA)
    endif

    
END




;Empty routine to compile
PRO mas_image_filters
END
