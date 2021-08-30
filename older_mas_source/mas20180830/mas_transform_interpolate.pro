function mas_read_transformation_matrix, file, ask=ask

    common scan_data
    
    I4 = diag_matrix(replicate(1D,4))
    
    if (keyword_set(ask)) then begin
        file = dialog_pickfile(path=project.current_path, /read, /file)
        if (file eq '') then return, I4
    endif
    
    out = dblarr(4,4)
    openr, lun, file, /get_lun
    readf, lun, out
    close, lun & free_lun, lun
    
    return, out

end

; Subroutine name: mas_interpolate
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:
; This procedure will interpolate the data by the factors set by the users, and chages the size of r/p/sdim

; Editing Information:
    ;Edited by HS 2006/10/02
    ;Fix spelling mistakes and commenting.
    ;Downsampling functionality added by GA 2008/09/25

function mas_interpolate, p_data

    COMPILE_OPT IDL2
    COMMON scan_data

    ;return, mas_interpolate_fft(p_data)
    
    ;;bring in the data that is necessary.
    freq_interp  = project.procPramArray[project.ci].freq_interp
    Phase_interp = project.procPramArray[project.ci].Phase_interp
    Slice_interp = project.procPramArray[project.ci].Slice_interp
    dim =  project.imndArray[project.ci].dimensions

    ;; tracks user cancelling the progress bar. only if this is still
    ;; zero at the end will the data get overwritten
    cancel_check = 0

    ;;if the image doesn't need to be
    ;;interpolated then don't do any thing
    if freq_interp eq 1 and phase_interp eq 1 and slice_interp eq 1 then return, p_data

    sz_scan = size((*p_DATA))

    update_status_bar, 'Interpolating'
    ;set the new sizes for the interpolated and down sampling cases
    if project.procPramArray[project.ci].down_sample eq 0 then begin
        temp_sx = freq_interp*sz_scan[1]
        new_sx = freq_interp*sz_scan[1]
        new_sy = Phase_interp *sz_scan[2]
        new_sz = Slice_interp * sz_scan[3]
    endif else begin
        temp_sx = freq_interp*sz_scan[1] ;bc downsampling occurs in the y-z direction first we have to make the empty matrix a little bigger in the freq. direction initially
        new_sx = sz_scan[1]/freq_interp + (sz_scan[1] mod freq_interp)/freq_interp
        new_sy = sz_scan[2]/Phase_interp + (sz_scan[2] mod Phase_interp)/Phase_interp
        new_sz = sz_scan[3]/Slice_interp + (sz_scan[3] mod Slice_interp)/Slice_interp
    endelse

        ;; single slice / 2D
        if sz_scan[0] eq 2 then begin
            (*p_DATA) = temporary(CONGRID((*p_DATA),new_sx ,new_sy, CUBIC=-0.5 ))
        endif

        ;; 3D scan
        if sz_scan[0] eq 3 then begin
            ;;depending upon what type the data that is passed in is we have
            ;; to make an array that is of the same type
            if sz_scan[4] eq 4 then begin
                ;;data must be floating point
                scan = ptr_new(fltarr(temp_sx, new_sy, new_sz ))
            endif else begin
                ;;data must be complex
                scan = ptr_new(complexarr(temp_sx, new_sy, new_sz))
            end

            progressbar = Obj_New('progressbar', Color='red', Text='Interpolating 1/2', /fast_loop)
            progressbar -> Start

            for ii=0, sz_scan[1]-1 do begin

                (*scan)[ii,*,*] = CONGRID(reform((*p_DATA)[ii,*,*]) , new_sy, new_sz, CUBIC=-0.5)
                progressBar -> Update, (float(ii)/float(sz_scan[3]-1))*100.0
                if (progressBar->checkCancel()) then begin
                    cancel_check = 1
                    break
                endif
            endfor

            progressbar -> Destroy

            if (cancel_check eq 0) then begin

                progressbar = Obj_New('progressbar', Color='red', Text='Interpolating 2/2')
                progressbar -> Start

                for ii=0, new_sz-1 do begin
                    (*scan)[0:new_sx-1,*,ii] = CONGRID(reform((*scan)[0:sz_scan[1]-1,*,ii]), new_sx , new_sy, CUBIC=-0.5)
                    progressBar -> Update, (float(ii)/float(new_sz-1))*100.0
                    if (progressBar->checkCancel()) then begin
                        cancel_check = 1
                        break
                    end
                end

                progressbar -> Destroy

                if (cancel_check eq 0) then begin
                    p_data = ptr_new((*scan)[0:new_sx-1,*,*], /no_copy)
                endif

                ptr_free, scan

            endif

        endif

        ;; 4D scan
        if sz_scan[0] eq 4 then begin

            ;;depending upon what type the data that is passed in is we have
            ;; to make an array that is of the same type
            if sz_scan[5] eq 4 then begin
                ;;data must be floating point
                scan = ptr_new(fltarr(temp_sx, new_sy, new_sz, sz_scan[4]))

            endif else begin
                ;;data must be complex
                scan = ptr_new(complexarr(temp_sx, new_sy, new_sz, sz_scan[4]))
            endelse

            progressbar = Obj_New('progressbar', Color='red', Text='Interpolating')
            progressbar -> Start

            ;;traverse the adim interpolation as you go.
            ;;*****************************************************
            ;;begins from folder 0 to folder 3
            for aa=0, sz_scan[4]-1 do begin

                ;; first congrid the y,z plane then congrid the x,y plane
                for ii=0, sz_scan[1]-1 do begin
                    (*scan)[ii,*,*,aa] = CONGRID(reform((*p_DATA)[ii,*,*,aa]) , new_sy, new_sz, sz_scan[4], CUBIC=-0.5)
                endfor

                if (progressbar->checkCancel()) then begin
                    cancel_check = 1
                    break
                endif

                for ii=0, new_sz-1 do begin
                    (*scan)[0:new_sx-1,*,ii,aa] = CONGRID(reform((*scan)[0:sz_scan[1]-1,*,ii,aa]), new_sx , new_sy, sz_scan[4], CUBIC=-0.5)
                endfor

                progressBar -> Update, (float(aa)/float(sz_scan[4]-1))*100.0

            endfor

            progressbar -> Destroy

            if (cancel_check eq 0) then begin
                ptr_free, p_data
                p_data = scan; ptr_new((*scan)[0:new_sx-1,*,*,*], /no_copy)
            endif

            ;ptr_free, scan

        endif

    ;; handle cancel by resetting everything to one and redrawing
    ;; the GUI
    if (cancel_check eq 1) then begin
        project.procPramArray[project.ci].freq_interp  = 1
        project.procPramArray[project.ci].Phase_interp = 1
        project.procPramArray[project.ci].Slice_interp = 1
        mas_redraw_gui
    endif

    update_status_bar, ''

    return, p_DATA

end

function mas_interpolate_fft, p_data

    common scan_data
    
    ci = project.ci
        ;;bring in the data that is necessary.
    freq_interp  = project.procPramArray[ci].freq_interp
    Phase_interp = project.procPramArray[ci].Phase_interp
    Slice_interp = project.procPramArray[ci].Slice_interp
    dim =  project.imndArray[ci].dimensions
    
    adim = project.imndarray[ci].adim
    
    ;; array of interpolation factors
    itp = [freq_interp,Phase_interp,Slice_interp]
    
    if (itp[0] eq 1 and itp[1] eq 1 and itp[2] eq 1) then return, p_data
    
    state1 = p_data 
    
    sz = size(*state1, /dimensions)

    ;; data structures to hold stage one of the interpolation
    fft_work1 = complexarr(sz[0]*itp[0], sz[1]*itp[1])
    sz_fft1 = size(fft_work1, /dimensions)
    out1 = complexarr(sz[0]*itp[0], sz[1]*itp[1], sz[2])
    sz_out1 = size(out1, /dimensions)
    
    ;; data structures to hold stage two of the interp.
    fft_work2 = complexarr(sz[1]*itp[1], sz[2]*itp[2])
    sz_fft2 = size(fft_work2, /dimensions)

    ;; final data array holding interpolated p_data
    out_final = fltarr([sz[0:2] * itp, adim])
    
    pbar = obj_new('progressbar', text='FFT Interpolating...', title='Interpolating')
    pbar->Start
    
    for a = 0, adim-1 do begin
    
        ;; midpoints
        st = (sz_fft1 - sz)/2
        en = (sz_fft1 + sz)/2
        
        for zz = 0, sz[2]-1 do begin
            fft_work1[st[0]:en[0]-1, st[1]:en[1]-1] = shift(reform(fft((*state1)[*,*,zz,a])), sz[0]/itp[0], sz[1]/itp[1])
            out1[*,*,zz] = fft(fft_work1, /inverse)
        endfor
        
        ;; midpoints
        st = ((sz_fft2 - sz_out1[1:2])/2)[1]
        en = ((sz_fft2 + sz_out1[1:2])/2)[1]
        
        for xx = 0, sz_out1[0]-1 do begin
            fft_work2[*, st:en-1] = shift(reform(fft(out1[xx,*,*])), sz_out1[1]/itp[1], sz_out1[2]/itp[2])
            out_final[xx,*,*,a] = abs(fft(fft_work2, /inverse))
        endfor

        pbar->Update, float(a+1)/float(adim) * 100.0
        
    endfor

    pbar->Destroy
    return, ptr_new(reform(out_final), /no_copy)

end

; Subroutine name: mas_3d_data_average
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:
; interpolate 4d datasets between pixels

; Editing Information:
    ;Edited by BT 2008/09/02
    ;Fix spelling mistakes and commenting.

function mas_3d_data_average, coord, data_ptr=data_ptr

    common scan_data
    
    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        return, -1
    endif

    data = data_ptr

    if (total(fix(coord)-coord) eq 0) then begin
        return, (*data)[coord[0],coord[1],coord[2]]
    endif
    
    base = floor(coord)
    sum  = 0
    lengths = fltarr(3,3)
    weights = fltarr(3,3,3)

    for i = 0,2 do begin
        if (round(coord[i]) eq base[i]) then begin
            lengths[i,0] = 0.5 + float(base[i])-coord[i]
            lengths[i,1] = coord[i] + 0.5-float(base[i])
            lengths[i,2] = 0.0
        endif else begin
            lengths[i,0] = 0.0
            lengths[i,1] = 1.5 + float(base[i])-coord[i]
            lengths[i,2] = coord[i] - 0.5-float(base[i])
        endelse
    endfor
    
    for x=0,2 do begin
        for y = 0,2 do begin
            for z = 0,2 do begin
                
                weights[x,y,z] = lengths[0, x]*lengths[1,y]*lengths[2,z]
                
                aa = base[0]-1+x
                bb = base[1]-1+y
                cc = base[2]-1+z

                sum += weights[x,y,z]*(*data)[aa,bb,cc]

            endfor
        endfor
    endfor

    return, sum
end


; Subroutine name: mas_4d_data_average
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:
; interpolate 4d datasets between pixels

; Editing Information:
    ;Edited by BT 2008/09/02
    ;Fix spelling mistakes and commenting.

function mas_4d_data_average, coord, data_ptr=data_ptr

    common scan_data

    if (ptr_valid(data_ptr)) then begin
        data = data_ptr
    endif else begin
        data = project.dataarray[project.ci].state1
    endelse
    
    nvolumes = (size(*data, /dimensions))[3]

    catch, error_state
    if (error_state ne 0) then begin
        catch, /cancel
        return, fltarr(nvolumes)-1
    endif
    
    base = floor(coord)
    sum  = fltarr(nvolumes)
    lengths = fltarr(3,3)
    weights = fltarr(3,3,3)

    for i = 0,2 do begin
        if (round(coord[i]) eq base[i]) then begin
            lengths[i,0] = 0.5 + float(base[i])-coord[i]
            lengths[i,1] = coord[i] + 0.5 - float(base[i])
            lengths[i,2] = 0.0
        endif else begin
            lengths[i,0] = 0.0
            lengths[i,1] = 1.5 + float(base[i])-coord[i]
            lengths[i,2] = coord[i] - 0.5 - float(base[i])
        endelse
    endfor
    
    for x=0,2 do begin
        for y = 0,2 do begin
            for z = 0,2 do begin
                
                weights[x,y,z] = lengths[0,x]*lengths[1,y]*lengths[2,z]
                
                aa = base[0]-1+x
                bb = base[1]-1+y
                cc = base[2]-1+z

                sum += weights[x,y,z]*(*data)[aa,bb,cc,*]

            endfor
        endfor
    endfor

    return, sum
end

; Subroutine name: mas_tricubic_4d
; Created by: Bill Triplett
; Calling Information:
; data = mas_tricubic([x, y, z])

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; Returns pointwise interpolated data for sub-voxel
; points in a 4D data set (ie, diffusion)

; Editing Information:

function mas_tricubic_4d, coord, data_ptr=data_ptr

    common scan_data

    if (keyword_set(data_ptr)) then begin
        data = data_ptr
    endif else begin
        data = project.dataarray[project.ci].state1
    endelse
        
    data_sz = size(*data, /dimensions)
    nvolumes = n_elements(data_sz) eq 4 ? data_sz[3] : 1

    sum = fltarr(nvolumes)

    fc = float(floor(coord))

    if (fc[0] le 0 or fc[1] le 0 or fc[2] le 0 OR $
        fc[0] ge data_sz[0]-2 or fc[1] ge data_sz[1]-2 or fc[2] ge data_sz[2]-2) then begin 
        fc = (fc > 0) < (data_sz[0:2]-1)
        return, reform((*data)[fc[0], fc[1], fc[2], *])
    endif

    x = float(coord[0])
    y = float(coord[1])
    z = float(coord[2])
    
    dx = (x - fc[0]) & dx2 = dx^2 & dx3 = dx^3
    dy = (y - fc[1]) & dy2 = dy^2 & dy3 = dy^3
    dz = (z - fc[2]) & dz2 = dz^2 & dz3 = dz^3

    u_x = [ -0.5 * dx3 +       dx2 - 0.5 * dx, $
             1.5 * dx3 - 2.5 * dx2 + 1.0     , $
            -1.5 * dx3 + 2.0 * dx2 + 0.5 * dx, $
             0.5 * dx3 - 0.5 * dx2 ]

    u_y = [ -0.5 * dy3 +       dy2 - 0.5 * dy, $
             1.5 * dy3 - 2.5 * dy2 + 1.0     , $
            -1.5 * dy3 + 2.0 * dy2 + 0.5 * dy, $
             0.5 * dy3 - 0.5 * dy2 ]

    u_z = [ -0.5 * dz3 +       dz2 - 0.5 * dz, $
             1.5 * dz3 - 2.5 * dz2 + 1.0     , $
            -1.5 * dz3 + 2.0 * dz2 + 0.5 * dz, $
             0.5 * dz3 - 0.5 * dz2 ]

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    uu_x = fltarr(4,4,4)

    uu_x[*,*,0] = transpose([ [replicate(u_x[0], 4)], [replicate(u_x[1], 4)], [replicate(u_x[2], 4)], [replicate(u_x[3], 4)] ])
    uu_x[*,*,1] = uu_x[*,*,0]
    uu_x[*,*,2] = uu_x[*,*,0]
    uu_x[*,*,3] = uu_x[*,*,0]
        
    uu_y  = transpose([ [replicate(u_y[0], 4)], [replicate(u_y[1], 4)], [replicate(u_y[2], 4)], [replicate(u_y[3], 4)] ])
    
    sum = fltarr(nvolumes)
    for i = 0, nvolumes-1 do begin
        block = (*data)[fc[0]-1:fc[0]+2, fc[1]-1:fc[1]+2, fc[2]-1:fc[2]+2, i]
        sum[i] = total( u_z * total( uu_y * total( uu_x * block , 1 ), 1 ) ) > 0.0
    endfor
    
    return, sum
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
;    q = fltarr(4, nvolumes)
;    r = fltarr(4, nvolumes)
;
;    for k = 0,3 do begin
;        q[k,*] = 0
;        for j = 0,3 do begin
;            r[j,*] = 0
;            for i = 0,3 do begin
;                r[j,*] += u_x[i] * (*data)[fc[0]+i-1, fc[1]+j-1, fc[2]+k-1, *]
;             ;   r[j,*] = total(u_x * (*data)[fc[0]-1:fc[0]+2, fc[1]+j-1, fc[2]+k-1, *], 1)
;            endfor
;            q[k,*] += u_y[j] * r[j,*]
;        endfor
;        sum += u_z[k] * q[k,*]
;    endfor
;
;    return, reform(sum)

end

; Subroutine name: mas_tricubic_3d
; Created by: Bill Triplett
; Calling Information:
; data = mas_tricubic([x, y, z])

; Bugs or Important Comments to Developers:

; Purpose of subroutine:
; Returns pointwise interpolated data for sub-voxel
; points in a 4D data set (ie, diffusion)

; Editing Information:

function mas_tricubic_3d, coord, data_ptr=data_ptr

    data = data_ptr

    data_sz = size(*data, /dimensions)

    fc = float(floor(coord))

    if (fc[0] eq 0 or fc[1] eq 0 or fc[2] eq 0 OR $
        fc[0] ge data_sz[0]-2 or fc[1] ge data_sz[1]-2 or fc[2] ge data_sz[2]-2) then begin 
        return, reform((*data)[fc[0], fc[1], fc[2]])
    endif

    x = float(coord[0])
    y = float(coord[1])
    z = float(coord[2])
    
    dx = (x - fc[0]) & dx2 = dx^2 & dx3 = dx^3
    dy = (y - fc[1]) & dy2 = dy^2 & dy3 = dy^3
    dz = (z - fc[2]) & dz2 = dz^2 & dz3 = dz^3
    
    u_x = [ -0.5 * dx3 +       dx2 - 0.5 * dx, $ ;; [ -1  1 -1 ]
             1.5 * dx3 - 2.5 * dx2 + 1.0     , $ ;; [  1 -1  0 ]
            -1.5 * dx3 + 2.0 * dx2 + 0.5 * dx, $ ;; [ -1  1  1 ]
             0.5 * dx3 - 0.5 * dx2 ]             ;;  [ 1  -1  0 ]

    u_y = [ -0.5 * dy3 +       dy2 - 0.5 * dy, $
             1.5 * dy3 - 2.5 * dy2 + 1.0     , $
            -1.5 * dy3 + 2.0 * dy2 + 0.5 * dy, $
             0.5 * dy3 - 0.5 * dy2 ]

    u_z = [ -0.5 * dz3 +       dz2 - 0.5 * dz, $
             1.5 * dz3 - 2.5 * dz2 + 1.0     , $
            -1.5 * dz3 + 2.0 * dz2 + 0.5 * dz, $
             0.5 * dz3 - 0.5 * dz2 ]
    
    uu_x = fltarr(4,4,4)

    uu_x[*,*,0] = transpose([ [replicate(u_x[0], 4)], [replicate(u_x[1], 4)], [replicate(u_x[2], 4)], [replicate(u_x[3], 4)] ])
    uu_x[*,*,1] = uu_x[*,*,0]
    uu_x[*,*,2] = uu_x[*,*,0]
    uu_x[*,*,3] = uu_x[*,*,0]
    
;    for i = 0, 3 do begin
;        uu_x[*,*,i] = transpose([ [replicate(u_x[0], 4)], [replicate(u_x[1], 4)], [replicate(u_x[2], 4)], [replicate(u_x[3], 4)] ])
;    endfor
    
    uu_y  = transpose([ [replicate(u_y[0], 4)], [replicate(u_y[1], 4)], [replicate(u_y[2], 4)], [replicate(u_y[3], 4)] ])
    
    block = (*data)[fc[0]-1:fc[0]+2, fc[1]-1:fc[1]+2, fc[2]-1:fc[2]+2]
    sum = total( u_z * total( uu_y * total( uu_x * block , 1 ), 1 ) ) > 0.0

    return, sum
    
end


; Subroutine name: mas_apply_transformation_matrix
; Created by: Bill Triplett
; Calling Information:
;
;    matrix: 4D transformation matrix
;    
;    data_ptr: optional user-provided data pointer. If not specified, state1 is used
;
;    new_size: 3 array containing the [x,y,z] dimensions to crop the result.  
;
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine:
; applies a homogenous coordinate transformation matrix to data_ptr (or state1 of project.ci).
;
; Editing Information:

pro mas_apply_transformation_matrix, $
        matrix, $
        progress=progress, $
        src_ptr=src_ptr, $
        dest_ptr=dest_ptr, $
        from_file=from_file, $
        premultiply=premultiply, $
        postmultiply=postmultiply,$
        interp=interp

    forward_function mas_tricubic_3d

    sz_src  = size(*src_ptr, /dimensions)
    sz_dest = size(*dest_ptr, /dimensions)

    if (keyword_set(from_file)) then begin
        matrix = dblarr(4,4)
        matfile = dialog_pickfile(/file, /read)
        openr, lun, matfile, /get_lun
        readf, lun, matrix
        close, lun
        free_lun, lun
    endif

    if (not keyword_set(interp)) then interp = 'nn'

    case interp of 
        'nn': interpolation = 0
        'tricubic': interpolation = 1
        else: interpolation = 0
    end
    
    if (not keyword_set(premultiply)) then premultiply = diag_matrix(replicate(1D, 4))
    if (not keyword_set(postmultiply)) then postmultiply = diag_matrix(replicate(1D, 4))
    
    tx = invert(premultiply # matrix # postmultiply)

    if (keyword_set(progress)) then begin
        pbar = obj_new('progressbar', title='Applying Transformation Matrix', $
                        text='Applying Transformation Matrix', color='Red')
        pbar->start
    endif else begin
        progress = 0
    endelse

    n = 0
    for x = 0L, sz_dest[0]-1 do begin
        for y = 0L, sz_dest[1]-1 do begin
            for z = 0L, sz_dest[2]-1 do begin
                inn = ([x,y,z,1.] # tx)[0:2]
                if (inn[0] lt 0 or inn[0] gt sz_src[0]-1) then continue
                if (inn[1] lt 0 or inn[1] gt sz_src[1]-1) then continue
                if (inn[2] lt 0 or inn[2] gt sz_src[2]-1) then continue
                if (interpolation eq 0) then begin
                    (*dest_ptr)[x,y,z] = (*src_ptr)[round(inn[0]), round(inn[1]), round(inn[2])]
                endif else begin
                    (*dest_ptr)[x,y,z] =  mas_tricubic_3d(inn, data_ptr=src_ptr)
                endelse
            endfor
            if (progress) then begin
                n += sz_dest[2]
                pbar->update, n/product(sz_dest) * 100.0
                if (pbar->checkCancel()) then begin
                    pbar->destroy
                    return
                endif
            endif
        endfor
    endfor

    if (progress) then pbar->destroy
    
end
