;; $Id$
;;
function mc_compute_mi_local, im_a, im_b, bs, normalized=normalized

 ;;   A = im_a/max(im_a)
 ;;   B = im_b/max(im_b)
    A = im_a
    B = im_b

    A = bytscl(A);;hist_equal(A))
    B = bytscl(B);;hist_equal(B))

    sz_a = size(A)
    sz_b = size(B)

    pab = float(hist_2d(A, B, bin1=bs, bin2=bs))

    hist_a = float(histogram(A, binsize=bs))
    hist_a = hist_a/total(hist_a)

    hist_b = float(histogram(B, binsize=bs))
    hist_b = hist_b/total(hist_b)        
    
    Ha = 0.0
    Hb = 1.0
    if (0 and keyword_set(normalized)) then begin
        n_save = where(hist_a ne 0)
        Ha = total((-hist_a(n_save)*alog(hist_a[n_save])))/alog(2)
        n_save = where(hist_b ne 0)
        Hb = total((-hist_b(n_save)*alog(hist_b[n_save])))/alog(2)
    endif

    pab = pab / total(pab)
    pa_pb = hist_a # transpose(hist_b)
    save = where(pa_pb gt 1e-12)
    u = pab(save)
    v = pa_pb(save)
    t = where(u gt 1e-12)
    y = u[t] * (alog(u[t]/v[t])/alog(2))

    return, total(y)/(Ha+Hb)

end

function array2string, array
    return, strcompress('['+string(array[0])+','+string(array[1])+']',/remove_all)
end

function mc_pad_image, image, pad_x=pad_x, pad_y=pad_y, out_x=out_x, out_y=out_y

    sz = size(image)
    sx = sz[1]
    sy = sz[2]

    amt_x = 0
    amt_y = 0

    if not keyword_set(pad_x) then pad_x = 0
    if not keyword_set(pad_y) then pad_y = 0

    if (sx gt sy) then begin
        amt_y = sx - sy
    endif else if (sy gt sx) then begin
        amt_x = sy - sx
    endif

    out_x = amt_x + 2*pad_x
    out_y = amt_y + 2*pad_y

    new_sx = sx + out_x
    new_sy = sy + out_y
    
    out = fltarr(new_sx, new_sy)

    out[pad_x:sx+pad_x-1,pad_y:sy+pad_y-1] = image

    out = shift(out, amt_x/2, amt_y/2)

    return, out

end

function mc_unpad_image, image, amt_x, amt_y

    sz_x = (size(image))[1]
    sz_y = (size(image))[2]
    
    new_sz_x = sz_x-amt_x
    new_sz_y = sz_y-amt_y
    
    out = fltarr(new_sz_x, new_sz_y)

    out = (shift(image, -amt_x/2, -amt_y/2))[0:new_sz_x-1, 0:new_sz_y-1]

    return, out

end

function mc_do_shear, im, amt, l_freq=l_freq, f_order=f_order

    compile_opt idl2

    sz_im = size(im)
    sx = sz_im[1]
    sy = sz_im[2]

    pad_x = 0
    pad_y = 0

;    if amt eq 0 then return, im
    if abs(amt) lt 1e-8 then return, im

    workspace = im 

;    if keyword_set(l_freq) then begin
;        if not keyword_set(f_order) then f_order = 1
;        bfilt = butterworth(sx, sy, cutoff=l_freq, /origin, order=f_order)
;    endif else begin
;        bfilt = fltarr(sx, sy)
;        bfilt[*,*] = 1.
;    endelse
        
    f_amt = amt*!PI/180/sx
    
    imft = shift(fft(workspace, 1), fix(sx/2),fix(sy/2))

;    kern = complexarr(sx,sy)
;    kern[*,*] = 1

    wub = findgen(sx) - fix(float(sx/2))

    C = 2*!PI*f_amt

    kern = exp( (C*wub # transpose(complex(0.0, findgen(sx)))) )

;; The above line is equivalent to the following:
;;    for k = 0, sy-1 do begin
;;        for l = 0,sx-1 do begin
;;            
;;            kern[l,k] = exp(C*float(wub[l])*complex(0.0,k))
;;
;;        endfor
;;    endfor

    res = complexarr(sx,sy)

    res = fft(imft, -1, dimension=2)
    res = res * transpose(kern)
    ;res = res * bfilt
    res = fft(res, -1, dimension=1)

    sh = abs(res)

    return, sh[0:sz_im[1]-1,0:sz_im[2]-1]

end

function mc_do_shear_x, im, amt, l_freq=l_freq, f_order=f_order
    compile_opt idl2
    return, mc_do_shear(im, amt, l_freq=l_freq, f_order=f_order)

end

function mc_do_shear_y, im, amt, l_freq=l_freq, f_order=f_order
    compile_opt idl2    
    return, transpose(mc_do_shear(transpose(im), amt, l_freq=l_freq,$
                            f_order=f_order))

end

function mc_do_scale, image, xscale, yscale

    if abs(xscale) lt 1e-8 and abs(yscale) lt 1e-8 then return, image
    factor = 2

    im_size = size(image)
    sx = im_size[1]
    sy = im_size[2]
    
    workspace = dblarr(sx*factor, sy*factor)
    
    mid_x = (sx*factor)/2
    mid_y = (sy*factor)/2

    ;;workspace[mid_x-sx:mid_x+sx-1, mid_y-sy:mid_y+sy-1] = image
    workspace[sx-sx/factor:sx+sx/factor-1, sy-sy/factor:sy+sy/factor-1] = image

    scaled = congrid(workspace, factor*sx+xscale,factor*sy+yscale, /center, cubic=-0.5)
    scaled = shift(scaled, -xscale/factor, -yscale/factor)
    
;    out = bytarr(sx,sy)
;    out = round(scaled[sx-sx/factor:sx+sx/factor-1, sy-sy/factor:sy+sy/factor-1])

    out = dblarr(sx,sy)
    out = scaled[sx-sx/factor:sx+sx/factor-1, sy-sy/factor:sy+sy/factor-1]

    return, out

end


;; Subroutine name: mc_do_fft_shift
;;
;; Created by: BT 2008-03-10
;;
;; Calling Information:
;;  image  - the image data to be shifted
;;  xshift - pixel shift amount in x direction
;;  yshift - pixel shift amount in y direction
;;  NO_FFT - set this keyword to indicate that image is 
;;           already in frequency domain
;;
;; Bugs or Important Comments to Developers:
;;
;; Purpose of subroutine:
;;  Shifts a single image data using fourier shift thm.
;;
;; Editing Information:
;;
function mc_do_fft_shift, image, xshift, yshift, NO_FFT=nofft

    compile_opt idl2
    im_size = size(image)
    sx = im_size[1]
    sy = im_size[2]

    ;; assume data is in freq domain unless told otherwise
    if keyword_set(nofft) then begin
        ;; yuck
    endif else begin
        nofft = 0
    endelse

    if nofft then begin
        ;; image data is already in frequency domain
        ;; NO_FFT = No FFT Required
        tmp_image = ptr_new(image)
        fft_image = *tmp_image
    endif else begin
        ;; image data is in magnitude mode
        tmp_image = ptr_new(double(image))
        fft_image = shift(fft(*tmp_image), fix(sx/2), fix(sy/2))
;        fft_image = shift(fft_image, fix(sx/2), fix(sy/2))
    endelse

    x_val = 2.0*!PI*xshift/sx
    y_val = 2.0*!PI*yshift/sy

    svec = complexarr(sx, sy)
    ;;svec = exp(-complex(0.0,findgen(sx)*y_val)) # exp(-transpose(complex(0.0, findgen(sy)*x_val)))
    svec = matrix_multiply( exp(-complex(0.0, findgen(sx)*y_val)), $
                            exp(-complex(0.0, findgen(sy)*x_val)), /btranspose)
;;;;
;; The above line is equivalent to the following:    
;;    for i = 0, sx-1 do begin 
;;        svec[i,*] = replicate(exp(-complex(0.0,i*y_val)), sy)
;;    endfor
;;    
;;    for i = 0, sy-1 do begin 
;;        svec[*,i] = svec[*,i] * replicate(exp(-complex(0.0,i*x_val)), sx)
;;    endfor
;;;;    
    new_fft_image = fft_image * svec

    if nofft then begin
        out_image = new_fft_image
    endif else begin
        out_image = shift(new_fft_image, fix(sx/2), fix(sy/2))
        out_image = abs(fft(out_image, /INVERSE))
    endelse

    ptr_free, tmp_image

    return, out_image

end

function mc_do_rotate, im, theta, f=f, o=o

    compile_opt idl2
    if (not keyword_set(f)) then f=5000
    if (not keyword_set(o)) then o=1

   if abs(theta) lt 1e-8 then return, im

    sz_im = size(im)
    sz_im_org = sz_im

    w_new = sz_im[1] * sin(theta*!PI/180) + sz_im[2] * cos(theta*!PI/180)
    h_new = sz_im[2] * sin(theta*!PI/180) + sz_im[1] * cos(theta*!PI/180)
    
    pd_x = abs(round((h_new-sz_im[1])/2.0))
    pd_y = abs(round((w_new-sz_im[2])/2.0))
    
    sz_im_org[1] += 2*pd_x
    sz_im_org[2] += 2*pd_y

    outx=0
    outy=0
    
    im2 = im

    sz_im = size(im2)
    sx = sz_im[1]
    sy = sz_im[2]

    x_dir = theta
    y_dir = theta

    sh = im2
    sh1 = mc_do_shear_x(sh,  float(x_dir/2))
    sh2 = mc_do_shear_y(sh1, float(-y_dir))
    sh3 = mc_do_shear_x(sh2, float(x_dir/2))

    out = sh3

    return, out

end

function mc_transform_image, image, arg_vector, scale_first=scale_first

    common _mc_working_set, working_set, _mc_working_set

    mn = min(image)
    mx = max(image)

    if (keyword_set(scale_first)) then begin
        if (scale_first ne 0) then begin
            scale_first = 1
        endif
    endif else begin
        scale_first = 0
    endelse

    temp = image

    if (n_elements(arg_vector) eq 7 and scale_first) then begin
        if (arg_vector[5] ne 0 or arg_vector[6] ne 0) then begin
            temp = mc_do_scale(temp, $
                               round(arg_vector[5]),$
                               round(arg_vector[6]))
        endif
    endif

    if (arg_vector[2] ne 0) then begin
        temp = mc_do_shear_x(temp,$
                             arg_vector[2])
    endif

    if (arg_vector[3] ne 0) then begin
        temp = mc_do_shear_y(temp,$
                             arg_vector[3])
    endif

    if (arg_vector[4] ne 0) then begin
        temp = mc_do_rotate(temp, arg_vector[4])
    endif

    if (arg_vector[0] ne 0 or arg_vector[1] ne 0) then begin
        temp = mc_do_fft_shift(temp, arg_vector[0], arg_vector[1])
    endif
    
    if (n_elements(arg_vector) eq 7 and not scale_first) then begin
        if (arg_vector[5] ne 0 or arg_vector[6] ne 0) then begin
            temp = mc_do_scale(temp, $
                               round(arg_vector[5]),$
                               round(arg_vector[6]))
        endif
    endif
    
    return, temp

end
