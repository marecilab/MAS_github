;; $Id$
;; Copyright 2003 University of Florida. All Rights Reserved


; Subroutine name: mas_make_sheet_display
; Created by:
; Calling Information:

; Bugs or Important Comments to Developers:


; Purpose of subroutine:
; the idea behind this is that the user would select multi and
; they would get a display that has multiple slices ready to display
; Input
;   data is pointer to the images,
;   x_fov is the x direction fov for 1 of the images
;   y_fov is the y direction fov for 1 of the images
; Returns
;   sheet.xfov
;   sheet.yfov
;   sheet.image

; Editing Information:
    ;Edited by HS 2006/10/04
    ;Fix spelling mistakes and commenting.

function mas_make_sheet_display, data, x_fov, y_fov
    common scan_data
    CI = project.ci

    ;get the size of the data to display
    sz_data = size((*data))
    x_dim = sz_data[1]
    y_dim = sz_data[2]
    frames = sz_data[3]



    x_frames = ceil(sqrt(frames))
    y_frames = ceil(float(frames)/float(x_frames))


    ;depending upon what type the data at the end of the
    ;pointer is, is what we have to initalize image as.
    if sz_data[4] eq 4 then begin
        ;data must be floating point
        image = fltarr(x_frames*x_dim,y_frames*y_dim)

    endif else begin
        ;data must be byte
        image = bytarr(x_frames*x_dim,y_frames*y_dim)
    end

    ;initalizing the structure
    sheet= {x_fov:x_frames*x_fov $
           ,y_fov:y_frames*y_fov $
           ,image:ptr_new()}

    ;copying the data into image
    for jj=0, y_frames-1 do begin
        for ii=0, x_frames-1 do begin
            if (ii+x_frames*jj ) le frames-1 then begin
                ;nr- image[ii*x_dim,jj*y_dim] = (*data)[*,*,ii+x_frames*jj]
                ;nr+ Flip image horizontally
                image[ii*x_dim,jj*y_dim] = reverse( (*data)[*,*,ii+x_frames*jj] , 2 )
            end
        end
    end

    ;reversing and assigning the image at the end of the pointer
    ;structure.
    sheet.image = ptr_new(reverse(image,2))

    return,sheet

end
