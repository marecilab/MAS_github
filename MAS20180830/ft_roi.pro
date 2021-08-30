;; $Id$
;;

; ------------------------------------------------------------------------------
; ft_roi.pro
;
; The FT_ROI function returns a pointer that points to two pointers first of
; which is a pointer array whose elements point to the coordinates of the seed
; points corresponding to the user defined ROI, where the second points to
; an n by 6 array carrying the normal vector and center point information
; for the seed planes.
;
; (That is just plain hilarious.) - BT
;
; SYNTAX
;	Result=ft_roi(filename)
;
; ARGUMENTS
;  filename
;	A string of filename that the user will use to specify the ROI.
;	It has to be in the same dimensions with the files that will be used
;	in fiber-tracking
;
; Created on October 18, 2001 by E. Ozarslan
;
; ------------------------------------------------------------------------------
; called by ft_xroi
pro ft_calculate_each, ev

  compile_opt idl2, hidden

  common share_ftroi, interp, sz_im2d, slice, new_dimx, new_dimy, fdim, pdim, sdim;;, roi_objects
  slice_temp=slice

  widget_control, ev.top, get_uvalue=pState
  og=(*pState).oROIGroup
  oSelROI = (*pState).oSelROI
  if (OBJ_VALID(oSelROI) ne 0) then begin
    result = (*pState).oROIModel->IsContained(oSelROI, POSITION=i_roi)
	roii=og->get(pos=i_roi)

	if total(obj_valid(roii)) then begin
	  mask=roii->ComputeMask(dimensions =[new_dimx, new_dimy], MASK_RULE=2)
	  ;mask = 255 * (mask gt 0)
	  if interp eq 2 then $
	    mask = temporary(congrid(mask, sz_im2d[1], sz_im2d[2]))
      if interp eq 1 then $
    	mask = temporary(rebin(mask, sz_im2d[1], sz_im2d[2]))
      if interp eq 0 then $
    	mask = temporary(rebin(mask, sz_im2d[1], sz_im2d[2], /sample))

      if slice[0] eq 1. then begin
        locs=where(mask, lcount)
        if lcount ne 0 then begin
          seeds=make_array(lcount, 3)
          for i=0, lcount-1 do begin
            zi=floor(locs[i]/pdim)
            yi=locs[i]-zi*pdim
            seeds[i,*]=[round(slice[3]), yi, zi]
          endfor
        endif else ok=dialog_message('message1')
      endif else begin

      if slice[1] eq 1. then begin
        locs=where(mask, lcount)
        if lcount ne 0 then begin
          seeds=make_array(lcount, 3)
          for i=0, lcount-1 do begin
            zi=floor(locs[i]/fdim)
            xi=locs[i]-zi*fdim
            seeds[i,*]=[xi, round(slice[4]), zi]
          endfor
        endif else ok=dialog_message('message1')
      endif else begin

      if slice[2] eq 1. then begin
        locs=where(mask, lcount)
        if lcount ne 0 then begin
          seeds=make_array(lcount, 3)
          for i=0, lcount-1 do begin
            yi=floor(locs[i]/fdim)
            xi=locs[i]-yi*fdim
            seeds[i,*]=[xi, yi, round(slice[5])]
          endfor
        endif else ok=dialog_message('message1')
      endif

      endelse
      endelse
    endif else ok=dialog_message('This feature is not for oblique slices.')

    planar=ptr_new(seeds)
    ptr_seeds_arr=ptrarr(1, /ALLOCATE_HEAP)
    ptr_seeds_arr[0]=planar
    ptr_ptr_seeds_arr=ptrarr(2, /ALLOCATE_HEAP)
    *ptr_ptr_seeds_arr[0]=ptr_seeds_arr
    ptr_temp=ptr_new(slice_temp, /no_copy)
    ptr_ptr_seeds_arr[1]=ptr_temp
    ft_roi_calculate, '', [1,1,1,0], ptr_ptr_seeds_arr

  endif

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function ft_roi_extracter, ptr_im, file=file, slice, local_resx, local_resy, local_resz, interp

    compile_opt idl2, hidden
    ;;common share_fibers_main, pTHR, pFA, pS_0, pAvD, pDIR, loadedthr, loadeddir
    common share_ftroi, interp_, sz_im2d, slice_, new_dimx, new_dimy, fdim, pdim, sdim;;, roi_objects
    common share_fibers_roiobj, roi_objects
    common share_fibers_main

    slice_=slice

;;Used to crash here because 'interp' is not defined. Now it is
;;defined as the global variable.
    if (n_elements(interpol) eq 0) then begin
        if (n_elements(interp) eq 0) then begin
            interpol = 0
        endif else begin
            interpol = interp
        endelse
    endif
        
    interp_=interpol

    DEVICE, GET_SCREEN_SIZE=screenSize
    winX = (floor(screenSize[0] / 1.25))  > 320 ; View window size.
    winY = (floor(screenSize[1] / 1.25))  > 320 ; View window size.
    winX = 8 * ((winX < winY) / 8)
    winY = winX
    
    ptr_temp=ptr_im
    im3d=*ptr_temp

    sz=size(im3d)
    fdim=sz[1]
    pdim=sz[2]
    sdim=sz[3]
    
    if slice[0] eq 1. then begin
        im2d=reform(im3d[round(slice[3]),*,*])
        tilesize=max([local_resy*pdim/winx, local_resz*sdim/winx])
        zoom1=local_resy/tilesize
        zoom2=local_resz/tilesize
    endif else begin

        if slice[1] eq 1. then begin
            im2d=reform(im3d[*,round(slice[4]),*])
            tilesize=max([local_resx*fdim/winx, local_resz*sdim/winx])
            zoom1=local_resx/tilesize
            zoom2=local_resz/tilesize
        endif else begin
            
            if slice[2] eq 1. then begin 
                ;; Mistake in this code if you chose the last slide (everything might be offset by 1 slice)
                im2d=reform(im3d[*,*,round(slice[5])])
                tilesize=max([local_resx*fdim/winx, local_resy*pdim/winx])
                zoom1=local_resx/tilesize
                zoom2=local_resy/tilesize
            endif else begin
                
                nx=slice[0]
                ny=slice[1]
                nz=slice[2]

                quad=sqrt(nx^2+ny^2)

                ip=[ny,-nx,0]/quad
                jp=[nx*nz,ny*nz,-quad^2]/quad
                kp=[nx,ny,nz]
                kp=kp/norm(kp)
                jp=crossp(kp,ip)

                A=[ [ ip[0], ip[1], ip[2] ],$
                    [ jp[0], jp[1], jp[2] ],$
                    [ kp[0], kp[1], kp[2] ] ] 

                Cen=[slice[3:5]]
                zref=round(nx*slice[3]+ny*slice[4]+nz*slice[5])
                maxdim=max([fdim,pdim,sdim])
                info=make_array(3*maxdim^2, 6)
                l=0

                for i=0,fdim-1 do begin
                    for j=0,pdim-1 do begin
                        for k=0,sdim-1 do begin

                            if round(nx*i+ny*j+nz*k) eq zref then begin
                                info[l,*]=[i,j,k,im3d[i,j,k], $
                                           round(ny*i/quad-nx*j/quad), $
                                           round(nx*nz*i/quad+ny*nz*j/quad-quad*k)]
                                l=l+1
                            endif

                        endfor
                    endfor
                endfor
                
                info1=make_array(l,6)

                for i=0,l-1 do info1[i,*]=info[i,*]

                info=0
                shiftx=min(info1[*,4])
                shifty=min(info1[*,5])
                info1[*,4]=info1[*,4]-shiftx
                info1[*,5]=info1[*,5]-shifty
                dims2dx=max(info1[*,4])+1
                dims2dy=max(info1[*,5])+1
                im2d=make_array(dims2dx,dims2dy)

                for i=0,l-1 do im2d[info1[i,4],info1[i,5]]=info1[i,3]
                for j=1,dims2dy-2 do for i=1,dims2dx-2 do begin
                    if im2d[i,j] eq 0. then begin
                        cc=0
                        aa=0.
                        for ii=-1,1 do begin
                            for jj=-1,1 do begin
                                if im2d[i+ii,j+jj] ne 0. then begin
                                    aa=aa+im2d[i+ii,j+jj]
                                    cc=cc+1
                                endif
                            endfor
                        endfor

                        im2d[i,j]=aa/cc

                    endif

                endfor

                im2d=temporary(reverse(reverse(im2d, 2)))
                sz_im2d=size(im2d)
                zoom1=round(winx/sz_im2d[1])
                zoom2=zoom1

            endelse
        endelse
    endelse
    
    sz_im2d=size(im2d)
    zoom2 = zoom1

    if (gfile eq 'Orientation') then begin
        
        if (slice[0] eq 1) then begin
            ;; noop
        endif else if (slice[1] eq 1) then begin
            ;; noop
        endif else if (slice[2] eq 1) then begin
            ;;im2d=reform(im3d[*,*,round(slice[5])])
            faimage = im2d
            evecs = *pDIR
            sz_adt = size(faimage)
            slice_axis = 0
            orientimage = fltarr(sz_adt[1],sz_adt[2],3)
            tempevecs =	fltarr(sz_adt[1],sz_adt[2], 3)
            
            tempfa = faimage
            tempevecs = REFORM(evecs[*,*,round(slice[5]),*])
            orientimage = orient(tempfa, tempevecs , sz_adt[1], sz_adt[2], sz_adt[1],sz_adt[2], slice_axis)
            im2d = orientimage
        endif
        
    endif
    
; This is where the interpolation type is defined.
; 	0 for linear
;  	1 for
;	2 for Congrid
    
    if (gfile ne 'Orientation') then begin ;Any other dataset except Orientation
        
        if interp eq 0 then begin
            new_dimx=floor(zoom1)*sz_im2d[1]
            new_dimy=floor(zoom2)*sz_im2d[2]
            im2d=temporary(rebin(im2d, new_dimx, new_dimy, /SAMPLE))
        endif
        if interp eq 1 then begin
            new_dimx=floor(zoom1)*sz_im2d[1]
            new_dimy=floor(zoom2)*sz_im2d[2]
            im2d=temporary(rebin(im2d, new_dimx, new_dimy))
        endif
        if interp eq 2 then begin
            new_dimx=round(sz_im2d[1]*zoom1)
            new_dimy=round(sz_im2d[2]*zoom2)
            new_dimx = new_dimx[0]
            new_dimy = new_dimy[0]
            
            im2d=temporary(congrid(im2d, new_dimx, new_dimy, cubic=-.5))
        endif
        
    endif else begin       ; the orientation file (added for 20070830)
        sz_im2d = size(im2d)
        
        if (interp eq 0) then begin
            
            new_dimx=floor(zoom1)*sz_im2d[1]
            new_dimy=floor(zoom2)*sz_im2d[2]
            temp_image2d = fltarr(new_dimx, new_dimy, 3)
            
            temp_image2d = rebin(im2d, new_dimx, new_dimy, 3)
            bytimage = temp_image2d
            
        endif else if (interp eq 1) then begin
            
            new_dimx=floor(zoom1)*sz_im2d[1]
            new_dimy=floor(zoom2)*sz_im2d[2]
                                ;temp_image2d = fltarr(new_dimx, new_dimy, 3)
            
            temp_image2d = rebin(im2d, new_dimx, new_dimy, 3)
            bytimage = temp_image2d
            
        endif else if (interp eq 2) then begin
            
            new_dimx=round(sz_im2d[1]*zoom1)
            new_dimy=round(sz_im2d[2]*zoom2)
            new_dimx = new_dimx[0]
            new_dimy = new_dimy[0]
            temp_image2d = fltarr(new_dimx, new_dimy, 3)
            
            
            temp_image2d=congrid(im2d, new_dimx, new_dimy, 3, cubic=-.5)
            bytimage = temp_image2d
            
        endif
        
        
    endelse                     ; For the interpolation IF statement
    
                                ; HS 20061115
                                ; Now I will check to see if file contains anything.
    if (N_ELEMENTS(file) EQ 0)  THEN BEGIN
        if (gfile ne 'Orientation') then begin
                                ; Doing the same thing as Evren did for the 'else' case on the 'Aver_D' search
            bytimage=bytscl(im2d, max=max(im3d))
        endif
        
    endif else begin
        
        if strmid(file, 5,6,/reverse_offset) eq 'Aver_D' then $
          bytimage=bytscl(im2d, max=.0008) else $
          bytimage=bytscl(im2d, max=max(im3d))
        
    endelse
    
                                ;im2d=temporary(congrid(im2d, round(sz_im2d[1]*zoom1), $
;					round(sz_im2d[2]*zoom2), cubic=-.5))
                                ;bytimage=bytscl(im2d, max=1.0)
    if slice[0] eq 1 or slice[1] eq 1 or slice[2] eq 1 then begin
        if interp eq 2 then newydim=round(zoom2/zoom1*sz_im2d[2]) $
        else newydim=(round(zoom2/zoom1)>1)*sz_im2d[2]
        
;;                                ; This is what launches the window with the 4 images.
;;        ft_display_adt, slice[1]+2*slice[2], $
;;          round(slice[3+slice[1]+2*slice[2]]) $
;;          ,new_dims=[sz_im2d[1], newydim]
        
    endif
    
; Launches the ROI selection window
    ft_xroi, bytimage, title='ROI for Fiber Tracking', $
      REGIONS_OUT=regions, renderer=0, /block ; /modal, group=0L

    if total(obj_valid(regions)) then begin

        sz_regions=size(regions)
        
        if sz_regions[0] then n_rois=sz_regions[1] else n_rois=1
        
;;        mask = regions[0]->ComputeMask(dimensions=[new_dimx, new_dimy], mask_rule=2)
        
        roimask = bytarr(n_rois, sz_im2d[1], sz_im2d[2])
;;;BOOKMARK
        for i=0, n_rois-1 do begin
            roi = regions[i]->ComputeMask(dimensions=[new_dimx, new_dimy], mask_rule=2)
            roimask[i,*,*] = congrid(roi, sz_im2d[1], sz_im2d[2])

            locs = where(roimask[i,*,*], lcount)
            
            if (lcount eq 0) then begin
               ;; return ptr_new()
            endif else begin

                roiseeds = make_array(lcount, 3)

                if (slice[0] eq 1.0) then begin
            
                    for pt=0, lcount-1 do begin
                        zi = floor(locs[pt]/pdim)
                        yi = locs[pt]-zi*pdim
                        roiseeds[pt,*] = [round(slice[3]), yi, zi]
                    endfor

                endif else if (slice[1] eq 1.0) then begin

                    for pt=0, lcount-1 do begin
                        zi = floor(locs[pt]/fdim)
                        xi = locs[pt]-zi*fdim
                        roiseeds[pt,*] = [xi, round(slice[4]), zi]
                    endfor

                endif else if (slice[2] eq 1.0) then begin

                    for pt=0, lcount-1 do begin
                        yi = floor(locs[pt]/fdim)
                        xi = locs[pt]-yi*fdim
                        roiseeds[pt,*] = [xi, yi, round(slice[5])]
                    endfor

                endif else begin
                    ;; oblique inverse not supported
                    continue
                endelse
                
                roi_objects = [roi_objects, ptr_new(transpose(roiseeds), /no_copy)]
                
            endelse
            
        endfor

        mask = make_array(new_dimx, new_dimy)

        for i=0, n_rois-1 do begin
            roi = regions[i]->ComputeMask(dimensions=[new_dimx, new_dimy], mask_rule=2)
            mask = mask + roi
            mask = 255 * (mask gt 0)
        endfor
        roimask = 255 * (roimask gt 0)

        ;;mask = temporary(congrid(mask, sz_im2d[1], sz_im2d[2], interp))
        if interp eq 2 then $
          mask = temporary(congrid(mask, sz_im2d[1], sz_im2d[2]))
        if interp eq 1 then $
          mask = temporary(rebin(mask, sz_im2d[1], sz_im2d[2]))
        if interp eq 0 then $
          mask = temporary(rebin(mask, sz_im2d[1], sz_im2d[2], /sample))
        
        if slice[0] eq 1. then begin
;;            mask = temporary(mask)
            locs=where(mask, lcount)
            if lcount ne 0 then begin
                seeds=make_array(lcount, 3)
                for i=0, lcount-1 do begin
                    zi=floor(locs[i]/pdim)
                    yi=locs[i]-zi*pdim
                    seeds[i,*]=[round(slice[3]), yi, zi]
                endfor
            endif else return, ptr_new()

        endif else begin
            
            if slice[1] eq 1. then begin
                locs=where(mask, lcount)
                if lcount ne 0 then begin
                    seeds=make_array(lcount, 3)
                    for i=0, lcount-1 do begin
                        zi=floor(locs[i]/fdim)
                        xi=locs[i]-zi*fdim
                        seeds[i,*]=[xi, round(slice[4]), zi]
                    endfor
                endif else return, ptr_new()

            endif else begin
                
                if slice[2] eq 1. then begin
                    locs=where(mask, lcount)
                    if lcount ne 0 then begin
                        seeds=make_array(lcount, 3)
                        for i=0, lcount-1 do begin
                            yi=floor(locs[i]/fdim)
                            xi=locs[i]-yi*fdim
                            seeds[i,*]=[xi, yi, round(slice[5])]
                        endfor
                    endif else return, ptr_new()

                endif else begin
                    
                    ;; oblique inverse
                    mask=temporary(reverse(reverse(mask, 2)))
                    locs=where(mask, lcount)
                    if lcount ne 0 then begin
                        sds=make_array(lcount*10,3)
                        ic=0
                        for i=0, lcount-1 do begin
                            vi=floor(locs[i]/dims2dx)
                            hi=locs[i]-vi*dims2dx
                            for j=0,l-1 do begin
                                if info1[j,4] eq hi and info1[j,5] eq vi then begin
                                    sds[ic,*]=[info1[j,0],info1[j,1],info1[j,2]]
                                    ic=ic+1
                                endif
                            endfor
                            seeds=make_array(ic,3)
                        endfor
                        for i=0,ic-1 do seeds[i,*]=sds[i,*]
                    endif else return, ptr_new()
                    
                endelse
            endelse
        endelse
        
        planar=ptr_new(seeds, /NO_COPY)
        
    endif else planar=ptr_new()
    
    return, planar
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;function ft_roi_seeder, ptr_im, file, resx, resy, resz, interp
function ft_roi_seeder, ptr_im, resx, resy, resz, interp, file=file

    compile_opt idl2, hidden

    ;; This is the comparison that determines how to call func_rd_slicer
    ;; When reading from MAS we don't need to call the function with the 'file' option.
    if N_ELEMENTS(file) eq 0 then begin
        print, 'Inside of ft roi seeder'
	    slices = func_rd_slicer(/ROI, aspect_ratio=resz/resx, POINTER=ptr_im)

    endif else begin

        slices=func_rd_slicer([file+'.flt'], /ROI, aspect_ratio=resz/resx, POINTER=ptr_im)
        ;; this is the endelse for the optional parameter "file" comparison.
    endelse

    sz_slices=size(slices)
    
    if sz_slices[0] ne 0 then begin
        if sz_slices[0] eq 2 then nnn=sz_slices[2] else nnn=1
        
        ptr_seeds_arr=ptrarr(nnn, /ALLOCATE_HEAP)
        
        for i=0, nnn-1 do begin
            ptr_seeds_arr[i]=ft_roi_extracter(ptr_im,file=file,slices[*,i],resx,resy,resz, interp)
        endfor
        
        ;;wdelete, 8
        ptr_temp = ptr_new(slices, /no_copy)
        ptr_ptr_seeds_arr = ptrarr(2, /ALLOCATE_HEAP)
        
        *ptr_ptr_seeds_arr[0] = ptr_seeds_arr
        ptr_ptr_seeds_arr[1] = ptr_temp
        
    endif else begin
        
        ptr_ptr_seeds_arr=ptr_new()
        
    endelse
    
    return, ptr_ptr_seeds_arr

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;function ft_roi, file, resx, resy, resz, interp
function ft_roi, local_resx, local_resy, local_resz, local_interp, file=file

  compile_opt idl2
  common scan_data
  ;common share_fibers_main, pTHR, pFA, pS_0, pAvD, pDIR, loadedthr, loadeddir
  common share_fibers_main

; need to call the appropriate function call when im reading data from MAS and not from a file.
IF N_ELEMENTS(loadedthr) EQ 0 THEN BEGIN

	; ptr_im= project.dataArray[project.CI].frac_Ani


	; HS - 20061114
	; From now on the file input is optional so that you can skip this if i'm loading data from MAS.
	; ptr_seeds=ft_roi_seeder(pTHR, file, resx, resy, resz, interp)

; Now we are going to decide what image type to use when selecting the ROI (parameter gfile).

	if gfile eq 'Fract_Anisotropy' then begin

		; Here I'm using the FA to guide the ROI
		ptr_seeds=ft_roi_seeder(pTHR, local_resx, local_resy, local_resz, local_interp)

	endif else if gfile eq 'Average Diffusion' then begin

		; Here I'm using the Average Diffusion to guide the ROI
		ptr_seeds=ft_roi_seeder(pAvD, local_resx, local_resy, local_resz, local_interp)

	endif else if gfile eq 'S0' then begin

		; Here I'm using the S0 to guide the ROI
		ptr_seeds=ft_roi_seeder(pS_0, local_resx, local_resy, local_resz, local_interp)

	endif else if gfile eq 'Orientation' then begin
;		print, 'Testing Orientation'
;
;		faimage = *pTHR
;		evecs = *pDIR
;		sz_adt = size(faimage)
;		slice_axis = 0
;		orientimage = fltarr(sz_adt[1],sz_adt[2], sz_adt[3],3)
;		tempevecs =	fltarr(sz_adt[1],sz_adt[2], 3)
;
;		; Have to do the orient image for every single slice in the data
;		FOR i=0,sz_adt[3]-1 DO BEGIN
;			tempfa = faimage[*,*,i]
;			tempevecs = REFORM(evecs[*,*,i,*])
;   			orientimage[*,*,i,*] = orient(tempfa, tempevecs , sz_adt[1], sz_adt[2], sz_adt[1],sz_adt[2], slice_axis)
;			;colorimage = color_mapper(bytscl(tempimage) , sz_adt[1], sz_adt[2], 0)
;		ENDFOR

		;ptr_seeds=ft_roi_seeder(ptr_new(orientimage), local_resx, local_resy, local_resz, local_interp)




		ptr_seeds=ft_roi_seeder(pTHR, local_resx, local_resy, local_resz, local_interp)
	endif

	return, ptr_seeds

endif else begin

  case file of
    'Fract_Anisotropy' 	: begin
       if loadedthr eq 0 then ptr_im = pTHR else begin
         im_temp=rflt('Fract_anisotropy.flt')
         ptr_im=ptr_new(im_temp, /no_copy)
       endelse
     end
    'S0'				: ptr_im = pS_0
    'Aver_D'			: ptr_im = pAvD
    'dir_coherence'		: begin
      if loadedthr eq 1 then ptr_im = pTHR else begin
         im_temp=rflt('dir_coherence.flt')
         ptr_im=ptr_new(im_temp, /no_copy)
       endelse
    end
    'fran_coh'		: begin
      if loadedthr eq 2 then ptr_im = pTHR else begin
         im_temp=rflt('fran_coh.flt')
         ptr_im=ptr_new(im_temp, /no_copy)
       endelse
    end
    else : begin
      if file_test(file+'.flt') eq 0 then return, 0
      im_temp=rflt(file+'.flt')
      ptr_im=ptr_new(im_temp, /no_copy)
    end
  endcase

	; Since 'file' is now optional we need to modify the following line
	;ptr_seeds=ft_roi_seeder(pTHR, file, resx, resy, resz, interp)
  ptr_seeds=ft_roi_seeder(ptr_im, resx, resy, resz, interp, file=file)
  return, ptr_seeds

 endelse
end
