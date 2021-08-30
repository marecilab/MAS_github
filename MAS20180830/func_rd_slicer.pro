;; $Id$
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; func_rd_slicer.pro
;
; The FUNC_RD_slicer function reads any number of *.flt files into
; SLICER3 for combined display. Each file must have the same size.
; It returns a 0 value if SLICER3 is run, and the 'planes' array
; specifying the planes chosen if ROI keyword is set.
;
; It now has support for slicing data passed from MAS.
;
;
;
; SYNTAX
; 	Func_rd_slicer, files [, /ROI][, POINTER=pointer]
;
; ARGUMENTS
;  files
;	A string array of file names to be displayed.  The
;	filenames should include the pathname (if they are not
;	in the current directory) and .flt extension.
;
; KEYWORDS
;  /roi
;	When this keyword is set, it reads a file into ROI3D_SLICER.
;
;  pointer
;  This keyword is set to a pointer pointing to the 3d image array.
;  In this case rflt routine is not called.
;
; Created on July 12,  2001 by E. Ozarslan
; by modifying the already existing code pro_rd_slicer.pro
;
; Last modified on September 25, 2003 by E. Ozarslan
; Modified by HS on 20061115.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function func_rd_slicer, files=files, ROI=roi, aspect_ratio=aspect_ratio, POINTER=pointer

compile_opt idl2
common share_fibers_main
; Now we check to see if the files parameter was included. If it wasn't
; then we are loading data from MAS.
if N_ELEMENTS(files) eq 0 then begin

	; This name will be passed to the roi3D_slicer function and it will
	; be displayed in that window.
        if (n_elements(gfile) eq 0) then begin
            gfile = 'UNKNOWN'
        endif
	names=gfile

		;print, 'Inside func_rd_slicer'

		; Reading the data and getting the different dimensions.
		im=*pointer
		sz=size(im)
		fdim=sz[1]
		pdim=sz[2]
		sdim=sz[3]

		; Creating a byte array with the dimensions of the data passed.
		images=bytarr(fdim,pdim,sdim)
		; Byte scaling the original data.
		images[*,*,*] = bytscl(im)

		; If user didn't specify an aspect ratio then we assume 1.
		if n_elements(aspect_ratio) eq 0 then aspect_ratio=1.

		; Creating a pointer array with 1 dimension
		ptarr=PTRARR(1, /ALLOCATE_HEAP)

		; Passing the data to slice to this pointer array.
		*ptarr[0] = images[*,*,*]

		; Calling the actual slicer function. Here we will select the planes
		planes=roi3d_slicer(ptarr, data_names=names, /modal, aspect_ratio=aspect_ratio)

		; Testing the returned data
		print, 'Inside func_rd_slicer'
		print, planes

		; Returning the planes selected with the previous tool.
		return, planes


; This is when we are loading data from files.
endif else begin

	num_files=(size(files))[1]

	if keyword_set(roi) then if (num_files ne 1) then begin
		ok=dialog_message('One file should be specified for ROI slicing!'$
			, /ERROR)
		return, 0
	endif

	if keyword_set(pointer) then im=*pointer else im=rflt(files[0])

	sz=size(im)
	fdim=sz[1]
	pdim=sz[2]
	sdim=sz[3]
	filecount=0
	names=['a']

	images=bytarr(fdim,pdim,sdim,1)


	for fn=0, num_files-1 do begin

	  if fn ne 0 then im=rflt(files[fn])

	  sz_im=size(im)

	  a2=Strpos(files[fn], '/', /Reverse_Search)

	  file_name=Strmid(files[fn], a2+1)

	  a3=Strpos(file_name, '.')

	  file_name=Strmid(file_name, 0, a3)

	  if sz_im[0] eq 3 or sz_im[0] eq 4 then $
	  	if sz_im[1] eq fdim and sz_im[2] eq pdim and sz_im[3] eq sdim then begin

	    	if sz_im[0] eq 3 then begin

	      		if fn eq 0 then images=bytscl(im) else begin
			        imagesnew=bytarr(fdim,pdim,sdim,filecount+1)
			        imagesnew[*,*,*,0:filecount-1]=images
			        imagesnew[*,*,*,filecount]=bytscl(im)
			        images=imagesnew
	      		endelse

		      filecount=filecount+1
		      names=[names, file_name]

	    	endif

		    if sz_im[0] eq 4 then begin

		      for ii=0,sz_im[4]-1 do begin

		        if fn+ii eq 0 then images=bytscl(reform(im[*,*,*,ii])) else begin
		          imagesnew=bytarr(fdim,pdim,sdim,filecount+1)
		          imagesnew[*,*,*,0:filecount-1]=images
		          imagesnew[*,*,*,filecount]=bytscl(reform(im[*,*,*,ii]))
		          images=imagesnew
		        endelse

		      filecount=filecount+1
		      names=[names, file_name+' component='+strtrim(string(ii),2)]

		      endfor

		    endif

	  endif

	endfor

	names=names[1:filecount]

	ptarr=PTRARR(filecount, /ALLOCATE_HEAP)

	for i=0, filecount-1 do *ptarr[i] = images[*,*,*,i]

	if n_elements(aspect_ratio) eq 0 then aspect_ratio=1.

		if KEYWORD_SET(roi) then begin

		  resolve_routine, 'roi3d_slicer', /is_function
		  planes=roi3d_slicer(ptarr, DATA_NAMES=files, /modal, aspect_ratio=aspect_ratio)

		endif else begin

		  resolve_routine, 'slicer3'
		  slicer3, ptarr, DATA_NAMES=files

		endelse

; This is the endelse for the 'file' comparison
endelse

return, planes

End
