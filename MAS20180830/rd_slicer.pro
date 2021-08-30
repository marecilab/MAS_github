;; $Id$
;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; rd_slicer.pro
;
; The RD_SLICER procedure reads any number of *.flt and *.trs files
; into slicer3 for combined display. Each file must have the same
; size. This procedure can be called by supplying a string array of
; file names to be displayed. If no argument is given the program
; prompts the user interactively until the user presses cancel.
;
; SYNTAX
; 	rd_slicer, files [, rf=rf] [, aspect_ratio=aspect_ratio]
;					 [, length_limit=length_limit]
;
; ARGUMENTS
;  files
;	A string array of file names to be displayed.  The
;	filenames should include the pathname (if they are not
;	in the current directory) and .flt extension.
;
;  rf
;   An array of resolution factors. The data will be interpolated
;   to rf[0]*fdim, rf[1]*pdim, rf[2]*sdim. High rf values make the
;   tracts look thinner.
;
;  aspect_ratio
;    A floating point number indicating the resolution along z
;    direction to that along x & y directions. When set, the default
;    z% value in the View mode of slicer3 is changed.
;
;  length_limit
;    Sensors out all tracts that are shorter than length_limit
;    voxels long.
;
; Last modified on September 4, 2003 by E. Ozarslan
; HS - 2006115
; Adding support to display data from MAS.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function decode_trs, filename, rf, dims, minsteps

  compile_opt idl2, hidden

  restore, filename
  if matrix[0] eq dims[0] and matrix[1] eq dims[1] and matrix[2] $
  										eq dims[2] then begin
    im=bytarr(dims[0]*rf[0],dims[1]*rf[1],dims[2]*rf[2])
    for i=0, n_elements(start_points)-2 do begin
      if start_points[i+1]-start_points[i] ge minsteps then begin
        ax=floor(rf[0]*rrx[start_points[i]:start_points[i+1]-1])
        ay=floor(rf[1]*rry[start_points[i]:start_points[i+1]-1])
        az=floor(rf[2]*rrz[start_points[i]:start_points[i+1]-1])
        for j=0, n_elements(ax)-1 do im[ax[j],ay[j],az[j]]=255
      endif
    endfor
    pt_im=ptr_new(im, /no_copy)
    return, pt_im
  endif else return, 0
end



pro rd_slicer, files, rf=rf, aspect_ratio=aspect_ratio, length_limit=length_limit

	common share_fibers_main
  	compile_opt idl2

  if n_elements(rf) eq 0 then rf=[1,1,1]

  if n_elements(length_limit) eq 0 then length_limit=0


  if (n_elements(files) eq 0) then begin
    files=dialog_pickfile(/must_exist, filter=['*.*'],/multiple_files,$
		path=Get_Path, Title='Read *.flt file(s).')
  endif
  if files[0] eq '' then return


; This is to read data from MAS. It will only receive the '.trs' file and therefore assume that the
; other data is coming from pointers in the common_share.
if N_ELEMENTS(files) eq 1 then begin
; Now I will display data coming from MAS.

		;to FUCKING STOP THE MOTHER****** CRASHING!!!
		file_exist = FILE_SEARCH(files[0])
		file_exist_size = size(file_exist)

			if file_exist_size[0] eq 0 then begin

				trs_file = DIALOG_PICKFILE( FILTER = '*.trs', Title = 'Select tracks file to open' ,$
                 path=path)
			endif

	restore, files[0]
    fdim=matrix[0]
    pdim=matrix[1]
    sdim=matrix[2]

  	images=bytarr(rf[0]*fdim,rf[1]*pdim,rf[2]*sdim,6)

	a3=Strpos(files[0], '.', /reverse_search)
    file_name=Strmid(files[0], 0, a3)
    if strlen(file_name) gt 35 then $
      file_name=strmid(file_name, 0, 10)+'...'+ strmid(file_name, 21, 22, /reverse_offset)

    pt_tract_im=decode_trs(files[0], rf, [fdim,pdim,sdim], length_limit)


	images[*,*,*,0]=*pt_tract_im
	ptr_free, pt_tract_im

	; pointers from MAS with the data
	; pTHR, pS_0, pAvD, pDIR
	im = *pTHR
	imagesnew=bytarr(rf[0]*fdim,rf[1]*pdim,rf[2]*sdim)
    imagesnew[*,*,*]=rebin(bytscl(im),rf[0]*fdim,rf[1]*pdim,rf[2]*sdim)
    images[*,*,*,1]=imagesnew
	im = *pS_0
	imagesnew=bytarr(rf[0]*fdim,rf[1]*pdim,rf[2]*sdim)
    imagesnew[*,*,*]=rebin(bytscl(im),rf[0]*fdim,rf[1]*pdim,rf[2]*sdim)
    images[*,*,*,2]=imagesnew
	im = *pAvD
	imagesnew=bytarr(rf[0]*fdim,rf[1]*pdim,rf[2]*sdim)
    imagesnew[*,*,*]=rebin(bytscl(im),rf[0]*fdim,rf[1]*pdim,rf[2]*sdim)
    images[*,*,*,3]=imagesnew
	im = *pDIR
	imagesnew=bytarr(rf[0]*fdim,rf[1]*pdim,rf[2]*sdim)
    imagesnew[*,*,*]=rebin(bytscl(im),rf[0]*fdim,rf[1]*pdim,rf[2]*sdim)
    images[*,*,*,4]=imagesnew

	; HS20070903
	; Adding a new dataset with the orientation images.
	; Using *pDIR (eigenvectors) to calculate the orientation image
	;im2d=reform(im3d[*,*,round(slice[5])])
;			faimage = im2d
;			evecs = *pDIR
;			sz_adt = size(faimage)
;			slice_axis = 0
;			orientimage = fltarr(sz_adt[1],sz_adt[2],3)
;			tempevecs =	fltarr(sz_adt[1],sz_adt[2], 3)
;
;				tempfa = faimage
;				tempevecs = REFORM(evecs[*,*,round(slice[5]),*])
;	   			orientimage = orient(tempfa, tempevecs , sz_adt[1], sz_adt[2], sz_adt[1],sz_adt[2], slice_axis)
;				im2d = orientimage


	;Adding a new dataset that contain the drawn ROIs in 3D space


	ROI_array = seeds_to_3Dspace(ROI_seeds_global, fdim,pdim,sdim)

	ROI_array = CONGRID(ROI_array, rf[0]*fdim,rf[1]*pdim,rf[2]*sdim)
	images[*,*,*,5] = ROI_array

	names=['Tracts', 'Fractional Anisotropy', 'S_0', 'Average Diffusivity', 'Eigenvectors', 'ROIs']

  	ptarr=PTRARR(6, /ALLOCATE_HEAP)

  	for i=0, 5 do *ptarr[i] = images[*,*,*,i]


endif else begin

  sz_files=size(files)
  num_files=sz_files[1]
  extension=strmid(files[0], strlen(files[0])-4, 4)

  if extension eq '.trs' then begin

    restore, files[0]
    fdim=matrix[0]
    pdim=matrix[1]
    sdim=matrix[2]

  endif else begin

  	if extension eq '.flt' then begin
    sz=size(rflt(files[0]))
    fdim=sz[1]
    pdim=sz[2]
    sdim=sz[3]
  	endif else return

  endelse



  filecount=0
  names=['a']

  images=bytarr(fdim,pdim,sdim,1)

  for fn=0, num_files-1 do begin

    extension=strmid(files[fn], strlen(files[fn])-4, 4)

    if extension eq '.trs' then begin

      a3=Strpos(files[fn], '.', /reverse_search)
      file_name=Strmid(files[fn], 0, a3)

      if strlen(file_name) gt 35 then $
        file_name=strmid(file_name, 0, 10)+'...' $
        +strmid(file_name, 21, 22, /reverse_offset)

      pt_tract_im=decode_trs(files[fn], rf, [fdim,pdim,sdim], length_limit)

      if ptr_valid(pt_tract_im) then begin

        if fn eq 0 then images=*pt_tract_im else begin
          imagesnew=bytarr(rf[0]*fdim,rf[1]*pdim,rf[2]*sdim,filecount+1)
          imagesnew[*,*,*,0:filecount-1]=images
          imagesnew[*,*,*,filecount]=*pt_tract_im
          images=imagesnew
        endelse

        ptr_free, pt_tract_im
        filecount=filecount+1
        names=[names, file_name]

      endif

    endif

    if extension eq '.flt' then begin

      im=rflt(files[fn])
      sz_im=size(im)

      a3=Strpos(files[fn], '.', /reverse_search)
      file_name=Strmid(files[fn], 0, a3)

      if strlen(file_name) gt 35 then $
        file_name=strmid(file_name, 0, 10)+'...' $
        +strmid(file_name, 21, 22, /reverse_offset)

      if sz_im[0] eq 3 or sz_im[0] eq 4 then $

      if sz_im[1] eq fdim and sz_im[2] eq pdim and sz_im[3] eq sdim then begin

        if sz_im[0] eq 3 then begin
          if fn eq 0 then images=rebin(bytscl(im),rf[0]*fdim,rf[1]*pdim,rf[2]*sdim) else begin
            imagesnew=bytarr(rf[0]*fdim,rf[1]*pdim,rf[2]*sdim,filecount+1)
            imagesnew[*,*,*,0:filecount-1]=images
            imagesnew[*,*,*,filecount]=rebin(bytscl(im),rf[0]*fdim,rf[1]*pdim,rf[2]*sdim)
            images=imagesnew
          endelse
          filecount=filecount+1
          names=[names, file_name]
        endif
        if sz_im[0] eq 4 then begin
          for ii=0,sz_im[4]-1 do begin
            if fn+ii eq 0 then images=rebin(bytscl(reform(im[*,*,*,ii])), $
            							rf[0]*fdim,rf[1]*pdim,rf[2]*sdim) else begin
              imagesnew=bytarr(rf[0]*fdim,rf[1]*pdim,rf[2]*sdim,filecount+1)
              imagesnew[*,*,*,0:filecount-1]=images
              imagesnew[*,*,*,filecount]=rebin(bytscl(reform(im[*,*,*,ii])), $
              								rf[0]*fdim,rf[1]*pdim,rf[2]*sdim)
              images=imagesnew
            endelse
            filecount=filecount+1
            names=[names, file_name+' component='+strtrim(string(ii),2)]
          endfor
        endif
      endif
    endif
  endfor

  names=names[1:filecount]

  ptarr=PTRARR(filecount, /ALLOCATE_HEAP)
  for i=0, filecount-1 do *ptarr[i] = images[*,*,*,i]


; This is the end else for the number of files comparison.
endelse



; HS - 20061031 (BOO!)
; Slicer3m is an edited version of slicer3 by Evren
; It takes the same parameters as the standard slicer3 and the aspect_ratio parameter.
; This aspect_ratio forces the Z axis slider to a specific number to keep the aspect ratio of the original image.
; See the source code for more information.

;These are here for testing only. Ignore.
;xvolume, *ptarr[0]
;xvolume, *ptarr[1]
;xvolume, *ptarr[2]

  resolve_routine, 'slicer3m'
  if n_elements(aspect_ratio) eq 0 then $
    slicer3m, ptarr, DATA_NAMES=names else $
    slicer3m, ptarr, DATA_NAMES=names, Aspect_ratio=aspect_ratio, SELECTED_SLICES=selectedslices


End
