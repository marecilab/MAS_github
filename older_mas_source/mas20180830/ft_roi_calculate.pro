;; $Id$
;;
; ------------------------------------------------------------------------------
;
; ft_roi_calculate.pro
;
;
; Created on October 10, 2002 by E. Ozarslan
; Last modified on September 9, 2003 by E. Ozarslan
;
; ------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Seeding function

function seeder, pv, w, sq

  compile_opt idl2, hidden

  sz_pv=size(*pv)
  n=sz_pv[1]		; number of voxels

  if w eq 'x' or w eq 'y' or w eq 'z' then begin

    case w of
  	    'x': begin
				ii=1
				jj=2
		 end
	    'y': begin
				ii=0
				jj=2
		 end
	    'z': begin
				ii=0
				jj=1
		 end
    endcase

    imin=min((*pv)[*,ii])
    imax=max((*pv)[*,ii])
    jmin=min((*pv)[*,jj])
    jmax=max((*pv)[*,jj])
    dimi=imax-imin+1
    dimj=jmax-jmin+1
    mask=bytarr(dimi, dimj)
    for i=0,dimi-1 do for j=0,dimj-1 do begin
      for t=0, n-1 do begin
	if (*pv)[t,ii] eq imin+i and (*pv)[t,jj] eq jmin+j then mask[i,j]=1
      endfor
    endfor
    lattice=make_array(ceil(dimi/sq)*ceil(dimj/sq),2)
    t=0
    for i=0,ceil(dimi/sq)-1 do for j=0,ceil(dimj/sq)-1 do begin
      si=i*sq
      sj=j*sq
      if mask[floor(si),floor(sj)] then begin
	lattice[t,*]=[imin+si,jmin+sj]
	t=t+1
      endif
    endfor
    seeds=make_array(t,3)

    case w of
	    'x':begin
      			seeds[*,0]=(*pv)[0,0]
  	        	seeds[*,1]=lattice[0:t-1,0]
      			seeds[*,2]=lattice[0:t-1,1]
		end
	    'y':begin
	      		seeds[*,1]=(*pv)[0,1]
			    seeds[*,0]=lattice[0:t-1,0]
	      		seeds[*,2]=lattice[0:t-1,1]
		end
	    'z':begin
	      		seeds[*,2]=(*pv)[0,2]
	  		    seeds[*,0]=lattice[0:t-1,0]
      			seeds[*,1]=lattice[0:t-1,1]
		end
    endcase
  endif else begin 							; oblique plane

    xmin=min((*pv)[*,0])
    ymin=min((*pv)[*,1])
    zmin=min((*pv)[*,2])
    xmax=max((*pv)[*,0])
    ymax=max((*pv)[*,1])
    zmax=max((*pv)[*,2])
    dimx=xmax-xmin+1
    dimy=ymax-ymin+1
    dimz=zmax-zmin+1
    mask=bytarr(dimx, dimy, dimz)
    for i=0,dimx-1 do for j=0,dimy-1 do for k=0,dimz-1 do begin
      for t=0,n-1 do begin
	if (*pv)[t,0] eq xmin+i and (*pv)[t,1] eq ymin+j and $
				(*pv)[t,2] eq zmin+k then mask[i,j,k]=1
      endfor
    endfor
    lattice=make_array(ceil(dimx/sq)*ceil(dimy/sq)*ceil(dimz/sq),3)
    t=0
    for i=0,ceil(dimx/sq)-1 do for j=0,ceil(dimy/sq)-1 do $
						for k=0,ceil(dimz/sq)-1 do begin
      si=i*sq
      sj=j*sq
      sk=k*sq
      if mask[floor(si),floor(sj),floor(sk)] then begin
	lattice[t,*]=[xmin+si,ymin+sj,zmin+sk]
	t=t+1
      endif
    endfor
    seeds=make_array(t,3)
    seeds[*,0]=lattice[0:t-1,0]
    seeds[*,1]=lattice[0:t-1,1]
    seeds[*,2]=lattice[0:t-1,2]
  endelse

  return, seeds
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Main function

pro ft_roi_calculate, localpath, files, ptr_seeds

compile_opt idl2

;common share_main, pFA, pS0, pAvD, pEiVec
;common share_fibers_main, pTHR, pFA,pS_0, pAvD, pDIR, loadedthr, loadeddir
common share_fibers_main

n=n_elements(files)
if n ne 0 then begin

  ; extracting the information pointed by the ptr_seeds
  pa=ptrarr(2, /ALLOCATE_HEAP)
  pa=*ptr_seeds[0]
  if total(ptr_valid(pa)) eq 0 then begin
    ok=dialog_message('No ROI selected.')
    return
  endif
  planes=*ptr_seeds[1]
  sz_planes=size(*ptr_seeds[1])
  if sz_planes[0] then num_planes=1 else num_planes=sz_planes[2]

  seeds=0
  cum_inc=0

  ; seed points are calculated and stored
  for p=0, num_planes-1 do begin
    if ptr_valid(pa[p]) then begin
      inc=0
      if planes[0,p] then w='x' else  $
  	  if planes[1,p] then w='y' else $
	  if planes[2,p] then w='z' else w='o'
      seedsp=seeder(pa[p], w, 1)
      sz_seedsp=size(seedsp)
      inc=sz_seedsp[1]
      if cum_inc eq 0 then seeds=seedsp else begin
	    seeds=[seeds, seedsp]
      endelse
      cum_inc=cum_inc+inc
    endif
  endfor

  sz_seeds=size(seeds)

 if sz_seeds[0] ne 0 then begin

  x0=seeds[*,0]
  y0=seeds[*,1]
  z0=seeds[*,2]
  npoints=n_elements(x0)

    names=['Fract_Anisotropy','Aver_D','S0','DC']
    names2=['Fractional Anisotropy','Average Diffusion','S0','Directional Coherence']
    if files[0] then begin
      out=fltarr(npoints)
      for i=0, npoints-1 do out[i]=(*pFA)[x0[i],y0[i],z0[i]]
      fa_mean=mean(out)
      fa_stdv=npoints gt 1 ? stdev(out) : 0.
    endif
    if files[1] then begin
      out=fltarr(npoints)
      for i=0, npoints-1 do out[i]=(*pAvD)[x0[i],y0[i],z0[i]]
      ad_mean=mean(out)*1.e4
      ad_stdv=npoints gt 1 ? (stdev(out)*1.e4): 0.
    endif
    if files[2] then begin
      out=fltarr(npoints)
      for i=0, npoints-1 do out[i]=(*pS_0)[x0[i],y0[i],z0[i]]
      s0_mean=mean(out)
      s0_stdv=npoints gt 1 ? stdev(out): 0.
    endif
    if files[3] then begin
      out=roi_dc(localpath+'Eigenvectors.flt', x0,y0,z0,npoints,sigma=1.0)
      dc_mean=out[0]
      dc_stdv=out[1]
    endif
 endif else return

text=[  localpath, ' ', ' Number of Voxels = '+string(npoints)]
if files[0] then text=[text, ' ','       Fractional Anisotropy',$
						string(fa_mean)+'    '+'+/-'+string(fa_stdv)]
if files[1] then text=[text, ' ','       Average Diffusion',$
						string(ad_mean)+'    '+'+/-'+string(ad_stdv)]
if files[2] then text=[text, ' ','       S0',$
						string(s0_mean)+'    '+'+/-'+string(s0_stdv)]
if files[3] then text=[text, ' ','       Directional Coherence',$
						string(dc_mean)+'    '+'+/-'+string(dc_stdv)]

r=dialog_message(text, /information, title='Results')

endif

end
