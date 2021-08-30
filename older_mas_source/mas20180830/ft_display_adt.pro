;; $Id$
;; Created by E. Ozarslan on 23 September 2003

; function that will calculate the orientation maps
function fibers_orient, fa, ev, fdim, pdim, ir, ip, congrid=congrid

  compile_opt idl2, hidden
  fa=smooth(fa,5)
  fa=1./(1.+20.*(1.-fa)^16)
  im=fltarr(ir, ip, 3)
  if keyword_set(congrid) then begin
    im[*,*,0]=congrid(fa*abs(ev[*,*,2]), ir, ip, cubic=-.5)
    im[*,*,1]=congrid(fa*abs(ev[*,*,1]), ir, ip, cubic=-.5)
    im[*,*,2]=congrid(fa*abs(ev[*,*,0]), ir, ip, cubic=-.5)
  endif else begin
    im[*,*,0]=rebin(fa*abs(ev[*,*,2]), ir, ip)
    im[*,*,1]=rebin(fa*abs(ev[*,*,1]), ir, ip)
    im[*,*,2]=rebin(fa*abs(ev[*,*,0]), ir, ip)
  endelse
  bim=bytscl(im)
  return, bim
end


; this function returns the 24-bit (3 dimensional) pixel values
; of a given grayscale slice according to the given color table.
function fibers_color_mapper, im, fdim, pdim, table

  compile_opt idl2, hidden
  color_im=bytarr(fdim,pdim,3)
  loadct, table, /silent
  tvlct, r, g, b, /get
  for j=0, pdim-1 do for i=0, fdim-1 do $
    color_im[i,j,*]=[r[im[i,j]],g[im[i,j]],b[im[i,j]]]
  return, color_im
end


; draws the orientation circle
function fibers_circle, fdim, pdim, sorient

  compile_opt idl2, hidden
  dx=floor(fdim/2)
  dy=floor(pdim/2)
  im=fltarr(fdim, pdim, 3)
  r=ceil(2./5.*min([fdim, pdim]))
  for j=0, r-1 do for i=-r+1, r-1 do begin
    if i^2+j^2 le r^2 then begin
      x=float(i)/float(r)
      y=float(j)/float(r)
      im[dx+i,dy+j,0]=sqrt(1.-x^2-y^2)
      im[dx+i,dy+j,1]=y
      im[dx+i,dy+j,2]=abs(x)
      im[dx+i,dy-j,0]=sqrt(1.-x^2-y^2)
      im[dx+i,dy-j,1]=y
      im[dx+i,dy-j,2]=abs(x)
    endif
  endfor
  im=bytscl(im)
  if sorient eq 0 then im=shift(im, 0,0,-1)
  if sorient eq 1 then im=reverse(shift(im, 0,0,1),3)

  return, im
end




pro ft_display_adt, sorient, sposition, new_dims=new_dims

  compile_opt idl2
  common share_fibers_main, pTHR, pFA, pS_0, pAvD, pDIR, loadedthr, loadeddir

;, roi_objects

  ;sposition=sposition-1
  case sorient of
  	0 : begin
  	  aver_d = reform((*pAvD)[sposition, *,*])
  	  s0 = reform((*pS_0)[sposition, *,*])
  	  fa = reform((*pTHR)[sposition, *,*])
      evecs=reform((*pDIR)[sposition, *,*,*])
  	end
  	1: begin
  	  aver_d = reform((*pAvD)[*, sposition, *])
  	  s0 = reform((*pS_0)[*, sposition, *])
  	  fa = reform((*pTHR)[*, sposition, *])
      evecs=reform((*pDIR)[*, sposition, *,*])
  	end
  	2: begin
  	  aver_d = reform((*pAvD)[*, *,sposition])
  	  s0 = reform((*pS_0)[*, *,sposition])
  	  fa = reform((*pTHR)[*, *,sposition])
      evecs=reform((*pDIR)[*, *,sposition,*])
  	end
  	else: return
  endcase

  sz_data = size(aver_d)
  xdim = sz_data[1]
  ydim = sz_data[2]
  ballsize=((xdim < ydim)/8)*2

  byt_aver_d = Bytscl(aver_d, max=0.0008)
  aver_d=0
  byt_S0 = Bytscl(S0)
  s0=0
  byt_fa = Bytscl(fa, max=1.0)


  if n_elements(new_dims) ne 0 then begin
    ir=new_dims[0]
    ip=new_dims[1]
    byt_im=bytarr(2*ir+ballsize, 2*ip+ballsize,3)
    byt_im[0:ir-1,ip+ballsize:2*ip+ballsize-1,*]= $
      			fibers_color_mapper(congrid(byt_S0, ir, ip, cubic=-.5), ir, ip, 8)
    byt_im[ir+ballsize:2*ir+ballsize-1,ip+ballsize:2*ip+ballsize-1,*]= $
  				fibers_color_mapper(congrid(byt_aver_d, ir, ip, cubic=-.5), ir, ip, 1)
    byt_im[ir+ballsize:2*ir+ballsize-1,0:ip-1,*]= $
      			fibers_color_mapper(congrid(byt_fa, ir, ip, cubic=-.5), ir, ip, 3)



    ;byt_im[0:ir-1,0:ip-1,*]=orient(fa, reform(evecs), xdim, ydim, ir, ip, /congrid)
  	byt_im[0:ir-1,0:ip-1,*]=fibers_orient(fa, reform(evecs), xdim, ydim, ir, ip, /congrid)

    byt_im[ir:ir+ballsize-1,ip:ip+ballsize-1,*]=255
    byt_im[ir:ir+ballsize-1,ip:ip+ballsize-1,*]=fibers_circle(ballsize, ballsize, sorient)
  endif else begin
    ir=xdim
    ip=ydim
    byt_im=bytarr(2*ir+ballsize, 2*ip+ballsize,3)
    byt_im[0:ir-1,ip+ballsize:2*ip+ballsize-1,*]=fibers_color_mapper(rebin(byt_S0, ir, ip), ir, ip, 8)
    byt_im[ir+ballsize:2*ir+ballsize-1,ip+ballsize:2*ip+ballsize-1,*]= $
  						fibers_color_mapper(rebin(byt_aver_d, ir, ip), ir, ip, 1)
    byt_im[ir+ballsize:2*ir+ballsize-1,0:ip-1,*]=fibers_color_mapper(rebin(byt_fa, ir, ip), ir, ip, 3)
    byt_im[0:ir-1,0:ip-1,*]=fibers_orient(fa, reform(evecs), xdim, ydim, ir, ip)
    byt_im[ir:ir+ballsize-1,ip:ip+ballsize-1,*]=255
    byt_im[ir:ir+ballsize-1,ip:ip+ballsize-1,*]=fibers_circle(ballsize, ballsize, sorient)
  endelse


device, decomposed=1
device, /bypass_translation
  window, 8, Xsize=4*(2*ir+ballsize), ysize=4*(2*ip+ballsize), Title=' S0, <D>, Orientation & FA'
  tv, rebin(byt_im,4*(2*ir+ballsize),4*(2*ip+ballsize),3), true=3
  ;tv, byt_im, true=3
  loadct, 0, /silent

End
