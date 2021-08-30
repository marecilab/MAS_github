;; $Id$
;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; given two symmetric tensor fields, returns the tensor dot product
function tensordotp, pa, pb

  pp=ptr_new(/allocate_heap)
  *pp=(*(pa[0]))*(*(pb[0])) + (*(pa[1]))*(*(pb[1]))+ $
	(*(pa[2]))*(*(pb[2])) + 2*((*(pa[3]))*(*(pb[3]))+ $
	(*(pa[4]))*(*(pb[4]))+(*(pa[5]))*(*(pb[5])))
  return, pp
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function coherence, pD, pDyad, pDdotD, fdim, pdim, sdim

  a=fltarr(fdim,pdim,sdim)
  s=where((*pDdotD) ne 0., c)
  if c ne 0 then begin
    pdot=tensordotp(pD,pDyad)
    a[s]=((*pdot)[s])/((*pDdotD)[s]) > 0.
    a=sqrt(a/max(a))
    pcoh=ptr_new(a)
  endif else pcoh=ptr_new()
  ptr_free, pdot
  return, pcoh
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Given a symmetric tensor field in vector form, returns
; the principal eigenvector field.
function eigen, pt, pDdotD, fdim,pdim,sdim

  pv=ptrarr(3, /allocate_heap)
  *(pv[0])=fltarr(fdim,pdim,sdim)
  *(pv[1])=fltarr(fdim,pdim,sdim)
  *(pv[2])=fltarr(fdim,pdim,sdim)
  for k=0,sdim-1 do for j=0,pdim-1 do for i=0, fdim-1 do begin
    if ((*pDdotD)[i,j,k] ne 0.) then begin
      a=make_array(3,3)
      a=[[(*(pt[0]))[i,j,k], (*(pt[3]))[i,j,k], (*(pt[4]))[i,j,k]], $
         [(*(pt[3]))[i,j,k], (*(pt[1]))[i,j,k], (*(pt[5]))[i,j,k]], $
         [(*(pt[4]))[i,j,k], (*(pt[5]))[i,j,k], (*(pt[2]))[i,j,k]]]
      evalsctrho=eigenql(a, eigenvectors=evec)
      (*(pv[0]))[i,j,k]=evec[0,0]
      (*(pv[1]))[i,j,k]=evec[1,0]
      (*(pv[2]))[i,j,k]=evec[2,0]
    endif
  endfor
  return, pv
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Given a symmetric tensor field, returns fractional anisotropy
function fa, pten, fdim, pdim, sdim

  pfa=ptr_new(/allocate_heap)
  *pfa=fltarr(fdim,pdim,sdim)
  I1=(*(pten[0]))+(*(pten[1]))+(*(pten[2]))
  I2=(*(pten[0]))*((*(pten[1]))+(*(pten[2])))+(*(pten[1]))*(*(pten[2])) $
	-(*(pten[3]))^2-(*(pten[4]))^2-(*(pten[5]))^2
  den=I1^2-2*I2
;  *pfa = (den ne 0.) * sqrt(1-I2/den)
  s=where(den ne 0., c)
  if c ne 0 then (*pfa)[s]=sqrt(1-I2[s]/den[s]) < 1.
  return, pfa
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; given a vector field, returns the dyadic tensor field
function dyad, pv

  pdyad=ptrarr(6, /allocate_heap)
  *(pdyad[0])=(*(pv[0]))^2
  *(pdyad[1])=(*(pv[1]))^2
  *(pdyad[2])=(*(pv[2]))^2
  *(pdyad[3])=(*(pv[0]))*(*(pv[1]))
  *(pdyad[4])=(*(pv[0]))*(*(pv[2]))
  *(pdyad[5])=(*(pv[1]))*(*(pv[2]))
  return, pdyad
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Constructs the coherence vector with correct signs
function cohvec, pDas, pDxy, pDxz, pDyz

  pcv=ptrarr(3, /allocate_heap)
  ker_x=fltarr(3,3,3)
  ker_x[0,1,1]=.5
  ker_x[2,1,1]=.5
  ker_y=fltarr(3,3,3)
  ker_y[1,0,1]=.5
  ker_y[1,2,1]=.5
  ker_z=fltarr(3,3,3)
  ker_z[1,1,0]=.5
  ker_z[1,1,2]=.5
  *(pcv[0])=convol(*(pDas[0]), ker_x, /center)
  *(pcv[1])=convol(*(pDas[1]), ker_y, /center)
  *(pcv[2])=convol(*(pDas[2]), ker_z, /center)
  t1=(*pDxy)*(*(pcv[0]))*(*(pcv[1]))
  t2=(*pDxz)*(*(pcv[0]))*(*(pcv[2]))
  t3=(*pDyz)*(*(pcv[1]))*(*(pcv[2]))
  t_ppm=t1-t2-t3
  t_pmp=-t1+t2-t3
  t_mpp=-t1-t2+t3
  t_ppp=t1+t2+t3
  c1=t_ppm gt t_pmp
  c2=t_ppm gt t_mpp
  c3=t_ppm gt t_ppp
  c4=t_pmp gt t_mpp
  c5=t_pmp gt t_ppp
  c6=t_mpp gt t_ppp
  s=where(c1+c2+c3 eq 3 , c)
  if c ne 0 then (*(pcv[2]))[s]=-((*(pcv[2]))[s])
  s=where(1-c1+c4+c5 eq 3 , c)
  if c ne 0 then (*(pcv[1]))[s]=-((*(pcv[1]))[s])
  s=where(2-c2-c4+c6 eq 3 , c)
  if c ne 0 then (*(pcv[0]))[s]=-((*(pcv[0]))[s])
  return, pcv
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; blurs a 3-d, n component image on a component by component
; basis with a specified standard deviation

function gaussconv, pDa, comp, sigma

  if sigma eq 0. then return, pDa
  n=ceil(4.*sigma)
  sf=max([2*n+1, 3])
  result=0.
  kernel=make_array(sf,sf,sf)
  den=2*sigma^2

  for k3=0,sf-1 do for k2=0,sf-1 do for k1=0,sf-1 do $
      kernel(k1,k2,k3)=exp(-((k1-n)^2+(k2-n)^2+(k3-n)^2)/den)

  pDasigma=ptrarr(comp, /allocate_heap)
  kernel=kernel/total(kernel)

  for i=0, comp-1 do $
    *(pDasigma[i])=convol(*(pDa[i]), kernel, /center)

  return, pDasigma
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Main procedure
; Calculates coherence and related parameters from ADT.flt file
function roi_dc, file, x0,y0,z0,nn, sigma=sigma

  common share_main, pFA, pS0, pAvD, pEiVec
;load

  if not n_elements(file) then $
    file=dialog_pickfile(/must_exist, filter='*.flt',path=Get_Path, $
	file='Eigenvectors.flt', Title='Read in the Eigenvectors.flt file')
  if file eq '' then return, 0
  ;ev=rflt(file)
  ;ev=ev[*,*,*,*,0]
  ev=*pEiVec
  sz_data=size(ev)
  fdim=sz_data[1]		; number of points in read
  pdim=sz_data[2]		; number of points in phase
  sdim=sz_data[3]		; number of points in slice

n=ceil(4.*sigma)+1
xmin=min(x0)
ymin=min(y0)
zmin=min(z0)
xmax=max(x0)
ymax=max(y0)
zmax=max(z0)
xdim=xmax-xmin+2*n+1
ydim=ymax-ymin+2*n+1
zdim=zmax-zmin+2*n+1

mask=bytarr(xdim,ydim,zdim,3)
for i=0, nn-1 do mask[x0[i]-xmin+n,y0[i]-ymin+n,z0[i]-zmin+n,*]=[1,1,1]
ev2=fltarr(xdim,ydim,zdim,3)
ev2=ev*mask
ev=ev2

ev2=0
  pevx=ptr_new(ev[*,*,*,0])
  pevy=ptr_new(ev[*,*,*,1])
  pevz=ptr_new(ev[*,*,*,2])
  pev=[pevx,pevy,pevz]
  ev=0
;  pabs=ptr_new(/allocate_heap)
;  *pabs=sqrt(*pevx^2+*pevy^2+*pezz^2)

;calculate the dyadic field for the principal eigenvector
  pdyad=dyad(pEv)

;blur with sigma
  if n_elements(sigma) ne 1 then sigma=1.
  pTsigma=gaussconv(pdyad,6,sigma)

;calculate FA(Tsigma)
  pFA=fa(pTsigma, xdim, ydim, zdim)

dc=fltarr(nn)
for i=0, nn-1 do dc[i]=(*pFA)[x0[i]-xmin+n,y0[i]-ymin+n,z0[i]-zmin+n]

  ptr_free, pevx
  ptr_free, pevy
  ptr_free, pevz
  ptr_free, pev
  ptr_free, pFA
  ptr_free, pTsigma
  ptr_free, pDyad

  return, [mean(dc), stdev(dc)]
end
