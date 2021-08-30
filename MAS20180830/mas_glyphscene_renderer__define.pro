;; $Id$
;;

pro mas_glyphscene_renderer::render,$
                           rendered_img=rendered_img,$
                           use_this_ptr=use_this_ptr

    vdim  = self.vdim
    vmult = self.vmult
    xdim  = self.xsize
    ydim  = self.ysize

    if (vdim eq 0 or vmult eq 0.0) then begin
        print, "Cannot render without vdim and vmult begin set!"
        return
    end

    obuf      = obj_new('idlgrbuffer', quality=2)
    omodel    = obj_new('idlgrmodel', depth_test_disable=2)
    surfmodel = obj_new('idlgrmodel', depth_test_disable=2)

    ;; now we're at the desired size.
    surfmodel -> scale, 1.0/vmult, 1.0/vmult, 1.0/vmult
    
    ;; now make sure we're not > 4096x4096
    rendered_img_dim = [xdim*vdim, ydim*vdim]
    maxdim = max(rendered_img_dim)
    maxind = where(rendered_img_dim eq maxdim, count)
    maxsz  = 4096.0  ;; The max size of IDLgrBuffer is 4096
    
    if (maxdim gt maxsz) then begin
        
        if (count eq 2) then begin
            maxind = reform(maxind[0])
        endif

        reduction = maxsz/rendered_img_dim[maxind]
        rendered_img_dim[1-maxind] = floor(rendered_img_dim[1-maxind] * reduction)
        vdim = (floor(vdim * reduction))[0]
        surfmodel->scale, reduction, reduction, reduction

        rendered_img_dim[maxind] = maxsz
        print, "Image dimensions capped to:", rendered_img_dim[0], rendered_img_dim[1]
        if (not self.skip_bufsize_warnings) then begin
            junk = dialog_message('Output image too large for render buffer. Final size will be: '+$
                                  strcompress(string(rendered_img_dim[0])+'wx'+string(rendered_img_dim[1])+'h.', $
                                              /remove_all), /center)
        endif
    endif
    
    oimg      = obj_new('idlgrimage', $ 
                        depth_test_disable=2, $
                        transform_mode=1,     $
                        interleave=self.interleave, $
                        alpha_channel=self.bg_image_opac, $
                        blend_function=[3,4], $
                        location=[0,0,0])
    
    oview     = obj_new('idlgrview', $
                        color=[0,0,0], $
                        zclip=[2.*vdim+1,-(2.*vdim+1)], $
                        units=0, $
                        depth_cue=[-vdim,2.*vdim], $
                        eye=vdim+5.*vdim);vdim+2.*vdim)
    
    olightmodel   = obj_new('idlgrmodel')
    olights       = objarr(3)
    li            = 0
    olights[li++] = obj_new('idlgrlight', type=2, intensity=self.light0_int, location=[1,1,1]) ; 0.65
    olights[li++] = obj_new('idlgrlight', type=0, intensity=self.light1_int, location=[1,1,1]) ; 0.5
    olights[li++] = obj_new('idlgrlight', type=1, intensity=self.light2_int, coneangle=90, location=[1,1,1]) ; 0.5
    
    omodel->add,      oimg
    omodel->add,      surfmodel
    
    olightmodel->add, olights
    oview->add,       omodel
    oview->add,       olightmodel
    
    tmp        = bytarr(xdim*vdim, ydim*vdim, 3)
    tmp[*,*,0] = congrid((*self.bg_image)[*,*,0], xdim*vdim, ydim*vdim, /center) 
    tmp[*,*,1] = congrid((*self.bg_image)[*,*,1], xdim*vdim, ydim*vdim, /center)
    tmp[*,*,2] = congrid((*self.bg_image)[*,*,2], xdim*vdim, ydim*vdim, /center)
    
    bg_image   = temporary(tmp)

    oview->setProperty, viewplane_rect=[0,0,vdim, ydim*vdim]
    obuf->setProperty, dimensions=[vdim, rendered_img_dim[1]]

    progressBar = obj_new('progressbar', color='red', $
                          text='Rendering glyph image...', $
                          title='Glyph Scene Renderer', /fast_loop)
    progressBar->start
    
    if ptr_valid(use_this_ptr) then begin
        prendered_img = use_this_ptr
    endif else begin
        prendered_img = ptr_new(bytarr(rendered_img_dim[0], rendered_img_dim[1], 3),$
                                /no_copy)
    endelse

    for x = 0, xdim-1 do begin

        scanmodel = obj_new('idlgrmodel', depth_test_disable=2)
        surfmodel->add, scanmodel

        oimg->setProperty, data=bg_image[x*vdim:(x+1)*vdim-1,*,*]
        
        good_glyphs = obj_valid(reform((*self.glyphs_xy)[x,*]))
        if (total(good_glyphs) ne 0) then begin
            scanmodel->add, (*self.glyphs_xy)[x,where(good_glyphs eq 1)];, /alias
        endif
        
        obuf->draw, oview
        obuf->getProperty, image_data=chunk

        obj_destroy, scanmodel

        (*prendered_img)[x*vdim:(x+1)*vdim-1, *, 0] = chunk[0,*,*]
        (*prendered_img)[x*vdim:(x+1)*vdim-1, *, 1] = chunk[1,*,*]
        (*prendered_img)[x*vdim:(x+1)*vdim-1, *, 2] = chunk[2,*,*]
        
        progressBar->update, float(x+1)/xdim * 100.0
        
        if (progressBar->checkcancel()) then begin
            obj_destroy, [omodel, obuf, oimg, olights, olightmodel, surfmodel, oview]
            obj_destroy, *self.glyphs_xy
            progressbar->Destroy
            ;if (not ptr_valid(use_this_ptr)) then begin
                ptr_free, prendered_img
            ;endif
            return
        endif
        
    endfor
    
    rendered_img = prendered_img

    obj_destroy, [omodel, obuf, oimg, olights, olightmodel, surfmodel, oview]
    
    progressBar->Destroy
   
end

function mas_glyphscene_renderer_getOptimalVdim, xdim, ydim, vdim

    vdim_in = vdim
    rendered_img_dim = [xdim*vdim_in, ydim*vdim_in]
    maxdim = max(rendered_img_dim)
    maxind = where(rendered_img_dim eq maxdim, count)
    maxsz  = 4096.0  ;; The max size of IDLgrBuffer is 4096

    if (maxdim gt maxsz) then begin
        
        if (count eq 2) then begin
            maxind = reform(maxind[0])
        endif

        reduction = maxsz/rendered_img_dim[maxind]
        rendered_img_dim[1-maxind] = floor(rendered_img_dim[1-maxind] * reduction)
        vdim_in = (floor(vdim_in * reduction))[0]
        rendered_img_dim[maxind] = maxsz

    endif

    return, vdim_in

end


pro mas_glyphscene_renderer::addGlyph, glyph, xpos, ypos, replace=replace

    if (obj_valid((*self.glyphs_xy)[xpos,ypos])) then begin
        if not keyword_set(replace) then return
    endif

    (*self.glyphs_xy)[xpos,ypos] = glyph

end


pro mas_glyphscene_renderer::getProperty, $
                           bg_image=bg_image, $
                           xsize=xsize, $
                           ysize=ysize, $
                           vdim=vdim, $
                           vmult=vmult, $
                           light0_int=light0_int, $
                           light1_int=light1_int, $
                           light2_int=light2_int, $
                           bg_img_opac=bg_img_opac, $
                           glyph_opac=glyph_opac, $
                           skip_bufsize_warnings=skip_bufsize_warnings

    if (arg_present(xsize)) then xsize=self.xsize
    if (arg_present(ysize)) then ysize=self.ysize
    if (arg_present(vdim)) then vdim=self.vdim
    if (arg_present(vmult)) then vmult=self.vmult
    if (arg_present(light0_int)) then light0_int=self.light0_int
    if (arg_present(light1_int)) then light1_int=self.light1_int
    if (arg_present(light2_int)) then light2_int=self.light2_int
    if (arg_present(bg_img_opac)) then bg_img_opac=self.bg_img_opac
    if (arg_present(glyph_opac)) then glyph_opac=self.glyph_opac
    if (arg_present(skip_bufsize_warnings)) then skip_bufsize_warnings=self.skip_bufsize_warnings
    if (arg_present(bg_image)) then bg_image=*self.bg_image

end

pro mas_glyphscene_renderer::setProperty, $
                           bg_image=bg_image, $
                           vdim=vdim, $
                           vmult=vmult, $
                           light0_int=light0_int, $
                           light1_int=light1_int, $
                           light2_int=light2_int, $
                           bg_img_opac=bg_img_opac, $
                           glyph_opac=glyph_opac, $
                           skip_bufsize_warnings=skip_bufsize_warnings

    if (arg_present(vdim)) then self.vdim=vdim
    if (arg_present(vmult)) then self.vmult=vmult
    if (arg_present(light0_int)) then self.light0_int=light0_int
    if (arg_present(light1_int)) then self.light1_int=light1_int
    if (arg_present(light2_int)) then self.light2_int=light2_int
    if (arg_present(bg_img_opac)) then self.bg_img_opac=bg_img_opac
    if (arg_present(glyph_opac)) then self.glyph_opac=glyph_opac
    if (arg_present(skip_bufsize_warnings)) then self.skip_bufsize_warnings=skip_bufsize_warnings
    if (arg_present(bg_image)) then self.bg_image=*bg_image


end

pro mas_glyphscene_renderer::setBGImage, bg_image

    if (n_elements(bg_image) eq 0) then return

    img_dims = size(bg_image, /dimensions)

    prev_bg_image = self.bg_image

    if n_elements(img_dims) eq 3 then begin
        if (img_dims[0] eq 3) then begin
            self.xsize = img_dims[1]
            self.ysize = img_dims[2]
            tmp = bytarr(img_dims[1], img_dims[2], 3)
            tmp[*,*,0] = bg_image[0,*,*]
            tmp[*,*,1] = bg_image[1,*,*]
            tmp[*,*,2] = bg_image[2,*,*]
            self.bg_image = ptr_new(tmp, /no_copy)
        endif else if (img_dims[2] eq 3) then begin
            self.xsize = img_dims[0]
            self.ysize = img_dims[1]
            self.bg_image = ptr_new(bg_image)
        endif else begin
            message, 'bg_image must be either an 2D array, 3 by M by N array, or an M by N by 3 array'
            return
            ;; ERROR WHAT KIND OF IMAGE IS THIS?
        endelse

    endif else if (n_elements(img_dims) eq 2) then begin
        tmp = bytarr(img_dims[0], img_dims[1],3)
        tmp[*,*,0] = bg_image
        tmp[*,*,1] = bg_image
        tmp[*,*,2] = bg_image
        self.bg_image = ptr_new(tmp, /no_copy)
        self.xsize = img_dims[0]
        self.ysize = img_dims[1]
    endif else begin
        message, 'bg_image must be either an 2D array, 3 by M by N array, or an M by N by 3 array'
        return
    endelse

    if (ptr_valid(prev_bg_image)) then ptr_free, prev_bg_image

end

pro mas_glyphscene_renderer::cleanup

    if (ptr_valid(self.bg_image)) then ptr_free, self.bg_image

    if (ptr_valid(self.glyphs_xy)) then begin
        obj_destroy, *self.glyphs_xy
        ptr_free, self.glyphs_xy
    endif

end

function mas_glyphscene_renderer::init, $
                                bg_image=bg_image, $
                                xsize=xsize, $
                                ysize=ysize, $
                                vdim=vdim, $
                                vmult=vmult, $
                                light0_int=light0_int, $
                                light1_int=light1_int, $
                                light2_int=light2_int, $
                                bg_img_opac=bg_img_opac, $
                                glyph_opac=glyph_opac, $
                                skip_bufsize_warnings=skip_bufsize_warnings

    if (n_elements(bg_image) ne 0) then begin
        if arg_present(xsize) or arg_present(ysize) then begin
            ;; shouldn't have both bg_image and {x,y}size
        endif
        self->setBGImage, bg_image
    endif else begin
        message, 'argument bg_image is required'
    endelse

    self.glyphs_xy = ptr_new(objarr(self.xsize, self.ysize))
    self.interleave = 2

    if (n_elements(vdim) eq 0) then begin

    endif else begin
        self.vdim = vdim
    endelse

    if (n_elements(vmult) eq 0) then begin

    endif else begin
        self.vmult = float(vmult)
    endelse

    if (n_elements(light0_int) eq 0) then begin

    endif else begin
        self.light0_int = float(light0_int)
    endelse

    if (n_elements(light1_int) eq 0) then begin

    endif else begin
        self.light1_int = float(light1_int)
    endelse

    if (n_elements(light2_int) eq 0) then begin

    endif else begin
        self.light2_int = float(light2_int)
    endelse

    if (n_elements(bg_img_opac) eq 0) then begin

    endif else begin
        self.bg_image_opac = float(bg_img_opac)
    endelse

    if (n_elements(glyph_opac) eq 0) then begin

    endif else begin
        self.glyph_opac = float(glyph_opac)
    endelse

    if (n_elements(skip_bufsize_warnings) eq 0) then begin

    endif else begin
        self.skip_bufsize_warnings = skip_bufsize_warnings
    endelse

    return, 1

end


pro mas_glyphscene_renderer__define

    struct = { MAS_GLYPHSCENE_RENDERER, $
               xsize: 0L, $
               ysize: 0L, $
               vdim: 0L, $
               vmult: 0.0, $
               bg_image: ptr_new(), $
               mask: ptr_new(), $
               interleave:0, $
               glyphs_xy: ptr_new(), $
               light0_int: 0.0, $
               light1_int: 0.0, $
               light2_int: 0.0, $
               bg_image_opac: 0.0, $
               glyph_opac: 0.0, $
               skip_bufsize_warnings: 0 }

end
