;$Id: //depot/idl/IDL_64/idldir/lib/utilities/idlexinscribingview__define.pro#1 $
;
;  Copyright (c) 1997-2007, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;+
; NAME:
;   IDLexInscribingView
;
; PURPOSE:
;   Provide a type of IDLgrView that has convenient methods for
;   defining the view volume and eyepoint, and getting viewport
;   dimensions.  By default, this class's method to define the
;   view volume, "SetViewVolume", defines the view volume such
;   that the volume inscribes the view's contents.  Hence the
;   class name: "IDLexInscribingView."
;
;-
function IDLexInscribingView::init, $
    _extra=e

self.isotropic = 1b ; Default.
self.do_aspect = 1b ; Default.

return, self->IDLgrView::init(_extra=e)
end
;--------------------------------------------------------------------
pro IDLexInscribingView::SetProperty, $
    isotropic=isotropic, $
    do_aspec=do_aspect, $
    _extra=e
compile_opt hidden

if n_elements(do_aspect) eq 1 then $
    self.do_aspect = do_aspect

if n_elements(isotropic) eq 1 then begin
    self.isotropic = isotropic
    self.do_aspect = 1b
    end

self->IDLgrView::SetProperty, _extra=e
end
;--------------------------------------------------------------------
pro IDLexInscribingView::ZoomViewplaneRect, factor
compile_opt hidden

self->GetProperty, viewplane_rect=viewplane_rect
old_center = viewplane_rect[0:1] + viewplane_rect[2:3] / 2
new_x_size = viewplane_rect[2] / factor
new_y_size = viewplane_rect[3] / factor
self->SetProperty, viewplane_rect=[ $
    new_x_size / (-2) + old_center[0], $
    new_y_size / (-2) + old_center[1], $
    new_x_size, $
    new_y_size $
    ]
end
;--------------------------------------------------------------------
function IDLexInscribingView::GetViewportDimensions, $
    oDestDevice, $      ; IN
    inches=inches, $    ; IN: (opt) boolean. Return result in inches.
    location=location, $; OUT: (opt)
    defaulted=defaulted ; OUT: (opt) boolean.
compile_opt hidden

defaulted = 0b
self->GetProperty, dimensions=dimensions, location=location, units=units
if max(dimensions) eq 0 then begin
    oDestDevice->GetProperty, units=orig_dest_units
    oDestDevice->SetProperty, units=units
    oDestDevice->GetProperty, dimensions=dimensions
    oDestDevice->SetProperty, units=orig_dest_units

    location = [0., 0.]
    defaulted = 1b
    end

if keyword_set(inches) then begin
    case units of
        0 : begin ; Pixels.
            oDestDevice->GetProperty, $
                resolution=cpd ; Centimeters per dot.
            dpi = 2.54 / float(cpd) ; Dots per inch.
            dimensions = dimensions / dpi
            location = location / (dpi - 1.)
            end
        1 :       ; Inches.
        2 : begin ; Centimeters.
            dimensions = dimensions / 2.54
            location = location / 2.54
            end
        3 : begin ; Normalized.
            oDestDevice->GetProperty, units=orig_dest_units
            oDestDevice->SetProperty, units=1 ; Inches.
            oDestDevice->GetProperty, dimensions=dest_dimensions
            oDestDevice->SetProperty, units=orig_dest_units

            dimensions = dimensions * dest_dimensions
            location = location * dest_dimensions
            end
        endcase
    end

return, dimensions
end
;--------------------------------------------------------------------
pro IDLexInscribingView::PadViewplaneRect, oDestDevice
;
;Purpose: Grow self's viewplane_rect in x or y so as to match
;the aspect ratio of self's viewport as it appears on the given
;destination device.  (Doing this prevents self's contents from
;stretching in x or y when projected onto the device.)
;
;Find viewport's aspect ratio as it appears on the destination
;device.
;
viewport_inch_dims = self->GetViewportDimensions(oDestDevice, /inches)
viewport_aspect = viewport_inch_dims[0] / double(viewport_inch_dims[1])
;
;Find viewplane_rect's aspect ratio.
;
self->GetProperty, viewplane_rect=viewplane_rect
vrect_aspect = viewplane_rect[2] / viewplane_rect[3]
;
;Compensate for the difference in aspect ratios.
;
old_center = viewplane_rect[0:1] + viewplane_rect[2:3] / 2

if vrect_aspect gt viewport_aspect then begin ; grow in y.
    new_y_size = viewplane_rect[3] * vrect_aspect / viewport_aspect
    viewplane_rect[1] = new_y_size / (-2) + old_center[1]
    viewplane_rect[3] = new_y_size
    end $
else begin ; grow in x.
    new_x_size = viewplane_rect[2] * viewport_aspect / vrect_aspect
    viewplane_rect[0] = new_x_size / (-2) + old_center[0]
    viewplane_rect[2] = new_x_size
    end

self->SetProperty, viewplane_rect=viewplane_rect
end
;--------------------------------------------------------------------
pro IDLexInscribingView::SetEye, $
    degrees ; IN: (opt)  Requested amount of perspective effect.
            ;   Default is 60.
;
;Purpose: For perspective projections, position self's eye
;to achieve the requested amount of perspective effect, or the maximum
;amount of perspective effect possible (whichever is less).
;
;DEGREES indicates an amount of perspective effect, and is analogous to
;selecting a "fish-eye" lens for a camera, where you can have a 60 degree
;lens, a 170 degree lens, etc.  The greater the amount of degrees,
;the greater the perspective effect.
;
;Since IDL's eye position cannot be less than zero (in z), method
;SetEye cannot achieve the requested amount of effect in all cases.
;
;DEGREES must be greater than 0 and less than 180.
;
on_error, 2

if n_elements(degrees) eq 0 then $
    deg = 60 $
else $
    deg = degrees

if deg le 0 then $
    message, 'Argument must be greater than zero.'

if deg ge 180 then $
    message, 'Argument must be less than 180.'

small_fraction = .01 ; Arbitrary.
self->GetProperty, zclip=zclip, viewplane_rect=viewplane_rect
self->SetProperty, $
    eye = ( $
        zclip[0] + (viewplane_rect[3] / 2.) / tan(!dtor * deg * .5) $
        ) > (zclip[0] + small_fraction * (zclip[0] - zclip[1]))
end
;--------------------------------------------------------------------
pro IDLexInscribingView::SetViewVolume, $
    oDestDevice, $      ; IN
    scale=scale, $      ; IN: (opt) Scale size of inscribing view volume.
    xrange=xrange, $    ; IN: (opt) Inscribe this range.
    yrange=yrange, $    ; IN: (opt) Inscribe this range.
    zrange=zrange, $    ; IN: (opt) Inscribe this range.
    quiet=quiet, $      ; IN: (opt)
    isotropic=isotropic_; IN: (opt)
;
;Purpose: Calculate and set self's view volume.  The view volume will
;   be calculated so as to enclose the contents of the view.
;
;INPUT:
;  oDestDevice: A reference to a Destination Object (i.e. an
;     IDLgrWindow or an IDLgrPrinter).
;
; MODIFICATION HISTORY:
;  Written by: PCS 1999
;       Influenced by procedure "set_view" (DD, February 1997),
;       and procedure "create_view" (Daniel Carr 1992)
;
if n_elements(scale) eq 0 then begin
    scl = 1./sqrt(3)
    end $
else begin
    scl = scale
    end

if scl lt 0 then $
    message, 'SCALE must be a positive value.'
if scl eq 0 then begin
    scl = 1./sqrt(3)
    end
;
;Determine the overall bounding box.
;
if not self->GetBounds(oDestDevice, xrng, yrng, zrng) $
then begin
    xrng = [-1., 1.]
    yrng = [-1., 1.]
    zrng = [-1., 1.]
    end

if n_elements(xrange) gt 0 then xrng = double(xrange)
if n_elements(yrange) gt 0 then yrng = double(yrange)
if n_elements(zrange) gt 0 then zrng = double(zrange)
;
;If isotropy is to be maintained, use largest of three dimensions.
;
if n_elements(isotropic_) eq 0 then begin
    isotropic = keyword_set(self.isotropic)
    end $
else begin
    isotropic = keyword_set(isotropic_)
    end

if keyword_set(isotropic) then begin
    xlen = xrng[1] - xrng[0]
    ylen = yrng[1] - yrng[0]
    zlen = zrng[1] - zrng[0]
    maxlen = xlen > ylen > zlen
    xrng[0] = xrng[0] - ((maxlen - xlen) / 2.0)
    xrng[1] = xrng[0] + maxlen
    yrng[0] = yrng[0] - ((maxlen - ylen) / 2.0)
    yrng[1] = yrng[0] + maxlen
    zrng[0] = zrng[0] - ((maxlen - zlen) / 2.0)
    zrng[1] = zrng[0] + maxlen
    end

xlen = xrng[1] - xrng[0]
ylen = yrng[1] - yrng[0]
zlen = zrng[1] - zrng[0]
;
;Scale the view volume.
;
xpad = (xlen * (1./scl) - xlen) / 2.
ypad = (ylen * (1./scl) - ylen) / 2.
zpad = (zlen * (1./scl) - zlen) / 2.
xrng[0] = xrng[0] - xpad
xrng[1] = xrng[1] + xpad
yrng[0] = yrng[0] - ypad
yrng[1] = yrng[1] + ypad
zrng[0] = zrng[0] - zpad
zrng[1] = zrng[1] + zpad

xlen = xrng[1] - xrng[0]
ylen = yrng[1] - yrng[0]
zlen = zrng[1] - zrng[0]
;
;Calculate Viewplane Rectangle.
;
viewplane_rect = [ $
    xrng[0], $
    yrng[0], $
    xlen, $
    ylen $
    ]

zclip = [zrng[1] + .1 * zlen, zrng[0] - .1 * zlen]

if keyword_set(quiet) then begin
;
;   If self's eyepoint would be inside the new view volume, IDL will
;   issue a warning message and move the eyepoint when we set the new
;   view volume.  To avoid the warning message, ensure that the
;   eyepoint will be outside the new view volume.
;
    self->GetProperty, eye=eye
    if eye le zclip[0] then begin
        self->SetProperty, eye=zclip[0] + .1 * zlen
        end
    end
;
;Set the new view volume.
;
self->SetProperty, viewplane_rect=viewplane_rect, zclip=zclip
if keyword_set(isotropic) or self.do_aspect then begin
    self->PadViewplaneRect, oDestDevice
    end

end
;--------------------------------------------------------------------
function IDLexInscribingView::GetBounds, $
    oDestDevice, $  ; IN: Supplies GetTextDimensions() if needed.
    xrange, $       ; OUT
    yrange, $       ; OUT
    zrange, $       ; OUT
    ignore=ignore, $; IN: (opt) classes.
    skip=skip       ; IN: (opt) classes. (non-recursive).
compile_opt hidden

;
;   Determine the overall bounding box for objects contained in self.
;
oModels = self->Get(/all, count=n_models)
for i=0L, n_models-1 do begin
    if get_obj_range(oModels[i], self, oDestDevice, $
                      model_range, ignore=ignore, skip=skip) then begin
        if n_elements(range) ne 0 then $
          range = [[range[0:2] < model_range[0:2]], $
                   [range[3:5] > model_range[3:5]]] $
        else range = model_range
    endif
endfor

if n_elements(range) ne 0 then begin ;Extract individual ranges
    xrange = range[[0,3]]
    yrange = range[[1,4]]
    zrange = range[[2,5]]
    return, 1
endif else return, 0
end
;--------------------------------------------------------------------
pro IDLexInscribingView__define
compile_opt idl2, hidden

struct_hide, {IDLexInscribingView, $
    inherits IDLgrView, $
    isotropic: 0b, $
    do_aspect: 0b $
    }
end
