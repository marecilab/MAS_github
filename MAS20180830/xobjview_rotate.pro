;$Id: //depot/idl/IDL_64/idldir/lib/utilities/xobjview_rotate.pro#1 $
;
;  Copyright (c) 1997-2007, ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
pro xobjview_rotate, $
    axis, $     ; IN: three-element vector about which to rotate.
    angle, $    ; IN: The amount (measured in degrees) of rotation
    tlb=tlb_, $ ; IN: (opt) Widget ID of an xobjview.  (Can be an array.)
    premultiply=premultiply ; IN: (opt) if set, do "data-centric" rotation.
;
;Purpose:  Provide an interface to the same behavior that occurs
;when a user rotates xobjview's graphic with the mouse.
;If input keyword TLB is not supplied, and if at least one instance of
;xobjview is currently running, then this routine will operate on
;the most recently created currently-running xobjview.
;
;Example:
;   xobjview, /test
;   for i=0,359 do begin
;       xobjview_rotate, [0,1,0], 1, /premult
;       xobjview_write_image, $
;           'img' + strcompress(i, /remove_all) + '.bmp', 'bmp'
;       end
;
on_error, 2
;
;Obtain valid TLB.
;
if n_elements(tlb_) gt 0 then begin
    if size(tlb_, /tname) ne 'LONG' then $
        message, 'TLB must be of type LONG.'
    if min(widget_info(tlb_, /valid_id)) eq 0 then $
        message, 'Invalid TLB.'
    for i=0,n_elements(tlb_)-1 do begin
        widget_control, tlb_[i], get_uvalue=pState
        if not ptr_valid(pState) then $
            message, 'Incorrect TLB.'
        if size(*pState, /tname) ne 'STRUCTURE' then $
            message, 'Incorrect TLB.'
        if max(tag_names(*pState) eq 'OOBJVIEWWID') eq 0 then $
            message, 'Incorrect TLB.'
        end
    tlb = tlb_
    end $
else begin
    if xregistered('xobjview') eq 0 then $
        message, 'No valid XOBJVIEW available.'
    tlb = LookupManagedWidget('xobjview')
    end
;
;Perform the desired rotation on each TLB.
;
for i=0,n_elements(tlb)-1 do begin
    widget_control, tlb[i], get_uvalue=pState
    (*pState).oObjViewWid->Rotate, axis, angle, premultiply=premultiply
    end

end
