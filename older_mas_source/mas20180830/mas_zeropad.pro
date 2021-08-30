; Subroutine name: mas_zeropad
; Created by: Magdoom Kulam
; Calling Information:
;
; data - State1 data 
;   
; Bugs or Important Comments to Developers:
;
; Purpose of subroutine: To zeropad the raw data by a factor
;
; Editing Information:


pro mas_zeropad, data

common scan_data

N = size(*data)
dim = project.imndarray[project.ci].dimensions      ; Type of experiment (2D or 3D)

; Zeropadding the raw data depending on the type of experiment
case N[0] of

1 :  begin 
     ; 1D experiment
     pad_data = make_array(2*N[1],/complex, value = 0) 
     pad_data[floor(N[1]/2):N[1]+floor(N[1]/2)-1] = *data
     end
     
2 :  begin 
     case dim of
      2   : begin
            ; 2D single slice experiment
            pad_data = make_array(2*N[1],2*N[2],/complex, value = 0)
            pad_data[floor(N[1]/2):N[1]+floor(N[1]/2)-1,floor(N[2]/2):N[2]+floor(N[2]/2)-1] = *data
            end
     else : begin                         
            ; 1D arrayed experiment
            pad_data = make_array(2*N[1],N[2],/complex, value = 0)
            pad_data[floor(N[1]/2):N[1]+floor(N[1]/2)-1,*] = *data
            end      
     endcase
     
     end
     
3 :  begin
     case dim of
      2 : begin
          ; 2D Multi-slice experiment
          pad_data = make_array(2*N[1],2*N[2],N[3],/complex, value = 0)
          pad_data[floor(N[1]/2):N[1]+floor(N[1]/2)-1,floor(N[2]/2):N[2]+floor(N[2]/2)-1,*] = *data
          end
      3  : begin
           ; 3D experiment 
           pad_data = make_array(2*N[1],2*N[2],2*N[3],/complex, value = 0)
           pad_data[floor(N[1]/2):N[1]+floor(N[1]/2)-1,floor(N[2]/2):N[2]+floor(N[2]/2)-1,floor(N[3]/2):N[3]+floor(N[3]/2)-1] = *data
           end    
     endcase
     end
     
4 :  begin
     case dim of
      2 : begin
          ; 2D Multi-slice arrayed experiment
          pad_data = make_array(2*N[1],2*N[2],N[3],N[4],/complex, value = 0)
          pad_data[floor(N[1]/2):N[1]+floor(N[1]/2)-1,floor(N[2]/2):N[2]+floor(N[2]/2)-1,*,*] = *data
          end
      3 : begin
          ; 3D arrayed experiment
          pad_data = make_array(2*N[1],2*N[2],2*N[3],N[4],/complex, value = 0)
          pad_data[floor(N[1]/2):N[1]+floor(N[1]/2)-1,floor(N[2]/2):N[2]+floor(N[2]/2)-1,floor(N[3]/2):N[3]+floor(N[3]/2)-1,*] = *data
          end   
     endcase      
     end
     
endcase

*data = pad_data 
end