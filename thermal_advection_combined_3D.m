%%The problem of 3D advection-diffusion case


function T_matrix=thermal_advection_combined_3D(Lx,Ly,Lz,xnode,ynode,znode,E_points,v_original,materials,Temp_Incoming,h,Temp_Conv_inside,...
       Rel_v_theta,Measure_Box_Sub,fileID_Temp_Red_Box,Tecplot_check,Graphic_profile,Grpahic_contour)

