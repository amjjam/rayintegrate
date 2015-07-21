function loadTest,file

openr,un,file,/get_lun

ny=0l & nz=0l
readu,un,ny,nz

d=dblarr(nz,ny)

readu,un,d

free_lun,un

return,{nz:nz,ny:ny,d:d}

end

