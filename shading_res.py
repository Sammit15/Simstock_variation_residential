import sys
import random as rd
from eppy import modeleditor
from eppy.modeleditor import IDF

iddfile = "C:/Users/Anmol/Desktop/UEM__6_amd_res/EnergyPlus/ep8.9_windows/Energy+.idd"

# fname1 = "C:/Users/Home/Desktop/session_9_downloads/idf_files/bi_1.idf"

IDF.setiddname(iddfile)
#idf1 = IDF(fname1)


residential = [ "35_00257_0017","37_00085_0003","37_00085_0004","37_00085_0006","37_00247_0002","37_00247_0003","37_00247_0004","37_00247_0005","37_00247_0006","37_00247_0007","37_00247_0008","37_00247_0010","37_00247_0011","37_00247_0012","37_00247_0013","37_00247_0014","37_00247_0015","37_00247_0016","37_00247_0018","37_00247_0019","37_00247_0020","37_00247_0022","37_00247_0023","37_00247_0025","37_00247_0026","37_00247_0028","37_00247_0030","37_00248_0002","37_00248_0003","37_00248_0004","37_00248_0005","37_00248_0006","37_00248_0007","37_00248_0008","37_00248_0010","37_00248_0011","37_00248_0012","37_00248_0013","37_00248_0014","37_00248_0015","37_00248_0016","37_00248_0018","37_00248_0020","37_00248_0021","37_00248_0022","37_00248_0023","37_00248_0114","37_00248_1016","37_00249_0001","37_00249_0004","37_00250_0001","37_00250_0002","37_00250_0003","37_00250_0004","37_00250_0005","37_00250_0006","37_00251_0002","37_00251_0003","37_00251_0006","37_00251_0007","37_00251_0008","37_00251_0009","37_00251_0010","37_00251_0011","37_00251_0012","37_00251_0013","37_00251_0014","37_00251_0016","37_00251_0017","37_00251_0018","37_00251_0020","37_00251_0022","37_00251_0024","37_00251_0025","37_00251_1025","37_00252_0001","37_00252_0002","37_00252_0003","37_00252_0004","37_00252_0005","37_00252_0006","37_00252_0007","37_00252_0008","37_00252_0009","37_00252_0011","37_00252_0012","37_00252_0013","37_00252_0014","37_00252_0015","37_00252_0016","37_00252_0017","37_00252_0018","37_00252_0019","37_00252_0020","37_00253_0002","37_00253_0003","37_00253_0004","37_00253_0005","37_00253_0006","37_00253_0007","37_00253_0009","37_00253_0010","37_00257_0001","37_00257_0002","37_00257_0003","37_00257_0004","37_00257_0012","37_00257_0017","37_00257_0018","37_00257_0019","38_00247_0028","38_00248_0007","38_00248_0013","38_00248_0015","38_00250_0006","38_00252_0001","38_00257_0001","38_00257_0003","38_00257_0018","38_00257_0019","39_00250_0006","39_00257_0001"]
 
education = [ "37_00247_0029","37_00085_0002","37_00085_0005","37_00085_0007","37_00247_0001","37_00247_0009","37_00248_0001","37_00248_1001","37_00249_0002","37_00249_0003","37_00249_0005","37_00251_0027","38_00085_0002","38_00251_0027","39_00085_0002"]

for i in range(0,137):
    fname1 = "C:/Users/Anmol/Desktop/UEM__6_amd_res/idf_files/bi_{}.idf".format(i+1)
    #print(fname1)
    idf1 = IDF(fname1)
    window_list = idf1.idfobjects['WINDOW']
    window_names = list()
    for item in window_list:
        window_names.append(item.Name)
    
    #print(window_names)

    residential_win = list()
    for item in window_names:
        if item[0:13] in residential:
            residential_win.append(item)
    #print(residential_win)
    
    for win in residential_win:
        idf1.newidfobject(
        "Shading:Overhang".upper(),
        Name = "Overhang_{}".format(win),
        Window_or_Door_Name = win, 
        Height_above_Window_or_Door = 0,
        Tilt_Angle_from_WindowDoor = 90,
        Left_extension_from_WindowDoor_Width = 0,
        Right_extension_from_WindowDoor_Width = 0,
        Depth = 0.1 * rd.randint(3,6))
		
    education_win = list()
    for item in window_names:
        if item[0:13] in education:
            education_win.append(item)
			
    for win in education_win:
        idf1.newidfobject(
        "Shading:Overhang".upper(),
        Name = "Overhang_{}".format(win),
        Window_or_Door_Name = win, 
        Height_above_Window_or_Door = 0,
        Tilt_Angle_from_WindowDoor = 90,
        Left_extension_from_WindowDoor_Width = 0,
        Right_extension_from_WindowDoor_Width = 0,
        Depth = 0.1 * rd.randint(3,6))
		
		
		
    idf1.saveas("C:/Users/Anmol/Desktop/UEM__6_amd_res/idf_files_overhang/bi_{}_overhang.idf".format(i+1))
