import argparse, sys, os
import glob
import shutil
import numpy as np
import configparser
from PIL import Image
from scipy.spatial.transform import Rotation as R
#import numpy as np
home_path = os.getenv("HOME")

parser = argparse.ArgumentParser(description='Predict depth, filter, and fuse')
parser.add_argument('--scene_folder', default='{}/stereo_pairs'.format(home_path), help='select folder with scene from blender')
parser.add_argument('--outdirr_img', default='vig_dka/images', help='select folder to store images')
parser.add_argument('--outdirr_cams', default='vig_dka/cams_1', help='select folder to store camera data')
parser.add_argument('--pair', default='vig_dka/pair.txt', help='select path to pairs.txt')
parser.add_argument('--number_of_vews', default='10', help='num vievs')
parser.add_argument('--depth_min', default='0.5', help='idn insert please something which represent min depth')
parser.add_argument('--depth_max', default='2.5', help='idn insert please something which represent max depth')
parser.add_argument('--some_strange_constant', default='100', help='type whatever you want')

args = parser.parse_args()


current_dir = os.getcwd()

List_of_files = sorted(glob.glob(args.scene_folder+'/'+'scene_*'))

for i in List_of_files:
    print(i)

# create an instance of configparser
config = configparser.ConfigParser()

# read the .cfg file
#config.read('sample.cfg')



# copy the images
for i, item in enumerate(List_of_files):

     file = glob.glob(item+'/'+'input_Cam000*')
     print(file)
     source_file = file[0]
     destination_file = current_dir+'/'+args.outdirr_img+'/'+'{:0>8}.jpg'.format(i)
     with Image.open(source_file) as img:
        img.save(destination_file)


# copy intrinsic and extrinsic parameters
T = np.array([[0, -1, 0],[0,0,1],[-1,0,0]])

#from.here

# for i, item in enumerate(List_of_files):
#     matrix=[]
#     tru_matrix=[]
#     my_dictionary = {}
#     quater={}
#
#     config.read('sample.cfg')
#     with open(item+'/'+'parameters.cfg') as f:
#         my_list = f.readlines()
#         for k, line in enumerate(my_list):
#
#             if "Matrix 4x4" in line:
#                 line = line.strip()
#                 name, value = line.split('=')
#                 value = value.replace('Matrix 4x4','').replace('<','').replace(')','').replace('(','').replace(',','')
#                 value = value.strip()
#                 matrix.append(value)
#
#                 for j in range(3):
#                     line = my_list[k+j+1].strip()
#                     value = line.replace('>','').replace('(','').replace(')','').replace(',','')
#                     value = value.strip()
#                     matrix.append(value)
#
#                 for j in range(4):
#                     temp = np.fromstring(matrix[j],sep=' ')
#                     tru_matrix.append(temp)
#                 tru_matrix=np.array(tru_matrix)
#
#                 tru_matrix[0:3,3]=np.dot(T,tru_matrix[0:3,3])*(-10)
#
#
#
#             #quaternion extraxtion
#             if "center_camera_matrix_world_decompose_rot" in line:
#
#                 line = line.strip()
#                 name, value = line.split('=',1)
#                 value = value.replace('Quaternion','').replace('<','').replace(')','').replace('(','').replace(',','').replace('>','')
#                 value=value.strip().split(' ')
#                 for item in value:
#                     axis, val = item.split('=')
#                     quater[axis]=float(val)
#
#
#             elif '=' in line:
#                 line = line.strip()
#                 #print(line)
#                 name, value = line.split('=',1)
#                 value = value.strip()
#                 my_dictionary[name]=value
#                 #print('AAAAAAAAAAA')
#


    #until here

    #wrighting data in proper format

for i, item in enumerate(List_of_files):


    config.read(item+'/'+'parameters.cfg')

    Tx=float(config['extrinsics']['center_cam_x_m'])
    Ty=float(config['extrinsics']['center_cam_y_m'])
    Tz=float(config['extrinsics']['center_cam_z_m'])

    Rx=float(config['extrinsics']['center_cam_rx_rad'])
    Ry=float(config['extrinsics']['center_cam_ry_rad'])
    Rz=float(config['extrinsics']['center_cam_rz_rad'])

    fx_px=config['intrinsics']['fx_px']
    cx=config['intrinsics']['image_resolution_x_px']
    cy=config['intrinsics']['image_resolution_y_px']

    location = current_dir+'/'+args.outdirr_cams+'/'+'{:0>8}_cam.txt'.format(i)
    #print(location)
    with open(location,'w') as file:




        #The rotation vector as provided by Blender was first transformed to a rotation matrix:
        r = R.from_euler('xyz', (57.2958*np.array([Rx,Ry,Rz])), degrees=True)
        matR = r.as_matrix()

        #Transpose the rotation matrix, to find matrix from the WORLD to BLENDER coordinate system:
        R_world2bcam = np.transpose(matR)

        #The matrix describing the transformation from BLENDER to CV/STANDARD coordinates is:
        R_bcam2cv = np.array([[1, 0, 0],
                              [0, -1, 0],
                              [0, 0, -1]])

        #Thus the representation from WORLD to CV/STANDARD coordinates is:
        R_world2cv = R_bcam2cv.dot(R_world2bcam)

        #The camera coordinate vector requires a similar transformation moving from BLENDER to WORLD coordinates:
        vectC=[Tx,Ty,Tz]
        T_world2bcam = -1 * R_world2bcam.dot(vectC)
        T_world2cv = R_bcam2cv.dot(T_world2bcam)

        T_world2cv=T_world2cv[:,np.newaxis]
        print(R_world2cv)
        print(T_world2cv)

        fin_mat=np.hstack((R_world2cv,T_world2cv))
        fin_mat=np.vstack((fin_mat,[0,0,0,1]))
        print(fin_mat)
        tru_matrix=fin_mat
        print('=========================================')


        file.write('extrinsic')
        file.write('\n')
        file.write(str(tru_matrix).replace(']','').replace('[','').strip())
        #for item in quter:
          #  file.write(str(item)+' ')


        file.write('\n\n')
        file.write('intrinsic')

        #intrinsic
        file.write('\n')
        file.write(fx_px+' 0.0 '+str(int(cx)/2))
        #file.write('1372.'+' 0.0 '+str(int(my_dictionary['image_resolution_x_px '])/2))
        file.write('\n')
        file.write('0.0 '+fx_px+' '+str(int(cy)/2))
        #file.write('0.0 '+'1372.'+' '+str(int(my_dictionary['image_resolution_x_px '])/2))
        file.write('\n')
        file.write('0.0 0.0 1.0')
        file.write('\n\n')
        file.write(args.depth_min+' '+args.depth_max)


 #create pair.txt

num_viewpoint =len(List_of_files)

with open(current_dir+'/'+args.pair,'w') as f:
    f.write(str(num_viewpoint))
    for i in range(num_viewpoint):
        f.write('\n')
        f.write(str(i)+'\n')
        f.write(args.number_of_vews+' ')
        for j in range(10):
            if (i+j+1<num_viewpoint):
                f.write(str(i+j+1)+' ')
                f.write(args.some_strange_constant+' ')
            else:
                f.write(str(i-j-1)+' ')
                f.write(args.some_strange_constant+' ')





