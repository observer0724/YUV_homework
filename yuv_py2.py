import cv2
import numpy as np
from matplotlib import pyplot as plt



def BGR2YUV_422(Image):
    Y = np.ndarray(shape = (np.shape(Image)[0],np.shape(Image)[1]), dtype = 'float')
    U = np.ndarray(shape = (np.shape(Image)[0],np.shape(Image)[1]), dtype = 'float')
    V = np.ndarray(shape = (np.shape(Image)[0],np.shape(Image)[1]), dtype = 'float')
    matrix = np.matrix([[0.114,0.587,0.299],[0.5,-0.331,-0.169],[-0.081,-0.419,0.5]])
    for I in range(np.shape(Image)[0]):
        for J in range(np.shape(Image)[1]):
            Image_matrix = np.matrix(Image[I,J])
            Y[I,J],U[I,J],V[I,J] = matrix*np.transpose(Image_matrix)+np.matrix([[0],[128],[128]])
    Half_U = np.ndarray(shape = (np.shape(Image)[0],np.shape(Image)[1]/2), dtype = 'float')
    Half_V = np.ndarray(shape = (np.shape(Image)[0],np.shape(Image)[1]/2), dtype = 'float')
    for I in range(np.shape(Half_U)[1]):
        Half_U[:,I] = U[:,I*2]
        Half_V[:,I] = V[:,I*2]
    return Y,Half_U,Half_V

def BGR2YUV_420(Image):
    Y = np.ndarray(shape = (np.shape(Image)[0],np.shape(Image)[1]), dtype = 'float')
    U = np.ndarray(shape = (np.shape(Image)[0],np.shape(Image)[1]), dtype = 'float')
    V = np.ndarray(shape = (np.shape(Image)[0],np.shape(Image)[1]), dtype = 'float')
    matrix = np.matrix([[0.114,0.587,0.299],[0.5,-0.331,-0.169],[-0.081,-0.419,0.5]])
    for I in range(np.shape(Image)[0]):
        for J in range(np.shape(Image)[1]):
            Image_matrix = np.matrix(Image[I,J])
            Y[I,J],U[I,J],V[I,J] = matrix*np.transpose(Image_matrix)+np.matrix([[0],[128],[128]])
    Half_U = np.ndarray(shape = (np.shape(Image)[0],np.shape(Image)[1]/2), dtype = 'float')
    Half_V = np.ndarray(shape = (np.shape(Image)[0],np.shape(Image)[1]/2), dtype = 'float')
    Qut_U = np.ndarray(shape = (np.shape(Image)[0]/2,np.shape(Image)[1]/2), dtype = 'float')
    Qut_V = np.ndarray(shape = (np.shape(Image)[0]/2,np.shape(Image)[1]/2), dtype = 'float')
    for I in range(np.shape(Half_U)[1]):
        Half_U[:,I] = U[:,I*2]
        Half_V[:,I] = V[:,I*2]
    for I in range(np.shape(Qut_U)[0]):
        Qut_U[I,:] = Half_U[I*2,:]
        Qut_V[I,:] = Half_V[I*2,:]

    return Y,Qut_U,Qut_V

def extend_UV_422(U):
    U_full = np.ndarray(shape = (np.shape(U)[0],np.shape(U)[1]*2), dtype = 'float')
    for I in range(np.shape(U)[1]):
        U_full[:,2*I] = U[:,I]
        U_full[:,2*I+1] = U[:,I]
    return U_full

def extend_UV_420(U):
    U_half = np.ndarray(shape = (np.shape(U)[0],np.shape(U)[1]*2), dtype = 'float')
    U_full = np.ndarray(shape = (np.shape(U)[0]*2,np.shape(U)[1]*2), dtype = 'float')
    for I in range(np.shape(U)[1]):
        U_half[:,2*I] = U[:,I]
        U_half[:,2*I+1] = U[:,I]
    for I in range(np.shape(U)[0]):
        U_full[2*I,:] = U_half[I,:]
        U_full[2*I+1,:] = U_half[I,:]
    return U_full

def limit(num):
    if num>255:
        num = 255
    elif num<0:
        num = 0
    else:
        pass
    return num

def YUV_4222BGR(Y,U,V):
    Image = np.ndarray(shape = (np.shape(Y)[0],np.shape(Y)[1],3), dtype = 'float')
    new_U = extend_UV_422(U)
    new_V = extend_UV_422(V)
    matrix = np.matrix([[1,1.77216,0.00099],[1,-0.3437,-0.71417],[1,-0.00093,1.401687]])
    for I in range(np.shape(Image)[0]):
        for J in range(np.shape(Image)[1]):
            point = np.matrix([[Y[I,J]],[new_U[I,J]-128],[new_V[I,J]-128]])
            point_rgb = matrix*point
            Image[I,J,0] = limit(point_rgb[0])
            Image[I,J,1] = limit(point_rgb[1])
            Image[I,J,2] = limit(point_rgb[2])
    return Image

def YUV_4202BGR(Y,U,V):
    Image = np.ndarray(shape = (np.shape(Y)[0],np.shape(Y)[1],3), dtype = 'float')
    new_U = extend_UV_420(U)
    new_V = extend_UV_420(V)
    matrix = np.matrix([[1,1.77216,0.00099],[1,-0.3437,-0.71417],[1,-0.00093,1.401687]])
    for I in range(np.shape(Image)[0]):
        for J in range(np.shape(Image)[1]):
            point = np.matrix([[Y[I,J]],[new_U[I,J]-128],[new_V[I,J]-128]])
            point_rgb = matrix*point
            Image[I,J,0] = limit(point_rgb[0])
            Image[I,J,1] = limit(point_rgb[1])
            Image[I,J,2] = limit(point_rgb[2])
    return Image
Image = cv2.imread('YUV1.png')[1:,:,:]

Y,U,V = BGR2YUV_420(Image)
recre_Image = YUV_4202BGR(Y,U,V)
cv2.imwrite('420.png', recre_Image)
