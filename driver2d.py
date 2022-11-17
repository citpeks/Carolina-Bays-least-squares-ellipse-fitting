# Python program to fit an ellipse to Carolina Bay rim coordinates.
# A file "driver2d-coords.txt" contains more than six lines, each with the
#   latitude and longitude of a point along the Carolina Bay rim.
# Antonio Zamora July 5, 2022
# 07/07/2022 - Added BOM handling
# 07/11/2022 - Corrected problem with lat./lon.
#              allowed comment lines starting with # or *
# DRIVER2D.PY reads 2D coordinates instead of latitude and longitude pairs
#  when a comment with *2D is included ahead of the data pairs
#  A title for the graph is specified with a line starting with *T= title
# 07/21/2022 - restricted azimuth calculation for 28>Lat<49 & -105<Lon<-66
# 08/10/2022 - calculated area
# 11/16/2022 - Updated azimuth calculation. 

import numpy as np
from ellipse import LsqEllipse
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# minlat and maxlon are global variables
minlat = ""  
maxlon = "" 
sw2d = 0   # 2D process indicator
calc_azimuth = 0  # =1 to calculate azimuth
title = ""

# * * * * * * *
# Read file with coordinates and return two lists containing the x and y data of the ellipse
# * * * * * * *
def read_coordinates():
    ellipse_x= []
    ellipse_y = []
    list1= []
    list2 = []
    global minlat
    global maxlon
    global sw2d
    global calc_azimuth
    global title

    file1 = open('driver2d-coords.txt', 'r', encoding='utf-8-sig') # remove BOM
    Lines = file1.readlines()
    for line in Lines:
      line = line.strip()    # remove trailing spaces and \r
      # print("line:", line)
      if line[0:3] == "*2D" or line[0:3] == "*2d" :   # check for 2D process indicator
        sw2d = 1
      if line[0:3] == "*T=" or line[0:3] == "*t=" :   # check for title line
        title = line[3:]
      # Skip blank lines and comments starting with * or #
      if not line.startswith('*') and not line.startswith('#') and not line == '':  
        line = line.split(",")
        # print("Lat: ", line[0], " Lon: ", line[1])  #print latitude and longitude
        if sw2d == 1:
          list1.append(line[1])   # plot (Y-coordinate)
          list2.append(line[0])   # plot (X-coordinate)
        else :
          list1.append(line[0])   # Latitude (Y-coordinate)
          list2.append(line[1])   # Longitude (X-coordinate)
      # end if not ...  
    file1.close()

    # Find minimum X and Y coordinates
    if sw2d == 0 :  # Lat./Lon. Coordinates were input    
      minlat = min(list1)  # minimum latitude (southmost)
      maxlat = max(list1)  # maximum latitude
      maxlon = max(list2)  # maximum longitude (westmost)
      minlon = min(list2)  # minimum longitude
      # print("minLat: ", minlat," maxlat: ",maxlat, " maxLon: ", maxlon, " minLon: ",minlon)
      # 28>Lat<49 & -105<Lon<-66  [coordinates are in the contiguous United States]
      if float(minlat) > 28.0 and float(maxlat) < 49.0 and \
        float(maxlon) > -105. and float(minlon)< -66.0 :
        calc_azimuth =1
      # end if float...
    else:  # plot coordinates were input
      minlat = min(list1,key=lambda x:float(x))  # minimum latitude (southmost)
      maxlat = max(list1,key=lambda x:float(x))  # maximum latitude
      maxlon = max(list2,key=lambda x:float(x))  # maximum longitude (westmost)
      minlon = min(list2,key=lambda x:float(x))  # minimum longitude
    # end else plot coordinates were input
        
    # print("cos(minlat)=", np.cos(np.radians(float(minlat))) )

    # Convert latitude and longitude to meters relative to minimum coordinates
    # This will place the ellipse in the first quadrant
    # One degree of latitude = 10,000,000 m/90 degrees = 111,111 meters/degree
    for j in list1:
      # print("latitude ", j)
      if sw2d == 0 :
        # Subtract latitude from minlat and convert to meters
        ellipse_y.append( (float(j) - float(minlat))*111111 ) 
        # print( (float(j) - float(minlat))*111111  )
      else :
        ellipse_y.append( (float(j) - float(minlat)) )

    # Process Longitude
    # the distance in meters between degrees of longitude depends 
    #    on the latitude: cos(minlat)*111111
    for j in list2:
      # print("longitude ", j)
      if sw2d == 0 :
        ellipse_x.append( (abs(float(maxlon)) - abs(float(j)))*111111*np.cos(np.radians(float(minlat)) ) ) 
      else :
        ellipse_x.append( (float(j) - float(minlon))  ) 
      # print( (abs(float(maxlon)) - abs(float(j)))*111111*np.cos(np.radians(float(minlat)) ) )

    return [ellipse_x, ellipse_y]

# * * * * * * *
# Make test ellipse
# * * * * * * *
def make_test_ellipse(center=[1, 1], width=1, height=.6, phi=3.14/5):
    """Generate Elliptical data with noise

    Parameters
    ----------
    center: list:float
        (<x_location>, <y_location>)
    width: float
        semimajor axis. Horizontal dimension of the ellipse (**)
    height: float
        semiminor axis. Vertical dimension of the ellipse (**)
    phi: float:radians
        tilt of the ellipse, the angle the semimajor axis
        makes with the x-axis

    Returns
    -------
    data:  list:list:float
        list of two lists containing the x and y data of the ellipse.
        of the form [[x1, x2, ..., xi],[y1, y2, ..., yi]]
    """
    t = np.linspace(0, 2*np.pi, 300)
    x_noise, y_noise = np.random.rand(2, len(t))

    ellipse_x = center[0] + width*np.cos(t)*np.cos(phi)-height*np.sin(t)*np.sin(phi) + x_noise/2.  # noqa: E501
    ellipse_y = center[1] + width*np.cos(t)*np.sin(phi)+height*np.sin(t)*np.cos(phi) + y_noise/2.  # noqa: E501

    return [ellipse_x, ellipse_y]

# * * * * * * *
# Main program
# * * * * * * *
if __name__ == '__main__':
    
    X1, Y1 = read_coordinates() 
#    X1, Y1 = make_test_ellipse()

    X = np.array(list(zip(X1, Y1)))
    reg = LsqEllipse().fit(X)
    center, semimajor, semiminor, phi = reg.as_parameters()

    # print("minlat=", minlat, " maxlon=", maxlon)
    # print(" MAIN sw2d=",sw2d)
    print(title)
    print("coefficients for ax**2 + 2bxy + cy**2 + 2dx + 2fy + g")
    a, b, c, d, f, g = reg.coeffs()
    print(f'a={a:.3f}, b={b:.3f}, c={c:.3f}, d={d:.3f}, f={f:.3f}, g={g:.3f}')
    
    print(f'center: {center[0]:.3f}, {center[1]:.3f}')
    if sw2d == 0 :
      # Calculate coordinates for center
      Latitude = float(minlat) + center[1]/111111
      Longitude = float(maxlon) + center[0]/(111111*np.cos(np.radians(float(minlat))))
      print(f'center Lat.,Lon. {Latitude:.6f}, {Longitude:.6f}')

    k1 = semimajor
    k2 = semiminor
    if k1 < k2:  # swap so major axis is larger than minor axis
      j = k2;
      k2 = k1
      k1 = j

    print(f'major: {2*k1:.3f} m')
    print(f'minor: {2*k2:.3f} m')
    print(f'area: {np.pi*k1*k2:,.1f} square meters')    
    print(f'phi: {phi:.3f} ({np.rad2deg(phi):.3f} degrees)')
    # calculate azimuth
    # phi = counterclockwise angle of rotation from the x-axis to the major-axis of the ellipse 
    if calc_azimuth == 1 :
      if phi < 0:       # phi is negative
        print("phi < 0")
        if semimajor < semiminor:     # case when major and minor axes were flipped.
          print("semimajor < semiminor")
          azrad = 3*np.pi/2 - abs(phi)         # 270 - |phi| (phi<0)
        else:
          print("semimajor > semiminor")  # untested
          azrad = np.pi - abs(phi)         # 180 - |phi| (phi<0)
      else:
        print("phi > 0")
        if semimajor < semiminor:     # case when major and minor axes were flipped.
          print("semimajor < semiminor")
          if phi > np.pi/2 :  # phi > 90
            azrad = 2*np.pi - phi      # 360 - phi
          else:
            azrad = np.pi - phi         # 180 - phi 
        else:
          print("semimajor > semiminor")  # untested
          azrad = np.pi + np.pi/2 - phi
      # end if phi < 0: 
      # convert radians to degrees 
      print(f'azimuth: {np.degrees(azrad):.2f}')
    # end if calc_azimuth == 1  

    fig = plt.figure(figsize=(6, 6))
    ax = plt.subplot()
    ax.axis('equal')
    ax.plot(X1, Y1, 'ro', zorder=1)
    ax.set_title(f'{title}\n major axis: {2*k1:.1f} m, minor axis: {2*k2:.1f} m')
    ax.plot(center[0], center[1], 'go', label='center')
    ellipse = Ellipse(
        xy=center, width=2*semimajor, height=2*semiminor, angle=np.rad2deg(phi),
        edgecolor='b', fc='None', lw=2, label='Fit', zorder=2
    )
    ax.add_patch(ellipse)

    plt.xlabel('$X_1$')
    plt.ylabel('$Y_1$')

    plt.legend()
    plt.show()
