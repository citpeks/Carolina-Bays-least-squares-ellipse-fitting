# Python program to fit an ellipse to Carolina Bay rim coordinates.
# A file "driver2d-coords.txt" contains more than four lines, each with the
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
# 08/30/2023 - Updated azimuth calculation phi<0 & semimajor < semiminor.
# 08/07/2024 - Added option to specify the input file name and performed some checks on input data
# 01/31/2025 - Calculate Goodness of Fit (Mean Squared Error MSE)
# 02/02/2025 - Added an option line *O= to activate DDT and GRID switches

import numpy as np
import math
from ellipse import LsqEllipse
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.optimize import fsolve
import os
import os.path

# These are global variables
minlat = ""  
maxlon = "" 
sw2d = 0   # 2D process indicator
calc_azimuth = 0  # =1 to calculate azimuth
title = ""

# Xc,Yc are the coordinates of the center of the ellipse, L = major axis, W = minor axis , phi= angle of rotation
def ellipse_parametric(t, Xc, Yc, L, W, phi):
    x = Xc + (L/2) * np.cos(t) * np.cos(phi) - (W/2) * np.sin(t) * np.sin(phi)
    y = Yc + (L/2) * np.cos(t) * np.sin(phi) + (W/2) * np.sin(t) * np.cos(phi)
    return x, y

def line_parametric(r, Xc, Yc, theta):
  x = Xc + r * np.cos(theta)
  y = Yc + r * np.sin(theta)
  return x, y


def equations(vars, Xc, Yc, L, W, phi, theta):
    t, r = vars
    x_ellipse, y_ellipse = ellipse_parametric(t, Xc, Yc, L, W, phi)
    x_line, y_line = line_parametric(r, Xc, Yc, theta)
    return [x_ellipse - x_line, y_ellipse - y_line]

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
    global ddt
    global displaygrid

    ddt = 0          # =1 for debugging, else =0
    displaygrid = 0  # =1 to display grid, else =0

    file_name = input("Input File Name: ")  # request file name from user
    if os.path.isfile(file_name) and os.access(file_name, os.R_OK):
      pass  # file exists and is readable
    else:  # File name supplied by user does not exist. 
      xstr = file_name.replace(" ","")  # delete all blank spaces
      if xstr == "":     # file name is null or blank
        # use default file name
        file_name = 'driver2d-coords.txt'
      else:  # print error message and quit.
        print("\n*** File not found: (" + file_name + ") ***\n")
        exit()    

    print("Processing file: " + file_name)
    errct = 0
    validct = 0
    valid = set("0123456789,.- ")  #list of valid characters for coordinates
    file1 = open(file_name, 'r', encoding='utf-8-sig') # remove BOM
    Lines = file1.readlines()
    for line in Lines:
      line = line.strip()    # remove leading and trailing spaces
      # print("line:", line)  ###
      if line[0:3] == "*2D" or line[0:3] == "*2d" :   # check for 2D process indicator
        sw2d = 1
      if line[0:3] == "*T=" or line[0:3] == "*t=" :   # check for title line
        title = line[3:]
      if line[0:3] == "*O=" or line[0:3] == "*o=" :   # check for options line
        if "DDT" in line or "ddt" in line :
          ddt = 1  # activate debugging 
        if "GRID" in line or "grid" in line :
          displaygrid = 1  # activate grid display in output image 
      # Skip blank lines and comments starting with * or #
      if not line.startswith('*') and not line.startswith('#') and not line == '':  
        lineset = set(line)
        if (lineset.issubset(valid) and "," in line):   # line contains only valid characters
          validct = validct + 1      
          line = line.split(",")
          # print("Lat: ", line[0], " Lon: ", line[1])  #print latitude and longitude
          if sw2d == 1:
            list1.append(line[1])   # plot (Y-coordinate)
            list2.append(line[0])   # plot (X-coordinate)
          else :
            list1.append(line[0])   # Latitude (Y-coordinate)
            list2.append(line[1])   # Longitude (X-coordinate)
        else:  # line contains invalid characters
          errct = errct + 1        
      # end if not * or # 
    file1.close()
    if errct > 0:
      print("Lines with errors=", errct)
    if validct < 5:
      print("Valid lines=",validct,"; five or more are required.")
      exit()
    print("")  # print blank line before main output

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
# Main program
# * * * * * * *
if __name__ == '__main__':
    
    X1, Y1 = read_coordinates() 

    X = np.array(list(zip(X1, Y1)))
    reg = LsqEllipse().fit(X)
    center, semimajor, semiminor, phi = reg.as_parameters()
    # phi = counterclockwise angle of rotation from the x-axis to the major-axis of the ellipse 
    # NOTE: LsqEllipse() does not reliably return samimajor > semiminor
    #       Any changes to LsqEllipse() would affect the display of the ellipse

    print(title)
    print(f'center: {center[0]:.3f}, {center[1]:.3f}')  # center of ellipse
    print(f'semimajor={semimajor:,.4f}, semiminor={semiminor:,.4f}, phi: {phi:.4f} ({np.rad2deg(phi):.3f} deg.)')    

    origphi = phi  # save original phi
    k1 = semimajor
    k2 = semiminor
    if k1 < k2:  # swap so major axis is larger than minor axis
      print("semimajor < semiminor")
      j = k2;
      k2 = k1
      k1 = j
      if phi > math.radians(90) : 
        phi -= math.radians(90)
      elif phi < 0 : 
        phi += math.radians(90)  
      else : 
        phi -= math.radians(90)  
    else : 
      print("semimajor > semiminor") 
      if phi > math.radians(90) : phi -= math.radians(180)

    if ddt > 0 : print(f' Adjusted phi = {phi:.4f} ({np.rad2deg(phi):.3f} deg.)')
    print("coefficients for  F(x,y) = ax**2 + 2bxy + cy**2 + 2dx + 2fy + g")
    a, b, c, d, f, g = reg.coeffs()
    print(f'a={a:.3f}, b={b:.3f}, c={c:.3f}, d={d:.3f}, f={f:.3f}, g={g:.3f}')

    X2 = np.array(X1)    # experimental points
    Y2 = np.array(Y1)
    xPred = 0; yPred = 0
    n1 = len(X2)
    calc_x = [0]*n1   # coordinates of points along the elliptical curve
    calc_y = [0]*n1
    print(f'Number of points = {n1}')

    n = 0
    MSE = 0   # Mean Squared Error [ 1/n * SUM(y_observed - y_predicted)^2 ]
    square_of_errors = 0   # (y_observed - y_predicted)^2 
    sumofsquares = 0 
    sumofabsoluteerrors = 0
    distance = 0
    sum_of_error_distances = 0
    while n < n1:
      x = X2[n]  # observed x
      y = Y2[n]  # observed y
      # Calculate theta (angle to the observed point) = arctan( (y-yc)/(x-xc) )
      theta = math.atan2( (y-center[1]), (x-center[0]) )
      if theta < 0 :
        theta = math.radians(360) + theta    # correct for quadrants III and IV
      if ddt > 0 :
        print(f'x[{n}]={x:,.4f} y[{n}]={y:,.4f} theta={theta:.4f} ({np.rad2deg(theta):.3f} deg.)')

      # calculate predicted intersection point on elliptical curve based on theta
      # Initial guess for t and r (important for convergence)
      # set initial guess: t=theta, r=semiminor
      initial_guess = [theta, k2]  
      # initial_guess = [0,1]   # override
      # Solve the system of equations
      t, r =  fsolve(equations, initial_guess, args=(center[0], center[1], k1*2, k2*2, phi, theta))
      # Calculate intersection point
      calc_x[n], calc_y[n] = ellipse_parametric(t, center[0], center[1], k1*2, k2*2, phi)
      if ddt > 0 :
        print(f'   calc_x[{n}]={calc_x[n]:.3f}, calc_y[{n}]={calc_y[n]:.3f}')

      square_of_errors = (y - calc_y[n])**2  # (y_observed - y_predicted)^2 
      sumofsquares += square_of_errors       # SUM(y_observed - y_predicted)^2
      sumofabsoluteerrors += abs(calc_y[n] - y)  # predicted minus observed value
      # distance between observed and predicted points
      distance = math.sqrt((x - calc_x[n])**2 + (square_of_errors))  
      sum_of_error_distances += distance

      n = n + 1

    print(f'center: {center[0]:.3f}, {center[1]:.3f}')
    print(f'Number of points = {n1}')
    print(f'Residual Sum of Squares (RSS) = {sumofsquares:,.4f}')
    MSE = sumofsquares/n
    print(f'Mean Squared Error (MSE) = {MSE:,.4f}')
    print(f'Mean Absolute Error (MAE) = {sumofabsoluteerrors/n:,.4f}')
    # print(f'Sum of error distances = {sum_of_error_distances:,.4f}')
    # the following varies according to the number of samples
    # print(f'Percent distance error relative to semiminor axis = {sum_of_error_distances*100/k2:.4f}%')
    average_error_distance = sum_of_error_distances/n
    print(f'Average error distance = {average_error_distance:.4f}')
    fitting_error = average_error_distance*100/k2
    # Percent average error distance relative to semiminor axis 
    print(f'Fitting error = {fitting_error:.4f}%')

    phi = origphi  # restore original phi
    if sw2d == 0 :
      # Calculate geographical coordinates for center
      Latitude = float(minlat) + center[1]/111111
      Longitude = float(maxlon) + center[0]/(111111*np.cos(np.radians(float(minlat))))
      print(f'center Lat.,Lon. {Latitude:.6f}, {Longitude:.6f}')

    print(f'major: {2*k1:.3f} m')
    print(f'minor: {2*k2:.3f} m')
    print(f'area: {np.pi*k1*k2:,.1f} square meters')    
    print(f'phi: {phi:.3f} ({np.rad2deg(phi):.3f} degrees)')
    # calculate azimuth
    if calc_azimuth == 1 :
      if phi < 0:       # phi is negative
        print("phi < 0")
        if semimajor < semiminor:     # case when major and minor axes were flipped.
          # print("semimajor < semiminor")
          azrad = np.pi + abs(phi)         # 180 + |phi| (phi<0)
        else:
          # print("semimajor > semiminor")  
          azrad = np.pi/2 + abs(phi)         # 90 + |phi| (phi<0)  #az 230307
      else:     # phi is positive
        print("phi > 0")
        if semimajor < semiminor:     # case when major and minor axes were flipped.
          # print("semimajor < semiminor")
          if phi > np.pi/2 :  # phi > 90
            azrad = 2*np.pi - phi      # 360 - phi
          else:
            azrad = np.pi - phi         # 180 - phi 
        else:
          # print("semimajor > semiminor")  # untested
          azrad = np.pi + np.pi/2 - phi
      # end if phi < 0: 
      # convert radians to degrees 
      print(f'azimuth: {np.degrees(azrad):.2f}')
    # end if calc_azimuth == 1  

    fig = plt.figure(figsize=(6, 6))
    ax = plt.subplot()
    ax.axis('equal')
    ax.plot(X1, Y1, 'ro', zorder=1)
    ax.plot(calc_x, calc_y, 'x')
    if ddt > 0 :
      ax.plot(calc_x[0], calc_y[0], 'go')  # green o for first calculated point
    ax.set_title(f'{title}\n major axis: {2*k1:.1f} m, minor axis: {2*k2:.1f} m, fitting error: {fitting_error:.3f}%')
    ax.plot(center[0], center[1], 'go', label='center')
    ellipse = Ellipse(
        xy=center, width=2*semimajor, height=2*semiminor, angle=np.rad2deg(phi),
        edgecolor='b', fc='None', lw=2, label='Fit', zorder=2
    )
    ax.add_patch(ellipse)

    plt.xlabel('$X_1$')
    plt.ylabel('$Y_1$')
    if displaygrid > 0:
      plt.grid()
    # plt.legend()
    plt.show()
