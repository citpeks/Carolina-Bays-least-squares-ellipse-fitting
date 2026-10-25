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
# 02/09/2025 - Used non-iterative method to calculate intersection point of line and ellipse
# 02/11/2025 - Allowed option to use EllipseModel() in skimage  *O= model
# 12/29/2025 - Calculated Scaled Mean Squared Error (SMSE) and MAD/semiminor
# 06/04/2026 - Implemented World Geodetic System WGS-84 ellipsoidal formula
# 10/25/2026 - Improved diagnostics for bad data

import numpy as np
import math
import re
# from skimage.measure import EllipseModel   # remove pound sign to activate this feature
from ellipse import LsqEllipse
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.optimize import fsolve
import os
import os.path

# These are global variables
minlat = ""  
maxlon = "" 
ddt = 0    # debugging switch
sw2d = 0   # 2D process indicator
calc_azimuth = 0  # =1 to calculate azimuth
title = ""

# * * * * * * *
# Calculate the coordinates of the intersection point of the line and ellipse
# * * * * * * *
def ellipse_line_intersection(Xc, Yc, L, W, phi, theta):
    # Calculate the parametric equations of the line
    dx = (L / 2) * np.cos(theta)
    dy = (L / 2) * np.sin(theta)
    
    # Transform the ellipse to the standard position
    a = L / 2
    b = W / 2
    
    # Rotate the line to align the ellipse with the coordinate axes
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    x1_rot = (dx * cos_phi + dy * sin_phi)
    y1_rot = (-dx * sin_phi + dy * cos_phi)
    
    # Coefficients of the quadratic equation
    A = (x1_rot / a) ** 2 + (y1_rot / b) ** 2
    B = 2 * ((Xc * cos_phi + Yc * sin_phi - Xc * cos_phi - Yc * sin_phi) * x1_rot / a ** 2 + 
             (Xc * sin_phi - Yc * cos_phi - Xc * sin_phi + Yc * cos_phi) * y1_rot / b ** 2)
    C = ((Xc * cos_phi + Yc * sin_phi - Xc * cos_phi - Yc * sin_phi) ** 2 / a ** 2 +
         (Xc * sin_phi - Yc * cos_phi - Xc * sin_phi + Yc * cos_phi) ** 2 / b ** 2) - 1
    
    # Solve the quadratic equation
    discriminant = B ** 2 - 4 * A * C
    if discriminant < 0:
        return []  # No intersection
    
    t1 = (-B + np.sqrt(discriminant)) / (2 * A)
    t2 = (-B - np.sqrt(discriminant)) / (2 * A)
    
    # Calculate the intersection points
    x1 = Xc + dx * t1  # positive root
    y1 = Yc + dy * t1
    x2 = Xc + dx * t2  # negative root
    y2 = Yc + dy * t2
    
    return [(x1, y1), (x2, y2)]  # return coordinates for positive and negative solutions

# * * * * * * *
# --- Ellipsoidal formulas for meters per degree (WGS-84) ---
# Latitude-dependent length of one degree of latitude:
#   M = 111132.92 - 559.82*cos(2*phi) + 1.175*cos(4*phi)
#
# Longitude-dependent length of one degree of longitude:
#   N = 111412.84*cos(phi) - 93.5*cos(3*phi)
# * * * * * * *

def meters_per_degree_lat(phi_deg):
    phi = np.radians(phi_deg)
    return (111132.92
            - 559.82 * np.cos(2 * phi)
            + 1.175 * np.cos(4 * phi))

def meters_per_degree_lon(phi_deg):
    phi = np.radians(phi_deg)
    return (111412.84 * np.cos(phi)
            - 93.5 * np.cos(3 * phi))

# * * * * * * *
# Read file with coordinates and return two lists containing the x and y data of the ellipse
# * * * * * * *
def read_coordinates():
    list1= []
    list2 = []
    global minlat
    global maxlon
    global sw2d
    global calc_azimuth
    global title
    global displaygrid
    global model

    ddt = 0          # =1 for debugging, else =0
    displaygrid = 0  # =1 to display grid, else =0
    model = 0        # = 1 to use EllipseModel() in skimage instead of lsqEllipse

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
    coords = 1  # assume imput is geographical coordinates
    valid = set("0123456789,.- ")  #list of valid characters for coordinates
    file1 = open(file_name, 'r', encoding='utf-8-sig') # remove BOM
    Lines = file1.readlines()
    for line in Lines:
      line = line.strip()    # remove leading and trailing spaces
      line = re.sub(r'\t', ' ', line) # change tabs to spaces
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
        if "MODEL" in line or "model" in line :
          model = 1  # activate EllipseModel() in skimage 
      # Skip blank lines and comments starting with * or #
      if not line.startswith('*') and not line.startswith('#') and not line == '':  
        lineset = set(line)
        if (lineset.issubset(valid) and "," in line):   # line contains only valid characters
          validct = validct + 1      
          line = line.split(",")
          # print("Lat: ", line[0], " Lon: ", line[1])  #print latitude and longitude
          # maximum value for latitude is 90 and for longitude 180
          if int(abs(float(line[0]))) > 90 or int(abs(float(line[1]))) > 180 :
            coords = 0  # not geographical coordinates
          if sw2d == 1:
            list1.append(line[1])   # plot (Y-coordinate)
            list2.append(line[0])   # plot (X-coordinate)
          else :
            list1.append(line[0])   # Latitude (Y-coordinate)
            list2.append(line[1])   # Longitude (X-coordinate)
        else:  # line contains invalid characters
          errct = errct + 1  
          print("***ERROR*** ", line)       
      # end if not * or # 
    file1.close()
    if errct > 0:
      print("Lines with errors=", errct)
      if sw2d == 1 and coords == 1 :
        print("***WARNING*** Invalid use of *2D option")
      if sw2d == 0 and coords == 0 :
        print("***ERROR*** Use *2D option")
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
    # Compute meters-per-degree at the reference latitude
    lat0 = float(minlat)
    M_lat = meters_per_degree_lat(lat0)   # meters per degree latitude
    N_lon = meters_per_degree_lon(lat0)   # meters per degree longitude

    # Convert latitude and longitude to meters relative to minimum coordinates
    # This will place the ellipse in the first quadrant

    # Convert latitude to meters
    ellipse_y = []
    for j in list1:
      if sw2d == 0:
        ellipse_y.append((float(j) - float(minlat)) * M_lat)
      else:
        ellipse_y.append(float(j) - float(minlat))

    # Convert longitude to meters
    ellipse_x = []
    for j in list2:
      if sw2d == 0:
        dx_deg = abs(float(maxlon)) - abs(float(j))
        ellipse_x.append(dx_deg * N_lon)
      else:
        ellipse_x.append(float(j) - float(minlon))

    return [ellipse_x, ellipse_y]


# * * * * * * *
# Main program
# * * * * * * *
if __name__ == '__main__':
    
    X1, Y1 = read_coordinates() 
    X = np.array(list(zip(X1, Y1)))
    # fit an ellipse to the set of points using a least squares method
    if model > 0 :  # use EllipseModel() in skimage 
      print(f'Using EllipseModel() in skimage.')
      reg = EllipseModel()  
      reg.estimate(X)
      xc, yc, semimajor, semiminor, phi = reg.params  
      center = (xc, yc)
    else :
      reg = LsqEllipse().fit(X)
      center, semimajor, semiminor, phi = reg.as_parameters()
      # phi = counterclockwise angle of rotation from the x-axis to the major-axis of the ellipse 
      print("coefficients for  F(x,y) = ax**2 + 2bxy + cy**2 + 2dx + 2fy + g")
      a, b, c, d, f, g = reg.coeffs()
      print(f'a={a:.3f}, b={b:.3f}, c={c:.3f}, d={d:.3f}, f={f:.3f}, g={g:.3f}')

    print(title)
    print(f'center: {center[0]:.3f}, {center[1]:.3f}')  # center of ellipse
    print(f'semimajor={semimajor:,.4f}, semiminor={semiminor:,.4f}, phi: {phi:.4f} ({np.rad2deg(phi):.3f} deg.)')    

    origphi = phi  # save original phi
    k1 = semimajor
    k2 = semiminor
    if k1 < k2:  # swap so major axis is larger than minor axis
      # NOTE: The previous version of LsqEllipse() did not guarantee semimajor > semiminor
      #       This program tests and corrects for this condition here and in the azimuth calculation
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
      if ddt > 0 : print("semimajor > semiminor") 
      if phi > math.radians(90) : phi -= math.radians(180)

    if ddt > 0 : print(f' Adjusted phi = {phi:.4f} ({np.rad2deg(phi):.3f} deg.)')

    X2 = np.array(X1)    # experimental points
    Y2 = np.array(Y1)
    xPred = 0; yPred = 0
    n1 = len(X2)
    calc_x = [0]*n1   # coordinates of points along the elliptical curve
    calc_y = [0]*n1
    scaledydiff = [0]*n1
    #  print(f'Number of points = {n1}')

    # Calculate scaling_factor relative to a semiminor axis of 1000
    scaling_factor = 1000/semiminor
    print(f'scaling factor = {scaling_factor}')

    n = 0
    MSE = 0   # Mean Squared Error [ 1/n * SUM(y_observed - y_predicted)^2 ]
    SMSE = 0  # Scaled Mean Squared Error
    square_of_errors = 0   # (y_observed - y_predicted)^2 
    scaled_square_of_errors = 0
    sumofsquares = 0 
    sumofscaledsquares = 0
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
      intersection_points = ellipse_line_intersection(center[0], center[1], k1*2, k2*2, phi, theta)
      calc_x[n], calc_y[n] = intersection_points[0]  # select only the positive solution
      if ddt > 0 :
        print(f'   calc_x[{n}]={calc_x[n]:.3f}, calc_y[{n}]={calc_y[n]:.3f}')

      square_of_errors = (y - calc_y[n])**2  # (y_observed - y_predicted)^2 
      # The scaling is applied to the linear distance before squaring
      scaledydiff[n] = scaling_factor * (y - calc_y[n])  # this is a signed number
      scaled_square_of_errors = (scaling_factor * abs(y - calc_y[n]))**2
      sumofsquares += square_of_errors       # SUM(y_observed - y_predicted)^2
      sumofscaledsquares += scaled_square_of_errors
      sumofabsoluteerrors += abs(calc_y[n] - y)  # predicted minus observed value
      # distance between observed and predicted points
      distance = math.sqrt((x - calc_x[n])**2 + (square_of_errors))  
      sum_of_error_distances += distance

      n = n + 1

    # print(f'center: {center[0]:.3f}, {center[1]:.3f}')
    print(f'Number of points = {n1}')
    print(f'Residual Sum of Squares (RSS) = {sumofsquares:,.4f}')
    MSE = sumofsquares/n
    print(f'Mean Squared Error (MSE) = {MSE:,.4f}')
    sigma = math.sqrt(MSE)  # standard deviation
    print(f'Standard Deviation = {sigma:,.4f}')
    print(f'Residual Sum of Scaled Squares (RSSS) = {sumofscaledsquares:,.4f}') 
    SMSE = sumofscaledsquares/n
    print(f'Scaled Mean Squared Error (SMSE) = {SMSE:,.4f}')
    print(f'Scaled Deviation = {math.sqrt(SMSE):,.4f}')
    MAD = sumofabsoluteerrors/n1
    print(f'Mean Absolute Deviation (MAD) = {MAD:,.4f}')
    print(f'MAD/semiminor = {(sumofabsoluteerrors*100)/(n*k2):,.4f}%')
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
      # Compute meters-per-degree at the reference latitude
      lat0 = float(minlat)
      M_lat = meters_per_degree_lat(lat0)   # meters per degree latitude
      N_lon = meters_per_degree_lon(lat0)   # meters per degree longitude
      Latitude = float(minlat) + center[1]/M_lat
      Longitude = float(maxlon) + center[0]/N_lon
      print(f'WGS center Lat.,Lon. {Latitude:.6f}, {Longitude:.6f}')


    print(f'major: {2*k1:.3f} m')
    print(f'minor: {2*k2:.3f} m')
    print(f'area: {np.pi*k1*k2:,.1f} square meters')    
    print(f'phi: {phi:.3f} ({np.rad2deg(phi):.3f} degrees)')
    # calculate azimuth
    if calc_azimuth == 1 :
      if phi < 0:       # phi is negative
        # print("phi < 0")
        if semimajor < semiminor:     # case when major and minor axes were flipped.
          # print("semimajor < semiminor")
          azrad = np.pi + abs(phi)         # 180 + |phi| (phi<0)
        else:
          # print("semimajor > semiminor")  
          azrad = np.pi/2 + abs(phi)         # 90 + |phi| (phi<0)  #az 230307
      else:     # phi is positive
        # print("phi > 0")
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

    # print(f'\n== Scaled difference ==')
    # for item in scaledydiff:
    #   print(item)

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
    if sw2d == 0 :
      ax.text(x=0.98,y=0.98,s='WGS',transform=ax.transAxes,    
      fontsize=10, va='top', ha='right' )  # World Geodetic System
    if displaygrid > 0:
      plt.grid()
    # plt.legend()
    plt.show()
