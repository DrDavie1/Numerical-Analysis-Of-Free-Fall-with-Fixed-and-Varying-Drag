#a) Write a Python program to plot the following analytical predictions for height 𝑦 and vertical speed 𝑣𝑦 as a function
#of time, for a free-falling object under constant gravity and constant drag factor, 𝑘:

import numpy as np
from math import *
import matplotlib.pyplot as plt

#________INIT PARAMETERS_________

g = 9.81 #gravitational constant
y0 = 1000 #inital height
rho_0 = 1.2 #air density
Cd = 1.2 #drag co-efficient
A = 1.8*0.5 #cross sectional area
m = 80 #mass
k = (Cd * rho_0 * A) / 2

#speed of sound:
M = 0.0289645 #molar mass dry air
gamma = 1.4
R = 8.31



#________THEORHETICAL EQUAIONS OF MOTION (A)_________

def y(y0,g,m,k,t):
    #Equation 4 (in report)
    yvals = np.array([])
    for time in t:
        y = y0 - (m/k)*log((cosh(sqrt((k*g)/m) * time)))
        yvals = np.append(yvals,y)
    return yvals


def vy(g,k,m,t):
    #Equation 5 (in report)
    vyvals = np.array([])
    for time in t:
        vy = - sqrt((m*g)/k) * tanh((sqrt((m*g)/k))*time)
        vyvals = np.append(vyvals,vy)
    return vyvals

def time_max(y0,g,m,k):
    #re-arrangement of equation 4 to find a max value of t

    t_max = sqrt(m/(k*g)) * acosh(exp(k*y0/m))

    return t_max

#__________DRAG FACTOR________

def drag_factor(y,Cd,A,rho_0,h):
    #Equation 2 (report)
    ky = Cd/2 * A * rho_0 * exp(-y/h)
    return ky


#__________TEMPERATURE AND SPEED OF SOUND________

def T_air(height):
    #how temperature changes with altitude
    if height <= 11000:
        T = 288 - 0.0065*height

    elif height > 11000 and height <= 25100:
        T = 216.5

    else:
        T = 141.3 + 0.0030*height

    return T


def v_sound(gamma,M,R,height):
    #how speed of sond changes with altitude(Equation 12)
    v = sqrt((gamma*R*T_air(height))/M)

    return v


#__________EULER METHOD EQUATIONS OF MOTION________

#CONSTANT DRAG FACTOR

def euler_free_fall_const(y0,g,m,k,delta_t):
    #Applying the Euler method for constant drag
    tvals = np.array([0])
    yvals = np.array([y0])
    vyvals = np.array([0])
    tn = 0
    yn = y0
    vyn = 0
    while yvals[len(yvals)-1] > 0: #while the most recent value of y is not 0
        tn += delta_t
        tvals = np.append(tvals,tn)
        yn += delta_t*vyn
        yvals = np.append(yvals,yn)
        vyn -= delta_t*( g + (k/m) * (abs(vyn) * vyn) )
        vyvals = np.append(vyvals,vyn)
    
    return tvals,yvals,vyvals

#CHANGING DRAG FACTOR

def euler_free_fall_dyn(y0,g,m,delta_t,Cd,rho_0,A,h,):
    #Applying the Euler method for varying drag.
    tvals = np.array([0])
    yvals = np.array([y0])
    vyvals = np.array([0])
    tn = 0
    yn = y0
    vyn = 0
    while yvals[len(yvals)-1] > 0: #while the most recent value of y is not 0
        tn += delta_t
        tvals = np.append(tvals,tn)
        yn += delta_t*vyn
        yvals = np.append(yvals,yn)

        k = drag_factor(yn,Cd,A,rho_0,h)

        vyn -= delta_t*( g + (k/m) * (abs(vyn) * vyn) )
        vyvals = np.append(vyvals,vyn)
    
    return tvals,yvals,vyvals

#________PLOTTING FUNCTION (A) and (B)_________

def plot_dis_speed(time_array,distance_array,velocity_array,title: str):
    #A function used to plot a Distance Time, Velocity Time figure 
    fig, (y_plot,vy_plot) = plt.subplots(1,2,figsize = (12,4))

    fig.suptitle(title)

    y_plot.set(xlabel='Time (s)', ylabel = 'Height (m)', title = 'Altitude')
    vy_plot.set(xlabel='Time (s)', ylabel = 'Speed (m/s)', title = 'Velocity')

    y_plot.plot(time_array,distance_array,color='red')
    vy_plot.plot(time_array,velocity_array,color='purple')

    plt.show()

#_________________(A)__________________

def section_A(NumPoints):
    #plotting the theorhetical model of free fall with constant drag
    tvals = np.linspace(0,time_max(y0,g,m,k),int(NumPoints))
    yvals = y(y0,g,m,k,tvals)
    vyvals = vy(g,m,k,tvals)
    
    plot_dis_speed(tvals,yvals,vyvals,'Height and speed data for a free-falling object under constant gravity and constant drag factor')

#_________________(B)__________________

def section_B(delta_t):
    #plotting the Euler model of free fall for constant drag
    tvals,yvals,vyvals = euler_free_fall_const(y0,g,m,k,delta_t)
    plot_dis_speed(tvals,yvals,vyvals,'Using Eulers numerical method to calculate height and speed data for a free-falling object ')


#_________________(C)__________________

def section_C(A,m,delta_t):
    #plooting the Euler model of free fall for varying drag 
    tvals,yvals,vyvals = euler_free_fall_dyn(39045,g,m,delta_t,Cd,rho_0,A,7640)
    plot_dis_speed(tvals,yvals,vyvals,'Using Eulers numerical method to calculate height and speed data for a free-falling object')

#_________________(D)__________________

def section_D(A,m,delta_t):
    #Finding Mach Number.
    tvals,yvals,vyvals = euler_free_fall_dyn(39045,g,m,delta_t,Cd,rho_0,A,7640)
    max_v = min(vyvals) #negative velocity
    index_max = np.where(vyvals == max_v)
    height_reached = yvals[index_max]

    sound_speed = v_sound(gamma,M,R,height_reached)

    mach_number = abs(max_v) / sound_speed

    fig, plot= plt.subplots(figsize = (6,4))

    fig.suptitle('The Sound Barrier')

    print('Mach Number:',mach_number)
    if mach_number > 1:
        print('Sound barrier would be broken')
        max = min(vyvals) - 20 
    else:
        print('Sound barrier would not be broken')
        max = -sound_speed - 20


    plot.plot(tvals,vyvals,color='purple')
    plot.set(xlabel='Time (s)', ylabel = 'Speed (m/s)', title = 'Velocity',ylim=(max,0))

    S_plot = plot.twinx()
    S_plot.set(ylim=(-max,0))
    S_plot.plot([tvals[0],tvals[int(len(tvals)-1)]],[sound_speed,sound_speed],color='black')



    plt.show()

#__________ADDITIONAL FUNCTIONS________


def changing_delta_T(NumPoints):

    #Comparing plots for a range of values of delta t.

    tvals = np.linspace(0,time_max(y0,g,m,k),NumPoints)
    yvals = y(y0,g,m,k,tvals)
    vyvals = vy(g,m,k,tvals)


    fig, (y_plot,vy_plot) = plt.subplots(1,2,figsize = (16,4))
    fig.suptitle('Comparing how changing delta effects numerical anaylsis')

   

    for i in range(0,5):
        temp_delta_t = 0.01 * 5**i 
        tvals_euler,yvals_euler,vyvals_euler = euler_free_fall_const(y0,g,m,k,temp_delta_t)
        current_yplot = y_plot.twinx()
        current_vyplot = vy_plot.twinx()
        current_yplot.set_axis_off()
        current_vyplot.set_axis_off()
        if i == 4:
            color = 'mistyrose'

        else:
            color = 'red'
        
        current_yplot.plot(tvals_euler,yvals_euler,color=color)
        current_vyplot.plot(tvals_euler,vyvals_euler,color=color)

    A_yplot = y_plot.twinx()
    A_vyplot = vy_plot.twinx()
    A_yplot.set_axis_off()
    A_vyplot.set_axis_off()

    A_yplot.plot(tvals,yvals,color='black')
    A_vyplot.plot(tvals,vyvals,color='black')


    y_plot.set(xlabel='Time (s)', ylabel = 'Height (m)', title = 'Altitude' , ylim= (0,max(yvals)))
    vy_plot.set(xlabel='Time (s)', ylabel = 'Speed (m/s)', title = 'Velocity', ylim= (min(vyvals),0))
    plt.show()


def changing_sound():
    #showing how speed of sound changes with height

    heightvals = np.linspace(0,40_000,1000)
    sound_vvals = np.array([])
    for h in heightvals:
        v_s = v_sound(gamma,M,R,h)
        sound_vvals = np.append(sound_vvals,v_s)
    
    fig, vs_plot= plt.subplots(1,figsize = (6,4))

    fig.suptitle('How speed of sound varies with altitude')

    vs_plot.set(xlabel='Height (m)', ylabel = 'Speed of Sound (m/s)')

    vs_plot.plot(heightvals,sound_vvals,color='red')


    plt.show()



#_________________Running the code__________________

MyInput = '0'

while MyInput != 'q':
    MyInput = input('Enter a choice, "a", "b", "c", "d" or "q" to quit: ')

    if MyInput == 'a':
        print('You have chosen part (a)')
        print('g = ',g,'air density = ',rho_0,'drag co-efficient = ',Cd,"start height = ",y0)
        int_ = 0
        while int_ == 0:
            NumPoints = input('Input number of points: ')
            try:
                NumPoints = abs(int(NumPoints))
                int_ += 1 
            except:
                print('input an interger')

        section_A(NumPoints)


    elif MyInput == 'b':
        print('You have chosen part (b)')
        print("================================================")
        print('this part users the Euler Method so you will need to chose a value for delta t')
        print("================================================")
        print('g = ',g,'air density = ',rho_0,'drag co-efficient = ',Cd,"start height = ",y0)
        float_ = 0
        choice = input('Would you like to chose a single value of delta t (input: "A") or show how delta t affects accuracy (input: "B"): ')
        if choice == 'A':
            while float_ == 0:
                delta_t = input('Input delta t: ')
                try:
                    delta_t = float(delta_t)
                    section_B(delta_t)
                    float_ += 1
                except:
                    print('input a floating point number')
        else:
            print('delta t = 0.01,0.05,0.25,1.25 and 6.25, black line is produced by theorhetical model')
            changing_delta_T(200)

    elif MyInput == 'c':
        print('You have chosen part (c)')
        print("This is a model of Felix's jump")
        float_ = 0
        while float_ == 0:
            m = input("Input a mass: ")
            A = input("Input cross sectional area: ")
            delta_t = input('Input delta t: ')

            try:
                delta_t = float(delta_t)
                m = float(m)
                A = float(A)
                section_C(A,m,delta_t)
                float_ += 1
            except:
                print('input a floating point number')

    elif MyInput == 'd':
        print('You have chosen part (c)')
        choice = input('Would you like to see how sound varies with altitude? ("yes" or "no"): ')
        if choice == 'yes':
            changing_sound()

        else:
            None
        print("================================================")
        print('Calculating Mach Number for a specific set of variables.')
        print("================================================")
        float_ = 0
        while float_ == 0:
            m = input("Input a mass: ")
            A = input("Input cross sectional area: ")
            delta_t = input('Input delta t: ')

            try:
                delta_t = float(delta_t)
                m = float(m)
                A = float(A)
                section_D(A,m,delta_t)
                float_ += 1

            except:
                print('input a floating point number')

    elif MyInput != 'q':
        print('This is not a valid choice')


print('You have chosen to finish - goodbye.')


            
        










    