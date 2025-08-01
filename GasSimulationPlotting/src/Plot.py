import numpy as np
import string as str
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.size': 14, 'font.family': 'sans-serif'})
#matplotlib.rcParams.update({'text.usetex': 'true'})

def makePlot():
    path = r'C:\Users\OrcaF\OneDrive\Documents\Github\Gas-Simulation'
    fid = open(path + '\Parameters.txt', 'r')

    linestring = fid.readline()
    linelist = linestring.split()
    nX = np.int64(linelist[1])

    linestring = fid.readline()
    linelist = linestring.split()
    nV = np.int64(linelist[1])
    
    linestring = fid.readline()
    linelist = linestring.split()
    aX = np.float64(linelist[1])
    
    linestring = fid.readline()
    linelist = linestring.split()
    bX = np.float64(linelist[1])
    
    linestring = fid.readline()
    linelist = linestring.split()
    aV = np.float64(linelist[1])
    
    linestring = fid.readline()
    linelist = linestring.split()
    bV = np.float64(linelist[1])
    
    linestring = fid.readline()
    linelist = linestring.split()
    knudson = np.float64(linelist[2])

    linestring = fid.readline()
    linelist = linestring.split()
    finalTime = np.float64(linelist[2])
    ftstr = '%4.2f'%finalTime

    fid.close()

    x = np.zeros(nX)
    v = np.zeros(nV)
    dX = (bX-aX)/nX
    dV = (bV-aV)/nV
    for i in range(0, nX):
        x[i] = aX + (i + 0.5)*dX
    for j in range(0, nV):
        v[j] = aV + (j + 0.5)*dV

    xl = np.zeros((nX+1))
    vl = np.zeros((nV+1))

    for i in range(0, nX+1):
        xl[i] = aX + i*dX

    for j in range(0, nV+1):
        vl[j] = aV + j*dV

    f0 = np.zeros((nX, nV))
    fid = open(path + '\OutputF0.txt', 'r')

    for i in range(0, nX):
        for j in range(0, nV):
            linestring = fid.readline()
            linelist = linestring.split()
            f0[i, j] = np.float64(linelist[0])
    fid.close()

    rho0 = np.zeros(nX)
    fid = open(path + '\OutputRho0.txt', 'r')

    for i in range(0, nX):
        linestring = fid.readline()
        linelist = linestring.split()
        rho0[i] = np.float64(linelist[0])
    fid.close()

    u0 = np.zeros(nX)
    fid = open(path + '\OutputU0.txt', 'r')

    for i in range(0, nX):
        linestring = fid.readline()
        linelist = linestring.split()
        u0[i] = np.float64(linelist[0])
    fid.close()

    T0 = np.zeros(nX)
    fid = open(path + '\OutputT0.txt', 'r')

    for i in range(0, nX):
        linestring = fid.readline()
        linelist = linestring.split()
        T0[i] = np.float64(linelist[0])
    fid.close()

    


    f1 = np.zeros((nX, nV))
    fid = open(path + '\OutputF1.txt', 'r')

    for i in range(0, nX):
        for j in range(0, nV):
            linestring = fid.readline()
            linelist = linestring.split()
            f1[i, j] = np.float64(linelist[0])
    fid.close()

    rho1 = np.zeros(nX)
    fid = open(path + '\OutputRho1.txt', 'r')

    for i in range(0, nX):
        linestring = fid.readline()
        linelist = linestring.split()
        rho1[i] = np.float64(linelist[0])
    fid.close()

    u1 = np.zeros(nX)
    fid = open(path + '\OutputU1.txt', 'r')

    for i in range(0, nX):
        linestring = fid.readline()
        linelist = linestring.split()
        u1[i] = np.float64(linelist[0])
    fid.close()

    T1 = np.zeros(nX)
    fid = open(path + '\OutputT1.txt', 'r')

    for i in range(0, nX):
        linestring = fid.readline()
        linelist = linestring.split()
        T1[i] = np.float64(linelist[0])
    fid.close()

    plt.figure(1)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([aX, bX])
    plt.plot(x, rho0, 'r-', linewidth=2.0)
    plt.grid(True)
    plt.title('Density at Time t=0')
    plt.xlabel('x')
    plt.ylabel('Density')
    plt.savefig('rho0.pdf', format='pdf', bbox_inches='tight')

    plt.figure(2)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([aX, bX])
    plt.plot(x, u0, 'r-', linewidth=2.0)
    plt.grid(True)
    plt.title('Velocity at Time t=0')
    plt.xlabel('x')
    plt.ylabel('Velocity')
    plt.savefig('u0.pdf', format='pdf', bbox_inches='tight')

    plt.figure(3)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([aX, bX])
    plt.plot(x, T0, 'r-', linewidth=2.0)
    plt.grid(True)
    plt.title('Temperature at Time t=0')
    plt.xlabel('x')
    plt.ylabel('Temperature')
    plt.savefig('T0.pdf', format='pdf', bbox_inches='tight')

    plt.figure(4)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([aX, bX])
    plt.gca().set_ylim([aV, bV])
    plt.pcolormesh(x, v, np.transpose(f0), linewidth = 0)
    plt.title('PDF at Time t=' + ftstr)
    plt.xlabel('x')
    plt.ylabel('v')
    plt.colorbar()
    plt.savefig('f0.png', format='png', bbox_inches='tight')

    ##########

    plt.figure(5)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([aX, bX])
    plt.plot(x, rho1, 'r-', linewidth=2.0)
    plt.grid(True)
    plt.title('Density at Time t=' + ftstr)
    plt.xlabel('x')
    plt.ylabel('Density')
    plt.savefig('rho1.pdf', format='pdf', bbox_inches='tight')

    plt.figure(6)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([aX, bX])
    plt.plot(x, u1, 'r-', linewidth=2.0)
    plt.grid(True)
    plt.title('Velocity at Time t=' + ftstr)
    plt.xlabel('x')
    plt.ylabel('Velocity')
    plt.savefig('u1.pdf', format='pdf', bbox_inches='tight')

    plt.figure(7)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([aX, bX])
    plt.plot(x, T1, 'r-', linewidth=2.0)
    plt.grid(True)
    plt.title('Temperature at Time t=' + ftstr)
    plt.xlabel('x')
    plt.ylabel('Temperature')
    plt.savefig('T1.pdf', format='pdf', bbox_inches='tight')

    print(x.shape)
    print(v.shape)
    print(f1.shape)
    plt.figure(8)
    plt.clf()
    plt.gca().set_aspect('auto')
    plt.gca().set_xlim([aX, bX])
    plt.gca().set_ylim([aV, bV])
    plt.pcolormesh(xl, vl, np.transpose(f1), linewidth = 0)
    plt.title('PDF at Time t=' + ftstr)
    plt.xlabel('x')
    plt.ylabel('v')
    plt.colorbar()
    plt.savefig('f1.png', format='png', bbox_inches='tight')

if __name__ == "__main__":
    makePlot()