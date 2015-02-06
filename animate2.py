__author__ = 'Kian'
# import numpy as np
# import matplotlib.pyplot as plt
# import mpl_toolkits.mplot3d.axes3d as p3
# import matplotlib.animation as animation
#
# class Animate(object):
#
#     def __init__(self, boxSize, sizeSimulation, dims, numParticles):
#         self.boxSize = boxSize
#         self.lineData = np.zeros((numParticles, dims, sizeSimulation))
#         fig = plt.figure()
#         ax = p3.Axes3D(fig)
#         ax.set_xlim3d([0.0, self.boxSize])
#         ax.set_xlabel('X')
#
#         ax.set_ylim3d([0.0, self.boxSize])
#         ax.set_ylabel('Y')
#
#         ax.set_zlim3d([0.0, self.boxSize])
#         ax.set_zlabel('Z')
#
#         ax.set_title('3D Test')
#         self.ax = ax
#         #
#         # self.parts, = ax.plot([], [], [], 'bo', ms=6)
#         # print type(self.parts)
#         # self.parts.set_data([], [], [])
#     def buildCoords(self, loopcount, positions):
#         self.lineData[:, :, loopcount] = positions
#         self.x = self.lineData[:, 0, :]
#         self.y = self.lineData[:, 1, :]
#         self.z = self.lineData[:, 2, :]
#
#     def anim(self):
#         self.parts, = self.ax.plot3D([], [], [])
#         print self.parts
# anim = Animate(100, 15, 3, 8)
# anim.anim()

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

class AnimatedScatter(object):
    def __init__(self, boxSize, positions,numParticles, dims, sizeSimulation):
        self.positions = positions
        self.lineData = np.zeros((numParticles, dims, sizeSimulation))
        print np.shape(self.positions)
        self.boxSize = boxSize
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111,projection = '3d')
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=1000,
                                           init_func=self.setup_plot, blit=True)
    def buildCoords(self, loopcount, positions):
        self.lineData[:, :, loopcount] = positions

    def setup_plot(self):
        # print np.shape(self.positions)
        X = self.positions
        self.scat = self.ax.scatter(X[:, 0], X[:, 1], X[:, 2], s=100, animated=True)

        self.ax.set_xlim3d(0, self.boxSize)
        self.ax.set_ylim3d(0, self.boxSize)
        self.ax.set_zlim3d(0, self.boxSize)
        return self.scat,

    def update(self, i):
        data = self.lineData[:, :, i]
        data = np.transpose(data)
        self.scat._offsets3d = ( np.ma.ravel(data[:,0]) , np.ma.ravel(data[:,0]) , np.ma.ravel(data[:,0]) )
        plt.draw()
        return self.scat,

    def show(self):
        plt.show()
#
# if __name__ == '__main__':
#     a = AnimatedScatter()
#     a.show()