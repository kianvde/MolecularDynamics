__author__ = 'Kian'
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation


class Animate(object):
    def __init__(self, boxSize, sizeSimulation, dims, numParticles):
        self.boxSize = boxSize
        self.lineData = np.zeros((numParticles, dims, sizeSimulation))
        # fig = plt.figure()
        # ax = p3.Axes3D(fig)
        # ax.set_xlim3d([0.0, self.boxSize])
        # ax.set_xlabel('X')
        #
        # ax.set_ylim3d([0.0, self.boxSize])
        # ax.set_ylabel('Y')
        #
        # ax.set_zlim3d([0.0, self.boxSize])
        # ax.set_zlabel('Z')
        #
        # ax.set_title('3D Test')
        #
        # self.parts, = ax.plot([], [], [], 'bo', ms=6)
        # self.parts.set_data([], [], [])

    def buildLines(self, loopcount, positions):
        self.lineData[:, :, loopcount] = positions
        print np.shape(self.lineData)

    def buildCoords(self, loopcount, positions):
        self.lineData[:, :, loopcount] = positions
        self.x = self.lineData[:, 0, :]
        self.y = self.lineData[:, 1, :]
        self.z = self.lineData[:, 2, :]

    def buildCoords2(self, loopcount, positions):
        self.parts.set_data[:, :, loopcount] = positions
        self.x = self.lineData[:, 0, :]
        self.y = self.lineData[:, 1, :]
        self.z = self.lineData[:, 2, :]

    def plotOne(self):
        fig2 = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter(self.x[0, :], self.y[0, :], self.z[0, :])
        plt.show()


    # def plot_anim2(self):
    #     fig3 = plt.figure()
    #     ax = plt.axes(projection='3d')
    #     scat = plt.scatter(self.x, self.y, self.z)



    def update_lines(self, num, dataLines, lines):
        for line, data in zip(lines, dataLines):
            # NOTE: there is no .set_data() for 3 dim data...
            line.set_data(data[0:2, :num])
            line.set_3d_properties(data[2, :num])
        return lines

    def plot_anim(self):
        # Attaching 3D axis to the figure
        fig = plt.figure()
        ax = p3.Axes3D(fig)
        lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in self.lineData]

        ax.set_xlim3d([0.0, self.boxSize])
        ax.set_xlabel('X')

        ax.set_ylim3d([0.0, self.boxSize])
        ax.set_ylabel('Y')

        ax.set_zlim3d([0.0, self.boxSize])
        ax.set_zlabel('Z')

        ax.set_title('3D Test')

        line_ani = animation.FuncAnimation(fig, self.update_lines, 25, fargs=(self.lineData, lines),
                                           interval=1000, blit=False)

        plt.show()