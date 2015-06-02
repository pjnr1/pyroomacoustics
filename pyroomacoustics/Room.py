# @version: 1.0  date: 05/06/2015 by Sidney Barthe
# @author: robin.scheibler@epfl.ch, ivan.dokmanic@epfl.ch, sidney.barthe@epfl.ch
# @copyright: EPFL-IC-LCAV 2015

import numpy as np

import beamforming as bf
from soundsource import SoundSource
from wall import Wall
from utilities import area, ccw3p
import constants



class Room(object):
    """
    This class represents a room instance.
    
    A room instance is formed by wall instances. Microphones and sound sources can be added.
    
    :attribute walls: (Wall array) list of walls forming the room
    :attribute fs: (int) sampling frequency
    :attribute t0: (float) time offset
    :attribute max_order: (int) the maximum computed order for images sources
    :attribute sigma2_awgn: (?) ambient additive white gaussian noise level
    :attribute sources: (SoundSource array) list of sound sources
    :attribute mics: (?) array of microphones
    :attribute normals: (np.array 2xN or 3xN, N=number of walls) array containing normal vector for each wall, used for calculations
    :attribute corners: (np.array 2xN or 3xN, N=number of walls) array containing a point belonging to each wall, used for calculations
    :attribute absorption: (np.array size N, N=number of walls)  array containing the absorption factor for each wall, used for calculations
    :attribute dim: (int) dimension of the room (2 or 3 meaning 2D or 3D)
    :attribute wallsId: (int dictionary) stores the mapping "wall name -> wall id (in the array walls)"
    """


    def __init__(
            self,
            walls,
            fs=8000,
            t0=0.,
            max_order=1,
            sigma2_awgn=None,
            sources=None,
            mics=None):

        self.walls = walls
        self.fs = fs
        self.t0 = t0
        self.max_order = max_order
        self.sigma2_awgn = sigma2_awgn
        
        if (sources is list):
            self.sources = sources
        else:
            self.sources = []

        if (mics is list):
            self.micArray = mics
        else:
            self.micArray = []
            
        self.normals = np.array([wall.normal for wall in self.walls]).T
        self.corners = np.array([wall.corners[:, 0] for wall in self.walls]).T
        self.absorption = np.array([wall.absorption for wall in self.walls])

        # Pre-compute RIR if needed
        if (len(self.sources) > 0 and self.micArray is not None):
            self.compute_RIR()
        else:
            self.rir = []
            
        self.dim = walls[0].dim
        self.wallsId = {}
        for i in range(len(walls)):
            if self.walls[i].name is not None:
                self.wallsId[self.walls[i].name] = i

    @classmethod
    def shoeBox2D(
            cls,
            p1,
            p2,
            absorption=1.,
            fs=8000,
            t0=0.,
            max_order=1,
            sigma2_awgn=None,
            sources=None,
            mics=None):
        """
        Creates a 2D "shoe box" room geometry (rectangle).
        
        :arg p1: (np.array dim 2) coordinates of the lower left corner of the room
        :arg p2: (np.array dim 2) coordinates the upper right corner of the room
        :arg absorption: (float) absorption factor reflection for all walls
        
        :returns: (Room) instance of a 2D shoe box room
        """

        walls = []
        walls.append(Wall(np.array([[p1[0], p2[0]], [p1[1], p1[1]]]), absorption, "south"))
        walls.append(Wall(np.array([[p2[0], p2[0]], [p1[1], p2[1]]]), absorption, "east"))
        walls.append(Wall(np.array([[p2[0], p1[0]], [p2[1], p2[1]]]), absorption, "north"))
        walls.append(Wall(np.array([[p1[0], p1[0]], [p2[1], p1[1]]]), absorption, "west"))

        return cls(walls, fs, t0, max_order, sigma2_awgn, sources, mics)
        
    @classmethod
    def fromCorners(cls, corners, absorption=1., **kwargs):
        """
        Creates a 2D room by giving an array of corners.
        
        :arg corners: (np.array dim 2xN, N>2) list of corners, must be antiClockwise oriented
        :arg absorption: (float array or float) list of absorption factor reflection for each wall or single value for all walls
        
        :returns: (Room) instance of a 2D room
        """
        
        corners = np.array(corners)
        if (corners.size[0] != 2 or len(corners) < 3):
            raise NameError('Room.fromCorners input error : corners must be more than two 2D points.')

        if (area(corners) >= 0):
            raise NameError('Room.fromCorners input error : corners must be anti-clockwise ordered.')

        self.corners = corners
        self.dim = corners.shape[0] 
            
        absorption = np.array(absorption, dtype='float64')
        if (absorption.ndim == 0):
            self.absorption = absorption * np.ones(corners.shape[1])
        elif (absorption.ndim > 1 and corners.shape[1] != len(absorption)):
            raise NameError('Room.fromCorners input error : absorption must be the same size as corners or must be a single value.')
        else:
            self.absorption = absorption
        
        walls = []
        for i in range(len(corners)):
            walls.append(Wall(np.array([corners[:, i], corners[:, (i+1)%len(corners)]]), absorption[i], "wall_"+str(i)))
            
        return cls(walls, **kwargs)

    def plot(self, img_order=None, freq=None, figsize=None, no_axis=False, mic_marker_size=10, **kwargs):

        import matplotlib
        from matplotlib.patches import Circle, Wedge, Polygon
        from matplotlib.collections import PatchCollection
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=figsize)


        if no_axis is True:
            ax = fig.add_axes([0, 0, 1, 1], aspect='equal', **kwargs)
            ax.axis('off')
            rect = fig.patch
            rect.set_facecolor('gray')
            rect.set_alpha(0.15)

        else:
            ax = fig.add_subplot(111, aspect='equal', **kwargs)

        # draw room
        polygons = [Polygon(self.corners.T, True)]
        # if no_axis is True:
        #   r = Rectangle((xlim[0],ylim[0]), xlim[1]-xlim[0], ylim[1]-ylim[0])
        #   polygons.append(r)
        p = PatchCollection(polygons, cmap=matplotlib.cm.jet,
                facecolor=np.array([1, 1, 1]), edgecolor=np.array([0, 0, 0]))
        ax.add_collection(p)

        # draw the microphones
        if (self.micArray is not None):
            for mic in self.micArray.R.T:
                ax.scatter(mic[0], mic[1],
                        marker='x', linewidth=0.5, s=mic_marker_size, c='k')

            # draw the beam pattern of the beamformer if requested (and available)
            if freq is not None \
                    and isinstance(self.micArray, bf.Beamformer) \
                    and (self.micArray.weights is not None
                            or self.micArray.filters is not None):

                freq = np.array(freq)
                if freq.ndim is 0:
                    freq = np.array([freq])

                # define a new set of colors for the beam patterns
                newmap = plt.get_cmap('autumn')
                desat = 0.7
                ax.set_color_cycle([newmap(k) for k in desat*np.linspace(0, 1, len(freq))])


                phis = np.arange(360) * 2 * np.pi / 360.
                newfreq = np.zeros(freq.shape)
                H = np.zeros((len(freq), len(phis)), dtype=complex)
                for i, f in enumerate(freq):
                    newfreq[i], H[i] = self.micArray.response(phis, f)

                # normalize max amplitude to one
                H = np.abs(H)**2/np.abs(H).max()**2

                # a normalization factor according to room size
                norm = np.linalg.norm((self.corners - self.micArray.center), axis=0).max()

                # plot all the beam patterns
                i = 0
                for f, h in zip(newfreq, H):
                    x = np.cos(phis) * h * norm + self.micArray.center[0, 0]
                    y = np.sin(phis) * h * norm + self.micArray.center[1, 0]
                    l = ax.plot(x, y, '-', linewidth=0.5)
                    # lbl = '%.2f' % f
                    # i0 = i*360/len(freq)
                    # ax.text(x[i0], y[i0], lbl, color=plt.getp(l[0], 'color'))
                    # i += 1

                #ax.legend(freq)

        # define some markers for different sources and colormap for damping
        markers = ['o', 's', 'v', '.']
        cmap = plt.get_cmap('YlGnBu')
        # draw the scatter of images
        for i, source in enumerate(self.sources):
            # draw source
            ax.scatter(
                source.position[0],
                source.position[1],
                c=cmap(1.),
                s=20,
                marker=markers[
                    i %
                    len(markers)],
                edgecolor=cmap(1.))
            #ax.text(source.position[0]+0.1, source.position[1]+0.1, str(i))

            # draw images
            if (img_order is None):
                img_order = self.max_order

            I = source.orders <= img_order

            val = (np.log2(source.damping[I]) + 10.) / 10.
            # plot the images
            ax.scatter(source.images[0, I], source.images[1, I], \
                       c=cmap(val), s=20,
                       marker=markers[i % len(markers)], edgecolor=cmap(val))

        return fig, ax

    def plotRIR(self, FD=False):

        if self.rir is None:
            self.compute_RIR()

        import matplotlib.pyplot as plt
        import utilities as u

        M = self.micArray.M
        S = len(self.sources)
        for r in xrange(M):
            for s in xrange(S):
                h = self.rir[r][s]
                plt.subplot(M, S, r*S + s + 1)
                if not FD:
                    plt.plot(np.arange(len(h)) / float(self.fs), h)
                else:
                    u.real_spectrum(h)
                plt.title('RIR: mic'+str(r)+' source'+str(s))
                if r == M-1:
                    if not FD:
                        plt.xlabel('Time [s]')
                    else:
                        plt.xlabel('Normalized frequency')

    def addMicrophoneArray(self, micArray):
        self.micArray = micArray

    def addSource(self, position, signal=None, delay=0):

        if (not self.isInside(np.array(position))):
            raise NameError('The source must be added inside the room.')

        # generate first order images
        i, d, w = self.firstOrderImages(np.array(position))
        images = [i]
        damping = [d]
        generators = [-np.ones(i.shape[1])]
        wall_indices = [w]

        # generate all higher order images up to max_order
        o = 1
        while o < self.max_order:
            # generate all images of images of previous order
            img = np.zeros((self.dim, 0))
            dmp = np.array([])
            gen = np.array([])
            wal = np.array([])
            for ind, si, sd in zip(xrange(len(images[o-1])), images[o - 1].T, damping[o - 1]):
                i, d, w = self.firstOrderImages(si)
                img = np.concatenate((img, i), axis=1)
                dmp = np.concatenate((dmp, d * sd))
                gen = np.concatenate((gen, ind*np.ones(i.shape[1])))
                wal = np.concatenate((wal, w))

            # sort
            ordering = np.lexsort(img)
            img = img[:, ordering]
            dmp = dmp[ordering]
            gen = gen[ordering]
            wal = wal[ordering]
                
            # add to array of images
            images.append(img)
            damping.append(dmp)
            generators.append(gen)
            wall_indices.append(wal)

            # next order
            o += 1

        # linearize the arrays
        images_lin = np.concatenate(images, axis=1)
        damping_lin = np.concatenate(damping)

        o_len = np.array([x.shape[0] for x in generators])
        # correct the pointers for linear structure
        for o in np.arange(2, len(generators)):
            generators[o] += np.sum(o_len[0:o-2])
            
        # linearize
        generators_lin = np.concatenate(generators)
        walls_lin = np.concatenate(wall_indices)

        # store the corresponding orders in another array
        ordlist = []
        for o in xrange(len(generators)):
            ordlist.append((o+1)*np.ones(o_len[o]))
        orders_lin = np.concatenate(ordlist)

        # add the direct source to the arrays
        images_lin = np.concatenate((np.array([position]).T, images_lin), axis=1)
        damping_lin = np.concatenate(([1], damping_lin))
        generators_lin = np.concatenate(([np.nan], generators_lin+1))
        walls_lin = np.concatenate(([np.nan], walls_lin))
        orders_lin = np.concatenate(([0], orders_lin))

        # add a new source to the source list
        self.sources.append(
            SoundSource(
                position,
                images=images_lin,
                damping=damping_lin,
                generators=generators_lin,
                walls=walls_lin,
                orders=orders_lin,
                signal=signal,
                delay=delay))

    def firstOrderImages(self, source_position):

        # projected length onto normal
        ip = np.sum(self.normals * (self.corners - source_position[:, np.newaxis]), axis=0)

        # projected vector from source to wall
        d = ip * self.normals

        # compute images points, positivity is to get only the reflections outside the room
        images = source_position[:, np.newaxis] + 2 * d[:, ip > 0]

        # collect absorption factors of reflecting walls
        damping = self.absorption[ip > 0]

        # collect the index of the wall corresponding to the new image
        wall_indices = np.arange(len(self.walls))[ip > 0]

        return images, damping, wall_indices

    def compute_RIR(self, c=constants.c):
        """
        Compute the room impulse response between every source and microphone
        """
        self.rir = []

        for mic in self.micArray.R.T:

            h = []

            for source in self.sources:

                h.append(source.getRIR(mic, self.fs, self.t0))

            self.rir.append(h)

    def simulate(self, recompute_rir=False):
        """Simulate the microphone signal at every microphone in the array"""

        # import convolution routine
        from scipy.signal import fftconvolve

        # Throw an error if we are missing some hardware in the room
        if (len(self.sources) is 0):
            raise NameError('There are no sound sources in the room.')
        if (self.micArray is None):
            raise NameError('There is no microphone in the room.')

        # compute RIR if necessary
        if len(self.rir) == 0 or recompute_rir:
            self.compute_RIR()

        # number of mics and sources
        M = self.micArray.M
        S = len(self.sources)

        # compute the maximum signal length
        from itertools import product
        max_len_rir = np.array([len(self.rir[i][j])
                                for i, j in product(xrange(M), xrange(S))]).max()
        f = lambda i: len(
            self.sources[i].signal) + np.floor(self.sources[i].delay * self.fs)
        max_sig_len = np.array([f(i) for i in xrange(S)]).max()
        L = max_len_rir + max_sig_len - 1
        if L % 2 == 1:
            L += 1

        # the array that will receive all the signals
        self.micArray.signals = np.zeros((M, L))

        # compute the signal at every microphone in the array
        for m in np.arange(M):
            rx = self.micArray.signals[m]
            for s in np.arange(S):
                sig = self.sources[s].signal
                if sig is None:
                    continue
                d = np.floor(self.sources[s].delay * self.fs)
                h = self.rir[m][s]
                rx[d:d + len(sig) + len(h) - 1] += fftconvolve(h, sig)

            # add white gaussian noise if necessary
            if self.sigma2_awgn is not None:
                rx += np.random.normal(0., np.sqrt(self.sigma2_awgn), rx.shape)

    def dSNR(self, x, source=0):
        """direct Signal-to-Noise Ratio"""

        if source >= len(self.sources):
            raise NameError('No such source')

        if self.sources[source].signal is None:
            raise NameError('No signal defined for source ' + str(source))

        if self.sigma2_awgn is None:
            return float('inf')

        x = np.array(x)

        sigma2_s = np.mean(self.sources[0].signal**2)

        d2 = np.sum((x - self.sources[source].position)**2)

        return sigma2_s/self.sigma2_awgn/(16*np.pi**2*d2)

    def getWallFromName(self, name):
        """
        Returns the instance of the wall by giving its name.
        
        :arg name: (string) name of the wall
        
        :returns: (Wall) instance of the wall with this name
        """
        
        if (name in self.wallsId):
            return self.walls[self.wallsId[name]]
        else:
            raise NameError('Room.getWallFromName : the wall '+name+' cannot be found.')
        
    def checkVisibilityForAllImages(self, source, p):
        """
        Checks visibility from a given point for all images of the given source.
        
        This function tests visibility for all images of the source and returns the results
        in an array.
        
        :arg source: (SoundSource) the sound source object (containing all its images)
        :arg p: (np.array size 2 or 3) coordinates of the point where we check visibility
        
        :returns: (int array) list of results of visibility for each image
            -1 : unchecked (only during execution of the function)
            0 (False) : not visible
            1 (True) : visible
        """
        
        visibilityCheck = np.zeros_like(source.images[0])-1
        
        for imageId in range(len(visibilityCheck)-1, -1, -1):
            visibilityCheck[imageId] = isVisible(source, p, imageId)
            
        return visibilityCheck
            
    def isVisible(self, source, p, imageId = 0):
        """
        Returns true if the given sound source (with image source id) is visible from point p.
        
        :arg source: (SoundSource) the sound source (containing all its images)
        :arg p: (np.array size 2 or 3) coordinates of the point where we check visibility
        :arg imageId: (int) id of the image within the SoundSource object
        
        :return: (bool)
            False (0) : not visible
            True (1) :  visible
        """

        p = np.array(p)

        if (source.order[imageId] > 0):
        
            # Check if the line of sight intersects the generating wall
            genWallId = source.walls[imageId]
            if self.walls[genWallId].intersects(p, np.array([source.images[:, imageId]])):

                    # Check if there is an obstruction
                    if(not self.isObstructed(source, p, imageId)):
                    
                        # Check visibility for the parent image by recursion
                        return isVisible(source, self.walls[genWallId].intersection(p, np.array([source.images[:, imageId]])), source.generators[imageId])
                    else:
                        return False
            else:
                return False
        else:
            return True
      
    def isObstructed(self, source, p, imageId = 0):
        """
        Checks if there is a wall obstructing the line of sight going from a source to a point.
        
        :arg source: (SoundSource) the sound source (containing all its images)
        :arg p: (np.array size 2 or 3) coordinates of the point where we check obstruction
        :arg imageId: (int) id of the image within the SoundSource object
        
        :returns: (bool)
            False (0) : not obstructed
            True (1) :  obstructed
        """
        
        #TODO optimization : use only walls inside the convex hull
        
        imageSide = self.walls[source.walls[imageId]].side(source.images[:, imageId])
        
        for wallId in range(len(self.walls)):
        
            # The generating wall can't be obstructive
            if(wallId != source.walls[imageId]):
            
                # Test if the line segment intersects the current wall
                if (self.walls[wallId].intersects(source.images[:, imageId], p)):
                    
                    # Only images with order > 0 have a generating wall. At this point, there is obstruction for source order 0.
                    if (source.orders[imageId] > 0):
                    
                        # Test if the intersection point and the image are at opposite sides of the generating wall
                        intersectionPoint = self.walls[wallId].intersection(source.images[:, imageId], p)
                        intersectionPointSide = self.walls[source.walls[imageId]].side(intersectionPoint)
                        if (intersectionPointSide != sourceSide):
                            return True
                    else:
                        return True
                
        return False

    def isInside(self, p, includeBorders = False):
        """
        Checks if the given point is inside the room.
        
        :arg p: (np.array dim 2 or 3) point to be tested
        :arg includeBorders: (bool) set true if a point on the wall must be considered inside the room
        
        :returns: (bool) True if the given point is inside the room, False otherwise.
        """
        
        p = np.array(p)
        if (self.dim != p.shape[0]):
            raise NameError('Room.isInside input error : dimension of room and p must match.')
        
        # Compute p0, which is a point outside the room at x coordinate xMin-1
        # (where xMin is the minimum x coordinate among the corners of the walls)
        p0 = np.array([np.amin(np.array([wall.corners[0, :] for wall in self.walls]).flatten())-1, p[0]])
        
        limitCase = False
        lastIntersection = 0
        count = 0
        for i in range(len(self.walls)):
            lastIntersection = self.walls[i].intersects(p0, p)
            if (lastIntersection == 2):
                limitCase = True
            if (lastIntersection > 0 and lastIntersection <= 3):
                count += 1
        if ((not limitCase and count % 2 == 1) or (limitCase and includeBorders)):
            return True
        else:
            return False