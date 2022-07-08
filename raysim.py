import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from time import sleep, time

print("hi")


def roots(coefficients):  # TODO: VECTORIZE
    roots = np.zeros([2, np.shape(coefficients)[1]])
    for n in range(np.shape(coefficients)[1]):
        roots[:, n] = np.roots(coefficients[:, n])
    return roots


def arrdot(a, b, size):
    p0 = np.reshape(a, [2 * size])[:size]
    p1 = np.reshape(a, [2 * size])[size:]
    p2 = np.reshape(b, [2 * size])[:size]
    p3 = np.reshape(b, [2 * size])[size:]
    return p0 * p2 + p1 * p3


def temp(z):
    if z < 11000:
        temperature = 288.14 - 0.00649 * z
    elif z > 25000:
        temperature = 141.89 + 0.00299 * z
    else:
        temperature = 216.64
    return temperature


def press(z):
    if z < 11000:
        pressure = 101290 * (temp(z) / 288.08) ** 5.256
    elif z > 25000:
        pressure = 2488 * (temp(z) / 216.6) ** (-11.388)
    else:
        pressure = 22650 * np.exp(1.73 - 0.000157 * z)
    return pressure


def compress(p, t, rh):  # work out compressibility factor
    t -= 273.1
    p /= 10 ** 6
    part1 = 1.00001 - 5.8057 * 10 ** (-3) * p + 2.6402 * 10 ** (-4) * p ** 2
    part2 = - 3.3297 * 10 ** (-7) * t + 1.2420 * 10 ** (-4) * p * t - 2.0158 * 10 ** (-6) * p ** 2 * t
    part3 = + 2.4925 * 10 ** (-9) * t ** 2 - 6.2873 * 10 ** (-7) * p * t ** 2 + 5.4174 * 10 ** (-9) * p ** 2 * t ** 2
    part4 = - 3.5 * 10 ** (-7) * rh - 5 * 10 ** (-9) * rh ** 2
    return part1 + part2 + part3 + part4


def edlen(z, rh, wlen):
    p = press(z)/100
    t = temp(z)
    wlen /= 1000
    sigma = 1/wlen
    part1 = (77.497 + 0.013) * (p / t) * (1 / compress(p, t, rh))
    part2 = (0.306007 + 88.2581 / (130 - sigma ** 2) + 0.5868 / (38.9 - sigma ** 2))
    part3 = part1 * part2
    return (part3 / 10 ** (6)) + 1


X = np.linspace(0, 85000, 200)
Y = X
for n in range(200):
    Y[n] = edlen(X[n], 0, 500)
plt.scatter(np.linspace(0, 85000, 200), Y)
plt.show()


earth_rad = 6731 * 1000
at_rad = 85 * 1000
tot_rad = earth_rad + at_rad


class Source:
    def __init__(self, angle, diameter, raynum, rings, refracs):
        source_vec = np.zeros([2, raynum])
        source_vec[0, :] = np.linspace(0, diameter, raynum) * np.cos(angle + np.pi / 2)
        source_vec[1, :] = np.linspace(0, diameter, raynum) * np.sin(angle + np.pi / 2)

        source_vec[1, :] -= 0.5 * diameter * np.cos(angle)  # positioning the source
        source_vec[1, :] -= 3 * diameter * np.sin(angle)
        source_vec[0, :] += 0.5 * diameter * np.sin(angle)
        source_vec[0, :] -= 3 * diameter * np.cos(angle)
        self.source_vec = source_vec
        self.angle = np.full(raynum, angle)
        self.raynum = raynum

        self.collided = np.zeros([raynum, rings], dtype='int')  # has the ray entered a circle yet
        self.entered = np.zeros([raynum, rings], dtype='int')

        self.refrac_arr = np.zeros([rings, 2, raynum])

        self.absorb = np.zeros([rings + 1])
        self.absorb[-1] = 1
        for n in range(rings):  # create array of each refractive index for each ring
            self.refrac_arr[n, :, :] = np.full((raynum, 2), [refracs[n], refracs[n + 1]]).T
            print(n, self.refrac_arr[n, :, :])

        self.absorbed = np.zeros([raynum])
        self.diameter = diameter
        self.rings = rings

    def calc_line(self, count, output, vis):
        self.grad = np.tan(self.angle)
        self.const = - self.grad * self.source_vec[0, :] + self.source_vec[1, :]
        if count < self.rings * 2:
            count = self.diameter * 1000
        if vis:
            #  visualise the rays
            if count > -1:
                for i in range(self.raynum):
                    if self.absorbed[i] == 0:
                        X = np.linspace(self.source_vec[0, i], count, 200)
                        #print("vis")
                        plt.plot(X, self.grad[i] * X + self.const[i])
        # const is the c in y = mx + c
        if output:
            return [self.grad, self.const, self.absorbed]

    def intersect(self, count, radius, vis):
        a = self.grad ** 2 + 1
        b = 2 * self.grad * self.const
        c = self.const ** 2 - radius ** 2
        coef = np.array([a, b, c])

        discriminant = b ** 2 - 4 * a * c
        intersected = np.array([discriminant > 0], dtype='int')[0]  # does a ray interesect?
        intersections = roots(coef)
        intersections = np.sort(intersections, axis=0)
        intersection = np.diag(intersections[self.collided[:, count]])  # 1st intersection if hasnt collided yet
        intersection_y = self.grad * intersection + self.const  # y values!
        intersection = np.array([intersection, intersection_y])
        self.source_vec += intersected * (intersection - self.source_vec)

        self.collided[:, count] = intersected  # TODO: make sure these are x and y coords

        if vis:
            print(radius, "radius")
            #plt.scatter(source0.source_vec[0], source0.source_vec[1], c='black', s=2)
            plt.scatter(np.linspace(-radius, radius, 1000), np.sqrt(radius ** 2 - np.linspace(-radius, radius, 1000) ** 2), s=1)
            plt.scatter(np.linspace(-radius, radius, 1000), - np.sqrt(radius ** 2 - np.linspace(-radius, radius, 1000) ** 2), s=1)

    def refract(self, count):
        n1 = np.diag(self.refrac_arr[count, :, :][self.entered][:, count])
        n2 = np.diag(self.refrac_arr[count, :, :][1 - self.entered][:, count])

        if bool(self.absorb[count + 1]):
            self.absorbed = self.entered[:, count]  # if it collided then it was absorbed

        is_negative = self.source_vec[0] / abs(self.source_vec[0])
        phi = np.arctan(self.source_vec[1]/self.source_vec[0]) + (np.pi) * ((1 - is_negative) / 2)

        theta_1 = phi - self.angle - np.pi * ((1 - is_negative) / 2)
        theta_2 = np.arcsin((n1 * np.sin(theta_1)) / n2)
        #print(theta_2, "theta_2")

        theta_o = self.angle + (theta_1 - theta_2) * (self.collided[:, count])
        #print(theta_o, "theta_o")
        self.angle = theta_o
        self.entered[:, count] = self.collided[:, count]


class Sensor:
    def __init__(self, x_pos, height):
        self.x_pos = x_pos
        self.height = height

    def image(self, in_grad, in_const, absorbed):
        #print(np.sum(absorbed), "absorbed")
        y_pos = in_grad * self.x_pos + in_const
        remaining = y_pos[np.array([absorbed == 0])[0]]
        less = remaining[np.array([remaining < self.height])[0]]
        final = less[np.array([less > - self.height])[0]]
        return final


shell_num = 5

ring_rad = at_rad / shell_num

in_refracs = np.ones([shell_num + 1])
for i in range(1, shell_num + 1):
    in_refracs[i] = edlen(at_rad - ring_rad * i, 0, 500)

source0 = Source(0.0, tot_rad * 2, 5000, shell_num, in_refracs)

t = time()
shell_ref = np.append(np.arange(shell_num), np.arange(shell_num)[::-1])
for n in range(shell_num * 2):
    source0.calc_line(n, False, vis=False)  # don't return these ray lines
    source0.intersect(shell_ref[n], tot_rad - ring_rad * shell_ref[n], vis=False)
    source0.refract(shell_ref[n])
rays = source0.calc_line(shell_num + 1, True, vis=False)  # final ray lines
print(time() - t)


fig = plt.figure()
ax = fig.gca()
def focus(i):
    if i == 0:
        sleep(1)
    ax.clear()
    sensor = Sensor(tot_rad + tot_rad * 1000 * (i + 900), tot_rad)
    image = sensor.image(rays[0], rays[1], rays[2])
    #image1 = sensor.image(rays1[0], rays1[1])
    ax.hist(image, bins=10)
    #ax.hist(image1, bins=50)
    print(tot_rad + tot_rad * 1000 * (i + 900), i)


anim = ani.FuncAnimation(fig, focus, interval=50)
plt.show()
'''

plt.show()
'''
