import numpy as np
import matplotlib.pyplot as plt

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


class source:
    def __init__(self, angle, diameter, raynum, rings):
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

        self.collided = np.zeros([raynum], dtype='int')  # has the ray entered a circle yet
        self.entered = np.zeros([raynum], dtype='int')

        self.refrac = 1.1
        self.refrac_arr = np.full((raynum, 2), [1, self.refrac]).T

    def calc_line(self, count):
        self.grad = np.tan(self.angle)
        self.const = - self.grad * self.source_vec[0, :] + self.source_vec[1, :]
        if count == 2:
            count = 15
        #'''  #visualise the rays
        if count > -1:
            for i in range(self.raynum):
                X = np.linspace(self.source_vec[0, i], count, 200)
                plt.plot(X, self.grad[i] * X + self.const[i])
        #'''
        # const is the c in y = mx + c

    def intersect(self, radius):
        a = self.grad ** 2 + 1
        b = 2 * self.grad * self.const
        c = self.const ** 2 - radius ** 2
        coef = np.array([a, b, c])

        discriminant = b ** 2 - 4 * a * c
        intersected = np.array([discriminant > 0], dtype='int')[0]  # does a ray interesect?
        intersections = roots(coef)
        intersections = np.sort(intersections, axis=0)
        intersection = np.diag(intersections[self.collided])  # 1st intersection if hasnt collided yet
        intersection_y = self.grad * intersection + self.const  # y values!
        intersection = np.array([intersection, intersection_y])
        self.source_vec += intersected * (intersection - self.source_vec)

        self.collided = intersected  # TODO: make sure these are x and y coords

    def refract(self):
        n1 = np.diag(self.refrac_arr[self.entered])
        n2 = np.diag(self.refrac_arr[1 - self.entered])  # TODO: FIX THIS (not generalisable for shells > 1)

        is_negative = self.source_vec[0] / abs(self.source_vec[0])
        print(is_negative)
        phi = np.arctan(self.source_vec[1]/self.source_vec[0]) + (np.pi) * ((1 - is_negative) / 2)
        print(np.degrees(phi)-180, "phi")

        theta_1 = phi - self.angle - np.pi * ((1 - is_negative) / 2)
        theta_2 = np.arcsin((n1 * np.sin(theta_1)) / n2)

        theta_o = self.angle + (theta_1 - theta_2) * self.collided
        self.angle = theta_o
        self.entered = self.collided


shell_num = 1
source0 = source(0.1, 1.9, 20, shell_num)

for n in range(shell_num + 1):
    source0.calc_line(n)
    source0.intersect(1)
    source0.refract()
source0.calc_line(2)

plt.scatter(source0.source_vec[0], source0.source_vec[1], c='black')
plt.scatter(np.linspace(-1, 1, 100), np.sqrt(1 - np.linspace(-1, 1, 100) ** 2))
plt.scatter(np.linspace(-1, 1, 100), - np.sqrt(1 - np.linspace(-1, 1, 100) ** 2))
plt.show()
