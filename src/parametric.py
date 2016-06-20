# Class Parametric
# Different implementations of
# Parametrization from unordered
# data
class Parametric:
    def __init__(self, data):
        self.data = data

    def get_points(self):
        raise NotImplementedError()

class UniformParametric (Parametric):
    def __init__ (self, data, k):
        super().__init__(data)
        self.k = k

        self.points = []
        for i in range(k+1):
            self.points.append(i/float(k))

    def get_points(self):
        return self.points

if __name__ == "__main__":
    data = [[2, 3], [5,5], [1,8], [9,6]]
    prm = UniformParametric(data, 20)

    print (prm.get_points())
