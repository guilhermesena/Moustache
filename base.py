# Class: SplineBase
# A spline base with k data points
# Has k+4 elements in the basis, 
# Which are 1, t, t^2, t^3 and
# (t-t_j)^3+ 

class SplineBase:
    def __init__(self, points):
        self.points = points
        points.sort()
    
    def apply(self, val, index):
        if index <= 3:
            return val**index

        return max((val - self.points[index - 4])**3, 0)
    

if __name__ == "__main__":
    points = []
    for i in range(11):
        points.append(i/float(10))

    print (points)
    s = SplineBase(points)

    for i in range(14):
        print (s.apply(0.3, i))

