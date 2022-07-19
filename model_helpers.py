import numpy as np

from sklearn.linear_model import LinearRegression


class LogitRegression(LinearRegression):
    def fit(self, x, p):
        eps = 1e-6
        p = np.asarray(p)
        p[p >= 1 - eps] = 1 - eps
        p[p <= eps] = eps
        y = np.log(p / (1 - p))
        return super().fit(x, y)

    def predict(self, x):
        y = super().predict(x)
        return 1 / (np.exp(-y) + 1)