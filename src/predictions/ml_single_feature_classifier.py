from sklearn import dummy
from src.predictions.ml_prepare_dataset import MISSING_CONSTANT

class SingleFeatureClassifier(dummy.DummyClassifier):

    def __init__(self, tool, direction, threshold, index):
        super(SingleFeatureClassifier, self).__init__()
        self.tool = tool
        self.direction = direction
        self.threshold = float(threshold)
        self.index = index

    def meaning(self,x):
        if x == MISSING_CONSTANT:
            return 3
        elif self.direction == ">":
            return True if x >= self.threshold else False
        else:
            return True if x < self.threshold else False

    def fit(self, X, y):
        return self

    def predict(self, X, y=None):
        #print(X.shape, self.index, "bb")
        return [self.meaning(x[self.index]) for x in X]
