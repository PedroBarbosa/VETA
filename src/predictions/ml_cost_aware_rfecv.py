import numpy as np

from sklearn.utils import check_X_y
from sklearn.feature_selection import RFE
from sklearn.model_selection import cross_val_score

from yellowbrick.features import RFECV as ybRFECV
from yellowbrick.base import ModelVisualizer
from yellowbrick.exceptions import YellowbrickValueError

class CostAwareRFECV(ybRFECV):
    """Expansion of RFECV from yellowbrick visualization library"""
    def __init__(self, model, ax=None, step=1, groups=None, cv=None,
                 scoring=None, cost=0, **kwargs):
        super(CostAwareRFECV, self).__init__(model, ax, step, groups, cv, scoring, **kwargs)
        self.cost_per_feature = cost
        
    def fit(self, X, y=None):
        check_X_y(X, y, "csr")
        n_features = X.shape[1]

        # This check is kind of unnecessary since RFE will do it, but it's
        # nice to get it out of the way ASAP and raise a meaningful error.
        if 0.0 < self.step < 1.0:
            step = int(max(1, self.step * n_features))
        else:
            step = int(self.step)

        if step < 0:
            raise YellowbrickValueError("step must be >0")

        # Create the RFE model
        rfe = RFE(self.estimator, step=step)
        n_feature_subsets = np.arange(1, n_features+1)

        # Create the cross validation params
        # TODO: handle random state
        cv_params = {
            key: self.get_params()[key]
            for key in ('groups', 'cv', 'scoring')
        }

        # Perform cross-validation for each feature subset
        scores = []
        for n_features_to_select in n_feature_subsets:
            rfe.set_params(n_features_to_select=n_features_to_select)
            scores.append(cross_val_score(rfe, X, y, **cv_params))

        # Convert scores to array
        self.cv_scores_ = np.array(scores)

        # Find the best RFE model
        best_ind = 0
        best_score = 0
        ongoing_cost = 0
        for ind, score in enumerate(self.cv_scores_.mean(axis=1)):
            if score > best_score * (1 + self.cost_per_feature):
                best_score = score
                best_ind = ind
                ongoing_cost = self.cost_per_feature
            else:
                ongoing_cost += self.cost_per_feature
        bestidx = best_ind
        #bestidx = self.cv_scores_.mean(axis=1).argmax()
        self.n_features_ = n_feature_subsets[bestidx]

        # Fit the final RFE model for the number of features
        self.rfe_estimator_ = rfe
        self.rfe_estimator_.set_params(n_features_to_select=self.n_features_)
        self.rfe_estimator_.fit(X, y)

        # Rewrap the visualizer to use the rfe estimator
        self._wrapped = self.rfe_estimator_

        # Hoist the RFE params to the visualizer
        self.support_ = self.rfe_estimator_.support_
        self.ranking_ = self.rfe_estimator_.ranking_

        self.draw()
        return self
