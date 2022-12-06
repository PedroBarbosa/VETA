import logging
from typing import Union
import numpy as np
from gplearn.genetic import SymbolicClassifier
from sklearn import tree
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.exceptions import NotFittedError, ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, cross_validate
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import Normalizer, StandardScaler
from sklearn.svm import LinearSVC, SVC
from sklearn.utils._testing import ignore_warnings


class Classifiers(object):
    """
    Class to test multiple classifiers using different hyperparameters
    """

    def __init__(self, x: np.array, y: np.array, scoring: dict, refit_metric: str):
        """
        Initialization with the scikit learn scoring metrics desired

        :param np.array x: Training data
        :param np.array y: Labels

        :param dict scoring: Dict with predefined metrics. Keys must
            be metric names, while values should be predefined metric
            string recognized by sklearn

        :param srt refit_metric: Refers to the metric to use for model
            selection after best parameter selection. Then, the model
            is refitted using all the data and the selected model. Value
            must be present in `scoring` keys.
        """
        self.x = x
        self.y = y.reshape(-1)
        self.scoring = scoring
        assert refit_metric in self.scoring.keys(), "Refit metric must exist in scoring keys."
        self.refit_metric = refit_metric

    def knn(self) -> Pipeline:
        """
        Creates a pipeline for a K-nearest neighbors classifier

        Since it uses distances to learn, a Normalization preprocessing
        step is also included.

        :return Pipeline: returns a Pipeline for the best estimator
        """

        pipeline = Pipeline(steps=[('normalization', Normalizer()),
                                   ('knn', KNeighborsClassifier())])
        params_grid = {'normalization__norm': ['l1', 'l2'],
                       'knn__n_neighbors': [3, 5, 10]}
        return self.do_grid_search("knn", pipeline, params_grid)

    def naive_bayes(self) -> Pipeline:
        """
        Creates a pipeline for a Gaussian Naive Bayes classifier

        Since it uses distances to learn, a Normalization preprocessing
        step is also included.

        :return Pipeline: returns a Pipeline for the fitted estimator
        """
        pipeline = Pipeline(steps=[('NaiveBayes', GaussianNB())])
        return self.do_grid_search("NaiveBayes", pipeline)

    @ignore_warnings(category=NotFittedError)
    def linear_discrimant_analysis(self) -> Pipeline:
        """
        Creates a pipeline for a Linear Discriminant analysis

        :return Pipeline: returns a Pipeline for the best estimator
        """

        pipeline = Pipeline(steps=[('lda', LinearDiscriminantAnalysis())])
        params_grid = {'lda__solver': ['svd', 'lsqr']}
        return self.do_grid_search("lda", pipeline, params_grid, njobs=1)

    @ignore_warnings(category=UserWarning)
    def quadratic_discriminant_analysis(self) -> Pipeline:
        """
        Creates a pipeline for a Quadratic Discriminant analysis

        :return Pipeline: returns a Pipeline for the best estimator
        """
        pipeline = Pipeline(steps=[('qda', QuadraticDiscriminantAnalysis())])
        params_grid = {'qda__reg_param': [0.0, 0.2, 0.4, 0.6, 0.8, 1],
                       'qda__store_covariance': [True]}

        return self.do_grid_search("qda", pipeline, params_grid, njobs=1)

    def decision_tree(self) -> Pipeline:
        """
        Creates a pipeline for a Decision Tree

        :return Pipeline: returns a Pipeline for the best estimator
        """
        pipeline = Pipeline(steps=[('dt', tree.DecisionTreeClassifier(random_state=0))])

        params_grid = {'dt__criterion': ['gini', 'entropy'],
                       'dt__splitter': ['best', 'random'],
                       'dt__max_depth': [3, 5, 7],
                       'dt__min_samples_split': [2, 4, 6],
                       'dt__min_samples_leaf': [1, 2, 3, 4],
                       'dt__max_features': [None, "auto", "sqrt"]}

        return self.do_grid_search("dt", pipeline, params_grid)

    def random_forest(self) -> Pipeline:
        """
        Creates a pipeline for a Random Forest

        :return Pipeline: returns a Pipeline for the best estimator
        """
        pipeline = Pipeline(steps=[('rf', RandomForestClassifier(random_state=0))])

        params_grid = {'rf__n_estimators': [10, 50, 100],
                       'rf__criterion': ['gini', 'entropy'],
                       'rf__max_depth': [5, 7],
                       'rf__min_samples_split': [2, 4],
                       'rf__min_samples_leaf': [2, 4],
                       'rf__max_features': ["auto"],
                       'rf__bootstrap': [True],
                       # 'rf__ccp_alpha': [0.0, 0.1],
                       'rf__max_samples': [None, 0.8, 0.9]}

        return self.do_grid_search("rf", pipeline, params_grid)

    def ada_boost(self) -> Pipeline:
        """
        Creates a pipeline for a AdaBoost

        :return Pipeline: returns a Pipeline for the best estimator
        """
        pipeline = Pipeline(steps=[('ab', AdaBoostClassifier(random_state=0))])

        params_grid = {'ab__n_estimators': [10, 50, 100],
                       'ab__learning_rate': [0.5, 1]
                       }

        return self.do_grid_search("svm", pipeline, params_grid)

    @ignore_warnings(category=ConvergenceWarning)
    def svm(self) -> Pipeline:
        """
        Creates a pipeline for Support Vector Machines

        :return Pipeline: returns a Pipeline for the best estimator
        """
        pipeline = Pipeline(steps=[('scaler', StandardScaler()),
                                   ('svm', LinearSVC(random_state=0))])

        # Pipeline that tests multiple SVM classifiers
        params_grid = [{'svm': (SVC(random_state=0),),
                        'svm__C': (0.001, 0.01, 0.1, 1, 10),
                        'svm__kernel': ('poly', 'rbf', 'sigmoid')},

                       {'svm': (LinearSVC(random_state=0),),
                        # 'svm__penalty': ('l1', 'l2'),
                        'svm__loss': ('hinge', 'squared_hinge'),
                        'svm__C': ([0.001, 0.01, 0.1, 1, 10]),
                        'svm__max_iter': (1500, 2000),
                        }
                       ]

        return self.do_grid_search("svm", pipeline, params_grid)

    def logistic_regression(self) -> Pipeline:
        """
        Creates a pipeline for Support Vector Machines

        :return Pipeline: returns a Pipeline for the best estimator
        """
        pipeline = Pipeline(steps=[('scaler', StandardScaler()),
                                   ('lr', LogisticRegression(random_state=0))])

        params_grid = {'lr__penalty': ['l1', 'l2'],
                       'lr__solver': ['liblinear']}

        return self.do_grid_search("lr", pipeline, params_grid)

    def gp(self) -> Pipeline:
        """
        Creates a pipeline for Genetic programming

        :return Pipeline: returns a Pipeline for the best estimator
        """
        pipeline = Pipeline(steps=[('scaler', StandardScaler()),
                                   ('gp', SymbolicClassifier(random_state=0))])

        params_grid = {'gp__generations': [10, 50, 100]}

        return self.do_grid_search("gp", pipeline, params_grid)

    def do_grid_search(self, name: str,
                       pipeline: Pipeline,
                       params_grid: Union[dict, list] = None,
                       njobs: int = -1) -> Pipeline:
        """
        Perform grid search for a given pipeline

        :param str name: Name of the estimator
        :param pipeline: Pipeline for a given estimator
        :param dict params_grid: Parameters to test in the
            grid search for the provided pipeline.
        :param int njobs: Number of jobs to run in parallel.
            Default: `-1`, use all available cpus.

        :return Pipeline: Returns a Pipeline for the best estimator
            so that it can be evaluated on test data (making predictions)
        """
        # If there are no parameters to tune, apply CV directly
        if params_grid is None:
            res = cross_validate(pipeline, self.x, self.y,
                                 cv=10,
                                 scoring=self.scoring,
                                 n_jobs=njobs)

            cv_score = np.mean(res['test_' + self.refit_metric])
            logging.info("Mean cross validated score ({}) for the {} estimator: {:.3f}".format(self.refit_metric,
                                                                                               name,
                                                                                               cv_score))
            # Fit all training data
            return pipeline.fit(self.x, self.y), None, cv_score

        else:
            # StratifiedKFold is performed for each of parameter combination
            # with ignore_warnings(category=[ConvergenceWarning, FitFailedWarning]):
            search = GridSearchCV(pipeline,
                                  params_grid,
                                  scoring=self.scoring,
                                  refit=self.refit_metric,
                                  cv=10,
                                  n_jobs=njobs,
                                  )

            search.fit(self.x, self.y)
            logging.info("Best parameters found for {} estimator: {}".format(name, search.best_params_))
            logging.info("Mean cross validated score ({}) of the best {} estimator: {:.3f}".format(self.refit_metric,
                                                                                                   name,
                                                                                                   search.best_score_))
            return search.best_estimator_, search.best_params_, search.best_score_

    def train(self, name: str) -> Pipeline:
        """
        Run pipeline for a given algorithm

        :param name:
        :return:
        """
        if name == "KNN":
            return self.knn()

        elif name == "Naive_Bayes":
            return self.naive_bayes()

        elif name == "LDA":
            return self.linear_discrimant_analysis()

        elif name == "QDA":
            return self.quadratic_discriminant_analysis()

        elif name == "Decision_Tree":
            return self.decision_tree()

        elif name == "Random_Forest":
            return self.random_forest()

        elif name == "Ada_Boost":
            return self.ada_boost()

        elif name == "SVM":
            return self.svm()

        elif name == "Logistic_Regression":
            return self.logistic_regression()

        elif name == "Genetic_Programming":
            return self.gp()
