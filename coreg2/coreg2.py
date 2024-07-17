import numpy as np
from time import time
from scipy.spatial.distance import minkowski
from sklearn.neighbors import KNeighborsRegressor
from sklearn.utils import shuffle
from sklearn.metrics import mean_squared_error
from data_utils import load_data
from sklearn.preprocessing import LabelEncoder
import pandas as pd
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split



class Coreg():
    """
    Instantiates a CoReg regressor.
    """
    def __init__(self, k1=3, k2=3, p1=2, p2=5, max_iters=100, pool_size=100):
        self.k1, self.k2 = k1, k2 # number of neighbors
        self.p1, self.p2 = p1, p2 # distance metrics
        self.max_iters = max_iters
        self.pool_size = pool_size
        self.h1 = KNeighborsRegressor(n_neighbors=self.k1, p=self.p1)
        self.h2 = KNeighborsRegressor(n_neighbors=self.k2, p=self.p2)
        self.h1_temp = KNeighborsRegressor(n_neighbors=self.k1, p=self.p1)
        self.h2_temp = KNeighborsRegressor(n_neighbors=self.k2, p=self.p2)
        self.labeled_data = None
        self.unlabeled_data = None

    def add_data(self, data_dir):
        """
        Adds data and splits into labeled and unlabeled.
        """
        df = pd.read_csv(data_dir)
        columns_to_drop = ["SDMVPSU", "SDMVSTRA", "WTMEC4YR", "SDDSRVYR", "RIDSTATR", "WTSCY4YR"]
        df = df.drop(columns=columns_to_drop)
        
        # Encode categorical variables
        label_encoders = {}
        for column in ["RIAGENDR", "RIDRETH2", "RIDEXPRG"]:
            le = LabelEncoder()
            df[column] = le.fit_transform(df[column].astype(str))
            label_encoders[column] = le
        
        # Separate labeled and unlabeled data
        self.labeled_data = df[df['simulated_mGFR'].notna()]
        self.unlabeled_data = df[df['simulated_mGFR'].isna()]
        self.unlabeled_data = self.unlabeled_data.drop(columns=['simulated_mGFR'])
        
        # to get rid of NaNs
        preprocessing_pipeline = Pipeline([
            ('imputer', SimpleImputer(strategy="mean")),
            ('std_scaler', StandardScaler())
        ])
        
        # Fit and transform labeled dataset
        self.labeled_data.iloc[:, :-1] = preprocessing_pipeline.fit_transform(self.labeled_data.iloc[:, :-1])
        self.labeled_data = pd.DataFrame(self.labeled_data)
        
        # Transform unlabeled dataset
        self.unlabeled_data = preprocessing_pipeline.transform(self.unlabeled_data)
        self.unlabeled_data = pd.DataFrame(self.unlabeled_data)

        print("labeled_data shape: ", self.labeled_data.shape)
        print("unlabeled_data shape: ", self.unlabeled_data.shape)

        self._split_data()

    def _split_data(self):
        """
        Splits the data into training, test, and unlabeled sets.
        """
        if self.labeled_data is None or self.unlabeled_data is None:
            raise ValueError("Data has not been added yet. Call add_data() first.")
        
        # Extract features and target
        y = self.labeled_data['simulated_mGFR'].values
        X = self.labeled_data.drop(columns=['simulated_mGFR']).values

        # Split labeled data into training and test sets
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, shuffle=True, random_state=42)

        # Initialize labeled sets for the two learners
        self.L1_X = X_train.copy()
        self.L1_y = y_train.copy()
        self.L2_X = X_train.copy()
        self.L2_y = y_train.copy()
        
        # Initialize test set and unlabeled set
        self.y_train = y_train
        self.X_test = X_test
        self.y_test = y_test
        self.X_train = X_train
        self.X_test = X_test
        self.U_X = self.unlabeled_data.values
        self.U_y = np.zeros((self.U_X.shape[0],))

    def run_trials(self, num_train=100, trials=10, verbose=False):
        """
        Runs multiple trials of training.
        """
        self.num_train = num_train
        self.num_trials = trials
        self._initialize_metrics()
        self.trial = 0
        while self.trial < self.num_trials:
            t0 = time()
            print(f'Starting trial {self.trial + 1}:')
            self.train(random_state=(self.trial + self.num_train), num_labeled=self.num_train, num_test=1000, verbose=verbose, store_results=True)
            print(f'Finished trial {self.trial + 1}: {time() - t0:.2f}s elapsed\n')
            self.trial += 1

    def train(self, random_state=-1, num_labeled=100, num_test=1000, verbose=False, store_results=False):
        """
        Trains the CoReg regressor.
        """
        t0 = time()
        self._split_data()
        self._fit_and_evaluate(verbose)
        if store_results: 
            self._store_results(0)
        self._get_pool()
        if verbose: 
            print(f'Initialized h1, h2: {time() - t0:.2f}s\n')
        for t in range(1, self.max_iters + 1):
            stop_training = self._run_iteration(t, t0, verbose, store_results)
            if stop_training:
                if verbose:
                    print(f'Done in {t} iterations: {time() - t0:.2f}s')
                break
        if verbose: 
            print(f'Finished {t} iterations: {time() - t0:.2f}s')

    def _run_iteration(self, t, t0, verbose=False, store_results=False):
        """
        Run t-th iteration of co-training, returns stop_training=True if no more unlabeled points are added to label sets.
        """
        stop_training = False
        if verbose: 
            print(f'Started iteration {t}: {time() - t0:.2f}s')
        self._find_points_to_add()
        added = self._add_points()
        if added:
            self._fit_and_evaluate(verbose)
            if store_results:
                self._store_results(t)
            self._remove_from_unlabeled()
            self._get_pool()
        else:
            stop_training = True
        return stop_training

    def _add_points(self):
        """
        Adds new examples to training sets.
        """
        added = False
        if self.to_add['x1'] is not None:
            print("self.L2_X shape: ", self.L2_X.shape)
            print("self.L2_y shape: ", self.L2_y.shape)
            self.L2_y = self.L2_y.reshape(-1, 1)
            self.L2_X = np.vstack((self.L2_X, self.to_add['x1']))
            self.L2_y = np.vstack((self.L2_y, self.to_add['y1']))
            added = True
        if self.to_add['x2'] is not None:
            self.L1_y = self.L1_y.reshape(-1, 1)
            self.L1_X = np.vstack((self.L1_X, self.to_add['x2']))
            self.L1_y = np.vstack((self.L1_y, self.to_add['y2']))
            added = True
        return added

    def _compute_delta(self, omega, L_X, L_y, h, h_temp):
        """
        Computes the improvement in MSE among the neighbors of the point being evaluated.
        """
        delta = 0
        for idx_o in omega:
            delta += (L_y[idx_o].reshape(1, -1) - h.predict(L_X[idx_o].reshape(1, -1))) ** 2
            delta -= (L_y[idx_o].reshape(1, -1) - h_temp.predict(L_X[idx_o].reshape(1, -1))) ** 2
        return delta

    def _compute_deltas(self, L_X, L_y, h, h_temp):
        """
        Computes the improvements in local MSE for all points in pool.
        """
        deltas = np.zeros((self.U_X_pool.shape[0],))
        for idx_u, x_u in enumerate(self.U_X_pool):
            x_u = x_u.reshape(1, -1)
            y_u_hat = h.predict(x_u).reshape(1, -1)
            omega = h.kneighbors(x_u, return_distance=False)[0]
            # print("shape of L_x: ", L_X.shape)
            # print("shape of x_u: ", x_u.shape)
            X_temp = np.vstack((L_X, x_u))

            L_y = L_y.reshape(-1, 1)
            # print("shape of L_y: ", L_y.shape)
            # print("shape of y_u_hat: ", y_u_hat.shape)
            y_temp = np.vstack((L_y, y_u_hat))
            h_temp.fit(X_temp, y_temp)
            delta = self._compute_delta(omega, L_X, L_y, h, h_temp)
            deltas[idx_u] = delta
        return deltas


    def _evaluate_metrics(self, verbose):
        """
        Evaluates KNN regressors on training and test data.
        """
        train1_hat = self.h1.predict(self.X_train)
        train2_hat = self.h2.predict(self.X_train)
        train_hat = 0.5 * (train1_hat + train2_hat)
        test1_hat = self.h1.predict(self.X_test)
        test2_hat = self.h2.predict(self.X_test)
        test_hat = 0.5 * (test1_hat + test2_hat)
        self.mse1_train = mean_squared_error(train1_hat, self.y_train)
        self.mse1_test = mean_squared_error(test1_hat, self.y_test)
        self.mse2_train = mean_squared_error(train2_hat, self.y_train)
        self.mse2_test = mean_squared_error(test2_hat, self.y_test)
        self.mse_train = mean_squared_error(train_hat, self.y_train)
        self.mse_test = mean_squared_error(test_hat, self.y_test)
        if verbose:
            print('MSEs:')
            print('  KNN1:')
            print('    Train: {:0.4f}'.format(self.mse1_train))
            print('    Test : {:0.4f}'.format(self.mse1_test))
            print('  KNN2:')
            print('    Train: {:0.4f}'.format(self.mse2_train))
            print('    Test : {:0.4f}'.format(self.mse2_test))
            print('  Combined:')
            print('    Train: {:0.4f}'.format(self.mse_train))
            print('    Test : {:0.4f}\n'.format(self.mse_test))

    def _find_points_to_add(self):
        """
        Finds unlabeled points (if any) to add to training sets.
        """
        self.to_add = {'x1': None, 'y1': None, 'idx1': None,
                       'x2': None, 'y2': None, 'idx2': None}
        # Keep track of added idxs
        added_idxs = []
        for idx_h in [1, 2]:
            if idx_h == 1:
                h = self.h1
                h_temp = self.h1_temp
                L_X, L_y = self.L1_X, self.L1_y
            elif idx_h == 2:
                h = self.h2
                h_temp = self.h2_temp
                L_X, L_y = self.L2_X, self.L2_y
            deltas = self._compute_deltas(L_X, L_y, h, h_temp)
            # Add largest delta (improvement)
            sort_idxs = np.argsort(deltas)[::-1] # max to min
            max_idx = sort_idxs[0]
            if max_idx in added_idxs: max_idx = sort_idxs[1]
            if deltas[max_idx] > 0:
                added_idxs.append(max_idx)
                x_u = self.U_X_pool[max_idx].reshape(1, -1)
                y_u_hat = h.predict(x_u).reshape(1, -1)
                self.to_add['x' + str(idx_h)] = x_u
                self.to_add['y' + str(idx_h)] = y_u_hat
                self.to_add['idx' + str(idx_h)] = self.U_idx_pool[max_idx]

    def _fit_and_evaluate(self, verbose):
        """
        Fits h1 and h2 and evaluates metrics.
        """
        self.h1.fit(self.L1_X, self.L1_y)
        self.h2.fit(self.L2_X, self.L2_y)
        self._evaluate_metrics(verbose)

    def _get_pool(self):
        """
        Gets unlabeled pool and indices of unlabeled.
        """
        self.U_X_pool, self.U_y_pool, self.U_idx_pool = shuffle(
            self.U_X, self.U_y, range(self.U_y.size))
        self.U_X_pool = self.U_X_pool[:self.pool_size]
        self.U_y_pool = self.U_y_pool[:self.pool_size]
        self.U_idx_pool = self.U_idx_pool[:self.pool_size]

    def _initialize_metrics(self):
        """
        Sets up metrics to be stored.
        """
        initial_metrics = np.full((self.num_trials, self.max_iters+1), np.inf)
        self.mses1_train = np.copy(initial_metrics)
        self.mses1_test = np.copy(initial_metrics)
        self.mses2_train = np.copy(initial_metrics)
        self.mses2_test = np.copy(initial_metrics)
        self.mses_train = np.copy(initial_metrics)
        self.mses_test = np.copy(initial_metrics)

    def _remove_from_unlabeled(self):
        # Remove added examples from unlabeled
        to_remove = []
        if self.to_add['idx1'] is not None:
            to_remove.append(self.to_add['idx1'])
        if self.to_add['idx2'] is not None:
            to_remove.append(self.to_add['idx2'])
        self.U_X = np.delete(self.U_X, to_remove, axis=0)
        self.U_y = np.delete(self.U_y, to_remove, axis=0)

    def _store_results(self, iteration):
        """
        Stores current MSEs.
        """
        self.mses1_train[self.trial,iteration] = self.mse1_train
        self.mses1_test[self.trial,iteration] = self.mse1_test
        self.mses2_train[self.trial,iteration] = self.mse2_train
        self.mses2_test[self.trial,iteration] = self.mse2_test
        self.mses_train[self.trial,iteration] = self.mse_train
        self.mses_test[self.trial,iteration] = self.mse_test

    