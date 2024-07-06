import numpy as np
import pandas as pd
import midasmlpy.date_functions as datef  # used to handle different frequencies of data and to create lags
import unittest


# load data from xlsx files and create a dataframe


class Test_TestDataTransformation(unittest.TestCase):
    def setUp(self):
        # Load data and set common variables
        self.Predictors = pd.read_excel('user_guide/predictors-monthly.xlsx').to_numpy()
        self.Target = pd.read_excel('user_guide/recessions-quarterly.xlsx').to_numpy()
        self.Y_date = self.Target[:, 0]
        self.Y = self.Target[:, 1]
        self.X_date = self.Predictors[:, 0]
        self.X = self.Predictors[:, 1:]

    def test_data_transform1(self):
        """Testing model"""
        # Expected output
        expected_Y = np.array([1, 0, 0, 0, 1, 1, 1, 1, 1, 0])
        expected_X_tilde = np.array([-0.60800157, -0.44238384, -0.0887135, 0.51130544, 0.38474615,
                                     1.17225514, 0.65501327, 0.46385537, 1.41947659,
                                     -0.6775583])  # test_1['X_tilde'][0:10,0]
        expected_Y_lagged = np.empty((174, 0))  # Adjust based on your function's processing and logic

        # Call the function
        output = datef.data_transform(self.Y, self.Y_date, self.X, self.X_date, x_lags=3, y_lags=0, horizon=0,
                                      degree=3, standardize=True)

        # Assertions to verify the output matches expected output
        np.testing.assert_array_equal(output['Y'][0:10], expected_Y)
        np.testing.assert_array_almost_equal(output['X_tilde'][0:10, 0], expected_X_tilde)
        np.testing.assert_array_equal(output['Y_lagged'], expected_Y_lagged)

    def test_data_transform2(self):
        """Testing y_lags"""
        # Expected output
        expected_Y = np.array([0, 0, 0, 1, 1, 1, 1, 1, 0, 0])
        expected_X_tilde = np.array([0.51130544, 0.38474615, 1.17225514, 0.58770945, 0.65501327,
                                     0.46385537, 1.41947659, 0.70855078, -0.6775583,
                                     -0.78863696])  # test_1['X_tilde'][0:10,0]
        expected_Y_lagged = np.array(
            [1, 0, 0, 0, 1, 1, 1, 1, 1, 0])  # Adjust based on your function's processing and logic

        # Call the function
        output = datef.data_transform(self.Y, self.Y_date, self.X, self.X_date, x_lags=3, y_lags=1, horizon=0,
                                      degree=4, standardize=True)

        # Assertions to verify the output matches expected output
        np.testing.assert_array_equal(output['Y'][0:10], expected_Y)
        np.testing.assert_array_almost_equal(output['X_tilde'][0:10, 0], expected_X_tilde)
        np.testing.assert_array_equal(output['Y_lagged'][0:10, 0], expected_Y_lagged)

    def test_data_transform3(self):
        """Testing horizon"""
        # Expected output
        expected_Y = np.array([1, 0, 0, 0, 1, 1, 1, 1, 1, 0])
        expected_X_tilde = np.array([-1.2609427, -0.41725772, -1.02339183, -0.63737169, 0.00529557,
                                     -1.2988975, -0.0322551, -1.9840987, 0.5949472,
                                     -0.8024663])  # test_1['X_tilde'][0:10,0]
        expected_Y_lagged = np.empty((174, 0))  # Adjust based on your function's processing and logic

        # Call the function
        output = datef.data_transform(self.Y, self.Y_date, self.X, self.X_date, x_lags=3, y_lags=0, horizon=2,
                                      degree=4, standardize=True)

        # Assertions to verify the output matches expected output
        np.testing.assert_array_equal(output['Y'][0:10], expected_Y)
        np.testing.assert_array_almost_equal(output['X_tilde'][0:10, 0], expected_X_tilde)
        np.testing.assert_array_equal(output['Y_lagged'], expected_Y_lagged)


if __name__ == '__main__':
    unittest.main()
